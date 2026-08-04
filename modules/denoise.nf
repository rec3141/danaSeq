// Amplicon denoising processes.
//
// Three-step workflow:
//   1. FILTER_TRIM    — per-sample quality filtering (uses --lang)
//   2. LEARN_ERRORS   — per-plate error model learning (uses --denoise_engine)
//   3. DENOISE        — per-plate denoising + pair merging (uses --denoise_engine)
//
// --denoise_engine controls which implementation to use:
//   'python' (default) = papa2 (byte-identical to R dada2, no R dependency)
//   'R' = R dada2 package (reference implementation)
//
// --lang controls filter_trim and downstream scripts independently.

// Resolve effective denoise engine: explicit param > lang fallback
def denoiseEngine() { return params.denoise_engine ?: params.lang }

// Output extensions: denoise steps always output .rds when engine=R
def errExt() { return denoiseEngine() == 'python' ? 'pkl' : 'rds' }
def seqExt() { return params.lang == 'python' ? 'pkl' : 'rds' }

// PER-SAMPLE auto-trim. Every sample is profiled on its own reads; the group's
// value is then chosen from those by TRUNC_POLICY.
//
// This used to run once per plate over every read in the group pooled into one
// directory. That did not produce a group consensus — the profiler samples the
// first `n_reads_sampled` reads it encounters, which all come from whichever
// sample sorts first, so the group's truncation was set by ONE arbitrary sample
// chosen by filename order. Measured on 492f42d0: the group value equalled the
// alphabetically-first sample's value in 12 of 14 groups, versus 6 of 14 for the
// group max and 6 of 14 for the median. On 2d9fa6af it looked like a deliberate
// "take the max" only because the first-sorting sample also happened to have the
// highest value — a coincidence that hid the real behaviour.
//
// Profiling per sample also restores the per-sample *_auto_trim.tsv files, which
// silently stopped existing when this went per-plate, taking with them the only
// record of what each sample would have chosen for itself.
process AUTO_TRIM {
    tag "${meta.id}"
    label 'process_low'
    conda params.lang == 'python' ? "${projectDir}/envs/python.yml" : "${projectDir}/envs/r.yml"
    publishDir "${params.outdir}/quality_check", mode: 'copy', pattern: "*_auto_trim.tsv"

    input:
    tuple val(meta), path(r1), path(r2)

    output:
    tuple val(meta), env(TRUNC_FWD), env(TRUNC_REV), emit: trim_params
    path("${meta.id}_auto_trim.tsv"), emit: params_tsv

    script:
    def plate_id = meta.id
    """
    microscape auto-trim "." \
        --min-quality ${params.auto_trim_min_quality} \
        --min-length ${params.auto_trim_min_length} \
        --output ${plate_id}_auto_trim.tsv \
        --verbose

    # Read the auto-detected values (exact key match, before we append anything).
    RAW_FWD=\$(awk -F'\\t' '\$1=="trunc_len_fwd"{print \$2}' ${plate_id}_auto_trim.tsv)
    RAW_REV=\$(awk -F'\\t' '\$1=="trunc_len_rev"{print \$2}' ${plate_id}_auto_trim.tsv)
    LEN_FWD=\$(awk -F'\\t' '\$1=="fwd_read_len"{print \$2}' ${plate_id}_auto_trim.tsv)
    LEN_REV=\$(awk -F'\\t' '\$1=="rev_read_len"{print \$2}' ${plate_id}_auto_trim.tsv)
    MIN_LEN=${params.auto_trim_min_length}

    # Enforce a minimum truncation length. A quality-driven short truncation
    # (e.g. ~30bp on a degraded library) leaves reads unable to overlap, so
    # mergePairs discards ~100% of reads and the sample silently vanishes at
    # DENOISE. Floor at MIN_LEN, capped at the actual read length. Issue #4.
    clamp() { v=\$1; lo=\$2; hi=\$3; [ "\$v" -lt "\$lo" ] && v=\$lo; [ "\$v" -gt "\$hi" ] && v=\$hi; echo "\$v"; }
    export TRUNC_FWD=\$(clamp "\${RAW_FWD:-0}" "\$MIN_LEN" "\${LEN_FWD:-\$MIN_LEN}")
    export TRUNC_REV=\$(clamp "\${RAW_REV:-0}" "\$MIN_LEN" "\${LEN_REV:-\$MIN_LEN}")

    # Record the applied values (and whether the floor kicked in) in the tsv so
    # quality_check reflects what was actually used, not just what was detected.
    if [ "\$TRUNC_FWD" != "\$RAW_FWD" ] || [ "\$TRUNC_REV" != "\$RAW_REV" ]; then
        printf 'floored\\ttrue\\ntrunc_len_fwd_applied\\t%s\\ntrunc_len_rev_applied\\t%s\\n' "\$TRUNC_FWD" "\$TRUNC_REV" >> ${plate_id}_auto_trim.tsv
        echo "[WARN] Auto-trim ${plate_id}: floored trunc fwd \$RAW_FWD->\$TRUNC_FWD rev \$RAW_REV->\$TRUNC_REV (min \$MIN_LEN) — degraded library, expect loss at the quality filter"
    fi
    echo "[INFO] Auto-trim ${plate_id}: fwd=\$TRUNC_FWD rev=\$TRUNC_REV"
    """
}

// Collapse a group's per-sample truncation lengths into the one pair the group
// will actually use, by an EXPLICIT policy, and record what was collapsed.
//
// There is no free choice here, only a trade. DADA2 discards any read shorter
// than truncLen, so a long value silently drops reads from every sample whose
// quality falls off earlier; a short one keeps them but spends overlap, and if
// overlap runs out mergePairs discards ~100% and the sample vanishes instead.
// Measured on 2d9fa6af (one group, 12 samples): min recovered +51% final reads
// and +26% ASVs over the longest value, at no cost in ASV length. Measured on
// 492f42d0 (14 groups): a few degraded samples bottom out near 22bp, so `min`
// there is set by the worst library in the group — which is why the default is a
// low quantile rather than the minimum, and why the floor below still applies.
//
// Whatever is chosen, the per-sample values and the resulting loss are written
// out. A truncation that discards half a sample's reads should never be silent.
process TRUNC_POLICY {
    tag "${plate_id}"
    label 'process_low'
    // Echo this process's stdout to the run log. Without it Nextflow files process
    // output in .command.out inside the work dir, where nobody looks — which is
    // why AUTO_TRIM's existing [WARN] lines have never appeared in a single run
    // log. A warning that a group is discarding half a sample's reads is worth
    // nothing if it is only readable by someone who already suspects it.
    debug true
    publishDir "${params.outdir}/quality_check", mode: 'copy', pattern: "*_trunc_policy.tsv"

    input:
    tuple val(plate_id), val(sample_ids), val(fwds), val(revs)

    output:
    tuple val(plate_id), env(TRUNC_FWD), env(TRUNC_REV), emit: trim_params
    path("${plate_id}_trunc_policy.tsv"), emit: policy_tsv

    script:
    def rows = [sample_ids, fwds, revs].transpose()
                   .collect { s, f, r -> "${s}\t${f}\t${r}" }.join('\n')
    """
    cat > samples.tsv <<'SAMPLES_EOF'
${rows}
SAMPLES_EOF

    python3 <<'PYEOF'
policy  = "${params.trunc_policy}"
min_len = int("${params.auto_trim_min_length}")
plate   = "${plate_id}"

rows = []
for line in open("samples.tsv"):
    p = line.rstrip("\\n").split("\\t")
    if len(p) >= 3 and p[1].isdigit() and p[2].isdigit():
        rows.append((p[0], int(p[1]), int(p[2])))
if not rows:
    raise SystemExit("TRUNC_POLICY: no per-sample truncation values for " + plate)

QUANTILES = {"q10": 0.10, "q25": 0.25, "median": 0.50, "q75": 0.75}

def pick(vals):
    v = sorted(vals)
    if policy == "min":
        return v[0]
    if policy == "max":
        return v[-1]
    if policy not in QUANTILES:
        raise SystemExit("TRUNC_POLICY: unknown --trunc_policy '" + policy +
                         "' (min, q10, q25, median, q75, max)")
    i = QUANTILES[policy] * (len(v) - 1)
    lo = int(i); hi = min(lo + 1, len(v) - 1)
    return int(round(v[lo] + (i - lo) * (v[hi] - v[lo])))

raw_fwd, raw_rev = pick([r[1] for r in rows]), pick([r[2] for r in rows])
# Same floor AUTO_TRIM applies per sample: a quality-driven short truncation
# leaves reads unable to overlap and the sample vanishes at DENOISE (issue #4).
fwd, rev = max(raw_fwd, min_len), max(raw_rev, min_len)

# How many samples are being truncated PAST their own quality cliff, and so will
# lose reads purely because of the group they landed in.
past = [r[0] for r in rows if r[1] < fwd or r[2] < rev]

with open(plate + "_trunc_policy.tsv", "w") as out:
    out.write("key\\tvalue\\n")
    out.write("policy\\t%s\\n" % policy)
    out.write("n_samples\\t%d\\n" % len(rows))
    out.write("trunc_len_fwd_applied\\t%d\\n" % fwd)
    out.write("trunc_len_rev_applied\\t%d\\n" % rev)
    out.write("floored\\t%s\\n" % str(fwd != raw_fwd or rev != raw_rev).lower())
    out.write("samples_truncated_past_own\\t%d\\n" % len(past))
    out.write("#\\tper-sample values follow\\n")
    out.write("#sample\\ttrunc_fwd\\ttrunc_rev\\ttruncated_past_own\\n")
    for s, f, r in sorted(rows):
        out.write("%s\\t%d\\t%d\\t%s\\n" % (s, f, r, str(f < fwd or r < rev).lower()))

print("[INFO] Trunc policy %s for %s: fwd=%d rev=%d from %d samples" %
      (policy, plate, fwd, rev, len(rows)))
if past:
    print("[WARN] Trunc policy %s for %s: %d of %d samples truncated past their own "
          "quality cliff and will lose reads at the filter: %s"
          % (policy, plate, len(past), len(rows), ", ".join(sorted(past)[:10])
             + (" ..." if len(past) > 10 else "")))
PYEOF

    export TRUNC_FWD=\$(awk -F'\\t' '\$1=="trunc_len_fwd_applied"{print \$2}' ${plate_id}_trunc_policy.tsv)
    export TRUNC_REV=\$(awk -F'\\t' '\$1=="trunc_len_rev_applied"{print \$2}' ${plate_id}_trunc_policy.tsv)
    """
}

process FILTER_TRIM {
    tag "${meta.id}"
    label 'process_low'
    conda params.lang == 'python' ? "${projectDir}/envs/python.yml" : "${projectDir}/envs/r.yml"
    publishDir "${params.outdir}/filtered", mode: 'copy', pattern: "*_filt_stats.tsv", enabled: !params.store_dir

    input:
    tuple val(meta), path(r1), path(r2), val(trunc_fwd), val(trunc_rev)

    output:
    tuple val(meta), path("${meta.id}_R1.filt.fastq.gz"), path("${meta.id}_R2.filt.fastq.gz"), emit: reads
    path("${meta.id}_filt_stats.tsv"), emit: stats

    script:
    def eff_trunc_fwd = params.truncLen_fwd > 0 ? params.truncLen_fwd : trunc_fwd
    def eff_trunc_rev = params.truncLen_rev > 0 ? params.truncLen_rev : trunc_rev
    if (params.lang == 'python')
    """
    papa2 filter-trim \
        "${r1}" "${meta.id}_R1.filt.fastq.gz" \
        "${r2}" "${meta.id}_R2.filt.fastq.gz" \
        --max-ee ${params.maxEE} \
        --trunc-q ${params.truncQ} \
        --max-n ${params.maxN} \
        --trunc-len-fwd ${eff_trunc_fwd} \
        --trunc-len-rev ${eff_trunc_rev} \
        --stats "${meta.id}_filt_stats.tsv" \
        --sample-id "${meta.id}"
    """
    else
    """
    dada2_filter_trim.R \
        "${meta.id}" "${r1}" "${r2}" \
        ${params.maxEE} ${params.truncQ} ${params.maxN} \
        ${params.truncLen_fwd} ${params.truncLen_rev} \
        ${task.cpus}
    """
}

// Learn error rates — uses denoise_engine (R or python/papa2)
process LEARN_ERRORS {
    tag "${meta.id}"
    label 'process_high'
    conda ((params.denoise_engine ?: params.lang) == 'python' ? "${projectDir}/envs/python.yml" : "${projectDir}/envs/r.yml")
    publishDir "${params.outdir}/error_models", mode: 'copy', enabled: !params.store_dir
    storeDir params.store_dir ? "${params.store_dir}/error_models" : null

    input:
    tuple val(meta), path(r1_files), path(r2_files)

    output:
    tuple val(meta.plate), val(meta), path("${meta.id}_errF.${errExt()}"), path("${meta.id}_errR.${errExt()}"), emit: error_models
    path("${meta.id}_error_rates.pdf"), emit: error_plots

    script:
    if (denoiseEngine() == 'python')
    """
    learn_errors.py "${meta.id}" ${task.cpus}
    """
    else
    """
    dada2_learn_errors.R "${meta.id}" ${task.cpus}
    """
}

// Per-plate denoising — uses denoise_engine (R or python/papa2)
// Output is always .pkl when lang=python (converted from .rds by wrapper if needed)
process DENOISE {
    tag "${meta.id}"
    label 'process_high'
    conda ((params.denoise_engine ?: params.lang) == 'python' ? "${projectDir}/envs/python.yml" : "${projectDir}/envs/r.yml")
    publishDir "${params.outdir}/seqtabs", mode: 'copy', enabled: !params.store_dir
    storeDir params.store_dir ? "${params.store_dir}/seqtabs" : null

    input:
    tuple val(meta), path(r1_files), path(r2_files), path(errF), path(errR)

    output:
    tuple val(meta), path("${meta.id}.seqtab.${seqExt()}"), emit: seqtab
    path("${meta.id}.seqtab.tsv"), emit: seqtab_tsv

    script:
    if (denoiseEngine() == 'python')
    """
    denoise.py \
        "${meta.id}" "${errF}" "${errR}" \
        ${params.min_overlap} ${task.cpus}
    """
    else if (params.lang == 'python')
    // R dada2 engine but Python downstream: run R, then convert .rds to .pkl
    """
    dada2_denoise.R \
        "${meta.id}" "${errF}" "${errR}" \
        ${params.min_overlap} ${task.cpus}

    rds_to_pkl.py "${meta.id}.seqtab.rds" "${meta.id}.seqtab.pkl"
    """
    else
    """
    dada2_denoise.R \
        "${meta.id}" "${errF}" "${errR}" \
        ${params.min_overlap} ${task.cpus}
    """
}
