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

// Per-sample auto-trim: each sample is profiled on its own reads and gets its own
// *_auto_trim.tsv. TRUNC_POLICY then chooses the group's value from these.
//
// Profile per sample, never over a group's reads pooled together: the profiler
// only reads the first `n_reads_sampled` reads it is given, so a pooled directory
// yields the first sample's answer rather than the group's.
process AUTO_TRIM {
    tag "${meta.id}"
    label 'process_low'
    conda params.lang == 'python' ? "${projectDir}/envs/python.yml" : "${projectDir}/envs/r.yml"
    publishDir "${params.outdir}/quality_check", mode: 'copy', pattern: "*_auto_trim.tsv"

    input:
    tuple val(meta), path(r1), path(r2)

    output:
    tuple val(meta), env(TRUNC_FWD), env(TRUNC_REV), env(READ_FWD), env(READ_REV), emit: trim_params
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
    export READ_FWD="\${LEN_FWD:-\$MIN_LEN}"
    export READ_REV="\${LEN_REV:-\$MIN_LEN}"
    export TRUNC_FWD=\$(clamp "\${RAW_FWD:-0}" "\$MIN_LEN" "\$READ_FWD")
    export TRUNC_REV=\$(clamp "\${RAW_REV:-0}" "\$MIN_LEN" "\$READ_REV")

    # Record the applied values (and whether the floor kicked in) in the tsv so
    # quality_check reflects what was actually used, not just what was detected.
    if [ "\$TRUNC_FWD" != "\$RAW_FWD" ] || [ "\$TRUNC_REV" != "\$RAW_REV" ]; then
        printf 'floored\\ttrue\\ntrunc_len_fwd_applied\\t%s\\ntrunc_len_rev_applied\\t%s\\n' "\$TRUNC_FWD" "\$TRUNC_REV" >> ${plate_id}_auto_trim.tsv
        echo "[WARN] Auto-trim ${plate_id}: floored trunc fwd \$RAW_FWD->\$TRUNC_FWD rev \$RAW_REV->\$TRUNC_REV (min \$MIN_LEN) — degraded library, expect loss at the quality filter"
    fi
    echo "[INFO] Auto-trim ${plate_id}: fwd=\$TRUNC_FWD rev=\$TRUNC_REV"
    """
}

// Collapse a group's per-sample truncation lengths into the single pair the group
// will use, by --trunc_policy, and record what was collapsed.
//
// The choice is a trade with a cliff on each side. DADA2 discards reads shorter
// than truncLen, so a long value drops reads from every sample whose quality
// falls off earlier. A short value keeps them but spends the overlap mergePairs
// needs: once fwd + rev drops below the amplicon length plus min_overlap, pairs
// stop merging and samples fall under min_reads instead. Both losses are silent
// in the final table, which is why the per-sample values and the resulting loss
// are written out here.
process TRUNC_POLICY {
    tag "${plate_id}"
    label 'process_low'
    // Echo stdout to the run log; without this the warnings below reach only
    // .command.out in the work dir.
    debug true
    publishDir "${params.outdir}/quality_check", mode: 'copy', pattern: "*_trunc_policy.tsv"

    input:
    tuple val(plate_id), val(sample_ids), val(fwds), val(revs), val(len_fwds), val(len_revs)
    path primer_files    // fwd/rev fasta when --primers_* were given, else empty

    output:
    tuple val(plate_id), env(TRUNC_FWD), env(TRUNC_REV), emit: trim_params
    path("${plate_id}_trunc_policy.tsv"), emit: policy_tsv

    script:
    def rows = [sample_ids, fwds, revs, len_fwds, len_revs].transpose()
                   .collect { s, f, r, lf, lr -> "${s}\t${f}\t${r}\t${lf}\t${lr}" }.join('\n')
    """
    cat > samples.tsv <<'SAMPLES_EOF'
${rows}
SAMPLES_EOF

    python3 <<'PYEOF'
policy  = "${params.trunc_policy}"
min_len = int("${params.auto_trim_min_length}")
plate   = "${plate_id}"
min_overlap = int("${params.min_overlap}")
min_seq_len = int("${params.min_seq_length}")
expected    = int("${params.expected_amplicon_length}")

# Fall back to the primer database when no length was given: expected_length in
# the vendored 16S table minus the two primers cutadapt removed. 18S and ITS
# tables carry no lengths, and de-novo inferred primers are in no table, so this
# resolves for some runs and not others — say which, never guess.
expected_src = "--expected_amplicon_length"
if expected <= 0:
    import glob as _glob, sys as _sys
    _sys.path.insert(0, "${projectDir}/bin")
    seqs = []
    for fa in sorted(_glob.glob("*.fa")) + sorted(_glob.glob("*.fasta")):
        cur = []
        for line in open(fa):
            line = line.strip()
            if line.startswith(">"):
                if cur: seqs.append("".join(cur)); cur = []
            elif line:
                cur.append(line.upper())
        if cur: seqs.append("".join(cur))
    try:
        from primer_db import insert_length_for
        for i in range(len(seqs)):
            for j in range(len(seqs)):
                if i == j: continue
                hit = insert_length_for(seqs[i], seqs[j])
                if hit:
                    expected, expected_src = hit, "primer database"
                    break
            if expected > 0: break
    except Exception as e:
        print("[INFO] Trunc policy: primer database unavailable (%s)" % e)

rows = []
for line in open("samples.tsv"):
    p = line.rstrip("\\n").split("\\t")
    if len(p) >= 5 and all(x.isdigit() for x in p[1:5]):
        rows.append((p[0], int(p[1]), int(p[2]), int(p[3]), int(p[4])))
if not rows:
    raise SystemExit("TRUNC_POLICY: no per-sample truncation values for " + plate)

# Percentiles of the group's per-sample truncation LENGTHS — not Phred scores.
PERCENTILES = {"p10": 0.10, "p25": 0.25, "p50": 0.50, "p75": 0.75}

def pick(vals):
    v = sorted(vals)
    if policy == "min":
        return v[0]
    if policy == "max":
        return v[-1]
    if policy not in PERCENTILES:
        hint = ""
        if policy.startswith("q") and policy[1:].isdigit():
            hint = (" — did you mean 'p" + policy[1:] + "'? These are percentiles of "
                    "read length, not Phred quality scores")
        elif policy == "median":
            hint = " — the median is 'p50'"
        raise SystemExit("TRUNC_POLICY: unknown --trunc_policy '" + policy +
                         "' (min, p10, p25, p50, p75, max)" + hint)
    i = PERCENTILES[policy] * (len(v) - 1)
    lo = int(i); hi = min(lo + 1, len(v) - 1)
    return int(round(v[lo] + (i - lo) * (v[hi] - v[lo])))

raw_fwd, raw_rev = pick([r[1] for r in rows]), pick([r[2] for r in rows])

# Same floor AUTO_TRIM applies per sample: a quality-driven short truncation
# leaves reads unable to overlap and the sample vanishes at DENOISE (issue #4).
# AUTO_TRIM also caps that floor at the read length, and so must this: a floor
# above the reads truncates every read out of existence, because dada2 discards
# reads shorter than truncLen. A 2x150 run trimmed to 148bp against a 150bp
# floor returns zero reads for every sample and reads as "no results".
# Cap by the same percentile over read lengths that picked the truncation, so
# the cap tracks the policy instead of letting one short library set it. Each
# sample's truncation is already <= its own read length, so this only ever
# constrains the floor -- never the quality-driven choice itself.
cap_fwd, cap_rev = pick([r[3] for r in rows]), pick([r[4] for r in rows])
fwd = min(max(raw_fwd, min_len), cap_fwd)
rev = min(max(raw_rev, min_len), cap_rev)

# How many samples are being truncated PAST their own quality cliff, and so will
# lose reads purely because of the group they landed in.
past = [r[0] for r in rows if r[1] < fwd or r[2] < rev]

# Truncating past a sample's own READ LENGTH is a different failure: it is not a
# trade, it is a guaranteed zero for that sample. Fail loudly when it would take
# the whole group, rather than emitting an empty table that looks like a primer
# or metadata problem downstream.
past_len = [r[0] for r in rows if r[3] < fwd or r[4] < rev]
if len(past_len) == len(rows):
    raise SystemExit(
        "TRUNC_POLICY: %s truncates fwd=%d rev=%d past the read length of all %d "
        "sample(s) (longest fwd=%d rev=%d) — every read would be discarded at the "
        "quality filter. Lower --auto_trim_min_length (currently %d) or set "
        "--truncLen_fwd/--truncLen_rev explicitly."
        % (plate, fwd, rev, len(rows), max(r[3] for r in rows),
           max(r[4] for r in rows), min_len))

# Can a merged fragment still exist? mergePairs needs fwd + rev to cover the
# fragment plus min_overlap. Too little does not raise — pairs quietly fail to
# merge, samples fall under min_reads, and the run returns a smaller table that
# looks fine. Check it here, before denoising, not by counting missing samples
# afterwards.
span = fwd + rev
floor_need = min_seq_len + min_overlap      # structural: shorter cannot merge at all
frag_need = (expected + min_overlap) if expected > 0 else 0
frag_warn = ""
if span < floor_need:
    frag_warn = ("span %d < min_seq_length %d + min_overlap %d = %d: no pair can merge "
                 "into a sequence this pipeline would keep"
                 % (span, min_seq_len, min_overlap, floor_need))
elif frag_need and span < frag_need:
    frag_warn = ("span %d < expected_amplicon_length %d + min_overlap %d = %d: pairs "
                 "will fail to merge and samples will drop out under min_reads"
                 % (span, expected, min_overlap, frag_need))

with open(plate + "_trunc_policy.tsv", "w") as out:
    out.write("key\\tvalue\\n")
    out.write("policy\\t%s\\n" % policy)
    out.write("n_samples\\t%d\\n" % len(rows))
    out.write("trunc_len_fwd_applied\\t%d\\n" % fwd)
    out.write("trunc_len_rev_applied\\t%d\\n" % rev)
    out.write("floored\\t%s\\n" % str(fwd != raw_fwd or rev != raw_rev).lower())
    out.write("read_len_cap_fwd\\t%d\\n" % cap_fwd)
    out.write("read_len_cap_rev\\t%d\\n" % cap_rev)
    out.write("capped_by_read_len\\t%s\\n"
              % str(max(raw_fwd, min_len) > cap_fwd or max(raw_rev, min_len) > cap_rev).lower())
    out.write("samples_truncated_past_own\\t%d\\n" % len(past))
    out.write("samples_truncated_past_read_len\\t%d\\n" % len(past_len))
    out.write("span_fwd_plus_rev\\t%d\\n" % span)
    out.write("span_required\\t%d\\n" % max(floor_need, frag_need))
    out.write("span_sufficient\\t%s\\n" % str(not frag_warn).lower())
    out.write("expected_amplicon_length\\t%s\\n" % (expected if expected > 0 else ""))
    out.write("expected_amplicon_source\\t%s\\n" % (expected_src if expected > 0 else "unknown"))
    out.write("#\\tper-sample values follow\\n")
    out.write("#sample\\ttrunc_fwd\\ttrunc_rev\\tread_len_fwd\\tread_len_rev"
              "\\ttruncated_past_own\\ttruncated_past_read_len\\n")
    for s, f, r, lf, lr in sorted(rows):
        out.write("%s\\t%d\\t%d\\t%d\\t%d\\t%s\\t%s\\n"
                  % (s, f, r, lf, lr, str(f < fwd or r < rev).lower(),
                     str(lf < fwd or lr < rev).lower()))

print("[INFO] Trunc policy %s for %s: fwd=%d rev=%d from %d samples (span %d)" %
      (policy, plate, fwd, rev, len(rows), span))
if frag_warn:
    print("[WARN] Trunc policy %s for %s: %s" % (policy, plate, frag_warn))
elif expected > 0:
    print("[INFO] Trunc policy %s for %s: span %d covers the %dbp fragment (%s) "
          "plus %dbp overlap" % (policy, plate, span, expected, expected_src, min_overlap))
else:
    print("[INFO] Trunc policy %s for %s: overlap unchecked — no amplicon length for "
          "these primers (set --expected_amplicon_length)" % (policy, plate))
if max(raw_fwd, min_len) > cap_fwd or max(raw_rev, min_len) > cap_rev:
    print("[INFO] Trunc policy %s for %s: floor --auto_trim_min_length %d capped to "
          "the reads actually present (fwd=%d rev=%d)"
          % (policy, plate, min_len, cap_fwd, cap_rev))
if past:
    print("[WARN] Trunc policy %s for %s: %d of %d samples truncated past their own "
          "quality cliff and will lose reads at the filter: %s"
          % (policy, plate, len(past), len(rows), ", ".join(sorted(past)[:10])
             + (" ..." if len(past) > 10 else "")))
if past_len:
    print("[WARN] Trunc policy %s for %s: %d of %d samples are shorter than fwd=%d "
          "rev=%d and will return ZERO reads at the filter: %s"
          % (policy, plate, len(past_len), len(rows), fwd, rev,
             ", ".join(sorted(past_len)[:10]) + (" ..." if len(past_len) > 10 else "")))
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
