// Preprocessing: optical deduplication, tile filtering, adapter trimming, artifact removal,
// human decontamination, and QC reporting

// ── Input-size-scaled resources ──────────────────────────────────────────────
// The per-sample preprocessing steps run under Nextflow's local executor inside
// one SLURM allocation, so cpus/memory here govern how many samples run at once.
// Reserving the full assembly allocation (24 CPU / 250 GB) for a 20 MB SRA file
// serialises the whole preprocessing stage and hands 100+ GB heaps to trivial
// jobs. Scale to the input instead: small base + growth per GB, capped at the
// assembly allocation so multi-GB samples still get real resources.
def _sizeGb(f)  { (f instanceof Collection ? (f.sum { it.size() } ?: 0L) : f.size()) / 1073741824.0 }
def cpusFor(f, base, perGb, cap) { Math.min((base + _sizeGb(f) * perGb) as int, cap as int) }
def memFor(f, base, perGb)       { def gb = (base + Math.ceil(_sizeGb(f) * perGb)) as int; "${gb} GB" }

process CLUMPIFY {
    tag "${meta.id}"
    cpus   { cpusFor([r1, r2], 4, 2, params.assembly_cpus) }
    memory { memFor([r1, r2], 8, 16) }
    conda "${projectDir}/conda-envs/dana-illumina-mag-bbmap"
    publishDir { "${params.outdir}/preprocess/${meta.id}" }, mode: 'copy', pattern: '*.fq.gz'
    storeDir { params.store_dir ? "${params.store_dir}/preprocess/${meta.id}" : null }

    input:
    tuple val(meta), path(r1), path(r2)

    output:
    tuple val(meta), path("${meta.id}.clumped.fq.gz"), emit: reads

    script:
    def xmx = task.memory ? "-Xmx${(task.memory.toGiga() * 0.85).intValue()}g" : ""
    """
    # Optical dedup requires native Illumina flowcell coordinates in the read
    # header. SRA reformats headers (accession.N prefix), which clumpify's
    # optical-coordinate parser chokes on and then hangs (IlluminaHeaderParser2
    # AssertionError in fetch threads → main thread waits until the timeout).
    # Detect the header style from the first read: a native Illumina header's
    # first whitespace-token is instr:run:fc:lane:tile:x:y (>= 6 colons); SRA's
    # is the bare accession (0 colons). Use optical only when native, else plain
    # dedupe (fast, correct, no crash). timeout remains a safety net.
    first_tok=\$(zcat "${r1}" 2>/dev/null | head -1 | awk '{print \$1}')
    ncolons=\$(printf '%s' "\$first_tok" | tr -cd ':' | wc -c | tr -d ' ')
    if [ "\${ncolons:-0}" -ge 6 ]; then
        dedupe_opt="optical"
        echo "[INFO] Native Illumina headers detected — clumpify optical dedup" >&2
    else
        dedupe_opt=""
        echo "[INFO] Non-native/SRA headers — clumpify plain dedupe (no optical)" >&2
    fi

    input_size=\$(stat -c%s "${r1}" 2>/dev/null || echo 0)
    timeout_sec=\$(( input_size / 1073741824 * 300 + 120 ))

    set +e
    timeout \${timeout_sec} clumpify.sh ${xmx} \\
        in1="${r1}" \\
        in2="${r2}" \\
        out="${meta.id}.clumped.fq.gz" \\
        dedupe \$dedupe_opt \\
        ow=t \\
        t=${task.cpus}
    clump_exit=\$?
    set -e

    if [ \$clump_exit -ne 0 ] || [ ! -s "${meta.id}.clumped.fq.gz" ]; then
        echo "[WARNING] Optical dedup failed (exit \$clump_exit) — falling back to interleave-only (SRA headers?)" >&2
        reformat.sh ${xmx} \\
            in1="${r1}" \\
            in2="${r2}" \\
            out="${meta.id}.clumped.fq.gz" \\
            ow=t
    fi

    if [ ! -s "${meta.id}.clumped.fq.gz" ]; then
        echo "[ERROR] Clumpify/reformat produced empty output for ${meta.id}" >&2
        exit 1
    fi
    """
}

process FILTER_BY_TILE {
    tag "${meta.id}"
    cpus   { cpusFor(reads, 4, 2, params.assembly_cpus) }
    memory { memFor(reads, 8, 16) }
    conda "${projectDir}/conda-envs/dana-illumina-mag-bbmap"
    publishDir { "${params.outdir}/preprocess/${meta.id}" }, mode: 'copy', pattern: '*.fq.gz'
    storeDir { params.store_dir ? "${params.store_dir}/preprocess/${meta.id}" : null }

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}.filtered_by_tile.fq.gz"), emit: reads

    script:
    def xmx = task.memory ? "-Xmx${(task.memory.toGiga() * 0.85).intValue()}g" : ""
    """
    # filterbytile requires native Illumina tile coordinates in read headers.
    # SRA-reformatted headers lack a parseable coordinate token in the first
    # field, so filterbytile crashes/hangs the same way clumpify optical does.
    # Detect header style (native = first token has >= 6 colons) and only run
    # filterbytile on native data; otherwise pass through immediately instead of
    # burning the timeout on a doomed attempt. timeout remains a safety net.
    first_tok=\$(zcat "${reads}" 2>/dev/null | head -1 | awk '{print \$1}')
    ncolons=\$(printf '%s' "\$first_tok" | tr -cd ':' | wc -c | tr -d ' ')

    if [ "\${ncolons:-0}" -ge 6 ]; then
        input_size=\$(stat -c%s "${reads}" 2>/dev/null || echo 0)
        timeout_sec=\$(( input_size / 1073741824 * 300 + 120 ))
        set +e
        timeout \${timeout_sec} filterbytile.sh ${xmx} \\
            in="${reads}" \\
            out="${meta.id}.filtered_by_tile.fq.gz" \\
            ow=t \\
            t=${task.cpus}
        fbt_exit=\$?
        set -e
        if [ \$fbt_exit -ne 0 ] || [ ! -s "${meta.id}.filtered_by_tile.fq.gz" ]; then
            echo "[WARNING] filterbytile failed (exit \$fbt_exit) — passing through" >&2
            cp "${reads}" "${meta.id}.filtered_by_tile.fq.gz"
        fi
    else
        echo "[INFO] Non-native/SRA headers — skipping tile filter (passthrough)" >&2
        cp "${reads}" "${meta.id}.filtered_by_tile.fq.gz"
    fi
    """
}

process BBDUK_TRIM {
    tag "${meta.id}"
    cpus   { cpusFor(reads, 4, 2, params.assembly_cpus) }
    memory { memFor(reads, 8, 8) }
    conda "${projectDir}/conda-envs/dana-illumina-mag-bbmap"
    publishDir { "${params.outdir}/preprocess/${meta.id}" }, mode: 'copy', pattern: '*.fq.gz'
    storeDir { params.store_dir ? "${params.store_dir}/preprocess/${meta.id}" : null }

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}.trimmed.fq.gz"), emit: reads

    script:
    def xmx = task.memory ? "-Xmx${(task.memory.toGiga() * 0.85).intValue()}g" : ""
    """
    bbduk.sh ${xmx} \\
        in="${reads}" \\
        out="${meta.id}.trimmed.fq.gz" \\
        ktrim=r k=23 mink=11 hdist=1 \\
        tbo tpe \\
        minlen=${params.min_readlen} \\
        ref=adapters \\
        ftm=5 ordered \\
        ow=t \\
        t=${task.cpus}

    if [ ! -s "${meta.id}.trimmed.fq.gz" ]; then
        echo "[ERROR] bbduk trim produced empty output for ${meta.id}" >&2
        exit 1
    fi
    """
}

process BBDUK_FILTER {
    tag "${meta.id}"
    cpus   { cpusFor(reads, 4, 2, params.assembly_cpus) }
    memory { memFor(reads, 8, 8) }
    conda "${projectDir}/conda-envs/dana-illumina-mag-bbmap"
    publishDir { "${params.outdir}/preprocess/${meta.id}" }, mode: 'copy', pattern: '*.fq.gz'
    storeDir { params.store_dir ? "${params.store_dir}/preprocess/${meta.id}" : null }

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}.filtered.fq.gz"), emit: reads

    script:
    def xmx = task.memory ? "-Xmx${(task.memory.toGiga() * 0.85).intValue()}g" : ""
    """
    bbduk.sh ${xmx} \\
        in="${reads}" \\
        out="${meta.id}.filtered.fq.gz" \\
        k=31 \\
        ref=artifacts,phix \\
        entropy=0.95 \\
        ordered cardinality \\
        ow=t \\
        t=${task.cpus}

    if [ ! -s "${meta.id}.filtered.fq.gz" ]; then
        echo "[ERROR] bbduk filter produced empty output for ${meta.id}" >&2
        exit 1
    fi
    """
}

process REMOVE_HUMAN {
    tag "${meta.id}"
    // Higher memory floor than the other steps: bbmap must load the masked
    // human reference index into memory before mapping.
    cpus   { cpusFor(reads, 8, 2, params.assembly_cpus) }
    memory { memFor(reads, 32, 8) }
    conda "${projectDir}/conda-envs/dana-illumina-mag-bbmap"
    publishDir { "${params.outdir}/preprocess/${meta.id}" }, mode: 'copy', pattern: '*.fq.gz'
    storeDir { params.store_dir ? "${params.store_dir}/preprocess/${meta.id}" : null }

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}.nohuman.fq.gz"), emit: reads
    path("${meta.id}.human_contam.fq.gz"),             emit: contam, optional: true

    script:
    def xmx = task.memory ? "-Xmx${(task.memory.toGiga() * 0.85).intValue()}g" : ""
    """
    # removehuman.sh uses bbmap internally against a masked human reference.
    # Input is interleaved; output keeps interleaved format.
    removehuman.sh ${xmx} \\
        in="${reads}" \\
        outu="${meta.id}.nohuman.fq.gz" \\
        outm="${meta.id}.human_contam.fq.gz" \\
        path="${params.human_ref}" \\
        t=${task.cpus}

    if [ ! -s "${meta.id}.nohuman.fq.gz" ]; then
        echo "[ERROR] Human removal produced empty output for ${meta.id}" >&2
        exit 1
    fi
    """
}

process FASTQC {
    tag "${meta.id}"
    label 'process_medium'
    conda "${projectDir}/conda-envs/dana-illumina-mag-bbmap"
    publishDir { "${params.outdir}/preprocess/${meta.id}/fastqc" }, mode: 'copy'
    storeDir { params.store_dir ? "${params.store_dir}/preprocess/${meta.id}/fastqc" : null }

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.html"), path("*.zip"), emit: reports

    script:
    """
    # FastQC expects separate files, not interleaved — deinterleave first
    reformat.sh in="${reads}" out1=r1.fq.gz out2=r2.fq.gz ow=t

    fastqc -t ${task.cpus} --noextract -o . r1.fq.gz r2.fq.gz

    # Rename outputs with sample ID prefix
    for f in *.html *.zip; do
        mv "\$f" "${meta.id}.\$f"
    done
    """
}
