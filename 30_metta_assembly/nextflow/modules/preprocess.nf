// Preprocessing: optical deduplication, tile filtering, adapter trimming, artifact removal

process CLUMPIFY {
    tag "${meta.id}"
    label 'process_high'
    conda "${projectDir}/conda-envs/dana-metta-bbmap"
    publishDir "${params.outdir}/preprocess/${meta.id}", mode: 'copy', pattern: '*.fq.gz'
    storeDir params.store_dir ? "${params.store_dir}/preprocess/${meta.id}" : null

    input:
    tuple val(meta), path(r1), path(r2)

    output:
    tuple val(meta), path("${meta.id}.clumped.fq.gz"), emit: reads

    script:
    def xmx = task.memory ? "-Xmx${(task.memory.toGiga() * 0.85).intValue()}g" : ""
    """
    # Optical dedup requires Illumina tile coordinates in read headers.
    # SRA-downloaded data has stripped headers, causing clumpify to crash
    # on both `dedupe optical` and `dedupe` modes (IlluminaHeaderParser).
    # Clumpify can also hang after thread errors, so use timeout (5 min per GB input, min 120s).
    # Fall back to simple interleaving (reformat.sh) on failure.
    input_size=\$(stat -c%s "${r1}" 2>/dev/null || echo 0)
    timeout_sec=\$(( input_size / 1073741824 * 300 + 120 ))

    set +e
    timeout \${timeout_sec} clumpify.sh ${xmx} \\
        in1="${r1}" \\
        in2="${r2}" \\
        out="${meta.id}.clumped.fq.gz" \\
        dedupe optical \\
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
    label 'process_high'
    conda "${projectDir}/conda-envs/dana-metta-bbmap"
    publishDir "${params.outdir}/preprocess/${meta.id}", mode: 'copy', pattern: '*.fq.gz'
    storeDir params.store_dir ? "${params.store_dir}/preprocess/${meta.id}" : null

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}.filtered_by_tile.fq.gz"), emit: reads

    script:
    def xmx = task.memory ? "-Xmx${(task.memory.toGiga() * 0.85).intValue()}g" : ""
    """
    # filterbytile requires Illumina tile coordinates in read headers.
    # SRA-downloaded data has stripped headers — can crash and hang.
    # Use timeout to kill hung Java processes, then fall back to passthrough.
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
        echo "[WARNING] filterbytile failed (exit \$fbt_exit) — passing through (SRA headers?)" >&2
        cp "${reads}" "${meta.id}.filtered_by_tile.fq.gz"
    fi
    """
}

process BBDUK_TRIM {
    tag "${meta.id}"
    label 'process_high'
    conda "${projectDir}/conda-envs/dana-metta-bbmap"
    publishDir "${params.outdir}/preprocess/${meta.id}", mode: 'copy', pattern: '*.fq.gz'
    storeDir params.store_dir ? "${params.store_dir}/preprocess/${meta.id}" : null

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
    label 'process_high'
    conda "${projectDir}/conda-envs/dana-metta-bbmap"
    publishDir "${params.outdir}/preprocess/${meta.id}", mode: 'copy', pattern: '*.fq.gz'
    storeDir params.store_dir ? "${params.store_dir}/preprocess/${meta.id}" : null

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
