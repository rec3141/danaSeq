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
    """
    clumpify.sh \\
        in1="${r1}" \\
        in2="${r2}" \\
        out="${meta.id}.clumped.fq.gz" \\
        dedupe optical \\
        ow=t \\
        t=${task.cpus}

    if [ ! -s "${meta.id}.clumped.fq.gz" ]; then
        echo "[ERROR] Clumpify produced empty output for ${meta.id}" >&2
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
    """
    filterbytile.sh \\
        in="${reads}" \\
        out="${meta.id}.filtered_by_tile.fq.gz" \\
        ow=t \\
        t=${task.cpus}

    if [ ! -s "${meta.id}.filtered_by_tile.fq.gz" ]; then
        echo "[ERROR] filterbytile produced empty output for ${meta.id}" >&2
        exit 1
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
    """
    bbduk.sh \\
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
    """
    bbduk.sh \\
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
