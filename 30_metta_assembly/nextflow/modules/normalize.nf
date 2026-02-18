// Normalization: coverage normalization with bbnorm

process NORMALIZE_READS {
    tag "${meta.id}"
    label 'process_high'
    conda "${projectDir}/conda-envs/dana-metta-bbmap"
    publishDir "${params.outdir}/normalize/${meta.id}", mode: 'copy', pattern: '*.{fq.gz,txt}'
    storeDir params.store_dir ? "${params.store_dir}/normalize/${meta.id}" : null

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}.normalized.fq.gz"), emit: reads
    path("${meta.id}.khist.txt"),                         emit: khist
    path("${meta.id}.peaks.txt"),                         emit: peaks

    script:
    """
    bbnorm.sh \\
        in="${reads}" \\
        out="${meta.id}.normalized.fq.gz" \\
        target=100 mindepth=2 \\
        prefilter=t \\
        hist="${meta.id}.khist.txt" \\
        peaks="${meta.id}.peaks.txt" \\
        ow=t \\
        t=${task.cpus}

    if [ ! -s "${meta.id}.normalized.fq.gz" ]; then
        echo "[ERROR] Normalization produced empty output for ${meta.id}" >&2
        exit 1
    fi
    """
}
