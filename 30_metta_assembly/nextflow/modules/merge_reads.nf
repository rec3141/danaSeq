// Read merging: merge overlapping paired reads, quality-trim unmerged

process MERGE_READS {
    tag "${meta.id}"
    label 'process_high'
    conda "${projectDir}/conda-envs/dana-metta-bbmap"
    publishDir "${params.outdir}/merge/${meta.id}", mode: 'copy', pattern: '*.fq.gz'
    storeDir params.store_dir ? "${params.store_dir}/merge/${meta.id}" : null

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}.merged.fq.gz"),   emit: merged
    tuple val(meta), path("${meta.id}.unmerged.fq.gz"), emit: unmerged

    script:
    def xmx = task.memory ? "-Xmx${(task.memory.toGiga() * 0.85).intValue()}g" : ""
    """
    bbmerge-auto.sh ${xmx} \\
        in="${reads}" \\
        out="${meta.id}.merged.fq.gz" \\
        outu="${meta.id}.unmerged.fq.gz" \\
        strict k=93 extend2=80 rem ordered \\
        ihist=ihist_merge.txt \\
        prefilter=2 \\
        ow=t \\
        t=${task.cpus}

    if [ ! -s "${meta.id}.merged.fq.gz" ]; then
        echo "[WARNING] No merged reads for ${meta.id}" >&2
        touch "${meta.id}.merged.fq.gz"
    fi
    if [ ! -s "${meta.id}.unmerged.fq.gz" ]; then
        echo "[WARNING] No unmerged reads for ${meta.id}" >&2
        touch "${meta.id}.unmerged.fq.gz"
    fi
    """
}

process QUALITY_TRIM {
    tag "${meta.id}"
    label 'process_medium'
    conda "${projectDir}/conda-envs/dana-metta-bbmap"
    publishDir "${params.outdir}/merge/${meta.id}", mode: 'copy', pattern: '*.fq.gz'
    storeDir params.store_dir ? "${params.store_dir}/merge/${meta.id}" : null

    input:
    tuple val(meta), path(unmerged)

    output:
    tuple val(meta), path("${meta.id}.qtrimmed.fq.gz"), emit: reads

    script:
    def xmx = task.memory ? "-Xmx${(task.memory.toGiga() * 0.85).intValue()}g" : ""
    """
    bbduk.sh ${xmx} \\
        in="${unmerged}" \\
        out="${meta.id}.qtrimmed.fq.gz" \\
        qtrim=r trimq=10 \\
        minlen=${params.min_readlen} \\
        ordered \\
        ow=t \\
        t=${task.cpus}

    if [ ! -s "${meta.id}.qtrimmed.fq.gz" ]; then
        echo "[WARNING] Quality trimming produced empty output for ${meta.id}" >&2
        touch "${meta.id}.qtrimmed.fq.gz"
    fi
    """
}
