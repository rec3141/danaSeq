// Error correction: three-phase correction (overlap, clump, kmer)

process ERROR_CORRECT_ECCO {
    tag "${meta.id}"
    label 'process_high'
    conda "${projectDir}/conda-envs/dana-metta-bbmap"
    publishDir "${params.outdir}/error_correct/${meta.id}", mode: 'copy', pattern: '*.fq.gz'
    storeDir params.store_dir ? "${params.store_dir}/error_correct/${meta.id}" : null

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}.ecco.fq.gz"), emit: reads

    script:
    """
    bbmerge.sh \\
        in="${reads}" \\
        out="${meta.id}.ecco.fq.gz" \\
        ecco mix vstrict ordered \\
        ihist=ihist_merge1.txt \\
        ow=t \\
        t=${task.cpus}

    if [ ! -s "${meta.id}.ecco.fq.gz" ]; then
        echo "[ERROR] Error correction (ecco) produced empty output for ${meta.id}" >&2
        exit 1
    fi
    """
}

process ERROR_CORRECT_ECC {
    tag "${meta.id}"
    label 'process_high'
    conda "${projectDir}/conda-envs/dana-metta-bbmap"
    publishDir "${params.outdir}/error_correct/${meta.id}", mode: 'copy', pattern: '*.fq.gz'
    storeDir params.store_dir ? "${params.store_dir}/error_correct/${meta.id}" : null

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}.eccc.fq.gz"), emit: reads

    script:
    """
    clumpify.sh \\
        in="${reads}" \\
        out="${meta.id}.eccc.fq.gz" \\
        ecc passes=4 reorder \\
        ow=t \\
        t=${task.cpus}

    if [ ! -s "${meta.id}.eccc.fq.gz" ]; then
        echo "[ERROR] Error correction (ecc) produced empty output for ${meta.id}" >&2
        exit 1
    fi
    """
}

process ERROR_CORRECT_TADPOLE {
    tag "${meta.id}"
    label 'process_high'
    conda "${projectDir}/conda-envs/dana-metta-bbmap"
    publishDir "${params.outdir}/error_correct/${meta.id}", mode: 'copy', pattern: '*.fq.gz'
    storeDir params.store_dir ? "${params.store_dir}/error_correct/${meta.id}" : null

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}.ecct.fq.gz"), emit: reads

    script:
    """
    tadpole.sh \\
        in="${reads}" \\
        out="${meta.id}.ecct.fq.gz" \\
        ecc k=62 ordered \\
        prefilter=2 prepasses=auto \\
        ow=t \\
        t=${task.cpus}

    if [ ! -s "${meta.id}.ecct.fq.gz" ]; then
        echo "[ERROR] Error correction (tadpole ecc) produced empty output for ${meta.id}" >&2
        exit 1
    fi
    """
}
