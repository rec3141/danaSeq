// Assembly: Tadpole, Megahit, SPAdes, metaSPAdes

process ASSEMBLE_TADPOLE {
    tag "${meta.id}"
    label 'process_high'
    conda "${projectDir}/conda-envs/dana-metta-bbmap"
    publishDir "${params.outdir}/assembly/${meta.id}", mode: 'copy', pattern: '*.fasta'
    storeDir params.store_dir ? "${params.store_dir}/assembly/${meta.id}" : null

    input:
    tuple val(meta), path(merged), path(qtrimmed)

    output:
    tuple val(meta), path("${meta.id}.tadpole.fasta"), emit: contigs

    script:
    """
    set +e
    tadpole.sh \\
        in="${merged},${qtrimmed}" \\
        out="${meta.id}.tadpole.fasta" \\
        k=124 \\
        prefilter=2 prepasses=auto \\
        ow=t \\
        t=${task.cpus}
    exit_code=\$?
    set -e

    if [ \$exit_code -ne 0 ] || [ ! -s "${meta.id}.tadpole.fasta" ]; then
        echo "[WARNING] Tadpole assembly failed or produced empty output for ${meta.id}" >&2
        touch "${meta.id}.tadpole.fasta"
    fi
    """
}

process ASSEMBLE_MEGAHIT {
    tag "${meta.id}"
    label 'process_high'
    conda "${projectDir}/conda-envs/dana-metta-megahit"
    publishDir "${params.outdir}/assembly/${meta.id}", mode: 'copy', pattern: '*.fasta'
    storeDir params.store_dir ? "${params.store_dir}/assembly/${meta.id}" : null

    input:
    tuple val(meta), path(merged), path(qtrimmed)

    output:
    tuple val(meta), path("${meta.id}.megahit.fasta"), emit: contigs

    script:
    """
    set +e
    megahit \\
        --k-min 45 --k-max 225 --k-step 26 \\
        --min-count 2 \\
        -r "${merged}" \\
        --12 "${qtrimmed}" \\
        -o megahit_out \\
        -t ${task.cpus}
    exit_code=\$?
    set -e

    if [ \$exit_code -ne 0 ] || [ ! -s megahit_out/final.contigs.fa ]; then
        echo "[WARNING] Megahit assembly failed or produced empty output for ${meta.id}" >&2
        touch "${meta.id}.megahit.fasta"
    else
        cp megahit_out/final.contigs.fa "${meta.id}.megahit.fasta"
    fi

    # Clean up intermediate files to save disk
    rm -rf megahit_out/intermediate_contigs
    """
}

process ASSEMBLE_SPADES {
    tag "${meta.id}"
    label 'process_high'
    conda "${projectDir}/conda-envs/dana-metta-spades"
    publishDir "${params.outdir}/assembly/${meta.id}", mode: 'copy', pattern: '*.fasta'
    storeDir params.store_dir ? "${params.store_dir}/assembly/${meta.id}" : null

    input:
    tuple val(meta), path(merged), path(qtrimmed)

    output:
    tuple val(meta), path("${meta.id}.spades.fasta"), emit: contigs

    script:
    def mem_gb = task.memory ? (task.memory.toGiga()) : 250
    """
    set +e
    spades.py \\
        -k 25,55,95,125 \\
        --phred-offset 33 \\
        -s "${merged}" \\
        --12 "${qtrimmed}" \\
        -o spades_out \\
        --only-assembler \\
        --mem ${mem_gb} \\
        -t ${task.cpus}
    exit_code=\$?
    set -e

    if [ \$exit_code -ne 0 ] || [ ! -s spades_out/contigs.fasta ]; then
        echo "[WARNING] SPAdes assembly failed or produced empty output for ${meta.id}" >&2
        touch "${meta.id}.spades.fasta"
    else
        cp spades_out/contigs.fasta "${meta.id}.spades.fasta"
    fi

    # Clean up intermediate k-mer directories
    rm -rf spades_out/K*
    """
}

process ASSEMBLE_METASPADES {
    tag "${meta.id}"
    label 'process_high'
    conda "${projectDir}/conda-envs/dana-metta-spades"
    publishDir "${params.outdir}/assembly/${meta.id}", mode: 'copy', pattern: '*.fasta'
    storeDir params.store_dir ? "${params.store_dir}/assembly/${meta.id}" : null

    input:
    tuple val(meta), path(normalized)

    output:
    tuple val(meta), path("${meta.id}.metaspades.fasta"), emit: contigs

    script:
    def mem_gb = task.memory ? (task.memory.toGiga()) : 250
    """
    set +e
    spades.py \\
        -k 25,55,77 \\
        --phred-offset 33 \\
        --12 "${normalized}" \\
        -o spadesmeta_out \\
        --only-assembler --meta \\
        --mem ${mem_gb} \\
        -t ${task.cpus}
    exit_code=\$?
    set -e

    if [ \$exit_code -ne 0 ] || [ ! -s spadesmeta_out/contigs.fasta ]; then
        echo "[WARNING] metaSPAdes assembly failed or produced empty output for ${meta.id}" >&2
        touch "${meta.id}.metaspades.fasta"
    else
        cp spadesmeta_out/contigs.fasta "${meta.id}.metaspades.fasta"
    fi

    # Clean up intermediate k-mer directories
    rm -rf spadesmeta_out/K*
    """
}
