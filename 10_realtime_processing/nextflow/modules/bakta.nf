// Bakta gene annotation in metagenome mode (alternative to Prokka)

process BAKTA_ANNOTATE {
    tag "${meta.id}"
    label 'process_medium'
    conda "${projectDir}/conda-envs/dana-bakta"
    publishDir "${params.outdir}/${meta.flowcell}/${meta.barcode}/bakta", mode: 'copy'

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${meta.id}/*.faa"),   emit: proteins
    tuple val(meta), path("${meta.id}/*.tsv"),   emit: tsv
    tuple val(meta), path("${meta.id}/*.gff3"),  emit: gff
    tuple val(meta), path("${meta.id}"),         emit: bakta_dir

    script:
    """
    bakta \
        --db "${params.bakta_db}" \
        --meta \
        --threads 1 \
        --output "${meta.id}" \
        --prefix "${meta.id}" \
        --force \
        "${fasta}"
    """
}
