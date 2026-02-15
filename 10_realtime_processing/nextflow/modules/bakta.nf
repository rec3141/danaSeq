// Bakta CDS-only annotation (fast path — minutes per sample)
// Skips ncRNA/tRNA/CRISPR/sORF scanning to avoid cmscan bottleneck.
// Produces .faa + .gff3 for HMM search and DB integration.

process BAKTA_CDS {
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
        --skip-trna --skip-tmrna --skip-rrna \
        --skip-ncrna --skip-ncrna-region \
        --skip-crispr --skip-sorf \
        --skip-gap --skip-ori --skip-plot \
        "${fasta}"
    """
}

// Bakta full annotation (slow path — complete with ncRNA, tRNA, CRISPR, etc.)
// Runs in parallel with downstream tools; does not block the pipeline.

process BAKTA_FULL {
    tag "${meta.id}"
    label 'process_medium'
    conda "${projectDir}/conda-envs/dana-bakta"
    publishDir "${params.outdir}/${meta.flowcell}/${meta.barcode}/bakta_full", mode: 'copy'

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
