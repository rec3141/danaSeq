// Bakta CDS-only annotation (fast path — minutes per sample)
// Skips ncRNA/tRNA/CRISPR/sORF scanning to avoid cmscan bottleneck.
// Produces .faa + .gff3 for HMM search and DB integration.
//
// Accepts one or many FASTA files per invocation. In batch mode, main.nf
// groups all files for a sample into a single call so DB loading overhead
// is paid once per barcode. In watch mode, files are passed individually.

process BAKTA_CDS {
    tag "${meta.id}"
    label 'process_medium'
    conda "${projectDir}/conda-envs/dana-bakta"
    publishDir "${params.outdir}/${meta.flowcell}/${meta.barcode}/bakta", mode: 'copy', enabled: !params.store_dir
    storeDir params.store_dir ? "${params.store_dir}/${meta.flowcell}/${meta.barcode}/bakta" : null

    input:
    tuple val(meta), path("fastas/*")

    output:
    tuple val(meta), path("${meta.id}/*.faa"),   emit: proteins
    tuple val(meta), path("${meta.id}/*.tsv"),   emit: tsv
    tuple val(meta), path("${meta.id}/*.gff3"),  emit: gff
    tuple val(meta), path("${meta.id}"),         emit: bakta_dir

    script:
    """
    export PATH="${projectDir}/bin:\$PATH"
    export TMPDIR="\$PWD/tmp"
    mkdir -p "\$TMPDIR"

    cat fastas/* > combined.fa

    bakta \
        --db "${params.bakta_db}" \
        --meta \
        --tmp-dir "\$TMPDIR" \
        --threads ${task.cpus} \
        --output "${meta.id}" \
        --prefix "${meta.id}" \
        --force \
        --skip-trna --skip-tmrna --skip-rrna \
        --skip-ncrna --skip-ncrna-region \
        --skip-crispr --skip-sorf \
        --skip-gap --skip-ori --skip-plot \
        combined.fa
    """
}

// Bakta full annotation (slow path — complete with ncRNA, tRNA, CRISPR, etc.)
// Runs in parallel with downstream tools; does not block the pipeline.

process BAKTA_FULL {
    tag "${meta.id}"
    label 'process_medium'
    conda "${projectDir}/conda-envs/dana-bakta"
    publishDir "${params.outdir}/${meta.flowcell}/${meta.barcode}/bakta_full", mode: 'copy', enabled: !params.store_dir
    storeDir params.store_dir ? "${params.store_dir}/${meta.flowcell}/${meta.barcode}/bakta_full" : null

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${meta.id}/*.faa"),   emit: proteins
    tuple val(meta), path("${meta.id}/*.tsv"),   emit: tsv
    tuple val(meta), path("${meta.id}/*.gff3"),  emit: gff
    tuple val(meta), path("${meta.id}"),         emit: bakta_dir

    script:
    """
    export PATH="${projectDir}/bin:\$PATH"

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
