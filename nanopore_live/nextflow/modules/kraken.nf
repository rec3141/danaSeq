// Kraken2 taxonomic classification
// maxForks=1 ensures only one instance runs at a time (database is loaded into RAM)
//
// Accepts one or many FASTA files per invocation. In batch mode, main.nf
// groups all files for a sample (flowcell+barcode) into a single call so the
// DB loads once per sample instead of once per file. In watch mode, files are
// passed individually for live streaming results.

process KRAKEN2_CLASSIFY {
    tag "${meta.id}"
    label 'process_kraken'
    conda "${projectDir}/conda-envs/dana-tools"
    maxForks 1
    publishDir "${params.outdir}/${meta.flowcell}/${meta.barcode}/kraken", mode: 'copy'

    input:
    tuple val(meta), path("fastas/*")

    output:
    tuple val(meta), path("${meta.id}.tsv"),    emit: parsed,  optional: true
    tuple val(meta), path("${meta.id}.report"), emit: report,  optional: true

    script:
    """
    kraken2 \
        --db ${params.kraken_db} \
        --use-names \
        --threads ${task.cpus} \
        --report "${meta.id}.report" \
        fastas/* \
    | gawk -f ${params.kraken_parse_awk} > "${meta.id}.tsv"
    """
}
