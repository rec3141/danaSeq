// Kraken2 taxonomic classification
// maxForks=1 ensures only one instance runs at a time (database is loaded into RAM)

process KRAKEN2_CLASSIFY {
    tag "${meta.id}"
    label 'process_kraken'
    conda "${projectDir}/conda-envs/dana-tools"
    maxForks 1
    publishDir "${params.outdir}/${meta.flowcell}/${meta.barcode}/kraken", mode: 'copy'

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${meta.id}.tsv"),    emit: parsed,  optional: true
    tuple val(meta), path("${meta.id}.report"), emit: report,  optional: true

    script:
    """
    kraken2 \
        --db ${params.kraken_db} \
        --use-names \
        --threads 1 \
        --report "${meta.id}.report" \
        "${fasta}" \
    | gawk -f ${params.kraken_parse_awk} > "${meta.id}.tsv"
    """
}
