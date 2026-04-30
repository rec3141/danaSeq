// Kraken2 taxonomic classification
// maxForks=1 ensures only one instance runs at a time within a pipeline.
// --memory-mapping makes kraken2 mmap() the *.k2d files instead of slurping
// them into its own heap; concurrent kraken invocations (e.g. multiple
// nanopore_live pipelines on the same host) then share the DB via the OS
// page cache, so RAM use is ~16 GB total regardless of how many run side by
// side. Cold-cache: the first task pays the disk read; subsequent tasks
// (this pipeline OR another) start near-instant.
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
    publishDir "${params.outdir}/${meta.flowcell}/${meta.barcode}/kraken", mode: 'copy', enabled: !params.store_dir
    storeDir params.store_dir ? "${params.store_dir}/${meta.flowcell}/${meta.barcode}/kraken" : null

    input:
    tuple val(meta), path("fastas/*")

    output:
    tuple val(meta), path("${meta.id}.tsv"),    emit: parsed,  optional: true
    tuple val(meta), path("${meta.id}.report"), emit: report,  optional: true

    script:
    """
    kraken2 \
        --db ${params.kraken_db} \
        --memory-mapping \
        --use-names \
        --threads ${task.cpus} \
        --report "${meta.id}.report" \
        fastas/* \
    | gawk -f ${params.kraken_parse_awk} > "${meta.id}.tsv"
    """
}
