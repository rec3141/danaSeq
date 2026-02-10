// Prokka gene annotation in metagenome mode

process PROKKA_ANNOTATE {
    tag "${meta.id}"
    label 'process_medium'
    publishDir "${params.outdir}/${meta.flowcell}/${meta.barcode}/prokka", mode: 'copy'

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${meta.id}/*.faa"),  emit: proteins
    tuple val(meta), path("${meta.id}/*.tsv"),  emit: tsv
    tuple val(meta), path("${meta.id}/*.gff"),  emit: gff
    tuple val(meta), path("${meta.id}"),        emit: prokka_dir

    script:
    """
    ${params.prokka_bin} \
        --metagenome \
        --fast \
        --cpus 1 \
        --evalue 1e-20 \
        --outdir "${meta.id}" \
        --force \
        --quiet \
        "${fasta}"

    # Clean up unnecessary Prokka outputs to save space
    rm -f "${meta.id}"/*.{err,fna,fsa,gbk,log,sqn,txt} 2>/dev/null || true
    """
}
