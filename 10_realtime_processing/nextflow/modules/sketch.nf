// Sendsketch taxonomic profiling via k-mer sketching against NCBI nt

process SENDSKETCH {
    tag "${meta.id}"
    label 'process_medium'
    conda 'bioconda::bbmap'
    publishDir "${params.outdir}/${meta.flowcell}/${meta.barcode}/sketch", mode: 'copy'

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${meta.id}.txt"), emit: sketch

    script:
    """
    sendsketch.sh \
        in="${fasta}" \
        address=nt \
        out="${meta.id}.txt" \
        format=3
    """
}
