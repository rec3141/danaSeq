// Tetranucleotide frequency analysis

process TETRAMER_FREQ {
    tag "${meta.id}"
    label 'process_low'
    publishDir "${params.outdir}/${meta.flowcell}/${meta.barcode}/tetra", mode: 'copy', enabled: !params.store_dir
    storeDir params.store_dir ? "${params.store_dir}/${meta.flowcell}/${meta.barcode}/tetra" : null

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${meta.id}.lrn"), emit: lrn
    tuple val(meta), path("${meta.id}.lengths"), emit: lengths

    script:
    """
    tetramer_freqs \
        -f "${fasta}" \
        -min ${params.min_readlen} \
        -max 10000000 \
        -o "${meta.id}.lrn" \
        --lengths "${meta.id}.lengths"
    """
}
