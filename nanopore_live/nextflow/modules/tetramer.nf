// Tetranucleotide frequency analysis
// Uses tetramer_freqs.py (Python replacement for tetramer_freqs_esom.pl)

process TETRAMER_FREQ {
    tag "${meta.id}"
    label 'process_low'
    conda "${projectDir}/conda-envs/dana-tools"
    publishDir "${params.outdir}/${meta.flowcell}/${meta.barcode}/tetra", mode: 'copy'

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${meta.id}.lrn"), emit: lrn
    tuple val(meta), path("${meta.id}.lengths"), emit: lengths

    script:
    """
    python3 ${projectDir}/bin/tetramer_freqs.py \
        -f "${fasta}" \
        -min ${params.min_readlen} \
        -max 10000000 \
        -o "${meta.id}.lrn" \
        --lengths "${meta.id}.lengths"
    """
}
