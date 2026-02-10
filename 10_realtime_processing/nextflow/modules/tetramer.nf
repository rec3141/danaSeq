// Tetranucleotide frequency analysis for ESOM clustering

process TETRAMER_FREQ {
    tag "${meta.id}"
    label 'process_low'
    publishDir "${params.outdir}/${meta.flowcell}/${meta.barcode}/tetra", mode: 'copy'

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${meta.id}.lrn"), emit: lrn

    script:
    """
    # Create annotation file from FASTA headers
    grep '>' "${fasta}" | sed 's/>//' | paste - - - > annotation.txt

    # Calculate tetranucleotide frequencies
    perl ${params.apps}/tetramer_freqs_esom.pl \
        -f "${fasta}" \
        -a annotation.txt \
        -min ${params.min_readlen} \
        -max 10000000

    # Combine names and frequency data into final output
    if ls Tetra_*.names >/dev/null 2>&1 && ls Tetra_*.lrn >/dev/null 2>&1; then
        paste <(awk '\$1!~/^%/' Tetra_*.names) <(awk '\$1!~/^%/' Tetra_*.lrn) \
            | cut -f3,5- > "${meta.id}.lrn"
    else
        echo "[ERROR] Tetramer output files not generated for ${meta.id}" >&2
        exit 1
    fi
    """
}
