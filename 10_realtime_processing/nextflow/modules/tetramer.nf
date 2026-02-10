// Tetranucleotide frequency analysis for ESOM clustering
// tetramer_freqs_esom.pl is not available via conda; downloaded from GitHub at build time

process TETRAMER_FREQ {
    tag "${meta.id}"
    label 'process_low'
    conda "${projectDir}/conda-envs/dana-tools"
    publishDir "${params.outdir}/${meta.flowcell}/${meta.barcode}/tetra", mode: 'copy'

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${meta.id}.lrn"), emit: lrn

    script:
    """
    # Resolve tetramer_freqs_esom.pl: use from PATH, or download as fallback
    TETRA_SCRIPT="\$(command -v tetramer_freqs_esom.pl 2>/dev/null || true)"
    if [ -z "\$TETRA_SCRIPT" ]; then
        curl -fsSL https://raw.githubusercontent.com/tetramerFreqs/Binning/master/tetramer_freqs_esom.pl \
            -o tetramer_freqs_esom.pl
        chmod +x tetramer_freqs_esom.pl
        TETRA_SCRIPT="./tetramer_freqs_esom.pl"
    fi

    # Create annotation file from FASTA headers
    grep '>' "${fasta}" | sed 's/>//' | paste - - - > annotation.txt

    # Calculate tetranucleotide frequencies
    perl "\$TETRA_SCRIPT" \
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
