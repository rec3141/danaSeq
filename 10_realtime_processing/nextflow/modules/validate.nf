// FASTQ validation and repair
// Checks gzip integrity; attempts repair with BBMap reformat.sh if corrupted

process VALIDATE_FASTQ {
    tag "${meta.id}"
    label 'process_low'
    conda 'bioconda::bbmap'
    errorStrategy 'ignore'

    input:
    tuple val(meta), path(fastq)

    output:
    tuple val(meta), path("${meta.id}.validated.fastq.gz"), emit: validated

    script:
    """
    if gzip -t "${fastq}" 2>/dev/null; then
        ln -s "${fastq}" "${meta.id}.validated.fastq.gz"
    else
        reformat.sh in="${fastq}" out="${meta.id}.validated.fastq.gz" ow 2>&1
        if ! gzip -t "${meta.id}.validated.fastq.gz" 2>/dev/null; then
            echo "[ERROR] Cannot repair corrupted FASTQ: ${fastq}" >&2
            exit 1
        fi
        echo "[INFO] Repaired corrupted FASTQ: ${fastq}" >&2
    fi
    """
}
