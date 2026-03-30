// FASTQ validation and repair
// Checks gzip integrity; attempts repair with BBMap reformat.sh if corrupted

process VALIDATE_FASTQ {
    tag "${meta.id}"
    label 'process_low'
    conda "${projectDir}/conda-envs/dana-bbmap"
    errorStrategy 'ignore'

    input:
    tuple val(meta), path(fastq)

    output:
    tuple val(meta), path("${meta.id}.validated.fastq.gz"), emit: validated

    script:
    def mem_mb = task.memory ? (task.memory.toMega() * 85 / 100).intValue() : 1700
    """
    if gzip -t "${fastq}" 2>/dev/null; then
        ln -s "${fastq}" "${meta.id}.validated.fastq.gz"
    else
        reformat.sh -Xmx${mem_mb}m in="${fastq}" out="${meta.id}.validated.fastq.gz" ow 2>&1
        if ! gzip -t "${meta.id}.validated.fastq.gz" 2>/dev/null; then
            echo "[ERROR] Cannot repair corrupted FASTQ: ${fastq}" >&2
            exit 1
        fi
        echo "[INFO] Repaired corrupted FASTQ: ${fastq}" >&2
    fi
    """
}
