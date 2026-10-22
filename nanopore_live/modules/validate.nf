// FASTQ validation and repair
// Checks gzip integrity; a damaged file is repaired by recovering every intact
// record (bin/fastq_gz_salvage.py), member by member.
//
// The repair used to be BBTools reformat.sh, which deadlocks instead of failing
// on some corrupt gzip streams (its multithreaded reader parks at 0% CPU; hit
// twice on a damaged archive file, 2026-10-20/21, stalling the run for a day)
// and in any case stops at the first bad byte. The salvage resynchronises on the
// next gzip member, so only the damaged members are lost, and it runs under a
// 2h timeout so a pathological file fails the task rather than hanging the run.

process VALIDATE_FASTQ {
    tag "${meta.id}"
    label 'process_low'
    conda "${projectDir}/conda-envs/dana-tools"
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
        set -o pipefail
        if ! timeout 2h fastq_gz_salvage.py "${fastq}" 2> salvage.log \
                | gzip -c > "${meta.id}.validated.fastq.gz" \
           || ! gzip -t "${meta.id}.validated.fastq.gz" 2>/dev/null; then
            echo "[ERROR] Cannot repair corrupted FASTQ: ${fastq} (\$(cat salvage.log))" >&2
            exit 1
        fi
        echo "[INFO] Repaired corrupted FASTQ: ${fastq}: \$(cat salvage.log)" >&2
    fi
    """
}
