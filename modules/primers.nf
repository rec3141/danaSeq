// Primer removal with cutadapt, and a record of what it actually matched.
//
// DETECT_PRIMERS samples reads and matches their 5' ends against the curated
// table in bin/primer_db.py. It used to run a full cutadapt pass per candidate
// primer file per sample — 17 passes over every sample, to /dev/null, just to
// count survivors — and reported only the winning filename.
//
// PRIMER_ASSIGNMENT turns cutadapt's own per-adapter counts into one row per
// sample, so which assay a sample really carries is reported rather than
// recomputed downstream from taxonomy. It states only what was observed — the
// primer name and read counts — because a primer FASTA carries nothing else.
//
// Two-pass approach:
//   1. DETECT_PRIMERS: runs all primer pairs on each sample, picks the best match
//   2. REMOVE_PRIMERS: runs the selected primer pair with --discard-untrimmed
//
// If metadata provides a primer_pair column, DETECT_PRIMERS is skipped.

process DETECT_PRIMERS {
    tag "${meta.id}"
    label 'process_low'
    conda "${projectDir}/envs/python.yml"
    publishDir "${params.outdir}/trimmed/detected", mode: 'copy', pattern: "*.json"

    input:
    tuple val(meta), path(r1), path(r2)

    output:
    tuple val(meta), path(r1), path(r2), path("detected_primers.fa"), emit: detected
    path("${meta.id}_detected.json"), emit: report, optional: true

    script:
    """
    detect_primers.py "${r1}" "${r2}" \
        -o detected_primers.fa \
        --json "${meta.id}_detected.json" \
        -n ${params.primer_detect_reads}
    """
}


process REMOVE_PRIMERS {
    tag "${meta.id}"
    cpus 1
    conda "${projectDir}/envs/python.yml"
    publishDir "${params.outdir}/trimmed", mode: 'copy', pattern: "*_cutadapt.log", enabled: !params.store_dir
    storeDir params.store_dir ? "${params.store_dir}/trimmed" : null

    input:
    tuple val(meta), path(r1), path(r2), path(primer_file)

    output:
    tuple val(meta), path("${meta.id}_R1.trimmed.fastq.gz"), path("${meta.id}_R2.trimmed.fastq.gz"), emit: reads
    path("${meta.id}_cutadapt.log"), emit: log

    script:
    """
    cutadapt \\
        -g file:${primer_file} \\
        -G file:${primer_file} \\
        --discard-untrimmed \\
        --pair-filter=any \\
        -j ${task.cpus} \\
        -e ${params.primer_error_rate} \\
        -o ${meta.id}_R1.trimmed.fastq.gz \\
        -p ${meta.id}_R2.trimmed.fastq.gz \\
        ${r1} ${r2} \\
        > ${meta.id}_cutadapt.log 2>&1

    echo "[INFO] ${meta.id}: primer removal complete (${primer_file})" >&2
    """
}

process PRIMER_ASSIGNMENT {
    tag "primer_assignment"
    label 'process_low'
    conda "${projectDir}/envs/python.yml"
    publishDir "${params.outdir}/trimmed", mode: 'copy'

    input:
    path(cutadapt_logs)

    output:
    path("primer_assignment.tsv"), emit: tsv

    script:
    """
    parse_primer_assignment.py primer_assignment.tsv ${cutadapt_logs}
    """
}
