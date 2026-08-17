// Primer removal with cutadapt, and a record of what it actually matched.
//
// PRIMER_ASSIGNMENT turns cutadapt's own per-adapter counts into one row per
// sample, so which assay a sample really carries is reported rather than
// recomputed downstream from taxonomy. It states only what was observed — the
// primer name and read counts — because a primer FASTA carries nothing else.
//
// Three stages:
//   1. DETECT_PRIMERS: samples reads once and matches their 5' ends against the
//      curated table in bin/primer_db.py, falling back to a de-novo consensus,
//      and reports what it found for that one sample
//   2. RESOLVE_PRIMER_SET: reduces those per-sample reports to the one set the
//      whole run is trimmed with
//   3. REMOVE_PRIMERS: one cutadapt pass, -g for R1 and -G for R2 from separate
//      files, so an end whose primer was stripped before submission is left
//      alone while the other end is still trimmed, with --discard-untrimmed
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
    // Evidence only. What a sample is trimmed with comes from RESOLVE_PRIMER_SET,
    // so a sample that yields no consensus reports nothing and still goes on to
    // REMOVE_PRIMERS with the set its assay-mates resolved.
    path("${meta.id}_detected.json"), emit: report, optional: true

    script:
    """
    detect_primers.py "${r1}" "${r2}" \
        -o detected_primers.fa \
        --json "${meta.id}_detected.json" \
        --sample "${meta.id}" \
        -n ${params.primer_detect_reads}
    """
}


// Reduce the per-sample detections to the one set the whole run is trimmed with.
// Handing cutadapt one adapter per sample lets it pick a different winner in each
// one, and plateFor() keys the AUTO_TRIM and LEARN_ERRORS groups on which it
// picked — so a single assay would be split into as many error models as there
// are samples. The set still holds several primers where a run really carries
// several assays; cutadapt chooses per read.
process RESOLVE_PRIMER_SET {
    label 'process_low'
    conda "${projectDir}/envs/python.yml"
    publishDir "${params.outdir}/trimmed", mode: 'copy',
        pattern: "{primer_set.json,assay_positions.tsv}"

    input:
    path(detections)

    output:
    path("primer_set_fwd.fa"),   emit: primers_fwd
    path("primer_set_rev.fa"),   emit: primers_rev
    path("primer_set.json"),     emit: report
    path("assay_positions.tsv"), emit: positions

    script:
    """
    resolve_primer_set.py ${detections} \
        --fwd-out primer_set_fwd.fa \
        --rev-out primer_set_rev.fa \
        --positions assay_positions.tsv \
        --json primer_set.json
    """
}


// Split a bundled pair file into its two ends, so every path into REMOVE_PRIMERS
// hands it one file per end. primers-<FWD>-<REV>.fa holds the forward record
// first and the reverse second, which is the only thing that says which is which
// — the headers carry primer names, not ends.
process SPLIT_PRIMER_PAIR {
    tag "${pair}"
    label 'process_low'
    conda "${projectDir}/envs/python.yml"

    input:
    tuple val(pair), path(primer_file)

    output:
    tuple val(pair), path("${pair}_fwd.fa"), path("${pair}_rev.fa"), emit: ends

    script:
    """
    split_primer_pair.py "${primer_file}" "${pair}_fwd.fa" "${pair}_rev.fa"
    """
}


process REMOVE_PRIMERS {
    tag "${meta.id}"
    cpus 1
    conda "${projectDir}/envs/python.yml"
    publishDir "${params.outdir}/trimmed", mode: 'copy', pattern: "*_cutadapt.log", enabled: !params.store_dir
    storeDir params.store_dir ? "${params.store_dir}/trimmed" : null

    input:
    tuple val(meta), path(r1), path(r2), path(primer_fwd), path(primer_rev)

    output:
    // ASSAY is the adapter cutadapt matched most often — the sample's amplicon,
    // as observed rather than declared. It rides along with the reads so error
    // models and truncation can be grouped per assay without a second pass over
    // the logs (issue #7).
    tuple val(meta), path("${meta.id}_R1.trimmed.fastq.gz"), path("${meta.id}_R2.trimmed.fastq.gz"), env('ASSAY'), emit: reads
    path("${meta.id}_cutadapt.log"), emit: log

    script:
    """
    # One adapter option per end, and none for an end with nothing to remove:
    # cutadapt given an empty adapter file matches nothing, and with
    # --discard-untrimmed that empties the run.
    # if-blocks, not `[ -s f ] && ...`: the script runs under `bash -ue`, where a
    # test that fails is a command that failed and kills the task.
    ARGS=""
    if [ -s "${primer_fwd}" ]; then ARGS="\$ARGS -g file:${primer_fwd}"; fi
    if [ -s "${primer_rev}" ]; then ARGS="\$ARGS -G file:${primer_rev}"; fi

    # --discard-untrimmed drops a pair when the filter applies to it. With both
    # ends trimmed that reads as "either mate carried no primer". With one end
    # trimmed the other mate never carries one, so requiring both is the only way
    # any pair survives.
    PAIR_FILTER=both
    if [ -s "${primer_fwd}" ] && [ -s "${primer_rev}" ]; then PAIR_FILTER=any; fi

    cutadapt \$ARGS \\
        --discard-untrimmed \\
        --pair-filter=\$PAIR_FILTER \\
        -j ${task.cpus} \\
        -e ${params.primer_error_rate} \\
        -o ${meta.id}_R1.trimmed.fastq.gz \\
        -p ${meta.id}_R2.trimmed.fastq.gz \\
        ${r1} ${r2} \\
        > ${meta.id}_cutadapt.log 2>&1

    # Winning 5' adapter, straight out of cutadapt's own per-adapter counts.
    # "none" when nothing matched — a sample with no assay must not silently
    # join whichever group happens to be first.
    ASSAY=\$(awk '/^=== First read: Adapter/{n=\$5}
                  /Trimmed: [0-9]+ times/{
                      if (match(\$0, /Trimmed: [0-9]+/)) {
                          v = substr(\$0, RSTART+9, RLENGTH-9) + 0
                          if (v > best) { best = v; bn = n }
                      }
                  }
                  END { print (bn == "" ? "none" : bn) }' "${meta.id}_cutadapt.log")
    export ASSAY

    echo "[INFO] ${meta.id}: primer removal complete (\$ARGS), assay=\$ASSAY" >&2
    """
}

process PRIMER_ASSIGNMENT {
    tag "primer_assignment"
    label 'process_low'
    conda "${projectDir}/envs/python.yml"
    publishDir "${params.outdir}/trimmed", mode: 'copy'

    input:
    path(cutadapt_logs)
    path(positions)

    output:
    path("primer_assignment.tsv"), emit: tsv

    script:
    // A pre-trimmed run produces no cutadapt log, so the coordinates are the
    // whole record. Both inputs are optional in practice and either may be empty.
    def pos_arg = positions.name != 'NO_POSITIONS' ? "--positions ${positions}" : ''
    """
    parse_primer_assignment.py primer_assignment.tsv ${pos_arg} ${cutadapt_logs}
    """
}
