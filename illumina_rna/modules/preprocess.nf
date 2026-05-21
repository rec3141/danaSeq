// Preprocessing: optical deduplication, tile filtering, adapter trimming, artifact removal,
// human decontamination, and QC reporting.
// Mirrors illumina_assembly's preprocess so RNA libraries get the same QA/QC treatment.

process CLUMPIFY {
    tag "${meta.id}"
    label 'process_high'
    conda "${projectDir}/conda-envs/dana-illumina-rna-bbmap"
    publishDir "${params.outdir}/preprocess/${meta.id}", mode: 'symlink', pattern: '*.fq.gz'
    storeDir params.store_dir ? "${params.store_dir}/preprocess/${meta.id}" : null

    input:
    tuple val(meta), path(r1), path(r2)

    output:
    tuple val(meta), path("${meta.id}.clumped.fq.gz"), emit: reads

    script:
    def xmx = task.memory ? "-Xmx${(task.memory.toGiga() * 0.85).intValue()}g" : ""
    """
    # Try optical dedup if --run_optical_dedup; only useful for Illumina libraries
    # with intact tile headers, and only feasible when memory >= ~600B/read in the
    # library. Otherwise (BGI/MGI, large libs on small boxes) skip straight to
    # interleave-only via reformat.sh.
    if [ "${params.run_optical_dedup}" = "true" ]; then
        input_size=\$(stat -c%s "${r1}" 2>/dev/null || echo 0)
        timeout_sec=\$(( input_size / 1073741824 * 300 + 120 ))
        set +e
        timeout \${timeout_sec} clumpify.sh ${xmx} \\
            in1="${r1}" in2="${r2}" \\
            out="${meta.id}.clumped.fq.gz" \\
            dedupe optical ow=t t=${task.cpus}
        clump_exit=\$?
        set -e
    else
        clump_exit=1   # force fallback
    fi

    if [ "\${clump_exit:-1}" -ne 0 ] || [ ! -s "${meta.id}.clumped.fq.gz" ] || [ \$(stat -c%s "${meta.id}.clumped.fq.gz" 2>/dev/null || echo 0) -lt 100 ]; then
        [ "${params.run_optical_dedup}" = "true" ] && \\
            echo "[INFO] Optical dedup skipped/failed (exit \${clump_exit:-skip}) — interleaving with reformat.sh" >&2
        rm -f "${meta.id}.clumped.fq.gz"
        reformat.sh ${xmx} in1="${r1}" in2="${r2}" out="${meta.id}.clumped.fq.gz" ow=t
    fi

    if [ ! -s "${meta.id}.clumped.fq.gz" ]; then
        echo "[ERROR] reformat produced empty output for ${meta.id}" >&2
        exit 1
    fi
    """
}

process FILTER_BY_TILE {
    tag "${meta.id}"
    label 'process_high'
    conda "${projectDir}/conda-envs/dana-illumina-rna-bbmap"
    publishDir "${params.outdir}/preprocess/${meta.id}", mode: 'symlink', pattern: '*.fq.gz'
    storeDir params.store_dir ? "${params.store_dir}/preprocess/${meta.id}" : null

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}.filtered_by_tile.fq.gz"), emit: reads

    script:
    def xmx = task.memory ? "-Xmx${(task.memory.toGiga() * 0.85).intValue()}g" : ""
    """
    input_size=\$(stat -c%s "${reads}" 2>/dev/null || echo 0)
    timeout_sec=\$(( input_size / 1073741824 * 300 + 120 ))

    set +e
    timeout \${timeout_sec} filterbytile.sh ${xmx} \\
        in="${reads}" \\
        out="${meta.id}.filtered_by_tile.fq.gz" \\
        ow=t \\
        t=${task.cpus}
    fbt_exit=\$?
    set -e

    if [ \$fbt_exit -ne 0 ] || [ ! -s "${meta.id}.filtered_by_tile.fq.gz" ]; then
        echo "[WARNING] filterbytile failed (exit \$fbt_exit) — passing through" >&2
        cp "${reads}" "${meta.id}.filtered_by_tile.fq.gz"
    fi
    """
}

process BBDUK_TRIM {
    tag "${meta.id}"
    label 'process_high'
    conda "${projectDir}/conda-envs/dana-illumina-rna-bbmap"
    publishDir "${params.outdir}/preprocess/${meta.id}", mode: 'symlink', pattern: '*.fq.gz'
    storeDir params.store_dir ? "${params.store_dir}/preprocess/${meta.id}" : null

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}.trimmed.fq.gz"), emit: reads

    script:
    def xmx = task.memory ? "-Xmx${(task.memory.toGiga() * 0.85).intValue()}g" : ""
    """
    bbduk.sh ${xmx} \\
        in="${reads}" \\
        out="${meta.id}.trimmed.fq.gz" \\
        ktrim=r k=23 mink=11 hdist=1 \\
        tbo tpe \\
        minlen=${params.min_readlen} \\
        ref=adapters \\
        ftm=5 ordered \\
        ow=t \\
        t=${task.cpus}

    if [ ! -s "${meta.id}.trimmed.fq.gz" ]; then
        echo "[ERROR] bbduk trim produced empty output for ${meta.id}" >&2
        exit 1
    fi
    """
}

process BBDUK_FILTER {
    tag "${meta.id}"
    label 'process_high'
    conda "${projectDir}/conda-envs/dana-illumina-rna-bbmap"
    publishDir "${params.outdir}/preprocess/${meta.id}", mode: 'symlink', pattern: '*.fq.gz'
    storeDir params.store_dir ? "${params.store_dir}/preprocess/${meta.id}" : null

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}.filtered.fq.gz"), emit: reads

    script:
    def xmx = task.memory ? "-Xmx${(task.memory.toGiga() * 0.85).intValue()}g" : ""
    """
    bbduk.sh ${xmx} \\
        in="${reads}" \\
        out="${meta.id}.filtered.fq.gz" \\
        k=31 \\
        ref=artifacts,phix \\
        entropy=0.95 \\
        ordered cardinality \\
        ow=t \\
        t=${task.cpus}

    if [ ! -s "${meta.id}.filtered.fq.gz" ]; then
        echo "[ERROR] bbduk filter produced empty output for ${meta.id}" >&2
        exit 1
    fi
    """
}

process REMOVE_HUMAN {
    tag "${meta.id}"
    label 'process_high'
    conda "${projectDir}/conda-envs/dana-illumina-rna-bbmap"
    publishDir "${params.outdir}/preprocess/${meta.id}", mode: 'symlink', pattern: '*.fq.gz'
    storeDir params.store_dir ? "${params.store_dir}/preprocess/${meta.id}" : null

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}.nohuman.fq.gz"), emit: reads
    path("${meta.id}.human_contam.fq.gz"),             emit: contam, optional: true

    script:
    def xmx = task.memory ? "-Xmx${(task.memory.toGiga() * 0.85).intValue()}g" : ""
    """
    removehuman.sh ${xmx} \\
        in="${reads}" \\
        outu="${meta.id}.nohuman.fq.gz" \\
        outm="${meta.id}.human_contam.fq.gz" \\
        path="${params.human_ref}" \\
        t=${task.cpus}

    if [ ! -s "${meta.id}.nohuman.fq.gz" ]; then
        echo "[ERROR] Human removal produced empty output for ${meta.id}" >&2
        exit 1
    fi
    """
}

process FASTQC {
    tag "${meta.id}"
    label 'process_medium'
    conda "${projectDir}/conda-envs/dana-illumina-rna-bbmap"
    publishDir "${params.outdir}/preprocess/${meta.id}/fastqc", mode: 'copy'
    storeDir params.store_dir ? "${params.store_dir}/preprocess/${meta.id}/fastqc" : null

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.html"), path("*.zip"), emit: reports

    script:
    """
    reformat.sh in="${reads}" out1=r1.fq.gz out2=r2.fq.gz ow=t

    fastqc -t ${task.cpus} --noextract -o . r1.fq.gz r2.fq.gz

    for f in *.html *.zip; do
        mv "\$f" "${meta.id}.\$f"
    done
    """
}

process COUNT_READS {
    // Lightweight stage-end read counter for the viz read-flow funnel.
    tag "${meta.id}:${stage}"
    label 'process_low'
    conda "${projectDir}/conda-envs/dana-illumina-rna-bbmap"

    input:
    tuple val(meta), val(stage), path(reads)

    output:
    tuple val(meta), val(stage), path("${meta.id}.${stage}.count.tsv"), emit: count

    script:
    """
    n=\$(reformat.sh in="${reads}" 2>&1 | awk '/Input:/ {print \$2; exit}')
    [ -z "\$n" ] && n=0
    printf '%s\\t%s\\t%s\\n' "${meta.id}" "${stage}" "\$n" > "${meta.id}.${stage}.count.tsv"
    """
}
