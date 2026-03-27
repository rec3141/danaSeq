// Quality control: BBDuk adapter/quality trimming, fastq_filter length/quality filtering, FASTA conversion

process QC_BBDUK {
    tag "${meta.id}"
    label 'process_medium'
    conda "${projectDir}/conda-envs/dana-bbmap"

    input:
    tuple val(meta), path(fastq)

    output:
    tuple val(meta), path("${meta.id}.bbduk.fastq.gz"), emit: trimmed

    script:
    def mem_mb = task.memory ? (task.memory.toMega() * 85 / 100).intValue() : 3400
    """
    bbduk.sh \
        -Xmx${mem_mb}m \
        in="${fastq}" \
        out="${meta.id}.bbduk.fastq.gz" \
        ref=adapters,artifacts,phix,lambda \
        qtrim=rl trimq=15 entropy=0.75 qin=33 \
        minlength=${params.min_readlen}
    """
}

process QC_FASTQ_FILTER {
    tag "${meta.id}"
    label 'process_low'

    input:
    tuple val(meta), path(fastq)

    output:
    tuple val(meta), path("${meta.id}.filtered.fastq"), emit: filtered

    script:
    """
    fastq_filter \
        --no_dedupe \
        --min_length ${params.min_readlen} \
        --keep_percent ${params.keep_percent} \
        -o "${meta.id}.filtered.fastq" \
        "${fastq}"

    # Fail if output is empty (all reads filtered out)
    if [ ! -s "${meta.id}.filtered.fastq" ]; then
        echo "[WARNING] All reads filtered out for ${meta.id}" >&2
        exit 1
    fi
    """
}

process CONVERT_TO_FASTA {
    tag "${meta.id}"
    label 'process_low'
    publishDir "${params.outdir}/${meta.flowcell}/${meta.barcode}/fa", mode: 'copy', enabled: !params.store_dir
    storeDir params.store_dir ? "${params.store_dir}/${meta.flowcell}/${meta.barcode}/fa" : null

    input:
    tuple val(meta), path(fastq)

    output:
    tuple val(meta), path("${meta.id}.fa"), emit: fasta

    script:
    """
    awk 'NR%4==1{sub(/^@/,">"); sub(/ .*/,""); print} NR%4==2{print}' \
        "${fastq}" > "${meta.id}.fa"

    if [ ! -s "${meta.id}.fa" ]; then
        echo "[WARNING] Empty FASTA output for ${meta.id}" >&2
        exit 1
    fi
    """
}
