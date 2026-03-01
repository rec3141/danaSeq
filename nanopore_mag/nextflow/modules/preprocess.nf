// Per-barcode read concatenation with optional deduplication

process CONCAT_READS {
    tag "${meta.id}"
    label 'process_low'
    time  = { params.dedupe ? 4.h : 1.h }
    conda "${projectDir}/conda-envs/dana-bbmap"
    publishDir "${params.outdir}/concat", mode: 'copy', enabled: !params.store_dir
    storeDir params.store_dir ? "${params.store_dir}/concat" : null

    input:
    tuple val(meta), path(fastqs)

    output:
    tuple val(meta), path("${meta.id}.fastq.gz"), emit: reads

    script:
    """
    # Concatenate all FASTQ files for this barcode
    cat ${fastqs} > ${meta.id}_raw.fastq.gz

    # Filter out tiny barcodes (< 1 KB after concat)
    filesize=\$(stat -c%s ${meta.id}_raw.fastq.gz)
    if [ "\$filesize" -lt 1024 ]; then
        echo "[WARNING] ${meta.id}: concat size \${filesize} bytes < 1 KB, skipping" >&2
        touch ${meta.id}.fastq.gz
        exit 0
    fi

    if [ "${params.dedupe}" = "true" ]; then
        # Deduplicate by read UUID only — same read basecalled twice shares a UUID
        # but may have different header metadata (field order, timezone format, runid).
        # BBMap dedupe uses full-header matching so misses these; awk on the first
        # whitespace-delimited field (the UUID) is exact. Uses ~200 MB regardless
        # of file size; pigz parallelises the output compression.
        zcat ${meta.id}_raw.fastq.gz \\
            | paste - - - - \\
            | awk -F'\\t' '{id=\$1; sub(/^@/,"",id); sub(/ .*/,"",id); if (!seen[id]++) print}' \\
            | tr '\\t' '\\n' \\
            | pigz -p ${task.cpus} > ${meta.id}.fastq.gz
        rm -f ${meta.id}_raw.fastq.gz
    else
        mv ${meta.id}_raw.fastq.gz ${meta.id}.fastq.gz
    fi
    """
}
