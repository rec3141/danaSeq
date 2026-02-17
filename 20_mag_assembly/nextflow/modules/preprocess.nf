// Per-barcode read concatenation with optional deduplication

process CONCAT_READS {
    tag "${meta.id}"
    label 'process_low'
    conda "${projectDir}/conda-envs/dana-bbmap"
    publishDir "${params.outdir}/concat", mode: 'copy'
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
        # Create empty output so Nextflow doesn't error, but downstream
        # processes should filter on file size
        touch ${meta.id}.fastq.gz
        exit 0
    fi

    if [ "${params.dedupe}" = "true" ]; then
        dedupe.sh in=${meta.id}_raw.fastq.gz out=${meta.id}.fastq.gz
        rm -f ${meta.id}_raw.fastq.gz
    else
        mv ${meta.id}_raw.fastq.gz ${meta.id}.fastq.gz
    fi
    """
}
