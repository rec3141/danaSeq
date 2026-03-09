// Per-barcode read concatenation with optional deduplication.
//
// When --input points to a nanopore run directory (fastq_pass/barcode*/),
// main.nf groups FASTQ files by flowcell+barcode and feeds each group here.
// CONCAT_READS concatenates the chunks into one file per barcode, optionally
// deduplicating by read UUID (--dedupe) to remove basecall duplicates.
// Tiny barcodes (< 1 KB) are skipped; downstream filters then drop empty outputs.

process CONCAT_READS {
    tag "${meta.id}"
    label 'process_medium'
    maxForks 32
    time { params.dedupe ? 4.h : 1.h }
    conda "${projectDir}/conda-envs/dana-bbmap"
    publishDir "${params.outdir}/concat", mode: 'copy', enabled: !params.store_dir
    storeDir params.store_dir ? "${params.store_dir}/concat" : null

    input:
    tuple val(meta), path(fastqs)

    output:
    tuple val(meta), path("${meta.id}.fastq.gz"), emit: reads

    script:
    """
    # Always produce output so downstream .collect() never deadlocks
    touch ${meta.id}.fastq.gz

    # Concatenate all FASTQ files for this barcode
    cat ${fastqs} > ${meta.id}_raw.fastq.gz

    # Filter out tiny barcodes (< 1 KB after concat)
    filesize=\$(stat -c%s ${meta.id}_raw.fastq.gz)
    if [ "\$filesize" -lt 1024 ]; then
        echo "[WARNING] ${meta.id}: concat size \${filesize} bytes < 1 KB, skipping" >&2
        exit 0
    fi

    if [ "${params.dedupe}" = "true" ]; then
        # Check if there are actually duplicate read IDs before doing the
        # expensive decompress-dedup-recompress cycle. Single-pass awk: count
        # total headers and unique IDs; exits early on first duplicate found.
        has_dupes=\$(zcat ${meta.id}_raw.fastq.gz \\
            | awk 'NR%4==1 {sub(/^@/,""); sub(/ .*/,""); if (seen[\$0]++) {print "yes"; exit}}')

        if [ "\$has_dupes" = "yes" ]; then
            echo "[INFO] ${meta.id}: duplicate read IDs found, deduplicating" >&2
            zcat ${meta.id}_raw.fastq.gz \\
                | paste - - - - \\
                | awk -F'\\t' '{id=\$1; sub(/^@/,"",id); sub(/ .*/,"",id); if (!seen[id]++) print}' \\
                | tr '\\t' '\\n' \\
                | pigz -p ${task.cpus} > ${meta.id}.fastq.gz
            rm -f ${meta.id}_raw.fastq.gz
        else
            echo "[INFO] ${meta.id}: no duplicate read IDs, skipping dedup" >&2
            mv ${meta.id}_raw.fastq.gz ${meta.id}.fastq.gz
        fi
    else
        mv ${meta.id}_raw.fastq.gz ${meta.id}.fastq.gz
    fi
    """
}
