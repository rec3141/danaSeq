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
    conda "${projectDir}/conda-envs/dana-mag-assembly"
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

    MIN_SIZE=${params.min_barcode_size ?: 10485760}  # default 10 MB

    # Filter input files: skip tiny and corrupt gzips before concatenation
    GOOD_FILES=""
    SKIPPED=0
    for f in ${fastqs}; do
        fsize=\$(stat -c%s "\$f" 2>/dev/null || echo 0)
        if [ "\$fsize" -lt "\$MIN_SIZE" ]; then
            SKIPPED=\$((SKIPPED + 1))
            continue
        fi
        if ! gzip -t "\$f" 2>/dev/null; then
            echo "[WARNING] ${meta.id}: corrupt gzip \$f (\${fsize} bytes), skipping" >&2
            SKIPPED=\$((SKIPPED + 1))
            continue
        fi
        GOOD_FILES="\$GOOD_FILES \$f"
    done

    if [ -z "\$GOOD_FILES" ]; then
        echo "[WARNING] ${meta.id}: no valid input files (\$SKIPPED skipped)" >&2
        exit 0
    fi
    [ "\$SKIPPED" -gt 0 ] && echo "[INFO] ${meta.id}: skipped \$SKIPPED files (< \${MIN_SIZE} bytes or corrupt)" >&2

    # Concatenate valid files and verify output
    cat \$GOOD_FILES > ${meta.id}_raw.fastq.gz
    if ! gzip -t ${meta.id}_raw.fastq.gz 2>/dev/null; then
        echo "[WARNING] ${meta.id}: concat gzip test failed, recompressing as single stream" >&2
        zcat \$GOOD_FILES | pigz -p ${task.cpus} > ${meta.id}_recomp.fastq.gz
        mv ${meta.id}_recomp.fastq.gz ${meta.id}_raw.fastq.gz
    fi

    # Filter out small barcodes after concat
    filesize=\$(stat -c%s ${meta.id}_raw.fastq.gz)
    if [ "\$filesize" -lt "\$MIN_SIZE" ]; then
        echo "[WARNING] ${meta.id}: concat size \${filesize} bytes < \$MIN_SIZE, skipping" >&2
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
