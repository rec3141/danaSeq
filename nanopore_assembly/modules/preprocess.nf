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
    conda "${projectDir}/conda-envs/dana-mag-assembly"
    // Per-barcode concatenations duplicate the input reads byte for byte
    // (~300 GB per large co-assembly), so they are neither published nor
    // stored unless --publish_concat is given.
    publishDir "${params.outdir}/concat", mode: 'copy', enabled: params.publish_concat

    input:
    tuple val(meta), path(fastqs)

    output:
    tuple val(meta), path("${meta.id}.fastq.gz"), emit: reads

    script:
    """
    # Always produce output so downstream .collect() never deadlocks
    touch ${meta.id}.fastq.gz

    MIN_SIZE=${params.min_barcode_size ?: 10485760}  # default 10 MB

    # Filter input files: skip corrupt gzips before concatenation
    GOOD_FILES=""
    SKIPPED=0
    for f in ${fastqs}; do
        if ! gzip -t "\$f" 2>/dev/null; then
            echo "[WARNING] ${meta.id}: corrupt gzip \$f (\$(stat -c%s "\$f" 2>/dev/null || echo 0) bytes), skipping" >&2
            SKIPPED=\$((SKIPPED + 1))
            continue
        fi
        GOOD_FILES="\$GOOD_FILES \$f"
    done

    if [ -z "\$GOOD_FILES" ]; then
        echo "[WARNING] ${meta.id}: no valid input files (\$SKIPPED skipped)" >&2
        exit 0
    fi
    [ "\$SKIPPED" -gt 0 ] && echo "[INFO] ${meta.id}: skipped \$SKIPPED corrupt files" >&2

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

// Concatenate all per-barcode reads, deduplicate, and optionally filter by quality/length.
// Pipes cat directly into fastq_filter (no intermediate file on disk).
// Intermediate read sets are large (~700 GB plain at 340 Gbp) and every
// downstream consumer pays a single-threaded gzip inflate per pass, so they
// are written as plain FASTQ in the work dir and never stored: regenerating
// them (1-2 h each at that scale) is cheaper than keeping them, and the work
// dir is meant to live on node-local disk (see run-nanopore-assembly.sh).
// Only the assembly itself goes to --store_dir.
def writeReads(String basename, int cpus) {
    return "cat > ${basename}"
}

// fastq_filter replaces both BBMap dedupe and filtlong in a single streaming pass.
process PREPARE_READS {
    tag "prepare-reads"
    label 'process_high'
    conda "${projectDir}/conda-envs/dana-mag-assembly"
    // No publishDir / storeDir — all_reads.fastq is a large transient
    // intermediate (see writeReads above), not a result.

    input:
    path(fastqs)

    output:
    path("all_reads.fastq*"),  emit: reads
    path("read_map.tsv.gz"),   emit: read_map, optional: true

    script:
    def filter_args = params.dedupe ? "" : "--no_dedupe"
    if (params.filtlong_size) {
        filter_args += " --target_bases ${params.filtlong_size}"
    }
    // read_id -> flowcell_barcode for every accepted read. fastq_filter already
    // parses each header, so this costs no extra pass over the reads. It lets a
    // pooled BAM (e.g. Flye's polishing alignment) be split per sample without
    // re-reading the fastqs.
    if (params.emit_read_map) {
        filter_args += " --read_map read_map.tsv"
    }
    """
    # Stream each input to fastq_filter. FASTA inputs (.fa/.fasta[.gz]) are
    # converted to FASTQ with a placeholder quality (Q40) on the fly — used for
    # pre-QC'd reads exported as fasta (e.g. nanopore_live fa/ store). zcat -f
    # transparently handles both plain and gzipped inputs. awk accumulates
    # multi-line fasta records so wrapped sequences are handled correctly.
    for f in ${fastqs}; do
        case "\$f" in
            *.fa|*.fasta|*.fa.gz|*.fasta.gz)
                zcat -f "\$f" | awk '
                    /^>/ { if (s) { print "@"n; print s; print "+"; q=s; gsub(/./,"I",q); print q } n=substr(\$0,2); s=""; next }
                    { s=s\$0 }
                    END { if (s) { print "@"n; print s; print "+"; q=s; gsub(/./,"I",q); print q } }' | gzip ;;
            *)
                cat "\$f" ;;
        esac
    done | fastq_filter ${filter_args} | ${writeReads('all_reads.fastq', task.cpus)}

    if [ -s read_map.tsv ]; then
        gzip -1 read_map.tsv
        echo "[INFO] read map: \$(zcat read_map.tsv.gz | wc -l) reads"
    fi
    """
}

// Map long reads to human reference with minimap2, keep only unmapped reads.
// Uses samtools -f 4 (unmapped) to extract non-human reads, then converts back to FASTQ.
process REMOVE_HUMAN {
    tag "remove-human"
    label 'process_high'
    conda "${projectDir}/conda-envs/dana-mag-assembly"
    // No storeDir — nohuman_reads.fastq is a large transient intermediate.

    input:
    path(reads)

    output:
    path("nohuman_reads.fastq*"), emit: reads

    script:
    """
    # Read counts are taken from the stream (primary records in = reads in,
    # FASTQ records out = reads out) instead of re-reading both files
    # afterwards, which cost two full decompressions of the read set.
    minimap2 -a -x map-ont --secondary=no -t ${task.cpus} \\
        "${params.human_ref}" "${reads}" \\
        | samtools view -b - \\
        | tee >(samtools view -c -F 0x900 - > input.count) \\
        | samtools view -b -f 4 - \\
        | samtools fastq -@ ${task.cpus} - \\
        | tee >(awk 'NR%4==1' | wc -l > output.count) \\
        | ${writeReads('nohuman_reads.fastq', task.cpus)}

    # process substitutions may still be flushing after the pipeline returns
    for i in \$(seq 1 120); do
        [ -s input.count ] && [ -s output.count ] && break
        sleep 1
    done
    if [ ! -s nohuman_reads.fastq* ]; then
        echo "[ERROR] Human removal produced empty output" >&2
        exit 1
    fi

    input_count=\$(cat input.count)
    output_count=\$(cat output.count)
    removed=\$((input_count - output_count))
    echo "[INFO] Human removal: \${input_count} reads in, \${output_count} out, \${removed} removed (\$(( removed * 100 / (input_count + 1) ))%)" >&2
    """
}
