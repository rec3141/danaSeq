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
    // Present only when something needs a human's attention; main.nf copies
    // each line into the Nextflow log, where it is seen, rather than leaving it
    // in .command.err, where it is not.
    path("prepare_reads.warnings.txt"), emit: warnings, optional: true

    script:
    def filter_args = params.dedupe ? "" : "--no_dedupe"
    if (params.filtlong_size) {
        filter_args += " --target_bases ${params.filtlong_size}"
        // fastq_filter selects the best --target_bases with a two-pass bucket
        // sort by default, so the kept reads depend only on their scores and
        // not on the order the samples happen to arrive in. Spill into the task
        // directory rather than the container's /tmp: it needs room for one
        // uncompressed copy of the input (~700 GB for a 250 Gbp co-assembly),
        // and the launcher already puts the work dir on node-local SLURM_TMPDIR.
        // --onepass restores the legacy order-dependent streaming threshold,
        // which is cheaper on disk and nothing else.
        filter_args += params.onepass_filter ? " --onepass" : " --spill_dir ."
    }
    // Flye documents --nano-hq for reads under ~5% error. Hold the reads to
    // that unless told not to, and measure the fit either way.
    if (params.read_type == 'nano-hq') {
        filter_args += params.nano_hq_filter ? " --nano_hq" : " --nano_hq_check"
    }
    """
    # Nextflow runs this with `bash -ue`, which does NOT include pipefail, so
    # the exit status of `... | fastq_filter | writeReads` is writeReads'.
    # fastq_filter aborted mid-stream on both grex co-assemblies on 2026-09-22
    # and this task still recorded .exitcode 0: marine reached Flye with 134 Gbp
    # of ~338, freshwater with 65 Gbp of ~400, and the freshwater run went on to
    # complete every downstream stage and report success on a 2.63 Gbp assembly.
    # A filter that dies must fail the task.
    set -o pipefail

    # Read map, if asked for: one file per sample, keyed on the staged filename.
    # Independent of the record stream, so it runs as its own pass.
    if [ "${params.emit_read_map}" = "true" ]; then
        for f in ${fastqs}; do
        # The staged files are CONCAT_READS outputs named <meta.id>.fastq.gz, so
        # the filename IS the sample id -- the same value main.nf derives from the
        # barcode directory. Do NOT parse the sample out of read headers. Three
        # header formats occur in practice and each defeats a different parser:
        #   MinKNOW key=value, but unbarcoded runs (filed under barcode00 by
        #     convention) carry no barcode= at all;
        #   SAM tags (PU:Z/SM:Z) from a Dorado -> BAM -> samtools fastq round
        #     trip, with no flow_cell_id=/barcode= anywhere;
        #   bare UUIDs, stripped of every field.
        # The directory is authoritative; the header is not.
        SAMPLE=\$(basename "\$f")
        SAMPLE=\${SAMPLE%.gz}; SAMPLE=\${SAMPLE%.fastq}; SAMPLE=\${SAMPLE%.fq}
        SAMPLE=\${SAMPLE%.fasta}; SAMPLE=\${SAMPLE%.fa}

        # Read map: one file per sample, written synchronously before the record
        # stream is handed to fastq_filter. Deliberately NOT `tee >(... >> shared)`:
        # bash does not wait for a process substitution before the next iteration,
        # so concurrent appends tear lines at ~4 KB buffer boundaries and invent
        # samples. Costs one extra decompress per file; correctness over speed,
        # and the whole thing is opt-in.
        if [ "${params.emit_read_map}" = "true" ]; then
            mkdir -p read_map.d
            case "\$f" in
                *.fa|*.fasta|*.fa.gz|*.fasta.gz)
                    zcat -f "\$f" | awk -v s="\$SAMPLE" \
                        '/^>/ { print substr(\$1, 2) "\\t" s }' > "read_map.d/\$SAMPLE.tsv" ;;
                *)
                    zcat -f "\$f" | awk -v s="\$SAMPLE" \
                        'NR % 4 == 1 { print substr(\$1, 2) "\\t" s }' > "read_map.d/\$SAMPLE.tsv" ;;
            esac
        fi

        done
    fi

    # Hand fastq_filter the files themselves rather than a pipe. Streaming them
    # through `for f; do cat "\$f"; done | fastq_filter` failed on three of
    # three grex co-assembly attempts on 2026-09-22/23: zlib's read() on the
    # pipe returned ENODATA (zlib=-1 errno=61) partway through, at a different
    # point each time on identical input, while a standalone replay of the same
    # 276 files read all 338.1 Gbp cleanly. A pipe read cannot return ENODATA
    # under POSIX, so whatever is happening is at the pipe layer inside the
    # container, and opening regular files sidesteps it entirely. It also gives
    # fastq_filter real file sizes, so its input estimate works again.
    #
    # FASTA inputs (.fa/.fasta[.gz], e.g. pre-QC'd reads from nanopore_live's
    # fa/ store) still need converting to FASTQ with a placeholder Q40 quality,
    # which only the stream can do, so that path keeps the loop.
    HAS_FASTA=0
    for f in ${fastqs}; do
        case "\$f" in *.fa|*.fasta|*.fa.gz|*.fasta.gz) HAS_FASTA=1 ;; esac
    done

    if [ "\$HAS_FASTA" = "0" ]; then
        fastq_filter ${filter_args} ${fastqs} 2> filter.log | ${writeReads('all_reads.fastq', task.cpus)}
    else
        echo "[INFO] FASTA inputs present: streaming through the conversion loop" >&2
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
        done | fastq_filter ${filter_args} 2> filter.log | ${writeReads('all_reads.fastq', task.cpus)}
    fi

    cat filter.log >&2

    # Warnings that need a human go to prepare_reads.warnings.txt, which main.nf
    # copies into the Nextflow log. Deliberately warnings, not exits -- see the
    # errorStrategy note in nextflow.config.
    warn() { echo "\$*" >> prepare_reads.warnings.txt; echo "[WARNING] \$*" >&2; }

    # 1. Did the whole read set arrive? Compare the bases fastq_filter read
    #    against the gzip it was given, which holds whatever the filtering does:
    #    ONT FASTQ runs ~1 base per gzipped byte (338.1 Gbp from 330.1 GB on
    #    marine) and FASTA runs higher. On 2026-09-22 truncated streams reached
    #    Flye as 134 Gbp of 338 and 65 of 327, and the freshwater run went on
    #    to assemble 2.63 Gbp and report success.
    IN_GZ=\$(stat -Lc%s ${fastqs} 2>/dev/null | awk '{s+=\$1} END{print s+0}')
    IN_BASES=\$(awk '/Input bases:/ {print \$NF}' filter.log)
    if [ -n "\$IN_BASES" ] && awk -v b="\$IN_BASES" -v g="\$IN_GZ" 'BEGIN{exit !(b < 0.6*g)}'; then
        warn "PREPARE_READS read only \$IN_BASES bases from \$IN_GZ bytes of gzip (expect ~1 base per byte). The read stream probably ended early; everything downstream will be built from a fraction of the data and still look successful. Check Flye's 'Total read length'."
    fi

    # 2. With --read_type nano-hq, how much of the data is outside the mode's
    #    envelope? Past the threshold the data are saying nano-raw may fit better.
    FIT=\$(grep 'nano-hq fit:' filter.log || true)
    if [ -n "\$FIT" ]; then
        PCT=\$(echo "\$FIT" | awk -F'[(%]' '{print \$2}')
        if awk -v p="\$PCT" -v t="${params.nano_hq_warn_pct}" 'BEGIN{exit !(p > t)}'; then
            if [ "${params.nano_hq_filter}" = "true" ]; then
                warn "--read_type nano-hq: \${PCT}% of input bases are in reads above 5% mean error, beyond nano-hq's envelope, and were DROPPED before assembly (threshold ${params.nano_hq_warn_pct}%). That data is lost to the assembly; consider --read_type nano-raw, which keeps it."
            else
                warn "--read_type nano-hq with --nano_hq_filter false: \${PCT}% of input bases are in reads above 5% mean error, beyond nano-hq's envelope, and were KEPT (threshold ${params.nano_hq_warn_pct}%). Flye will treat them as under 5% error; consider --read_type nano-raw."
            fi
        fi
        echo "[INFO] \$FIT" >&2
    fi

    if [ -d read_map.d ]; then
        cat read_map.d/*.tsv | gzip -1 > read_map.tsv.gz
        n_parts=\$(ls read_map.d/*.tsv | wc -l)
        rm -rf read_map.d
        echo "[INFO] read map: \$(zcat read_map.tsv.gz | wc -l) reads, \$(zcat read_map.tsv.gz | cut -f2 | sort -u | wc -l) samples from \$n_parts files" >&2
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
    # A failing stage in a pipe must fail the task (see PREPARE_READS).
    set -o pipefail

    # Find the human reads, then drop them -- rather than pushing every read
    # through SAM. The old pipe emitted a SAM record for all ~97 M reads and
    # ran `samtools view -b` twice, each on one core, so ~700 GB went through
    # single-threaded BGZF deflate twice: on 2026-09-23 two samtools pinned at
    # a core each while minimap2 sat blocked on its output at ~6-8 of 128
    # threads, and the stage took 2.8 h at 338 Gbp.
    #
    # PAF prints a line only for reads that hit the reference, which for
    # environmental samples is a sliver of the input. -c keeps base-level
    # alignment so the calls match what the SAM path (-a) made; --secondary=no
    # as before. A read with any line here is one the old -f 4 filter dropped.
    minimap2 -c -x map-ont --secondary=no -t ${task.cpus} \\
        "${params.human_ref}" "${reads}" > human.paf
    cut -f1 human.paf | sort -u > human.ids

    # fastq_filter streams the file once and drops the listed IDs. It opens the
    # file itself, never a pipe (see PREPARE_READS), and reports its counts.
    # Dedup already happened upstream.
    fastq_filter --no_dedupe --exclude_ids human.ids "${reads}" 2> filter.log \\
        | ${writeReads('nohuman_reads.fastq', task.cpus)}
    cat filter.log >&2

    if [ ! -s nohuman_reads.fastq* ]; then
        echo "[ERROR] Human removal produced empty output" >&2
        exit 1
    fi

    input_count=\$(awk '/Total reads:/ {print \$NF}' filter.log)
    removed=\$(awk '/Excluded:/ {print \$NF}' filter.log)
    output_count=\$((input_count - removed))
    echo "[INFO] Human removal: \${input_count} reads in, \${output_count} out, \${removed} removed (\$(( removed * 100 / (input_count + 1) ))%), \$(wc -l < human.paf) alignments" >&2
    """
}
