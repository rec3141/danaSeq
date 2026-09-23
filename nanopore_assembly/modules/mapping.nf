// Read mapping: align each sample to the co-assembly, calculate coverage depths.
//
// Processes:
//   MAP_READS            — minimap2 map-ont per sample, samtools sort + index.
//                          Drops unmapped + secondary (-F 0x104), keeps supplementary
//                          for read-bridged adjacency. CoverM ignores supplementary
//                          when computing depths so they don't inflate coverage.
//   CALCULATE_DEPTHS     — CoverM metabat-mode depth table across all BAMs.
//                          Replaces jgi_summarize_bam_contig_depths (overflow bug).
//
// Note: CALCULATE_GENE_DEPTHS lives in mag_analysis (depends on annotation).

process MAP_READS {
    tag "${meta.id}"
    label 'process_medium'
    // minimap2's peak RSS is set by the INDEX BATCH, not the whole reference:
    // with --split-prefix it builds the index in -I sized chunks (4 Gbp by
    // default), so an 8.9 Gbp assembly still peaks at ~48 GB, not ~90 GB.
    // Measured: 48.05 GB peak on a 8.9 Gbp / 1.18M-contig marine assembly.
    //
    // process_medium declares 16 GB. Nextflow's local executor derives
    // concurrency from the DECLARED value, so it packed ~3x too many tasks onto
    // the node, they exhausted 750 GB between them, and samtools sort -- the
    // last stage to allocate -- died with exit 1. 116 task failures across 39
    // of 276 samples, each retried once and then silently ignored.
    memory { def refGb  = (assembly.size() / (1024L**3)) as double
             def batch  = Math.min(refGb, 4.0d)          // minimap2 -I default
             def needed = Math.ceil(batch * 10.0d + 16.0d) as int
             (Math.max(needed, 24) * task.attempt).GB }
    conda "${projectDir}/conda-envs/dana-mag-assembly"
    publishDir "${params.outdir}/mapping", mode: 'copy', enabled: !params.store_dir, pattern: '*.{bam,bai}'
    storeDir params.store_dir ? "${params.store_dir}/mapping" : null

    input:
    tuple val(meta), path(fastq), path(assembly)

    output:
    tuple val(meta), path("${meta.id}.sorted.bam"), path("${meta.id}.sorted.bam.bai"), emit: bam

    script:
    """
    # Refuse to map against an empty-reference fasta. A header-only fasta would
    # produce a 0-@SQ BAM that minimap2/samtools still exit 0 on (the original
    # silent-failure mode), with downstream coverm panicking on the result.
    if [ ! -s "${assembly}" ]; then
        echo "[ERROR] Reference assembly is empty: ${assembly}" >&2
        exit 1
    fi
    REF_CONTIGS=\$(grep -c '^>' "${assembly}" || true)
    if [ "\${REF_CONTIGS:-0}" -lt 1 ]; then
        echo "[ERROR] Reference assembly contains 0 contigs: ${assembly}" >&2
        exit 1
    fi

    # -F 0x104: drop unmapped (0x4) and secondary (0x100), keep supplementary (0x800)
    # Supplementary alignments are kept for read-bridged adjacency (cross-contig links)
    # CoverM's metabat method ignores supplementary alignments, so depths are unchanged
    #
    # --split-prefix: assemblies >4GB (default -I) trigger a multi-part minimap2
    # index. WITHOUT --split-prefix, minimap2 silently drops @SQ records from the
    # SAM output ("For a multi-part index, no \\@SQ lines will be outputted"),
    # samtools view then errors with "no SQ lines present", and the resulting
    # BAM has 0 references / 0 reads while exiting 0 — the exact silent failure
    # that wasted 6.5 days of the original myloasm run. Always pass it.
    minimap2 -a -x map-ont --secondary=no -t ${task.cpus} \\
        --split-prefix "${meta.id}_split" \\
        "${assembly}" "${fastq}" \\
        | samtools view -u -F 0x104 \\
        | samtools sort -@ ${task.cpus} -o "${meta.id}.sorted.bam" -

    samtools index -@ ${task.cpus} "${meta.id}.sorted.bam"

    # A BAM with no alignments is a silent failure downstream: coverm emits a
    # zero column and the sample vanishes from the depth matrix without comment.
    N_ALN=\$(samtools view -c "${meta.id}.sorted.bam")
    if [ "\${N_ALN:-0}" -eq 0 ]; then
        echo "[ERROR] ${meta.id}: 0 alignments in output BAM" >&2
        exit 1
    fi
    echo "[INFO] ${meta.id}: \${N_ALN} alignments" >&2

    # Validate BAM: must contain @SQ records matching the reference, else minimap2
    # mapped against an empty index (the bug we're guarding against).
    BAM_SQ=\$(samtools view -H "${meta.id}.sorted.bam" | grep -c '^@SQ' || true)
    if [ "\${BAM_SQ:-0}" -lt 1 ]; then
        echo "[ERROR] BAM for ${meta.id} has 0 \\@SQ records (reference had \${REF_CONTIGS} contigs); minimap2 saw an empty reference" >&2
        exit 1
    fi
    """
}

process CALCULATE_DEPTHS {
    tag "depths"
    label 'process_low'
    conda "${projectDir}/conda-envs/dana-mag-assembly"
    publishDir "${params.outdir}/mapping", mode: 'copy', enabled: !params.store_dir
    storeDir params.store_dir ? "${params.store_dir}/mapping" : null

    input:
    path(bams)
    path(assembly)
    val(expected_ids)

    output:
    path("depths.txt"),            emit: jgi_depth
    path("unmapped_samples.txt"),  emit: unmapped, optional: true

    script:
    // Built here, not inside the script: a '\n' in an interpolated expression
    // is unescaped by the GString lexer before the expression is parsed.
    def expected_list = expected_ids.join(System.lineSeparator())
    """
    # MAP_READS runs under errorStrategy 'ignore' (a failed mapping must not
    # take down a multi-day assembly), which means a sample can disappear from
    # the depth matrix with no error anywhere. coverm then builds the matrix
    # from whatever BAMs exist and MetaBAT/SemiBin bin on missing columns.
    # Reconcile what was asked for against what arrived, and say so loudly.
    cat > expected.txt <<'EXPECTED_EOF'
${expected_list}
EXPECTED_EOF
    ls *.sorted.bam 2>/dev/null | sed 's/\\.sorted\\.bam\$//' | sort -u > got.txt
    sort -u expected.txt > exp.txt
    comm -23 exp.txt got.txt > unmapped_samples.txt || true
    N_EXP=\$(wc -l < exp.txt); N_GOT=\$(wc -l < got.txt); N_MISS=\$(wc -l < unmapped_samples.txt)
    if [ "\$N_MISS" -gt 0 ]; then
        echo "========================================================" >&2
        echo "[WARNING] \$N_MISS of \$N_EXP samples have NO BAM and are ABSENT" >&2
        echo "          from the depth matrix. Binning will run without their" >&2
        echo "          differential coverage. See mapping/unmapped_samples.txt" >&2
        sed 's/^/          /' unmapped_samples.txt >&2
        echo "========================================================" >&2
    else
        rm -f unmapped_samples.txt
        echo "[INFO] all \$N_GOT expected samples present in the depth matrix" >&2
    fi

    # CoverM handles supplementary alignments correctly and avoids the integer
    # overflow bug in jgi_summarize_bam_contig_depths (MetaBAT2 <=2.17)
    coverm contig \\
        -b *.sorted.bam \\
        --methods metabat \\
        --min-read-percent-identity 80 \\
        --min-read-aligned-percent 0 \\
        --threads ${task.cpus} \\
        --output-file depths.txt

    if [ ! -s depths.txt ]; then
        echo "[ERROR] Depth calculation produced empty output" >&2
        exit 1
    fi
    """
}

