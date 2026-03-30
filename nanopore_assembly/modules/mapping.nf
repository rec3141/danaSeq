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
    conda "${projectDir}/conda-envs/dana-mag-assembly"
    publishDir "${params.outdir}/mapping", mode: 'copy', enabled: !params.store_dir, pattern: '*.{bam,bai}'
    storeDir params.store_dir ? "${params.store_dir}/mapping" : null

    input:
    tuple val(meta), path(fastq), path(assembly)

    output:
    tuple val(meta), path("${meta.id}.sorted.bam"), path("${meta.id}.sorted.bam.bai"), emit: bam

    script:
    """
    # -F 0x104: drop unmapped (0x4) and secondary (0x100), keep supplementary (0x800)
    # Supplementary alignments are kept for read-bridged adjacency (cross-contig links)
    # CoverM's metabat method ignores supplementary alignments, so depths are unchanged
    minimap2 -a -x map-ont --secondary=no -t ${task.cpus} \\
        "${assembly}" "${fastq}" \\
        | samtools view -b -F 0x104 \\
        | samtools sort -@ ${task.cpus} -o "${meta.id}.sorted.bam" -

    samtools index -@ ${task.cpus} "${meta.id}.sorted.bam"

    # Validate BAM
    if [ ! -s "${meta.id}.sorted.bam" ]; then
        echo "[ERROR] Mapping produced empty BAM for ${meta.id}" >&2
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

    output:
    path("depths.txt"), emit: jgi_depth

    script:
    """
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

