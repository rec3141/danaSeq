// Read mapping: align each sample to the co-assembly, calculate coverage depths

process MAP_READS {
    tag "${meta.id}"
    label 'process_medium'
    conda "${projectDir}/conda-envs/dana-mag-mapping"
    publishDir "${params.outdir}/mapping", mode: 'link', enabled: !params.store_dir, pattern: '*.{bam,bai}'
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
    conda "${projectDir}/conda-envs/dana-mag-mapping"
    publishDir "${params.outdir}/mapping", mode: 'link', enabled: !params.store_dir
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

process CALCULATE_GENE_DEPTHS {
    tag "gene_depths"
    label 'process_low'
    conda "${projectDir}/conda-envs/dana-mag-mapping"
    publishDir "${params.outdir}/mapping", mode: 'link', enabled: !params.store_dir
    storeDir params.store_dir ? "${params.store_dir}/mapping" : null

    input:
    path(bams)

    output:
    path("gene_depths.tsv"), emit: gene_depths

    script:
    // Discover annotation TSV using same priority as viz.nf
    def storeRoot = params.store_dir ?: ''
    """
    #!/bin/bash
    set -euo pipefail

    # Discover annotation TSV: BAKTA_EXTRA > BAKTA_BASIC > PROKKA
    STORE="${storeRoot}"
    OUT="${params.outdir}"
    find_first() { for f in "\$@"; do [ -f "\${f}" ] && echo "\${f}" && return; done; }
    ANNOT_TSV=\$(find_first \\
        "\${STORE:+\${STORE}/annotation/bakta/extra/annotation.tsv}" \\
        "\${OUT}/annotation/bakta/extra/annotation.tsv" \\
        "\${STORE:+\${STORE}/annotation/bakta/basic/annotation.tsv}" \\
        "\${OUT}/annotation/bakta/basic/annotation.tsv" \\
        "\${STORE:+\${STORE}/annotation/prokka/annotation.tsv}" \\
        "\${OUT}/annotation/prokka/annotation.tsv")

    if [ -z "\${ANNOT_TSV}" ]; then
        echo "[WARNING] No annotation TSV found — writing empty gene_depths.tsv" >&2
        echo -e "locus_tag\\tcontig\\tstart\\tend\\tmean_depth" > gene_depths.tsv
        exit 0
    fi

    # Convert annotation TSV to BED4 (0-based start): contig start-1 end locus_tag
    awk -F'\\t' 'NR>1 && \$2!="region" && \$6!="" {
        start = \$3 - 1; if (start < 0) start = 0;
        print \$1 "\\t" start "\\t" \$4 "\\t" \$6
    }' "\${ANNOT_TSV}" | sort -k1,1 -k2,2n > genes.bed

    if [ ! -s genes.bed ]; then
        echo "[WARNING] No genes with locus tags in annotation — writing empty gene_depths.tsv" >&2
        echo -e "locus_tag\\tcontig\\tstart\\tend\\tmean_depth" > gene_depths.tsv
        exit 0
    fi

    # samtools bedcov: sum of per-base depths per region per BAM
    samtools bedcov genes.bed *.sorted.bam > bedcov_raw.tsv

    # Compute mean depth per gene (average across samples)
    python3 -c "
import sys
header = ['locus_tag', 'contig', 'start', 'end', 'mean_depth']
print('\\t'.join(header))
for line in open('bedcov_raw.tsv'):
    parts = line.rstrip().split('\\t')
    contig, start, end, locus_tag = parts[0], int(parts[1]), int(parts[2]), parts[3]
    region_len = end - start
    if region_len <= 0:
        continue
    # Columns 4+ are depth sums, one per BAM
    depth_sums = [int(x) for x in parts[4:]]
    n_samples = len(depth_sums)
    if n_samples == 0:
        continue
    # Mean across samples: (sum_per_sample / region_len) averaged over samples
    mean_depth = sum(s / region_len for s in depth_sums) / n_samples
    # Report 1-based start for consistency with annotation TSV
    print(f'{locus_tag}\\t{contig}\\t{start + 1}\\t{end}\\t{mean_depth:.2f}')
    "  > gene_depths.tsv

    echo "Gene depths: \$(tail -n+2 gene_depths.tsv | wc -l) genes"
    """
}
