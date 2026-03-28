// Per-gene depth calculation from BAMs + annotation GFF.
// Moved here from assembly pipeline because it depends on annotation output.
//
// Processes:
//   CALCULATE_GENE_DEPTHS — samtools bedcov per-gene depths from annotation GFF.

process CALCULATE_GENE_DEPTHS {
    tag "gene_depths"
    label 'process_low'
    conda "${projectDir}/conda-envs/dana-mag-binning"
    publishDir "${params.outdir}/mapping", mode: 'copy', enabled: !params.store_dir
    storeDir params.store_dir ? "${params.store_dir}/mapping" : null

    input:
    path(bams)
    path(gff)

    output:
    path("gene_depths.tsv"), emit: gene_depths

    script:
    """
    #!/bin/bash
    set -euo pipefail

    GFF="${gff}"

    if [ ! -s "\${GFF}" ]; then
        echo "[WARNING] GFF is empty — writing empty gene_depths.tsv" >&2
        echo -e "locus_tag\\tcontig\\tstart\\tend\\tmean_depth" > gene_depths.tsv
        exit 0
    fi

    # Convert GFF to BED4 (0-based start): contig start-1 end locus_tag
    # GFF columns: seqid source type start end score strand phase attributes
    # Extract locus_tag from attributes field (;locus_tag=XXX or ID=XXX)
    awk -F'\\t' '!/^#/ && \$3 == "CDS" {
        match(\$9, /locus_tag=([^;]+)/, lt);
        if (lt[1] == "") { match(\$9, /ID=([^;]+)/, lt); }
        if (lt[1] != "") {
            start = \$4 - 1; if (start < 0) start = 0;
            print \$1 "\\t" start "\\t" \$5 "\\t" lt[1]
        }
    }' "\${GFF}" | sort -k1,1 -k2,2n > genes.bed

    if [ ! -s genes.bed ]; then
        echo "[WARNING] No CDS features with locus tags in GFF — writing empty gene_depths.tsv" >&2
        echo -e "locus_tag\\tcontig\\tstart\\tend\\tmean_depth" > gene_depths.tsv
        exit 0
    fi

    echo "[INFO] Found \$(wc -l < genes.bed) CDS features for depth calculation"

    # samtools bedcov: sum of per-base depths per region per BAM
    # bedcov is single-threaded, so split BED and run in parallel
    CPUS=${task.cpus}

    mkdir -p bed_chunks
    awk -v n="\$CPUS" -v dir="bed_chunks" '
        { print > (dir "/chunk_" (NR-1)%n) }
    ' genes.bed

    ls bed_chunks/chunk_* | xargs -P "\$CPUS" -I{} sh -c '
        samtools bedcov "{}" *.sorted.bam > "{}.cov"
    '

    cat bed_chunks/chunk_*.cov > bedcov_raw.tsv

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
    depth_sums = [int(x) for x in parts[4:]]
    n_samples = len(depth_sums)
    if n_samples == 0:
        continue
    mean_depth = sum(s / region_len for s in depth_sums) / n_samples
    # Report 1-based start for consistency with GFF
    print(f'{locus_tag}\\t{contig}\\t{start + 1}\\t{end}\\t{mean_depth:.2f}')
    "  > gene_depths.tsv

    echo "[INFO] Gene depths: \$(tail -n+2 gene_depths.tsv | wc -l) genes"
    """
}
