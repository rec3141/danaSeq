# Assembler Comparison: Flye vs metaMDBG vs myloasm

Test data: 2.7 GB nanopore reads, 8 barcodes (Ebb & Flow 2025-10-28)

## Assembly Statistics

| Metric | metaMDBG | Flye | myloasm |
|--------|----------|------|---------|
| Contigs | 11,091 | 15,222 | 21,981 |
| Total size | 94.4 Mbp | 175.5 Mbp | 286.7 Mbp |
| Largest contig | 1.73 Mbp | 830 Kbp | 2.30 Mbp |
| N50 | 10,384 | 16,010 | 13,658 |
| Circular contigs | 28 | 252 | 106 |
| Assembly time | 31 min | 1h 10m | 45 min |

## DAS Tool Binning Results

| Metric | metaMDBG | Flye | myloasm |
|--------|----------|------|---------|
| Total bins | 26 | 42 | 56 |
| High-quality (>=90% comp, <5% redun) | 5 | 15 | 10 |
| Medium-quality (>=50% comp, <10% redun) | 18 | 26 | 44 |
| Total >=50% complete | 26 | 42 | 56 |

Five binners were used in all runs: SemiBin2, MetaBAT2, MaxBin2, LorBin, COMEBin.

### Top Binner per Assembler

| Assembler | Top binner | Bins contributed |
|-----------|-----------|-----------------|
| metaMDBG | LorBin | 8 |
| Flye | COMEBin | 16 |
| myloasm | COMEBin | 29 |

## Total Pipeline Runtime

| Assembler | Time | CPUs / RAM |
|-----------|------|------------|
| metaMDBG | 1h 7m | 4 CPU / 8 GB |
| myloasm | 2h 11m | 4 CPU / 8 GB |
| Flye | 1h 58m | 16 CPU / 60 GB |

metaMDBG and myloasm ran with reduced resources (test profile); Flye required standard resources (crashed at 8 GB during polishing).

## Pairwise Assembly Overlap

Assemblies were aligned pairwise with minimap2 (asm5 mode, primary alignments only, contigs >= 5 Kbp). The table shows what fraction of the query assembly aligns to the reference.

| Query -> Reference | Contigs | Mbp aligned / total | % aligned |
|----|----|----|-----|
| flye -> metamdbg | 8,771 / 8,771 (100%) | 99.2 / 145.2 | 68.3% |
| metamdbg -> flye | 6,310 / 6,310 (100%) | 73.6 / 77.1 | 95.4% |
| flye -> myloasm | 9,646 / 9,646 (100%) | 132.9 / 153.0 | 86.8% |
| myloasm -> flye | 17,555 / 17,555 (100%) | 201.3 / 250.8 | 80.3% |
| metamdbg -> myloasm | 6,292 / 6,292 (100%) | 70.2 / 77.2 | 90.9% |
| myloasm -> metamdbg | 16,343 / 16,343 (100%) | 153.0 / 238.7 | 64.1% |

metaMDBG is nearly entirely contained in both Flye (95.4%) and myloasm (90.9%). myloasm has the most unique sequence: ~86 Mbp not in metaMDBG and ~50 Mbp not in Flye.

## t-SNE + k-means Clustering

t-SNE on tetranucleotide frequencies (37,811 contigs >= 5 Kbp) with k=128 k-means clustering shows:

- **37 clusters** dominated (>66%) by myloasm, totaling 9,944 contigs / 120.6 Mbp
- **0 clusters** dominated by Flye or metaMDBG alone
- **47 clusters** shared roughly equally (no assembler >50%), totaling 228.6 Mbp
- **1 near-exclusive cluster** (>90% myloasm): 246 contigs / 2.5 Mbp

sendsketch against NCBI nt on a myloasm-dominated cluster (136 contigs near myloasm_contig_06013) returned zero hits for 135/136 contigs — this cluster represents novel, low-coverage organisms that only myloasm recovered.

## Takeaways

- **Flye** produces the most high-quality bins (15 HQ) and best N50 (16 Kbp), with the most circular contigs (252). Best for quality-focused analysis. Requires the most memory.
- **myloasm** recovers the most total bins (56) and the largest assembly (287 Mbp) with the largest single contig (2.3 Mbp). Best for capturing diversity.
- **metaMDBG** is fastest (31 min assembly) and most memory-efficient, but produces the fewest bins (26). Good for quick surveys or resource-constrained environments.
