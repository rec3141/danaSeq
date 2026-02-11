# MAG Assembly Pipeline

Metagenome-assembled genome (MAG) reconstruction from Oxford Nanopore long reads. Co-assembles reads with Flye, maps back with minimap2, runs three binning algorithms in parallel (SemiBin2, MetaBAT2, MaxBin2), and integrates results with DAS Tool consensus.

## Quick Start

### Nextflow (recommended)

```bash
cd nextflow

# Install conda environments (~15 min first time)
./install.sh
./install.sh --check

# Run the pipeline
mamba run -p conda-envs/dana-mag-flye \
    nextflow run main.nf --input /path/to/reads -resume

# Show all options
mamba run -p conda-envs/dana-mag-flye \
    nextflow run main.nf --help
```

### Launcher script

```bash
cd nextflow

# Local conda (default)
./run-mag.sh --input /path/to/reads --outdir /path/to/output

# Docker mode
docker build -t danaseq-mag .
./run-mag.sh --docker --input /path/to/reads --outdir /path/to/output
```

## Pipeline Overview

```
Sample FASTQs (N files)
         |  collect()
   ASSEMBLY_FLYE           All reads -> 1 co-assembly (Flye --meta)
         |
   MAP_READS (xN)          Per-sample alignment (minimap2, -F 0x904)
         |  collect()
   CALCULATE_DEPTHS         Coverage depth table (CoverM)
         |
    +---------+---------+
    |         |         |
 SemiBin2  MetaBAT2  MaxBin2   Three binners in parallel
    |         |         |
    +---------+---------+
         |
   DASTOOL_CONSENSUS        Best bin per contig
         |
   results/binning/consensus/
```

## Output

```
results/
├── assembly/
│   └── assembly.fasta            Co-assembly
├── mapping/
│   ├── *.sorted.bam              Per-sample alignments
│   ├── *.sorted.bam.bai
│   └── depths.txt                CoverM depth table
├── binning/
│   ├── semibin/contig_bins.tsv   SemiBin2 assignments
│   ├── metabat/contig_bins.tsv   MetaBAT2 assignments
│   ├── maxbin/contig_bins.tsv    MaxBin2 assignments
│   └── consensus/
│       ├── bins/*.fa             Final consensus MAG FASTAs
│       ├── contig2bin.tsv        Contig-to-bin assignments
│       ├── allbins.fa            All bins concatenated
│       └── scores.tsv            Per-bin quality scores
└── pipeline_info/
    ├── timeline.html
    ├── report.html
    └── trace.txt
```

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--input` | (required) | Directory containing `*.fastq.gz` files |
| `--outdir` | `results` | Output directory |
| `--min_overlap` | `1000` | Flye `--min-overlap` |
| `--dedupe` | `false` | BBDuk deduplication before assembly |
| `--filtlong_size` | (skip) | Filtlong target bases (e.g. `40000000000`) |
| `--run_maxbin` | `true` | Include MaxBin2 in consensus |
| `--metabat_min_cls` | `50000` | MetaBAT2 minimum cluster size |
| `--assembly_cpus` | `24` | CPUs for assembly |
| `--assembly_memory` | `64 GB` | Memory for assembly |

### Profiles

| Profile | Use case |
|---------|----------|
| `standard` | Local execution (default) |
| `test` | Small test data, reduced resources |
| `shipboard` | Production: 32 CPUs, 256 GB RAM |

## Design Notes

**CoverM for depth calculation.** Replaces `jgi_summarize_bam_contig_depths`, which has an integer overflow bug in MetaBAT2 <=2.17 when processing long-read BAMs. CoverM handles supplementary alignments correctly.

**Supplementary alignment filtering.** The mapping step uses `samtools view -F 0x904` to drop unmapped, secondary, and supplementary alignments. Chimeric long reads produce supplementary records that cause massive depth overcounting.

**Dynamic binner architecture.** Each binner emits `[label, file]` tuples that are mixed and collected for DAS Tool. Adding a new binner requires only a process definition and one line in `main.nf`.

**GPU-accelerated SemiBin2.** The local conda env includes `pytorch-gpu` for GPU acceleration. The Docker image uses CPU-only PyTorch to keep the image small (~7 GB vs ~12 GB).

**Graceful failure handling.** SemiBin2 and DAS Tool handle edge cases (0 bins, no bins above score threshold) without crashing the pipeline.

## Remaining Bash Scripts

These scripts have not yet been ported to Nextflow and remain as standalone tools:

### Polishing (40s)

| Script | Description |
|--------|-------------|
| `40_polish_assemblies.sh` | Racon (2 rounds) + Medaka polishing |
| `41_medaka_split.sh` | Parallelized Medaka for speed |

### Bin Characterization (50s)

| Script | Description |
|--------|-------------|
| `50_run_kraken_all.sh` | Kraken2 taxonomic classification of bins |
| `51_sketch_bins.sh` | Sendsketch species-level identification |
| `52_tetra_frequency.sh` | Tetranucleotide composition profiles |

### Visualization (70s-80s)

| Script | Description |
|--------|-------------|
| `70_kraken2_to_bandage.sh` | Convert Kraken output for Bandage |
| `71_circular_gfa.awk` | Handle circular contigs in GFA format |
| `80_plot_bins.R` | PCA, t-SNE, heatmaps of MAG composition |
| `80_plot_bins_save.R` | Variant with file output |
| `81_plot_mapping.R` | Coverage plots and mapping statistics |
| `82_inter_binning_analysis.r` | t-SNE/UMAP projections, graph clustering |

### Integration (90s)

| Script | Description |
|--------|-------------|
| `90_edna_schema.r` | Database schema for eDNA integration |
| `91_foam_ecoserv.r` | FOAM ecosystem services analysis |

## MAG Quality Standards (MIMAG)

| Tier | Completeness | Contamination | Additional |
|------|-------------|---------------|------------|
| High quality | >90% | <5% | 23S, 16S, 5S rRNA + tRNAs |
| Medium quality | >50% | <10% | -- |
| Low quality | <50% | <10% | -- |

## Directory Structure

```
20_mag_assembly/
├── nextflow/              Nextflow pipeline (primary interface)
│   ├── main.nf
│   ├── modules/
│   ├── envs/
│   ├── nextflow.config
│   ├── Dockerfile
│   └── install.sh
├── 40_polish_assemblies.sh    Not yet ported
├── ...                        (see table above)
├── archive/                   Replaced bash scripts (reference only)
│   ├── 10_assembly_flye.sh
│   ├── 20_mapping.sh
│   ├── 30_binning_semibin.sh
│   ├── 61_map_and_bin_optimized.sh
│   └── ...
├── CLAUDE.md
└── README.md
```

## References

- Flye: Kolmogorov et al., Nature Biotechnology 2019
- SemiBin2: Pan et al., Nature Communications 2023
- MetaBAT2: Kang et al., PeerJ 2019
- DAS Tool: Sieber et al., Nature Microbiology 2018
- CheckM2: Chklovski et al., Nature Methods 2023
- MIMAG: Bowers et al., Nature Biotechnology 2017
- CoverM: [github.com/wwood/CoverM](https://github.com/wwood/CoverM)
