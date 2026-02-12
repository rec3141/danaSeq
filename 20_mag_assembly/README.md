# MAG Assembly Pipeline

Metagenome-assembled genome (MAG) reconstruction from Oxford Nanopore long reads. Co-assembles reads with Flye, maps back with minimap2, runs five binning algorithms (SemiBin2, MetaBAT2, MaxBin2, LorBin, COMEBin), and integrates results with DAS Tool consensus.

## Quick Start

```bash
cd nextflow

# Install conda environments (~15 min first time)
./install.sh
./install.sh --check

# Run (local conda, handles activation automatically)
./run-mag.sh --input /path/to/reads --outdir /path/to/output

# Or with Docker
docker build -t danaseq-mag .
./run-mag.sh --docker --input /path/to/reads --outdir /path/to/output

# Show all options
./run-mag.sh --help
```

### Kitchen sink — all options with defaults

```bash
cd nextflow
./run-mag.sh --input /data/reads --outdir /data/output \
    --dedupe \
    --filtlong_size 40000000000 \
    --min_overlap 1000 \
    --run_maxbin true \
    --run_lorbin true \
    --run_comebin true \
    --lorbin_min_length 80000 \
    --metabat_min_cls 50000 \
    --checkm2_db /path/to/checkm2_db \
    --assembly_cpus 24 \
    --assembly_memory '64 GB'
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
    +---------+---------+---------+---------+
    |         |         |         |         |
 SemiBin2  MetaBAT2  MaxBin2  LorBin   COMEBin   Five binners (serial)
    |         |         |         |         |
    +---------+---------+---------+---------+
         |
   DASTOOL_CONSENSUS        Best bin per contig
         |
   CHECKM2                  Quality assessment (optional, needs --checkm2_db)
         |
   results/binning/
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
│   ├── lorbin/contig_bins.tsv    LorBin assignments
│   ├── comebin/contig_bins.tsv   COMEBin assignments
│   ├── dastool/
│   │   ├── bins/*.fa             Final consensus MAG FASTAs
│   │   ├── contig2bin.tsv        Contig-to-bin assignments
│   │   ├── allbins.fa            All bins concatenated
│   │   ├── bin_quality.tsv       DAS Tool SCG-based quality scores
│   │   └── summary.tsv           Consensus winners with scores
│   └── checkm2/
│       └── quality_report.tsv    CheckM2 completeness/contamination (if --checkm2_db)
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
| `--run_lorbin` | `true` | Include LorBin in consensus (deep learning, long-read) |
| `--run_comebin` | `true` | Include COMEBin in consensus (contrastive learning) |
| `--lorbin_min_length` | `80000` | LorBin `--bin_length` minimum (bp) |
| `--metabat_min_cls` | `50000` | MetaBAT2 minimum cluster size |
| `--checkm2_db` | (skip) | Path to CheckM2 DIAMOND database |
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

**GPU-accelerated ML binners.** The local conda envs include `pytorch-gpu` for GPU-accelerated SemiBin2, LorBin, and COMEBin. The Docker image uses CPU-only PyTorch to keep the image small.

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
- LorBin: Gao et al., Briefings in Bioinformatics 2024
- COMEBin: Xie et al., Nature Communications 2024
- DAS Tool: Sieber et al., Nature Microbiology 2018
- CheckM2: Chklovski et al., Nature Methods 2023
- MIMAG: Bowers et al., Nature Biotechnology 2017
- CoverM: [github.com/wwood/CoverM](https://github.com/wwood/CoverM)
