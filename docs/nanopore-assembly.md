# Nanopore Assembly Pipeline

Long-read metagenomic assembly pipeline for Oxford Nanopore data. Preprocesses reads, runs assembly with selectable assemblers, maps reads back to contigs, and produces depth tables for downstream MAG analysis.

## Quick Start

```bash
cd nanopore_assembly
./install.sh && ./install.sh --check

# Local (conda)
./run-nanopore-assembly.sh --input /path/to/reads --outdir /path/to/output

# Choose assembler
./run-nanopore-assembly.sh --input /path/to/reads --outdir /path/to/output --assembler metamdbg

# Apptainer (HPC)
./run-nanopore-assembly.sh --apptainer --input /path/to/reads --outdir /path/to/output

# Docker
./run-nanopore-assembly.sh --docker --input /path/to/reads --outdir /path/to/output
```

## Pipeline Stages

| Step | Process | Tool | Description |
|------|---------|------|-------------|
| 1 | `CONCAT_READS` | cat | Concatenate per-barcode FASTQs into one file per sample |
| 2 | `PREPARE_READS` | BBDuk / Filtlong | Deduplication, size-selection, quality filtering |
| 3 | `REMOVE_HUMAN` | minimap2 | Remove human-mapped reads (optional) |
| 4a | `FLYE_ASSEMBLE` | Flye --meta | Long-read metagenome assembly (default) |
| 4b | `ASSEMBLY_METAMDBG` | metaMDBG | Minimizer-space de Bruijn graph assembly (alternative) |
| 4c | `ASSEMBLY_MYLOASM` | myloasm | Hybrid overlap/DBG assembly (alternative) |
| 5 | `FLYE_POLISH` | Flye | Assembly polishing iterations (optional, default for Flye) |
| 6 | `MAP_READS` | minimap2 + samtools | Per-sample alignment to assembly (flag `-F 0x904`) |
| 7 | `CALCULATE_DEPTHS` | CoverM | Coverage depth table in MetaBAT2 format |
| 8 | `CALCULATE_TNF` | tetramer_freqs (C) | 136-feature tetranucleotide frequency profiles |

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--input` | (required) | Directory containing `*.fastq.gz` or `fastq_pass/barcode*/` structure |
| `--outdir` | `results` | Output directory |
| `--assembler` | `flye` | Assembler: `flye`, `metamdbg`, or `myloasm` |
| `--read_type` | `auto` | Nanopore read type: `auto`, `nano-hq`, `nano-raw`, `nano-corr` |
| `--min_overlap` | `1000` | Flye `--min-overlap` parameter |
| `--polish` | (auto) | Run Flye polisher; auto-enabled for Flye, disabled for others |
| `--dedupe` | `true` | BBDuk deduplication before assembly |
| `--dedupe_memory` | `24 GB` | Memory allocated to deduplication |
| `--filtlong_size` | (skip) | Filtlong target bases (e.g. `40000000000`); skip if not set |
| `--min_barcode_size` | `10485760` | Skip barcodes with fewer than 10 MB of reads |
| `--run_remove_human` | `true` | Remove human reads via minimap2 |
| `--human_ref` | `databases/human_ref/GRCh38_noalt_as.fa.gz` | Path to human reference FASTA or `.mmi` |
| `--assembly_cpus` | `16` | CPUs for assembly processes |
| `--assembly_memory` | `60 GB` | Memory for assembly processes |
| `--store_dir` | (none) | Persistent cache directory (storeDir) |

## Outputs

```
results/
├── assembly/
│   ├── assembly.fasta          Co-assembly
│   ├── assembly_info.txt       Contig metadata (length, coverage, circularity)
│   ├── assembly_graph.gfa      Assembly graph
│   ├── tnf.tsv                 Tetranucleotide frequencies (136 features)
│   └── gc.tsv                  Per-contig GC content
├── mapping/
│   ├── <sample>.sorted.bam     Per-sample alignments
│   ├── <sample>.sorted.bam.bai BAM indices
│   └── depths.txt              CoverM depth matrix (MetaBAT2 format)
└── pipeline_info/
    ├── timeline.html           Execution timeline
    ├── report.html             Resource usage report
    ├── trace.txt               Per-process trace
    ├── dag.html                Pipeline DAG
    └── versions.yml            Tool versions and parameters
```

Feed these outputs into `mag_analysis`:

```bash
cd ../mag_analysis
./run-mag-analysis.sh \
    --assembly /path/to/results/assembly/assembly.fasta \
    --depths /path/to/results/mapping/depths.txt \
    --bam_dir /path/to/results/mapping/ \
    --outdir /path/to/analysis --db_dir /path/to/databases
```

## Profiles

| Profile | Use case |
|---------|----------|
| `standard` | Local execution (default) |
| `test` | Small test data, reduced resources (4 CPUs, 8 GB) |

## Resource Requirements

| Component | CPUs | RAM | Notes |
|-----------|------|-----|-------|
| Flye assembly | 16 | 60 GB | Scales with metagenome complexity |
| metaMDBG assembly | 16 | 60 GB | Lower memory than Flye for large datasets |
| myloasm assembly | 16 | 60 GB | Hybrid overlap/DBG approach |
| minimap2 mapping | 8 | 16 GB | Per-sample; runs in parallel |
| CoverM depths | 2 | 4 GB | Operates on all BAMs |
| Human decontamination | 8 | 16 GB | minimap2 against GRCh38 |

## Input

Nanopore barcode directory structure (MinKNOW output):

```
input_dir/
├── fastq_pass/
│   ├── barcode01/*.fastq.gz
│   ├── barcode02/*.fastq.gz
│   └── ...
```

Or a flat directory of FASTQ files:

```
input_dir/
├── sample1.fastq.gz
├── sample2.fastq.gz
└── ...
```

Barcodes with fewer than `--min_barcode_size` bytes (default 10 MB) are automatically skipped.

## Design Notes

- **Multiple assemblers.** Flye (default), metaMDBG, and myloasm each handle different read quality profiles. Select with `--assembler`.
- **CoverM depth calculation.** Uses CoverM instead of MetaBAT2's `jgi_summarize_bam_contig_depths` to avoid integer overflow with supplementary alignments.
- **Human decontamination.** Enabled by default via minimap2 against GRCh38. Disable with `--run_remove_human false`.
- **Tetranucleotide frequencies.** 136-feature TNF profiles computed by a compiled C binary for downstream binning and visualization.
- **Persistent caching.** Use `--store_dir` to skip completed processes across runs, even after `work/` cleanup.
- **Auto-polishing.** Polishing is automatically enabled for Flye assemblies and disabled for metaMDBG/myloasm. Override with `--polish true/false`.
