# Real-Time Nanopore Pipeline

Shipboard metagenomic analysis for Oxford Nanopore sequencing data. Processes reads as they stream from MinKNOW, providing live taxonomic classification, gene annotation, and functional profiling.

## Quick Start

```bash
cd nanopore_live/nextflow

# Install conda environments (~15 min first time)
./install.sh
./install.sh --check

# Run with all modules
./run-realtime.sh --input /path/to/nanopore/run --outdir /path/to/output \
    --run_kraken --kraken_db /path/to/krakendb \
    --run_prokka --run_sketch --run_tetra

# Docker
docker build -t danaseq-realtime .
./run-realtime.sh --docker --input /path/to/nanopore/run --outdir /path/to/output \
    --run_kraken --kraken_db /path/to/krakendb --run_prokka

# Show all options
./run-realtime.sh --help
```

## Pipeline Overview

```
MinKNOW FASTQ
    |
VALIDATE_FASTQ          Gzip integrity + BBMap repair
    |
QC_FASTQ_FILTER         Length + quality filtering (C, streaming)
    |
CONVERT_TO_FASTA        Header cleanup
    |
    +---> KRAKEN2_CLASSIFY     Taxonomic classification (batched per sample)
    +---> PROKKA_ANNOTATE      ORF prediction + annotation
    +---> HMM_SEARCH           HMMER3 against user-supplied databases
    +---> SENDSKETCH           Rapid taxonomic sketching
    +---> TETRAMER_FREQ        Tetranucleotide composition (C)
              |
         DB_INTEGRATION        Load results into DuckDB
              |
         CLEANUP               Compress/delete source files (optional)
```

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--input` | (required) | Directory **containing** `fastq_pass/` |
| `--outdir` | `results` | Output directory |
| `--run_kraken` | `false` | Enable Kraken2 classification |
| `--kraken_db` | (required if kraken) | Path to Kraken2 database |
| `--run_prokka` | `false` | Enable Prokka annotation |
| `--hmm_databases` | (skip) | Path to HMM file(s) for functional profiling |
| `--run_sketch` | `false` | Enable Sendsketch profiling |
| `--run_tetra` | `false` | Enable tetranucleotide frequency |
| `--watch` | `false` | Monitor for new files during live sequencing |
| `--watch_glob` | `*/fastq_pass/barcode*/*.fastq.gz` | Glob pattern for watch mode |
| `--db_sync_minutes` | `10` | DuckDB sync interval in watch mode |
| `--run_db_integration` | `false` | Load results into DuckDB |
| `--cleanup` | `false` | Compress/delete files after DB import |

### Profiles

| Profile | Use case |
|---------|----------|
| `standard` | Local execution, auto-detect resources |
| `test` | Small test files, reduced resources |
| `shipboard` | Production: 32 CPUs, 256 GB RAM, 100 GB Kraken |

## Watch Mode

Monitor a directory for new FASTQ files during active sequencing:

```bash
nextflow run main.nf --input /path/to/runs \
    --watch --db_sync_minutes 10 \
    --run_kraken --kraken_db /path/to/db \
    --run_prokka --run_db_integration
```

`DB_SYNC` runs as a long-lived process that periodically loads new results into DuckDB. R scripts are idempotent and track imports via `import_log`.

## Output

```
results/
├── FLOWCELL/
│   ├── barcode01/
│   │   ├── fa/           FASTA (quality-filtered sequences)
│   │   ├── fq/           Intermediate FASTQ (BBDuk output)
│   │   ├── kraken/       *.tsv (per-read), *.report (summary)
│   │   ├── prokka/       SAMPLE/PROKKA_* (GFF, FAA, FFN, TSV)
│   │   ├── hmm/          *.DBNAME.tsv, *.DBNAME.tbl
│   │   ├── sketch/       Sendsketch profiles
│   │   ├── tetra/        Tetranucleotide frequencies
│   │   └── log.txt
│   └── ...
├── dana.duckdb           Integrated database
└── pipeline_info/        Nextflow reports (timeline, trace, DAG)
```

## DuckDB Integration

R scripts in `nextflow/bin/` load results into DuckDB:

| Script | Data |
|--------|------|
| `40_kraken_db.r` | Kraken2 classifications |
| `41_krakenreport_db.r` | Kraken2 summary reports |
| `42_prokka_db.r` | Prokka annotations |
| `43_sketch_db.r` | Sendsketch profiles |
| `44_tetra_db.r` | Tetranucleotide frequencies |
| `45_stats_db.r` | Assembly statistics |
| `46_log_db.r` | Processing logs |
| `47_merge_db.r` | Merge per-run databases |

## Design Notes

- **Kraken2 batching.** In batch mode, FASTAs are grouped per sample so the DB loads once per barcode. In watch mode, files stream individually. Uses `maxForks = 1` since kraken2 loads 50-100 GB into RAM.
- **Compiled C tools.** `fastq_filter` (QC) and `tetramer_freqs` (TNF) are compiled C binaries, replacing filtlong and a Python script.
- **Watch mode.** Uses `Channel.watchPath()` with a flat glob pattern (Java WatchService limitation -- no recursive `**`).
- **Conda environments.** Four isolated environments avoid dependency conflicts: dana-bbmap, dana-prokka, dana-bakta, dana-tools.

## Input

Oxford Nanopore directory structure with multiplexed barcodes:

```
input_dir/fastq_pass/
├── barcode01/*.fastq.gz
├── barcode02/*.fastq.gz
└── ...
```
