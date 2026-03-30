# Real-Time Nanopore Pipeline

Real-time metagenomic analysis for Oxford Nanopore sequencing data. Processes reads as they stream from MinKNOW, providing live taxonomic classification, gene annotation, and functional profiling.

## Quick Start

```bash
cd nanopore_live
./install.sh && ./install.sh --check

./run-realtime.sh --input /path/to/nanopore/run --outdir /path/to/output \
    --run_kraken --kraken_db /path/to/krakendb \
    --run_prokka --run_sketch --run_tetra
```

## Pipeline Stages

| Step | Process | Tool | Description |
|------|---------|------|-------------|
| 1 | `VALIDATE_FASTQ` | BBMap | Gzip integrity check + read repair |
| 2 | `QC_BBDUK` | BBDuk | Adapter and quality trimming |
| 3 | `QC_FASTQ_FILTER` | fastq_filter (C) | Length + quality filtering (streaming) |
| 4 | `CONVERT_TO_FASTA` | awk | Header cleanup and format conversion |
| 5a | `KRAKEN2_CLASSIFY` | Kraken2 | Taxonomic classification (batched per sample) |
| 5b | `PROKKA_ANNOTATE` | Prokka | ORF prediction + functional annotation |
| 5c | `BAKTA_CDS` / `BAKTA_FULL` | Bakta | CDS-only or full annotation (alternative to Prokka) |
| 5d | `HMM_SEARCH` | HMMER3 | Profile HMM search against user-supplied databases |
| 5e | `SENDSKETCH` | BBTools | Rapid taxonomic sketching via MinHash |
| 5f | `TETRAMER_FREQ` | tetramer_freqs (C) | Tetranucleotide composition profiles |
| 6 | `DB_INTEGRATION` | R + DuckDB | Load all results into DuckDB |
| 7 | `DB_SYNC` | R + DuckDB | Periodic sync during watch mode |
| 8 | `CLEANUP` | bash | Compress/delete source files after DB import |

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--input` | (required) | Directory **containing** `fastq_pass/` |
| `--outdir` | `results` | Output directory |
| `--run_kraken` | `false` | Enable Kraken2 classification |
| `--kraken_db` | (required if kraken) | Path to Kraken2 database |
| `--run_prokka` | `false` | Enable Prokka annotation |
| `--annotator` | `bakta` | Gene annotator: `prokka`, `bakta`, or `none` |
| `--bakta_db` | (required if bakta) | Path to Bakta database |
| `--bakta_full` | `false` | Run full Bakta annotation (ncRNA/tRNA/CRISPR) |
| `--hmm_databases` | (skip) | Comma-delimited paths to HMM files |
| `--run_sketch` | `false` | Enable Sendsketch profiling |
| `--run_tetra` | `false` | Enable tetranucleotide frequency |
| `--watch` | `false` | Monitor for new files during live sequencing |
| `--watch_glob` | `*/fastq_pass/barcode*/*.fastq.gz` | Glob pattern for watch mode |
| `--db_sync_minutes` | `10` | DuckDB sync interval in watch mode |
| `--run_db_integration` | `false` | Load results into DuckDB |
| `--cleanup` | `false` | Compress/delete files after DB import |
| `--min_readlen` | `1500` | Minimum read length after filtering |
| `--keep_percent` | `80` | Percent of reads to keep (by quality) |
| `--min_file_size` | `1000000` | Minimum FASTQ file size in bytes (1 MB) |
| `--store_dir` | (none) | Persistent cache directory (storeDir) |

## Outputs

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

## Profiles

| Profile | Use case |
|---------|----------|
| `standard` | Local execution |
| `test` | Small test files, reduced resources |

## Resource Requirements

| Component | CPUs | RAM | Notes |
|-----------|------|-----|-------|
| Kraken2 | 8 | 50-100 GB | Depends on database size; serialized via `maxForks=1` |
| Prokka/Bakta | 4 | 8 GB | Per-file annotation |
| HMMER3 | 4 | 4 GB | Per-file search |
| DuckDB integration | 2 | 4 GB | R-based import scripts |

## Watch Mode

Monitor a directory for new FASTQ files during active sequencing:

```bash
./run-realtime.sh --input /path/to/runs --outdir /path/to/output \
    --watch --db_sync_minutes 10 \
    --run_kraken --kraken_db /path/to/db \
    --run_prokka --run_db_integration
```

`DB_SYNC` runs as a long-lived process that periodically loads new results into DuckDB. R scripts are idempotent and track imports via `import_log`.

## Post-DB Cleanup

The `--cleanup` flag compresses or deletes source files after DuckDB import:

| Directory/Files | Action |
|----------------|--------|
| `fa/*.fa` | Gzip in place (kept as compressed backup) |
| `kraken/`, `sketch/`, `tetra/` | Delete (data lives in DuckDB) |
| `prokka/*/PROKKA_*.tsv` | Delete (loaded into DuckDB) |
| `prokka/*/PROKKA_*.gff`, `.faa`, `.ffn` | Gzip in place |
| `hmm/`, `dana.duckdb`, `log.txt` | Kept (not cleaned) |

Safe for watch mode -- operates per-file and checks `import_log` before deleting.

## DuckDB Integration

R scripts in `bin/` load results into DuckDB:

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

## Input

Oxford Nanopore directory structure with multiplexed barcodes:

```
input_dir/fastq_pass/
├── barcode01/*.fastq.gz
├── barcode02/*.fastq.gz
└── ...
```

## Design Notes

- **Kraken2 batching.** In batch mode, FASTAs are grouped per sample (`groupTuple`) so the DB loads once per barcode. In watch mode, files stream individually. Uses `maxForks = 1` since Kraken2 loads 50-100 GB into RAM.
- **Compiled C tools.** `fastq_filter` (QC) and `tetramer_freqs` (TNF) are compiled C binaries in `bin/`, replacing filtlong and a Python script.
- **Watch mode.** Uses `Channel.watchPath()` with a flat glob pattern (Java WatchService limitation -- no recursive `**`).
- **Conda environments.** Four isolated environments avoid dependency conflicts: dana-bbmap, dana-prokka, dana-bakta, dana-tools.
