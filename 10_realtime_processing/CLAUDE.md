# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

**Important**: If user poses a question, agent responds with an answer, not a codebase change.

## Project Overview

**dānaSeq** is a real-time metagenomic analysis pipeline for Oxford Nanopore sequencing on oceanographic expeditions. It processes DNA reads as they stream from the sequencer, providing live taxonomic classification, gene annotation, and interactive geographic visualization.

The pipeline is transitioning from standalone bash scripts to a **Nextflow DSL2** implementation. The Nextflow pipeline in `nextflow/` is the primary interface; legacy bash scripts are preserved in `archive/` for reference.

## Repository Structure

```
dānaSeq/
├── 10_realtime_processing/       # This directory
│   ├── nextflow/                 # Primary pipeline (Nextflow DSL2)
│   │   ├── main.nf              # Pipeline entry point
│   │   ├── modules/             # Process definitions
│   │   │   ├── validate.nf      #   FASTQ integrity checks
│   │   │   ├── qc.nf            #   BBDuk, Filtlong, FASTA conversion
│   │   │   ├── kraken.nf        #   Kraken2 classification
│   │   │   ├── prokka.nf        #   Gene annotation
│   │   │   ├── hmm.nf           #   HMM functional gene search
│   │   │   ├── sketch.nf        #   Sendsketch profiling
│   │   │   ├── tetramer.nf      #   Tetranucleotide frequency
│   │   │   └── db_integration.nf #  DuckDB loading + cleanup
│   │   ├── bin/                 # Helper scripts called by processes
│   │   │   ├── 30_kraken_parse.awk
│   │   │   ├── 4X_*_db.r       #   DuckDB integration R scripts
│   │   │   └── 60_edna_mapping_viz.r
│   │   ├── nextflow.config      # Params, profiles, resources
│   │   ├── conda-envs/          # Pre-built conda environments
│   │   ├── envs/                # Conda YAML specs
│   │   ├── Dockerfile           # Self-contained Docker image
│   │   ├── entrypoint.sh        # Docker entrypoint
│   │   ├── run-realtime.sh      # Pipeline launcher (local/Docker)
│   │   ├── install.sh           # Conda env installer
│   │   └── test-data/           # Bundled test data
│   ├── archive/                 # Legacy bash scripts (reference only)
│   │   ├── 24_process_reads_optimized.sh  # Former production script
│   │   └── ...
│   ├── CLAUDE.md                # This file
│   └── README.md                # User-facing documentation
├── 20_mag_assembly/             # Post-expedition MAG reconstruction (Nextflow)
├── 30_archive/                  # Archived root-level scripts and docs
├── tests/                       # Pipeline tests
├── README.md                    # Project overview and quick start
├── CITATION.bib                 # References
└── LICENSE
```

## Running the Pipeline

### Nextflow (primary)

```bash
# Activate conda environment
conda activate nextflow/conda-envs/dana-tools

# Basic QC only
nextflow run nextflow/main.nf --input /path/to/data -resume

# Full pipeline
nextflow run nextflow/main.nf --input /path/to/data \
    --run_kraken --kraken_db /path/to/db \
    --run_prokka --run_sketch --run_tetra \
    -resume

# With HMM functional gene profiling
nextflow run nextflow/main.nf --input /path/to/data \
    --run_prokka --hmm_databases /path/to/CANT-HYD.hmm \
    -resume

# Watch mode for live sequencing
nextflow run nextflow/main.nf --input /path/to/runs \
    --watch --db_sync_minutes 10 \
    --run_kraken --kraken_db /path/to/db \
    --run_prokka --run_db_integration

# Quick test with bundled data
nextflow run nextflow/main.nf --input nextflow/test-data -profile test -resume

# Show all options
nextflow run nextflow/main.nf --help
```

### Launcher script

```bash
cd nextflow

# Local conda (default)
./run-realtime.sh --input /path/to/data --outdir /path/to/output \
    --run_kraken --kraken_db /path/to/krakendb \
    --run_prokka

# Docker mode
docker build -t danaseq-realtime .
./run-realtime.sh --docker --input /path/to/data --outdir /path/to/output \
    --run_kraken --kraken_db /path/to/krakendb \
    --run_prokka

# Manual docker run
docker run --user $(id -u):$(id -g) \
    -v /path/to/data:/data/input:ro \
    -v /path/to/output:/data/output \
    -v /path/to/krakendb:/kraken_db:ro \
    danaseq-realtime \
    run /pipeline/main.nf \
        --input /data/input --outdir /data/output \
        --kraken_db /kraken_db \
        --run_kraken --run_prokka -resume
```

## Pipeline Architecture

### Processing Stages

1. **VALIDATE_FASTQ** -- Gzip integrity check, BBMap repair if corrupted
2. **QC_BBDUK** -- Adapter/quality trimming
3. **QC_FILTLONG** -- Length filtering (>=1500bp), quality filtering (Q7+, top 80%)
4. **CONVERT_TO_FASTA** -- FASTQ to FASTA with header cleanup
5. **KRAKEN2_CLASSIFY** -- Taxonomic classification (maxForks=1, serialized)
6. **PROKKA_ANNOTATE** -- ORF prediction and functional annotation
7. **HMM_SEARCH** -- HMMER3 search against user-supplied HMM databases
8. **SENDSKETCH** -- Rapid taxonomic sketching
9. **TETRAMER_FREQ** -- Tetranucleotide composition for binning
10. **DB_INTEGRATION / DB_SYNC** -- Load results into DuckDB
11. **CLEANUP** -- Compress/delete source files after DB import

### Key Design Decisions

**Kraken2 serialization:** The `KRAKEN2_CLASSIFY` process uses `maxForks = 1` in its module definition. Nextflow handles this natively -- no manual semaphores needed. Kraken2 loads 50-100GB into RAM, so only one instance can run at a time.

**Resume:** Nextflow's built-in `-resume` uses task hashing. No manual checkpoint logic is needed. If a run is interrupted, rerunning with `-resume` picks up where it left off.

**Watch mode:** Uses `Channel.watchPath()` with a flat glob pattern (Java WatchService limitation -- no recursive `**`). `DB_SYNC` runs as a long-lived process with an internal sleep loop that periodically loads new results into DuckDB.

**Post-DB cleanup:** The `--cleanup` flag compresses or deletes source files after DuckDB import. Safe for watch mode since it operates per-file and checks `import_log` before deleting.

**Atomic operations:** Nextflow publishes outputs via symlinks or copy. Process scripts write to temp files then rename. Interrupted tasks leave no partial outputs in the results directory.

### Conda Environments

Three isolated environments avoid dependency conflicts:
- **dana-bbmap** -- BBMap (samtools conflicts with R)
- **dana-prokka** -- Prokka (BioPerl pins perl 5.26)
- **dana-tools** -- Filtlong, Kraken2, HMMER, R/DuckDB, Nextflow, OpenJDK

Each env has its own JDK; wrapper scripts prepend the correct `bin/` to PATH so BBMap, Prokka, and Nextflow each find their own java.

### Nextflow Config Profiles

| Profile | Use case |
|---------|----------|
| `standard` | Local execution, auto-detect resources |
| `test` | Small test files (1KB min, 500bp reads) |
| `shipboard` | Production: 32 CPUs, 256GB RAM, 100GB Kraken |
| `local_tools` | Use system-installed tools instead of conda |

## Expected Input/Output

### Input

```
input_dir/
└── fastq_pass/
    ├── barcode01/
    │   └── FLOWCELL_pass_barcode01_*.fastq.gz
    └── ...
```

`--input` must point to the directory **containing** `fastq_pass/`, not `fastq_pass/` itself.

### Output

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
- `40_kraken_db.r` -- Kraken2 classifications
- `41_krakenreport_db.r` -- Kraken2 summary reports
- `42_prokka_db.r` -- Prokka annotations
- `43_sketch_db.r` -- Sendsketch profiles
- `44_tetra_db.r` -- Tetranucleotide frequencies
- `45_stats_db.r` -- Assembly statistics
- `46_log_db.r` -- Processing logs
- `47_merge_db.r` -- Merge per-run databases
- `48_merge_all_db.r` -- Consolidate expedition data
- `49_kraken_table.r` -- Summary tables

Scripts are idempotent and track imports via `import_log`. Failures are wrapped in `|| true` so they don't crash the pipeline.

## Development Notes

### Adding a New Analysis Module

1. Create `nextflow/modules/new_analysis.nf` with a process definition
2. Add process to the import list and workflow in `main.nf`
3. Add a `--run_new_analysis` param in `nextflow.config`
4. If it needs a new tool, add to the appropriate conda env YAML in `nextflow/envs/`
5. Test with: `nextflow run nextflow/main.nf --input nextflow/test-data -profile test --run_new_analysis -resume`

### Modifying an Existing Module

Module files are in `nextflow/modules/*.nf`. Each defines one or more processes with:
- `conda` directive pointing to the env YAML
- `publishDir` for output routing
- `maxForks` for resource-constrained tools (Kraken2)
- Input/output channel declarations

### Legacy Bash Scripts

The `archive/` directory contains the original bash pipeline scripts for reference. Key files:
- `24_process_reads_optimized.sh` -- Former production script with GNU parallel, semaphores, and manual resume logic
- `22_process_reads_fast.sh` -- Rapid screening variant
- `RESUME_LOGIC.md`, `CRASH_SAFETY.md`, etc. -- Historical design documentation

These are not actively maintained. Use the Nextflow pipeline for all new work.

### MAG Assembly (20_mag_assembly/)

The MAG assembly stage is a separate Nextflow DSL2 pipeline in `20_mag_assembly/nextflow/`. It runs post-expedition (not real-time) and includes:
- Flye metagenomic co-assembly
- CoverM depth calculation (replaces jgi_summarize_bam_contig_depths)
- Consensus binning (SemiBin2 + MetaBAT2 + MaxBin2 + DAS Tool)
- Dynamic binner architecture for future extensibility

See `20_mag_assembly/CLAUDE.md` for full details.

## Common Issues

### No FASTQ files found
`--input` must point at the directory **containing** `fastq_pass/`. The pipeline gives targeted advice when files aren't found.

### Kraken2 out of memory
Nextflow serializes Kraken2 via `maxForks = 1`. If running outside Nextflow, never launch multiple Kraken2 instances. The database requires 50-100GB RAM.

### Conda environment build fails
Run `./install.sh --check` to diagnose. Ensure mamba/conda is on PATH. Build requires internet access for package downloads.

### Watch mode not picking up files
Java WatchService requires a flat glob pattern (no `**`). Adjust `--watch_glob` to match your directory depth. Default: `*/fastq_pass/barcode*/*.fastq.gz`.
