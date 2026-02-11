# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

**Important**: If user poses a question, agent responds with an answer, not a codebase change.

## Project Overview

This is the **Dana Pipeline** - a real-time metagenomic analysis system for Oxford Nanopore sequencing on oceanographic expeditions. It processes DNA reads as they stream from the sequencer, providing live taxonomic classification, gene annotation, and interactive geographic visualization.

### Expected Input Structure

The pipeline expects Oxford Nanopore directory structure:
```
input_dir/
└── fastq_pass/
    ├── barcode01/
    │   └── FLOWCELL_pass_barcode01_xxx.fastq.gz
    ├── barcode02/
    │   └── FLOWCELL_pass_barcode02_xxx.fastq.gz
    └── ...
```

### Output Structure

```
output/
├── FLOWCELL/
│   ├── barcode01/
│   │   ├── fa/           # FASTA files (final QC'd sequences)
│   │   ├── fq/           # Intermediate FASTQ (BBDuk output)
│   │   ├── kraken/       # Kraken2 classifications (.tsv, .report)
│   │   ├── prokka/       # Prokka annotations (per-sample directories)
│   │   ├── hmm/          # HMM search results (.tsv, .tbl)
│   │   ├── sketch/       # Sendsketch profiles
│   │   ├── tetra/        # Tetranucleotide frequencies
│   │   └── log.txt       # Processing log
│   └── barcode02/
│       └── ...
```

## Environment Variables & Dependencies

### Required Tools
- **BBMap** (`BBMAP=/work/apps/bbmap`)
- **Filtlong** (`FILTLONG=/work/apps/Filtlong/bin/filtlong`)
- **GNU parallel** (must be GNU version, not moreutils)

### Optional Tools (flag-dependent)
- **Kraken2** (`KRAKEN2=/usr/bin/kraken2`, `KRAKEN_DB=/path/to/db`) - for `-K`
- **Prokka** (`PROKKA_BIN=/work/apps/prokka/bin/prokka`) - for `-P`
- **HMMER** (`HMMSEARCH=/usr/bin/hmmsearch`) - for `--hmm`
- **Perl** with tetramer scripts (`APPS=/work/apps`) - for `-T`

### DuckDB Integration (Optional)
R scripts (`4X_*_db.r`) integrate results into DuckDB:
- `40_kraken_db.r` - Kraken classifications
- `42_prokka_db.r` - Prokka annotations
- `43_sketch_db.r` - Sendsketch profiles
- `44_tetra_db.r` - Tetranucleotide frequencies

**Note**: These are wrapped in `|| true` - silent failures won't crash pipeline.

### Nextflow Pipeline

A DSL2 Nextflow implementation lives in `nextflow/`. It uses three conda
environments under `nextflow/conda-envs/` (built from YAML specs in `nextflow/envs/`):
- **dana-bbmap** - BBMap (samtools conflicts with R)
- **dana-prokka** - Prokka (BioPerl pins perl 5.26)
- **dana-tools** - Filtlong, Kraken2, HMMER, R/DuckDB, Nextflow, OpenJDK

Nextflow and Java are bundled in the dana-tools env. Users activate it to run:
```bash
conda activate nextflow/conda-envs/dana-tools
nextflow run nextflow/main.nf --input /path/to/data -resume
```

Key features beyond the bash pipeline:
- **Watch mode** (`--watch`): Monitors for new FASTQ files during live sequencing;
  `DB_SYNC` periodically loads results into DuckDB
- **Post-DB cleanup** (`--cleanup`): After confirming DuckDB import via `import_log`,
  gzips fa/ in place, deletes kraken/sketch/tetra files (data in DB), deletes
  prokka TSVs, gzips prokka .gff/.faa/.ffn. Per-file operation is safe for watch mode
- Native `-resume` via Nextflow's caching (no manual checkpoint logic)
- **Docker**: Self-contained image; supports `--user` for non-root execution.
  `bin/30_kraken_parse.awk` is bundled (the `../` parent path doesn't exist in containers)

### Docker Usage

```bash
# Build
cd nextflow && docker build -t danaseq-realtime .

# Run as current user (output files owned by you, not root)
docker run --user $(id -u):$(id -g) \
    -v /path/to/data:/data/input:ro \
    -v /path/to/output:/data/output \
    -v /path/to/krakendb:/kraken_db:ro \
    danaseq-realtime \
    run /pipeline/main.nf \
        --input /data/input --outdir /data/output \
        --kraken_db /kraken_db \
        --run_kraken --run_prokka --run_sketch --run_tetra \
        --cleanup -resume
```

The container uses `/home/dana` as a writable HOME for Nextflow metadata.
Each conda env has its own JDK; wrapper scripts prepend the correct env's
`bin/` to PATH so BBMap, Prokka, and Nextflow each find their own java.

## Common Issues & Solutions

## Documentation Files

- `README.md` - User-facing guide with workflow overview
- `RESUME_LOGIC.md` - Stage-aware resume behavior with scenarios
- `HMM_SEARCH_GUIDE.md` - HMM integration details, CANT-HYD usage
- `OUTPUT_FORMATS.md` - Progress line formats, failure interpretation
- `DEPLOYMENT_ISSUES.md` - Hardcoded paths, dependency checks, portability
- `CRASH_SAFETY.md` - Atomic operations, resume guarantees
- `BUGFIX_*.md` - Historical fixes (reference for similar issues)

When helping users, **always check relevant .md files first** - they contain detailed implementation notes and known issues.
