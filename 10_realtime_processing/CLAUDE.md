# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

**Important**: If user poses a question, agent responds with an answer, not a codebase change.

## Project Overview

This is the **Dana Pipeline** - a real-time metagenomic analysis system for Oxford Nanopore sequencing on oceanographic expeditions. It processes DNA reads as they stream from the sequencer, providing live taxonomic classification, gene annotation, and interactive geographic visualization.

## Primary Script: 24_process_reads_optimized.sh

This is the production pipeline script. **All other 2X_*.sh scripts are legacy/experimental** - direct users to `24_process_reads_optimized.sh` for any processing tasks.

### Running the Pipeline

```bash
# Basic QC only (BBDuk → Filtlong → FASTA)
./24_process_reads_optimized.sh -i /path/to/nanopore/data

# With Kraken2 classification
./24_process_reads_optimized.sh -i /path/to/data -K

# Full pipeline: QC + Kraken + Prokka annotation
./24_process_reads_optimized.sh -i /path/to/data -K -P

# Add HMM search for functional genes (e.g., CANT-HYD hydrocarbon genes)
./24_process_reads_optimized.sh -i /path/to/data -P --hmm /path/to/CANT-HYD.hmm

# Force re-run of optional stages even if output exists
./24_process_reads_optimized.sh -i /path/to/data -P --hmm /path/to/HMM.hmm --force

# Multiple HMM databases (comma-delimited full paths)
./24_process_reads_optimized.sh -i data -P --hmm /path/to/CANT-HYD.hmm,/path/to/FOAM.hmm
```

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
└── failed_files.txt      # Failure log with diagnostics
```

## Architecture & Key Concepts

### 1. Parallelization Strategy (32×1 Model)

The pipeline uses **GNU parallel** with 32 workers processing 1 file each:
- Most steps (BBDuk, Filtlong, FASTA conversion): Fully parallel
- **Kraken2**: Serialized with semaphore (`sem --id kraken_db_lock`) to prevent RAM exhaustion (50-100GB per instance)
- **DuckDB writes**: Serialized with semaphore (`sem --id duckdb_lock`) to prevent lock conflicts

### 2. Resume Logic & Stage-Aware Processing

The pipeline has intelligent resume capability:

**SKIP_BASIC flag**: If `.fa` file exists, skip BBDuk → Filtlong → FASTA
**Optional stages** (Prokka, HMM, Kraken, Sketch, Tetra): Each has its own resume check

```bash
# Example: Add HMM to already-processed files
./24_process_reads_optimized.sh -i data -P              # Initial run
./24_process_reads_optimized.sh -i data -P --hmm HMM.hmm  # Only HMM runs
```

**Important**: Each optional stage checks:
```bash
if (( FORCE )) || [[ ! -s "$output_file" ]]; then
  # Run stage
fi
```

See `RESUME_LOGIC.md` for detailed scenarios.

### 3. HMM Search Integration

HMM searches run **after Prokka** on predicted proteins (`.faa` files):
- Uses `hmmsearch --cut_tc` (trusted cutoffs from HMM files)
- Supports multiple HMM databases (comma-delimited paths)
- Each database creates separate output: `sample.DBNAME.tsv`, `sample.DBNAME.tbl`
- Resume-aware: Skips if `.tsv` exists (unless `--force`)

**HMM database requirements**:
- HMMER3 format
- Uncompressed (not gzipped)
- Contains trusted cutoffs (TC scores)

See `HMM_SEARCH_GUIDE.md` for details.

### 4. FASTQ Validation & Caching

Nanopore sequencers can produce corrupted gzip files. The pipeline:
1. Validates with `gzip -t`
2. Attempts repair with BBMap `reformat.sh`
3. Caches validated files in `/data/.fastq_pass` (configurable via `CACHE_FASTQ`)

Cache enables fast resume on re-runs (validation is I/O intensive).

### 5. Crash Safety & Atomic Operations

**Atomic rename pattern** for FASTA generation:
```bash
# Write to temporary → cut headers → atomic rename
reformat.sh ... > fafile.tmp.fa
cut -f1 -d' ' fafile.tmp.fa > fafile.partial
mv fafile.partial fafile  # Atomic
```

If interrupted, `.partial` or `.tmp.fa` exist but `fafile` doesn't → Resume will retry.

### 6. Failure Handling

Failures are logged but don't halt the pipeline:
```bash
log_failure "$file" "BBDuk" "corrupted gzip - validation should have caught this"
return 0  # Continue processing other files
```

End-of-run summary shows:
- Failures by stage
- File-specific diagnostics
- Common fixes

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

See `DEPLOYMENT_ISSUES.md` for comprehensive setup guidance.

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

### Prokka Not Detecting Existing Output
Check that `PROKKA_*.tsv` exists in `output/FLOWCELL/barcode/prokka/SAMPLENAME/`.

The resume check uses:
```bash
shopt -s nullglob
local prokfiles=("$prokdir/PROKKA_"*.tsv)
if (( ${#prokfiles[@]} == 0 )); then
  # Run Prokka
fi
```

### HMM Not Running on Existing Prokka Data
Ensure HMM file paths are **absolute** (not relative) and **uncompressed**.

Check for `.faa` files:
```bash
ls output/*/barcode*/prokka/*/PROKKA_*.faa
```

### Kraken RAM Exhaustion
Only 1 Kraken runs at a time (semaphore). If multiple instances exist, check for orphaned semaphores:
```bash
sem --id kraken_db_lock --wait  # Wait for semaphore to clear
```

### ASCII Banner Not Rendering
Must use `echo -e` for ANSI escape codes:
```bash
echo -e "\033[0;36mText\033[0m"  # Correct
cat << 'EOF'                      # Wrong - heredoc doesn't interpret escapes
```

## Progress Display

The pipeline uses `parallel --eta --line-buffer`:
- **--eta**: Shows completion percentage and ETA (not filename)
- **--line-buffer**: Ensures RUN lines output cleanly

Output format:
```
ETA: 27m10s 0avg 41% 4234:5931=27m10s
RUN: sample_001 : BBDUK FILTLONG FASTA KRAKEN PROKKA HMM DONE
RUN: sample_002 : BBDUK FILTLONG FASTA DONE
```

See `OUTPUT_FORMATS.md` for complete progress line reference.

## Testing & Validation

Use debug mode to see internal state:
```bash
./24_process_reads_optimized.sh -i data -d
```

Shows:
- File existence checks
- Resume logic decisions
- Command execution details
- Semaphore acquisition

## Documentation Files

- `README.md` - User-facing guide with workflow overview
- `RESUME_LOGIC.md` - Stage-aware resume behavior with scenarios
- `HMM_SEARCH_GUIDE.md` - HMM integration details, CANT-HYD usage
- `OUTPUT_FORMATS.md` - Progress line formats, failure interpretation
- `DEPLOYMENT_ISSUES.md` - Hardcoded paths, dependency checks, portability
- `CRASH_SAFETY.md` - Atomic operations, resume guarantees
- `BUGFIX_*.md` - Historical fixes (reference for similar issues)

When helping users, **always check relevant .md files first** - they contain detailed implementation notes and known issues.
