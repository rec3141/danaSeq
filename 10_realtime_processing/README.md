# Real-Time Processing Pipeline

Shipboard metagenomic analysis for Oxford Nanopore sequencing data.

---

## Overview

This directory contains scripts for real-time taxonomic classification, gene annotation, and functional profiling during active sequencing runs. Results are integrated into DuckDB for immediate SQL-queryable analysis.

**Pipeline Stages:**
```
MinKNOW Output → Validation → QC → Classification → Annotation → Database → Visualization
```

---

## Script Organization

Scripts are organized by function in numbered groups:

### 10s: Preprocessing
MinKNOW output validation and format standardization
- `10_preprocess_fastq.py` — Python preprocessing utilities
- `11_minknow_copy.sh` — Copy files from MinKNOW output directory
- `12_minknow_rename.sh` — Standardize file naming

### 20s: Read Processing
Primary analysis workflows
- `20_process_reads_basic.sh` — Single-sample processing
- `21_process_reads_barcode.sh` — Multiplexed barcode processing
- `22_process_reads_fast.sh` — Rapid screening mode
- `23_process_reads_fast2.sh` — Alternative fast processing
- **`24_process_reads_optimized.sh`** — **Recommended production pipeline**
- `25_process_reads_catchup.sh` — Batch processing of accumulated data
- `26_process_reads_test.sh` — Development testing
- `27_process_reads_fast_viz.r` — Quick visualization for fast mode

### 30s: Parsing Utilities
Data transformation and formatting
- `30_kraken_parse.awk` — Kraken2 output parsing
- `31_taxonomy_colors.awk` — Taxonomy color scheme generation

### 40s: Database Integration
DuckDB data loading and management
- `40_kraken_db.r` — Load Kraken2 classifications
- `41_krakenreport_db.r` — Load Kraken2 summary reports
- `42_prokka_db.r` — Load Prokka annotations
- `43_sketch_db.r` — Load sendsketch profiles
- `44_tetra_db.r` — Load tetranucleotide frequencies
- `45_stats_db.r` — Load assembly statistics
- `46_log_db.r` — Load processing logs
- `47_merge_db.r` — Merge run databases
- `48_merge_all_db.r` — Consolidate all expedition data
- `49_kraken_table.r` — Generate summary tables

### 50s: Visualization
Quality metrics and exploratory plots
- `50_plot_nano_qscores.r` — Nanopore quality score distributions
- `51_plot_general.r` — General plotting utilities

### 60s: Interactive Dashboards
Real-time geographic visualization
- `60_edna_mapping_viz.r` — Interactive mapping with contamination filtering

---

## Primary Workflow: 24_process_reads_optimized.sh

### Usage

**Standard Pipeline:**
```bash
./24_process_reads_optimized.sh -i <input_dir> -K -P -S
```

**With Functional Gene Profiling:**
```bash
./24_process_reads_optimized.sh -i <input_dir> -K -P --hmm /path/to/CANT-HYD.hmm
```

**Multiple HMM Databases:**
```bash
./24_process_reads_optimized.sh -i <input_dir> -P \
  --hmm /path/to/CANT-HYD.hmm,/path/to/FOAM.hmm,/path/to/NCycDB.hmm
```

### Flags

| Flag | Function | Notes |
|------|----------|-------|
| `-i` | Input directory | Required; expects `fastq_pass/` subdirectory |
| `-K` | Kraken2 classification | Memory-intensive (50-100GB RAM) |
| `-P` | Prokka annotation | ORF prediction and functional annotation |
| `-S` | Sendsketch profiling | Fast taxonomic sketching |
| `-T` | Tetranucleotide freq | Composition analysis for binning |
| `--hmm` | HMM search | Comma-delimited list of HMM database paths |
| `--force` | Force rerun | Override resume logic |
| `-d` | Debug mode | Verbose output for troubleshooting |

### Processing Steps

1. **FASTQ Validation**
   - Verify gzip integrity (`gzip -t`)
   - Repair corrupted files (BBMap `reformat.sh`)
   - Cache validated files for resume performance

2. **Quality Control**
   - Adapter removal (BBDuk)
   - Length filtering (≥1kb)
   - Quality filtering (≥Q7 mean)
   - Conversion to FASTA format

3. **Taxonomic Classification** (if `-K`)
   - Kraken2 k-mer classification
   - Serialized execution to prevent memory exhaustion
   - Output: `.tsv` (per-read), `.report` (summary)

4. **Gene Annotation** (if `-P`)
   - Prokka ORF prediction
   - Functional annotation
   - Output: GFF, FAA, FFN, TSV

5. **Functional Profiling** (if `--hmm`)
   - HMMER3 search against specified databases
   - Trusted cutoff filtering (`--cut_tc`)
   - Output: `.tsv` (hit table), `.tbl` (domain table)

6. **Optional Profiling**
   - Sendsketch taxonomic sketching (if `-S`)
   - Tetranucleotide frequency calculation (if `-T`)

### Critical: Kraken2 Memory Management

**The `-K` flag MUST only be used with `24_process_reads_optimized.sh`.**

Kraken2 loads entire database into memory (50-100GB). This script uses GNU semaphores to serialize Kraken2 calls while maintaining parallelism for other steps. Other scripts will spawn multiple instances simultaneously and crash the system.

Implementation:
```bash
sem --id kraken_db_lock kraken2 [options]
```

### Resume Logic

The pipeline implements stage-aware resume capability:

- **Basic QC**: Skips if `.fa` file exists
- **Kraken2**: Skips if `.tsv` exists (unless `--force`)
- **Prokka**: Skips if `PROKKA_*.tsv` exists
- **HMM**: Skips if `.tsv` exists for each database

Atomic operations prevent corruption on interrupt:
```bash
process > output.tmp
mv output.tmp output  # Atomic rename
```

---

## Script Comparison

| Script | Parallelism | Kraken2 Safe | Use Case |
|--------|-------------|--------------|----------|
| 20_basic | Serial | ❌ | Single sample, maximum control |
| 21_barcode | Parallel | ❌ | Standard workflow (no Kraken) |
| 22_fast | Parallel | ❌ | Rapid screening |
| **24_optimized** | **Parallel+Serial** | **✅** | **Production (Kraken OK)** |
| 25_catchup | Parallel | ❌ | Batch accumulated data |

---

## Nextflow Pipeline

A Nextflow DSL2 implementation of the same workflow with DAG-based parallelism,
native `-resume`, and resource-aware scheduling.

### Installation

```bash
cd nextflow
./install.sh          # Creates 3 conda environments under conda-envs/
./install.sh --check  # Verify all tools installed
```

Nextflow and Java are included in the `dana-tools` conda environment. Activate
it before running:

```bash
conda activate nextflow/conda-envs/dana-tools
```

### Batch Mode

Process all files then exit:

```bash
nextflow run nextflow/main.nf \
  --input /path/to/data \
  --run_kraken --kraken_db /path/to/db \
  --run_prokka --run_sketch --run_tetra \
  --run_db_integration --danadir /path/to/r_scripts \
  -resume
```

### Watch Mode

Monitor a directory for new FASTQ files during live sequencing:

```bash
nextflow run nextflow/main.nf \
  --input /path/to/runs \
  --watch --db_sync_minutes 10 \
  --run_kraken --kraken_db /path/to/db \
  --run_prokka \
  --run_db_integration --danadir /path/to/r_scripts
```

In watch mode, `DB_SYNC` runs as a long-lived process that periodically scans
output directories and loads new results into DuckDB. The R scripts are
idempotent (they track imported files via `import_log`).

### Post-DB Cleanup

On long sequencing runs, source files become redundant after loading into
DuckDB. The `--cleanup` flag compresses or deletes them after DB import:

```bash
nextflow run nextflow/main.nf --input /data/run \
  --run_db_integration --danadir /path/to/r_scripts \
  --cleanup
```

**What gets cleaned per barcode directory:**

| Directory/Files | Action |
|----------------|--------|
| `fa/*.fa` | Gzip in place (not in DuckDB, kept as compressed backup) |
| `kraken/`, `sketch/`, `tetra/` | Delete files (data lives in DuckDB) |
| `prokka/*/PROKKA_*.tsv` | Delete (already loaded into DuckDB) |
| `prokka/*/PROKKA_*.gff`, `.faa`, `.ffn` | Compress in place (gzip) |
| `hmm/`, `dana.duckdb`, `log.txt` | Kept (not cleaned) |

Cleanup operates on individual files, not whole directories, so it's safe for
watch mode where new files arrive between sync cycles. Only runs after verifying
`import_log` has entries in the barcode's DuckDB.

### Docker

A self-contained Docker image bundles all three conda environments. No host
dependencies required beyond Docker itself.

**Build:**
```bash
cd nextflow
docker build -t danaseq-realtime .
```

**Run (as current user):**
```bash
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

**Notes:**
- Always use `--user $(id -u):$(id -g)` so output files are owned by your user
- Mount the Kraken2 database directory as a read-only volume (`/kraken_db:ro`)
- The container uses `/home/dana` as a writable HOME for Nextflow metadata
- Each tool wrapper prepends its conda env to PATH, ensuring correct JDK isolation

---

## DuckDB Integration

Real-time SQL queries on growing datasets:

```r
library(DuckDB)
con <- dbConnect(duckdb(), "expedition.duckdb")

# Query cyanobacteria abundance
dbGetQuery(con, "
  SELECT sample_id, COUNT(*) as cyano_reads
  FROM kraken_results
  WHERE taxonomy LIKE '%Cyanobacteria%'
  GROUP BY sample_id
  ORDER BY cyano_reads DESC
")
```

---

## Interactive Dashboard

Launch geographic visualization:
```bash
Rscript 60_edna_mapping_viz.r
```

**Features:**
- Real-time updates as data accumulates
- GPS-tagged sample locations
- Taxonomic composition per site
- Contamination filtering (human, plant, reagents)
- Time-series visualization
- Export capabilities

---

## Input Structure

Expected directory organization:
```
input_dir/
└── fastq_pass/
    ├── barcode01/
    │   └── *.fastq.gz
    ├── barcode02/
    │   └── *.fastq.gz
    └── ...
```

---

## Output Structure

```
output/
├── FLOWCELL_ID/
│   ├── barcode01/
│   │   ├── fa/              Quality-filtered FASTA sequences
│   │   ├── fq/              Intermediate FASTQ (BBDuk output)
│   │   ├── kraken/          Taxonomic classifications
│   │   │   ├── sample.tsv       Per-read assignments
│   │   │   └── sample.report    Summary statistics
│   │   ├── prokka/          Gene annotations (per-sample dirs)
│   │   │   └── SAMPLE/
│   │   │       ├── PROKKA_*.gff
│   │   │       ├── PROKKA_*.faa    Protein sequences
│   │   │       └── PROKKA_*.tsv    Feature table
│   │   ├── hmm/             Functional gene search results
│   │   │   ├── sample.DBNAME.tsv   Hit table
│   │   │   └── sample.DBNAME.tbl   Domain table
│   │   ├── sketch/          Sendsketch profiles
│   │   ├── tetra/           Tetranucleotide frequencies
│   │   └── log.txt          Processing log
│   └── barcode02/
│       └── ...
├── expedition.duckdb        Integrated database
└── failed_files.txt         Failure diagnostics
```

---

## Troubleshooting

### No New Data Processed

**Cause:** All files already processed (resume logic working correctly)

**Solution:** Use `--force` to override, or process new data

### Kraken2 Memory Error

**Cause:** Multiple Kraken2 instances running simultaneously

**Solution:** Only use `-K` with `24_process_reads_optimized.sh`

### Prokka Not Detecting Existing Output

**Cause:** Missing `PROKKA_*.tsv` file in expected location

**Check:**
```bash
ls output/*/barcode*/prokka/*/PROKKA_*.tsv
```

### HMM Search Not Running

**Cause:** Relative path to HMM file, or gzipped HMM

**Solution:**
- Use absolute paths: `/path/to/file.hmm`
- Decompress: `gunzip file.hmm.gz`
- Verify format: `head -1 file.hmm` (should show "HMMER3")

---

## Performance Optimization

### Parallelism

Default: 32 parallel workers processing 1 file each
- Adjust with GNU parallel `-j` flag
- Monitor CPU usage: `htop`
- Monitor memory: `free -h`

### Storage

- Enable FASTQ caching for faster resume: `CACHE_FASTQ=/data/.fastq_pass`
- Use SSD for output directory
- Monitor space: `df -h`

### Resource Requirements

**Minimum (no Kraken):**
- 16 cores
- 32GB RAM
- 500GB storage

**Recommended (with Kraken):**
- 32 cores
- 128GB RAM
- 1TB storage

---

## Technical Documentation

Detailed implementation notes:
- `CLAUDE.md` — Architecture and design decisions
- `RESUME_LOGIC.md` — Stage-aware checkpoint system
- `HMM_SEARCH_GUIDE.md` — Functional gene profiling details
- `OUTPUT_FORMATS.md` — Progress line interpretation
- `DEPLOYMENT_ISSUES.md` — Portability and hardcoded paths
- `CRASH_SAFETY.md` — Atomic operations and data integrity
- `BUGFIX_*.md` — Historical issue resolutions

---

## References

See main repository `METHODS.md` and `CITATION.bib` for complete methodology and references.
