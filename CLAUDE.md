# CLAUDE.md

Guidance for Claude Code when working in this repository.

---

## Critical Safety Rules

### File Operations

**ALWAYS use `trash-cli` instead of `rm`:**
```bash
# Correct
trash old_file.txt
trash-put deprecated_directory/

# Incorrect - permanent deletion
rm -rf directory/  # NEVER use this
```

**Pre-deletion checklist:**
1. Create manifest: `ls -laR path/ > manifest_$(date +%Y%m%d_%H%M%S).txt`
2. Confirm with user if ambiguous
3. Use `trash` for recoverability
4. Document what was removed in commit message

**Clarify ambiguous terms:**
- "Clean up" → Ask: organize files or delete files?
- "Skip" → Ask: ignore for now or permanently delete?
- "Remove" → Ask: trash (recoverable) or delete (permanent)?

### Recovery Operations

```bash
# List trashed files
trash-list

# Restore specific file
trash-restore

# Empty trash (use with extreme caution)
trash-empty
```

---

## Project Overview

**dānaSeq** is a production bioinformatics pipeline for real-time metagenomic analysis on oceanographic expeditions. It processes Oxford Nanopore sequencing data for taxonomic classification, functional gene profiling, and MAG assembly.

**Use Cases:**
- Shipboard monitoring during research cruises
- Real-time HAB (Harmful Algal Bloom) detection
- Pathogen surveillance in aquatic environments
- Publication-quality genome reconstruction

---

## Repository Architecture

```
dānaSeq/
├── 10_realtime_processing/   Real-time analysis (26 scripts)
├── 20_mag_assembly/           MAG reconstruction (26 scripts)
├── 30_archive/                Deprecated code
├── agents/                    Expert advisor system
├── METHODS.md                 Scientific methodology
├── CONTRIBUTING.md            Development guidelines
└── CITATION.bib               Reference database
```

### Numbered Naming Convention

Scripts use incremental numbering (10, 20, 30...) with intentional gaps for future insertion:
- **10s-20s:** Core processing
- **30s-40s:** Secondary analysis
- **50s-60s:** Integration pipelines
- **70s+:** Visualization

**Example:** `24_process_reads_optimized.sh` ← step 24, allows insertion of 21-23, 25-29

---

## Essential Commands

### Real-Time Processing

```bash
cd 10_realtime_processing

# Production pipeline (recommended)
./24_process_reads_optimized.sh -i /path/to/data -K -P

# With functional gene profiling
./24_process_reads_optimized.sh -i /path/to/data -P \
  --hmm /path/to/CANT-HYD.hmm,/path/to/FOAM.hmm

# Rapid screening mode
./22_process_reads_fast.sh -i /path/to/data -P -S

# Interactive dashboard
Rscript 60_edna_mapping_viz.r
```

**Flag Reference:**
- `-K` Kraken2 classification (ONLY with `24_process_reads_optimized.sh`)
- `-P` Prokka gene annotation
- `-S` Sendsketch taxonomic profiling
- `-T` Tetranucleotide frequency
- `--hmm` HMM functional gene search (comma-separated paths)
- `--force` Override resume logic

### MAG Assembly

```bash
cd 20_mag_assembly

# Complete automated pipeline
./61_map_and_bin_optimized.sh

# Visualization
Rscript 80_plot_bins.R
Rscript 82_inter_binning_analysis.r
```

### Diagnostics

```bash
./status.sh              # Dependency verification
./banner.sh              # Pipeline information
./agents.sh              # Expert advisors
```

---

## Core Technical Concepts

### 1. Critical: Kraken2 Memory Management

**Problem:** Kraken2 loads 50-100GB database into RAM. Parallel execution crashes systems.

**Solution:** Only `24_process_reads_optimized.sh` uses semaphores to serialize Kraken2:
```bash
sem --id kraken_db_lock kraken2 [options]
```

**Rule:** `-K` flag ONLY with `24_process_reads_optimized.sh`

Documentation: `10_realtime_processing/CRITICAL_KRAKEN_BUG.md`

### 2. Parallelization Model (32×1)

- **32 workers** processing **1 file each** (GNU parallel)
- **Parallel:** BBDuk, Filtlong, FASTA conversion, Prokka, HMM
- **Serial:** Kraken2 (memory), DuckDB writes (locks)

### 3. Resume Logic

All scripts implement stage-aware checkpointing:

```bash
# Check for existing output
if [[ -f "$OUTPUT" ]] && (( ! FORCE )); then
    echo "Output exists, skipping"
    return 0
fi

# Atomic write
process > "$OUTPUT.tmp"
mv "$OUTPUT.tmp" "$OUTPUT"  # Atomic rename prevents corruption
```

**Stages checked independently:**
- Basic QC → `.fa` file
- Kraken2 → `.tsv` file
- Prokka → `PROKKA_*.tsv` files
- HMM → `.tsv` per database

### 4. MAG Binning Consensus

Three complementary algorithms:
- **SemiBin2:** Deep learning (best for complex communities)
- **MetaBAT2:** TNF + coverage (fast, reliable)
- **MaxBin2:** Marker genes (good for known taxa)

**DAS Tool** integrates all three, selecting best bin for each contig.

**Never use single binner for publication.**

### 5. Data Integrity

**FASTQ Validation:**
1. `gzip -t` integrity check
2. BBMap `reformat.sh` repair if corrupted
3. Cache validated files (`/data/.fastq_pass`)

**Atomic Operations:**
- Write to `.tmp` files
- Atomic `mv` on completion
- Prevents partial writes on interrupt

---

## Software Engineering Best Practices

### Code Style

**Bash:**
```bash
#!/usr/bin/env bash
set -euo pipefail  # Exit on error, undefined vars, pipe failures

# Variables: UPPER_CASE
INPUT_DIR="/path/to/input"
OUTPUT_FILE="${INPUT_DIR}/results.txt"

# Functions: lowercase_with_underscores
process_sample() {
    local sample=$1
    local output=$2

    # Check inputs
    [[ -f "$sample" ]] || { echo "Error: $sample not found" >&2; return 1; }

    # Process with error handling
    if ! tool --input "$sample" > "$output" 2>&1 | tee -a "$LOG"; then
        log_error "Processing failed for $sample"
        return 1
    fi

    return 0
}

# Prefer [[ ]] over [ ]
[[ -f "$file" ]] && process_file "$file"

# Quote all variables
cp "$source_file" "$destination_dir/"
```

**R:**
```r
library(tidyverse)

# snake_case for variables and functions
process_coverage_data <- function(input_file, min_coverage = 10) {
  # Explicit handling of missing data
  data <- read_tsv(input_file, col_types = cols()) %>%
    filter(!is.na(coverage), coverage >= min_coverage)

  return(data)
}
```

**Python:**
```python
from pathlib import Path
from typing import List, Optional

def validate_fastq(
    input_path: Path,
    min_length: int = 1000,
    min_quality: float = 7.0
) -> bool:
    """
    Validate FASTQ file integrity and quality.

    Args:
        input_path: Path to FASTQ file
        min_length: Minimum read length in bp
        min_quality: Minimum mean quality score

    Returns:
        True if valid, False otherwise
    """
    if not input_path.exists():
        raise FileNotFoundError(f"FASTQ not found: {input_path}")

    # Implementation
    pass
```

### Error Handling

**Comprehensive error reporting:**
```bash
log_error() {
    local message=$1
    local file=$2
    local stage=$3

    echo "ERROR [$(date +'%Y-%m-%d %H:%M:%S')] $stage failed for $file: $message" >&2
    echo "$file,$stage,$message" >> failed_files.txt
}

# Usage
if ! run_kraken "$sample"; then
    log_error "Kraken classification failed" "$sample" "kraken2"
    return 1
fi
```

**User-friendly messages:**
```bash
# Bad
echo "Error 127"

# Good
echo "ERROR: Kraken2 not found. Install with: conda install -c bioconda kraken2"
echo "See DEPLOYMENT_ISSUES.md for full setup instructions"
```

### Testing Strategy

**Before committing:**
1. Test with small dataset (1K reads)
2. Verify resume logic works
3. Test failure recovery
4. Check resource usage
5. Validate output format

**Test datasets:**
```
/data/test/
├── small/    1K reads   (30s runtime)
├── medium/   10K reads  (5min runtime)
└── large/    100K reads (1hr runtime)
```

**Testing commands:**
```bash
# Dry run to see commands
./24_process_reads_optimized.sh -i test_data --dry-run

# Debug mode
./24_process_reads_optimized.sh -i test_data -d

# Check outputs
ls -lh output/
checksum validation here
```

### Performance Optimization

**Profiling:**
```bash
# CPU usage
time ./script.sh -i data
htop  # Monitor during execution

# Memory usage
/usr/bin/time -v ./script.sh -i data

# Disk I/O
iotop  # Monitor during execution
```

**Optimization checklist:**
- [ ] Minimize disk I/O (use tmpfs for intermediate files)
- [ ] Avoid redundant processing (implement resume logic)
- [ ] Parallelize independent operations
- [ ] Serialize resource-intensive operations (Kraken2, DuckDB)
- [ ] Use streaming where possible (avoid loading entire files)

### Security Considerations

**Never commit:**
- API keys (use environment variables)
- Passwords or tokens
- Hardcoded credentials
- Private data or results

**Check before committing:**
```bash
# Scan for secrets
git diff --cached | grep -iE '(password|api[_-]?key|token|secret)'

# Remove if found
git reset HEAD file_with_secret.txt
# Edit to remove secret
git add file_with_secret.txt
```

**Secure practices:**
```bash
# Use environment variables
API_KEY="${OPENAI_API_KEY:-}"
[[ -z "$API_KEY" ]] && { echo "Set OPENAI_API_KEY" >&2; exit 1; }

# Never log secrets
echo "Processing with API key: ${API_KEY:0:8}..." # Show only prefix
```

---

## Code Review Checklist

Before approving changes:

**Functionality:**
- [ ] Implements requested feature correctly
- [ ] Handles edge cases (empty input, missing files, corrupted data)
- [ ] Resume logic works correctly
- [ ] Error messages are informative

**Code Quality:**
- [ ] Follows project naming conventions
- [ ] Consistent indentation and style
- [ ] No hardcoded paths (or documented in DEPLOYMENT_ISSUES.md)
- [ ] Proper error handling with useful messages
- [ ] Comments explain complex logic

**Testing:**
- [ ] Tested with small dataset
- [ ] Tested resume from interruption
- [ ] Resource usage is reasonable
- [ ] No memory leaks or excessive disk usage

**Documentation:**
- [ ] README.md updated if user-facing change
- [ ] METHODS.md updated if methodology changed
- [ ] CLAUDE.md updated if architecture changed
- [ ] Inline comments for complex sections
- [ ] Script header complete and accurate

**Safety:**
- [ ] Uses `trash` instead of `rm`
- [ ] No secrets committed
- [ ] Atomic operations for data integrity
- [ ] Appropriate use of `set -euo pipefail`

---

## Input/Output Specifications

### Real-Time Processing Input

```
input_dir/
└── fastq_pass/
    ├── barcode01/
    │   └── FLOWCELL_pass_barcode01_*.fastq.gz
    └── barcode02/
        └── FLOWCELL_pass_barcode02_*.fastq.gz
```

### Real-Time Processing Output

```
output/
├── FLOWCELL_ID/
│   ├── barcode01/
│   │   ├── fa/              FASTA (final QC'd sequences)
│   │   ├── fq/              FASTQ (BBDuk intermediate)
│   │   ├── kraken/          *.tsv, *.report
│   │   ├── prokka/          SAMPLE/PROKKA_*
│   │   ├── hmm/             *.DBNAME.tsv, *.DBNAME.tbl
│   │   └── log.txt
│   └── ...
├── expedition.duckdb        Integrated database
└── failed_files.txt         Failure log
```

### MAG Assembly Output

```
mag_output/
├── assembly/assembly.fasta
├── bins/das_tool_consensus/    ← Use these
│   ├── MAG_00001.fa
│   └── ...
├── polished/MAG_*.polished.fa
├── checkm2/quality_report.tsv
└── taxonomy/assignments.txt
```

---

## Common Development Patterns

### Adding New Analysis Step

```bash
# 1. Choose number (e.g., 35 between 30 and 40)
# 2. Create script
cat > 35_new_analysis.sh << 'EOF'
#!/usr/bin/env bash
set -euo pipefail

# Script header (see existing scripts for template)

# Parse arguments
# Implement resume logic
# Run analysis
# Log results
EOF

chmod +x 35_new_analysis.sh

# 3. Test
./35_new_analysis.sh -i test_data/

# 4. Update documentation
vim README.md
```

### Implementing Resume Logic

```bash
run_analysis() {
    local input=$1
    local output=$2

    # Resume check
    if [[ -f "$output" ]] && (( ! FORCE )); then
        log_info "Output exists: $output (skipping, use --force to override)"
        return 0
    fi

    # Process to temporary file
    if ! tool --input "$input" --output "$output.tmp" 2>&1 | tee -a "$LOG"; then
        log_error "Tool failed" "$input" "tool"
        rm -f "$output.tmp"  # Clean up on failure
        return 1
    fi

    # Atomic rename
    mv "$output.tmp" "$output"
    log_info "Completed: $output"
    return 0
}
```

### Database Integration Pattern

```r
library(DuckDB)
library(tidyverse)

load_results_to_db <- function(results_dir, db_path) {
  # Connect
  con <- dbConnect(duckdb(), db_path)
  on.exit(dbDisconnect(con, shutdown = TRUE))

  # Create table if needed
  if (!dbExistsTable(con, "results")) {
    dbExecute(con, "
      CREATE TABLE results (
        sample_id VARCHAR,
        metric DOUBLE,
        timestamp TIMESTAMP DEFAULT CURRENT_TIMESTAMP
      )
    ")
  }

  # Load data
  files <- list.files(results_dir, pattern = "\\.tsv$", full.names = TRUE)
  for (f in files) {
    data <- read_tsv(f, col_types = cols()) %>%
      mutate(sample_id = tools::file_path_sans_ext(basename(f)))

    dbWriteTable(con, "results", data, append = TRUE)
  }

  message("Loaded ", length(files), " files to database")
}
```

---

## Troubleshooting Guide

### Script Not Finding Dependencies

**Symptom:** `command not found: kraken2`

**Solution:**
1. Check paths: `./status.sh`
2. Verify installation: `which kraken2`
3. Update script variables: Edit `KRAKEN2` path in script header
4. See `DEPLOYMENT_ISSUES.md`

### Resume Logic Not Working

**Symptom:** Script re-processes existing files

**Solution:**
1. Check output file exists: `ls -lh output/sample.fa`
2. Check file is non-empty: `[[ -s output/sample.fa ]] && echo "exists"`
3. Use `--force` to intentionally override
4. Check for `.tmp` or `.partial` files indicating interrupted run

### Memory Exhaustion

**Symptom:** System freeze, OOM killer, swap thrashing

**Solution:**
1. Kraken2: Only use `-K` with `24_process_reads_optimized.sh`
2. Check parallel workers: Reduce `-j` parameter in script
3. Monitor: `watch -n1 free -h`
4. Increase swap: `sudo fallocate -l 128G /swapfile` (temporary fix)

### DuckDB Lock Errors

**Symptom:** `database is locked`

**Solution:**
1. Only one write process at a time
2. Check for orphaned connections: `lsof expedition.duckdb`
3. Scripts use semaphore: `sem --id duckdb_lock`
4. If corrupted: Rebuild from source files

---

## Reference Documentation

**Repository:**
- `README.md` — Project overview
- `METHODS.md` — Scientific methodology
- `CONTRIBUTING.md` — Development guidelines
- `CITATION.bib` — Reference database

**Real-Time Processing:**
- `10_realtime_processing/README.md` — Workflow guide
- `10_realtime_processing/CLAUDE.md` — Architecture details
- `RESUME_LOGIC.md` — Checkpoint system
- `HMM_SEARCH_GUIDE.md` — Functional profiling
- `CRASH_SAFETY.md` — Data integrity
- `CRITICAL_KRAKEN_BUG.md` — Memory management

**MAG Assembly:**
- `20_mag_assembly/README.md` — Assembly methodology

**Always consult relevant documentation before making changes.**

---

## Quick Decision Tree

**User asks to "clean up" →** Ask: "Organize files or delete files?"

**User asks to "skip" →** Ask: "Ignore for now or permanently delete?"

**Adding new script →** Use next available number with gap (25, 35, 45...)

**Script fails →** Check: dependencies (`./status.sh`), logs (`tail -f log.txt`), failed files (`cat failed_files.txt`)

**Kraken2 OOM →** Verify: Only `24_process_reads_optimized.sh` uses `-K`

**Need to delete files →** Use: `trash` (not `rm`)

**Committing changes →** Check: No secrets, documentation updated, tested on small dataset

---

**Repository:** https://github.com/rec3141/danaSeq
**License:** MIT (see LICENSE file)
**Contact:** rec3141@gmail.com
