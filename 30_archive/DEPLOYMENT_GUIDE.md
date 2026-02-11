# dānaSeq Deployment Guide

**Version:** 2.0 (Post-Security Audit)
**Last Updated:** 2026-02-06
**Target Audience:** Human operators + LLM assistants
**Deployment Type:** Shipboard production environment

---

## Table of Contents

1. [Quick Start](#quick-start)
2. [Pre-Deployment Checklist](#pre-deployment-checklist)
3. [System Requirements](#system-requirements)
4. [Installation](#installation)
5. [Configuration](#configuration)
6. [Verification](#verification)
7. [First Run](#first-run)
8. [Production Operations](#production-operations)
9. [Monitoring](#monitoring)
10. [Troubleshooting](#troubleshooting)
11. [Maintenance](#maintenance)
12. [Rollback Procedures](#rollback-procedures)
13. [For LLM Assistants](#for-llm-assistants)

---

## Quick Start

**For experienced operators:**

```bash
# 1. Clone and verify
git clone https://github.com/rec3141/danaSeq.git
cd danaSeq
git log --oneline -4  # Should show security fix commits

# 2. Check dependencies
./status.sh

# 3. Test run (small dataset)
cd 10_realtime_processing
./24_process_reads_optimized.sh -i /path/to/test/data -P

# 4. Verify output
ls -lh out_*/*/barcode*/fa/*.fa  # Should see FASTA files
```

**Estimated time:** 30 minutes (small test dataset)

---

## Pre-Deployment Checklist

### For Human Operators

Use this checklist before deploying to production:

```
Hardware & Environment:
[ ] Server has 32+ CPU cores
[ ] Server has 128+ GB RAM (for Kraken2)
[ ] Storage has 1+ TB available space
[ ] Network connectivity verified (for initial setup)
[ ] UPS/power backup configured
[ ] Temperature monitoring in place

Software Dependencies:
[ ] Linux OS (Ubuntu 20.04+ or CentOS 7+)
[ ] Bash 4.0+
[ ] Python 3.9+
[ ] R 4.2+
[ ] All bioinformatics tools installed (see status.sh)
[ ] GNU parallel installed (not moreutils parallel)
[ ] trash-cli installed (for safe file deletion)

Security & Access:
[ ] User has sudo access for initial setup
[ ] Data directories have correct permissions
[ ] Backup system configured
[ ] Log rotation configured
[ ] Monitoring alerts configured

Data Preparation:
[ ] Test dataset available (~10K reads)
[ ] Reference databases downloaded
[ ] HMM databases prepared (if using)
[ ] Metadata files created

Documentation:
[ ] Team trained on basic operation
[ ] Emergency contacts documented
[ ] Backup procedures documented
```

### For LLM Assistants

**Critical verification points:**

1. **Security status:** Verify repository is post-audit (commit 3a8964d or later)
2. **Dependencies:** Parse `./status.sh` output for missing tools
3. **Permissions:** Check write access to `/data/` directory
4. **Resources:** Verify RAM ≥128GB if using Kraken2 (`-K` flag)
5. **Disk space:** Ensure ≥1TB available in output directory

**Automated checks:**
```bash
# Run as preliminary diagnostic
./status.sh 2>&1 | tee deployment-check.log
df -h /data
free -h
```

---

## System Requirements

### Minimum (No Kraken2)

```yaml
CPU: 16 cores
RAM: 32 GB
Storage: 500 GB
OS: Linux (kernel 3.10+)
Network: Not required (fully offline capable)
```

### Recommended (With Kraken2)

```yaml
CPU: 32+ cores (64 optimal)
RAM: 128 GB (256 GB optimal)
Storage: 1-2 TB SSD/NVMe
OS: Ubuntu 20.04+ / CentOS 7+
Network: 1 Gbps (for initial DB download)
```

### Storage Breakdown

```
Component                Size    Purpose
─────────────────────────────────────────────────────
Kraken2 database         50-100 GB   Taxonomic classification
Reference databases      20-50 GB    GTDB, Kaiju, etc.
HMM databases           1-5 GB      Functional gene profiling
Input data (per run)    10-100 GB   Raw FASTQ files
Output data (per run)   20-200 GB   Processed results + MAGs
Temporary files         50-100 GB   Intermediate processing
Working space           100+ GB     Buffer for operations
─────────────────────────────────────────────────────
TOTAL RECOMMENDED       1-2 TB
```

---

## Installation

### Step 1: Clone Repository

```bash
# Choose installation directory (must be on large partition)
cd /data
git clone https://github.com/rec3141/danaSeq.git
cd danaSeq

# Verify security fixes are present
git log --oneline --grep="security\|CRITICAL\|HIGH\|MEDIUM\|LOW" | head -4
# Expected: 4 security fix commits from 2026-02-06
```

**Verification:**
```bash
# Check for critical security documentation
ls -1 SECURITY_*.md
# Expected output:
# SECURITY_AUDIT.md
# SECURITY_PATCH_2026-02-06.md
```

### Step 2: Install System Dependencies

#### Ubuntu/Debian

```bash
# Update package manager
sudo apt update

# Core utilities
sudo apt install -y \
  build-essential \
  git \
  parallel \
  python3 python3-pip \
  r-base \
  gzip pigz \
  trash-cli

# Verify GNU parallel (not moreutils)
parallel --version | head -1
# Expected: GNU parallel 20XXXXXX
```

#### CentOS/RHEL

```bash
# Enable EPEL repository
sudo yum install -y epel-release

# Core utilities
sudo yum install -y \
  gcc gcc-c++ make \
  git \
  parallel \
  python3 python3-pip \
  R \
  gzip pigz \
  trash-cli

# Verify GNU parallel
parallel --version | head -1
```

### Step 3: Install Bioinformatics Tools

**Option A: Use Conda (Recommended)**

```bash
# Install Miniconda
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh -b -p $HOME/miniconda3
source $HOME/miniconda3/bin/activate

# Create danaSeq environment
conda create -n danaseq -c bioconda -c conda-forge \
  bbmap filtlong kraken2 prokka hmmer \
  flye minimap2 samtools \
  semibin2 metabat2 maxbin2 das_tool \
  racon medaka checkm2 \
  kaiju gtdbtk

# Activate environment
conda activate danaseq
```

**Option B: Manual Installation**

See `INSTALLATION.md` for tool-by-tool instructions.

### Step 4: Download Reference Databases

**Critical:** Databases are large (50-100GB). Download before expedition.

```bash
# Create database directory
mkdir -p /data/databases
cd /data/databases

# Kraken2 database (choose one)
# Option 1: Standard Plus PFP (8GB, fast)
kraken2-build --download-library bacteria --db kraken_standard_8gb
kraken2-build --download-library viral --db kraken_standard_8gb
kraken2-build --build --db kraken_standard_8gb

# Option 2: Pre-built database (recommended)
wget https://genome-idx.s3.amazonaws.com/kraken/k2_pluspfp_08gb_20230605.tar.gz
tar -xzf k2_pluspfp_08gb_20230605.tar.gz -C kraken_db/

# GTDB-Tk database
wget https://data.gtdb.ecogenomic.org/releases/latest/auxillary_files/gtdbtk_data.tar.gz
tar -xzf gtdbtk_data.tar.gz -C gtdbtk_db/

# Kaiju database
kaiju-makedb -s progenomes -t /data/databases/kaiju_db
```

**Verification:**
```bash
# Check database sizes
du -sh /data/databases/*
# Expected:
# 8G    kraken_db
# 70G   gtdbtk_db
# 50G   kaiju_db
```

### Step 5: Install R Packages

```bash
# Launch R
R

# Install required packages
install.packages(c("dplyr", "tidyr", "ggplot2", "readr"))
install.packages("DBI")
install.packages("duckdb")
install.packages("Rtsne")
install.packages("umap")

# Verify installation
library(duckdb)
library(dplyr)
# Should load without errors
```

---

## Configuration

### Step 1: Set Environment Variables

Create `/etc/profile.d/danaseq.sh`:

```bash
# dānaSeq Environment Configuration
export DANASEQ_ROOT="/data/danaSeq"
export THREADS=32
export BBMAP="/path/to/bbmap"
export DANADIR="/data/danaSeq/10_realtime_processing"
export PROKKA_BIN="/path/to/prokka/bin/prokka"
export FILTLONG="/path/to/filtlong/bin/filtlong"
export KRAKEN2="/usr/bin/kraken2"
export KRAKEN_DB="/data/databases/kraken_db"
export HMMSEARCH="/usr/bin/hmmsearch"
export APPS="/data/apps"

# Add to PATH
export PATH="${DANASEQ_ROOT}/10_realtime_processing:${PATH}"
export PATH="${DANASEQ_ROOT}/20_mag_assembly:${PATH}"
```

**Load configuration:**
```bash
source /etc/profile.d/danaseq.sh
```

### Step 2: Configure Logging

Create `/etc/logrotate.d/danaseq`:

```
/data/danaSeq/logs/*.log {
    daily
    rotate 30
    compress
    delaycompress
    notifempty
    create 0640 $USER $GROUP
    sharedscripts
    postrotate
        # Restart monitoring if needed
    endscript
}
```

### Step 3: Set Up Monitoring

**Basic monitoring script** (`/usr/local/bin/danaseq-monitor.sh`):

```bash
#!/usr/bin/env bash
# dānaSeq Monitoring Script

ALERT_EMAIL="your-email@example.com"
ALERT_RAM_THRESHOLD=90  # Percent
ALERT_DISK_THRESHOLD=85  # Percent

# Check RAM usage
RAM_USAGE=$(free | grep Mem | awk '{print int($3/$2 * 100)}')
if (( RAM_USAGE > ALERT_RAM_THRESHOLD )); then
    echo "HIGH RAM USAGE: ${RAM_USAGE}%" | mail -s "dānaSeq Alert" "$ALERT_EMAIL"
fi

# Check disk usage
DISK_USAGE=$(df /data | tail -1 | awk '{print int($3/$2 * 100)}')
if (( DISK_USAGE > ALERT_DISK_THRESHOLD )); then
    echo "HIGH DISK USAGE: ${DISK_USAGE}%" | mail -s "dānaSeq Alert" "$ALERT_EMAIL"
fi

# Check for error messages in recent logs
ERROR_COUNT=$(find /data/danaSeq -name "*.log" -mtime -1 -exec grep -c "^\[ERROR\]" {} \; | awk '{s+=$1} END {print s}')
if (( ERROR_COUNT > 10 )); then
    echo "HIGH ERROR COUNT: ${ERROR_COUNT} errors in last 24h" | mail -s "dānaSeq Alert" "$ALERT_EMAIL"
fi
```

**Add to crontab:**
```bash
# Run monitoring every hour
0 * * * * /usr/local/bin/danaseq-monitor.sh
```

---

## Verification

### Step 1: Dependency Check

```bash
cd /data/danaSeq
./status.sh > deployment-verification.txt 2>&1

# Review output
cat deployment-verification.txt
```

**Expected output:**
```
✓ BBMap found
✓ Filtlong found
✓ Kraken2 found
✓ Prokka found
✓ HMMER found
✓ Flye found
✓ minimap2 found
✓ samtools found
...all tools present
```

**If any tools are missing:**
```bash
# Check which tools are required for your workflow
# Minimum (basic QC): BBMap, Filtlong
# Standard: + Prokka, Kraken2
# Full: + all MAG assembly tools
```

### Step 2: Security Verification

```bash
# Verify security patches are present
cd /data/danaSeq

# Check for atomic lock implementation (CRITICAL fix)
grep -A5 "mkdir.*lockfile" 10_realtime_processing/24_process_reads_optimized.sh
# Should show atomic lock acquisition

# Check for bash -c (not eval)
grep "eval" 10_realtime_processing/24_process_reads_optimized.sh
# Should ONLY appear in comments, not in code

# Check for signal handlers
grep "trap.*EXIT" 10_realtime_processing/24_process_reads_optimized.sh
# Should show cleanup handlers

# Check for error message standardization
grep '\[ERR\]' 10_realtime_processing/24_process_reads_optimized.sh
# Should return NO matches (all should be [ERROR])
```

### Step 3: Test Run

**Small test dataset (recommended first run):**

```bash
cd /data/danaSeq/10_realtime_processing

# Create test data directory
mkdir -p /data/test/fastq_pass/barcode01
cp /path/to/sample/*.fastq.gz /data/test/fastq_pass/barcode01/

# Run basic QC (no Kraken, no Prokka)
./24_process_reads_optimized.sh -i /data/test -v

# Expected runtime: 5-10 minutes for 10K reads
```

**Verify test output:**
```bash
# Check directory structure
tree -L 3 out_test_*/

# Expected structure:
# out_test_YYYYMMDD_HHMMSS/
# ├── FLOWCELL/
# │   └── barcode01/
# │       ├── fa/          # FASTA files (critical)
# │       ├── fq/          # Intermediate FASTQ
# │       └── log.txt      # Processing log
# └── failed_files.txt     # Should be empty or not exist

# Verify FASTA files were created
ls -lh out_test_*/*/barcode01/fa/*.fa
# Should show non-zero file sizes

# Check log for errors
grep '^\[ERROR\]' out_test_*/*/barcode01/log.txt
# Should show NO errors (or only expected errors like missing tools)

# Check log for completion
grep '^\[INFO\]' out_test_*/*/barcode01/log.txt | tail -5
# Should show completion messages
```

### Step 4: Full Pipeline Test

**With all features enabled:**

```bash
# Run with Kraken2, Prokka, and HMM search
./24_process_reads_optimized.sh \
  -i /data/test \
  -K \
  -P \
  --hmm /data/databases/CANT-HYD.hmm \
  -v

# Expected runtime: 30-60 minutes for 10K reads
```

**Verify full output:**
```bash
OUTDIR="out_test_$(date +%Y%m%d)_"*

# Check all output directories exist
ls -ld "$OUTDIR"/*/barcode01/{fa,kraken,prokka,hmm}
# All should exist

# Check Kraken output
head -5 "$OUTDIR"/*/barcode01/kraken/*.tsv
# Should show taxonomic assignments

# Check Prokka output
ls "$OUTDIR"/*/barcode01/prokka/*/PROKKA_*.tsv
# Should exist

# Check HMM output
head -5 "$OUTDIR"/*/barcode01/hmm/*.tsv
# Should show gene matches (if present in data)
```

---

## First Run

### Production Data Structure

**Expected MinKNOW output:**
```
/data/minknow/
└── RUN_2025XXXX/
    └── fastq_pass/
        ├── barcode01/
        │   └── FLOWCELL_pass_barcode01_*.fastq.gz
        ├── barcode02/
        │   └── FLOWCELL_pass_barcode02_*.fastq.gz
        └── ...
```

### Command Construction

**Basic workflow (QC only):**
```bash
./24_process_reads_optimized.sh \
  -i /data/minknow/RUN_2025XXXX \
  -o /data/results/RUN_2025XXXX_$(date +%Y%m%d_%H%M%S)
```

**Standard workflow (QC + classification + annotation):**
```bash
./24_process_reads_optimized.sh \
  -i /data/minknow/RUN_2025XXXX \
  -o /data/results/RUN_2025XXXX_$(date +%Y%m%d_%H%M%S) \
  -K \
  -P \
  -v
```

**Full workflow (with functional profiling):**
```bash
./24_process_reads_optimized.sh \
  -i /data/minknow/RUN_2025XXXX \
  -o /data/results/RUN_2025XXXX_$(date +%Y%m%d_%H%M%S) \
  -K \
  -P \
  --hmm /data/databases/CANT-HYD.hmm,/data/databases/FOAM.hmm \
  -v
```

**Resume interrupted run:**
```bash
# Simply re-run the same command
# Resume logic automatically skips completed files
./24_process_reads_optimized.sh \
  -i /data/minknow/RUN_2025XXXX \
  -o /data/results/RUN_2025XXXX_TIMESTAMP \
  -K -P -v
```

**Force re-run (override resume):**
```bash
./24_process_reads_optimized.sh \
  -i /data/minknow/RUN_2025XXXX \
  -o /data/results/RUN_2025XXXX_TIMESTAMP \
  -K -P \
  --force
```

### Real-Time Processing

**Monitor progress:**
```bash
# In separate terminal
watch -n 30 'find /data/results/LATEST -name "*.fa" | wc -l'

# Or watch log
tail -f /data/results/LATEST/*/barcode01/log.txt
```

**Check for errors:**
```bash
# Real-time error monitoring
tail -f /data/results/LATEST/failed_files.txt

# Or check all logs
find /data/results/LATEST -name "log.txt" -exec grep '^\[ERROR\]' {} +
```

---

## Production Operations

### Daily Checklist

**Morning startup:**
```
[ ] Check system resources: df -h /data && free -h
[ ] Review overnight logs: grep ERROR /data/results/*/log.txt
[ ] Verify MinKNOW is running: systemctl status minknow
[ ] Check available disk space: > 500GB free
[ ] Start processing if data accumulated overnight
```

**During run:**
```
[ ] Monitor progress: watch output directory
[ ] Check for errors: tail -f failed_files.txt
[ ] Monitor resources: htop
[ ] Verify data integrity: check file sizes increasing
```

**End of day:**
```
[ ] Verify all data processed: check file counts
[ ] Review error log: failed_files.txt
[ ] Backup results: rsync to backup server
[ ] Check disk space: df -h
[ ] Update run log: record samples processed
```

### Resource Management

**Monitor during operation:**

```bash
# CPU usage
htop

# RAM usage (especially important for Kraken2)
free -h
watch -n 5 free -h

# Disk usage
df -h /data
du -sh /data/results/*

# I/O performance
iostat -x 5

# Process status
ps aux | grep 24_process_reads
```

**Safe stop procedures:**

```bash
# Graceful stop (Ctrl+C in terminal)
# - Signal handlers will clean up
# - Temp files removed
# - Locks released
# - Can resume later

# If process is backgrounded
kill -TERM <PID>  # Graceful termination
# Do NOT use: kill -9  (skips cleanup)
```

### Batch Processing

**Process accumulated data:**

```bash
# Process all runs from past week
for run in /data/minknow/RUN_202501*/; do
  echo "Processing: $run"
  ./24_process_reads_optimized.sh \
    -i "$run" \
    -o "/data/results/$(basename "$run")_$(date +%Y%m%d_%H%M%S)" \
    -K -P -v
done
```

**Time-limited batches:**

```bash
# Process for 2 hours, then stop gracefully
./24_process_reads_optimized.sh \
  -i /data/minknow/RUN_2025XXXX \
  -K -P \
  --max-duration 7200
```

---

## Monitoring

### Key Metrics

**Performance metrics:**
```bash
# Files processed per hour
OUTDIR="/data/results/LATEST"
FILES_PROCESSED=$(find "$OUTDIR" -name "*.fa" | wc -l)
RUNTIME_HOURS=$(echo "scale=2; $(date +%s) - $(stat -c%Y "$OUTDIR")" | bc)
echo "Rate: $(echo "$FILES_PROCESSED / $RUNTIME_HOURS" | bc) files/hour"

# Failure rate
TOTAL_FILES=$(find /data/minknow/LATEST/fastq_pass -name "*.fastq.gz" | wc -l)
FAILED_FILES=$(wc -l < "$OUTDIR/failed_files.txt")
echo "Failure rate: $(echo "scale=2; $FAILED_FILES / $TOTAL_FILES * 100" | bc)%"
```

**Quality metrics:**
```bash
# Average read count per file
find "$OUTDIR" -name "*.fa" -exec grep -c '^>' {} \; | \
  awk '{s+=$1; n++} END {print "Avg reads:", s/n}'

# File size distribution
find "$OUTDIR" -name "*.fa" -exec stat -c%s {} \; | \
  awk '{sum+=$1; n++} END {print "Avg size:", sum/n/1024/1024, "MB"}'
```

### Log Analysis

**Parse structured logs:**

```bash
# Count errors by type
grep '^\[ERROR\]' "$OUTDIR"/*/barcode*/log.txt | \
  cut -d: -f2 | sort | uniq -c | sort -rn

# Count warnings
grep '^\[WARNING\]' "$OUTDIR"/*/barcode*/log.txt | wc -l

# Extract processing times
grep 'RUN:' "$OUTDIR"/*/barcode*/log.txt | \
  grep -o 'DONE' | wc -l
```

**Generate daily report:**

```bash
#!/usr/bin/env bash
# daily-report.sh

OUTDIR="/data/results/$(date +%Y%m%d)_"*

cat << EOF
dānaSeq Daily Report - $(date +%Y-%m-%d)
========================================

Files Processed: $(find "$OUTDIR" -name "*.fa" 2>/dev/null | wc -l)
Failed Files: $(wc -l < "$OUTDIR/failed_files.txt" 2>/dev/null || echo 0)
Errors: $(find "$OUTDIR" -name "log.txt" -exec grep -c '^\[ERROR\]' {} + | awk '{s+=$1} END {print s}')
Warnings: $(find "$OUTDIR" -name "log.txt" -exec grep -c '^\[WARNING\]' {} + | awk '{s+=$1} END {print s}')

Disk Usage: $(df -h /data | tail -1 | awk '{print $3 "/" $2 " (" $5 " used)"}')
RAM Usage: $(free -h | grep Mem | awk '{print $3 "/" $2}')

Top Errors:
$(grep '^\[ERROR\]' "$OUTDIR"/*/barcode*/log.txt 2>/dev/null | \
  cut -d: -f2 | sort | uniq -c | sort -rn | head -5)
EOF
```

---

## Troubleshooting

### Common Issues

#### 1. "Cannot create output directory"

**Symptoms:**
```
[ERROR] Cannot create output directory: /data/results/run_XXXXX
```

**Diagnosis:**
```bash
# Check permissions
ls -ld /data/results
# Should show write permissions for your user

# Check disk space
df -h /data
# Should have >500GB free
```

**Solution:**
```bash
# Fix permissions
sudo chown -R $USER:$USER /data/results

# Or create directory manually
mkdir -p /data/results/run_XXXXX
```

#### 2. "Kraken2 failed" or Out of Memory

**Symptoms:**
```
[ERROR] Kraken2 classification failed
# Or system becomes unresponsive during processing
```

**Diagnosis:**
```bash
# Check RAM
free -h
# Total should be >128GB for Kraken2

# Check if multiple Kraken2 instances running
ps aux | grep kraken2
# Should only see 1 instance (semaphore serialization)
```

**Solution:**
```bash
# If insufficient RAM, run without Kraken2
./24_process_reads_optimized.sh -i DATA -P  # No -K flag

# If semaphore issue, clear locks
sem --id kraken_db_lock --wait
```

#### 3. "No FASTQ files found"

**Symptoms:**
```
[WARNING] No FASTQ files found matching criteria
```

**Diagnosis:**
```bash
# Check directory structure
find /data/input -name "*.fastq.gz"

# Expected structure
ls -R /data/input/fastq_pass/
```

**Solution:**
```bash
# Verify input path contains fastq_pass/ subdirectory
# Correct: -i /data/minknow/RUN_XXXX
# Wrong: -i /data/minknow/RUN_XXXX/fastq_pass
```

#### 4. "Validation failed" errors

**Symptoms:**
```
[ERROR] Validation failed: /path/to/file.fastq.gz
```

**Diagnosis:**
```bash
# Test gzip integrity
gzip -t /path/to/file.fastq.gz
# Will show "unexpected end of file" if corrupt

# Check file size
ls -lh /path/to/file.fastq.gz
# Very small files (<1KB) are likely incomplete
```

**Solution:**
```bash
# File is corrupt - re-copy from sequencer
# Or if actively sequencing, wait for file to complete

# Check for incomplete writes
lsof /path/to/file.fastq.gz
# If file is open, wait for write to complete
```

#### 5. Pipeline hangs or is very slow

**Diagnosis:**
```bash
# Check I/O wait
iostat -x 5
# High %util or await indicates disk bottleneck

# Check RAM
free -h
# If swap is being used, RAM is insufficient

# Check process status
ps aux | grep 24_process_reads
# Should show active processes
```

**Solution:**
```bash
# Reduce parallelism if I/O bound
export THREADS=16  # Instead of 32
./24_process_reads_optimized.sh -i DATA -t 16

# Or use fewer analysis flags
./24_process_reads_optimized.sh -i DATA  # Basic QC only
```

### Error Log Interpretation

**Standard error patterns:**

| Message | Severity | Action |
|---------|----------|--------|
| `[ERROR] Cannot create directory` | HIGH | Check permissions, disk space |
| `[ERROR] Validation failed` | MEDIUM | Check file integrity, wait if copying |
| `[WARNING] Time limit reached` | LOW | Normal - resume will continue |
| `[ERROR] Kraken2 failed` | HIGH | Check RAM, verify database |
| `[ERROR] Cannot repair corrupted FASTQ` | MEDIUM | File is unrecoverable, check sequencer |
| `[WARNING] No bins produced` | LOW | Normal for low-quality samples |

### Debug Mode

**Enable detailed logging:**

```bash
./24_process_reads_optimized.sh -i DATA -K -P -d

# Output shows:
# [DEBUG] Checking: /path/to/file.fastq.gz
# [DEBUG] Processing sample_001 - input size: 12345678 bytes
# [DEBUG] BBduk output size: 11234567 bytes
# [DEBUG] Running: /path/to/command
```

**Useful for diagnosing:**
- File discovery issues
- Command execution problems
- Size/quality issues
- Pipeline stage failures

---

## Maintenance

### Regular Tasks

**Daily:**
```bash
# Check logs for errors
grep '^\[ERROR\]' /data/results/*/*/log.txt | tail -20

# Verify disk space
df -h /data

# Check backup status
rsync -av --dry-run /data/results/ /backup/results/
```

**Weekly:**
```bash
# Clean old temporary files
find /data/.fastq_pass -name ".lock.*" -type d -mtime +7 -exec rmdir {} \;
find /data/results -name "*.tmp" -mtime +7 -delete
find /data/results -name "*.partial" -mtime +7 -delete

# Update databases (if connected to network)
cd /data/databases
kraken2-build --update-db kraken_db

# Review and archive old results
du -sh /data/results/202501*
# Move to long-term storage if needed
```

**Monthly:**
```bash
# Update dānaSeq
cd /data/danaSeq
git fetch origin
git log HEAD..origin/main  # Review new commits
git pull  # If updates available

# Update bioinformatics tools
conda update -n danaseq --all

# Review and clean logs
logrotate -f /etc/logrotate.d/danaseq

# Generate monthly statistics
/usr/local/bin/monthly-stats.sh
```

### Backup Procedures

**Critical data to backup:**

```bash
# Results (priority 1)
rsync -av --progress \
  /data/results/ \
  /backup/results/

# Raw data (priority 2)
rsync -av --progress \
  /data/minknow/ \
  /backup/minknow/

# Databases (priority 3 - can re-download)
rsync -av --progress \
  /data/databases/ \
  /backup/databases/

# Pipeline code (priority 4 - in git)
cd /data/danaSeq
git bundle create /backup/danaseq-$(date +%Y%m%d).bundle --all
```

**Automated backup script:**

```bash
#!/usr/bin/env bash
# /usr/local/bin/danaseq-backup.sh

BACKUP_ROOT="/backup/danaseq"
DATE=$(date +%Y%m%d)

# Backup results (incremental)
rsync -av --delete --link-dest="${BACKUP_ROOT}/latest" \
  /data/results/ \
  "${BACKUP_ROOT}/${DATE}/"

# Update latest symlink
ln -snf "${BACKUP_ROOT}/${DATE}" "${BACKUP_ROOT}/latest"

# Clean old backups (keep 30 days)
find "${BACKUP_ROOT}" -maxdepth 1 -type d -mtime +30 -exec rm -rf {} \;

echo "Backup completed: ${DATE}"
```

### Updates and Patches

**Check for updates:**

```bash
cd /data/danaSeq
git fetch origin

# View available updates
git log HEAD..origin/main --oneline

# Review changes
git diff HEAD..origin/main
```

**Apply updates (safe method):**

```bash
# 1. Backup current version
cp -r /data/danaSeq /data/danaSeq.backup.$(date +%Y%m%d)

# 2. Test updates on small dataset first
git pull origin main

# 3. Run test
cd 10_realtime_processing
./24_process_reads_optimized.sh -i /data/test -v

# 4. If test passes, proceed with production
# 5. If test fails, rollback
cd /data
rm -rf danaSeq
mv danaSeq.backup.YYYYMMDD danaSeq
```

---

## Rollback Procedures

### Emergency Rollback

**If new version causes issues:**

```bash
# 1. Stop all running pipelines
pkill -f 24_process_reads_optimized.sh

# 2. Restore previous version
cd /data
rm -rf danaSeq
mv danaSeq.backup.YYYYMMDD danaSeq

# 3. Verify rollback
cd danaSeq
git log --oneline -1
# Should show previous version

# 4. Resume processing
cd 10_realtime_processing
./24_process_reads_optimized.sh -i /data/minknow/LATEST -K -P
```

### Rollback to Specific Security Patch

**Pre-security audit version (NOT RECOMMENDED):**
```bash
git checkout <commit-before-a4d9ec2>
# WARNING: This removes all security fixes!
```

**Specific security level:**
```bash
# After CRITICAL fixes only
git checkout a4d9ec2

# After HIGH fixes
git checkout e66a6ef

# After MEDIUM fixes
git checkout 23fcf9c

# Latest (all fixes)
git checkout 3a8964d
```

---

## For LLM Assistants

### Quick Reference

**When assisting with dānaSeq deployment, prioritize:**

1. **Verify security status** - Confirm post-audit version (commit 3a8964d+)
2. **Check dependencies** - Parse `./status.sh` output
3. **Validate resources** - RAM ≥128GB for Kraken2, disk ≥1TB
4. **Test before production** - Small dataset first
5. **Monitor for errors** - Check `[ERROR]` patterns in logs

### Common User Requests

**"How do I run the pipeline?"**
```bash
# Standard workflow
cd /data/danaSeq/10_realtime_processing
./24_process_reads_optimized.sh -i /data/minknow/RUN_XXX -K -P -v
```

**"The pipeline failed, what's wrong?"**
```bash
# Diagnostic sequence
1. Check last 20 errors: grep '^\[ERROR\]' out_*/*/log.txt | tail -20
2. Check disk space: df -h /data
3. Check RAM: free -h
4. Check dependencies: ./status.sh
5. Review failed_files.txt
```

**"How do I resume an interrupted run?"**
```bash
# Simply re-run the same command
# Resume logic automatically skips completed files
./24_process_reads_optimized.sh -i INPUT -o EXISTING_OUTPUT -K -P
```

**"How do I process just QC without Kraken/Prokka?"**
```bash
# Omit -K and -P flags
./24_process_reads_optimized.sh -i INPUT
```

**"How do I add HMM search to existing results?"**
```bash
# Re-run with --hmm flag - will skip QC/Kraken/Prokka (already done)
./24_process_reads_optimized.sh -i INPUT -o EXISTING_OUTPUT -P --hmm FILE.hmm
```

### Diagnostic Prompts

**When user reports "not working":**
1. Ask for specific error message
2. Request: `tail -50 failed_files.txt`
3. Request: `grep '^\[ERROR\]' log.txt | tail -10`
4. Check: `df -h /data` and `free -h`

**When user reports "slow performance":**
1. Check parallelism: `ps aux | grep 24_process_reads | wc -l`
2. Check I/O: `iostat -x 5`
3. Check Kraken serialization: `ps aux | grep kraken2`
4. Suggest: reduce threads or disable Kraken2

**When user reports "corrupt files":**
1. Test integrity: `gzip -t FILE.fastq.gz`
2. Check if file is still being written: `lsof FILE.fastq.gz`
3. Check file size: `ls -lh FILE.fastq.gz`
4. Suggest: wait if actively sequencing, re-copy if complete

### Key Files to Check

```
Critical files for diagnostics:
- ./status.sh              (dependency verification)
- failed_files.txt         (failure tracking)
- */barcode*/log.txt       (detailed processing logs)
- SECURITY_AUDIT.md        (vulnerability documentation)
- CLAUDE.md                (architecture and concepts)
```

### Safety Reminders for LLMs

1. **NEVER** suggest `rm -rf` - always use `trash` or `trash-put`
2. **NEVER** suggest `git reset --hard` without backup confirmation
3. **NEVER** suggest `--force` flag without explaining resume behavior
4. **ALWAYS** check security patch status before modifications
5. **ALWAYS** suggest testing on small dataset first
6. **ALWAYS** verify disk space before starting large jobs

---

## Quick Command Reference

```bash
# Essential commands at a glance

# Basic run (QC only)
./24_process_reads_optimized.sh -i INPUT

# Standard run (QC + Kraken + Prokka)
./24_process_reads_optimized.sh -i INPUT -K -P -v

# Full run (add HMM search)
./24_process_reads_optimized.sh -i INPUT -K -P --hmm FILE.hmm -v

# Resume interrupted run
./24_process_reads_optimized.sh -i INPUT -o EXISTING_OUTPUT -K -P

# Force re-run
./24_process_reads_optimized.sh -i INPUT -K -P --force

# Check dependencies
./status.sh

# Monitor progress
watch -n 30 'find OUTPUT -name "*.fa" | wc -l'

# Check for errors
tail -f OUTPUT/failed_files.txt

# Debug mode
./24_process_reads_optimized.sh -i INPUT -d

# MAG assembly
cd ../20_mag_assembly
./61_map_and_bin_optimized.sh OUTPUT_DIR
```

---

## Support and Resources

### Documentation

- `README.md` - Project overview and quick start
- `METHODS.md` - Scientific methodology
- `SECURITY_AUDIT.md` - Complete security analysis
- `CLAUDE.md` - Architecture for LLM assistants
- `10_realtime_processing/README.md` - Real-time pipeline details
- `20_mag_assembly/README.md` - MAG assembly details

### Troubleshooting

- Check `failed_files.txt` in output directory
- Review `log.txt` files in each barcode directory
- Enable debug mode with `-d` flag
- Consult troubleshooting section in this guide

### Contact

- **GitHub Issues:** https://github.com/rec3141/danaSeq/issues
- **Email:** rec3141@gmail.com
- **Documentation:** All `.md` files in repository

---

## Appendix: Security Verification

### Post-Deployment Security Check

**Verify all security fixes are present:**

```bash
cd /data/danaSeq

# Check 1: Atomic lock implementation (CRITICAL)
if grep -q "mkdir.*lockfile" 10_realtime_processing/24_process_reads_optimized.sh; then
    echo "✓ CRITICAL: Race condition fix present"
else
    echo "✗ CRITICAL: Race condition fix MISSING"
fi

# Check 2: No eval usage (CRITICAL)
if ! grep -q '^[^#]*eval ' 10_realtime_processing/24_process_reads_optimized.sh; then
    echo "✓ CRITICAL: Command injection fix present"
else
    echo "✗ CRITICAL: eval still in use"
fi

# Check 3: Signal handlers (CRITICAL)
if grep -q "trap.*cleanup_on_exit.*EXIT" 10_realtime_processing/24_process_reads_optimized.sh; then
    echo "✓ CRITICAL: Signal handlers present"
else
    echo "✗ CRITICAL: Signal handlers MISSING"
fi

# Check 4: Path validation (CRITICAL)
if grep -q "validate_output_dir" 20_mag_assembly/61_map_and_bin_optimized.sh; then
    echo "✓ CRITICAL: Path traversal fix present"
else
    echo "✗ CRITICAL: Path validation MISSING"
fi

# Check 5: Error message standardization (LOW)
if ! grep -q '\[ERR\]' 10_realtime_processing/24_process_reads_optimized.sh; then
    echo "✓ LOW: Error messages standardized"
else
    echo "✗ LOW: Old error format present"
fi

echo ""
echo "Security verification complete"
```

**Expected output:** All checks should show ✓

---

**Document Version:** 2.0
**Last Updated:** 2026-02-06
**Deployment Status:** Production-Ready
**Security Status:** ✅ All vulnerabilities fixed

**Ready for shipboard deployment!** 🚢
