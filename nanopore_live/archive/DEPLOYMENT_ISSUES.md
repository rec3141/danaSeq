# 🚨 Deployment Issues & Known Gotchas

## **Common Failure Points When Sharing This Script**

This document lists issues that will break when users run this script on their systems for the first time, ordered by likelihood.

---

## 💥 **1. Hardcoded Paths (100% Guaranteed Break)**

### Problem
Lines 41-47 assume your specific HPC cluster layout:

```bash
BBMAP=${BBMAP:-/work/apps/bbmap}
DANADIR=${DANADIR:-/work/apps/dana}
PROKKA_BIN=${PROKKA_BIN:-/work/apps/prokka/bin/prokka}
FILTLONG=${FILTLONG:-/work/apps/Filtlong/bin/filtlong}
KRAKEN2=${KRAKEN2:-/usr/bin/kraken2}
KRAKEN_DB=${KRAKEN_DB:-/data/scratch/refdbs/krakendb/pluspfp_08gb}
APPS=${APPS:-/work/apps}
```

### Error Message
```
/work/apps/bbmap/bbduk.sh: No such file or directory
```

### Fix Options

**Option A: Environment variables (current approach)**
```bash
# User must set before running:
export BBMAP=/usr/local/bbmap
export DANADIR=/home/user/dana
export FILTLONG=/usr/local/bin/filtlong
./24_process_reads_optimized.sh -i data
```

**Option B: Config file**
```bash
# Create dana.conf:
BBMAP=/usr/local/bbmap
DANADIR=/home/user/dana
FILTLONG=/usr/local/bin/filtlong
KRAKEN_DB=/data/kraken2/pluspfp

# Script sources it:
[[ -f dana.conf ]] && source dana.conf
[[ -f ~/.dana.conf ]] && source ~/.dana.conf
[[ -f /etc/dana.conf ]] && source /etc/dana.conf
```

**Option C: Auto-detection**
```bash
# Try to find tools automatically
BBMAP=${BBMAP:-$(command -v bbmap.sh 2>/dev/null | xargs dirname 2>/dev/null)}
FILTLONG=${FILTLONG:-$(command -v filtlong 2>/dev/null)}
```

**Recommended:** Combine B + C, with A as override.

---

## 💥 **2. Hardcoded Cache Directory (95% Will Break)**

### Problem
Line 114:
```bash
CACHE_FASTQ="/data/.fastq_pass"
mkdir -p "${CACHE_FASTQ}"
```

- `/data/` doesn't exist on many systems (especially macOS)
- User might not have write permissions
- Could fill up a shared partition

### Error Message
```
mkdir: cannot create directory '/data/.fastq_pass': Permission denied
```

### Fix
```bash
# Use system temp dir with fallback
CACHE_FASTQ="${CACHE_FASTQ:-${TMPDIR:-/tmp}/dana_fastq_cache}"
mkdir -p "${CACHE_FASTQ}" || {
  echo "[ERROR] Cannot create cache directory: ${CACHE_FASTQ}"
  echo "[TIP] Set CACHE_FASTQ to a writable directory"
  exit 1
}
```

---

## 💥 **3. Missing DuckDB R Scripts (90% Will Break)**

### Problem
Script calls R scripts that might not exist:
```bash
Rscript ${DANADIR}/sketch-db.r $bcdir
Rscript ${DANADIR}/kraken-db.r $bcdir
Rscript ${DANADIR}/krakenreport-db.r $bcdir
Rscript ${DANADIR}/prokka-db.r $bcdir
Rscript ${DANADIR}/tetra-db.r $bcdir
```

These are wrapped in `|| true`, so they won't crash the pipeline, but:
- **Silent failures** - no DuckDB integration
- **No visualization** - dashboard won't work
- **Confusing** - users don't know why it's not working

### Error (in log files)
```
Rscript: /work/apps/dana/kraken-db.r: No such file or directory
```

### Fix
```bash
# Check once at startup
if [[ ! -f "${DANADIR}/kraken-db.r" ]]; then
  echo "[WARN] DuckDB integration disabled: R scripts not found in ${DANADIR}"
  echo "[INFO] Results will be saved but not integrated into database"
  DUCKDB_INTEGRATION=0
else
  DUCKDB_INTEGRATION=1
fi

# Later in code:
if (( DUCKDB_INTEGRATION )); then
  run_cmd "Rscript ${DANADIR}/kraken-db.r $bcdir" "$logfile" || true
fi
```

---

## 💥 **4. Kraken Database (85% Will Break)**

### Problem
Default database path won't exist:
```bash
KRAKEN_DB=${KRAKEN_DB:-/data/scratch/refdbs/krakendb/pluspfp_08gb}
```

Issues:
- Path doesn't exist
- Database format version mismatch
- Database is 100GB+ (hours to download)
- Multiple database types (Standard, PlusPF, PlusPFP, MiniKraken)

### Error Message
```
kraken2: database ("/data/scratch/refdbs/krakendb/pluspfp_08gb") does not contain necessary file taxo.k2d
```

### Fix
```bash
if (( RUN_KRAKEN )); then
  if [[ ! -d "${KRAKEN_DB}" ]]; then
    echo "[ERROR] Kraken2 database not found: ${KRAKEN_DB}"
    echo ""
    echo "Download databases from:"
    echo "  https://benlangmead.github.io/aws-indexes/k2"
    echo ""
    echo "Recommended databases:"
    echo "  - k2_viral (400MB) - Viruses only"
    echo "  - MiniKraken (8GB) - Bacteria/archaea/viruses"
    echo "  - k2_standard (50GB) - Bacteria/archaea/viruses/protozoa"
    echo "  - k2_pluspfp (75GB) - Standard + plants/fungi"
    echo ""
    echo "Then set: export KRAKEN_DB=/path/to/database"
    exit 1
  fi

  # Check database integrity
  if [[ ! -f "${KRAKEN_DB}/hash.k2d" ]] || [[ ! -f "${KRAKEN_DB}/taxo.k2d" ]]; then
    echo "[ERROR] Kraken2 database appears incomplete: ${KRAKEN_DB}"
    echo "[INFO] Database should contain: hash.k2d, opts.k2d, taxo.k2d"
    exit 1
  fi
fi
```

---

## 💥 **5. Missing GNU Parallel (70% Will Break)**

### Problem
Script requires **GNU parallel**, not the "moreutils" parallel:
- Ubuntu/Debian have two different `parallel` packages
- macOS doesn't include it by default
- Version must support `--bar` flag

### Error Messages
```
parallel: Error: --bar is not a valid option
```
Or:
```
bash: parallel: command not found
```

### Fix
```bash
# Check for GNU parallel
if ! command -v parallel &>/dev/null; then
  echo "[ERROR] GNU parallel not found"
  echo ""
  echo "Install on Ubuntu/Debian:"
  echo "  sudo apt install parallel"
  echo ""
  echo "Install on macOS:"
  echo "  brew install parallel"
  echo ""
  echo "Install on CentOS/RHEL:"
  echo "  sudo yum install parallel"
  exit 1
fi

# Check it's GNU parallel (has --bar flag)
if ! parallel --help 2>&1 | grep -q "\-\-bar"; then
  echo "[ERROR] Wrong version of parallel installed"
  echo "[INFO] You have 'moreutils parallel', need 'GNU parallel'"
  echo ""
  echo "Fix on Ubuntu/Debian:"
  echo "  sudo apt remove moreutils"
  echo "  sudo apt install parallel"
  exit 1
fi
```

---

## 💥 **6. Missing Bioinformatics Tools (60% Will Break)**

### Problem
If BBTools, Filtlong, Kraken2, Prokka aren't installed or in wrong location.

### Error Messages
```
/work/apps/bbmap/bbduk.sh: No such file or directory
filtlong: command not found
```

### Fix
```bash
check_tool() {
  local tool="$1"
  local path="$2"
  local install_hint="$3"

  if [[ ! -x "$path" ]] && [[ ! -f "$path" ]]; then
    echo "[ERROR] $tool not found: $path"
    [[ -n "$install_hint" ]] && echo "[INFO] $install_hint"
    return 1
  fi
  return 0
}

# Check required tools
echo "[INFO] Checking dependencies..."
MISSING=0

check_tool "BBMap" "${BBMAP}/bbduk.sh" "Download from: https://sourceforge.net/projects/bbmap/" || MISSING=1
check_tool "Filtlong" "${FILTLONG}" "Install: conda install -c bioconda filtlong" || MISSING=1

if (( RUN_KRAKEN )); then
  check_tool "Kraken2" "${KRAKEN2}" "Install: conda install -c bioconda kraken2" || MISSING=1
fi

if (( RUN_PROKKA )); then
  check_tool "Prokka" "${PROKKA_BIN}" "Install: conda install -c bioconda prokka" || MISSING=1
fi

if (( MISSING )); then
  echo ""
  echo "[TIP] Run ./status.sh to check all dependencies"
  exit 1
fi
```

---

## 💥 **7. Input Structure Assumptions (40% Will Break)**

### Problem
Script expects specific Nanopore directory structure:
- `fastq_pass/` subdirectories
- Filenames: `FLOWCELL_pass_barcodeXX_RUNID_FILEID.fastq.gz`
- Barcode format: `barcodeXX` (where XX = digits)

If data is from different source or renamed:
```
[INFO] Found 0 FASTQ files
[WARN] No FASTQ files found
```

### Fix Options

**Option A: Relaxed filename requirements**
```bash
# Current (strict):
find "${INPUT}" -regex '.*/fastq_pass/.*' -path '*barcode*'

# Relaxed:
find "${INPUT}" -type f -name '*.fastq.gz' -size +"${MIN_SIZE}"
```

**Option B: Better error message**
```bash
if (( NUMFILES == 0 )); then
  echo "[WARN] No FASTQ files found matching expected structure"
  echo ""
  echo "Expected directory structure:"
  echo "  input_dir/"
  echo "    └── fastq_pass/"
  echo "        └── barcodeXX/"
  echo "            └── FLOWCELL_pass_barcodeXX_xxx.fastq.gz"
  echo ""
  echo "Files in your input directory:"
  find "${INPUT}" -type f -name '*.fastq.gz' | head -5
  exit 1
fi
```

---

## 💥 **8. R/Python Dependencies (30% Will Break)**

### Problem
Script calls external interpreters:
- `Rscript` (R might not be installed)
- `python3` (usually present, but dependencies might be missing)
- R packages: tidyverse, DuckDB, etc.

### Error Messages
```
Rscript: command not found
python3: No module named 'hashlib'  # (rare, but possible)
```

### Fix
```bash
# Check interpreters exist
if (( RUN_SKETCH )) || (( RUN_KRAKEN )) || (( RUN_PROKKA )) || (( RUN_TETRA )); then
  if ! command -v Rscript &>/dev/null; then
    echo "[WARN] Rscript not found - DuckDB integration disabled"
    echo "[INFO] Install R to enable visualization features"
  fi
fi

if (( SHUFFLE )); then
  if ! command -v python3 &>/dev/null; then
    echo "[ERROR] python3 required for --shuffle mode"
    exit 1
  fi
fi
```

---

## 🛡️ **Recommended Pre-flight Check Function**

Add this comprehensive check at script startup:

```bash
################################################################################
# Pre-flight dependency checks
################################################################################
preflight_checks() {
  local errors=0
  local warnings=0

  echo "[INFO] Running pre-flight checks..."

  # Critical: GNU Parallel
  if ! command -v parallel &>/dev/null; then
    echo "[ERROR] GNU parallel not found"
    echo "        Install: apt install parallel (Ubuntu) or brew install parallel (macOS)"
    ((errors++))
  elif ! parallel --help 2>&1 | grep -q "\-\-bar"; then
    echo "[ERROR] Wrong parallel version (need GNU parallel, not moreutils)"
    ((errors++))
  fi

  # Critical: Core tools
  if [[ ! -d "${BBMAP}" ]]; then
    echo "[ERROR] BBMap not found: ${BBMAP}"
    echo "        Set: export BBMAP=/path/to/bbmap"
    ((errors++))
  fi

  if [[ ! -x "${FILTLONG}" ]]; then
    echo "[ERROR] Filtlong not found: ${FILTLONG}"
    echo "        Set: export FILTLONG=/path/to/filtlong"
    ((errors++))
  fi

  # Optional tools (only check if flags set)
  if (( RUN_KRAKEN )); then
    if [[ ! -x "${KRAKEN2}" ]]; then
      echo "[ERROR] Kraken2 not found: ${KRAKEN2}"
      ((errors++))
    fi
    if [[ ! -d "${KRAKEN_DB}" ]]; then
      echo "[ERROR] Kraken database not found: ${KRAKEN_DB}"
      echo "        Download from: https://benlangmead.github.io/aws-indexes/k2"
      ((errors++))
    fi
  fi

  if (( RUN_PROKKA )); then
    if [[ ! -x "${PROKKA_BIN}" ]]; then
      echo "[ERROR] Prokka not found: ${PROKKA_BIN}"
      ((errors++))
    fi
  fi

  # Warnings for DuckDB integration
  if [[ ! -d "${DANADIR}" ]]; then
    echo "[WARN] Dana scripts not found: ${DANADIR}"
    echo "       DuckDB integration and visualization will be disabled"
    ((warnings++))
  fi

  if ! command -v Rscript &>/dev/null; then
    echo "[WARN] Rscript not found - DuckDB integration disabled"
    ((warnings++))
  fi

  # Check cache directory is writable
  if ! mkdir -p "${CACHE_FASTQ}" 2>/dev/null; then
    echo "[ERROR] Cannot create cache directory: ${CACHE_FASTQ}"
    echo "        Set: export CACHE_FASTQ=/writable/path"
    ((errors++))
  fi

  # Summary
  if (( errors > 0 )); then
    echo ""
    echo "[ERROR] ${errors} critical issue(s) found. Cannot proceed."
    echo "[TIP] Run ./status.sh for detailed dependency information"
    return 1
  fi

  if (( warnings > 0 )); then
    echo "[WARN] ${warnings} warning(s) - some features will be disabled"
  fi

  echo "[OK] Pre-flight checks passed!"
  return 0
}

# Run checks before processing
preflight_checks || exit 1
```

---

## 📋 **Quick Setup Guide for New Users**

### Minimal Setup (QC only, no classification)

```bash
# 1. Install dependencies
sudo apt install parallel  # or: brew install parallel

# 2. Set environment variables
export BBMAP=/path/to/bbmap
export FILTLONG=/path/to/filtlong
export CACHE_FASTQ=/tmp/dana_cache

# 3. Run basic QC
./24_process_reads_optimized.sh -i /path/to/nanopore/data
```

### Full Setup (with Kraken classification)

```bash
# 1. Install all tools
sudo apt install parallel
# ... install BBMap, Filtlong, Kraken2 ...

# 2. Download Kraken database (one-time, large!)
# MiniKraken (8GB) - fastest
wget https://genome-idx.s3.amazonaws.com/kraken/k2_minusb_20220926.tar.gz
tar -xzf k2_minusb_20220926.tar.gz

# 3. Set environment
export BBMAP=/usr/local/bbmap
export FILTLONG=/usr/local/bin/filtlong
export KRAKEN2=/usr/local/bin/kraken2
export KRAKEN_DB=/data/kraken2/k2_minusb_20220926
export DANADIR=/path/to/dana
export CACHE_FASTQ=/tmp/dana_cache

# 4. Run with Kraken
./24_process_reads_optimized.sh -i /path/to/data -K
```

---

## 🔧 **Future Improvements**

### High Priority
1. ✅ **Add pre-flight checks** (see function above)
2. ✅ **Create config file support** (dana.conf)
3. ✅ **Better error messages** with installation hints
4. ✅ **Graceful degradation** (disable features if tools missing)

### Medium Priority
5. ⏳ **Auto-detect tool locations** (check PATH, common locations)
6. ⏳ **Conda environment** (environment.yml for easy setup)
7. ⏳ **Docker container** (eliminate all path issues!)
8. ⏳ **Setup wizard** (interactive config generator)

### Low Priority
9. ⏳ **Homebrew formula** (for macOS users)
10. ⏳ **Package for apt/yum** (system-wide installation)

---

## 📞 **Support Checklist**

When a user reports "it doesn't work," ask them to provide:

```bash
# 1. System info
uname -a
cat /etc/os-release  # Linux
sw_vers              # macOS

# 2. Tool versions
parallel --version
bbmap.sh --version
filtlong --version
kraken2 --version
prokka --version
Rscript --version

# 3. Environment variables
echo $BBMAP
echo $FILTLONG
echo $KRAKEN_DB
echo $DANADIR

# 4. Run status check
./status.sh

# 5. Try with verbose mode
./24_process_reads_optimized.sh -i data -v -d
```

---

**Last updated:** 2025-12-01
**Maintainer:** Dana Pipeline Team
**Related docs:** `README.md`, `CRITICAL_KRAKEN_BUG.md`, `status.sh`
