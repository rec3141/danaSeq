# 🔄 Resume Logic Guide

## **How Resume Works**

The pipeline now has **smart stage-aware resume logic** that allows you to:
1. Resume interrupted runs
2. Add new analyses to already-processed files
3. Rerun specific stages without reprocessing everything

---

## **Resume Behavior**

### **Scenario 1: First Run**
```bash
./24_process_reads_optimized.sh -i data -P --hmm /path/to/CANT-HYD.hmm
```

**What happens:**
- Basic processing (BBDuk → Filtlong → FASTA): ✅ Runs
- Prokka: ✅ Runs
- HMM search: ✅ Runs

**Output:** `RUN: file : BBDUK FILTLONG FASTA PROKKA HMM DONE`

---

### **Scenario 2: Interrupted Run (Resume)**
```bash
./24_process_reads_optimized.sh -i data -P --hmm /path/to/CANT-HYD.hmm
^C  # Interrupted after 100 files

# Resume
./24_process_reads_optimized.sh -i data -P --hmm /path/to/CANT-HYD.hmm
```

**What happens:**
- Files 1-100: ✅ Skipped (all stages complete)
- Files 101+: ✅ Processed

**Output:**
```
[DEBUG] Skipping file_001 (all requested stages complete)
[DEBUG] Skipping file_002 (all requested stages complete)
...
RUN: file_101 : BBDUK FILTLONG FASTA PROKKA HMM DONE
```

---

### **Scenario 3: Add HMM to Already-Processed Files** ⭐
```bash
# Initial run without HMM
./24_process_reads_optimized.sh -i data -P

# Later: Add HMM search to existing files
./24_process_reads_optimized.sh -i data -P --hmm /path/to/CANT-HYD.hmm
```

**What happens:**
- Basic processing (BBDuk → Filtlong → FASTA): ⏭️ **Skipped** (already exists)
- Prokka: ⏭️ **Skipped** (already exists)
- HMM search: ✅ **Runs** (doesn't exist)

**Output:** `RUN: file : HMM DONE`

---

### **Scenario 4: Add Second HMM Database**
```bash
# Initial run with CANT-HYD
./24_process_reads_optimized.sh -i data -P --hmm /path/to/CANT-HYD.hmm

# Later: Add FOAM database
./24_process_reads_optimized.sh -i data -P --hmm /path/to/CANT-HYD.hmm,/path/to/FOAM.hmm
```

**What happens:**
- Basic processing: ⏭️ Skipped
- Prokka: ⏭️ Skipped
- CANT-HYD HMM: ⏭️ **Skipped** (already exists)
- FOAM HMM: ✅ **Runs** (doesn't exist)

**Output:** `RUN: file : HMM DONE`

---

### **Scenario 5: Run Just Prokka (No HMM)**
```bash
# Initial run without Prokka
./24_process_reads_optimized.sh -i data -K

# Later: Add Prokka
./24_process_reads_optimized.sh -i data -K -P
```

**What happens:**
- Basic processing: ⏭️ Skipped
- Kraken: ⏭️ Skipped (already exists)
- Prokka: ✅ **Runs** (doesn't exist)

**Output:** `RUN: file : PROKKA DONE`

---

## **How It Works Internally**

### **Stage Detection Logic**

```bash
# 1. Check if FASTA exists
if [[ -s "$fafile" ]]; then
  SKIP_BASIC=1  # Basic processing done

  # 2. Check if optional stages are needed
  if Prokka requested AND Prokka output doesn't exist; then
    NEED_OPTIONAL=1
  fi

  if HMM requested AND any HMM output doesn't exist; then
    NEED_OPTIONAL=1
  fi

  # 3. Skip entirely if nothing to do
  if NEED_OPTIONAL == 0; then
    return 0  # Skip this file
  fi
fi

# 4. Run basic processing if needed
if SKIP_BASIC == 0; then
  # BBDuk → Filtlong → FASTA
fi

# 5. Run optional stages (each has its own resume check)
if Prokka requested; then
  if Prokka output doesn't exist; then
    run Prokka
  fi
fi

if HMM requested; then
  for each HMM database; do
    if HMM output doesn't exist; then
      run HMM search
    fi
  done
fi
```

---

## **File Existence Checks**

### **Basic Processing**
- **Check:** Does `$bcdir/fa/$base.fa` exist?
- **If yes:** Skip BBDuk, Filtlong, FASTA conversion

### **Prokka**
- **Check:** Does `$bcdir/prokka/$base/PROKKA_*.tsv` exist?
- **If yes:** Skip Prokka annotation

### **HMM Search**
- **Check:** Does `$bcdir/hmm/$base.DBNAME.tsv` exist for EACH database?
- **If yes:** Skip that specific HMM database

### **Kraken**
- **Check:** Does `$bcdir/kraken/$base.kraken.tsv` exist?
- **If yes:** Skip Kraken classification

### **Tetramer**
- **Check:** Does `$bcdir/tetra/$base.lrn` exist?
- **If yes:** Skip tetramer analysis

---

## **Important Notes**

### ✅ **Advantages**

1. **Efficient reruns** - Only processes what's needed
2. **Incremental analysis** - Add new stages without reprocessing
3. **Crash-safe** - Can interrupt and resume anytime
4. **Mix-and-match** - Different flags per run

### ⚠️ **Considerations**

1. **Parameters changed?** - Pipeline won't know if you changed `--min-readlen` or other parameters. If parameters change, use `--force` or delete output and rerun.

2. **Want to force reprocessing?** - Use the `--force` flag:
   ```bash
   # Force all optional stages (Prokka, HMM, Kraken, etc.) to rerun
   ./24_process_reads_optimized.sh -i data -P --hmm /path/to/CANT-HYD.hmm --force
   ```

   Or manually delete the stage's output files:
   ```bash
   # Force Prokka to rerun
   rm -rf output/FC*/barcode*/prokka

   # Force HMM to rerun
   rm -rf output/FC*/barcode*/hmm

   # Force everything to rerun
   rm -rf output
   ```

3. **Partial Prokka runs** - If Prokka failed partway, its output directory might exist but be incomplete. The resume logic checks for `.tsv` files, so incomplete runs should be detected.

---

## **Examples**

### **Force HMM Rerun on Existing Data**

```bash
# Remove HMM results
find output -name "*.CANT-HYD.tsv" -delete
find output -name "*.CANT-HYD.tbl" -delete

# Rerun with HMM
./24_process_reads_optimized.sh -i data -P --hmm /path/to/CANT-HYD.hmm
# Only HMM stage will run
```

### **Add Multiple New Stages**

```bash
# Initial run: Just basic QC
./24_process_reads_optimized.sh -i data

# Later: Add Kraken
./24_process_reads_optimized.sh -i data -K

# Later: Add Prokka + HMM
./24_process_reads_optimized.sh -i data -K -P --hmm /path/to/CANT-HYD.hmm
```

### **Debug What Will Run**

```bash
# Use debug mode to see what will be skipped
./24_process_reads_optimized.sh -i data -P --hmm /path/to/CANT-HYD.hmm -d

# Look for:
# [DEBUG] Skipping file_001 (all requested stages complete)
# [DEBUG] Found Prokka proteins: ...
# [DEBUG] HMM results exist for CANT-HYD: ...
```

---

## **Troubleshooting**

### **HMM not running on existing Prokka output**

**Problem:** You ran with `-P`, then reran with `-P --hmm`, but HMM didn't run.

**Cause:** Likely a bug (should be fixed now!)

**Solution:**
```bash
# Check if Prokka .faa exists
ls -lh output/FC*/barcode*/prokka/*/PROKKA_*.faa

# Check if HMM output exists (shouldn't if this is the issue)
ls -lh output/FC*/barcode*/hmm/*.CANT-HYD.tsv

# Force HMM to run
rm -rf output/FC*/barcode*/hmm
./24_process_reads_optimized.sh -i data -P --hmm /path/to/CANT-HYD.hmm
```

### **Pipeline reruns everything unnecessarily**

**Problem:** Basic processing reruns even though `.fa` files exist.

**Possible causes:**
- `.fa` files are empty (size 0)
- `.fa` files are gzipped but script checks for uncompressed
- Wrong output directory specified

**Solution:**
```bash
# Check FASTA files exist and have content
find output -name "*.fa" -type f -exec ls -lh {} \; | head

# Verify not empty
find output -name "*.fa" -type f -size 0
# Should return nothing

# Check correct output directory
./24_process_reads_optimized.sh -i data -o output -P --hmm /path/to/HMM.hmm
```

---

**Created:** 2025-12-02
**Status:** Production-ready
**Related:** `24_process_reads_optimized.sh`, `HMM_SEARCH_GUIDE.md`
