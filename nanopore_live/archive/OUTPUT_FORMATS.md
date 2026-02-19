# 📊 Output Format Reference

## **Progress Line Formats**

### ✅ **Successful Processing**
```
RUN: FBA73618_pass_barcode13_f90ff888_76b7e129_172 : BBDUK FILTLONG FASTA KRAKEN PROKKA DONE
```
- All stages completed successfully
- File processed completely

---

### ❌ **Stage Failure**
```
RUN: FBA73618_pass_barcode13_f90ff888_76b7e129_172 : BBDUK (FAILED)
```
- Shows which stage failed
- Processing stops at failure point
- Details logged to `failed_files.txt`

**With verbose mode (`-v`):**
```
RUN: FBA73618_pass_barcode13_f90ff888_76b7e129_172 : BBDUK (FAILED)
  └─ BBDuk: corrupted gzip - validation should have caught this
```
- Shows failure reason inline
- Helps with immediate debugging

---

### ⚠️ **No Data After Filtering**
```
RUN: FBA73618_pass_barcode13_f90ff888_76b7e129_172 : BBDUK (NO DATA - all reads filtered)
```
**Or:**
```
RUN: FBA73618_pass_barcode13_f90ff888_76b7e129_172 : BBDUK FILTLONG (NO DATA - all reads too short)
```
**Or:**
```
RUN: FBA73618_pass_barcode13_f90ff888_76b7e129_172 : BBDUK FILTLONG FASTA (NO DATA - empty after processing)
```

**What this means:**
- Not a failure - file processed successfully
- But no data remained after quality filtering
- Common for:
  - Low-quality samples
  - Very short reads
  - Empty barcodes (no DNA assigned to this barcode)

---

### 🔄 **Skipped (Resume)**
No output - file already processed in previous run.

**With debug mode (`-d`):**
```
[DEBUG] Skipping FBA73618_pass_barcode13_f90ff888_76b7e129_172 (already complete)
```

---

## **Stage Abbreviations**

| Output | Stage | What It Does |
|--------|-------|--------------|
| `BBDUK` | BBDuk QC | Remove adapters, low-quality bases, artifacts |
| `FILTLONG` | Filtlong filtering | Length/quality filtering |
| `FASTA` | Format conversion | FASTQ → FASTA |
| `SKETCH` | Sendsketch | Fast taxonomic sketch (if `-S` flag) |
| `KRAKEN` | Kraken2 | Precise taxonomic classification (if `-K` flag) |
| `PROKKA` | Prokka | Gene annotation (if `-P` flag) |
| `HMM` | HMM search | Custom HMM searches on genes (if `--hmm` flag) |
| `TETRA` | Tetramer analysis | 4-mer frequencies for binning (if `-T` flag) |
| `COMPRESS` | Compression | Gzip outputs (if `--compress` flag) |

---

## **Example Run Outputs**

### Small Dataset (Mixed Results)
```bash
./24_process_reads_optimized.sh -i test_data -K -t 4

[INFO] 14:06:15 Found 10 FASTQ files
[INFO] 14:06:15 Processing 10 files with 4 parallel workers

RUN: file_001 : BBDUK FILTLONG FASTA KRAKEN DONE
RUN: file_002 : BBDUK FILTLONG FASTA KRAKEN DONE
RUN: file_003 : BBDUK (FAILED)
RUN: file_004 : BBDUK FILTLONG (NO DATA - all reads too short)
RUN: file_005 : BBDUK FILTLONG FASTA KRAKEN DONE
RUN: file_006 : BBDUK FILTLONG FASTA KRAKEN DONE
RUN: file_007 : BBDUK FILTLONG FASTA KRAKEN DONE
RUN: file_008 : BBDUK (NO DATA - all reads filtered)
RUN: file_009 : BBDUK FILTLONG FASTA KRAKEN DONE
RUN: file_010 : BBDUK FILTLONG FASTA KRAKEN DONE

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
⚠️  WARNING: 1 file(s) failed processing
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

Failures by stage:
  1 BBDuk

Failed files:
  ❌ file_003.fastq.gz
     Stage: BBDuk
     Reason: corrupted gzip - validation should have caught this

[DONE] Output directory: out_test_data_20251201
```

**Result Summary:**
- ✅ 7 files processed successfully
- ⚠️ 2 files had no data (normal for some barcodes)
- ❌ 1 file failed (corrupted)

---

### Large Dataset (Resume from Interruption)
```bash
# First run (interrupted after 50 files)
./24_process_reads_optimized.sh -i data -K -t 16
^C  # Ctrl+C

# Resume
./24_process_reads_optimized.sh -i data -K -t 16

[INFO] 14:06:15 Found 1000 FASTQ files
[INFO] 14:06:15 Validating FASTQ files
[INFO] 14:06:20 Processing 1000 files with 16 parallel workers

# First 50 files: Silently skipped (already complete)
# Progress bar shows fast movement through these

RUN: file_051 : BBDUK FILTLONG FASTA KRAKEN DONE
RUN: file_052 : BBDUK FILTLONG FASTA KRAKEN DONE
...
RUN: file_1000 : BBDUK FILTLONG FASTA KRAKEN DONE

✅ All files processed successfully!

[DONE] Output directory: out_data_20251201
```

---

## **Verbosity Levels**

### Normal Mode (default)
```
RUN: file_001 : BBDUK FILTLONG FASTA KRAKEN DONE
RUN: file_002 : BBDUK (FAILED)
```
- One line per file
- Stage names shown
- Failures marked

### Verbose Mode (`-v`)
```
RUN: file_001 : BBDUK FILTLONG FASTA KRAKEN DONE
[VERBOSE] BBduk completed for file_001
[VERBOSE] Filtlong completed for file_001
...

RUN: file_002 : BBDUK (FAILED)
  └─ BBDuk: input file empty or missing
```
- Shows completion messages
- Shows failure reasons inline
- Displays command outputs

### Debug Mode (`-d`)
```
[DEBUG] Processing file_001 - input size: 45123456 bytes
RUN: file_001 : BBDUK [DEBUG] Running bbduk on file_001
[DEBUG] Running: /work/apps/bbmap/bbduk.sh in=...
[... full command output ...]
[DEBUG] BBduk output size: 42567890 bytes
FILTLONG [DEBUG] Running filtlong on file_001
[DEBUG] Filtlong output size: 38234567 bytes
FASTA [DEBUG] Converting to FASTA: file_001
[DEBUG] Final FASTA size: 38234567 bytes
KRAKEN [DEBUG] Running kraken on file_001 (acquiring Kraken lock...)
[DEBUG] Kraken DB: /data/scratch/refdbs/krakendb/pluspfp_08gb
DONE

[DEBUG] Skipping file_002 (already complete)
```
- Full diagnostic output
- All commands shown
- File sizes tracked
- Resume logic visible

---

## **Progress Bar**

GNU Parallel provides a live progress bar:

```
RUN: file_172 : BBDUK FILTLONG FASTA KRAKEN DONE
100% 1000:0=2h15m ./CMO2025/20251124/.../barcode15/file_999.fastq.gz
```

**Format:**
- `100%` - Percentage complete
- `1000:0` - Completed:Running
- `2h15m` - Estimated time to completion
- Shows current file being processed

---

## **End-of-Run Summary**

### All Success
```
✅ All files processed successfully!

[DONE] Output directory: out_data_20251201
```

### With Failures
```
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
⚠️  WARNING: 15 file(s) failed processing
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

Failures by stage:
  10 BBDuk
   3 Filtlong
   2 Reformat

Failed files (showing first 10):

  ❌ file_023.fastq.gz
     Stage: BBDuk
     Reason: corrupted gzip - validation should have caught this

  ❌ file_045.fastq.gz
     Stage: BBDuk
     Reason: input file empty or missing

  ... (full list in failed_files.txt)

Common fixes:
  • Corrupted files: Re-download or re-extract original data
  • Tool errors: Check log files in output directories
  • Empty files: These are skipped automatically (not errors)
  • Path issues: Verify BBMAP, FILTLONG paths are correct

Full failure log: out_data_20251201/failed_files.txt
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

[DONE] Output directory: out_data_20251201
```

---

## **Quick Visual Guide**

```
✅ DONE                  = Success
❌ (FAILED)              = Error, check logs
⚠️ (NO DATA - ...)       = Filtered to empty (not an error)
🔄 Silent skip           = Already processed (resume)
```

---

## **Interpreting Results**

### High Success Rate (>95%)
```
✅ All files processed successfully!
```
**Interpretation:** Pipeline working well, data quality good

### Some Empty Files (5-10%)
```
RUN: file_X : BBDUK (NO DATA - all reads filtered)
```
**Interpretation:** Normal - some barcodes empty or low quality

### Many Empty Files (>20%)
```
Many files: (NO DATA - all reads too short)
```
**Interpretation:** Check sequencing parameters:
- Read length too short for MinION?
- Wrong quality threshold?
- Wrong barcode demultiplexing?

### Few Failures (<5%)
```
⚠️  WARNING: 3 file(s) failed processing
```
**Interpretation:** Normal - probably corrupted transfers

### Many Failures (>5%)
```
⚠️  WARNING: 50 file(s) failed processing
Failures by stage:
  50 BBDuk
```
**Interpretation:** Systematic issue:
- BBDuk failures: Check BBMap installation/paths
- Filtlong failures: Check Filtlong installation
- Kraken failures: Check database path/integrity

---

## **Monitoring Long Runs**

### Terminal 1: Run Pipeline
```bash
./24_process_reads_optimized.sh -i data -K -t 16
```

### Terminal 2: Watch Failures
```bash
tail -f out_data_*/failed_files.txt
```

### Terminal 3: Count Progress
```bash
watch -n 5 'find out_data_*/ -name "*.fa" | wc -l'
```

---

**Created:** 2025-12-01
**Related:** `24_process_reads_optimized.sh`, `BUGFIX_RESILIENT_PROCESSING.md`
