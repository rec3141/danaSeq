# 🔍 FASTA Corruption Diagnosis & Fix

## **Symptom: NULL Bytes in FASTA Headers**

```
^@^@^@^@^@^@^@^@^@^@...  (many NULL bytes)
GATTCTAATGAGGCTATATGGC... (sequence without header)
>85e2a360-2cad-43d8-968b-0621d697cfb0
ACAAGAAAGTTGTCGGTGTCTTTGTG
```

**What's wrong:**
- NULL bytes (`^@` = 0x00) at start of file
- Sequence data without a header
- Missing or corrupted first FASTA record

---

## **Likely Causes**

### 1. **BBMap reformat.sh Bug or Corruption** (Most Likely)

The corruption happens BEFORE our processing:
```bash
# This step might have produced corrupted output:
reformat.sh in=$ftfile out=$fafile.tmp.fa

# We then process it:
cut -f1 -d' ' "$fafile.tmp.fa" > "$fafile.partial"
```

**If reformat.sh crashed or had a bug:**
- `.tmp.fa` file already contains NULL bytes
- We faithfully copy the corruption to final file
- Result: Corrupted `.fa` file

**Evidence:**
- NULL bytes at file start (not typical of interrupted writes)
- Sequence without header (suggests upstream corruption)
- Pattern doesn't match interrupted shell redirect behavior

### 2. **Disk I/O Error**

```bash
cut ... > "$fafile.partial"
# Disk error here → partial write with NULL padding
```

**Check for disk errors:**
```bash
dmesg | grep -i "I/O error"
dmesg | grep -i "disk"
```

### 3. **Sparse Files / Seek Operations**

NULL bytes can appear when:
- File created with `fseek()` or `lseek()` beyond EOF
- Filesystem stores "holes" as metadata
- Holes appear as NULL bytes when read

**From research:**
[Sparse files](https://en.wikipedia.org/wiki/Sparse_file) store consecutive NULL bytes as metadata, not actual data blocks. If BBMap or another tool used seek operations incorrectly, this could create sparse regions.

**Check if file is sparse:**
```bash
du -h file.fa    # Actual disk usage
ls -lh file.fa   # Apparent size
# If du << ls, it's sparse!
```

### 4. **NFS / Network Storage Issues**

NULL bytes are common with:
- [DFSR replication problems](https://serverfault.com/questions/679322/why-do-files-end-up-being-filled-entirely-with-null-bytes-after-replication-usin)
- Stale file handles
- Network interruptions during write

**Check if on NFS:**
```bash
df -T /path/to/output_dir
# If shows 'nfs' or 'nfs4', that could be the issue
```

---

## **What We Learned from Research** 📚

**MYTH:** Interrupted `cut` command (Ctrl+C) leaves NULL bytes
**REALITY:** [Interrupted shell redirects](https://unix.stackexchange.com/questions/348921/redirect-shell-output-to-file-even-when-i-kill-the-program) result in:
- Truncated file (ends where buffer was flushed)
- Empty file (if nothing flushed)
- **NOT NULL bytes**

**NULL bytes come from:**
1. Sparse files / seek operations
2. Disk/filesystem corruption
3. NFS/network issues
4. Application bugs in upstream tools

---

## **Diagnosis Steps**

### Step 1: How Many Files Are Affected?

```bash
cd /Users/ericcollins/Desktop/claude-code/dana/nanopore_live
./check_fasta_corruption.sh /path/to/output_dir

# Or search manually:
find output_dir -name "*.fa" -exec grep -l $'\x00' {} \;
```

**If 1-2 files:** Likely random corruption (I/O glitch, Ctrl+C timing)
**If many files:** Systematic issue (disk problem, pipeline bug)

### Step 2: Check Source Files

For the corrupted FASTA, trace back to source:
```bash
# Corrupted file: output/FC/barcode/fa/FBA73618_pass_barcode13_xxx.fa

# Check BBDuk output:
ls -lh output/FC/barcode/fq/FBA73618_pass_barcode13_xxx.fastq.gz

# Check if .partial or .tmp.fa still exist:
ls -lh output/FC/barcode/fa/FBA73618_pass_barcode13_xxx.fa.*
```

### Step 3: Check Original Input

```bash
# Find original input file:
find input_dir -name "FBA73618_pass_barcode13_xxx*.fastq.gz"

# Validate it:
gzip -t input_dir/.../FBA73618_pass_barcode13_xxx.fastq.gz
echo $?  # Should be 0 if OK
```

### Step 4: Check Disk Health

```bash
# Linux:
dmesg | tail -100
df -h  # Check disk space
iostat -x 1 5  # Check I/O errors

# macOS:
diskutil verifyVolume /
```

---

## **Immediate Fix: Remove Corrupted File**

```bash
# Find the corrupted file
CORRUPTED_FILE="output/FC/barcode/fa/FBA73618_pass_barcode13_xxx.fa"

# Remove it
rm "$CORRUPTED_FILE"

# Re-run pipeline (will reprocess just this file)
./24_process_reads_optimized.sh -i input_dir -K -t 16

# The file will be reprocessed from scratch
```

---

## **Prevention: Add Validation**

Let me add a validation step to catch this before the file is considered "done".

### Option A: Validate FASTA After Creation

Add to script after line 400 (after mv):

```bash
# After atomic rename
mv "$fafile.partial" "$fafile"

# Validate: Check first line is a header
first_line=$(head -1 "$fafile")
if [[ ! "$first_line" =~ ^\> ]]; then
  echo "ERROR: Corrupted FASTA - first line not a header!"
  echo "First line: ${first_line:0:80}"
  rm "$fafile"  # Remove corrupted file
  log_failure "$file" "FASTA validation" "corrupted header"
  return 0
fi

# Validate: Check no NULL bytes
if grep -q $'\x00' "$fafile"; then
  echo "ERROR: FASTA contains NULL bytes!"
  rm "$fafile"  # Remove corrupted file
  log_failure "$file" "FASTA validation" "NULL bytes detected"
  return 0
fi
```

### Option B: Stricter Atomic Write

Instead of shell redirect, use a more robust approach:

```bash
# Instead of:
cut -f1 -d' ' "$fafile.tmp.fa" > "$fafile.partial"

# Use:
if ! cut -f1 -d' ' "$fafile.tmp.fa" > "$fafile.partial"; then
  log_failure "$file" "Header cleaning" "cut command failed"
  return 0
fi

# Verify partial file exists and has content
if [[ ! -s "$fafile.partial" ]]; then
  log_failure "$file" "Header cleaning" "partial file empty"
  return 0
fi

# Verify first line is a header
if ! head -1 "$fafile.partial" | grep -q '^>'; then
  log_failure "$file" "Header cleaning" "invalid first line"
  rm "$fafile.partial"
  return 0
fi

# Now safe to rename
mv "$fafile.partial" "$fafile"
```

---

## **Root Cause Investigation**

### Was This During a Ctrl+C?

Check your shell history:
```bash
history | grep -B5 -A5 "24_process_reads_optimized"
```

**If you interrupted the pipeline:**
- The `cut` command was likely mid-execution
- Signal caught during buffer flush
- Partial write + atomic rename = corrupted file

### Is Disk Almost Full?

```bash
df -h output_directory
```

**If disk >95% full:**
- Writes may fail silently
- Files get partial data + NULL padding
- **Fix:** Free up space before re-running

### Is This NFS/Network Storage?

Network filesystems can have:
- Write buffering issues
- Stale file handles
- Incomplete writes during network glitches

**Fix:** Use local storage for processing, copy to NFS after

---

## **Long-Term Fix: Improved Robustness**

### 1. Add Checksum Validation

```bash
# After creating final FASTA:
md5sum "$fafile" > "$fafile.md5"

# On resume, verify:
if [[ -f "$fafile.md5" ]]; then
  if md5sum -c "$fafile.md5" --quiet; then
    return 0  # Valid, skip
  else
    rm "$fafile" "$fafile.md5"  # Corrupted, redo
  fi
fi
```

### 2. Write-and-Verify Pattern

```bash
# Write with verification
cut -f1 -d' ' "$fafile.tmp.fa" > "$fafile.partial" || {
  log_failure "$file" "Header cleaning" "cut failed"
  return 0
}

# Sync to disk (ensure buffers flushed)
sync "$fafile.partial"

# Verify size matches input (roughly)
input_size=$(stat -c%s "$fafile.tmp.fa")
output_size=$(stat -c%s "$fafile.partial")

if (( output_size == 0 )); then
  log_failure "$file" "Header cleaning" "output empty"
  return 0
fi

# Basic sanity: output shouldn't be bigger than input
if (( output_size > input_size )); then
  log_failure "$file" "Header cleaning" "output larger than input (corruption?)"
  return 0
fi
```

### 3. Trap Signals

Add signal handler at top of script:

```bash
cleanup_on_interrupt() {
  echo ""
  echo "[WARN] Interrupted! Cleaning up partial files..."

  # Find and remove any .partial files
  find "${OUTPUT}" -name "*.partial" -delete

  exit 130
}

trap cleanup_on_interrupt SIGINT SIGTERM
```

---

## **Quick Workaround Script**

```bash
#!/bin/bash
# fix_corrupted_fastas.sh - Remove corrupted FASTA files for re-processing

OUTPUT_DIR="$1"

if [[ -z "$OUTPUT_DIR" ]]; then
  echo "Usage: $0 <output_directory>"
  exit 1
fi

echo "Finding corrupted FASTA files in: $OUTPUT_DIR"

REMOVED=0

# Find files with NULL bytes
find "$OUTPUT_DIR" -name "*.fa" -type f | while read fafile; do
  if grep -q $'\x00' "$fafile" 2>/dev/null; then
    echo "Removing corrupted: $fafile"
    rm "$fafile"
    ((REMOVED++))
  fi
done

# Find files where first line isn't a header
find "$OUTPUT_DIR" -name "*.fa" -type f | while read fafile; do
  first_line=$(head -1 "$fafile" 2>/dev/null)
  if [[ ! "$first_line" =~ ^\> ]]; then
    echo "Removing malformed: $fafile (first line: ${first_line:0:60})"
    rm "$fafile"
    ((REMOVED++))
  fi
done

echo ""
echo "Removed $REMOVED corrupted files"
echo ""
echo "Re-run pipeline to reprocess:"
echo "  ./24_process_reads_optimized.sh -i input_dir -K -t 16"
```

---

## **Testing for Corruption**

### Quick Test

```bash
# Check specific file
file="output/FC/barcode/fa/suspect.fa"

# Look for NULL bytes
if hexdump -C "$file" | grep -q "00 00 00 00"; then
  echo "CORRUPTED: Contains NULL bytes"
fi

# Check first line
head -1 "$file" | cat -v
# Should start with ">" not "^@^@^@"
```

### Comprehensive Test

```bash
cd /Users/ericcollins/Desktop/claude-code/dana/nanopore_live
./check_fasta_corruption.sh /path/to/output_dir > corruption_report.txt
cat corruption_report.txt
```

---

## **Your Specific Case**

Based on your output:
```
^@^@^@^@... (NULL bytes)
GATTCTAATGAGGC... (sequence)
>85e2a360-2cad-43d8-968b-0621d697cfb0 (header)
```

**Diagnosis:**
1. First record completely missing (just NULL bytes)
2. Second record missing header (just sequence)
3. Third record OK (has header and sequence)

**Most likely:** File was interrupted during write, buffer had:
- NULL padding where header should be
- Sequence from record 2
- Proper record 3

**Recommended action:**
```bash
# Find the file
find output_dir -name "*.fa" -exec grep -l $'\x00' {} \;

# Remove it
rm /path/to/corrupted.fa

# Reprocess
./24_process_reads_optimized.sh -i input_dir -K -t 16
```

---

## **When to Worry**

### Don't Worry (Probably):
- 1-2 corrupted files out of 10,000
- Happened during Ctrl+C
- Can be reprocessed successfully

### Do Worry:
- Many corrupted files (>1%)
- Corruption after full pipeline completion
- Re-processing produces same corruption
- Disk errors in dmesg

**If worried:** Check disk health, free up space, use local storage

---

**Created:** 2025-12-01
**Status:** Diagnostic guide
**Next step:** Run `check_fasta_corruption.sh` to assess scope of issue
