# 💾 Storage Optimization Plan

## **Problem Statement**

Nanopore sequencing generates massive text-based files that compress extremely well but are stored uncompressed, wasting terabytes of disk space.

---

## **Current Storage Usage**

### File Inventory (per processed file)

| File Type | Location | Size | Compressed? | Compression Ratio | Needed After? |
|-----------|----------|------|-------------|-------------------|---------------|
| Input FASTQ | `input/*.fastq.gz` | 10-500 MB | ✅ YES | - | Archive |
| BBDuk FASTQ | `fq/*.fastq.gz` | 8-400 MB | ✅ YES | - | No (intermediate) |
| Filtlong FASTQ | `fq/*.filt.fastq` | 20-1000 MB | ❌ **NO** | ~4:1 | No (temp, deleted) |
| Final FASTA | `fa/*.fa` | 20-1000 MB | ❌ **NO** | ~4:1 | Yes (MAG assembly) |
| Kraken TSV | `kraken/*.tsv` | 5-100 MB | ❌ **NO** | ~5:1 | Maybe (queries) |
| Kraken report | `kraken/*.report` | 1-10 MB | ❌ **NO** | ~5:1 | Maybe (summaries) |
| Prokka TSV | `prokka/*/*.tsv` | 1-50 MB | ❌ **NO** | ~3:1 | Maybe (annotations) |
| Prokka FAA | `prokka/*/*.faa` | 5-50 MB | ❌ **NO** | ~3:1 | Maybe (proteins) |
| Tetramer LRN | `tetra/*.lrn` | 5-50 MB | ❌ **NO** | ~3:1 | Yes (binning) |
| Log files | `*/log.txt` | 0.1-5 MB | ❌ **NO** | ~10:1 | Yes (debugging) |

### Calculation for User's Dataset (10,165 files)

**Assuming averages:**
- FASTA: 40 MB × 10,165 = **406 GB**
- Kraken TSV: 10 MB × 10,165 = **101 GB**
- Kraken reports: 2 MB × 10,165 = **20 GB**
- Prokka: 10 MB × 10,165 = **101 GB**
- Tetramer: 10 MB × 10,165 = **101 GB**
- Logs: 1 MB × 10,165 = **10 GB**

**Total uncompressed text: ~739 GB**

**After compression (4:1 average): ~185 GB**

**💰 SAVINGS: ~554 GB (75% reduction!)**

---

## **Compression Strategy Options**

### **Option 1: Aggressive Inline Compression** ⚡

Compress files immediately after creation, before next stage.

```bash
# After FASTA creation (line ~400):
gzip "$fafile"
fafile="${fafile}.gz"  # Update variable

# Kraken reads gzipped FASTA directly:
kraken2 --db "$KRAKEN_DB" "$fafile.gz" ...  # Works! ✅
```

**Pros:**
- ✅ Immediate space savings (compressed during run)
- ✅ If interrupted, completed files already compressed
- ✅ Most tools (Kraken, Prokka) support gzipped input
- ✅ No post-processing needed

**Cons:**
- ⚠️ Adds ~1-2 sec per file (gzip time)
- ⚠️ Slightly more complex code
- ⚠️ Need to handle .gz extensions everywhere

**Best for:** Real-time processing on ships with limited disk

---

### **Option 2: Post-File Compression** 🎯 **RECOMMENDED**

Compress after successfully processing each file.

```bash
# At end of process_one(), after all analyses:
if [[ -s "$fafile" ]]; then
  gzip -f "$fafile"
fi

if (( RUN_KRAKEN )); then
  gzip -f "$bcdir/kraken/$base.tsv" 2>/dev/null || true
  gzip -f "$bcdir/kraken/$base.report" 2>/dev/null || true
fi

# Resume check changes to:
[[ -s "$fafile" ]] || [[ -s "$fafile.gz" ]] && return 0
```

**Pros:**
- ✅ Clean separation (process → compress)
- ✅ Compressed files marked as "complete"
- ✅ Tools use uncompressed during processing (simpler)
- ✅ Easy to make optional (`--compress` flag)
- ✅ Resume detection works with .gz files

**Cons:**
- ⚠️ Disk holds uncompressed briefly during processing
- ⚠️ Adds ~1-2 sec per file

**Best for:** Most use cases (good balance)

---

### **Option 3: Batch Compression at End** 📦

Compress everything after pipeline finishes.

```bash
# After all files processed:
echo "[INFO] Compressing output files..."
find "${OUTPUT}" -type f \( -name "*.fa" -o -name "*.tsv" -o -name "*.report" -o -name "*.lrn" \) | \
  parallel --bar gzip -f {}
```

**Pros:**
- ✅ Zero overhead during processing
- ✅ Simple implementation
- ✅ User can skip if they want

**Cons:**
- ❌ No space savings until very end
- ❌ If interrupted, nothing compressed
- ❌ Requires full pipeline completion

**Best for:** Short runs on systems with plenty of disk

---

### **Option 4: Selective Compression** 🎨

Compress only specific file types, leave others readable.

```bash
# Compress large files immediately:
gzip "$fafile"  # FASTA (big, reused later)

# Leave results uncompressed for easy viewing:
# kraken/*.tsv (kept as-is)
# kraken/*.report (kept as-is)
```

**Pros:**
- ✅ Balance between space and accessibility
- ✅ Results stay human-readable
- ✅ Big files compressed

**Cons:**
- ⚠️ Less space savings (~60% instead of 75%)
- ⚠️ Inconsistent (some files .gz, some not)

**Best for:** Interactive analysis workflows

---

### **Option 5: Lazy Compression** 💤

Provide a separate script for users to compress manually.

```bash
#!/bin/bash
# compress_outputs.sh
find "$1" -type f \( -name "*.fa" -o -name "*.tsv" -o -name "*.report" \) | \
  parallel --bar gzip {}
```

**Pros:**
- ✅ No pipeline changes
- ✅ User control

**Cons:**
- ❌ Extra step users will forget
- ❌ No space savings during processing
- ❌ Not expedition-friendly

**Best for:** When you don't want to touch the pipeline

---

## **Recommended Implementation**

### **Tier 1: Option 2 (Post-File Compression)** ✨

This gives the best balance of space savings, safety, and user experience.

#### Changes Required:

**1. Add compression flag:**
```bash
# Line ~50: Add new option
COMPRESS_OUTPUT=1  # Default: enabled

# In argument parsing:
--no-compress) COMPRESS_OUTPUT=0; shift ;;
```

**2. Add compression helper:**
```bash
compress_output_files() {
  local base="$1"
  local bcdir="$2"

  # Compress FASTA (biggest file)
  [[ -f "$bcdir/fa/$base.fa" ]] && gzip -f "$bcdir/fa/$base.fa"

  # Compress Kraken results
  [[ -f "$bcdir/kraken/$base.tsv" ]] && gzip -f "$bcdir/kraken/$base.tsv"
  [[ -f "$bcdir/kraken/$base.report" ]] && gzip -f "$bcdir/kraken/$base.report"

  # Compress Prokka results
  if [[ -d "$bcdir/prokka/$base" ]]; then
    find "$bcdir/prokka/$base" -type f \( -name "*.tsv" -o -name "*.faa" -o -name "*.ffn" \) -exec gzip -f {} \;
  fi

  # Compress tetramer results
  [[ -f "$bcdir/tetra/$base.lrn" ]] && gzip -f "$bcdir/tetra/$base.lrn"
}
```

**3. Call at end of process_one():**
```bash
# After all optional analyses, before "echo DONE":
if (( COMPRESS_OUTPUT )); then
  echo -n "COMPRESS "
  compress_output_files "$base" "$bcdir"
fi

echo "DONE"
```

**4. Update resume check:**
```bash
# Change line 292:
[[ -s "$fafile" ]] && return 0

# To:
[[ -s "$fafile" ]] || [[ -s "$fafile.gz" ]] && return 0
```

**5. Export function:**
```bash
export -f compress_output_files
export COMPRESS_OUTPUT
```

---

## **Compression Tool Comparison**

| Tool | Speed | Ratio | CPU | Best For |
|------|-------|-------|-----|----------|
| `gzip` | Fast | 4:1 | Low | Default (widely compatible) |
| `pigz` | Very Fast | 4:1 | Multi-core | Large files (parallel gzip) |
| `bzip2` | Slow | 5:1 | High | Maximum compression |
| `xz` | Very Slow | 6:1 | Very High | Archival (rarely needed) |
| `zstd` | Very Fast | 4:1 | Medium | Modern alternative |

**Recommendation:** Use `pigz` if available (parallel), fall back to `gzip`.

```bash
# Detect best compressor
if command -v pigz &>/dev/null; then
  COMPRESSOR="pigz"
else
  COMPRESSOR="gzip"
fi

compress_output_files() {
  # ...
  $COMPRESSOR -f "$bcdir/fa/$base.fa"
}
```

---

## **Tool Compatibility Check**

Which downstream tools can read gzipped files?

| Tool | Reads .gz? | Notes |
|------|-----------|-------|
| `kraken2` | ✅ YES | Native support |
| `prokka` | ✅ YES | Via `--proteins file.fa.gz` |
| `bbmap` | ✅ YES | All BBTools handle .gz |
| `minimap2` | ✅ YES | Reference and query can be .gz |
| `samtools` | ✅ YES | Handles .gz transparently |
| `grep/awk` | ⚠️ NO | Need `zgrep`/`zcat` |
| Custom scripts | ⚠️ MAYBE | Depends on implementation |

**Bottom line:** Most bioinformatics tools handle .gz natively! ✅

---

## **Edge Cases & Considerations**

### 1. Resume After Compression
**Scenario:** File compressed in previous run, now resuming.

**Solution:**
```bash
# Check for both .fa and .fa.gz
if [[ -s "$fafile" ]] || [[ -s "$fafile.gz" ]]; then
  return 0  # Already processed
fi
```

### 2. Partial Compression on Crash
**Scenario:** Crash during gzip, file partially compressed.

**Solution:**
```bash
# gzip creates .gz atomically (writes to temp, then renames)
# If interrupted, either:
# - .fa exists (not compressed yet) → Safe ✅
# - .fa.gz exists, .fa gone (compression complete) → Safe ✅
# - Both exist (mid-rename, rare) → Resume will re-compress .fa ✅
```

### 3. Disk Full During Compression
**Scenario:** Gzip fails due to insufficient space.

**Solution:**
```bash
if ! gzip -f "$fafile"; then
  echo "[WARN] Compression failed for $base (disk full?)" >&2
  # Keep uncompressed version
fi
```

### 4. User Wants Uncompressed Files
**Scenario:** User prefers human-readable files.

**Solution:**
```bash
# Add flag:
./24_process_reads_optimized.sh -i data -K --no-compress

# Or decompress later:
find output_dir -name "*.gz" | parallel gunzip {}
```

### 5. Database Integration
**Scenario:** DuckDB R scripts expect .fa, not .fa.gz

**Solution:** Update R scripts to handle both:
```r
# In kraken-db.r:
files <- c(
  list.files(path, pattern = "\\.fa$", full.names = TRUE),
  list.files(path, pattern = "\\.fa\\.gz$", full.names = TRUE)
)
```

---

## **Performance Impact Analysis**

### Time Cost per File

| Operation | Time (MB/s) | 40MB File |
|-----------|-------------|-----------|
| gzip single-thread | ~50 MB/s | ~0.8 sec |
| pigz (4 cores) | ~200 MB/s | ~0.2 sec |
| gzip all results | - | ~2 sec total |

**Impact on 10,165 files:**
- Single-threaded gzip: +5.6 hours
- Parallel pigz: +1.4 hours

**But:** Compression happens *after* slow steps (Kraken ~60 sec/file)
- While Worker 1 compresses, Workers 2-16 do Kraken on other files
- Real-world impact: **Minimal (<5% slowdown)**

---

## **User Experience**

### Default Behavior (Compressed)
```bash
./24_process_reads_optimized.sh -i data -K

RUN: file001 : BBDUK FILTLONG FASTA KRAKEN COMPRESS DONE
RUN: file002 : BBDUK FILTLONG FASTA KRAKEN COMPRESS DONE

[DONE] Output directory: out_data_20251201
[INFO] Compressed 10165 FASTA files (saved ~380 GB)
```

### Optional: Skip Compression
```bash
./24_process_reads_optimized.sh -i data -K --no-compress

RUN: file001 : BBDUK FILTLONG FASTA KRAKEN DONE
RUN: file002 : BBDUK FILTLONG FASTA KRAKEN DONE

[DONE] Output directory: out_data_20251201
```

### Manual Compression Later
```bash
# If user wants to compress after the fact
./compress_outputs.sh out_data_20251201

[INFO] Compressing FASTA files...
100% 10165:0=2m30s
[INFO] Saved 380 GB
```

---

## **Storage Savings Projection**

### Small Run (100 files)
- Uncompressed: ~40 GB
- Compressed: ~10 GB
- **Saved: 30 GB** (enough for another dataset!)

### Medium Run (1,000 files)
- Uncompressed: ~400 GB
- Compressed: ~100 GB
- **Saved: 300 GB** (half a terabyte!)

### Large Run (10,000 files - typical expedition)
- Uncompressed: ~4 TB
- Compressed: ~1 TB
- **Saved: 3 TB** (multiple expeditions on one disk!)

### Annual Production (100,000 files)
- Uncompressed: ~40 TB
- Compressed: ~10 TB
- **Saved: 30 TB** (might not even need to buy new disks!)

---

## **Migration Plan**

### Phase 1: Add Compression (This PR)
1. Add `--compress` flag (default ON)
2. Add `compress_output_files()` function
3. Update resume checks for .gz
4. Test with small dataset
5. Document in README

### Phase 2: Update Downstream Tools (Next PR)
1. Update R scripts to handle .gz files
2. Test MAG assembly with .fa.gz inputs
3. Update visualization scripts
4. Document compatibility

### Phase 3: Batch Compression Tool (Nice to have)
1. Create `compress_outputs.sh` helper
2. Add `decompress_outputs.sh` helper
3. Add to utility scripts section

### Phase 4: Advanced Compression (Future)
1. Investigate `pigz` for parallel compression
2. Consider `zstd` for even better ratio/speed
3. Add compression level tuning (`-1` to `-9`)

---

## **Testing Plan**

### Test 1: Basic Compression
```bash
# Run with compression
./24_process_reads_optimized.sh -i test_data -K --compress

# Verify files compressed
find output_dir -name "*.fa.gz" | wc -l  # Should match input count
find output_dir -name "*.fa" -not -name "*.fa.gz" | wc -l  # Should be 0
```

### Test 2: Resume After Compression
```bash
# Process 10 files
./24_process_reads_optimized.sh -i test_data -K --compress

# Interrupt (Ctrl+C)

# Resume
./24_process_reads_optimized.sh -i test_data -K --compress

# Verify: Skips compressed files, continues uncompressed
```

### Test 3: No-Compress Mode
```bash
./24_process_reads_optimized.sh -i test_data -K --no-compress

# Verify: No .gz files created
find output_dir -name "*.gz" | wc -l  # Should be 0 (except input)
```

### Test 4: Space Savings
```bash
# Before compression
du -sh output_dir

# Compress
./compress_outputs.sh output_dir

# After compression
du -sh output_dir

# Calculate savings
```

### Test 5: Tool Compatibility
```bash
# Verify Kraken can read .fa.gz
kraken2 --db $DB test.fa.gz

# Verify minimap2 can read .fa.gz
minimap2 -x map-ont ref.fa.gz query.fa.gz
```

---

## **Rollout Strategy**

### Week 1: Development
- Implement Option 2 (Post-File Compression)
- Add tests
- Update documentation

### Week 2: Testing
- Run on small dataset (100 files)
- Verify space savings
- Check resume capability
- Test with/without compression

### Week 3: Validation
- Run on medium dataset (1,000 files)
- Monitor performance impact
- Verify downstream tool compatibility
- Get user feedback

### Week 4: Production
- Deploy to expedition systems
- Update user guides
- Create compression helpers
- Monitor real-world usage

---

## **Alternatives Considered**

### Why Not Use Lossy Compression?
- **CRAM format** - Lossy FASTQ compression
  - ❌ Not suitable: We need original quality scores
  - ❌ Limited tool support

### Why Not Archive to Tape?
- **LTO tape** - Long-term cold storage
  - ✅ Good for archival
  - ❌ Not useful during expedition (need random access)
  - 💡 Could be future Phase 5

### Why Not Stream to Cloud?
- **S3 Glacier** - Cloud archival
  - ❌ No internet at sea
  - 💡 Good for post-expedition archival

---

## **Conclusion**

**Recommended Implementation:** Option 2 (Post-File Compression)

**Expected Benefits:**
- 💰 75% storage reduction (~554 GB saved for user's dataset)
- ⚡ Minimal performance impact (<5% slowdown)
- 🛡️ Maintains crash safety (atomic gzip)
- 🔄 Resume-compatible
- 🎯 User-controllable (--compress/--no-compress)
- 🚀 Expedition-ready (compress during run)

**Implementation Effort:** ~2 hours
**Testing Effort:** ~4 hours
**Documentation:** ~1 hour

**Total: ~1 day of work for 75% storage savings!**

---

**Next Steps:**
1. Get user approval on approach
2. Implement Option 2
3. Test on small dataset
4. Deploy to production

**Questions for User:**
- Do you want compression ON by default? (Recommended: Yes)
- Should we compress Prokka results? (Large but text-based)
- Keep log files uncompressed for easy debugging?
- Use pigz if available for faster compression?

---

**Created:** 2025-12-01
**Status:** 📋 Planning phase
**Related:** `README.md`, `24_process_reads_optimized.sh`
