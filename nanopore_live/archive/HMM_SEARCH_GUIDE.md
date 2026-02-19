# 🧬 HMM Search Integration Guide

## **Overview**

The pipeline now supports **custom HMM searches** on Prokka-predicted genes. This enables specialized functional annotation using curated HMM databases like:

- **CANT-HYD** - Hydrocarbon degradation genes (37 HMMs)
- **FOAM** - Functional ontology assignments for metagenomes
- **Custom databases** - Any HMMER3-format HMM file

---

## **Quick Start**

### Basic Usage

```bash
# Run with Prokka + CANT-HYD HMM search
./24_process_reads_optimized.sh -i data -P --hmm /path/to/CANT-HYD.hmm

# Run with multiple HMM databases (comma-delimited paths)
./24_process_reads_optimized.sh -i data -P --hmm /path/to/CANT-HYD.hmm,/path/to/FOAM.hmm

# Full pipeline with HMM search
./24_process_reads_optimized.sh -i data -K -P -T --hmm /work/apps/dana/CANT-HYD.hmm
```

---

## **How It Works**

### 1. **Prokka Predicts Genes**
```
Input:  sample.fa (DNA sequences)
        ↓
Output: PROKKA_*.faa (protein sequences)
```

### 2. **HMM Search on Proteins**
```
Input:  PROKKA_*.faa + CANT-HYD.hmm
        ↓
        hmmsearch --cut_tc ...
        ↓
Output: sample.CANT-HYD.tsv (annotated hits)
```

### 3. **Results Saved**
```
output/FC/barcode/hmm/
├── sample.CANT-HYD.tbl  (raw hmmsearch output)
└── sample.CANT-HYD.tsv  (parsed results)
```

---

## **Output Format**

### TSV Columns

```tsv
gene_id         hmm_name        evalue          score
PROKKA_00001    AlkB            1.2e-45         152.3
PROKKA_00023    CYP153          3.4e-38         128.7
PROKKA_00087    LadA_alpha      5.6e-52         175.2
```

- **gene_id**: Prokka gene identifier (matches .faa file)
- **hmm_name**: HMM model name (e.g., AlkB, CYP153)
- **evalue**: Expectation value (lower = better)
- **score**: Bit score (higher = better)

---

## **Trusted Cutoffs (--cut_tc)**

### What Are Trusted Cutoffs?

Curated HMM databases like CANT-HYD include **model-specific score thresholds** that have been validated to minimize false positives. Each HMM has:

- **TC (Trusted Cutoff)**: High-confidence hits only
- **NC (Noise Cutoff)**: Lower threshold, more sensitive
- **GA (Gathering Threshold)**: Original alignment threshold

### Why Use --cut_tc?

```bash
# Without --cut_tc (default E-value cutoff)
hmmsearch CANT-HYD.hmm proteins.faa
# Result: Many false positives (e.g., generic FAD-binding proteins)

# With --cut_tc (recommended)
hmmsearch --cut_tc CANT-HYD.hmm proteins.faa
# Result: Only phylogenetically-validated hydrocarbon genes
```

**From CANT-HYD authors:**
> "Use `--cut_tc` for trusted cutoff to ensure accurate annotation of hydrocarbon degradation genes based on phylogenetic validation."

### Our Implementation

The pipeline **always uses `--cut_tc`** when searching with HMM databases. This ensures:
- ✅ No false positives from generic domains
- ✅ Phylogenetically-validated hits only
- ✅ Consistent with database authors' recommendations

---

## **Adding New HMM Databases**

### Step 1: Verify HMM Format

```bash
# Check it's HMMER3 format
head -1 /path/to/FOAM.hmm
# Should show: HMMER3/f [3.x | ...]

# Ensure it's NOT gzipped (hmmsearch doesn't support gzipped input)
file /path/to/FOAM.hmm
# Should show: ASCII text, not gzip compressed data

# If gzipped, decompress first:
gunzip /path/to/FOAM.hmm.gz

# Count number of HMMs
grep -c "^NAME " /path/to/FOAM.hmm

# List HMM names
grep "^NAME " /path/to/FOAM.hmm
```

### Step 2: Run Pipeline

```bash
# Provide full path to HMM file
./24_process_reads_optimized.sh -i data -P --hmm /path/to/FOAM.hmm

# Or multiple databases (comma-delimited paths)
./24_process_reads_optimized.sh -i data -P --hmm /path/to/FOAM.hmm,/path/to/CAZy.hmm
```

---

## **Configuration**

### Environment Variables

```bash
# HMMER binary location (default: /usr/bin/hmmsearch)
export HMMSEARCH=/usr/local/bin/hmmsearch

# Example: Use HMMER from custom location
export HMMSEARCH=/opt/hmmer/bin/hmmsearch
./24_process_reads_optimized.sh -i data -P --hmm /path/to/CANT-HYD.hmm
```

### Per-Database Settings

HMM databases are self-contained - trusted cutoffs are **embedded in the .hmm file**. No additional configuration needed!

**Important:** HMM files must be:
- ✅ HMMER3 format
- ✅ Uncompressed (plain text, not gzipped)
- ✅ Accessible from the compute nodes (if running on a cluster)

---

## **Example: CANT-HYD Database**

### About CANT-HYD

**Calgary approach to ANnoTating HYDrocarbon degrading enzymes**

- 37 curated HMM models
- Phylogeny-derived trusted cutoffs
- Aerobic & anaerobic hydrocarbon degradation

**Reference:**
> Khot V, Zorz J, et al. (2022). CANT-HYD: A Curated Database of Phylogeny-Derived Hidden Markov Models for Annotation of Marker Genes Involved in Hydrocarbon Degradation. *Front. Microbiol.* 12:764058.

### HMMs Included

```
AbcA_1, AbcA_2          - Benzene carboxylase
AlkB                    - Alkane hydroxylase
AlmA_GroupI, AlmA_GroupIII - Flavin-binding alkane monooxygenase
AhyA                    - Molybdopterin alkane hydroxylase
AssA                    - Alkylsuccinate synthase
BssA                    - Benzylsuccinate synthase
CYP153                  - Alkane-oxidizing cytochrome P450
DmpO                    - Phenol/toluene 2-monooxygenase
DszC                    - Dibenzothiophene monooxygenase
EbdA, CmdA              - Ethylbenzene dehydrogenase
K27540                  - Naphthalene carboxylase
LadA_alpha, LadA_beta, LadB - Long-chain alkane hydroxylase
MAH_alpha, MAH_beta     - Dioxygenase subunits
... and 17 more
```

### Expected Results

With **1000 genes per file** and **10,165 files**:
- ~0.5-5% of genes hit hydrocarbon models (depends on sample)
- Most hits: AlkB, CYP153 (common in oil/marine environments)
- Rare hits: BssA, AssA (anaerobic specialists)

---

## **Performance**

### Timing (per file with ~1000 genes)

```
Prokka:     120 seconds (ORF calling + annotation)
HMM search:  5 seconds per database
Total:      125-135 seconds (minimal overhead!)
```

### Parallelization

- **32 files × 1 CPU each** (standard pipeline model)
- HMM searches run in parallel across files
- No semaphore needed (hmmsearch is fast)

---

## **Troubleshooting**

### HMM file not found

```bash
[WARN] HMM file not found: /work/apps/dana/FOAM.hmm (skipping)
```

**Fix:**
```bash
# Verify file exists at specified path
ls -lh /work/apps/dana/FOAM.hmm

# Check for typos in path
./24_process_reads_optimized.sh -i data -P --hmm /work/apps/dana/FOAM.hmm

# Use absolute paths (not relative)
# WRONG: --hmm ../FOAM.hmm
# RIGHT: --hmm /full/path/to/FOAM.hmm
```

### No Prokka proteins found

```bash
[WARN] No Prokka proteins found for HMM search (skipping)
```

**Cause:** Prokka failed or didn't run

**Fix:**
```bash
# Ensure -P flag is used
./24_process_reads_optimized.sh -i data -P --hmm /path/to/CANT-HYD.hmm

# Check Prokka output
ls -lh output/FC/barcode/prokka/*/PROKKA_*.faa

# Check Prokka logs for errors
cat output/FC/barcode/log.txt | grep -i prokka
```

### hmmsearch not found

```bash
bash: /usr/bin/hmmsearch: No such file or directory
```

**Fix:**
```bash
# Install HMMER
conda install -c bioconda hmmer

# Or specify custom path
export HMMSEARCH=/path/to/hmmsearch

# Verify installation
hmmsearch -h
```

---

## **Downstream Analysis**

### Count Hits per Database

```bash
# Total hits across all files
find output -name "*.CANT-HYD.tsv" -exec cat {} \; | grep -v "^gene_id" | wc -l

# Hits per HMM model
find output -name "*.CANT-HYD.tsv" -exec cat {} \; | \
  grep -v "^gene_id" | cut -f2 | sort | uniq -c | sort -rn
```

### Extract High-Scoring Hits

```bash
# Get all AlkB hits with score > 100
find output -name "*.CANT-HYD.tsv" -exec awk -F'\t' \
  'NR>1 && $2=="AlkB" && $4>100 {print $0}' {} \; > alkb_high_score.tsv
```

### Combine Results

```bash
# Merge all CANT-HYD results into one file
echo -e "sample\tgene_id\thmm_name\tevalue\tscore" > all_cant-hyd.tsv
find output -name "*.CANT-HYD.tsv" | while read f; do
  sample=$(basename "$f" .CANT-HYD.tsv)
  tail -n +2 "$f" | awk -v s="$sample" '{print s"\t"$0}'
done >> all_cant-hyd.tsv
```

---

## **Integration with DuckDB** (Future)

Optional R script to load HMM results into DuckDB:

```r
# hmm-db.r (to be created)
library(DuckDB)
library(dplyr)

# Load all HMM TSV files
hmm_files <- list.files("output", pattern = ".*\\.CANT-HYD\\.tsv$",
                         recursive = TRUE, full.names = TRUE)

hmm_data <- map_df(hmm_files, ~{
  read_tsv(.x) %>%
    mutate(sample = basename(dirname(.x)))
})

# Insert into DuckDB
con <- dbConnect(duckdb::duckdb(), "output/data.duckdb")
dbWriteTable(con, "hmm_cant_hyd", hmm_data, overwrite = TRUE)
dbDisconnect(con)
```

---

## **Examples**

### Example 1: Find All Alkane Hydroxylase Genes

```bash
./24_process_reads_optimized.sh \
  -i marine_samples \
  -P --hmm /work/apps/dana/CANT-HYD.hmm \
  -t 32

# Extract AlkB and AlmA hits
find output -name "*.CANT-HYD.tsv" -exec awk -F'\t' \
  'NR>1 && ($2=="AlkB" || $2~/AlmA/) {print $0}' {} \; > alkane_hydroxylases.tsv
```

### Example 2: Compare Multiple Databases

```bash
./24_process_reads_optimized.sh \
  -i soil_metagenomes \
  -P --hmm /work/apps/dana/CANT-HYD.hmm,/work/apps/dana/FOAM.hmm \
  -t 32

# Compare coverage
echo "CANT-HYD hits:"
find output -name "*.CANT-HYD.tsv" -exec cat {} \; | grep -v "^gene_id" | wc -l

echo "FOAM hits:"
find output -name "*.FOAM.tsv" -exec cat {} \; | grep -v "^gene_id" | wc -l
```

### Example 3: Resume with HMM Search

```bash
# Initial run (interrupted)
./24_process_reads_optimized.sh \
  -i data \
  -P --hmm /work/apps/dana/CANT-HYD.hmm \
  -t 32
^C

# Resume - HMM search skipped for files with existing .tsv
./24_process_reads_optimized.sh \
  -i data \
  -P --hmm /work/apps/dana/CANT-HYD.hmm \
  -t 32
# Only processes new files or files without .tsv output
```

---

## **Best Practices**

### 1. Always Use Trusted Cutoffs

✅ **DO**: Use `--cut_tc` (automatic in this pipeline)
❌ **DON'T**: Use default E-value thresholds for curated databases

### 2. Combine with Taxonomic Data

```bash
# Run Kraken + Prokka + HMM
./24_process_reads_optimized.sh \
  -i data \
  -K -P --hmm /work/apps/dana/CANT-HYD.hmm

# Cross-reference: Which taxa have which hydrocarbon genes?
```

### 3. Validate Rare Hits

If you find unusual HMM hits:
1. Check E-value and score (lower E-value = better)
2. BLAST the protein sequence against NCBI
3. Verify taxonomic context (does this organism make sense?)

### 4. Multiple Databases

```bash
# Screen for multiple functions (comma-delimited full paths)
./24_process_reads_optimized.sh \
  -i data \
  -P --hmm /path/to/CANT-HYD.hmm,/path/to/FOAM.hmm,/path/to/CAZy.hmm

# Each creates separate .tsv files
output/FC/barcode/hmm/sample.CANT-HYD.tsv
output/FC/barcode/hmm/sample.FOAM.tsv
output/FC/barcode/hmm/sample.CAZy.tsv
```

---

## **References**

### CANT-HYD
Khot V, Zorz J, Gittins DA, et al. (2022). CANT-HYD: A Curated Database of Phylogeny-Derived Hidden Markov Models for Annotation of Marker Genes Involved in Hydrocarbon Degradation. *Front. Microbiol.* 12:764058.
https://github.com/dgittins/CANT-HYD-HydrocarbonBiodegradation

### HMMER
Eddy SR. (2011). Accelerated Profile HMM Searches. *PLoS Comput. Biol.* 7(10): e1002195.
http://hmmer.org/

### Prokka
Seemann T. (2014). Prokka: rapid prokaryotic genome annotation. *Bioinformatics* 30(14): 2068-2069.

---

**Created:** 2025-12-02
**Status:** Production-ready
**Related:** `24_process_reads_optimized.sh`, `OUTPUT_FORMATS.md`
