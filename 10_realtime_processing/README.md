```
╔═══════════════════════════════════════════════════════════════════════════╗
║                                                                           ║
║     🌊  REAL-TIME PROCESSING PIPELINE  🌊                                ║
║                                                                           ║
║            ⚡ LIVE ANALYSIS AT SEA ⚡                                     ║
║                                                                           ║
║     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~                   ║
║        🐟     🦑      🐡     🦈      🐙     🐠                            ║
║     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~                   ║
║                                                                           ║
║        Process DNA sequences as they stream from the sequencer!          ║
║                                                                           ║
╚═══════════════════════════════════════════════════════════════════════════╝
```

## 🎯 What's This Directory?

This is your **MISSION CONTROL** for real-time metagenomic analysis during expeditions! 🚢

While your sequencer is churning out reads, these scripts are:
- 🔍 Classifying taxonomy in real-time
- 💾 Storing results in lightning-fast DuckDB
- 🗺️ Updating live dashboards
- 🚨 Detecting harmful algae blooms
- 📊 Generating quality metrics

**Think of it as a metagenomic news ticker!** 📰

---

## 🌊 The Processing Pipeline

```
┌─────────────────────────────────────────────────────────────────────┐
│                                                                     │
│  MinKNOW Sequencer → 📦 → 🧹 → 🏷️ → 📝 → 💾 → 📊 → 🗺️              │
│                                                                     │
│  Raw FASTQ      →  Preprocess  →  QC  →  Classify  →  Store  →  Viz│
│                                                                     │
└─────────────────────────────────────────────────────────────────────┘
```

---

## 📚 Workflow Overview

### 📦 **10s: Preprocessing** - Wrangle that data!

```
    ╔════════════════════════════════╗
    ║   📥 MINKNOW OUTPUT            ║
    ║   ├─ FASTQ.gz files            ║
    ║   ├─ Sometimes corrupted! 😱   ║
    ║   └─ Need standardization      ║
    ╚════════════════════════════════╝
              ↓
    ╔════════════════════════════════╗
    ║   🔧 PREPROCESSING             ║
    ╚════════════════════════════════╝
```

**Scripts:**
- `10_preprocess_fastq.py` - 🐍 Python preprocessing utilities
- `11_minknow_copy.sh` - 📋 Copy files from MinKNOW output
- `12_minknow_rename.sh` - ✏️ Rename to standard format

---

### 🔄 **20s: Read Processing** - The main event!

```
╔═══════════════════════════════════════════════════════════════╗
║                                                               ║
║  🧬 FASTQ READS → [MAGIC HAPPENS HERE] → 🏷️ TAXONOMY + 📝 GENES ║
║                                                               ║
╚═══════════════════════════════════════════════════════════════╝
```

**Choose your fighter:**

| Script | Speed | Accuracy | Use Case | Emoji |
|--------|-------|----------|----------|-------|
| `20_process_reads_basic.sh` | 🐌 | ⭐⭐⭐⭐⭐ | Single samples, max quality | 🎯 |
| `21_process_reads_barcode.sh` | 🐌 | ⭐⭐⭐⭐⭐ | Standard multiplexed | 📊 |
| `22_process_reads_fast.sh` | 🚀 | ⭐⭐⭐ | URGENT! Need results NOW | ⚡ |
| `23_process_reads_fast2.sh` | 🚀 | ⭐⭐⭐ | Alternative fast mode | ⚡⚡ |
| `24_process_reads_optimized.sh` | 🏎️ | ⭐⭐⭐⭐⭐ | **RECOMMENDED!** AI-enhanced | 🏆 |
| `25_process_reads_catchup.sh` | 🏎️ | ⭐⭐⭐⭐ | Batch accumulated samples | 📦 |
| `26_process_reads_test.sh` | 🧪 | ⭐⭐⭐ | Testing new features | 🔬 |
| `27_process_reads_fast_viz.r` | - | - | Quick plots for fast mode | 📈 |

**🏆 TL;DR: Use `24_process_reads_optimized.sh` for most work!**

---

#### 🔬 What Happens During Processing?

```
┌──────────────────────────────────────────────────────────────────┐
│                                                                  │
│  Step 1: 🔍 FASTQ Validation & Repair                           │
│          └─ Fix corrupted files, verify integrity               │
│                                                                  │
│  Step 2: 🧹 Quality Filtering                                    │
│          ├─ BBDuk: Remove adapters, contaminants                │
│          └─ Filtlong: Length & quality filtering                │
│                                                                  │
│  Step 3: 🏷️ Taxonomic Classification (Kraken2)                   │
│          └─ What organisms are in here?                         │
│                                                                  │
│  Step 4: 📝 Gene Annotation (Prokka)                             │
│          └─ What genes do we have?                              │
│                                                                  │
│  Step 5: 🎯 Taxonomic Profiling (sendsketch)                     │
│          └─ Fast sketch-based classification                    │
│                                                                  │
│  Step 6: 🧮 Tetranucleotide Frequency (TNF)                      │
│          └─ Composition signatures for binning                  │
│                                                                  │
└──────────────────────────────────────────────────────────────────┘
```

---

### 🔧 **30s: Parsing Utilities** - AWK magic!

```
    ┌──────────────────┐
    │  📄 Raw Output   │
    └────────┬─────────┘
             │
       awk magic ✨
             │
             ▼
    ┌──────────────────┐
    │  📊 Clean Data   │
    └──────────────────┘
```

- `30_kraken_parse.awk` - 🔍 Parse Kraken2 output elegantly
- `31_taxonomy_colors.awk` - 🎨 Generate taxonomy color schemes

---

### 💾 **40s: Database Integration** - DuckDB FTW!

```
╔═══════════════════════════════════════════════════════════════╗
║                                                               ║
║  🦆 DUCKDB - The Secret Weapon for Real-Time Analysis 🦆      ║
║                                                               ║
║  Why DuckDB?                                                  ║
║  ✅ FAST analytical queries (OLAP)                            ║
║  ✅ Embedded database (no server needed!)                     ║
║  ✅ Works great on laptops & ships                            ║
║  ✅ SQL queries while sequencing! 🎉                          ║
║                                                               ║
╚═══════════════════════════════════════════════════════════════╝
```

**Database scripts:**

```
40_kraken_db.r           →  🏷️ Taxonomy classifications
41_krakenreport_db.r     →  📊 Summary reports
42_prokka_db.r           →  🧬 Gene annotations
43_sketch_db.r           →  🎯 Sketch profiles
44_tetra_db.r            →  🧮 TNF signatures
45_stats_db.r            →  📈 Assembly stats
46_log_db.r              →  📝 Processing logs
47_merge_db.r            →  🔗 Merge run databases
48_merge_all_db.r        →  🌍 Merge EVERYTHING
49_kraken_table.r        →  📋 Summary tables
```

**Example workflow:**
```r
# While at sea, query your growing database!
library(duckdb)
con <- dbConnect(duckdb(), "expedition_data.duckdb")

# How much Cyanobacteria do we have?
dbGetQuery(con, "SELECT sample_id, COUNT(*) as cyano_reads
                 FROM kraken_results
                 WHERE taxonomy LIKE '%Cyanobacteria%'
                 GROUP BY sample_id")

# 🔥 Real-time science! 🔥
```

---

### 📊 **50s: Visualization** - Pretty plots!

```
    📊 ┌─────────────────┐
       │  Quality Scores │
       │    ╱╲    ╱╲     │  ← Your beautiful data
       │   ╱  ╲  ╱  ╲    │
       │  ╱    ╲╱    ╲   │
       └─────────────────┘
```

- `50_plot_nano_qscores.r` - 📈 Quality score distributions
- `51_plot_general.r` - 🎨 General plotting utilities

---

### 🗺️ **60s: Interactive Dashboards** - The showstopper!

```
╔═══════════════════════════════════════════════════════════════════╗
║                                                                   ║
║             🗺️ EDNA GEOGRAPHIC MAPPING DASHBOARD 🗺️               ║
║                                                                   ║
║   ┌────────────────────────────────────────────────────────┐    ║
║   │  🌍 Interactive Map                                     │    ║
║   │                                                         │    ║
║   │    🔴 Station 1: High Cyanobacteria                    │    ║
║   │    🟡 Station 2: Moderate diversity                    │    ║
║   │    🟢 Station 3: Clean water                           │    ║
║   │                                                         │    ║
║   │  Click markers for:                                    │    ║
║   │  • Taxonomic breakdown                                 │    ║
║   │  • Read counts                                         │    ║
║   │  • Contamination flags                                 │    ║
║   │  • Time series plots                                   │    ║
║   │                                                         │    ║
║   └────────────────────────────────────────────────────────┘    ║
║                                                                   ║
╚═══════════════════════════════════════════════════════════════════╝
```

**The crown jewel:**
- `60_edna_mapping_viz.r` - 🗺️ Geographic mapping with contamination filtering

**Features:**
- 🌍 Interactive leaflet maps
- 📍 GPS-tagged samples
- 🔴 Real-time updates as data streams in
- 🚨 Contamination alerts (human, plant, lab reagents)
- 📊 Taxonomic pie charts per location
- 🎨 Beautiful color-coded taxonomy
- ⏱️ Time-series animations

**Run it:**
```r
# Fire up the dashboard!
Rscript 60_edna_mapping_viz.r

# Opens in browser: http://localhost:8080
# Watch your expedition data come alive! 🎉
```

---

## 🔥 Key Features

```
┌──────────────────────────────────────────────────────────────┐
│                                                              │
│  ⚡ REAL-TIME PROCESSING                                     │
│     └─ Process reads as they stream from sequencer          │
│                                                              │
│  🔄 INCREMENTAL UPDATES                                      │
│     └─ Smart resume: only process new data                  │
│                                                              │
│  💾 DUCKDB BACKEND                                           │
│     └─ Lightning-fast SQL queries during expeditions        │
│                                                              │
│  🎚️ MULTIPLE SPEED MODES                                     │
│     └─ Trade speed vs accuracy based on your needs          │
│                                                              │
│  🧹 CONTAMINATION DETECTION                                  │
│     └─ Flag human, plant, lab reagent DNA                   │
│                                                              │
│  📊 LIVE DASHBOARDS                                          │
│     └─ Watch taxonomy appear in real-time on maps!          │
│                                                              │
└──────────────────────────────────────────────────────────────┘
```

---

## 🚀 Quick Start Guide

### Basic Usage

```bash
# 1. 🏆 The recommended way (optimized, reliable)
./24_process_reads_optimized.sh /path/to/barcode_dir

# 2. ⚡ Need results in a hurry?
./22_process_reads_fast.sh /path/to/barcode_dir

# 3. 📊 Launch the dashboard (in another terminal)
Rscript 60_edna_mapping_viz.r
```

### At-Sea Workflow

```
Day 1: 🚢 Board the vessel
       └─ Set up MinKNOW sequencer
       └─ Start dashboard: Rscript 60_edna_mapping_viz.r

Day 2: 🧪 Collect first samples
       └─ Load into sequencer
       └─ Run: ./24_process_reads_optimized.sh barcode01
       └─ Watch magic happen on dashboard! ✨

Day 3+: 🔄 Continuous sampling
       └─ Script auto-detects new data
       └─ Dashboard updates automatically
       └─ Monitor for HABs, pathogens, biodiversity

End:   📊 Export DuckDB for post-expedition analysis
       └─ Generate publication-ready figures
```

---

## 💡 Pro Tips

### 🎯 Script Selection Guide

**Use `24_process_reads_optimized.sh` when:**
- ✅ You want the best quality + speed balance
- ✅ Standard expedition workflow
- ✅ You want AI-enhanced error handling
- ✅ Resume capability is important

**Use `22_process_reads_fast.sh` when:**
- ⚡ You need results in <1 hour
- ⚡ Initial screening for pathogens/HABs
- ⚡ Real-time decision-making required
- ⚡ Can re-run with full analysis later

**Use `21_process_reads_barcode.sh` when:**
- 🔬 Maximum quality is critical
- 🔬 Publication-grade analysis
- 🔬 Time is not a constraint

### 🎚️ Speed vs Accuracy Trade-offs

```
🐌 ←─────────────────────────────────────────────────→ 🚀

20_basic              24_optimized         22_fast
   ⭐⭐⭐⭐⭐                  ⭐⭐⭐⭐⭐                ⭐⭐⭐
   Slowest               BEST               Fastest
   Most accurate        Great accuracy      Good enough
```

### 🧹 Contamination Filtering

The scripts automatically flag:
- 🧍 Human DNA (from handling)
- 🌱 Plant DNA (pollen, terrestrial runoff)
- 🧪 Lab reagents (from extraction kits)

Check the dashboard for red flags! 🚩

### 💾 Managing DuckDB Files

```bash
# Check database size
ls -lh *.duckdb

# Query from command line
duckdb expedition.duckdb "SELECT COUNT(*) FROM kraken_results"

# Export to CSV
duckdb expedition.duckdb "COPY (SELECT * FROM kraken_results) TO 'results.csv'"
```

---

## 🎓 Understanding the Output

### Directory Structure

```
output/
├── 📁 barcode01/
│   ├── 📄 filtered_reads.fastq.gz        ← Quality-filtered reads
│   ├── 📄 kraken_output.txt              ← Taxonomic classification
│   ├── 📄 prokka/                        ← Gene annotations
│   │   ├── annotations.gff
│   │   ├── proteins.faa
│   │   └── ...
│   ├── 📄 sendsketch_results.txt         ← Taxonomic sketches
│   └── 📄 tetra_frequencies.txt          ← TNF profiles
│
├── 📁 barcode02/
│   └── ...
│
├── 💾 expedition.duckdb                  ← All results in one DB!
└── 📊 logs/                              ← Processing logs
```

---

## 🐛 Troubleshooting

### Problem: Script says "No new data"
```
🔍 Cause: All data already processed
✅ Solution: This is normal! Script is smart about skipping done work.
```

### Problem: MinKNOW files have weird names
```
🔍 Cause: MinKNOW output naming varies by version
✅ Solution: Run `12_minknow_rename.sh` first
```

### Problem: Dashboard won't load
```
🔍 Cause: Missing R packages or port conflict
✅ Solution:
   1. Install packages: install.packages(c("tidyverse", "leaflet", "DuckDB"))
   2. Check port: Try http://localhost:8081 instead
```

### Problem: Kraken2 says "database not found"
```
🔍 Cause: Database path incorrect
✅ Solution: Edit script header to point to your Kraken2 DB
```

---

## 📚 Further Reading

- 📖 Kraken2 manual: [github.com/DerrickWood/kraken2](https://github.com/DerrickWood/kraken2)
- 📖 Prokka docs: [github.com/tseemann/prokka](https://github.com/tseemann/prokka)
- 📖 DuckDB docs: [duckdb.org](https://duckdb.org)
- 📖 Oxford Nanopore: [nanoporetech.com](https://nanoporetech.com)

---

```
╔═══════════════════════════════════════════════════════════════╗
║                                                               ║
║           🌊 RIDE THE WAVE OF REAL-TIME GENOMICS 🌊           ║
║                                                               ║
║     ⚡ Fast processing • 💾 Smart storage • 🗺️ Live viz ⚡       ║
║                                                               ║
║                  Now go catch some microbes! 🦠               ║
║                                                               ║
╚═══════════════════════════════════════════════════════════════╝
```
