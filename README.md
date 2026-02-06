# dānaSeq

**Genomics-Based Ecosystem Service Assessment Tool**

"Ecosystem Services" are the benefits that accrue to humans from having a healthy ecosystem. When natural forces or human activities disturb the ecosystem, these benefits may diminish or change. This software takes its name from *dāna*, the Buddhist concept of selfless giving and generosity in service of others.

dānaSeq conducts real-time analysis of Nanopore sequencing data to identify microbial communities, taxa and genes that affect Ecosystem Services using automated binning, genome assembly, and annotation of marker genes using a custom database of hidden markov models (HMMs) that build upon FOAM (Prestat et al. 2014), CANT-HYD (Khot et al. 2022), NCycDB (Tu et al. 2019), HADEG (Rojas-Vargas et al. 2023), HMDB (Wang et al. 2023), TASmania (Akarsu 2019), IDOPS (Díaz-Valerio et al. 2021).

---

```
╔═══════════════════════════════════════════════════════════════════════════════╗
║                                                                               ║
║   ███╗   ███╗███████╗████████╗ █████╗  ██████╗ ███████╗███╗   ██╗ ██████╗   ║
║   ████╗ ████║██╔════╝╚══██╔══╝██╔══██╗██╔════╝ ██╔════╝████╗  ██║██╔═══██╗  ║
║   ██╔████╔██║█████╗     ██║   ███████║██║  ███╗█████╗  ██╔██╗ ██║██║   ██║  ║
║   ██║╚██╔╝██║██╔══╝     ██║   ██╔══██║██║   ██║██╔══╝  ██║╚██╗██║██║   ██║  ║
║   ██║ ╚═╝ ██║███████╗   ██║   ██║  ██║╚██████╔╝███████╗██║ ╚████║╚██████╔╝  ║
║   ╚═╝     ╚═╝╚══════╝   ╚═╝   ╚═╝  ╚═╝ ╚═════╝ ╚══════╝╚═╝  ╚═══╝ ╚═════╝   ║
║                                                                               ║
║              🌊 OXFORD NANOPORE EDNA ANALYSIS PIPELINE 🌊                    ║
║                                                                               ║
║          Real-time metagenomic sequencing for oceanographic research         ║
║                                                                               ║
╚═══════════════════════════════════════════════════════════════════════════════╝
```

```
                     🦠      🧬      🦠      🧬      🦠
                  ╔═══════════════════════════════════════╗
                  ║  ATCGATCGATCGATCGATCGATCGATCGATCGATCG ║
                  ║  TAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAG ║
                  ╚═══════════════════════════════════════╝
                     🔬   DECODE THE OCEANS   🔬
```

## 🎯 What Is This?

A **LEGENDARY** bioinformatics pipeline for analyzing environmental DNA from the ocean in **REAL-TIME** during research expeditions. We're talking shipboard sequencing, live taxonomic classification, and MAG assembly while you're still at sea! 🚢

---

## 📁 Project Structure

```
🌊 DANA METAGENOMIC PIPELINE
│
├─ 🏃 10_realtime_processing/     ⚡ LIVE ANALYSIS AT SEA
│   ├─ 10s: 📦 Preprocessing (MinKNOW wrangling)
│   ├─ 20s: 🔄 Read processing (6 flavors of awesome)
│   ├─ 30s: 🔧 Parsing utilities (AWK magic)
│   ├─ 40s: 💾 DuckDB integration (FAST queries)
│   ├─ 50s: 📊 Visualization (plots galore)
│   └─ 60s: 🗺️  Interactive dashboards (where's the cyano?)
│
├─ 🧬 20_mag_assembly/            🔬 MAG RECONSTRUCTION
│   ├─ 10s: 🏗️  Assembly (Flye power)
│   ├─ 20s: 📍 Mapping (minimap2 precision)
│   ├─ 30s: 🗂️  Binning (3-tool consensus)
│   ├─ 40s: ✨ Polishing (Racon + Medaka shine)
│   ├─ 50s: 🏷️  Characterization (what IS this?)
│   ├─ 60s: 🎼 Complete pipelines (one-stop shop)
│   ├─ 70s: 🔄 Format conversion (Bandage ready)
│   ├─ 80s: 📈 Visualization (t-SNE dreams)
│   └─ 90s: 🌐 Integration (schema + ecosystem)
│
└─ 📦 30_archive/                 💀 THE GRAVEYARD
    └─ Old experiments & deprecated code (RIP)
```

---

## 🚀 Quick Start

### ⚡ Real-Time Processing (shipboard/field)

```bash
cd 10_realtime_processing

# 🎯 The OPTIMIZED version (AI-enhanced, recommended!)
./24_process_reads_optimized.sh -i <barcode_dir> -K -P -S

# ⚡ Need SPEED? Go fast mode!
./22_process_reads_fast.sh -i <barcode_dir> -P -S

# 📊 Watch it live on the dashboard
Rscript 60_edna_mapping_viz.r
```

> **⚠️ CRITICAL WARNING:** When using Kraken2 (`-K` flag), ONLY use `24_process_reads_optimized.sh`!
> Kraken loads 50-100GB database into RAM. Other scripts will try to run multiple instances in parallel and **crash your system!** 💥
> The optimized script uses semaphores to serialize only Kraken calls while keeping other steps parallel. See `CRITICAL_KRAKEN_BUG.md` for details.

### 🧬 MAG Assembly (post-expedition)

```bash
cd 20_mag_assembly

# 🎼 Full orchestra - Assembly → Mapping → Binning → Polish
./60_map_and_bin_optimized.sh

# 🎨 Visualize those beautiful bins
Rscript 80_plot_bins.R
```

---

## 🛠️ Utility Scripts

Before diving in, check out these helpful utilities:

### 🎨 `banner.sh` - Welcome Banner
Display a beautiful introduction to the pipeline!
```bash
./banner.sh
```
Perfect for showing off to collaborators or starting your terminal session with style! 🌊

### 📊 `status.sh` - Dependency Checker
Quickly check which tools are installed and what's missing:
```bash
./status.sh
```
Shows:
- ✅ Installed tools (with versions!)
- ❌ Missing dependencies
- 💡 Installation recommendations
- 📁 Directory structure verification

**Pro tip:** Run `./status.sh` before starting an expedition to make sure everything's ready!

### 🎭 `agents.sh` - Expert Advisors
Get specialized guidance from domain experts:
```bash
./agents.sh
```

Meet your team:
- 🌊 **The Oceanographer** - Sampling strategy, water masses, HABs
- 💻 **The Bioinformatician** - Pipeline optimization, debugging, HPC
- 🌊 **The Ocean** - Deep wisdom from the waters themselves
- 🦠 **The Microbial Ecologist** - Community ecology, metabolic guilds

Each agent provides:
- 📚 Domain-specific knowledge
- 💡 Best practices and tips
- ⚠️ Common pitfalls to avoid
- 🎯 Interpretation guidance
- 📖 Recommended reading

**Pro tip:** Consult the Oceanographer BEFORE sampling, the Bioinformatician DURING analysis!

---

## 🎨 File Naming Convention

```
  10_step_name.sh       ← Nice round numbers
     ↓
  20_next_step.sh       ← Plenty of gaps
     ↓
  30_another.sh         ← Room to grow!
     ↓
  ...future steps...    ← Add whenever!
```

**The system:**
- 🔟 **10s-20s:** Core processing (the essentials)
- 🔧 **30s-40s:** Secondary analysis (getting fancy)
- 🎯 **50s-60s:** Integration & pipelines (the big guns)
- 📊 **70s+:** Visualization & reporting (make it pretty!)

---

## ✨ Key Features

```
    ┌─────────────────────────────────────────────────────┐
    │  ⚡ REAL-TIME      Process as sequencer streams!   │
    │  🎯 MULTI-TOOL     3-way binning consensus         │
    │  ✅ QUALITY        QC at every step                │
    │  🚢 EXPEDITION     Shipboard-ready deployment      │
    │  📊 VISUALIZATION  Interactive maps & clustering   │
    │  🧬 MAG ASSEMBLY   High-quality genome bins        │
    │  🌊 EDNA FOCUS     Marine & freshwater optimized   │
    └─────────────────────────────────────────────────────┘
```

---

## 🎯 Target Applications

### 🌊 Marine & Freshwater Ecology
Track microbial communities across oceanographic gradients

### 🦠 Harmful Algal Blooms
Real-time monitoring of toxic cyanobacteria

### 🚨 Pathogen Surveillance
Waterborne pathogens & fecal indicators (E. coli, Vibrio, etc.)

### 🌍 Biodiversity Assessments
Who's out there? Full taxonomic profiling

### 🧬 Environmental DNA
Complete eDNA metabarcoding workflow

---

## 🔬 The Pipeline Flow

```
┌─────────────────┐
│ 🧪 SAMPLE       │  Collect water, filter, extract DNA
│ COLLECTION      │
└────────┬────────┘
         │
         ▼
┌─────────────────┐
│ 🧬 NANOPORE     │  Load into MinION/GridION/PromethION
│ SEQUENCING      │
└────────┬────────┘
         │
         ▼
┌─────────────────────────────────────────────────────┐
│ ⚡ REAL-TIME PROCESSING (10_realtime_processing/)   │
├─────────────────────────────────────────────────────┤
│  1. 📦 Validate & repair FASTQ                      │
│  2. 🧹 Quality filter (BBDuk + Filtlong)            │
│  3. 🏷️  Taxonomic classify (Kraken2)                │
│  4. 📝 Annotate genes (Prokka)                      │
│  5. 💾 Store in DuckDB                              │
│  6. 🗺️  Update live dashboard                       │
└────────┬────────────────────────────────────────────┘
         │
         ▼
┌─────────────────────────────────────────────────────┐
│ 🧬 MAG ASSEMBLY (20_mag_assembly/)                  │
├─────────────────────────────────────────────────────┤
│  1. 🏗️  Co-assemble with Flye                       │
│  2. 📍 Map reads back (minimap2)                    │
│  3. 🗂️  Bin contigs (SemiBin + MetaBAT + MaxBin)   │
│  4. 🤝 Consensus binning (DAS Tool)                 │
│  5. ✨ Polish bins (Racon + Medaka)                 │
│  6. ✅ Quality check (CheckM2)                      │
│  7. 🏷️  Taxonomic assign (Kaiju + GTDB)            │
└────────┬────────────────────────────────────────────┘
         │
         ▼
┌─────────────────┐
│ 📊 RESULTS      │  High-quality MAGs, taxonomy, abundance
│ & PUBLICATION   │  Interactive visualizations, reports
└─────────────────┘
```

---

## 🛠️ Dependencies

### Core Bioinformatics Tools

```
┌──────────────────────┬─────────────────────────────────────┐
│ 🧬 Sequencing        │ Oxford Nanopore MinKNOW             │
│ 🏗️  Assembly          │ Flye (metagenomic mode)            │
│ 📍 Mapping           │ minimap2                            │
│ 🏷️  Taxonomy          │ Kraken2, Kaiju                      │
│ 📝 Annotation        │ Prokka                              │
│ 🗂️  Binning           │ SemiBin2, MetaBAT2, MaxBin2        │
│ 🤝 Consensus         │ DAS Tool                            │
│ ✅ Quality           │ CheckM2                             │
│ 🧹 Preprocessing     │ BBTools, Filtlong                   │
│ ✨ Polishing         │ Racon, Medaka                       │
└──────────────────────┴─────────────────────────────────────┘
```

### Analysis & Visualization

```
┌──────────────────────┬─────────────────────────────────────┐
│ 📊 R packages        │ tidyverse, DuckDB, leaflet          │
│ 🐍 Python            │ Python 3.x                          │
│ 📈 Clustering        │ t-SNE, UMAP, graph-based            │
└──────────────────────┴─────────────────────────────────────┘
```

---

## 🌍 Active Expeditions

### 🚢 CMO2025
California to Mexico Oceanographic Survey

### 🧊 QEI2025
Queen Elizabeth Islands Arctic Expedition

> **Note:** Scripts contain hardcoded paths for these expeditions.
> Update paths before running on new projects!

---

## 📚 Documentation

Each directory contains detailed READMEs:
- 📖 `10_realtime_processing/README.md` - Real-time workflow guide
- 📖 `20_mag_assembly/README.md` - MAG assembly deep dive
- 📖 `30_archive/README.md` - What's deprecated

---

## 🎓 Citation

If this pipeline helps your research, buy the developer a coffee ☕ and cite appropriately!

---

## 💪 Power Tips

### 🔥 Optimized Scripts
Scripts ending in `_optimized.sh` are AI-enhanced versions with:
- Better error handling
- Smarter resource usage
- Progress tracking
- Resume capability

### ⚡ Fast Mode
When you need results NOW:
- Use `22_process_reads_fast.sh`
- Trade accuracy for speed
- Perfect for initial screening

### 🎯 Consensus > Single Tool
Always use multi-tool binning:
- SemiBin2 (deep learning)
- MetaBAT2 (TNF + coverage)
- MaxBin2 (marker genes)
- DAS Tool (combines all three)

### 📊 Live Dashboards
Keep `60_edna_mapping_viz.r` running during expeditions to watch taxonomy appear in real-time on a map! 🗺️

---

```
╔══════════════════════════════════════════════════════════════╗
║                                                              ║
║           🌊 DECODE THE OCEANS, ONE READ AT A TIME 🌊        ║
║                                                              ║
║                    🦠 → 🧬 → 💻 → 📊 → 🌍                    ║
║                                                              ║
╚══════════════════════════════════════════════════════════════╝
```

**Now go forth and sequence! 🚀🔬🌊**
