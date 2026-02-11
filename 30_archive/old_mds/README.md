# dānaSeq

**Real-time Metagenomic Analysis for Ecosystem Service Assessment**

Named after the Buddhist concept of *dāna* (selfless giving), this pipeline conducts shipboard analysis of Oxford Nanopore sequencing data to assess microbial ecosystem services in aquatic environments.

---

## Overview

dānaSeq provides a complete workflow for real-time metagenomic analysis during oceanographic expeditions, from raw nanopore reads to high-quality metagenome-assembled genomes (MAGs) with functional annotation.

**Core Capabilities:**
- Real-time taxonomic classification during sequencing
- Automated quality control and preprocessing
- Consensus genome binning (SemiBin2, MetaBAT2, MaxBin2)
- Functional gene annotation using custom HMM databases
- Interactive geographic visualization
- MAG assembly and polishing to MIMAG standards

**Functional Gene Databases:**
FOAM (Prestat et al. 2014), CANT-HYD (Khot et al. 2022), NCycDB (Tu et al. 2019), HADEG (Rojas-Vargas et al. 2023), HMDB (Wang et al. 2023), TASmania (Akarsu 2019), IDOPS (Díaz-Valerio et al. 2021)

---

## Architecture

```
dānaSeq Pipeline Structure

10_realtime_processing/     Real-time analysis during sequencing
├─ 10s  Preprocessing       MinKNOW output validation
├─ 20s  Read processing     QC, classification, annotation
├─ 30s  Parsing utilities   Data transformation
├─ 40s  Database ops        DuckDB integration
├─ 50s  Visualization       Quality metrics
└─ 60s  Dashboards          Interactive mapping

20_mag_assembly/            Post-expedition genome assembly
├─ 10s  Assembly            Flye metagenomic assembly
├─ 20s  Mapping             Read alignment and coverage
├─ 30s  Binning             Multi-tool consensus binning
├─ 40s  Polishing           Racon + Medaka refinement
├─ 50s  Characterization    Taxonomy and quality assessment
├─ 60s  Pipelines           End-to-end workflows
├─ 70s  Format conversion   Interoperability
├─ 80s  Visualization       Ordination and clustering
└─ 90s  Integration         Ecosystem service prediction

30_archive/                 Deprecated code
```

**Numbering Convention:**
Scripts are numbered in increments of 10 (10, 20, 30...) to allow insertion of intermediate steps without renumbering the entire pipeline.

---

## Quick Start

### Prerequisites

Check installed dependencies:
```bash
./status.sh
```

### Real-Time Processing

Recommended workflow for shipboard analysis:

```bash
cd 10_realtime_processing

# Standard pipeline (QC + taxonomy + annotation)
./24_process_reads_optimized.sh -i <input_dir> -K -P -S

# With functional gene profiling
./24_process_reads_optimized.sh -i <input_dir> -K -P --hmm /path/to/CANT-HYD.hmm

# Launch interactive dashboard (separate terminal)
Rscript 60_edna_mapping_viz.r
```

**Critical:** When using Kraken2 (`-K` flag), always use `24_process_reads_optimized.sh`. This script serializes Kraken2 calls to prevent memory exhaustion (Kraken2 loads 50-100GB into RAM). Other scripts will spawn multiple instances and crash the system.

### MAG Assembly

Complete assembly workflow:

```bash
cd 20_mag_assembly

# Automated pipeline: assembly → mapping → binning → polishing
./61_map_and_bin_optimized.sh

# Visualization
Rscript 80_plot_bins.R
Rscript 82_inter_binning_analysis.r
```

---

## Input/Output

**Input:**
Oxford Nanopore directory structure with multiplexed barcodes:
```
input/fastq_pass/
├── barcode01/*.fastq.gz
├── barcode02/*.fastq.gz
└── ...
```

**Output:**
```
output/
├── <flowcell>/
│   ├── <barcode>/
│   │   ├── fa/              Final quality-filtered sequences
│   │   ├── kraken/          Taxonomic classifications
│   │   ├── prokka/          Gene annotations
│   │   ├── hmm/             Functional gene matches
│   │   └── log.txt          Processing log
│   └── ...
├── assembly/                Assembled contigs
├── bins/                    MAG consensus bins
├── polished/                Refined genomes
├── checkm2/                 Quality metrics
└── <database>.duckdb        Integrated results
```

---

## Methodology

### Real-Time Processing Workflow

1. **Validation:** FASTQ integrity checking with automated repair (BBMap)
2. **Quality Control:** Adapter removal (BBDuk) and length/quality filtering (Filtlong Q7+, >1kb)
3. **Taxonomic Classification:** Kraken2 with custom marine databases
4. **Gene Annotation:** Prokka for rapid ORF prediction
5. **Functional Profiling:** HMMER search against ecosystem service gene databases
6. **Data Integration:** DuckDB for SQL-queryable results during expeditions

### MAG Assembly Workflow

1. **Assembly:** Flye metagenomic mode with overlap-based consensus
2. **Read Mapping:** minimap2 alignment for coverage calculation
3. **Binning:** Consensus approach using:
   - SemiBin2 (deep learning, trained on thousands of genomes)
   - MetaBAT2 (tetranucleotide frequency + differential coverage)
   - MaxBin2 (marker gene phylogeny)
   - DAS Tool (automated selection of best bins)
4. **Polishing:** Two rounds Racon + Medaka neural network correction
5. **Quality Assessment:** CheckM2 for completeness/contamination (MIMAG standards)
6. **Taxonomic Assignment:** Kaiju + GTDB-Tk classification

### Quality Standards

MAGs are classified according to MIMAG (Minimum Information about MAGs) standards:

| Quality | Completeness | Contamination | Criteria |
|---------|--------------|---------------|----------|
| High | >90% | <5% | + 23S, 16S, 5S rRNA & ≥18 tRNAs |
| Medium | >50% | <10% | - |
| Low | <50% | <10% | - |

---

## Dependencies

### Core Tools

**Sequencing & Assembly:**
Oxford Nanopore MinKNOW, Flye, minimap2

**Preprocessing:**
BBTools (BBDuk, BBMap), Filtlong

**Classification:**
Kraken2, Kaiju, GTDB-Tk

**Annotation:**
Prokka, HMMER3

**Binning:**
SemiBin2, MetaBAT2, MaxBin2, DAS Tool

**Quality Control:**
CheckM2, Racon, Medaka

### Analysis Environment

**Database:**
DuckDB (embedded OLAP)

**Statistics & Visualization:**
R (tidyverse, leaflet, plotly), Python 3.9+

**System:**
GNU parallel, trash-cli (safer file operations)

---

## Utility Scripts

**status.sh** — Dependency verification and version checking
**banner.sh** — Pipeline information display
**agents.sh** — Expert advisor system for domain-specific guidance

---

## Active Deployments

**CMO2025** — California-Mexico Oceanographic Survey
**QEI2025** — Queen Elizabeth Islands Arctic Expedition

*Note:* Scripts contain deployment-specific paths. Update before use in new projects.

---

## Documentation

**Repository Root:**
- `CLAUDE.md` — AI assistant instructions and architecture notes
- `CHANGELOG.md` — Version history and transformation log
- `EPIC_TRANSFORMATION.md` — Development narrative

**Subdirectories:**
- `10_realtime_processing/CLAUDE.md` — Detailed real-time pipeline architecture
- `10_realtime_processing/README.md` — Real-time processing guide
- `20_mag_assembly/README.md` — MAG assembly methodology

**Technical Notes:**
- `RESUME_LOGIC.md` — Stage-aware checkpoint system
- `HMM_SEARCH_GUIDE.md` — Functional gene profiling
- `CRASH_SAFETY.md` — Atomic operations and data integrity
- `DEPLOYMENT_ISSUES.md` — Portability and configuration

---

## Contributing

This is research software under active development. For bug reports or feature requests, please open an issue on GitHub.

**Development Guidelines:**
- Use `trash` instead of `rm` for all file operations
- Follow numbered naming convention (10s, 20s, 30s...)
- Include resume logic for all long-running operations
- Test on small datasets before expedition deployment
- Document hardcoded paths in DEPLOYMENT_ISSUES.md

---

## Citation

If this pipeline contributes to your research, please cite appropriately and consider contributing improvements back to the project.

**Key References:**

Prestat E, et al. (2014) FOAM: Functional Ontology Assignments for Metagenomes. *Nucleic Acids Research*.

Khot V, et al. (2022) CANT-HYD: A Curated Database of Phylogeny-Derived Hidden Markov Models for Annotation of Marker Genes Involved in Hydrocarbon Degradation.

Bowers RM, et al. (2017) Minimum information about a single amplified genome (MISAG) and a metagenome-assembled genome (MIMAG). *Nature Biotechnology* 35:725-731.

---

**Repository:** https://github.com/rec3141/danaSeq
**License:** [To be specified]
**Contact:** rec3141@gmail.com

---

*Decode the oceans, one read at a time.*
