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
dānaSeq/
├── 10_realtime_processing/        Real-time analysis during sequencing
│   ├── nextflow/                  Primary pipeline (Nextflow DSL2)
│   │   ├── main.nf               Pipeline entry point
│   │   ├── modules/              Process definitions (8 modules)
│   │   ├── bin/                  Helper scripts (AWK, R)
│   │   ├── conda-envs/           Pre-built conda environments
│   │   ├── envs/                 Conda YAML specs
│   │   ├── Dockerfile            Self-contained Docker image
│   │   ├── run-docker.sh         Docker convenience wrapper
│   │   ├── install.sh            Conda env installer
│   │   └── test-data/            Bundled test data
│   └── archive/                  Legacy bash scripts (reference only)
│
├── 20_mag_assembly/               Post-expedition genome assembly
│   ├── 10s  Assembly              Flye metagenomic assembly
│   ├── 20s  Mapping               Read alignment and coverage
│   ├── 30s  Binning               Multi-tool consensus binning
│   ├── 40s  Polishing             Racon + Medaka refinement
│   ├── 50s  Characterization      Taxonomy and quality assessment
│   ├── 60s  Pipelines             End-to-end workflows
│   ├── 70s  Format conversion     Interoperability
│   ├── 80s  Visualization         Ordination and clustering
│   └── 90s  Integration           Ecosystem service prediction
│
├── 30_archive/                    Archived scripts and documentation
└── tests/                         Pipeline tests
```

---

## Quick Start

### Option 1: Docker (recommended for new users)

No dependencies required beyond Docker. Build once, runs anywhere:

```bash
git clone https://github.com/rec3141/danaSeq.git
cd danaSeq/10_realtime_processing/nextflow

# Build the container (~15-20 min first time)
docker build -t danaseq-realtime .

# Quick test with bundled test data
./run-docker.sh --input test-data --outdir /tmp/test-output -profile test

# Run on your own data with Kraken2
./run-docker.sh --input /path/to/nanopore/run --outdir /path/to/output \
    --run_kraken --kraken_db /path/to/krakendb \
    --run_prokka --run_sketch --run_tetra

# With HMM functional gene profiling (e.g. CANT-HYD hydrocarbon degradation)
./run-docker.sh --input /path/to/nanopore/run --outdir /path/to/output \
    --run_kraken --kraken_db /path/to/krakendb \
    --run_prokka --hmm_databases /path/to/CANT-HYD.hmm \
    --run_sketch --run_tetra
```

The `run-docker.sh` helper validates all paths, creates the output directory,
sets up volume mounts, and runs as your user automatically. Run `./run-docker.sh --help`
for all options, or `./run-docker.sh --help-pipeline` for Nextflow flags.

<details>
<summary>Manual docker run (without helper script)</summary>

```bash
mkdir -p /path/to/output
docker run --user $(id -u):$(id -g) \
    -v /path/to/nanopore/run:/data/input:ro \
    -v /path/to/output:/data/output \
    -v /path/to/krakendb:/kraken_db:ro \
    danaseq-realtime \
    run /pipeline/main.nf \
        --input /data/input --outdir /data/output \
        --run_kraken --kraken_db /kraken_db \
        --run_prokka --run_sketch --run_tetra \
        -resume
```

**Notes for manual usage:**
- `--input` must point to the directory **containing** `fastq_pass/`, not `fastq_pass/` itself
- Always use `--user $(id -u):$(id -g)` so output files are owned by your user
- Create the output directory on the host first (`mkdir -p`) to avoid root-owned dirs
- Mount HMM files individually: `-v /path/to/A.hmm:/hmm/A.hmm:ro`

</details>

### Option 2: Conda

Install the pipeline tools into isolated conda environments:

```bash
git clone https://github.com/rec3141/danaSeq.git
cd danaSeq/10_realtime_processing/nextflow

# Install conda environments (~15-20 min first time)
./install.sh

# Verify installation
./install.sh --check

# Activate and run
conda activate conda-envs/dana-tools
nextflow run main.nf --input /path/to/data \
    --run_kraken --kraken_db /path/to/krakendb \
    --run_prokka --hmm_databases /path/to/CANT-HYD.hmm \
    --run_sketch --run_tetra \
    -resume
```

Requires: [Miniforge](https://github.com/conda-forge/miniforge) (conda/mamba) on PATH.

### Watch Mode (live sequencing)

Monitor a directory for new FASTQ files during active sequencing:

```bash
nextflow run main.nf --input /path/to/runs \
    --watch --db_sync_minutes 10 \
    --run_kraken --kraken_db /path/to/db \
    --run_prokka \
    --run_db_integration --danadir /path/to/r_scripts
```

In watch mode, `DB_SYNC` runs as a long-lived process that periodically scans
output directories and loads new results into DuckDB.

### MAG Assembly (post-expedition)

```bash
cd 20_mag_assembly

# Automated pipeline: assembly -> mapping -> binning -> polishing
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
results/
├── <flowcell>/
│   ├── <barcode>/
│   │   ├── fa/              Final quality-filtered sequences
│   │   ├── kraken/          Taxonomic classifications
│   │   ├── prokka/          Gene annotations
│   │   ├── hmm/             Functional gene matches
│   │   └── log.txt          Processing log
│   └── ...
├── dana.duckdb             Integrated results database
├── assembly/                Assembled contigs
├── bins/                    MAG consensus bins
├── polished/                Refined genomes
└── checkm2/                 Quality metrics
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
| High | >90% | <5% | + 23S, 16S, 5S rRNA & >=18 tRNAs |
| Medium | >50% | <10% | - |
| Low | <50% | <10% | - |

---

## Dependencies

All dependencies are managed automatically via conda environments or Docker. Manual installation is not required for most users.

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

---

## Active Deployments

**CMO2025** -- Churchill Marine Observatory Mesocosms
**QEI2025** -- Queen Elizabeth Islands Arctic Expedition

*Note:* Some scripts contain deployment-specific paths. Update before use in new projects.

---

## Documentation

- `10_realtime_processing/CLAUDE.md` -- Pipeline architecture and development guide
- `10_realtime_processing/README.md` -- Real-time processing details
- `20_mag_assembly/README.md` -- MAG assembly methodology
- `CITATION.bib` -- Reference database

Historical documentation (archived in `30_archive/`):
- `METHODS.md` -- Scientific methodology
- `CONTRIBUTING.md` -- Development guidelines
- `CHANGELOG.md` -- Version history
- `SECURITY_AUDIT.md` -- Security review notes

---

## Contributing

This is research software under active development. For bug reports or feature requests, please open an issue on GitHub.

**Development Guidelines:**
- Use `trash` instead of `rm` for all file operations
- New real-time analysis steps should be Nextflow modules (see `10_realtime_processing/CLAUDE.md`)
- Include resume logic for all long-running operations
- Test on small datasets before expedition deployment

---

## Citation

If this pipeline contributes to your research, please cite appropriately and consider contributing improvements back to the project.

**Key References:**

Prestat E, et al. (2014) FOAM: Functional Ontology Assignments for Metagenomes. *Nucleic Acids Research*.

Khot V, et al. (2022) CANT-HYD: A Curated Database of Phylogeny-Derived Hidden Markov Models for Annotation of Marker Genes Involved in Hydrocarbon Degradation.

Bowers RM, et al. (2017) Minimum information about a single amplified genome (MISAG) and a metagenome-assembled genome (MIMAG). *Nature Biotechnology* 35:725-731.

---

**Repository:** https://github.com/rec3141/danaSeq
**License:** MIT (see LICENSE)
**Contact:** rec3141@gmail.com

---

*Decode the oceans, one read at a time.*
