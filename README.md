# dānaSeq

**Real-time metagenomic analysis for Oxford Nanopore sequencing on oceanographic expeditions.**

Named after the Buddhist concept of *dāna* (selfless giving), this pipeline processes DNA reads as they stream from the sequencer, providing live taxonomic classification, gene annotation, and functional profiling. Post-expedition, a separate pipeline reconstructs metagenome-assembled genomes (MAGs) from the collected data.

A companion Illumina short-read pipeline (METTA) handles metagenomic assembly with multi-assembler consensus and MetaBAT2 binning.

All pipelines are implemented in Nextflow DSL2 and run via conda, Docker, or Apptainer with no hardcoded paths.

## Getting Started

### 1. Get the code

```bash
# Clone the repository
git clone https://github.com/rec3141/danaSeq.git
cd danaSeq

# Or download a specific release
gh release download v0.1.0-alpha --archive=tar.gz
```

### 2. Choose a runtime

| Method | Best for | Setup |
|--------|----------|-------|
| **Conda** | Local/laptop development | `./install.sh && ./install.sh --check` |
| **Docker** | Reproducible runs, CI | `docker pull ghcr.io/rec3141/danaseq-mag:latest` |
| **Apptainer** | HPC clusters (auto-detected) | `apptainer pull danaseq-mag.sif docker://ghcr.io/rec3141/danaseq-mag:latest` |

### 3. Download databases (MAG pipeline)

Most MAG pipeline modules require reference databases. Download them before first use:

```bash
cd nanopore_mag/nextflow

# Interactive menu (shows sizes and descriptions)
./download-databases.sh

# Or download specific databases
./download-databases.sh --genomad --checkv --checkm2 --kaiju

# With Apptainer on HPC (no local conda needed, auto-pulls SIF)
./download-databases.sh --apptainer --genomad --checkv --checkm2 --kaiju

# Or with Docker
./download-databases.sh --docker --genomad --checkv --checkm2 --kaiju

# All databases at once (~150+ GB)
./download-databases.sh --all

# Custom download location (default: ./databases)
./download-databases.sh --dir /path/to/databases --genomad --checkv
```

See [`nanopore_mag/README.md`](nanopore_mag/README.md#databases) for the full database list with sizes.

## Quick Start

### Real-time processing ([details](nanopore_live/README.md))

```bash
cd danaSeq/nanopore_live/nextflow
./install.sh && ./install.sh --check

./run-realtime.sh --input /path/to/nanopore/run --outdir /path/to/output \
    --run_kraken --kraken_db /path/to/krakendb \
    --run_prokka --run_sketch --run_tetra
```

### MAG assembly ([details](nanopore_mag/README.md))

```bash
cd danaSeq/nanopore_mag/nextflow

# Local conda
./install.sh && ./install.sh --check
./run-mag.sh --input /path/to/reads --outdir /path/to/output

# Docker
./run-mag.sh --docker --input /path/to/reads --outdir /path/to/output

# Apptainer (HPC) — auto-pulls SIF on first run
./run-mag.sh --apptainer --input /path/to/reads --outdir /path/to/output \
    --db_dir /path/to/databases

# Or pre-pull the SIF and pass it explicitly
apptainer pull danaseq-mag.sif docker://ghcr.io/rec3141/danaseq-mag:latest
./run-mag.sh --apptainer --sif ./danaseq-mag.sif \
    --input /path/to/reads --outdir /path/to/output \
    --db_dir /path/to/databases
```

### METTA assembly — Illumina short-read ([details](illumina_mag/README.md))

```bash
cd danaSeq/illumina_mag/nextflow
./install.sh && ./install.sh --check

./run-metta.sh --input /path/to/reads --outdir /path/to/output

# SLURM profile (Compute Canada)
./run-metta.sh --input /path/to/reads --outdir /path/to/output \
    -profile slurm --slurm_account def-myaccount
```

### Test with bundled data

```bash
# Real-time pipeline
cd nanopore_live/nextflow
nextflow run main.nf --input test-data -profile test -resume

# MAG pipeline
cd nanopore_mag/nextflow
mamba run -p conda-envs/dana-mag-flye \
    nextflow run main.nf --input test-data -profile test -resume
```

## Architecture

```
dānaSeq/
├── nanopore_live/       Real-time analysis during sequencing
│   ├── nextflow/                 Nextflow DSL2 pipeline
│   │   ├── main.nf              Entry point (11 processing stages)
│   │   ├── modules/             validate, qc, kraken, prokka, hmm, sketch, tetramer, db
│   │   ├── bin/                 AWK parsers, R/DuckDB integration scripts
│   │   ├── envs/                Conda YAML specs (3 environments)
│   │   ├── Dockerfile           Self-contained Docker image
│   │   └── test-data/           Bundled test data
│   └── archive/                 Legacy bash scripts (reference only)
│
├── nanopore_mag/              Metagenome-assembled genome reconstruction
│   ├── nextflow/                 Nextflow DSL2 pipeline
│   │   ├── main.nf              Entry point (40+ processes)
│   │   ├── modules/             12 modules: assembly, mapping, binning, annotation,
│   │   │                          taxonomy, mge, metabolism, eukaryotic, rrna, refinement,
│   │   │                          preprocess, viz
│   │   ├── bin/                  Pipeline scripts (merge, map, MinPath, KEGG-Decoder, etc.)
│   │   ├── viz/                  Interactive dashboard (Svelte + Vite + Plotly)
│   │   ├── envs/                Conda YAML specs (27 environments)
│   │   ├── Dockerfile           Pipeline image (thin layer on base)
│   │   ├── Dockerfile.base      Base image (all conda envs + wrapper scripts)
│   │   └── test-data/           Bundled test data
│   └── archive/                 Replaced bash scripts (reference only)
│
├── illumina_mag/                  Illumina metagenomic assembly (METTA)
│   └── nextflow/                 Nextflow DSL2 pipeline
│       ├── main.nf              Entry point (18 processes)
│       ├── modules/             8 modules: preprocess, error_correct, normalize,
│       │                          merge_reads, assembly, dedupe, mapping, binning
│       ├── envs/                Conda YAML specs (4 environments)
│       └── run-metta.sh         Pipeline launcher
│
├── archive/                      Archived root-level scripts and documentation
├── tests/                        Pipeline tests
├── CITATION.bib                  References
└── LICENSE                       MIT
```

## Pipelines

### [Real-time processing](nanopore_live/README.md)

Processes FASTQ files as they arrive from Oxford Nanopore MinKNOW:

```
MinKNOW FASTQ → Validate → BBDuk QC → Filtlong → FASTA
                                                    ├── Kraken2 (taxonomy)
                                                    ├── Prokka (gene annotation)
                                                    ├── HMMER3 (functional genes)
                                                    ├── Sendsketch (profiling)
                                                    └── Tetranucleotide freq
                                                          ↓
                                                    DuckDB integration
```

Key features:
- **Watch mode** for live sequencing (`--watch`)
- Kraken2 serialized to one instance (50-100 GB database)
- DuckDB integration for real-time SQL queries
- HMM search against functional gene databases (FOAM, CANT-HYD, NCycDB, etc.)
- Post-DB cleanup to compress/delete source files after import

### [MAG assembly](nanopore_mag/README.md)

Reconstructs metagenome-assembled genomes alongside real-time processing:

```
Sample FASTQs → Flye co-assembly → minimap2 mapping → CoverM depths
                        │                                   │
                        │                            ┌──────┼──────┬──────┬──────┐
                        │                         SemiBin2 MetaBAT2 MaxBin2 LorBin COMEBin
                        │                            └──────┼──────┴──────┴──────┘
                        │                              DAS Tool consensus → CheckM2
                        │                                   │
                        │                              NCLB bin refinement (optional)
                        │                                   │
                        │                              VIZ_PREPROCESS (dashboard + static site)
                        │
                   Parallel annotation & classification:
                        ├── Prokka/Bakta gene annotation
                        │     ├── Kaiju protein-level taxonomy
                        │     ├── KofamScan + eggNOG-mapper + dbCAN → merge → bin mapping
                        │     │     ├── KEGG module completeness
                        │     │     ├── MinPath pathway reconstruction
                        │     │     └── KEGG-Decoder biogeochemical functions
                        │     ├── IslandPath-DIMOB genomic islands
                        │     ├── MacSyFinder secretion/conjugation systems
                        │     └── DefenseFinder anti-phage defense
                        ├── Kraken2 k-mer taxonomy
                        ├── sendsketch GTDB MinHash taxonomy
                        ├── barrnap + vsearch rRNA classification (SILVA)
                        ├── geNomad virus/plasmid → CheckV quality
                        ├── IntegronFinder integron detection
                        ├── Tiara + Whokaryote eukaryotic classification
                        │     └── MetaEuk eukaryotic gene prediction
                        │           └── MarFERReT marine eukaryotic taxonomy + Pfam
                        └── Tetranucleotide frequency profiles
```

Key features:
- **Five-binner consensus**: SemiBin2, MetaBAT2, MaxBin2, LorBin, COMEBin → DAS Tool
- **Four taxonomy classifiers**: Kaiju (protein), Kraken2 (k-mer), sendsketch (GTDB MinHash), rRNA (SILVA)
- **Metabolic profiling**: KofamScan + eggNOG-mapper + dbCAN → KEGG modules, MinPath, KEGG-Decoder
- **Mobile genetic elements**: geNomad, CheckV, IntegronFinder, IslandPath, MacSyFinder, DefenseFinder
- **Eukaryotic analysis**: Tiara + Whokaryote classification, MetaEuk gene prediction
- **NCLB bin refinement**: LLM-guided iterative bin refinement (optional)
- **Interactive dashboard**: Auto-generated Svelte viz with Plotly charts (`--run_viz`)
- CheckM2 quality assessment (completeness/contamination per MIMAG standards)
- CoverM for depth calculation (avoids MetaBAT2 integer overflow bug)
- GPU-accelerated ML binners (local), CPU-only in Docker

### [METTA assembly — Illumina short-read](illumina_mag/README.md)

Processes Illumina paired-end metagenomic reads through multi-assembler consensus:

```
Paired-end FASTQs → BBTools QC (clumpify, filterbytile, bbduk)
                          → 3-phase error correction (ecco, ecc, tadpole)
                          → bbnorm normalization → bbmerge read merging
                               │
                    ┌──────────┼──────────┬──────────┐
                 Tadpole    Megahit    SPAdes   metaSPAdes
                    └──────────┼──────────┴──────────┘
                          Cascade deduplication (100% → 99% → 98%)
                               │
                          bbmap read mapping → jgi_summarize_bam_contig_depths
                               │
                          MetaBAT2 binning
```

Key features:
- **Four-assembler consensus**: Tadpole, Megahit, SPAdes, metaSPAdes with cascade deduplication
- **BBTools-based QC**: Optical deduplication, tile filtering, adapter/artifact removal
- **Three-phase error correction**: Overlap, clump, and k-mer-based correction
- **Per-sample and co-assembly modes**: `--coassembly` pools all samples
- **SRA-safe**: Automatic fallbacks for SRA-stripped read headers
- **SLURM support**: `-profile slurm` for Compute Canada HPC clusters

## Input

Oxford Nanopore directory structure with multiplexed barcodes (real-time pipeline):
```
input_dir/fastq_pass/
├── barcode01/*.fastq.gz
├── barcode02/*.fastq.gz
└── ...
```

Directory of FASTQ files (MAG pipeline):
```
input_dir/
├── sample1.fastq.gz
├── sample2.fastq.gz
└── ...
```

Illumina paired-end reads (METTA pipeline):
```
input_dir/
├── sampleA_S1_L001_R1_001.fastq.gz
├── sampleA_S1_L001_R2_001.fastq.gz
└── ...
```

## Quality Standards

MAGs are classified per MIMAG (Bowers et al. 2017):

| Quality | Completeness | Contamination | Additional |
|---------|-------------|---------------|------------|
| High | >90% | <5% | 23S, 16S, 5S rRNA + tRNAs |
| Medium | >50% | <10% | -- |
| Low | <50% | <10% | -- |

## Resource Requirements

| Configuration | CPUs | RAM | Storage | Notes |
|--------------|------|-----|---------|-------|
| Minimum (no Kraken) | 16 | 32 GB | 500 GB | QC and annotation only |
| Recommended | 32 | 128 GB | 1 TB | Full pipeline with Kraken2 |
| Shipboard production | 32 | 256 GB | 2 TB | Real-time + MAG assembly |

## References

**Assembly & Binning:**
- BBTools/BBMap: Bushnell, [sourceforge.net/projects/bbmap](https://sourceforge.net/projects/bbmap/)
- Megahit: Li et al., Bioinformatics 2015
- SPAdes/metaSPAdes: Nurk et al., Genome Research 2017
- Flye: Kolmogorov et al., Nature Biotechnology 2019
- SemiBin2: Pan et al., Nature Communications 2023
- MetaBAT2: Kang et al., PeerJ 2019
- MaxBin2: Wu et al., Bioinformatics 2016
- LorBin: Gao et al., Briefings in Bioinformatics 2024
- COMEBin: Xie et al., Nature Communications 2024
- DAS Tool: Sieber et al., Nature Microbiology 2018
- CheckM2: Chklovski et al., Nature Methods 2023
- CoverM: [github.com/wwood/CoverM](https://github.com/wwood/CoverM)

**Annotation:**
- Prokka: Seemann, Bioinformatics 2014
- Bakta: Schwengers et al., Microbial Genomics 2021

**Taxonomy:**
- Kaiju: Menzel et al., Nature Communications 2016
- Kraken2: Wood et al., Genome Biology 2019
- BBSketch/sendsketch: Bushnell, [sourceforge.net/projects/bbmap](https://sourceforge.net/projects/bbmap/)
- barrnap: Seemann, [github.com/tseemann/barrnap](https://github.com/tseemann/barrnap)
- vsearch: Rognes et al., PeerJ 2016
- SILVA: Quast et al., Nucleic Acids Research 2013

**Mobile Genetic Elements:**
- geNomad: Camargo et al., Nature Biotechnology 2024
- CheckV: Nayfach et al., Nature Biotechnology 2021
- IntegronFinder: Cury et al., Nucleic Acids Research 2016
- IslandPath-DIMOB: Bertelli & Brinkman, Bioinformatics 2018
- MacSyFinder: Abby et al., PLoS ONE 2014
- DefenseFinder: Tesson et al., Nature Communications 2022

**Metabolic Profiling:**
- KofamScan: Aramaki et al., Bioinformatics 2020
- eggNOG-mapper: Cantalapiedra et al., Molecular Biology and Evolution 2021
- dbCAN3: Zheng et al., Nucleic Acids Research 2023
- MinPath: Ye & Doak, PLoS Computational Biology 2009
- KEGG-Decoder: Graham et al., bioRxiv 2018

**Eukaryotic Analysis:**
- Tiara: Karlicki et al., Bioinformatics 2022
- Whokaryote: Pronk et al., Microbial Genomics 2022
- MetaEuk: Levy Karin et al., Microbiome 2020
- MarFERReT: Carradec et al., Scientific Data 2023

**Standards:**
- MIMAG: Bowers et al., Nature Biotechnology 2017
- FOAM: Prestat et al., Nucleic Acids Research 2014
- CANT-HYD: Khot et al., 2022

See `CITATION.bib` for the complete reference list.

## License

MIT. See [LICENSE](LICENSE).

**Repository:** https://github.com/rec3141/danaSeq
**Contact:** rec3141@gmail.com
