# danaSeq

**Real-time metagenomic analysis for Oxford Nanopore sequencing on oceanographic expeditions.**

Named after the Buddhist concept of *dana* (selfless giving), danaSeq processes DNA reads as they stream from the sequencer, providing live taxonomic classification, gene annotation, and functional profiling. Post-expedition, separate pipelines handle assembly and downstream MAG analysis.

The platform consists of four independent Nextflow DSL2 pipelines. Assembly pipelines produce an assembly FASTA + depth table that feeds into `mag_analysis`. All pipelines run via conda, Docker, or Apptainer with no hardcoded paths.

## Pipelines

| Pipeline | Purpose | Key tools |
|----------|---------|-----------|
| [**nanopore_live**](nanopore-live.md) | Real-time analysis during sequencing | Kraken2, Prokka, HMMER3, DuckDB |
| [**nanopore_assembly**](nanopore-assembly.md) | Nanopore assembly + mapping + depth | Flye, metaMDBG, minimap2, CoverM |
| [**illumina_assembly**](illumina-assembly.md) | Illumina multi-assembler consensus | Tadpole, Megahit, SPAdes, metaSPAdes, BBMap |
| [**mag_analysis**](mag-analysis.md) | Technology-agnostic downstream analysis | 7-binner consensus, DAS Tool, Binette, 50+ processes |

## Architecture

```
danaSeq/
├── nanopore_live/          Real-time analysis during sequencing
│   └── nextflow/           11 processing stages -> DuckDB
│
├── nanopore_assembly/      Nanopore assembly + mapping + depth
│   └── nextflow/           Flye/metaMDBG/myloasm -> minimap2 -> CoverM
│                           Output: assembly.fasta + depths.txt + BAMs
│
├── illumina_assembly/      Illumina multi-assembler + mapping + depth
│   └── nextflow/           4 assemblers -> cascade dedupe -> BBMap
│                           Output: assembly.fasta + depths.txt + BAMs
│
├── mag_analysis/           Technology-agnostic downstream analysis
│   └── nextflow/           Input: assembly + depths (from any assembler)
│       ├── modules/        binning, annotation, taxonomy, mge, eukaryotic,
│       │                    metabolism, rrna, phylogeny, viz, gene_depths
│       ├── viz/            Interactive Svelte dashboard
│       └── bin/            Pipeline scripts
│
└── tests/                  Pipeline tests
```

## Getting Started

### 1. Clone the repository

```bash
git clone https://github.com/rec3141/danaSeq.git
cd danaSeq
```

### 2. Choose a runtime

| Method | Best for | Setup |
|--------|----------|-------|
| **Conda** | Local/laptop development | `cd <pipeline>/nextflow && ./install.sh` |
| **Docker** | Reproducible runs, CI | `docker pull ghcr.io/rec3141/danaseq-mag-analysis:latest` |
| **Apptainer** | HPC clusters | `./run-*.sh --apptainer` (auto-pulls SIF) |

### 3. Download databases

```bash
# Interactive menu (shows sizes and descriptions)
./download-databases.sh

# Human reference (~4 GB, required by both assembly pipelines)
./download-databases.sh --human

# MAG analysis databases
./download-databases.sh --genomad --checkv --checkm2 --kaiju
```

## Quick Start

### Real-time processing

```bash
cd nanopore_live/nextflow
./install.sh && ./install.sh --check

./run-realtime.sh --input /path/to/nanopore/run --outdir /path/to/output \
    --run_kraken --kraken_db /path/to/krakendb \
    --run_prokka --run_sketch --run_tetra
```

### Nanopore assembly

```bash
cd nanopore_assembly/nextflow
./install.sh && ./install.sh --check

./run-nanopore-assembly.sh --input /path/to/reads --outdir /path/to/output
```

### Illumina assembly

```bash
cd illumina_assembly/nextflow
./install.sh && ./install.sh --check

./run-illumina-assembly.sh --input /path/to/reads --outdir /path/to/output
```

### MAG analysis

```bash
cd mag_analysis/nextflow
./install.sh && ./install.sh --check

./run-mag-analysis.sh \
    --assembly /path/to/assembly.fasta \
    --depths /path/to/depths.txt \
    --outdir /path/to/output \
    --db_dir /path/to/databases --all
```

## Quality Standards

MAGs are classified per MIMAG (Bowers et al. 2017):

| Quality | Completeness | Contamination | Additional |
|---------|-------------|---------------|------------|
| High | >90% | <5% | 23S, 16S, 5S rRNA + tRNAs |
| Medium | >50% | <10% | -- |
| Low | <50% | <10% | -- |

## Resource Requirements

| Configuration | CPUs | RAM | Storage |
|--------------|------|-----|---------|
| Minimum (no Kraken) | 16 | 32 GB | 500 GB |
| Recommended | 32 | 128 GB | 1 TB |
| Shipboard production | 32 | 256 GB | 2 TB |

## License

MIT. See [LICENSE](https://github.com/rec3141/danaSeq/blob/main/LICENSE).
