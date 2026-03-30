# danaSeq

[![Docs](https://img.shields.io/badge/docs-GitHub%20Pages-blue)](https://rec3141.github.io/danaSeq/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

**Real-time metagenomic analysis for Oxford Nanopore sequencing on oceanographic expeditions.**

Named after the Buddhist concept of *dana* (selfless giving), danaSeq processes DNA reads as they stream from the sequencer, providing live taxonomic classification, gene annotation, and functional profiling. Post-expedition, separate pipelines handle assembly and downstream MAG analysis.

The platform consists of four independent Nextflow DSL2 pipelines. Assembly pipelines produce an assembly FASTA + depth table that feeds into `mag_analysis`. All pipelines run via conda, Docker, or Apptainer with no hardcoded paths.

**[Full documentation](https://rec3141.github.io/danaSeq/)**

## Pipelines

| Pipeline | Purpose | Key tools |
|----------|---------|-----------|
| **[nanopore_live](nanopore_live/README.md)** | Real-time analysis during sequencing | Kraken2, Prokka, HMMER3, DuckDB |
| **[nanopore_assembly](nanopore_assembly/README.md)** | Nanopore assembly + mapping + depth | Flye, metaMDBG, minimap2, CoverM |
| **[illumina_assembly](illumina_assembly/README.md)** | Illumina multi-assembler consensus | Tadpole, Megahit, SPAdes, metaSPAdes, BBMap |
| **[mag_analysis](mag_analysis/README.md)** | Technology-agnostic downstream analysis | 7-binner consensus, DAS Tool, Binette, 50+ processes |

## Quick Start

```bash
# Clone
git clone https://github.com/rec3141/danaSeq.git
cd danaSeq

# Real-time processing
cd nanopore_live/nextflow
./install.sh && ./install.sh --check
./run-realtime.sh --input /path/to/nanopore/run --outdir /path/to/output \
    --run_kraken --kraken_db /path/to/krakendb \
    --run_prokka --run_sketch --run_tetra

# Nanopore assembly
cd ../../nanopore_assembly/nextflow
./install.sh && ./install.sh --check
./run-nanopore-assembly.sh --input /path/to/reads --outdir /path/to/output

# Illumina assembly
cd ../../illumina_assembly/nextflow
./install.sh && ./install.sh --check
./run-illumina-assembly.sh --input /path/to/reads --outdir /path/to/output

# MAG analysis (using assembly outputs)
cd ../../mag_analysis/nextflow
./install.sh && ./install.sh --check
./run-mag-analysis.sh \
    --assembly /path/to/assembly.fasta \
    --depths /path/to/depths.txt \
    --outdir /path/to/output \
    --db_dir /path/to/databases --all
```

## Runtime Options

| Method | Best for | Setup |
|--------|----------|-------|
| **Conda** | Local/laptop development | `cd <pipeline>/nextflow && ./install.sh` |
| **Docker** | Reproducible runs, CI | `docker pull ghcr.io/rec3141/danaseq-mag-analysis:latest` |
| **Apptainer** | HPC clusters | `./run-*.sh --apptainer` (auto-pulls SIF) |

## Databases

```bash
# Interactive menu (shows sizes and descriptions)
./download-databases.sh

# Human reference (~4 GB, required by both assembly pipelines)
./download-databases.sh --human

# MAG analysis databases
./download-databases.sh --genomad --checkv --checkm2 --kaiju
```

## Resource Requirements

| Configuration | CPUs | RAM | Storage |
|--------------|------|-----|---------|
| Minimum (no Kraken) | 16 | 32 GB | 500 GB |
| Recommended | 32 | 128 GB | 1 TB |
| Shipboard production | 32 | 256 GB | 2 TB |

## License

MIT. See [LICENSE](LICENSE).

**Repository:** <https://github.com/rec3141/danaSeq>
**Documentation:** <https://rec3141.github.io/danaSeq/>
