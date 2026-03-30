# danaSeq

Metagenomic analysis pipelines for Oxford Nanopore and Illumina sequencing data. Four independent Nextflow DSL2 pipelines cover real-time read classification, long-read and short-read assembly, and downstream MAG analysis including binning, annotation, taxonomy, and metabolic profiling. Assembly pipelines produce a FASTA + depth table that feeds directly into `mag_analysis`.

## Architecture

```
danaSeq/
├── nanopore_live/          Real-time analysis during sequencing
│   │                       9 modules, 14 processes -> DuckDB
│
├── nanopore_assembly/      Long-read assembly + mapping + depth
│   │                       3 modules: preprocess, assembly, mapping
│   │                       Flye/metaMDBG/myloasm -> minimap2 -> CoverM
│   │                       Output: assembly.fasta + depths.txt + BAMs
│
├── illumina_assembly/      Multi-assembler consensus + mapping + depth
│   │                       7 modules: preprocess, error_correct, normalize,
│   │                       merge_reads, assembly, dedupe, mapping
│   │                       Output: assembly.fasta + depths.txt + BAMs
│
├── mag_analysis/           Technology-agnostic downstream analysis
│   │                       10 modules: binning, annotation, taxonomy, rrna,
│   │                       metabolism, mge, eukaryotic, gene_depths,
│   │                       phylogeny, viz
│   ├── viz/                Interactive Svelte dashboard
│   └── bin/                Pipeline scripts
│
└── tests/                  Pipeline tests
```

## Pipelines

| Pipeline | Purpose | Key tools |
|----------|---------|-----------|
| [**nanopore_live**](nanopore-live.md) | Real-time analysis during sequencing | Kraken2, Prokka/Bakta, HMMER3, DuckDB |
| [**nanopore_assembly**](nanopore-assembly.md) | Long-read assembly + mapping + depth | Flye, metaMDBG, myloasm, minimap2, CoverM |
| [**illumina_assembly**](illumina-assembly.md) | Multi-assembler consensus assembly | Tadpole, Megahit, SPAdes, metaSPAdes, BBMap |
| [**mag_analysis**](mag-analysis.md) | Technology-agnostic downstream analysis | 7-binner consensus, DAS Tool, Binette, 50+ processes |

## Getting Started

### 1. Clone the repository

```bash
git clone https://github.com/rec3141/danaSeq.git
cd danaSeq
```

### 2. Choose a runtime

| Method | Best for | Setup |
|--------|----------|-------|
| **Conda** | Local/laptop development | `cd <pipeline> && ./install.sh` |
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
| Full analysis | 32 | 256 GB | 2 TB |

## License

MIT. See [LICENSE](https://github.com/rec3141/danaSeq/blob/main/LICENSE).
