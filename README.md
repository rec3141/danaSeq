# danaSeq

[![Docs](https://img.shields.io/badge/docs-GitHub%20Pages-blue)](https://rec3141.github.io/danaSeq/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

Metagenomic analysis pipelines for Oxford Nanopore and Illumina sequencing data. Four Nextflow DSL2 pipelines cover real-time read processing, long-read assembly, short-read assembly, and downstream MAG analysis.

**[Full documentation](https://rec3141.github.io/danaSeq/)**

## Pipelines

| Pipeline | Purpose | Key tools |
|----------|---------|-----------|
| **[nanopore_live](nanopore_live/)** | Real-time analysis during sequencing | Kraken2, Prokka/Bakta, HMMER3, DuckDB |
| **[nanopore_assembly](nanopore_assembly/)** | Long-read assembly + mapping + depth | Flye, metaMDBG, myloasm, minimap2, CoverM |
| **[illumina_assembly](illumina_assembly/)** | Multi-assembler consensus assembly | Tadpole, Megahit, SPAdes, metaSPAdes, BBMap |
| **[mag_analysis](mag_analysis/)** | Downstream MAG analysis | 7-binner consensus, DAS Tool, Binette, 50+ processes |

## Quick Start

```bash
git clone https://github.com/rec3141/danaSeq.git
cd danaSeq

# Real-time processing
cd nanopore_live
./install.sh && ./install.sh --check
./run-realtime.sh --input /path/to/nanopore/run --outdir /path/to/output \
    --run_kraken --kraken_db /path/to/krakendb \
    --run_prokka --run_sketch --run_tetra

# MAG analysis
cd ../mag_analysis
./install.sh && ./install.sh --check
./run-mag-analysis.sh \
    --assembly /path/to/assembly.fasta \
    --depths /path/to/depths.txt \
    --outdir /path/to/output \
    --db_dir /path/to/databases --all
```

## Databases

```bash
./download-databases.sh              # Interactive menu
./download-databases.sh --human      # Human reference (~4 GB)
./download-databases.sh --genomad --checkv --checkm2 --kaiju
```

## License

MIT. See [LICENSE](LICENSE).

**Repository:** <https://github.com/rec3141/danaSeq>
**Documentation:** <https://rec3141.github.io/danaSeq/>
