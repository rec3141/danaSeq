# dānaSeq

**Real-time metagenomic analysis for Oxford Nanopore sequencing on oceanographic expeditions.**

Named after the Buddhist concept of *dāna* (selfless giving), this pipeline processes DNA reads as they stream from the sequencer, providing live taxonomic classification, gene annotation, and functional profiling. Post-expedition, a separate pipeline reconstructs metagenome-assembled genomes (MAGs) from the collected data.

Both pipelines are implemented in Nextflow DSL2 and run via conda or Docker with no hardcoded paths.

## Quick Start

### Real-time processing

```bash
git clone https://github.com/rec3141/danaSeq.git
cd danaSeq/10_realtime_processing/nextflow

# Option A: Conda
./install.sh && ./install.sh --check
conda activate conda-envs/dana-tools
nextflow run main.nf --input /path/to/nanopore/run \
    --run_kraken --kraken_db /path/to/krakendb \
    --run_prokka --run_sketch --run_tetra \
    -resume

# Option B: Docker
docker build -t danaseq-realtime .
./run-realtime.sh --docker --input /path/to/nanopore/run --outdir /path/to/output \
    --run_kraken --kraken_db /path/to/krakendb \
    --run_prokka
```

### MAG assembly (post-expedition)

```bash
cd danaSeq/20_mag_assembly/nextflow

# Option A: Conda
./install.sh && ./install.sh --check
mamba run -p conda-envs/dana-mag-flye \
    nextflow run main.nf --input /path/to/reads -resume

# Option B: Docker
docker build -t danaseq-mag .
./run-mag.sh --docker --input /path/to/reads --outdir /path/to/output
```

### Kitchen sink — all modules (real-time)

```bash
cd danaSeq/10_realtime_processing/nextflow
./run-realtime.sh --input /data/run1 --outdir /data/output \
    --run_kraken --kraken_db /path/to/krakendb \
    --run_prokka \
    --hmm_databases /path/to/CANT-HYD.hmm,/path/to/FOAM.hmm \
    --run_sketch \
    --run_tetra \
    --run_db_integration \
    --cleanup \
    --min_readlen 1500 \
    --keep_percent 80 \
    --min_file_size 1000000
```

### Kitchen sink — all options (MAG)

```bash
cd danaSeq/20_mag_assembly/nextflow
./run-mag.sh --input /data/reads --outdir /data/output \
    --dedupe \
    --filtlong_size 40000000000 \
    --min_overlap 1000 \
    --run_maxbin true \
    --metabat_min_cls 50000 \
    --assembly_cpus 24 \
    --assembly_memory '64 GB'
```

### Test with bundled data

```bash
# Real-time pipeline
cd 10_realtime_processing/nextflow
nextflow run main.nf --input test-data -profile test -resume

# MAG pipeline
cd 20_mag_assembly/nextflow
mamba run -p conda-envs/dana-mag-flye \
    nextflow run main.nf --input test-data -profile test -resume
```

## Architecture

```
dānaSeq/
├── 10_realtime_processing/       Real-time analysis during sequencing
│   ├── nextflow/                 Nextflow DSL2 pipeline
│   │   ├── main.nf              Entry point (11 processing stages)
│   │   ├── modules/             validate, qc, kraken, prokka, hmm, sketch, tetramer, db
│   │   ├── bin/                 AWK parsers, R/DuckDB integration scripts
│   │   ├── envs/                Conda YAML specs (3 environments)
│   │   ├── Dockerfile           Self-contained Docker image
│   │   └── test-data/           Bundled test data
│   └── archive/                 Legacy bash scripts (reference only)
│
├── 20_mag_assembly/              Post-expedition genome reconstruction
│   ├── nextflow/                 Nextflow DSL2 pipeline
│   │   ├── main.nf              Entry point (7 processes)
│   │   ├── modules/             assembly, mapping, binning
│   │   ├── envs/                Conda YAML specs (5 environments)
│   │   ├── Dockerfile           Self-contained Docker image
│   │   └── test-data/           Bundled test data
│   ├── 40s-90s scripts          Polishing, taxonomy, visualization (not yet ported)
│   └── archive/                 Replaced bash scripts (reference only)
│
├── 30_archive/                   Archived root-level scripts and documentation
├── tests/                        Pipeline tests
├── CITATION.bib                  References
└── LICENSE                       MIT
```

## Pipelines

### Real-time processing (10_realtime_processing/)

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

See `10_realtime_processing/README.md` for full details.

### MAG assembly (20_mag_assembly/)

Reconstructs metagenome-assembled genomes from collected reads:

```
Sample FASTQs → Flye co-assembly → minimap2 mapping → CoverM depths
                                                          ├── SemiBin2
                                                          ├── MetaBAT2
                                                          └── MaxBin2
                                                                ↓
                                                          DAS Tool consensus
```

Key features:
- CoverM for depth calculation (avoids MetaBAT2 integer overflow bug)
- Supplementary alignment filtering (`-F 0x904`) for long reads
- Dynamic binner architecture (add new binners with one line)
- Graceful failure handling for small/empty datasets
- GPU-accelerated SemiBin2 (local), CPU-only in Docker

See `20_mag_assembly/README.md` for full details.

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

- Flye: Kolmogorov et al., Nature Biotechnology 2019
- SemiBin2: Pan et al., Nature Communications 2023
- MetaBAT2: Kang et al., PeerJ 2019
- DAS Tool: Sieber et al., Nature Microbiology 2018
- CoverM: [github.com/wwood/CoverM](https://github.com/wwood/CoverM)
- MIMAG: Bowers et al., Nature Biotechnology 2017
- FOAM: Prestat et al., Nucleic Acids Research 2014
- CANT-HYD: Khot et al., 2022

See `CITATION.bib` for the complete reference list.

## License

MIT. See [LICENSE](LICENSE).

**Repository:** https://github.com/rec3141/danaSeq
**Contact:** rec3141@gmail.com
