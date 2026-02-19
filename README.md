# dānaSeq

**Real-time metagenomic analysis for Oxford Nanopore sequencing on oceanographic expeditions.**

Named after the Buddhist concept of *dāna* (selfless giving), this pipeline processes DNA reads as they stream from the sequencer, providing live taxonomic classification, gene annotation, and functional profiling. Post-expedition, a separate pipeline reconstructs metagenome-assembled genomes (MAGs) from the collected data.

A companion Illumina short-read pipeline (METTA) handles metagenomic assembly with multi-assembler consensus and MetaBAT2 binning.

All pipelines are implemented in Nextflow DSL2 and run via conda or Docker with no hardcoded paths.

## Quick Start

### Real-time processing

```bash
git clone https://github.com/rec3141/danaSeq.git
cd danaSeq/nanopore_live/nextflow
./install.sh && ./install.sh --check

# Run (local conda, handles activation automatically)
./run-realtime.sh --input /path/to/nanopore/run --outdir /path/to/output \
    --run_kraken --kraken_db /path/to/krakendb \
    --run_prokka --run_sketch --run_tetra

# Or with Docker
docker build -t danaseq-realtime .
./run-realtime.sh --docker --input /path/to/nanopore/run --outdir /path/to/output \
    --run_kraken --kraken_db /path/to/krakendb \
    --run_prokka
```

### MAG assembly

```bash
cd danaSeq/nanopore_mag/nextflow
./install.sh && ./install.sh --check

# Run (local conda, handles activation automatically)
./run-mag.sh --input /path/to/reads --outdir /path/to/output

# Or with Docker
docker build -t danaseq-mag .
./run-mag.sh --docker --input /path/to/reads --outdir /path/to/output
```

### METTA assembly (Illumina short-read)

```bash
cd danaSeq/illumina_mag/nextflow
./install.sh && ./install.sh --check

# Run (local conda, handles activation automatically)
./run-metta.sh --input /path/to/reads --outdir /path/to/output

# Co-assembly mode
./run-metta.sh --input /path/to/reads --outdir /path/to/output --coassembly

# SLURM profile (Compute Canada)
./run-metta.sh --input /path/to/reads --outdir /path/to/output \
    -profile slurm --slurm_account def-myaccount \
    --conda_path ~/scratch/miniforge3/bin
```

### Kitchen sink — all modules (real-time)

```bash
cd danaSeq/nanopore_live/nextflow
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
cd danaSeq/nanopore_mag/nextflow
./run-mag.sh --input /data/reads --outdir /data/output \
    --dedupe \
    --filtlong_size 40000000000 \
    --min_overlap 1000 \
    --annotator prokka \
    --run_maxbin true \
    --run_lorbin true \
    --run_comebin true \
    --lorbin_min_length 80000 \
    --metabat_min_cls 50000 \
    --checkm2_db /path/to/checkm2_db \
    --genomad_db /path/to/genomad_db \
    --checkv_db /path/to/checkv_db \
    --kaiju_db /path/to/kaiju_db \
    --run_kraken2 true --kraken2_db /path/to/kraken2_db \
    --run_sendsketch true --sendsketch_address http://host:3068/sketch \
    --run_rrna true --silva_ssu_db /path/to/SILVA_SSU.fasta \
    --run_metabolism true \
    --kofam_db /path/to/kofam_db \
    --eggnog_db /path/to/eggnog_db \
    --dbcan_db /path/to/dbcan_db \
    --run_eukaryotic true --run_metaeuk true --metaeuk_db /path/to/metaeuk_db \
    --run_marferret true --marferret_db /path/to/marferret_db \
    --macsyfinder_models /path/to/macsyfinder_models \
    --defensefinder_models /path/to/defensefinder_models \
    --assembly_cpus 24 \
    --assembly_memory '64 GB'
```

### Kitchen sink — all options (METTA)

```bash
cd danaSeq/illumina_mag/nextflow
./run-metta.sh --input /data/reads --outdir /data/output \
    --coassembly \
    --min_readlen 70 \
    --run_normalize true \
    --run_tadpole true \
    --run_megahit true \
    --run_spades true \
    --run_metaspades true \
    --dedupe_identity 98 \
    --metabat_min_cls 2000 \
    --store_dir /scratch/metta_store \
    --assembly_cpus 24 \
    --assembly_memory '250 GB'
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
│   │   ├── Dockerfile           Self-contained Docker image
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

### Real-time processing (nanopore_live/)

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

See `nanopore_live/README.md` for full details.

### MAG assembly (nanopore_mag/)

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

See `nanopore_mag/README.md` for full details.

### METTA assembly (illumina_mag/)

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

See `illumina_mag/README.md` for full details.

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
