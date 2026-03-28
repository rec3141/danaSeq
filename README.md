# dānaSeq

**Real-time metagenomic analysis for Oxford Nanopore sequencing on oceanographic expeditions.**

Named after the Buddhist concept of *dāna* (selfless giving), this pipeline processes DNA reads as they stream from the sequencer, providing live taxonomic classification, gene annotation, and functional profiling. Post-expedition, separate pipelines handle assembly and downstream MAG analysis.

The platform consists of four independent Nextflow DSL2 pipelines:
1. **nanopore_live** — Real-time streaming analysis during sequencing
2. **nanopore_assembly** — Nanopore assembly + mapping + depth calculation
3. **illumina_assembly** — Illumina multi-assembler consensus + mapping + depth calculation
4. **mag_analysis** — Technology-agnostic downstream: binning, annotation, taxonomy, metabolism, MGEs, eukaryotes, ecosystem services, viz

Assembly pipelines produce an assembly FASTA + depth table that feeds into `mag_analysis`.
All pipelines run via conda, Docker, or Apptainer with no hardcoded paths.

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

### 3. Download databases

Both pipelines require reference databases. Download them before first use:

```bash
# Human reference (required by both pipelines, ~4 GB)
./download-databases.sh --human

# Nanopore MAG pipeline databases
./download-databases.sh --genomad --checkv --checkm2 --kaiju

# Interactive menu (shows sizes and descriptions)
./download-databases.sh

# With Apptainer on HPC (no local conda needed, auto-pulls SIF if --sif omitted)
./download-databases.sh --apptainer --sif /path/to/danaseq-mag.sif \
    --all --dir /path/to/databases

# Or with Docker
./download-databases.sh --docker --genomad --checkv --checkm2 --kaiju

# Custom download location (default: ./databases)
./download-databases.sh --dir /path/to/databases --genomad --checkv
```

See [`mag_analysis/nextflow/download-databases.sh`](mag_analysis/nextflow/download-databases.sh) for the full database list with sizes.

## Quick Start

### Real-time processing ([details](nanopore_live/README.md))

```bash
cd danaSeq/nanopore_live/nextflow
./install.sh && ./install.sh --check

./run-realtime.sh --input /path/to/nanopore/run --outdir /path/to/output \
    --run_kraken --kraken_db /path/to/krakendb \
    --run_prokka --run_sketch --run_tetra
```

### Nanopore assembly

```bash
cd danaSeq/nanopore_assembly/nextflow
./install.sh && ./install.sh --check

./run-nanopore-assembly.sh --input /path/to/reads --outdir /path/to/output
```

### Illumina assembly

```bash
cd danaSeq/illumina_assembly/nextflow
./install.sh && ./install.sh --check

./run-illumina-assembly.sh --input /path/to/reads --outdir /path/to/output

# SLURM profile (Compute Canada)
./run-illumina-assembly.sh --input /path/to/reads --outdir /path/to/output \
    -profile slurm --slurm_account def-myaccount
```

### MAG analysis (downstream of any assembly)

```bash
cd danaSeq/mag_analysis/nextflow
./install.sh && ./install.sh --check

# Using outputs from nanopore_assembly
./run-mag-analysis.sh \
    --assembly /path/to/assembly/assembly.fasta \
    --depths /path/to/mapping/depths.txt \
    --bam_dir /path/to/mapping/ \
    --outdir /path/to/output \
    --db_dir /path/to/databases \
    --all
```

### Test with bundled data

```bash
# Real-time pipeline
cd nanopore_live/nextflow
nextflow run main.nf --input test-data -profile test -resume
```

## Architecture

```
dānaSeq/
├── nanopore_live/          Real-time analysis during sequencing
│   └── nextflow/           11 processing stages → DuckDB
│
├── nanopore_assembly/      Nanopore assembly + mapping + depth
│   └── nextflow/           Flye/metaMDBG/myloasm → minimap2 → CoverM
│                           Output: assembly.fasta + depths.txt + BAMs
│
├── illumina_assembly/      Illumina multi-assembler + mapping + depth
│   └── nextflow/           4 assemblers → cascade dedupe → BBMap
│                           Output: assembly.fasta + depths.txt + BAMs
│
├── mag_analysis/           Technology-agnostic downstream analysis
│   └── nextflow/           Input: assembly + depths (from any assembler)
│       ├── modules/        binning, annotation, taxonomy, mge, eukaryotic,
│       │                    metabolism, rrna, phylogeny, viz, gene_depths
│       ├── viz/            Interactive Svelte dashboard
│       ├── ecossdb/        Ecosystem services (git submodule)
│       └── bin/            Pipeline scripts
│
├── archive/                Archived scripts and documentation
├── tests/                  Pipeline tests
├── CITATION.bib            References
└── LICENSE                 MIT
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

### Nanopore assembly (nanopore_assembly)

Preprocesses nanopore reads and produces a co-assembly:

```
Sample FASTQs → Concat per barcode → Dedupe + Filtlong → Remove human
                    → Flye/metaMDBG/myloasm co-assembly → Polish (optional)
                    → minimap2 mapping → CoverM depths
                    → Tetranucleotide frequencies

Output: assembly.fasta + depths.txt + BAMs + tnf.tsv
```

### Illumina assembly (illumina_assembly)

Processes Illumina paired-end reads through multi-assembler consensus:

```
Paired-end FASTQs → BBTools QC → Human decontamination → FastQC
                  → 3-phase error correction → bbnorm → bbmerge
                       │
            ┌──────────┼──────────┬──────────┐
         Tadpole    Megahit    SPAdes   metaSPAdes
            └──────────┼──────────┴──────────┘
                  Cascade deduplication (100% → 99% → 98%)
                  → bbmap mapping → jgi_summarize_bam_contig_depths

Output: assembly.fasta + depths.txt + BAMs
```

### MAG analysis (mag_analysis)

Technology-agnostic downstream analysis — accepts assembly + depths from any source:

```
assembly.fasta + depths.txt + BAMs (optional)
         │
    ┌────┼────┬────┬────┬────┬────┐
 SemiBin2 MetaBAT2 MaxBin2 LorBin COMEBin VAMB
    └────┼────┴────┴────┴────┴────┘
   DAS Tool / Binette / MAGScoT consensus → CheckM2 → GTDB-Tk
         │
    Parallel annotation & classification:
         ├── Prokka/Bakta → KofamScan + eggNOG + dbCAN → KEGG modules
         ├── Kaiju, Kraken2, sendsketch, rRNA (SILVA)
         ├── geNomad → CheckV, IntegronFinder, IslandPath, MacSyFinder, DefenseFinder
         ├── Tiara + Whokaryote → MetaEuk → MarFERReT
         ├── ECOSSDB ecosystem services
         └── Interactive Svelte dashboard (--run_viz)
```

Key features:
- **Seven-binner consensus**: SemiBin2, MetaBAT2, MaxBin2, LorBin, COMEBin, VAMB, VAMB-tax → DAS Tool + Binette + MAGScoT
- **Four taxonomy classifiers**: Kaiju (protein), Kraken2 (k-mer), sendsketch (GTDB MinHash), rRNA (SILVA)
- **Metabolic profiling**: KofamScan + eggNOG-mapper + dbCAN → KEGG modules, MinPath, KEGG-Decoder
- **Mobile genetic elements**: geNomad, CheckV, IntegronFinder, IslandPath, MacSyFinder, DefenseFinder
- **Eukaryotic analysis**: Tiara + Whokaryote classification, MetaEuk gene prediction, MarFERReT
- **Ecosystem services**: ECOSSDB mapping to CICES 5.2 + UN SDG targets
- **Interactive dashboard**: Svelte + Plotly + Phylocanvas.gl (`--run_viz`)
- **Works with any assembler**: Nanopore, Illumina, HiFi, or external assemblies

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

Illumina paired-end reads (illumina_mag pipeline):
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
