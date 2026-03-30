# MAG Analysis Pipeline

Technology-agnostic downstream analysis for metagenome-assembled genomes (MAGs). Accepts a pre-computed assembly + depth table from any assembler ([nanopore_assembly](nanopore-assembly.md), [illumina_assembly](illumina-assembly.md), or external) and runs binning, annotation, taxonomy, metabolic profiling, mobile genetic element detection, eukaryotic analysis, ecosystem services mapping, phylogenetics, and interactive visualization.

## Quick Start

```bash
cd mag_analysis/nextflow
./install.sh && ./install.sh --check

# Basic run
./run-mag-analysis.sh \
    --assembly /path/to/assembly.fasta \
    --depths /path/to/depths.txt \
    --outdir /path/to/output \
    --annotator bakta --db_dir /path/to/databases

# Kitchen sink (all modules)
./run-mag-analysis.sh \
    --assembly /path/to/assembly.fasta \
    --depths /path/to/depths.txt \
    --bam_dir /path/to/mapping/ \
    --outdir /path/to/output \
    --all --db_dir /path/to/databases

# Apptainer (HPC)
./run-mag-analysis.sh --apptainer \
    --assembly /path/to/assembly.fasta \
    --depths /path/to/depths.txt \
    --outdir /path/to/output \
    --all --db_dir /path/to/databases
```

## Required Inputs

| Flag | Format | Source |
|------|--------|--------|
| `--assembly` | FASTA | Any assembler output |
| `--depths` | MetaBAT2-format TSV | CoverM or jgi_summarize_bam_contig_depths |
| `--bam_dir` (optional) | Dir with `*.sorted.bam` + `.bai` | For BAM-based binners (SemiBin2, LorBin, COMEBin) |

## Pipeline Overview

```
assembly.fasta + depths.txt + BAMs (optional)
         |
    +----+----+----+----+----+----+
 SemiBin2 MetaBAT2 MaxBin2 LorBin COMEBin VAMB
    +----+----+----+----+----+----+
   DAS Tool / Binette / MAGScoT consensus -> CheckM2 -> GTDB-Tk
         |
    Parallel annotation & classification:
         +-- Prokka/Bakta -> KofamScan + eggNOG + dbCAN -> KEGG modules
         +-- Kaiju, Kraken2, sendsketch, rRNA (SILVA)
         +-- geNomad -> CheckV, IntegronFinder, IslandPath, MacSyFinder, DefenseFinder
         +-- Tiara + Whokaryote -> MetaEuk -> MarFERReT
         +-- ECOSSDB ecosystem services
         +-- antiSMASH biosynthetic gene clusters
         +-- Interactive Svelte dashboard (--run_viz)
```

## Key Features

- **Seven-binner consensus**: SemiBin2, MetaBAT2, MaxBin2, LorBin, COMEBin, VAMB, VAMB-tax -> DAS Tool + Binette + MAGScoT
- **Four taxonomy classifiers**: Kaiju (protein), Kraken2 (k-mer), sendsketch (GTDB MinHash), rRNA (SILVA)
- **Metabolic profiling**: KofamScan + eggNOG-mapper + dbCAN -> KEGG modules, MinPath, KEGG-Decoder
- **Mobile genetic elements**: geNomad, CheckV, IntegronFinder, IslandPath, MacSyFinder, DefenseFinder
- **Eukaryotic analysis**: Tiara + Whokaryote classification, MetaEuk gene prediction, MarFERReT taxonomy
- **Ecosystem services**: ECOSSDB mapping to CICES 5.2 + UN SDG targets
- **Biosynthetic gene clusters**: antiSMASH detection and classification
- **Interactive dashboard**: Svelte + Plotly + Phylocanvas.gl (`--run_viz`)
- **Works with any assembler**: Nanopore, Illumina, HiFi, or external assemblies
- CheckM2 quality assessment (completeness/contamination per MIMAG standards)

## Databases

Download all databases:

```bash
./download-databases.sh --all --dir /path/to/databases
```

Or use `--db_dir` to auto-resolve paths from a standard layout.

## Parameters

### Assembly & Binning

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--assembly` | (required) | Assembly FASTA file |
| `--depths` | (required) | Depth matrix (MetaBAT2 format) |
| `--bam_dir` | (optional) | Directory with sorted BAMs for BAM-based binners |

### Annotation

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--annotator` | `prokka` | Gene annotator (`prokka` or `bakta`) |

### Taxonomy

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--kaiju_db` | (required) | Kaiju database path |
| `--run_kraken2` | `true` | Run Kraken2 classification |
| `--run_sendsketch` | `true` | Run sendsketch GTDB taxonomy |
| `--run_rrna` | `true` | Run rRNA classification (SILVA) |

### Mobile Genetic Elements

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--genomad_db` | (required) | geNomad database path |
| `--checkv_db` | (required) | CheckV database path |

### Metabolism

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--run_metabolism` | `true` | Enable metabolic profiling |
| `--kofam_db` | (required if metabolism) | KofamScan database path |
| `--eggnog_db` | (required if metabolism) | eggNOG database path |
| `--dbcan_db` | (required if metabolism) | dbCAN database path |

### Eukaryotic Analysis

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--run_eukaryotic` | `true` | Enable eukaryotic classification |
| `--run_metaeuk` | `true` | Run MetaEuk gene prediction |
| `--run_marferret` | `true` | Run MarFERReT taxonomy |

### Convenience

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--all` | `false` | Enable all analysis modules |
| `--db_dir` | (none) | Auto-resolve database paths from standard layout |
