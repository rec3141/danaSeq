# MAG Analysis Pipeline

Technology-agnostic downstream analysis for metagenome-assembled genomes (MAGs). Accepts a pre-computed assembly + depth table from any assembler ([nanopore_assembly](../nanopore_assembly/), [illumina_assembly](../illumina_assembly/), or external) and runs binning, annotation, taxonomy, metabolic profiling, mobile genetic element detection, eukaryotic analysis, ecosystem services mapping, phylogenetics, and interactive visualization.

## Quick Start

```bash
cd nextflow
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
    --all --db_dir /path/to/databases \
    --sendsketch_address http://host:3068/sketch

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

## Databases

Download all databases:
```bash
./download-databases.sh --all --dir /path/to/databases
```

Or use `--db_dir` to auto-resolve paths from a standard layout:
```bash
./run-mag-analysis.sh --db_dir /path/to/databases ...
```

See `./download-databases.sh --help` for individual database options and sizes.

## Output

```
results/
├── assembly/                   Staged inputs (assembly_info, gc.tsv, depths)
├── annotation/                 Prokka or Bakta gene annotation
├── binning/
│   ├── metabat/                MetaBAT2 bins
│   ├── semibin/                SemiBin2 bins (optional)
│   ├── maxbin/                 MaxBin2 bins (optional)
│   ├── lorbin/                 LorBin bins (optional)
│   ├── comebin/                COMEBin bins (optional)
│   ├── vamb/                   VAMB bins (optional)
│   ├── dastool/                DAS Tool consensus bins
│   ├── binette/                Binette consensus (optional)
│   ├── magscot/                MAGScoT consensus (optional)
│   └── checkm2/               Bin quality assessment
├── taxonomy/
│   ├── kaiju/                  Protein-level taxonomy
│   ├── kraken2/                K-mer taxonomy
│   ├── sendsketch/             GTDB MinHash taxonomy
│   └── rrna/                   rRNA + tRNA gene classification
├── eukaryotic/
│   ├── tiara/                  Deep learning eukaryotic classification
│   └── whokaryote/             Gene structure classification
├── mge/
│   ├── genomad/                Virus + plasmid detection
│   ├── checkv/                 Viral quality assessment
│   ├── integrons/              IntegronFinder
│   ├── islandpath/             Genomic islands
│   ├── macsyfinder/            Secretion/conjugation systems
│   └── defensefinder/          Anti-phage defense systems
├── metabolism/
│   ├── kofamscan/              KEGG Orthology
│   ├── emapper/                eggNOG annotations
│   ├── dbcan/                  CAZyme annotations
│   ├── merged/                 Unified annotation table
│   ├── per_mag/                Per-bin annotation tables
│   ├── modules/                KEGG module completeness
│   ├── minpath/                MinPath pathway reconstruction
│   ├── kegg_decoder/           Biogeochemical function scoring
│   └── ecossdb/                Ecosystem services (CICES 5.2 + SDG)
├── mapping/
│   └── gene_depths.tsv         Per-gene depth profiles
├── viz/
│   ├── data/                   Dashboard JSON files
│   └── site/                   Static Svelte dashboard
└── pipeline_info/
```

## Modules

All modules are optional and controlled via `--run_*` flags. Use `--all` to enable everything.

| Category | Modules | Key flags |
|----------|---------|-----------|
| **Annotation** | Prokka, Bakta (+Extra) | `--annotator bakta` |
| **Binning** | SemiBin2, MetaBAT2, MaxBin2, LorBin, COMEBin, VAMB, VAMB-tax | `--run_semibin`, etc. |
| **Consensus** | DAS Tool (always), Binette, MAGScoT | `--run_binette`, `--run_magscot` |
| **Taxonomy** | Kaiju, Kraken2, Sendsketch, rRNA (SILVA) | `--run_kaiju`, etc. |
| **Eukaryotic** | Tiara, Whokaryote, MetaEuk, MarFERReT | `--run_eukaryotic` |
| **MGE** | geNomad, CheckV, IntegronFinder, IslandPath, MacSyFinder, DefenseFinder | `--run_genomad`, etc. |
| **Metabolism** | KofamScan, eggNOG, dbCAN, KEGG modules, MinPath, KEGG-Decoder, antiSMASH | `--run_metabolism` |
| **Ecosystem** | ECOSSDB (CICES 5.2 + UN SDG targets) | `--run_ecossdb` |
| **Quality** | CheckM2 | `--checkm2_db` |
| **Phylogenetics** | GTDB-Tk | `--run_gtdbtk` |
| **Visualization** | Interactive Svelte dashboard | `--run_viz` |

## Pipeline

```
assembly.fasta + depths.txt + BAMs (optional)
     │
     ├── Annotation (Prokka/Bakta)
     ├── Taxonomy (Kaiju, Kraken2, Sendsketch, rRNA)
     ├── Eukaryotic (Tiara, Whokaryote)
     ├── MGE (geNomad, CheckV, IntegronFinder, IslandPath, MacSyFinder, DefenseFinder)
     ├── Metabolism (KofamScan, eggNOG, dbCAN → merge → per-MAG → KEGG/MinPath/Decoder)
     │
     ├── Binning (7 binners → DAS Tool + Binette + MAGScoT → CheckM2)
     ├── Phylogenetics (GTDB-Tk)
     ├── Ecosystem services (ECOSSDB → CICES 5.2 → UN SDG)
     └── Visualization (4-stage incremental Svelte dashboard)
```

40+ processes, 14+ conda environments.
