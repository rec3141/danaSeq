# Nanopore Assembly Pipeline

Preprocesses Oxford Nanopore reads and produces a co-assembly with depth tables and BAMs for downstream analysis by [mag_analysis](mag-analysis.md).

## Quick Start

```bash
cd nanopore_assembly/nextflow
./install.sh && ./install.sh --check

# Local (conda)
./run-nanopore-assembly.sh --input /path/to/reads --outdir /path/to/output

# Apptainer (HPC)
./run-nanopore-assembly.sh --apptainer --input /path/to/reads --outdir /path/to/output

# Docker
./run-nanopore-assembly.sh --docker --input /path/to/reads --outdir /path/to/output
```

## Pipeline Overview

```
Sample FASTQs (N barcodes)
         |
   CONCAT_BARCODES         Concatenate per-barcode FASTQs
         |
   DEDUPE_READS            Deduplicate reads (optional)
         |
   FILTLONG                Size-select and quality-filter
         |
   REMOVE_HUMAN            Remove human reads (optional)
         |
   ASSEMBLY                Co-assembly (Flye --meta / metaMDBG / myloasm)
         |
   POLISH                  Assembly polishing (optional)
         |
   MAP_READS (xN)          Per-sample alignment (minimap2, -F 0x904)
         |
   CALCULATE_DEPTHS        Coverage depth table (CoverM)
         |
   TETRAMER_FREQ           Tetranucleotide frequencies

Output: assembly.fasta + depths.txt + BAMs + tnf.tsv
```

## Input

Nanopore barcode directory structure (MinKNOW output):

```
input_dir/
├── fastq_pass/
│   ├── barcode01/*.fastq.gz
│   ├── barcode02/*.fastq.gz
│   └── ...
```

Or a flat directory of FASTQ files:

```
input_dir/
├── sample1.fastq.gz
├── sample2.fastq.gz
└── ...
```

## Output

```
results/
├── assembly/
│   ├── assembly.fasta          Co-assembly
│   ├── assembly_info.txt       Contig metadata (length, coverage, circularity)
│   ├── assembly_graph.gfa      Assembly graph
│   ├── tnf.tsv                 Tetranucleotide frequencies
│   └── gc.tsv                  Per-contig GC%
├── mapping/
│   ├── <sample>.sorted.bam     Per-sample alignments
│   ├── <sample>.sorted.bam.bai BAM indices
│   └── depths.txt              CoverM depth matrix (MetaBAT2 format)
└── pipeline_info/              Nextflow reports
```

## Key Features

- **Multiple assemblers**: Flye, metaMDBG, myloasm (configurable)
- **CoverM depth calculation**: Avoids MetaBAT2 integer overflow with supplementary alignments
- **Human decontamination**: Optional removal via BBTools
- **Tetranucleotide frequencies**: 136-feature TNF profiles for downstream binning
- Assembly + depths output feeds directly into `mag_analysis`
