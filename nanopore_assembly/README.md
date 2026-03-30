# Nanopore Assembly Pipeline

Preprocesses Oxford Nanopore reads and produces a co-assembly with depth tables and BAMs for downstream analysis by [mag_analysis](../mag_analysis/).

## Quick Start

```bash
cd nextflow
./install.sh && ./install.sh --check

# Local (conda)
./run-nanopore-assembly.sh --input /path/to/reads --outdir /path/to/output

# Apptainer (HPC)
./run-nanopore-assembly.sh --apptainer --input /path/to/reads --outdir /path/to/output

# Docker
./run-nanopore-assembly.sh --docker --input /path/to/reads --outdir /path/to/output
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
│   ├── tnf.tsv                 Tetranucleotide frequencies (136 features)
│   └── gc.tsv                  Per-contig GC%
├── mapping/
│   ├── *.sorted.bam            Per-sample alignments
│   ├── *.sorted.bam.bai        BAM indices
│   └── depths.txt              CoverM depth table (MetaBAT2 format)
└── pipeline_info/
```

Feed these into mag_analysis:
```bash
../mag_analysis/run-mag-analysis.sh \
    --assembly results/assembly/assembly.fasta \
    --depths results/mapping/depths.txt \
    --bam_dir results/mapping/ \
    --outdir /path/to/analysis --db_dir /path/to/databases
```

## Options

| Flag | Default | Description |
|------|---------|-------------|
| `--assembler` | `flye` | `flye`, `metamdbg`, or `myloasm` |
| `--polish` | auto | Flye polishing (auto: true for flye, false for others) |
| `--dedupe` | true | BBDuk deduplication before assembly |
| `--filtlong_size` | null | Filtlong target bases (e.g. `40000000000`) |
| `--run_remove_human` | true | Remove human reads via minimap2 |
| `--assembly_cpus` | 16 | CPUs for assembly |
| `--assembly_memory` | `60 GB` | Memory for assembly |

## Pipeline

```
Sample FASTQs → Concat per barcode → Dedupe + Filtlong → Remove human
    → Flye/metaMDBG/myloasm co-assembly → Polish (optional)
    → minimap2 mapping → CoverM depths
    → Tetranucleotide frequencies
```

10 processes, 1 conda environment (`dana-mag-assembly`).
