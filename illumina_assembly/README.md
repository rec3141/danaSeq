# Illumina Assembly Pipeline

Processes Illumina paired-end metagenomic reads through multi-assembler consensus and produces assemblies with depth tables and BAMs for downstream analysis by [mag_analysis](../mag_analysis/).

## Quick Start

```bash
cd nextflow
./install.sh && ./install.sh --check

# Local (conda)
./run-illumina-assembly.sh --input /path/to/reads --outdir /path/to/output

# Apptainer (HPC)
./run-illumina-assembly.sh --apptainer --input /path/to/reads --outdir /path/to/output

# SLURM profile (Compute Canada)
./run-illumina-assembly.sh --input /path/to/reads --outdir /path/to/output \
    -profile slurm --slurm_account def-myaccount
```

## Input

Directory of paired-end Illumina reads (`*_R1_*.fastq.gz`, R2 auto-detected):
```
input_dir/
├── sampleA_S1_L001_R1_001.fastq.gz
├── sampleA_S1_L001_R2_001.fastq.gz
├── sampleB_S2_L001_R1_001.fastq.gz
├── sampleB_S2_L001_R2_001.fastq.gz
└── ...
```

## Output

```
results/
├── preprocess/<sample>/      QC'd reads + FastQC reports
├── error_correct/<sample>/   3-phase error-corrected reads
├── normalize/<sample>/       Coverage-normalized reads
├── merge/<sample>/           Merged + quality-trimmed reads
├── assembly/<sample>/
│   ├── <sample>.dedupe.fasta    Final deduplicated assembly
│   ├── <sample>.tadpole.fasta   Per-assembler outputs
│   ├── <sample>.megahit.fasta
│   ├── <sample>.spades.fasta
│   ├── <sample>.metaspades.fasta
│   └── <sample>.assembly_stats.txt
├── mapping/<sample>/
│   ├── <sample>.sorted.bam      Alignments
│   ├── <sample>.sorted.bam.bai  BAM index
│   ├── <sample>.depths.txt      Depth table (MetaBAT2 format)
│   └── <sample>.covstats.txt    Coverage statistics
└── pipeline_info/
```

Feed these into mag_analysis:
```bash
../mag_analysis/nextflow/run-mag-analysis.sh \
    --assembly results/assembly/<sample>/<sample>.dedupe.fasta \
    --depths results/mapping/<sample>/<sample>.depths.txt \
    --bam_dir results/mapping/<sample>/ \
    --outdir /path/to/analysis --db_dir /path/to/databases
```

## Options

| Flag | Default | Description |
|------|---------|-------------|
| `--coassembly` | false | Co-assemble all samples (default: per-sample) |
| `--run_remove_human` | true | Human decontamination (requires `--human_ref`) |
| `--human_ref` | databases/human_ref | BBTools human reference index |
| `--run_fastqc` | true | FastQC on preprocessed reads |
| `--run_normalize` | true | bbnorm coverage normalization |
| `--run_tadpole` | true | Run Tadpole assembler |
| `--run_megahit` | true | Run Megahit assembler |
| `--run_spades` | true | Run SPAdes assembler |
| `--run_metaspades` | true | Run metaSPAdes assembler |
| `--dedupe_identity` | 98 | Cascade deduplication threshold |
| `--min_contig_len` | 500 | Minimum contig length after dedup |
| `--assembly_cpus` | 24 | CPUs for assembly |
| `--assembly_memory` | `250 GB` | Memory for assembly |

## Pipeline

```
Paired-end FASTQs → BBTools QC (clumpify, filterbytile, bbduk)
    → Human decontamination → FastQC
    → 3-phase error correction (ecco, ecc, tadpole)
    → bbnorm normalization → bbmerge read merging
         │
    ┌────┼────┬────┐
 Tadpole Megahit SPAdes metaSPAdes
    └────┼────┴────┘
    Cascade deduplication (100% → 99% → 98%)
    → BBMap read mapping → jgi_summarize_bam_contig_depths
```

19 processes, 4 conda environments.
