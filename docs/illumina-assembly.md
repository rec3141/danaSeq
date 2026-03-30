# Illumina Assembly Pipeline

Processes Illumina paired-end metagenomic reads through multi-assembler consensus and produces assemblies with depth tables and BAMs for downstream analysis by [mag_analysis](mag-analysis.md).

## Quick Start

```bash
cd illumina_assembly/nextflow
./install.sh && ./install.sh --check

# Local (conda)
./run-illumina-assembly.sh --input /path/to/reads --outdir /path/to/output

# Co-assembly mode
./run-illumina-assembly.sh --input /path/to/reads --outdir /path/to/output --coassembly

# Skip human removal
./run-illumina-assembly.sh --input /path/to/reads --outdir /path/to/output \
    --run_remove_human false

# Apptainer (HPC) -- auto-pulls SIF on first run
./run-illumina-assembly.sh --apptainer --pull --input /path/to/reads --outdir /path/to/output

# SLURM profile (Compute Canada)
./run-illumina-assembly.sh --input /path/to/reads --outdir /path/to/output \
    -profile slurm --slurm_account def-myaccount \
    --conda_path ~/scratch/miniforge3/bin

# Show all options
./run-illumina-assembly.sh --help
```

## Pipeline Overview

```
Paired-end FASTQs (*_R1_*.fastq.gz, N samples)
         |
   CLUMPIFY (xN)              Optical deduplication
         |
   FILTER_BY_TILE (xN)        Remove low-quality tiles
         |
   BBDUK_TRIM (xN)            Adapter trimming (k=23, mink=11, minlen=70)
         |
   BBDUK_FILTER (xN)          Artifact + PhiX removal (k=31, entropy=0.95)
         |
   REMOVE_HUMAN (xN)          Human read removal (optional)
         |
   FASTQC (xN)                QC report on final preprocessed reads (optional)
         |
   ERROR_CORRECT_ECCO (xN)    Phase 1: overlap-based correction
         |
   ERROR_CORRECT_ECC (xN)     Phase 2: clump-based correction
         |
   ERROR_CORRECT_TADPOLE (xN) Phase 3: k-mer correction (k=62)
         |
   NORMALIZE_READS (xN)       Coverage normalization (target=100, mindepth=2)
         |
   MERGE_READS (xN)           Merge overlapping pairs (bbmerge-auto k=93)
   |              |
   merged     QUALITY_TRIM     Quality-trim unmerged reads
   |              |
   +---------+---+
   |         |         |              |
 Tadpole  Megahit   SPAdes      metaSPAdes      Four parallel assemblers
 (k=124)  (k=45-225) (k=25-125) (k=25-77,--meta)
   |         |         |              |
   +---------+---------+--------------+
         |
   DEDUPE_ASSEMBLIES           Cascade deduplication (100% -> 99% -> 98%)
         |
   MAP_READS_BBMAP (xN)       Map QC'd reads to assembly (minid=90)
         |
   CALCULATE_DEPTHS            jgi_summarize_bam_contig_depths

Output: assembly.fasta + depths.txt + BAMs
```

## Key Features

- **Four-assembler consensus**: Tadpole, Megahit, SPAdes, metaSPAdes with cascade deduplication
- **BBTools-based QC**: Optical deduplication, tile filtering, adapter/artifact removal, human decontamination
- **FastQC reports** on final preprocessed reads
- **Three-phase error correction**: Overlap, clump, and k-mer-based correction
- **Per-sample and co-assembly modes**: `--coassembly` pools all samples
- **SRA-safe**: Automatic fallbacks for SRA-stripped read headers
- **SLURM support**: `-profile slurm` for Compute Canada HPC clusters
- Assembly + depths output feeds directly into `mag_analysis`

## Input

`--input` must point to a directory containing paired-end Illumina reads:

```
input_dir/
├── sampleA_S1_L001_R1_001.fastq.gz
├── sampleA_S1_L001_R2_001.fastq.gz
├── sampleB_S2_L001_R1_001.fastq.gz
├── sampleB_S2_L001_R2_001.fastq.gz
└── ...
```

R2 files are auto-detected by replacing `_R1_` with `_R2_`. The sample ID is extracted as the first field before `_` in the filename.

## Output

```
results/
├── preprocess/<sample>/        QC'd reads + FastQC reports
├── error_correct/<sample>/     Three-phase error-corrected reads
├── normalize/<sample>/         Coverage-normalized reads + k-mer histograms
├── merge/<sample>/             Merged + quality-trimmed reads
├── assembly/<sample>/
│   ├── <sample>.dedupe.fasta   Final deduplicated assembly
│   ├── <sample>.tadpole.fasta  Per-assembler outputs
│   ├── <sample>.megahit.fasta
│   ├── <sample>.spades.fasta
│   ├── <sample>.metaspades.fasta
│   └── <sample>.assembly_stats.txt
├── mapping/<sample>/
│   ├── <sample>.sorted.bam     Alignments
│   ├── <sample>.sorted.bam.bai BAM index
│   ├── <sample>.depths.txt     Depth matrix (jgi format)
│   └── <sample>.covstats.txt   Per-contig coverage
└── pipeline_info/              Nextflow reports
```

In co-assembly mode, `<sample>` is replaced with `coassembly` for assembly and mapping directories.

## Parameters

### Required

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--input` | (required) | Directory containing paired-end `*_R1_*.fastq.gz` files |
| `--outdir` | `results` | Output directory |

### Mode

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--coassembly` | `false` | Co-assemble all samples together |

### Preprocessing

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--min_readlen` | `70` | Minimum read length after trimming |
| `--run_remove_human` | `true` | Remove human reads via removehuman.sh |
| `--human_ref` | `databases/human_ref` | Path to BBTools human reference index |
| `--run_fastqc` | `true` | Run FastQC on final preprocessed reads |
| `--run_normalize` | `true` | Enable bbnorm coverage normalization |

### Assembly

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--run_tadpole` | `true` | Run Tadpole assembler |
| `--run_megahit` | `true` | Run Megahit assembler |
| `--run_spades` | `true` | Run SPAdes assembler |
| `--run_metaspades` | `true` | Run metaSPAdes assembler |
| `--dedupe_identity` | `98` | Final deduplication identity threshold |

### Resources

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--assembly_cpus` | `24` | CPUs for assembly processes |
| `--assembly_memory` | `250 GB` | Memory for assembly processes |

### SLURM

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--slurm_account` | `def-rec3141` | SLURM `--account` for job submission |
| `--conda_path` | (none) | Path to conda/mamba `bin/` for SLURM jobs |

### Profiles

| Profile | Use case |
|---------|----------|
| `standard` | Local execution (default) |
| `test` | Small test data, reduced resources |
| `slurm` | SLURM cluster execution |

## Human Reference

Download the human reference index from the danaSeq root:

```bash
./download-databases.sh --human
```

## Design Notes

- **Multi-assembler approach.** Four assemblers capture different k-mer ranges. Cascade deduplication (100% -> 99% -> 98%) merges them into a single non-redundant contig set.
- **Three-phase error correction.** Overlap-based (bbmerge ecco), clump-based (clumpify ecc, 4 passes), and k-mer-based (tadpole ecc k=62) correction each catch different error types.
- **Human decontamination.** Enabled by default via BBTools' `removehuman.sh`. Disable with `--run_remove_human false`.
- **Interleaved read handling.** BBTools operates on interleaved FASTQ throughout, simplifying channel plumbing.
- **Graceful assembler failure.** All assemblers use `set +e` and produce an empty FASTA on error. The deduplication step skips empty files.
