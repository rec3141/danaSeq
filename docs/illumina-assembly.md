# Illumina Assembly Pipeline

Multi-assembler consensus pipeline for Illumina paired-end metagenomic reads. Produces deduplicated assemblies with depth tables and BAMs for downstream MAG analysis.

## Quick Start

```bash
cd illumina_assembly
./install.sh && ./install.sh --check

# Local (conda)
./run-illumina-assembly.sh --input /path/to/reads --outdir /path/to/output

# Co-assembly mode
./run-illumina-assembly.sh --input /path/to/reads --outdir /path/to/output --coassembly

# Apptainer (HPC) -- auto-pulls SIF on first run
./run-illumina-assembly.sh --apptainer --pull --input /path/to/reads --outdir /path/to/output

# SLURM profile
./run-illumina-assembly.sh --input /path/to/reads --outdir /path/to/output \
    -profile slurm --slurm_account def-myaccount \
    --conda_path ~/scratch/miniforge3/bin
```

## Pipeline Stages

| Step | Process | Tool | Description |
|------|---------|------|-------------|
| 1 | `CLUMPIFY` | BBTools clumpify | Optical deduplication |
| 2 | `FILTER_BY_TILE` | BBTools filterbytile | Remove low-quality tiles |
| 3 | `BBDUK_TRIM` | BBDuk | Adapter trimming (k=23, mink=11, minlen=70) |
| 4 | `BBDUK_FILTER` | BBDuk | Artifact + PhiX removal (k=31, entropy=0.95) |
| 5 | `REMOVE_HUMAN` | BBTools removehuman | Human read removal (optional) |
| 6 | `FASTQC` | FastQC | QC report on preprocessed reads (optional) |
| 7 | `ERROR_CORRECT_ECCO` | BBMerge ecco | Phase 1: overlap-based error correction |
| 8 | `ERROR_CORRECT_ECC` | BBTools clumpify ecc | Phase 2: clump-based error correction (4 passes) |
| 9 | `ERROR_CORRECT_TADPOLE` | Tadpole | Phase 3: k-mer-based error correction (k=62) |
| 10 | `NORMALIZE_READS` | bbnorm | Coverage normalization (target=100, mindepth=2) |
| 11 | `MERGE_READS` | BBMerge | Merge overlapping pairs (k=93) |
| 12 | `QUALITY_TRIM` | BBDuk | Quality-trim unmerged reads |
| 13a | `ASSEMBLE_TADPOLE` | Tadpole | Assembly (k=124) |
| 13b | `ASSEMBLE_MEGAHIT` | Megahit | Assembly (k=45-225) |
| 13c | `ASSEMBLE_SPADES` | SPAdes | Assembly (k=25-125) |
| 13d | `ASSEMBLE_METASPADES` | metaSPAdes | Assembly (k=25-77, --meta) |
| 14 | `DEDUPE_ASSEMBLIES` | BBTools dedupe | Cascade deduplication (100% -> 99% -> 98%) |
| 15 | `MAP_READS_BBMAP` | BBMap | Map QC'd reads to assembly (minid=90) |
| 16 | `CALCULATE_DEPTHS` | jgi_summarize_bam_contig_depths | Depth matrix in MetaBAT2 format |

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
| `--min_contig_len` | `500` | Minimum contig length after deduplication |

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

## Outputs

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
│   ├── <sample>.depths.txt     Depth matrix (MetaBAT2 format)
│   └── <sample>.covstats.txt   Per-contig coverage statistics
└── pipeline_info/              Nextflow reports (timeline, trace, DAG)
```

In co-assembly mode, `<sample>` is replaced with `coassembly` for assembly and mapping directories.

Feed these outputs into `mag_analysis`:

```bash
cd ../mag_analysis
./run-mag-analysis.sh \
    --assembly results/assembly/<sample>/<sample>.dedupe.fasta \
    --depths results/mapping/<sample>/<sample>.depths.txt \
    --bam_dir results/mapping/<sample>/ \
    --outdir /path/to/analysis --db_dir /path/to/databases
```

## Profiles

| Profile | Use case |
|---------|----------|
| `standard` | Local execution (default) |
| `test` | Small test data, reduced resources (4 CPUs, 8 GB) |
| `slurm` | SLURM cluster execution |

## Resource Requirements

| Component | CPUs | RAM | Notes |
|-----------|------|-----|-------|
| Preprocessing (BBTools) | 8 | 16 GB | Per-sample; parallelized |
| Error correction | 8 | 16 GB | Three sequential phases per sample |
| SPAdes / metaSPAdes | 24 | 250 GB | Most memory-intensive step |
| Megahit | 24 | 250 GB | Lower peak memory than SPAdes |
| Tadpole | 24 | 250 GB | Fast, single k-mer assembly |
| BBMap read mapping | 8 | 16 GB | Per-sample mapping |
| Human decontamination | 8 | 16 GB | BBTools removehuman.sh |

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

R2 files are matched by replacing `_R1_` with `_R2_`. The sample ID is extracted as the first field before `_` in the filename.

## Design Notes

- **Multi-assembler approach.** Four assemblers capture different k-mer ranges. Cascade deduplication (100% -> 99% -> 98%) merges them into a single non-redundant contig set.
- **Three-phase error correction.** Overlap-based (bbmerge ecco), clump-based (clumpify ecc, 4 passes), and k-mer-based (tadpole ecc k=62) correction each catch different error types.
- **Human decontamination.** Enabled by default via BBTools' `removehuman.sh`. Disable with `--run_remove_human false`.
- **Interleaved read handling.** BBTools operates on interleaved FASTQ throughout, simplifying channel plumbing.
- **Graceful assembler failure.** All assemblers use `set +e` and produce an empty FASTA on error. The deduplication step skips empty files.
- **SRA-safe.** Automatic fallbacks for SRA-stripped read headers that lack optical duplicate coordinates.
