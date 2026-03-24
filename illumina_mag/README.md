# Illumina MAG Pipeline

Illumina paired-end metagenomic assembly and binning pipeline. Processes reads through BBTools quality control (optical deduplication, tile filtering, adapter trimming, artifact/PhiX removal, human decontamination), FastQC reporting, three-phase error correction (overlap, clump, k-mer), coverage normalization, read merging, four parallel assemblers (Tadpole, Megahit, SPAdes, metaSPAdes), cascade deduplication, read mapping, depth calculation (jgi_summarize_bam_contig_depths), and MetaBAT2 binning. Supports per-sample (default) and co-assembly modes.

## Quick Start

```bash
cd nextflow

# Install conda environments (~15 min first time)
./install.sh
./install.sh --check

# Run (local conda, handles activation automatically)
./run-illumina-mag.sh --input /path/to/reads --outdir /path/to/output \
    --human_ref /path/to/hg19_index

# Co-assembly mode
./run-illumina-mag.sh --input /path/to/reads --outdir /path/to/output \
    --human_ref /path/to/hg19_index --coassembly

# Skip human removal (e.g. for non-clinical samples)
./run-illumina-mag.sh --input /path/to/reads --outdir /path/to/output \
    --run_remove_human false

# SLURM profile
./run-illumina-mag.sh --input /path/to/reads --outdir /path/to/output \
    --human_ref /path/to/hg19_index \
    -profile slurm --slurm_account def-rec3141 \
    --conda_path ~/scratch/miniforge3/bin

# Docker mode
./run-illumina-mag.sh --docker --input /path/to/reads --outdir /path/to/output

# Apptainer (HPC) — auto-pulls SIF on first run
./run-illumina-mag.sh --apptainer --pull --input /path/to/reads --outdir /path/to/output

# Show all options
./run-illumina-mag.sh --help
```

### Kitchen sink -- all options with defaults

```bash
cd nextflow
./run-illumina-mag.sh --input /data/reads --outdir /data/output \
    --coassembly \
    --min_readlen 70 \
    --run_remove_human true \
    --human_ref /data/hg19_index \
    --run_fastqc true \
    --run_normalize true \
    --run_tadpole true \
    --run_megahit true \
    --run_spades true \
    --run_metaspades true \
    --dedupe_identity 98 \
    --metabat_min_cls 2000 \
    --store_dir /scratch/illumina_mag_store \
    --assembly_cpus 24 \
    --assembly_memory '250 GB'
```

## Human Reference Index

Human decontamination is enabled by default. Download the reference once from the danaSeq root:

```bash
# From the danaSeq root directory
./download-databases.sh --human
```

This creates `databases/human_ref/` with the BBTools masked hg19 index. The pipeline default `--human_ref` already points there, so no extra flags needed.

To skip human removal entirely, use `--run_remove_human false`.

## Pipeline Overview

```
Paired-end FASTQs (*_R1_*.fastq.gz, N samples)
         |
   CLUMPIFY (xN)              Optical deduplication (clumpify dedupe optical)
         |
   FILTER_BY_TILE (xN)        Remove low-quality tiles (filterbytile)
         |
   BBDUK_TRIM (xN)            Adapter trimming (k=23, mink=11, minlen=70)
         |
   BBDUK_FILTER (xN)          Artifact + PhiX removal (k=31, entropy=0.95)
         |
   REMOVE_HUMAN (xN)          Human read removal (removehuman.sh, optional)
         |
   FASTQC (xN)                QC report on final preprocessed reads (optional)
         |
   ERROR_CORRECT_ECCO (xN)    Phase 1: overlap-based correction (bbmerge ecco)
         |
   ERROR_CORRECT_ECC (xN)     Phase 2: clump-based correction (clumpify ecc)
         |
   ERROR_CORRECT_TADPOLE (xN) Phase 3: k-mer correction (tadpole ecc k=62)
         |
   NORMALIZE_READS (xN)       Coverage normalization (bbnorm target=100, mindepth=2)
         |
   MERGE_READS (xN)           Merge overlapping pairs (bbmerge-auto k=93)
   |              |
   merged     QUALITY_TRIM     Quality-trim unmerged reads (qtrim=r trimq=10)
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
   MAP_READS_BBMAP (xN)       Map QC'd reads to assembly (bbmap, minid=90)
         |
   CALCULATE_DEPTHS            jgi_summarize_bam_contig_depths
         |
   BIN_METABAT2                MetaBAT2 binning (--saveCls)
```

### Per-sample vs co-assembly

**Per-sample mode** (default): each sample is independently preprocessed, assembled, deduplicated, mapped, and binned. Each sample produces its own assembly and bins.

**Co-assembly mode** (`--coassembly`): all samples are preprocessed independently, then their merged/quality-trimmed reads are pooled for a single combined assembly. All samples are mapped back to the co-assembly, and depths from all BAMs are used together for binning.

## Input

`--input` must point to a directory containing paired-end Illumina reads with the naming convention:

```
*_R1_*.fastq.gz
```

R2 files are auto-detected by replacing `_R1_` with `_R2_`. The sample ID is extracted as the first field before `_` in the filename.

Example: `4Qbjj-2017000184_S123_L001_R1_001.fastq.gz` produces sample ID `4Qbjj-2017000184`.

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
├── preprocess/<sample>/
│   ├── <sample>.clumped.fq.gz           Optically deduplicated reads
│   ├── <sample>.filtered_by_tile.fq.gz  Tile-filtered reads
│   ├── <sample>.trimmed.fq.gz           Adapter-trimmed reads
│   ├── <sample>.filtered.fq.gz          Artifact/PhiX-filtered reads
│   ├── <sample>.nohuman.fq.gz           Human-decontaminated reads
│   └── fastqc/
│       ├── <sample>.r1_fastqc.html      FastQC report (R1)
│       └── <sample>.r2_fastqc.html      FastQC report (R2)
├── error_correct/<sample>/
│   ├── <sample>.ecco.fq.gz             Overlap-corrected reads
│   ├── <sample>.eccc.fq.gz             Clump-corrected reads
│   └── <sample>.ecct.fq.gz             K-mer-corrected reads
├── normalize/<sample>/
│   ├── <sample>.normalized.fq.gz       Coverage-normalized reads
│   ├── <sample>.khist.txt              K-mer frequency histogram
│   └── <sample>.peaks.txt              Coverage peaks
├── merge/<sample>/
│   ├── <sample>.merged.fq.gz           Merged overlapping pairs
│   ├── <sample>.unmerged.fq.gz         Non-overlapping pairs
│   └── <sample>.qtrimmed.fq.gz         Quality-trimmed unmerged reads
├── assembly/<sample>/
│   ├── <sample>.tadpole.fasta          Tadpole contigs
│   ├── <sample>.megahit.fasta          Megahit contigs
│   ├── <sample>.spades.fasta           SPAdes contigs
│   ├── <sample>.metaspades.fasta       metaSPAdes contigs
│   ├── <sample>.dedupe.fasta           Deduplicated combined assembly
│   └── <sample>.assembly_stats.txt     Per-assembler + combined statistics
├── mapping/<sample>/
│   ├── <sample>.sorted.bam             Sorted alignment
│   ├── <sample>.sorted.bam.bai         BAM index
│   ├── <sample>.covhist.txt            Coverage histogram
│   ├── <sample>.covstats.txt           Per-contig coverage statistics
│   └── <sample>.depths.txt             Depth matrix (jgi format)
├── binning/<sample>/
│   ├── <sample>.metabat_bins.tsv       Contig-to-bin assignments (contig\tbin)
│   └── bins/
│       ├── metabat_001.fa              Bin FASTA files
│       ├── metabat_002.fa
│       └── ...
└── pipeline_info/
    ├── run_command.sh                  Exact re-runnable command (for -resume)
    ├── timeline.html                   Process execution timeline
    ├── report.html                     Nextflow execution report
    ├── trace.txt                       Per-task resource usage
    └── dag.html                        Pipeline DAG visualization
```

In co-assembly mode, `<sample>` is replaced with `coassembly` for the assembly, mapping, and binning directories.

## Parameters

### Required

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--input` | (required) | Directory containing paired-end `*_R1_*.fastq.gz` files |
| `--outdir` | `results` | Output directory |

### Mode

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--coassembly` | `false` | Co-assemble all samples together (default: per-sample assembly) |

### Caching

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--store_dir` | (off) | Persistent cache directory (storeDir); completed processes are skipped across runs even after `work/` cleanup |

### Preprocessing

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--min_readlen` | `70` | Minimum read length after adapter trimming and quality trimming |
| `--run_remove_human` | `true` | Remove human reads via removehuman.sh |
| `--human_ref` | `databases/human_ref` | Path to BBTools human reference index (run `./download-databases.sh --human`) |
| `--run_fastqc` | `true` | Run FastQC on final preprocessed reads |

### Normalization

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--run_normalize` | `true` | Enable bbnorm coverage normalization (target=100, mindepth=2) |

### Assembly

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--run_tadpole` | `true` | Run Tadpole assembler (k=124) |
| `--run_megahit` | `true` | Run Megahit assembler (k=45-225, step=26) |
| `--run_spades` | `true` | Run SPAdes assembler (k=25,55,95,125) |
| `--run_metaspades` | `true` | Run metaSPAdes assembler (k=25,55,77, `--meta`) |

### Deduplication

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--dedupe_identity` | `98` | Final deduplication identity threshold (cascade: 100 -> 99 -> this value) |

### Binning

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--metabat_min_cls` | `2000` | MetaBAT2 `--minClsSize` minimum cluster size (bp) |

### Resources

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--assembly_cpus` | `24` | CPUs for assembly and preprocessing (`process_high` label) |
| `--assembly_memory` | `250 GB` | Memory for assembly and preprocessing (`process_high` label) |

### SLURM

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--slurm_account` | `def-rec3141` | SLURM `--account` for job submission |
| `--conda_path` | (none) | Path to conda/mamba `bin/` directory for SLURM jobs (e.g. `~/scratch/miniforge3/bin`). Required when using `-profile slurm` because SLURM jobs do not source `.bashrc`. |

### Profiles

| Profile | Use case |
|---------|----------|
| `standard` | Local execution (default) |
| `test` | Small test data, reduced resources (4 CPUs, 8 GB) |
| `slurm` | SLURM cluster execution (requires `--slurm_account` and `--conda_path`) |

## SLURM Profile

On HPC clusters with SLURM, use `-profile slurm` to submit each pipeline process as a separate SLURM job. Two additional flags are needed:

- `--slurm_account` sets the SLURM `--account` for billing (default: `def-rec3141`).
- `--conda_path` must point to the `bin/` directory of your conda or mamba installation. SLURM jobs run in a minimal shell without `.bashrc`, so conda is not on PATH by default. The pipeline injects `export PATH=<conda_path>:$PATH` as a `beforeScript` for every process.

```bash
./run-illumina-mag.sh --input /project/reads --outdir /scratch/output \
    --human_ref /project/hg19_index \
    -profile slurm \
    --slurm_account def-rec3141 \
    --conda_path ~/scratch/miniforge3/bin
```

Resource allocation per process label:

| Label | CPUs | Memory | Time |
|-------|------|--------|------|
| `process_low` | 2 | 4 GB | 1 h |
| `process_medium` | 8 | 24 GB | 4 h |
| `process_high` | `--assembly_cpus` (24) | `--assembly_memory` (250 GB) | 24 h |

Processes automatically retry once on failure (`errorStrategy = 'retry'`, `maxRetries = 1`), then are ignored to allow the rest of the pipeline to complete.

## Conda Environments

Four isolated environments avoid dependency conflicts between the assemblers and the BBTools suite:

| Environment | Tools | Key binaries |
|-------------|-------|--------------|
| `dana-illumina-mag-bbmap` | BBTools suite + FastQC + samtools + Nextflow runtime | bbduk.sh, clumpify.sh, filterbytile.sh, bbmerge.sh, removehuman.sh, tadpole.sh, bbnorm.sh, dedupe.sh, bbmap.sh, fastqc, samtools, nextflow |
| `dana-illumina-mag-megahit` | Megahit assembler | megahit |
| `dana-illumina-mag-spades` | SPAdes (includes metaSPAdes) | spades.py |
| `dana-illumina-mag-binning` | MetaBAT2 + depth tools + samtools | metabat2, jgi_summarize_bam_contig_depths, samtools |

The `dana-illumina-mag-bbmap` environment also hosts the Nextflow runtime, so `run-illumina-mag.sh` uses `mamba run -p conda-envs/dana-illumina-mag-bbmap` to launch the pipeline.

### Installing environments

```bash
cd nextflow

# Install all 4 environments
./install.sh

# Install to a custom prefix
./install.sh --prefix /custom/path

# Check installation status
./install.sh --check

# Remove all environments and rebuild
./install.sh --clean
./install.sh
```

YAML specs are in `nextflow/envs/` (`bbmap.yml`, `megahit.yml`, `spades.yml`, `binning.yml`). Environments are prefix-installed under `nextflow/conda-envs/`. Requires conda or mamba on `$PATH`.

## Resuming Runs

`run-illumina-mag.sh` automatically records the Nextflow session ID in:

```
<outdir>/pipeline_info/run_command.sh
```

To re-run (e.g. after adding samples or fixing a failure):

```bash
cd nextflow
tail -1 <outdir>/pipeline_info/run_command.sh | bash
```

Session handling works three ways:

1. **Auto-detect (default):** reads the last session ID from `run_command.sh` if it exists.
2. **Explicit:** `--session <uuid>` to resume from a specific session.
3. **Post-run capture:** after each run, extracts the actual session UUID from `.nextflow.log` and writes it into `run_command.sh` for future runs.

For persistent caching that survives `work/` cleanup, use `--store_dir`:

```bash
./run-illumina-mag.sh --input /data/reads --outdir /data/output \
    --human_ref /data/hg19_index \
    --store_dir /scratch/illumina_mag_store
```

## Design Notes

**Multi-assembler approach.** Running four assemblers (Tadpole, Megahit, SPAdes, metaSPAdes) on the same reads captures different k-mer ranges and assembly algorithms. The cascade deduplication (100% -> 99% -> 98%) merges them into a single non-redundant contig set, keeping the longest representative of each near-duplicate cluster.

**Three-phase error correction.** Reads pass through overlap-based correction (bbmerge ecco), clump-based correction (clumpify ecc, 4 passes), and k-mer-based correction (tadpole ecc k=62). Each phase catches different error types with complementary strategies.

**Human decontamination.** Enabled by default via BBTools' `removehuman.sh`. Uses a masked human genome reference to identify and remove human-derived reads before error correction and assembly. Can be disabled with `--run_remove_human false` for datasets where human contamination is not a concern.

**Interleaved read handling.** BBTools operates on interleaved FASTQ throughout (R1 and R2 in the same file). This simplifies channel plumbing and avoids paired-read synchronization issues across processes.

**metaSPAdes uses normalized reads.** Unlike the other three assemblers which use merged + quality-trimmed reads, metaSPAdes receives the normalized interleaved reads directly. This matches the pipeline design where metaSPAdes benefits from normalized coverage.

**Graceful assembler failure.** All assembler processes use `set +e` to catch failures and produce an empty FASTA on error. The deduplication step skips empty files. This means a single assembler crash does not halt the pipeline.

**jgi_summarize_bam_contig_depths.** Unlike the long-read MAG pipeline (which uses CoverM to avoid integer overflow with supplementary alignments), this short-read pipeline uses `jgi_summarize_bam_contig_depths` directly. Short reads produce standard primary alignments without the supplementary alignment inflation seen in long-read data.

**Resume.** Nextflow's built-in `-resume` uses task hashing. No manual checkpoint logic needed. Processes whose script block has not changed are skipped automatically.

## Pipeline Structure

```
illumina_mag/
├── README.md                    This file
└── nextflow/                    Pipeline implementation (Nextflow DSL2)
    ├── main.nf                  Pipeline entry point + workflow logic
    ├── nextflow.config          Parameters, profiles, resource labels
    ├── run-illumina-mag.sh      Pipeline launcher (conda, Docker, Apptainer)
    ├── install.sh               Conda environment builder (install/check/clean)
    ├── Dockerfile               Container image build (all 4 envs + wrappers)
    ├── entrypoint.sh            Container entrypoint (outdir validation + Nextflow)
    ├── .dockerignore            Files excluded from container build
    ├── modules/
    │   ├── preprocess.nf        CLUMPIFY, FILTER_BY_TILE, BBDUK_TRIM, BBDUK_FILTER, REMOVE_HUMAN, FASTQC
    │   ├── error_correct.nf     ERROR_CORRECT_ECCO, ERROR_CORRECT_ECC, ERROR_CORRECT_TADPOLE
    │   ├── normalize.nf         NORMALIZE_READS
    │   ├── merge_reads.nf       MERGE_READS, QUALITY_TRIM
    │   ├── assembly.nf          ASSEMBLE_TADPOLE, ASSEMBLE_MEGAHIT, ASSEMBLE_SPADES, ASSEMBLE_METASPADES
    │   ├── dedupe.nf            DEDUPE_ASSEMBLIES
    │   ├── mapping.nf           MAP_READS_BBMAP, CALCULATE_DEPTHS
    │   └── binning.nf           BIN_METABAT2
    ├── envs/                    Conda YAML specs
    │   ├── bbmap.yml            BBTools, FastQC, samtools, Nextflow, OpenJDK
    │   ├── megahit.yml          Megahit
    │   ├── spades.yml           SPAdes
    │   └── binning.yml          MetaBAT2, jgi_summarize_bam_contig_depths, samtools
    └── conda-envs/              Pre-built conda environments (created by install.sh)
        ├── dana-illumina-mag-bbmap/
        ├── dana-illumina-mag-megahit/
        ├── dana-illumina-mag-spades/
        └── dana-illumina-mag-binning/
```

## References

**Quality Control & Error Correction:**
- BBTools: Bushnell, [sourceforge.net/projects/bbmap](https://sourceforge.net/projects/bbmap/)
- FastQC: Andrews, [bioinformatics.babraham.ac.uk/projects/fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

**Assemblers:**
- SPAdes/metaSPAdes: Nurk et al., *Genome Research* 2017
- Megahit: Li et al., *Bioinformatics* 2015
- Tadpole: Bushnell, [sourceforge.net/projects/bbmap](https://sourceforge.net/projects/bbmap/)

**Binning:**
- MetaBAT2: Kang et al., *PeerJ* 2019
