# CLAUDE.md

## Project Overview

**dānaSeq Illumina Assembly** processes Illumina paired-end reads through multi-assembler consensus and produces assemblies with depth tables and BAMs for downstream analysis by `mag_analysis`.

Implemented in **Nextflow DSL2** in `nextflow/`.

## Repository Structure

```
illumina_assembly/
└── nextflow/
    ├── main.nf              Pipeline entry point (19 processes)
    ├── nextflow.config      Params, profiles, resources
    ├── run-illumina-assembly.sh  Pipeline launcher (local/Docker/Apptainer)
    ├── install.sh           Conda env builder
    ├── modules/
    │   ├── preprocess.nf    CLUMPIFY, FILTER_BY_TILE, BBDUK_TRIM,
    │   │                    BBDUK_FILTER, REMOVE_HUMAN, FASTQC
    │   ├── error_correct.nf ERROR_CORRECT_ECCO, ECC, TADPOLE
    │   ├── normalize.nf     NORMALIZE_READS
    │   ├── merge_reads.nf   MERGE_READS, QUALITY_TRIM
    │   ├── assembly.nf      ASSEMBLE_TADPOLE, MEGAHIT, SPADES, METASPADES
    │   ├── dedupe.nf        DEDUPE_ASSEMBLIES (cascade 100%->99%->98%)
    │   └── mapping.nf       MAP_READS_BBMAP, CALCULATE_DEPTHS
    ├── envs/                Conda YAML specs (4 environments)
    └── Dockerfile           Thin layer on danaseq-assembly-base
```

## Running

```bash
cd nextflow
./install.sh && ./install.sh --check

./run-illumina-assembly.sh --input /path/to/reads --outdir /path/to/output

# Then feed into mag_analysis:
../mag_analysis/nextflow/run-mag-analysis.sh \
    --assembly /path/to/output/assembly/<sample>/<sample>.dedupe.fasta \
    --depths /path/to/output/mapping/<sample>/<sample>.depths.txt \
    --bam_dir /path/to/output/mapping/<sample>/
```

## Modes

- **Per-sample** (default): each sample assembled independently
- **Co-assembly** (`--coassembly`): all samples pooled into one assembly

## Output

```
results/
├── preprocess/<sample>/     QC'd reads + FastQC reports
├── error_correct/<sample>/  3-phase error-corrected reads
├── normalize/<sample>/      Coverage-normalized reads
├── merge/<sample>/          Merged + quality-trimmed reads
├── assembly/<sample>/       Per-assembler + deduplicated contigs
├── mapping/<sample>/        BAM + coverage stats + depth table
└── pipeline_info/
```

## Key Design Decisions

- **Four-assembler consensus**: Tadpole + Megahit + SPAdes + metaSPAdes, cascade deduplicated
- **Three-phase error correction**: overlap (ecco), clump (ecc), k-mer (tadpole)
- **SRA-safe**: automatic fallbacks for SRA-stripped read headers
- **SLURM support**: `-profile slurm` for HPC clusters
