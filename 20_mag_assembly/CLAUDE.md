# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

**Important**: If user poses a question, agent responds with an answer, not a codebase change.

## Project Overview

**dānaSeq MAG Assembly** is a post-expedition metagenome-assembled genome (MAG) reconstruction pipeline. It co-assembles nanopore reads with Flye, maps reads back, runs three binning algorithms in parallel, and integrates results with DAS Tool consensus.

The pipeline is implemented in **Nextflow DSL2** in `nextflow/`. Legacy bash scripts are preserved in the parent directory for reference but are not actively maintained.

## Repository Structure

```
20_mag_assembly/
├── nextflow/                    Primary pipeline (Nextflow DSL2)
│   ├── main.nf                 Pipeline entry point
│   ├── nextflow.config         Params, profiles, resources
│   ├── modules/
│   │   ├── assembly.nf         ASSEMBLY_FLYE
│   │   ├── mapping.nf          MAP_READS, CALCULATE_DEPTHS
│   │   └── binning.nf          BIN_SEMIBIN2, BIN_METABAT2, BIN_MAXBIN2,
│   │                           DASTOOL_CONSENSUS
│   ├── envs/                   Conda YAML specs
│   │   ├── flye.yml            Flye, Filtlong, Nextflow, OpenJDK
│   │   ├── mapping.yml         minimap2, samtools, CoverM
│   │   ├── semibin.yml         SemiBin2, PyTorch GPU
│   │   ├── semibin-cpu.yml     SemiBin2, PyTorch CPU (for Docker)
│   │   ├── binning.yml         MetaBAT2, MaxBin2, DAS_Tool
│   │   └── bbmap.yml           BBMap (optional dedupe)
│   ├── conda-envs/             Pre-built envs (created by install.sh)
│   ├── install.sh              Conda environment builder
│   ├── Dockerfile              Docker image (CPU-only SemiBin2)
│   ├── entrypoint.sh           Docker entrypoint
│   ├── run-mag.sh              Pipeline launcher (local/Docker)
│   ├── .dockerignore           Excludes conda-envs/, work/ from image
│   └── .gitignore              Excludes runtime artifacts
├── CLAUDE.md                   This file
├── 61_map_and_bin_optimized.sh Legacy orchestrator (reference only)
└── ...                         Other legacy scripts
```

## Running the Pipeline

```bash
cd nextflow

# Install conda environments
./install.sh
./install.sh --check

# Basic run
./run-mag.sh --input /path/to/reads --outdir /path/to/output

# With Filtlong pre-filtering, no MaxBin2
./run-mag.sh --input /path/to/reads --outdir /path/to/output \
    --filtlong_size 40000000000 --run_maxbin false

# Kitchen sink — all options with defaults
./run-mag.sh --input /path/to/reads --outdir /path/to/output \
    --dedupe \
    --filtlong_size 40000000000 \
    --min_overlap 1000 \
    --run_maxbin true \
    --metabat_min_cls 50000 \
    --assembly_cpus 24 \
    --assembly_memory '64 GB'

# Docker mode
docker build -t danaseq-mag .
./run-mag.sh --docker --input /path/to/reads --outdir /path/to/output

# Quick test with bundled data (direct nextflow)
mamba run -p conda-envs/dana-mag-flye \
    nextflow run main.nf --input test-data -profile test -resume

# Show all options
./run-mag.sh --help
```

## Pipeline Architecture

### Processing DAG (fan-in → fan-out → fan-in)

```
Sample FASTQs (N files)
         │ collect()
   ASSEMBLY_FLYE          Fan-in: all reads → 1 co-assembly
         │
   MAP_READS (×N)         Fan-out: assembly × each sample → N BAMs
         │ collect()
   CALCULATE_DEPTHS       Fan-in: all BAMs → depth table (CoverM)
         │
    ┌────┼────┐
 SemiBin2 MetaBAT2 MaxBin2   Parallel binning
    └────┼────┘
   DASTOOL_CONSENSUS      Consensus integration
```

### Key Design Decisions

**CoverM for depth calculation:** Replaces `jgi_summarize_bam_contig_depths` which has an integer overflow bug in MetaBAT2 <=2.17 on long reads. CoverM handles supplementary alignments correctly and outputs MetaBAT2-compatible depth tables.

**Supplementary alignment filtering:** MAP_READS uses `samtools view -F 0x904` to drop unmapped, secondary, and supplementary alignments before sorting. Long reads produce chimeric alignments that cause massive depth overcounting.

**Dynamic binner architecture:** Binners emit `[label, file]` tuples that are mixed into a single channel and collected for DAS_Tool. New binners can be added by appending to `ch_binner_results` in `main.nf` -- no changes needed in the DAS_Tool process.

**Graceful failure handling:** SemiBin2 catches crashes on small datasets (0 bins → empty ORFs → hmmsearch fail) and produces an empty output file. DAS_Tool filters out empty binner inputs and handles the "no bins above score threshold" case. The pipeline completes successfully even when individual binners fail.

**GPU vs CPU PyTorch:** The local conda env (`semibin.yml`) includes `pytorch-gpu` for GPU-accelerated SemiBin2. The Docker image uses `semibin-cpu.yml` (~500MB vs ~5GB) since most deployments won't have `--gpus`.

**Resume:** Nextflow's built-in `-resume` uses task hashing. No manual checkpoint logic needed.

### Conda Environments

Five isolated environments avoid dependency conflicts:

| Environment | Tools | Rationale |
|-------------|-------|-----------|
| `dana-mag-flye` | Flye, Filtlong, Nextflow, OpenJDK | Python version conflicts; also hosts Nextflow runtime |
| `dana-mag-mapping` | minimap2, samtools, CoverM | Universal mapping tools |
| `dana-mag-semibin` | SemiBin2, PyTorch GPU | ML dependencies isolated |
| `dana-mag-binning` | MetaBAT2, MaxBin2, DAS_Tool | Binning suite |
| `dana-bbmap` | BBMap | Optional dedupe (only if `params.dedupe`) |

### Nextflow Config Profiles

| Profile | Use case |
|---------|----------|
| `standard` | Local execution |
| `test` | Small test data, reduced resources (4 CPUs, 8GB) |
| `shipboard` | Production: 32 CPUs, 256GB RAM |

## Expected Input/Output

### Input

```
input_dir/
├── sample1.fastq.gz
├── sample2.fastq.gz
└── ...
```

`--input` must point to a directory containing `*.fastq.gz` files. All files are co-assembled.

### Output

```
results/
├── assembly/
│   └── assembly.fasta         Co-assembly
├── mapping/
│   ├── *.sorted.bam           Per-sample alignments
│   ├── *.sorted.bam.bai       BAM indices
│   └── depths.txt             CoverM depth table
├── binning/
│   ├── semibin/contig_bins.tsv
│   ├── metabat/contig_bins.tsv
│   ├── maxbin/contig_bins.tsv
│   └── consensus/
│       ├── bins/*.fa           Final consensus MAG FASTAs
│       ├── contig2bin.tsv      Contig-to-bin assignments
│       ├── allbins.fa          All bins concatenated
│       └── scores.tsv          Per-bin quality scores
└── pipeline_info/
    ├── timeline.html
    ├── report.html
    └── trace.txt
```

## Development Notes

### Adding a New Binner

1. Add a process to `modules/binning.nf` that outputs `path("LABEL_bins.tsv")` in DAS_Tool format (contig\tbin, tab-separated)
2. In `main.nf`, invoke the process and mix its output into `ch_binner_results`:
   ```groovy
   BIN_NEWTOOL(ASSEMBLY_FLYE.out.assembly, CALCULATE_DEPTHS.out.jgi_depth)
   ch_binner_results = ch_binner_results.mix(
       BIN_NEWTOOL.out.bins.map { ['newtool', it] }
   )
   ```
3. Add a conda env YAML if needed, or add the tool to an existing env
4. Add a `--run_newtool` param to `nextflow.config` if it should be optional

### Modifying an Existing Module

Module files are in `nextflow/modules/*.nf`. Each process has:
- `conda` directive pointing to the pre-built env in `conda-envs/`
- `publishDir` for output routing (with `saveAs` to normalize filenames)
- Resource labels (`process_low`, `process_medium`, `process_high`)
- Graceful error handling for edge cases (empty output, tool crashes)

### Known Limitations

- **Small test data:** The bundled test data (~17MB) is too small for meaningful binning. SemiBin2 produces 0 bins, DAS_Tool finds no bins above the 0.5 score threshold. All handled gracefully.
- **DAS_Tool hanging:** On very small assemblies, DAS_Tool's ruby SCG annotation script can hang. The `set +e` error handling catches this on retry.
- **Docker GPU:** The Docker image uses CPU-only PyTorch. For GPU SemiBin2 in Docker, build with `semibin.yml` instead of `semibin-cpu.yml` and use `--gpus all`.

## Common Issues

### No FASTQ files found
`--input` must point to a directory containing `*.fastq.gz` files directly (not subdirectories).

### Conda environment build fails
Run `./install.sh --check` to diagnose. Ensure mamba is on PATH. Build requires internet access.

### Depth values are wrong (overflow/negative)
This pipeline uses CoverM instead of `jgi_summarize_bam_contig_depths` to avoid the MetaBAT2 integer overflow bug. If you see overflow values, ensure the CALCULATE_DEPTHS process is using CoverM (check `modules/mapping.nf`).

---

**Repository:** https://github.com/rec3141/danaSeq
**License:** MIT (see LICENSE file)
