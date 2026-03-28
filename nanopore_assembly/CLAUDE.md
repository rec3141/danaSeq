# CLAUDE.md

## Project Overview

**dānaSeq Nanopore Assembly** preprocesses Oxford Nanopore reads and produces a co-assembly with depth tables and BAMs for downstream analysis by `mag_analysis`.

Implemented in **Nextflow DSL2** in `nextflow/`.

## Repository Structure

```
nanopore_assembly/
└── nextflow/
    ├── main.nf              Pipeline entry point (10 processes)
    ├── nextflow.config      Params, profiles, resources
    ├── run-nanopore-assembly.sh  Pipeline launcher (local/Docker/Apptainer)
    ├── install.sh           Conda env builder + C binary compiler
    ├── modules/
    │   ├── preprocess.nf    CONCAT_READS, PREPARE_READS, REMOVE_HUMAN
    │   ├── assembly.nf      FLYE_ASSEMBLE, ASSEMBLY_METAMDBG, ASSEMBLY_MYLOASM,
    │   │                    FLYE_POLISH, CALCULATE_TNF
    │   └── mapping.nf       MAP_READS, CALCULATE_DEPTHS
    ├── bin/
    │   ├── fastq_filter     Compiled C++ streaming QC filter
    │   ├── fastq_filter.cpp Source
    │   ├── tetramer_freqs   Compiled C tetranucleotide calculator
    │   └── tetramer_freqs.c Source
    ├── envs/                Conda YAML specs
    ├── Dockerfile           Thin layer on danaseq-assembly-base
    └── Dockerfile.base      Shared base image (nanopore + illumina envs)
```

## Running

```bash
cd nextflow
./install.sh && ./install.sh --check

./run-nanopore-assembly.sh --input /path/to/reads --outdir /path/to/output

# Then feed into mag_analysis:
../mag_analysis/nextflow/run-mag-analysis.sh \
    --assembly /path/to/output/assembly/assembly.fasta \
    --depths /path/to/output/mapping/depths.txt \
    --bam_dir /path/to/output/mapping/
```

## Output

```
results/
├── assembly/
│   ├── assembly.fasta       Co-assembly
│   ├── assembly_info.txt    Contig metadata (length, coverage, circularity)
│   ├── assembly_graph.gfa   Assembly graph
│   ├── tnf.tsv              Tetranucleotide frequencies (136 features)
│   └── gc.tsv               Per-contig GC%
├── mapping/
│   ├── *.sorted.bam         Per-sample alignments
│   ├── *.sorted.bam.bai     BAM indices
│   └── depths.txt           CoverM depth table (MetaBAT2 format)
└── pipeline_info/
```

## Key Design Decisions

- **CoverM** replaces jgi_summarize_bam_contig_depths (integer overflow bug fix)
- **Compiled C binaries** for fastq_filter and tetramer_freqs (31x faster than Python)
- **Three assembler options**: Flye (default), metaMDBG, myloasm via `--assembler`
- **Optional Flye polishing**: auto for Flye, disabled for others (built-in polishing)
