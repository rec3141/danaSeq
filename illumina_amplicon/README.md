# danaSeq illumina_amplicon

**Amplicon sequencing analysis pipeline — from raw reads to interactive visualization**

[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A523.04-brightgreen)](https://www.nextflow.io/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

A Nextflow DSL2 pipeline for amplicon sequencing analysis. Takes demultiplexed
paired-end FASTQ files and produces ASV tables, taxonomy, phylogeny, ordinations,
and correlation networks.

## Quick Start

```bash
git clone https://github.com/rec3141/danaSeq.git
cd danaSeq/illumina_amplicon

# Launcher: picks Apptainer > Singularity > Docker > local conda, appends -resume
./run-illumina-amplicon.sh \
    --input /path/to/reads \
    --ref_databases "silva:/path/to/silva_train_set.fasta:Domain,Phylum,Class,Order,Family,Genus"

# Or run Nextflow directly (conda environments are created automatically)
nextflow run main.nf \
    --input /path/to/reads \
    --ref_databases "silva:/path/to/silva_train_set.fasta:Domain,Phylum,Class,Order,Family,Genus" \
    -resume
```

The pipeline lives in this subdirectory; there is no root `main.nf`, so it is
run from a local clone (or from the container image, see
[the docs](https://rec3141.github.io/danaSeq/amplicon/container/)).

## Dependencies

The pipeline pulls its tools from [conda-forge](https://conda-forge.org/) and
[bioconda](https://bioconda.github.io/) (`envs/python.yml`), except papa2,
which is pip-installed from GitHub at a pinned tag:

- **[papa2](https://github.com/rec3141/papa2)** — DADA2 denoising and taxonomy (Python, `pip install git+https://github.com/rec3141/papa2.git@v1.0.1`)
- **[cutadapt](https://cutadapt.readthedocs.io/)** — Primer removal
- **[MAFFT](https://mafft.cbrc.jp/alignment/software/)** — Multiple sequence alignment
- **[FastTree](http://www.microbesonline.org/fasttree/)** — Phylogenetic tree inference

## Pipeline Stages

| Stage | Process | Description |
|-------|---------|-------------|
| 1 | REMOVE_PRIMERS | Primer trimming with cutadapt |
| 2 | FILTER_TRIM | Quality filtering (maxEE, truncQ, PhiX removal) |
| 3 | LEARN_ERRORS | Per-plate error model learning |
| 4 | DENOISE | Denoising, pair merging, per-plate chimera removal |
| 5 | MERGE_SEQTABS | Merge per-plate sequence tables |
| 6 | REMOVE_CHIMERAS | Chimera removal on merged data (pooled by default, `--chimera_method`) |
| 7 | FILTER_SEQTAB | Length, prevalence, abundance, and depth filtering |
| 8 | ASSIGN_TAXONOMY | Naive Bayesian classification (parallel per ref DB) |
| 9 | BUILD_PHYLOGENY | MAFFT alignment + FastTree tree (optional) |
| 10 | RENORMALIZE | Taxonomic grouping and within-group proportions |
| 11 | CLUSTER_TSNE | t-SNE ordination of samples and ASVs |
| 12 | NETWORK_SPARCC | CLR + Pearson co-occurrence network (not the SparCC algorithm) |

## Parameters

```bash
nextflow run main.nf --help
```

### Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--input` | required | Directory of paired-end `*.fastq.gz` files |
| `--ref_databases` | required | Reference DBs (`"name:path:Levels;..."`) |
| `--outdir` | `results` | Output directory |
| `--maxEE` | `2` | Max expected errors per read |
| `--truncQ` | `2` | Truncate at first base with quality <= Q |
| `--min_overlap` | `12` | Min overlap for pair merging |
| `--run_phylogeny` | `false` | Build phylogenetic tree |
| `--denoise_cpus` | `8` | CPUs for the heavy (`process_high`) steps |
| `--denoise_memory` | `16 GB` | Memory for the heavy (`process_high`) steps |

## Profiles

```bash
-profile standard   # Local executor (default)
-profile slurm      # Submit each task to Slurm
-profile test       # Reduced resources for testing
```

Conda environments are created automatically (`conda.enabled = true` in
`nextflow.config`); there is no `conda` profile. There are no `docker` or
`singularity` profiles either: container execution goes through
`./run-illumina-amplicon.sh` (`--runtime docker|apptainer|singularity|conda`),
which runs the image's entrypoint with conda disabled.

## Related

- **[papa2](https://github.com/rec3141/papa2)** — Python DADA2 port
- **[danaSeq](https://github.com/rec3141/danaSeq)** — this pipeline's home, alongside the metagenomics stages
- **[microscape.app](https://microscape.app)** — hosting and sharing the interactive viz

## Citation

> Callahan BJ, McMurdie PJ, Rosen MJ, Han AW, Johnson AJA, Holmes SP (2016).
> DADA2: High-resolution sample inference from Illumina amplicon data.
> *Nature Methods*, 13, 581-583.

## License

MIT, the same licence as the rest of danaSeq; see [LICENSE](LICENSE).
