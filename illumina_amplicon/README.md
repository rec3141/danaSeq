# microscape-nf

**Amplicon sequencing analysis pipeline — from raw reads to interactive visualization**

[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A523.04-brightgreen)](https://www.nextflow.io/)
[![License: BSD-3-Clause](https://img.shields.io/badge/License-BSD_3--Clause-blue.svg)](LICENSE)

A Nextflow DSL2 pipeline for amplicon sequencing analysis. Takes demultiplexed
paired-end FASTQ files and produces ASV tables, taxonomy, phylogeny, ordinations,
and correlation networks.

## Quick Start

```bash
nextflow run rec3141/microscape-nf \
    --input /path/to/reads \
    --ref_databases "silva:/path/to/silva_train_set.fasta:Domain,Phylum,Class,Order,Family,Genus" \
    -profile conda \
    -resume
```

## Dependencies

The pipeline pulls its tools from [bioconda](https://bioconda.github.io/):

- **[papa2](https://github.com/rec3141/papa2)** — DADA2 denoising (Python, `conda install -c bioconda papa2`)
- **[microscape](https://github.com/rec3141/microscape)** — Downstream analysis (Python, `conda install -c bioconda microscape`)
- **[cutadapt](https://cutadapt.readthedocs.io/)** — Primer removal
- **[MAFFT](https://mafft.cbrc.jp/alignment/software/)** — Multiple sequence alignment

For R users, the pipeline also supports `--lang R` which uses:
- **[dada2](https://bioconductor.org/packages/dada2/)** — R/Bioconductor
- **[microscapeR](https://github.com/rec3141/microscapeR)** — R companion package

## Pipeline Stages

| Stage | Process | Description |
|-------|---------|-------------|
| 1 | REMOVE_PRIMERS | Primer trimming with cutadapt |
| 2 | DADA2_FILTER_TRIM | Quality filtering (maxEE, truncQ, PhiX removal) |
| 3 | DADA2_LEARN_ERRORS | Per-plate error model learning |
| 4 | DADA2_DENOISE | Denoising, pair merging, per-plate chimera removal |
| 5 | MERGE_SEQTABS | Merge per-plate sequence tables |
| 6 | REMOVE_CHIMERAS | Consensus chimera removal on merged data |
| 7 | FILTER_SEQTAB | Length, prevalence, abundance, and depth filtering |
| 8 | ASSIGN_TAXONOMY | Naive Bayesian classification (parallel per ref DB) |
| 9 | BUILD_PHYLOGENY | MSA + neighbor-joining tree (optional) |
| 10 | RENORMALIZE | Taxonomic grouping and within-group proportions |
| 11 | ORDINATE | t-SNE ordination of samples and ASVs |
| 12 | NETWORK | SparCC correlation networks |

## Parameters

```bash
nextflow run rec3141/microscape-nf --help
```

### Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--input` | required | Directory of paired-end `*.fastq.gz` files |
| `--ref_databases` | required | Reference DBs (`"name:path:Levels;..."`) |
| `--outdir` | `results` | Output directory |
| `--lang` | `python` | Language: `python` (papa2) or `R` (dada2) |
| `--maxEE` | `2` | Max expected errors per read |
| `--truncQ` | `11` | Truncate at first base with quality <= Q |
| `--min_overlap` | `10` | Min overlap for pair merging |
| `--run_phylogeny` | `false` | Build phylogenetic tree |
| `--threads` | `8` | CPU threads |

## Profiles

```bash
-profile conda      # Auto-create conda environments
-profile docker     # Use Docker container
-profile singularity # Use Singularity/Apptainer
-profile slurm      # Submit to SLURM cluster
-profile test       # Reduced resources for testing
```

## Related Packages

- **[papa2](https://github.com/rec3141/papa2)** — Python DADA2 port (bioconda)
- **[microscape](https://github.com/rec3141/microscape)** — Python downstream analysis (bioconda)
- **[microscapeR](https://github.com/rec3141/microscapeR)** — R downstream analysis (Bioconductor)

## Citation

> Callahan BJ, McMurdie PJ, Rosen MJ, Han AW, Johnson AJA, Holmes SP (2016).
> DADA2: High-resolution sample inference from Illumina amplicon data.
> *Nature Methods*, 13, 581-583.

## License

BSD-3-Clause — see [LICENSE](LICENSE).
