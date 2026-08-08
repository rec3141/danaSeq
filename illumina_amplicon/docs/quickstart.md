# Quickstart

## Prerequisites

- [Nextflow](https://www.nextflow.io/) >= 23.04
- [Conda](https://docs.conda.io/) or [Mamba](https://mamba.readthedocs.io/) (for `-profile conda`)
- Or: Docker / Apptainer (for container execution)

## Run from GitHub

No installation needed — Nextflow pulls the pipeline automatically:

```bash
nextflow run rec3141/microscape-nf \
    --input /path/to/reads \
    --ref_databases "silva:/path/to/silva_train_set.fasta:Domain,Phylum,Class,Order,Family,Genus" \
    -profile conda \
    -resume
```

## Run with Container

```bash
# Docker
nextflow run rec3141/microscape-nf \
    --input /path/to/reads \
    --ref_databases "silva:/path/to/silva.fasta:Domain,Phylum,Class,Order,Family,Genus" \
    -with-docker ghcr.io/rec3141/microscape-nf:latest \
    -resume

# Apptainer (HPC)
nextflow run rec3141/microscape-nf \
    --input /path/to/reads \
    --ref_databases "silva:/path/to/silva.fasta:Domain,Phylum,Class,Order,Family,Genus" \
    -with-singularity docker://ghcr.io/rec3141/microscape-nf:latest \
    -resume
```

## Run from Local Clone

```bash
git clone https://github.com/rec3141/microscape-nf.git
cd microscape-nf

nextflow run main.nf \
    --input /path/to/reads \
    --ref_databases "silva:/path/to/silva.fasta:Domain,Phylum,Class,Order,Family,Genus" \
    -resume
```

## Multi-Database Taxonomy

Classify against multiple reference databases in parallel:

```bash
nextflow run rec3141/microscape-nf \
    --input /path/to/reads \
    --ref_databases "silva:/db/silva.fasta:Domain,Phylum,Class,Order,Family,Genus;pr2:/db/pr2.fasta:Domain,Supergroup,Division,Class,Order,Family,Genus,Species" \
    --run_phylogeny \
    -profile conda -resume
```

## Persistent Cache

For large projects, use `--store_dir` to skip completed steps across runs:

```bash
nextflow run rec3141/microscape-nf \
    --input /path/to/reads \
    --ref_databases "silva:/db/silva.fasta:Domain,Phylum,Class,Order,Family,Genus" \
    --store_dir /scratch/microscape_cache \
    -profile conda -resume
```

## Outputs

Results are written to `--outdir` (default: `results/`):

```
results/
├── filtered/              Filtered FASTQ files
├── errors/                Error model plots
├── seqtab/                Per-plate sequence tables
├── merged/                Merged sequence table
├── chimeras/              Chimera-filtered table
├── filtered_seqtab/       QC-filtered final table
├── taxonomy/              Per-database classifications
├── phylogeny/             MSA + NJ tree (if --run_phylogeny)
├── renormalized/          Per-group proportional tables
├── ordination/            t-SNE/PCA coordinates
├── network/               Correlation edge lists
├── viz/                   JSON files for web viewer
└── pipeline_info/         Nextflow timeline, trace, DAG
```
