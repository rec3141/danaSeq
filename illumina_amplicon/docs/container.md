# Container

Pre-built container images include Nextflow, papa2, microscape, cutadapt,
MAFFT, and R/dada2 — everything needed to run the full pipeline without
installing dependencies.

## Docker

```bash
# Pull the image
docker pull ghcr.io/rec3141/microscape-nf:latest

# Run the pipeline
docker run --rm -v $(pwd):/data ghcr.io/rec3141/microscape-nf \
    run /pipeline/main.nf \
    --input /data/reads \
    --ref_databases "silva:/data/silva.fasta:Domain,Phylum,Class,Order,Family,Genus" \
    --outdir /data/results \
    -resume
```

## Apptainer (HPC)

```bash
# Pull and convert to SIF
apptainer pull microscape-nf.sif docker://ghcr.io/rec3141/microscape-nf:latest

# Run the pipeline
apptainer exec microscape-nf.sif nextflow run /pipeline/main.nf \
    --input /path/to/reads \
    --ref_databases "silva:/path/to/silva.fasta:Domain,Phylum,Class,Order,Family,Genus" \
    --outdir results \
    -resume
```

## Nextflow with Container Profile

Nextflow can pull and manage the container automatically:

```bash
# Docker
nextflow run rec3141/microscape-nf \
    --input /path/to/reads \
    --ref_databases "silva:/path/to/silva.fasta:Domain,Phylum,Class,Order,Family,Genus" \
    -with-docker ghcr.io/rec3141/microscape-nf:latest \
    -resume

# Singularity/Apptainer
nextflow run rec3141/microscape-nf \
    --input /path/to/reads \
    --ref_databases "silva:/path/to/silva.fasta:Domain,Phylum,Class,Order,Family,Genus" \
    -with-singularity docker://ghcr.io/rec3141/microscape-nf:latest \
    -resume
```

## What's Inside

The container includes both Python and R environments:

| Component | Version | Source |
|-----------|---------|--------|
| Nextflow | latest | nextflow.io |
| papa2 | 0.1.0 | bioconda |
| microscape | 0.1.0 | bioconda |
| cutadapt | latest | bioconda |
| MAFFT | latest | bioconda |
| R + dada2 | 4.x + 1.36 | bioconda |
| microscapeR | 0.99.0 | GitHub |

## Building Locally

```bash
git clone https://github.com/rec3141/microscape-nf.git
cd microscape-nf
docker build -t microscape-nf .
```
