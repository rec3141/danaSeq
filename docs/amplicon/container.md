# Container

Pre-built container images include Nextflow, papa2, cutadapt, MAFFT,
and FastTree — everything needed to run the full pipeline without
installing dependencies.

## Docker

```bash
# Pull the image
docker pull ghcr.io/rec3141/danaseq-illumina-amplicon:latest

# Run the pipeline
docker run --rm -v $(pwd):/data ghcr.io/rec3141/danaseq-illumina-amplicon \
    run /pipeline/main.nf \
    --input /data/reads \
    --ref_databases "silva:/data/silva.fasta:Domain,Phylum,Class,Order,Family,Genus" \
    --outdir /data/results \
    -resume
```

## Apptainer (HPC)

```bash
# Pull and convert to SIF
apptainer pull danaseq-illumina-amplicon.sif docker://ghcr.io/rec3141/danaseq-illumina-amplicon:latest

# Run the pipeline through the image entrypoint (`run`, not `exec`: the
# entrypoint loads /pipeline/docker.config, which turns conda off so the
# tools already on PATH are used)
apptainer run danaseq-illumina-amplicon.sif run /pipeline/main.nf \
    --input /path/to/reads \
    --ref_databases "silva:/path/to/silva.fasta:Domain,Phylum,Class,Order,Family,Genus" \
    --outdir results \
    -resume
```

## Launcher

`run-illumina-amplicon.sh` does the above for you: it picks Apptainer >
Singularity > Docker > local conda, pulls the image if needed (`--pull`, or
`--sif PATH` for a pre-built SIF), binds `--outdir` and `--input`, forwards
the host CA bundle, and appends `-resume`:

```bash
cd danaSeq/illumina_amplicon
./run-illumina-amplicon.sh --runtime apptainer \
    --input /path/to/reads \
    --ref_databases "silva:/path/to/silva.fasta:Domain,Phylum,Class,Order,Family,Genus"
```

There are no `docker`/`singularity` Nextflow profiles in this pipeline.

## What's Inside

The container includes:

| Component | Version | Source |
|-----------|---------|--------|
| Nextflow | latest | nextflow.io |
| papa2 | v1.0.1 (pinned tag in `envs/python.yml`, bumped deliberately) | [GitHub release](https://github.com/rec3141/papa2/releases) via pip |
| cutadapt | latest | bioconda |
| MAFFT + FastTree | latest | bioconda |

There is no R in the image: papa2 is byte-identical to R dada2 1.40
(see its [parity documentation](https://rec3141.github.io/papa2/parity/)),
so the pipeline runs on one engine.

## Building Locally

```bash
git clone https://github.com/rec3141/danaSeq.git
cd danaSeq
docker build -t danaseq-illumina-amplicon illumina_amplicon
```
