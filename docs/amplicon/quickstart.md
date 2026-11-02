# Quickstart

## Prerequisites

- [Nextflow](https://www.nextflow.io/) >= 23.04
- [Conda](https://docs.conda.io/) or [Mamba](https://mamba.readthedocs.io/) (environments are created automatically when running outside a container)
- Or: Docker / Apptainer (for container execution)

The pipeline lives in `illumina_amplicon/` of the repository; there is no
root `main.nf`, so it is run from a local clone.

## Run from Local Clone

```bash
git clone https://github.com/rec3141/danaSeq.git
cd danaSeq/illumina_amplicon

nextflow run main.nf \
    --input /path/to/reads \
    --ref_databases "silva:/path/to/silva.fasta:Domain,Phylum,Class,Order,Family,Genus" \
    -resume
```

## Run with Container

The launcher picks Apptainer > Singularity > Docker > local conda, binds
`--outdir` and `--input`, and appends `-resume`:

```bash
cd danaSeq/illumina_amplicon

./run-illumina-amplicon.sh \
    --input /path/to/reads \
    --ref_databases "silva:/path/to/silva.fasta:Domain,Phylum,Class,Order,Family,Genus"

# Force a runtime
./run-illumina-amplicon.sh --runtime apptainer --input /path/to/reads ...
```

See [Container](container.md) for the image and for running it without the
launcher.

## Multi-Database Taxonomy

Classify against multiple reference databases in parallel:

```bash
nextflow run main.nf \
    --input /path/to/reads \
    --ref_databases "silva:/db/silva.fasta:Domain,Phylum,Class,Order,Family,Genus;pr2:/db/pr2.fasta:Domain,Supergroup,Division,Subdivision,Class,Order,Family,Genus,Species" \
    --run_phylogeny \
    -resume
```

The first database named is the primary one: renormalization groups, the
per-ASV lineage strings and the viewer's default labels come from it.

## Persistent Cache

For large projects, use `--store_dir` to skip completed steps across runs:

```bash
nextflow run main.nf \
    --input /path/to/reads \
    --ref_databases "silva:/db/silva.fasta:Domain,Phylum,Class,Order,Family,Genus" \
    --store_dir /scratch/danaseq_amplicon_cache \
    -resume
```

## Outputs

Results are written to `--outdir` (default: `results/`):

```
results/
├── demultiplexed/         Demultiplexed reads (--run_demultiplex only)
├── trimmed/               cutadapt logs, detected primers (detected/), primer assignment
├── quality_check/         Per-sample auto-trim and per-group truncation policy tables
├── filtered/              Filter/trim statistics
├── error_models/          Per-group error models and plots
├── seqtabs/               Per-group sequence tables
├── seqtab_final/          Merged, chimera-filtered and QC-filtered tables and stats
├── taxonomy/              Per-database classifications
├── phylogeny/             Alignment + FastTree tree (--run_phylogeny)
├── renormalized/          Per-group proportional tables
├── clustering/            t-SNE coordinates (samples and ASVs)
├── network/               Correlation edge list and stats
├── viz/                   JSON data for the web viewer (+ tree.nwk, run_manifest.json)
├── site/                  Built static viewer (--build_viz_site)
└── pipeline_info/         Nextflow timeline, report, trace, DAG, pipeline_timing.tsv
```

With `--store_dir`, the per-sample intermediates (`demultiplexed/`,
`filtered/`, `error_models/`, `seqtabs/`, cutadapt logs) go to the cache
directory instead of `--outdir`.
