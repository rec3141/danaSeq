# danaSeq · illumina_rna

Paired-end Illumina RNA-seq pipeline tailored for metatranscriptomic data mapped
against references from other danaSeq pipelines (assemblies from
`illumina_assembly`/`nanopore_assembly`, or annotated MAGs from `mag_analysis`).

## Stages

1. **QA/QC** — BBTools (clumpify → filterbytile → bbduk_trim → bbduk_filter →
   removehuman) + FastQC. Identical to `illumina_assembly`'s preprocess.
2. **rRNA depletion** — SortMeRNA against user-supplied rRNA FASTAs. Default-on
   for metatranscriptomic libraries; disable with `--run_remove_rrna false`.
3. **Mapping** — BBmap to each `<name>.fasta` in `--references`. Per-(sample, reference)
   sorted BAM, idxstats, flagstat, covstats.
4. **Quantification** — `featureCounts` from subread when a `<name>.gff` sits
   alongside the reference FASTA (e.g. from Bakta/Prokka annotation). Always emits
   `samtools idxstats` for contig-level abundance.
5. **Summarize** — Per-reference `gene_counts.tsv` (genes × samples) and JSONs
   for the viz.

## Quick start

```bash
# Install conda envs
./install.sh

# Run pipeline
./run-illumina-rna.sh \
    --input /path/to/reads \
    --references /path/to/refs_dir \
    --outdir /path/to/output \
    --human_ref /path/to/bbtools_human_index \
    --sortmerna_refs /path/to/sortmerna_fastas
```

`--references` should be a directory containing one or more `<name>.fasta`
files. If a matching `<name>.gff` is present (Bakta/Prokka output), gene-level
counts are emitted; otherwise quantification is contig-level only.

## Outputs

```
results/
├── preprocess/<sample>/                Cleaned reads + FastQC + SortMeRNA log
├── mapping/<sample>/<ref>/             sorted.bam(.bai), idxstats, flagstat, covstats
├── quantify/<sample>/<ref>/            featureCounts per-sample counts + summary
├── expression/<ref>/                   <ref>.gene_counts.tsv (merged, genes × samples)
├── viz/                                JSONs for the Svelte viz
└── pipeline_info/                      Nextflow timeline/report/trace/dag
```

## Container

```bash
# Build locally
docker build -t danaseq-illumina-rna .

# Run via the launcher
./run-illumina-rna.sh --container --input ... --references ... --outdir ...
```

## Viz

```bash
cd viz
npm install
VIZ_DATA_DIR=/path/to/output/viz npm run dev
# open http://localhost:5173
```

Pages: Overview (mapping-rate matrix, read-flow funnel), Samples (per-sample
detail), Expression (gene × sample heatmap, top-N by variance, log₂ transform),
References (per-contig read totals).
