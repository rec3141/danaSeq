# Illumina RNA Pipeline

Paired-end Illumina RNA-seq pipeline tailored for metatranscriptomic data mapped against references from other danaSeq pipelines (assemblies from `illumina_assembly`/`nanopore_assembly`, or annotated MAGs from `mag_analysis`).

## Quick Start

```bash
cd illumina_rna
./install.sh && ./install.sh --check

# Local (conda)
./run-illumina-rna.sh \
    --input /path/to/reads \
    --references /path/to/refs_dir \
    --outdir /path/to/output \
    --human_ref /path/to/bbtools_human_index \
    --sortmerna_refs /path/to/sortmerna_fastas

# Apptainer (HPC) -- auto-pulls SIF on first run
./run-illumina-rna.sh --apptainer --pull \
    --input /path/to/reads --references /path/to/refs --outdir /path/to/output

# SLURM profile
./run-illumina-rna.sh --input /path/to/reads --references /path/to/refs --outdir /path/to/output \
    -profile slurm --slurm_account def-myaccount \
    --conda_path ~/scratch/miniforge3/bin
```

## Inputs

- **`--input`** — directory of paired-end reads. Supports two naming conventions:
    - Illumina: `*_R1_*.fastq.gz` / `*_R2_*.fastq.gz`
    - BGI/MGI:  `*_1.fq.gz` / `*_2.fq.gz` (and `*.fastq.gz` variants)
- **`--references`** — directory of `<name>.fasta` files. If a sibling `<name>.gff` is present (typically the output of `bakta`/`prokka` annotation), gene-level counts are emitted via `featureCounts`; otherwise quantification is contig-level only via `samtools idxstats`.

## Pipeline Stages

| Step | Process | Tool | Description |
|------|---------|------|-------------|
| 1 | `CLUMPIFY` | BBTools clumpify | Optical deduplication |
| 2 | `FILTER_BY_TILE` | BBTools filterbytile | Bad-tile filtering |
| 3 | `BBDUK_TRIM` | BBTools bbduk | Adapter trimming + length filter |
| 4 | `BBDUK_FILTER` | BBTools bbduk | Artifact + PhiX removal |
| 5 | `REMOVE_HUMAN` | BBTools removehuman | Host decontamination (optional) |
| 6 | `FASTQC` | FastQC | Per-sample QC report (optional) |
| 7 | `REMOVE_RRNA` | SortMeRNA | rRNA depletion (optional, on by default) |
| 8 | `BBMAP_INDEX` | BBMap | Reference index build |
| 9 | `MAP_READS_BBMAP` | BBMap + samtools | Per-(sample, reference) alignment |
| 10 | `FEATURECOUNTS` | subread featureCounts | Gene-level counts (when GFF provided) |
| 11 | `MERGE_GENE_COUNTS` | python/pandas | Genes × samples matrix per reference |
| 12 | `VIZ_PREPROCESS` | python | JSONs for the Svelte viz |

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

## Viz

```bash
cd viz
npm install
VIZ_DATA_DIR=/path/to/output/viz npm run dev
# open http://localhost:5173
```

Pages: **Overview** (sample × reference mapping-rate matrix, read-flow funnel), **Samples** (per-sample read flow + per-reference mapping detail), **Expression** (top-variance gene × sample heatmap with log₂ transform; per-reference dropdown), **References** (per-contig read totals).
