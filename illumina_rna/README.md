# danaSeq · illumina_rna

Paired-end Illumina/BGI RNA-seq for metatranscriptomes, mapped against
references from other danaSeq pipelines (or any FASTA + GFF). Sized for
co-assemblies: an 11 Gbp, 1.5M-contig Flye assembly with 13M Bakta genes and
130 transcriptomes runs overnight on a handful of nodes.

## Stages

1. **QC** — one streaming pass per sample: bbduk adapter trim → bbduk
   artifacts/phiX/entropy → bbduk rRNA split (k=31 against e.g. SILVA SSU+LSU).
   Read counts come from each stage's summary, so nothing is re-read.
2. **Human removal** — strobealign against a masked human FASTA; a pair goes
   when a mate aligns over ≥80% of its length at ≥95% identity (the
   stringency of BBTools `removehuman.sh`).
3. **Mapping** — strobealign (default; the index is built once per reference)
   or BBmap. Name-grouped BAMs, no coordinate sort.
4. **Quantification** — featureCounts on genes (from the reference GFF; CDS,
   rRNA, tRNA, tmRNA, ncRNA) and on whole contigs.
5. **Summarize** — genes/contigs × samples sparse matrices, per-MAG and
   per-function (KO / COG / CAZy) TPM, read fate per sample, viz JSONs.

## Chaining to other danaSeq outputs

| Reference source | Flag | What is used |
|---|---|---|
| `mag_analysis` | `--from_mag DIR` | `assembly/assembly.fasta`, Bakta GFF (`annotation/bakta/{extra,basic}`), bins (`binning/<bin_set>/`, default `dastool`), CheckM2, GTDB-Tk, `metabolism/merged/merged_annotations.tsv` |
| `nanopore_assembly` / `illumina_assembly` | `--from_assembly DIR` | the assembly FASTA(s); contig-level counts |
| anything | `--references DIR` | `<name>.fasta`, optional `<name>.gff` and `<name>.bins.tsv` (contig⇥bin) |

## Quick start

```bash
./run-illumina-rna.sh --container \
    --input /path/to/reads \
    --from_mag /path/to/mag_analysis_out \
    --human_ref /path/to/hg19_masked.fa.gz \
    --rrna_refs /path/to/silva_ssu_lsu_dir \
    --outdir /path/to/output
```

Reads are found recursively under `--input` (`*_R1_*.fastq.gz`, `*_1.fq.gz`,
`*_1.fastq.gz`), or listed in `--samplesheet` (`sample,r1,r2`).

In container mode every path argument is bound into the container at the same
path; add `--bind PATH` for anything else the inputs point into (symlink
targets, samplesheet files). The pipeline code runs from this checkout.

## Large cohorts

- `--max_read_pairs N` uses the first N pairs of each sample.
- `--prepare_only` builds the shared references (strobealign indexes, gene
  SAFs, rRNA k-mers) and stops.
- `--chunk I/N --summarize false` processes one of N fixed sample chunks; run
  the chunks as separate jobs with one shared `--store_dir`, then one run
  without `--chunk` merges everything (every per-sample task is a store hit).
- Inside a SLURM job the local executor is capped at the job's CPUs and memory.

## Outputs

```
<outdir or store_dir>/
├── references/<ref>/        index, genes/contigs SAF, genes.tsv.gz (annotation), contig lengths
├── preprocess/<sample>/     <sample>.qc.fq.gz, .nohuman.fq.gz, QC logs and counts
├── mapping/<sample>/<ref>/  name-grouped BAM, flagstat, strobealign log
├── quantify/<sample>/<ref>/ gene and contig counts, featureCounts summary
<outdir>/
├── expression/<ref>/        gene_counts.mtx.gz, contig_counts.mtx.gz (+ samples.tsv),
│                            sample_summary.tsv, bin_counts.tsv, bin_tpm.tsv,
│                            ko_tpm.tsv / cog_tpm.tsv / cazy_tpm.tsv
└── viz/                     index.json + <ref>.{samples,bins,functions,top_genes}.json.gz
```

Matrix rows follow `references/<ref>/<ref>.genes.tsv.gz` (genes) and
`<ref>.contig_lengths.tsv` (contigs); a dense `gene_counts.tsv.gz` is also
written when the matrix is small.

## Viz

```bash
cd viz
npm install
VIZ_DATA_DIR=/path/to/output/viz npm run dev     # http://localhost:5173
./deploy.sh --preprocess-dir /path/to/output/viz --slug my-run --name "My run"
```

Pages: Overview (where each sample's reads went), MAGs (per-MAG TPM with
CheckM2 / GTDB-Tk), Functions (KO / COG / CAZy), Genes (most expressed genes),
Samples (QC and mapping table). `--viz_mag_url` links the header to the
mag_analysis viz of the same reference.
