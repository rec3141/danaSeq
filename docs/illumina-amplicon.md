# danaSeq illumina_amplicon

**Amplicon sequencing analysis pipeline — from raw reads to interactive visualization**

danaSeq illumina_amplicon is a Nextflow DSL2 pipeline for amplicon sequencing analysis. It takes
demultiplexed paired-end FASTQ files and produces ASV tables, multi-database taxonomy,
phylogenies, ordinations, and correlation networks.

## Architecture

Denoising is [papa2](https://rec3141.github.io/papa2/) (a DADA2 port: filter,
dereplicate, denoise, merge, chimera removal, taxonomy). Everything downstream —
QC filtering, ordination, networks, visualization — is carried by this stage's
own scripts in `bin/`.

```
Raw FASTQ files
    │
    ├── REMOVE_PRIMERS (cutadapt)
    ├── DADA2_FILTER_TRIM (papa2)
    ├── DADA2_LEARN_ERRORS (papa2, per-plate)
    ├── DADA2_DENOISE (papa2, per-plate)
    ├── MERGE_SEQTABS
    ├── REMOVE_CHIMERAS (papa2)
    ├── FILTER_SEQTAB (bin/)
    │
    ├── ASSIGN_TAXONOMY (papa2, parallel per DB)
    ├── BUILD_PHYLOGENY (bin/ + MAFFT)
    ├── RENORMALIZE (bin/)
    │
    ├── ORDINATE (bin/, t-SNE/PCA)
    ├── NETWORK (bin/, SparCC)
    └── EXPORT_VIZ (bin/, JSON → Svelte)
```

## Related

| Project | Purpose |
|---|---|
| [papa2](https://github.com/rec3141/papa2) | DADA2 denoising (pinned GitHub release, byte-identical to R dada2) |
| [danaSeq](https://github.com/rec3141/danaSeq) | This pipeline's home, alongside the metagenomics stages |
| [microscape.app](https://microscape.app) | Hosting and sharing the interactive viz |

The `microscape` (Python) and `microscapeR` (R) packages are retired. They
mirrored this pipeline's downstream steps in two more places, which is how the
truncation logic came to disagree with itself; that code now lives here in
`bin/`, and the archived packages are reference snapshots only.
