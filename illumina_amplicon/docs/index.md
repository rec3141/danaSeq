# microscape-nf

**Amplicon sequencing analysis pipeline — from raw reads to interactive visualization**

microscape-nf is a Nextflow DSL2 pipeline for amplicon sequencing analysis. It takes
demultiplexed paired-end FASTQ files and produces ASV tables, multi-database taxonomy,
phylogenies, ordinations, and correlation networks.

## Architecture

The pipeline uses two bioconda packages for all computation:

- **[papa2](https://rec3141.github.io/papa2/)** — DADA2 denoising (filter, dereplicate, denoise, merge, chimera removal, taxonomy)
- **[microscape](https://rec3141.github.io/microscape/)** — Downstream analysis (QC filtering, ordination, networks, visualization)

An R alternative is available via `--lang R`, which uses Bioconductor's dada2 and
[microscapeR](https://github.com/rec3141/microscapeR).

```
Raw FASTQ files
    │
    ├── REMOVE_PRIMERS (cutadapt)
    ├── DADA2_FILTER_TRIM (papa2)
    ├── DADA2_LEARN_ERRORS (papa2, per-plate)
    ├── DADA2_DENOISE (papa2, per-plate)
    ├── MERGE_SEQTABS
    ├── REMOVE_CHIMERAS (papa2)
    ├── FILTER_SEQTAB (microscape)
    │
    ├── ASSIGN_TAXONOMY (papa2, parallel per DB)
    ├── BUILD_PHYLOGENY (microscape + MAFFT)
    ├── RENORMALIZE (microscape)
    │
    ├── ORDINATE (microscape, t-SNE/PCA)
    ├── NETWORK (microscape, SparCC)
    └── EXPORT_VIZ (microscape, JSON → Svelte)
```

## Related Packages

| Package | Language | Channel | Purpose |
|---|---|---|---|
| [papa2](https://github.com/rec3141/papa2) | Python | bioconda | DADA2 denoising |
| [microscape](https://github.com/rec3141/microscape) | Python | bioconda | Downstream analysis |
| [microscapeR](https://github.com/rec3141/microscapeR) | R | Bioconductor | R companion |
| [microscape-nf](https://github.com/rec3141/microscape-nf) | Nextflow | GitHub | This pipeline |
