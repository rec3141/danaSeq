# Pipeline Stages

## Stage A: Preprocessing and DADA2

| Step | Process | Tool | Description |
|------|---------|------|-------------|
| 1 | `DEMULTIPLEX` | Mr_Demuxy | Optional inner-barcode demultiplexing |
| 2 | `REMOVE_PRIMERS` | cutadapt | Primer trimming (auto-selects by 16S/18S/ITS prefix) |
| 3 | `DADA2_FILTER_TRIM` | papa2 | Per-sample quality filtering (maxEE, truncQ, PhiX) |
| 4 | `DADA2_LEARN_ERRORS` | papa2 | Per-plate error model learning (LOESS) |
| 5 | `DADA2_DENOISE` | papa2 | Denoising, pair merging, per-plate chimera removal |
| 6 | `MERGE_SEQTABS` | papa2 | Merge per-plate tables (long-format) |
| 7 | `REMOVE_CHIMERAS` | papa2 | Consensus chimera removal on merged data |
| 8 | `FILTER_SEQTAB` | `filter_seqtab.py` | Length, prevalence, abundance, and depth filtering |

!!! note "Per-plate processing"
    Steps 3–5 run independently per plate (samples sharing PCR history).
    This allows per-plate error learning while scaling to hundreds of plates.

---

## Stage B: Taxonomy, Phylogeny, and Normalization

| Step | Process | Tool | Description |
|------|---------|------|-------------|
| 9 | `ASSIGN_TAXONOMY` | papa2 | Naive Bayesian classification (parallel per ref DB) |
| 10 | `BUILD_PHYLOGENY` | `build_phylogeny.py` + MAFFT | MSA + NJ tree (optional, `--run_phylogeny`) |
| 11 | `RENORMALIZE` | `renormalize.py` | Group ASVs by taxonomy, normalize within groups |

!!! tip "Multiple databases"
    Supply multiple reference databases with `--ref_databases` and taxonomy
    is assigned against each in parallel. Results are merged in the output.

---

## Stage C: Ordination, Networks, and Visualization

| Step | Process | Tool | Description |
|------|---------|------|-------------|
| 12 | `LOAD_METADATA` | `load_metadata.py` | Sample metadata integration (MIMARKS) |
| 13 | `CLUSTER_TSNE` | `cluster_tsne.py` | Bray-Curtis + t-SNE ordination (samples and ASVs) |
| 14 | `NETWORK_SPARCC` | `network_sparcc.py` | SparCC-style CLR correlation networks |
| 15 | `BUILD_VIZ` | `build_viz.py` | JSON export for Svelte web viewer |

---

## Data Flow

The pipeline uses **long-format DataFrames** as its canonical representation
from the merge step onward:

```
sample          sequence            count
plate1_A01      TACGGAGGATGCGA...   1523
plate1_A01      TACGGAGGATCCGA...   847
plate1_A02      TACGGAGGATGCGA...   2041
```

This avoids materializing a dense matrix (samples × ASVs), which can exceed
memory for large datasets (4K+ samples, 100K+ ASVs).

---

## Profiles

| Profile | Description |
|---------|-------------|
| `conda` | Auto-create conda environments from `envs/*.yml` |
| `docker` | Use Docker container |
| `singularity` | Use Singularity/Apptainer container |
| `slurm` | Submit jobs to SLURM cluster |
| `test` | Reduced resources for CI/testing |
