# Pipeline Runtime Reference

Resource usage from the Ebb Flow 20260218 run on a 16-CPU / 64 GB local workstation with GPU (NVIDIA A2000).

Runtimes estimated from storeDir file timestamps and process monitoring. Estimates marked with `~` may include minor queuing overhead from CPU contention.

## Run Context

| Metric | Value |
|--------|-------|
| Input | 313 GB, 15 barcodes (nanopore long reads) |
| Assembly | 1.39 Gb, 183,002 contigs, N50 = 4,738 bp, largest = 1.48 Mb |
| Proteins | ~1.1M (Bakta CDS prediction) |
| CPUs | 16 available (Nextflow `fair = true`) |
| RAM | 64 GB |
| GPU | NVIDIA A2000 (used by SemiBin2, LorBin, COMEBin, Tiara) |

## Target HPC Partitions

| Partition | Time Limit | Use For |
|-----------|------------|---------|
| `cpubase_bynode_b1` | 3h | Fast processes: mapping, depths, small tools, QC |
| `cpubase_bynode_b2` | 12h | Medium processes: binning, taxonomy, most annotation |
| `cpubase_bynode_b3` | 1d | Heavy processes: assembly, geNomad, MetaEuk, ML binners |
| `cpubase_bynode_b4` | 3d | Bakta full annotation, eggNOG-mapper, very large assemblies |
| `cpubase_bynode_b5` | 7d | Reserved for exceptional cases |

## GPU vs CPU-Only Note

The runtimes below were measured with GPU acceleration available. On CPU-only HPC nodes, the following deep learning tools will be significantly slower:

| Process | GPU Runtime | CPU-Only Estimate | Notes |
|---------|-------------|-------------------|-------|
| BIN_SEMIBIN2 | ~7h | ~24-48h | PyTorch semi-supervised; consider b4 partition |
| BIN_LORBIN | ~7.5h | ~24-48h | PyTorch; shares model architecture with SemiBin2 |
| BIN_COMEBIN | ~30 min | ~2-4h | Contrastive learning; failed on this dataset |
| TIARA_CLASSIFY | ~73 min | ~3-6h | TensorFlow k-mer DNN |

If running on CPU-only nodes, use the `semibin-cpu.yml` and `comebin-cpu.yml` conda environments and increase time limits accordingly. Consider requesting GPU partitions for these processes if available.

## Process Runtimes and Partition Mapping

### Assembly and Mapping

| Process | CPUs | Mem | Actual Runtime | Partition | Notes |
|---------|------|-----|---------------|-----------|-------|
| CONCAT_READS (x15) | 2 | 4 GB | < 1 min | b1 (3h) | Per-barcode cat |
| ASSEMBLY_FLYE | 16 | 60 GB | ~12h | b3 (1d) | Scales with read volume |
| CALCULATE_TNF | 2 | 4 GB | ~90 min | b1 (3h) | 136 tetranucleotide features, 183K contigs |
| MAP_READS (x14) | 8 | 24 GB | 5-20 min each | b1 (3h) | Per-sample minimap2 + samtools sort |
| CALCULATE_DEPTHS | 8 | 24 GB | ~34 min | b1 (3h) | CoverM on 14 BAMs |

### Binning

| Process | CPUs | Mem | Actual Runtime | Partition | Notes |
|---------|------|-----|---------------|-----------|-------|
| BIN_SEMIBIN2 | 16 | 60 GB | ~7h | b2 (12h) | ML binner, 393 bins; GPU accelerated |
| BIN_LORBIN | 16 | 60 GB | ~7.5h | b2 (12h) | ML binner, 1230 bins; GPU accelerated |
| BIN_COMEBIN | 16 | 60 GB | ~30 min | b2 (12h) | Failed on this dataset (0 bins); GPU accelerated |
| BIN_METABAT2 | 8 | 24 GB | ~2-6h | b2 (12h) | 3298 bins; variable runtime |
| BIN_MAXBIN2 | 8 | 24 GB | ~4h | b2 (12h) | EM algorithm, 278 bins |
| DASTOOL_CONSENSUS | 8 | 24 GB | ~2.5h | b1 (3h) | SCG evaluation + consensus |
| CHECKM2 | 8 | 24 GB | ~2h | b1 (3h) | DIAMOND-based; 8 CPUs sufficient |

### Taxonomy

| Process | CPUs | Mem | Actual Runtime | Partition | Notes |
|---------|------|-----|---------------|-----------|-------|
| KAIJU_CONTIG_CLASSIFY | 16 | 60 GB | ~5h | b2 (12h) | 6-frame translation, greedy mode; ~40 GB DB |
| KRAKEN2_CLASSIFY | 8 | 60 GB | ~9h | b2 (12h) | Loads ~50 GB DB into RAM; maxForks 1 |
| SENDSKETCH_CLASSIFY | 8 | 24 GB | ~70 min | b1 (3h) | Remote GTDB TaxServer |
| RRNA_CLASSIFY | 8 | 24 GB | ~5h | b2 (12h) | barrnap 3x + vsearch 6x + aragorn |

### Annotation

| Process | CPUs | Mem | Actual Runtime | Partition | Notes |
|---------|------|-----|---------------|-----------|-------|
| BAKTA_BASIC | 8 | 24 GB | ~4h | b2 (12h) | Light DB, CDS-only (skips tRNA/rRNA/ncRNA/CRISPR) |
| **BAKTA_EXTRA** | **8** | **24 GB** | **~35.5h** | **b4 (3d)** | **Full DB; DIAMOND PSC (12h) + pseudo-gene (8.7h) + BLASTX (0.6h)** |
| KOFAMSCAN | 8 | 24 GB | ~2-4h (est.) | b2 (12h) | HMM search, ~1.1M proteins |
| EMAPPER | 8 | 24 GB | 5h+ (still running) | b2 (12h) | DIAMOND `--sensitive --iterate` on eggNOG DB |
| DBCAN | 8 | 24 GB | ~1-2h (est.) | b1 (3h) | 3-method CAZyme consensus |
| MERGE_ANNOTATIONS | 2 | 4 GB | ~10s | b1 (3h) | |
| MAP_TO_BINS | 2 | 4 GB | ~10s | b1 (3h) | |
| KEGG_MODULES | 2 | 4 GB | ~5s | b1 (3h) | |
| KEGG_DECODER | 2 | 4 GB | ~20s | b1 (3h) | |
| MINPATH | 2 | 4 GB | ~43 min | b1 (3h) | |

### MGE Detection

| Process | CPUs | Mem | Actual Runtime | Partition | Notes |
|---------|------|-----|---------------|-----------|-------|
| GENOMAD_CLASSIFY | 16 | 60 GB | ~11h | b2 (12h) | Neural network marker-gene detection |
| CHECKV_QUALITY | 8 | 24 GB | ~5-8h | b2 (12h) | AAI + HMM viral QA |
| INTEGRONFINDER | 8 | 24 GB | ~30 min | b1 (3h) | |
| DEFENSEFINDER | 8 | 24 GB | ~5h | b2 (12h) | Parallel HMM on ~1.1M proteins |
| MACSYFINDER | 8 | 24 GB | ~8h | b2 (12h) | TXSScan + CONJScan |

### Eukaryotic

| Process | CPUs | Mem | Actual Runtime | Partition | Notes |
|---------|------|-----|---------------|-----------|-------|
| TIARA_CLASSIFY | 8 | 24 GB | ~73 min | b1 (3h) | DL k-mer classification; GPU accelerated |
| WHOKARYOTE_CLASSIFY | 8 | 24 GB | ~9h | b2 (12h) | Gene-structure RF; needs Prodigal + GFF |
| METAEUK_PREDICT | 16 | 60 GB | ~11h | b2 (12h) | Multi-exon eukaryotic gene prediction |

### Visualization

| Process | CPUs | Mem | Actual Runtime | Partition | Notes |
|---------|------|-----|---------------|-----------|-------|
| VIZ_PREPROCESS | 8 | 24 GB | ~10 min | b1 (3h) | TSV-to-JSON + PCA/t-SNE/UMAP + Vite build |

## Recommended Nextflow Config for HPC

### Resource Labels

```groovy
withLabel: process_low {
    cpus   = 2
    memory = 4.GB
    time   = 3.h        // b1 partition
}
withLabel: process_medium {
    cpus   = 8
    memory = 24.GB
    time   = 12.h       // b2 partition
}
withLabel: process_high {
    cpus   = 16
    memory = 60.GB
    time   = 1.d        // b3 partition
}
withLabel: process_kraken {
    cpus   = 8
    memory = 60.GB
    time   = 12.h       // b2 partition (needs high memory)
}
```

### Per-Process Overrides

```groovy
// BAKTA_EXTRA is the critical bottleneck (17h+ on full DB)
withName: 'BAKTA_EXTRA|BAKTA_FULL' {
    time   = 3.d        // b4 partition
}

// CHECKM2 is overprovisioned at process_high — demote to save slots
withName: 'CHECKM2' {
    cpus   = 8
    memory = 24.GB
    time   = 3.h        // b1 partition; ~2h actual
}

// EMAPPER DIAMOND can run very long on large protein sets
withName: 'EMAPPER' {
    time   = 1.d        // b3 partition; 5h+ observed, may reach 12h+
}

// CPU-only overrides for ML binners (if no GPU partition available)
withName: 'BIN_SEMIBIN2|BIN_LORBIN' {
    time   = 3.d        // b4 partition when running CPU-only
}
```

## Critical Path Analysis

On a multi-node HPC where all processes can run in parallel:

```
Assembly (12h) ──┬── Bakta Basic (4h) ──┬── DefenseFinder (5h)
                 │                      ├── MacSyFinder (8h)
                 │                      ├── Whokaryote (9h)
                 │                      ├── eggNOG-mapper (12h+)
                 │                      ├── KofamScan (4h)
                 │                      └── KEGG modules (<1h)
                 │
                 ├── Bakta Extra (24h+) ── IslandPath (<1h)
                 │
                 ├── Binning (7h) ── DAS Tool (2.5h) ── CheckM2 (2h)
                 │
                 ├── geNomad (11h) ── CheckV (8h)
                 │
                 ├── Kraken2 (9h)
                 ├── Kaiju (5h)
                 ├── MetaEuk (11h)
                 └── rRNA (5h)
```

**Minimum wall time** (unlimited parallelism): ~48h
- Critical path: Assembly (12h) → Bakta Extra (35.5h)

**Without Bakta Extra**: ~29h
- Critical path: Assembly (12h) → Bakta Basic (4h) → eggNOG-mapper (12h+)

## Scaling Notes

- Runtimes scale roughly linearly with assembly size (contigs x length)
- This assembly (1.39 Gb, 183K contigs) is from 15 nanopore barcodes over ~3 weeks
- Smaller runs (1-4 barcodes, <200 Mb assembly) typically complete within b2 limits
- Memory is dominated by database loading: Kraken2 ~50 GB, Kaiju ~40 GB, Bakta full ~30 GB
- DIAMOND steps (Bakta pseudo-genes, eggNOG-mapper, MarFERReT) are CPU-bound and benefit from more threads
