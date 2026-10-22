# MAG Analysis Pipeline

Technology-agnostic downstream analysis for metagenome-assembled genomes (MAGs). Accepts a pre-computed assembly + depth table from any assembler and runs binning, annotation, taxonomy, metabolic profiling, mobile genetic element detection, eukaryotic analysis, ecosystem services mapping, phylogenetics, and interactive visualization.

## Quick Start

```bash
cd mag_analysis
./install.sh && ./install.sh --check

# Basic run
./run-mag-analysis.sh \
    --assembly /path/to/assembly.fasta \
    --depths /path/to/depths.txt \
    --outdir /path/to/output \
    --annotator bakta --db_dir /path/to/databases

# All modules
./run-mag-analysis.sh \
    --assembly /path/to/assembly.fasta \
    --depths /path/to/depths.txt \
    --bam_dir /path/to/mapping/ \
    --outdir /path/to/output \
    --all --db_dir /path/to/databases

# Apptainer (HPC)
./run-mag-analysis.sh --apptainer \
    --assembly /path/to/assembly.fasta \
    --depths /path/to/depths.txt \
    --outdir /path/to/output \
    --all --db_dir /path/to/databases
```

## Pipeline Stages

### Binning

| Step | Process | Tool | Description |
|------|---------|------|-------------|
| 1a | `BIN_SEMIBIN2` | SemiBin2 | Self-supervised contig binning (needs BAMs) |
| 1b | `BIN_METABAT2` | MetaBAT2 | Depth + TNF binning (always runs) |
| 1c | `BIN_MAXBIN2` | MaxBin2 | EM-based abundance binning |
| 1d | `BIN_LORBIN` | LorBin | Long-read-aware binning (needs BAMs) |
| 1e | `BIN_COMEBIN` | COMEBin | Contrastive learning binning (needs BAMs) |
| 1f | `BIN_VAMB` | VAMB | Variational autoencoder binning |
| 1g | `BIN_VAMB_TAX` | VAMB | Taxonomy-guided variational autoencoder binning |
| 2a | `DASTOOL_CONSENSUS` | DAS Tool | Score-based consensus of multiple binners |
| 2b | `BINETTE_CONSENSUS` | Binette | CheckM2-guided bin refinement |
| 2c | `MAGSCOT_CONSENSUS` | MAGScoT | Marker-gene-based consensus |
| 3 | `CHECKM2` | CheckM2 | Completeness and contamination assessment |
| 3b | `DEREP_GALAH` | galah + skani | Species-level dereplication of the consensus bins, ranked by CheckM2 |

### Annotation

| Step | Process | Tool | Description |
|------|---------|------|-------------|
| 4a | `PROKKA_ANNOTATE` | Prokka | ORF prediction + functional annotation |
| 4b | `BAKTA_BASIC` | Bakta | CDS annotation (lightweight) |
| 4c | `BAKTA_EXTRA` | Bakta | Full annotation (ncRNA, tRNA, CRISPR, sORFs) |

### Taxonomy

| Step | Process | Tool | Description |
|------|---------|------|-------------|
| 5a | `KAIJU_CONTIG_CLASSIFY` | Kaiju | Protein-level contig classification |
| 5b | `KAIJU_CLASSIFY` | Kaiju | Protein-level MAG classification |
| 5c | `KRAKEN2_CLASSIFY` | Kraken2 | k-mer-based classification |
| 5d | `SENDSKETCH_CLASSIFY` | BBTools sendsketch | MinHash taxonomy (GTDB) |
| 5e | `RNA_CLASSIFY` | barrnap + BLAST | rRNA extraction + SILVA classification |
| 5f | `GTDBTK_CLASSIFY` | GTDB-Tk | Genome-based taxonomy (requires ~120 GB RAM) |

### Metabolism

| Step | Process | Tool | Description |
|------|---------|------|-------------|
| 6a | `KOFAMSCAN` | KofamScan | KEGG ortholog assignment via HMM profiles |
| 6b | `EMAPPER` | eggNOG-mapper | Ortholog annotation (COG, GO, KEGG, CAZy) |
| 6c | `DBCAN` | dbCAN | CAZyme annotation via HMM + diamond + HMMER |
| 6d | `MERGE_ANNOTATIONS` | R | Merge KofamScan + eggNOG + dbCAN per gene |
| 6e | `MAP_TO_BINS` | R | Assign annotations to MAGs |
| 6f | `KEGG_MODULES` | R | KEGG module completeness per MAG |
| 6g | `MINPATH` | MinPath | Pathway parsimony analysis |
| 6h | `KEGG_DECODER` | KEGG-Decoder | Metabolic pathway heatmaps |
| 6i | `ANTISMASH` | antiSMASH | Biosynthetic gene cluster detection |
| 6j | `ECOSSDB_MAP` | R | Ecosystem services mapping (CICES 5.2) |
| 6k | `ECOSSDB_SCORE` | R | Ecosystem service scoring per MAG |
| 6l | `ECOSSDB_SDG` | R | UN Sustainable Development Goal mapping |
| 6m | `ECOSSDB_VIZ` | R | Ecosystem services visualization |

### Mobile Genetic Elements

| Step | Process | Tool | Description |
|------|---------|------|-------------|
| 7a | `GENOMAD_CLASSIFY` | geNomad | Virus/plasmid identification |
| 7b | `CHECKV_QUALITY` | CheckV | Viral genome quality assessment |
| 7c | `INTEGRONFINDER` | IntegronFinder | Integron detection |
| 7d | `ISLANDPATH_DIMOB` | IslandPath-DIMOB | Genomic island prediction |
| 7e | `MACSYFINDER` | MacSyFinder | Secretion system detection |
| 7f | `DEFENSEFINDER` | DefenseFinder | Anti-phage defense system detection |

### Eukaryotic Analysis

| Step | Process | Tool | Description |
|------|---------|------|-------------|
| 8a | `TIARA_CLASSIFY` | Tiara | Domain-level classification (eukaryote/prokaryote) |
| 8b | `WHOKARYOTE_CLASSIFY` | Whokaryote | Eukaryote/prokaryote classifier |
| 8c | `METAEUK_PREDICT` | MetaEuk | Eukaryotic gene prediction |
| 8d | `MARFERRET_CLASSIFY` | MarFERReT | Marine eukaryote taxonomy |

### Other

| Step | Process | Tool | Description |
|------|---------|------|-------------|
| 9 | `CALCULATE_GENE_DEPTHS` | samtools + R | Per-gene coverage from BAMs |
| 10 | `VIZ_PREPROCESS` | R + Svelte | Interactive dashboard generation |

## Parameters

### Required Inputs

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--assembly` | (required) | Assembly FASTA file |
| `--depths` | (required) | Depth matrix (MetaBAT2 format) |
| `--bam_dir` | (optional) | Directory with sorted BAMs for BAM-based binners |

### Binning

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--run_semibin` | `false` | Include SemiBin2 (needs `--bam_dir`) |
| `--run_maxbin` | `false` | Include MaxBin2 |
| `--run_lorbin` | `false` | Include LorBin (needs `--bam_dir`) |
| `--run_comebin` | `false` | Include COMEBin (needs `--bam_dir`). Without a visible GPU it is skipped, unless `--comebin_on_cpu` is set or a COMEBin result is already in `--store_dir` |
| `--comebin_on_cpu` | `false` | Run COMEBin even without a GPU; on a co-assembly that takes days |
| `--binners_only` | `false` | Run the selected binners and stop, with no DAS Tool, Binette, MAGScoT or CheckM2. Lets a GPU node run the binners ahead of a full run sharing `--store_dir` |
| `--gpu` | `auto` | Train SemiBin2, LorBin, COMEBin and VAMB on a GPU: `auto` uses one when the tool's PyTorch sees it, `true` also warns when none is usable, `false` forces CPU. Needs the GPU image (`--gpu` in `run-mag-analysis.sh`) or a CUDA conda install |
| `--run_vamb` | `false` | Include VAMB |
| `--run_vamb_tax` | `false` | Include taxonomy-guided VAMB |
| `--run_binette` | `false` | Run Binette consensus (needs `--checkm2_db`) |
| `--run_magscot` | `false` | Run MAGScoT consensus |
| `--metabat_min_cls` | `50000` | MetaBAT2 minimum cluster size |
| `--lorbin_min_length` | `80000` | LorBin minimum contig length |

### Dereplication

`--run_derep` clusters the consensus bins (DAS Tool, plus Binette and MAGScoT
when they run) at species-level ANI with galah, and keeps one representative
per cluster: the genome with the highest completeness − 5 × contamination, less
small penalties for contig count and ambiguous bases (galah's
`Parks2020_reduced`). It needs `--checkm2_db`, and clusters only bins that
CheckM2 scored. Within one run, it picks the best of the refiners' versions of
each genome; it does not dereplicate across runs or co-assemblies.

`binning/derep/derep_clusters.tsv` lists every consensus bin. Its `status` is
`representative`, `member` (with its skani ANI and aligned fractions against
the representative; `NA` if the pair falls below skani's 80% ANI screen),
`failed_quality` (outside the completeness/contamination thresholds), or
`no_checkm2` (no CheckM2 estimate, e.g. a bin above `--checkm2_max_bin_bp`).

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--run_derep` | `false` | Dereplicate the consensus bins with galah (needs `--checkm2_db`) |
| `--derep_ani` | `95` | ANI (%) at which genomes are clustered |
| `--derep_min_af` | `15` | Minimum aligned fraction (%) for two genomes to cluster |
| `--derep_min_completeness` | `50` | CheckM2 completeness (%) below which a bin is not clustered |
| `--derep_max_contamination` | `10` | CheckM2 contamination (%) above which a bin is not clustered |
| `--derep_min_size_bp` | `500000` | Bins shorter than this (bp) are not clustered, so a fragment cannot become a representative |

### Annotation

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--annotator` | `bakta` | Gene annotator: `prokka`, `bakta`, or `none` |
| `--bakta_db` | (required if bakta) | Path to Bakta database |
| `--bakta_extra` | `false` | Run full Bakta annotation |

### Taxonomy

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--run_kaiju` | `false` | Run Kaiju classification |
| `--kaiju_db` | (required if kaiju) | Kaiju database path |
| `--run_kraken2` | `false` | Run Kraken2 classification |
| `--kraken2_db` | (required if kraken2) | Kraken2 database path |
| `--run_sendsketch` | `false` | Run sendsketch GTDB taxonomy |
| `--run_rrna` | `false` | Run rRNA classification (SILVA) |
| `--silva_ssu_db` | (required if rrna) | SILVA SSU database path |
| `--rrna_min_identity` | `0.80` | Minimum identity for rRNA BLAST |
| `--run_gtdbtk` | `false` | Run GTDB-Tk classification |
| `--gtdbtk_db` | (required if gtdbtk) | GTDB-Tk database path |

### Metabolism

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--run_metabolism` | `false` | Enable metabolic profiling |
| `--kofam_db` | (required if metabolism) | KofamScan database path |
| `--eggnog_db` | (required if metabolism) | eggNOG database path |
| `--dbcan_db` | (required if metabolism) | dbCAN database path |
| `--emapper_batch_size` | `50000` | eggNOG-mapper batch size (proteins) |
| `--run_antismash` | `false` | Run antiSMASH BGC detection |
| `--antismash_db` | (required if antismash) | antiSMASH database path |
| `--run_ecossdb` | `true` | Enable ecosystem services mapping |

### Mobile Genetic Elements

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--run_genomad` | `false` | Run geNomad virus/plasmid detection |
| `--genomad_db` | (required if genomad) | geNomad database path |
| `--run_checkv` | `false` | Run CheckV quality assessment |
| `--checkv_db` | (required if checkv) | CheckV database path |
| `--run_integronfinder` | `false` | Run IntegronFinder |
| `--integron_chunk_size` | `2000` | Contigs per IntegronFinder task. IntegronFinder keeps two files per contig until its input is done, so this, times `--integron_max_forks`, bounds the files it holds at once |
| `--integron_max_forks` | `16` | IntegronFinder chunks run at once |
| `--run_islandpath` | `false` | Run IslandPath-DIMOB |
| `--run_macsyfinder` | `false` | Run MacSyFinder |
| `--run_defensefinder` | `false` | Run DefenseFinder |

### Eukaryotic Analysis

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--run_eukaryotic` | `false` | Enable eukaryotic classification |
| `--tiara_min_len` | `3000` | Tiara minimum contig length |
| `--whokaryote_min_len` | `5000` | Whokaryote minimum contig length |
| `--run_metaeuk` | `false` | Run MetaEuk gene prediction |
| `--metaeuk_db` | (required if metaeuk) | MetaEuk database path |
| `--run_marferret` | `false` | Run MarFERReT taxonomy |
| `--marferret_db` | (required if marferret) | MarFERReT database path |

### Convenience

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--all` | `false` | Enable all analysis modules |
| `--db_dir` | (none) | Auto-resolve database paths from standard layout |
| `--store_dir` | (none) | Persistent cache directory (storeDir) |
| `--run_viz` | `false` | Build interactive Svelte dashboard |
| `--viz_port` | `5174` | Dashboard dev server port |

## Outputs

```
results/
├── binning/
│   ├── semibin/                SemiBin2 bins
│   ├── metabat/                MetaBAT2 bins
│   ├── maxbin/                 MaxBin2 bins
│   ├── lorbin/                 LorBin bins
│   ├── comebin/                COMEBin bins
│   ├── vamb/                   VAMB bins
│   ├── dastool/                DAS Tool consensus bins
│   ├── binette/                Binette refined bins
│   ├── magscot/                MAGScoT consensus bins
│   ├── checkm2/                Quality assessment (completeness/contamination)
│   └── derep/                  Dereplicated consensus bins (--run_derep):
│                                 representatives/ (one FASTA per species cluster) and
│                                 derep_clusters.tsv (genome, representative, ani,
│                                 af_genome, af_representative, completeness,
│                                 contamination, status)
├── annotation/
│   ├── prokka/ or bakta/       Gene annotations (GFF, FAA, FFN)
│   └── bakta_extra/            Full Bakta output (if enabled)
├── taxonomy/
│   ├── kaiju/                  Protein-level classification
│   ├── kraken2/                k-mer classification
│   ├── sendsketch/             GTDB MinHash taxonomy
│   ├── rrna/                   rRNA SILVA classification
│   └── gtdbtk/                 GTDB-Tk genome taxonomy
├── metabolism/
│   ├── kofamscan/              KEGG ortholog assignments
│   ├── emapper/                eggNOG annotations
│   ├── dbcan/                  CAZyme annotations
│   ├── merged/                 Merged annotation table
│   ├── per_mag/                Per-MAG annotation summaries
│   ├── modules/                KEGG module completeness
│   ├── minpath/                MinPath pathway analysis
│   ├── kegg_decoder/           Pathway heatmaps
│   ├── antismash/              Biosynthetic gene clusters
│   └── ecossdb/                Ecosystem services (CICES + SDG)
├── mge/
│   ├── genomad/                Virus/plasmid predictions
│   ├── checkv/                 Viral quality assessment
│   ├── integrons/              Integron predictions
│   ├── islandpath/             Genomic islands
│   ├── macsyfinder/            Secretion systems
│   └── defensefinder/          Defense systems
├── eukaryotic/
│   ├── tiara/                  Domain classification
│   ├── whokaryote/             Eukaryote/prokaryote predictions
│   ├── metaeuk/                Eukaryotic gene predictions
│   └── marferret/              Marine eukaryote taxonomy
├── viz/                        Interactive Svelte dashboard
└── pipeline_info/              Nextflow reports (timeline, trace, DAG)
```

The contig explorer's t-SNE uses [fastsne](https://github.com/rec3141/fastsne) when it is installed (the
containers ship it; the `-gpu` image also has the CUDA build): seconds instead of minutes, and `--dedup 2`
keeps over-represented genomes from stretching into worms. The viz preprocess picks the engine with
`--tsne-engine auto|gpu|cpu|sklearn` (auto: GPU, then CPU, then sklearn) and also writes
`viz/fastsne/viewer.html?data=contigs.fsv`, a WebGL view of the t-SNE, UMAP and PCA maps coloured by bin,
taxonomy, GC, depth and length (`--no-fastsne-viewer` to skip).

## Profiles

| Profile | Use case |
|---------|----------|
| `standard` | Local execution (default) |
| `test` | Small test data, reduced resources (4 CPUs, 8 GB) |

## Resource Requirements

| Component | CPUs | RAM | Notes |
|-----------|------|-----|-------|
| MetaBAT2 / MaxBin2 | 16 | 60 GB | Default `process_high` label |
| SemiBin2 | 8 | 16 GB | GPU optional |
| COMEBin | 8 | 16 GB | GPU; skipped on CPU unless `--comebin_on_cpu` |
| Prokka / Bakta | 16 | 60 GB | Per-bin annotation |
| Bakta (extra) | 16 | 24 GB | Full annotation mode |
| eggNOG-mapper | 16 | 48 GB | Batched protein annotation |
| MetaEuk | 16 | 32 GB | Eukaryotic gene prediction |
| GTDB-Tk | 16 | 120 GB | Most memory-intensive step |
| Kraken2 | 8 | 24 GB | Depends on database |
| CheckM2 | 16 | 60 GB | Quality assessment |

## Databases

Download all databases with the interactive menu:

```bash
./download-databases.sh --all --dir /path/to/databases
```

Or download individually:

```bash
./download-databases.sh --checkm2   # CheckM2 (~3 GB)
./download-databases.sh --genomad   # geNomad (~3 GB)
./download-databases.sh --checkv    # CheckV (~2 GB)
./download-databases.sh --kaiju     # Kaiju nr_euk (~60 GB)
./download-databases.sh --gtdbtk    # GTDB-Tk (~85 GB)
```

Use `--db_dir` to auto-resolve paths from a standard layout:

```bash
./run-mag-analysis.sh --db_dir /path/to/databases ...
```

## Design Notes

- **Seven-binner consensus.** SemiBin2, MetaBAT2, MaxBin2, LorBin, COMEBin, VAMB, VAMB-tax feed into DAS Tool + Binette + MAGScoT for consensus refinement.
- **Modular activation.** Every analysis module is off by default; use `--run_*` flags or `--all` to enable selectively.
- **Technology-agnostic.** Accepts assemblies from Nanopore, Illumina, HiFi, or any external assembler. Only requires FASTA + MetaBAT2-format depth table.
- **Persistent caching.** Use `--store_dir` to skip completed processes across runs, even after `work/` cleanup.
- **MIMAG quality standards.** CheckM2 classifies MAGs as high (>90% completeness, <5% contamination), medium (>50%, <10%), or low quality.
