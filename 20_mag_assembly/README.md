# MAG Assembly Pipeline

Metagenome-assembled genome (MAG) reconstruction from Oxford Nanopore long reads. Co-assembles reads with Flye, maps back with minimap2, runs five binning algorithms (SemiBin2, MetaBAT2, MaxBin2, LorBin, COMEBin), and integrates results with DAS Tool consensus. Includes mobile genetic element detection, gene annotation, contig-level taxonomy, anti-phage defense system detection, and metabolic profiling (KofamScan, eggNOG-mapper, dbCAN3).

## Quick Start

```bash
cd nextflow

# Install conda environments (~15 min first time)
./install.sh
./install.sh --check

# Download databases (interactive menu or specify --genomad, --checkv, --checkm2, --kaiju, etc.)
./download-databases.sh

# Run (local conda, handles activation automatically)
./run-mag.sh --input /path/to/reads --outdir /path/to/output

# Or with Docker
docker build -t danaseq-mag .
./run-mag.sh --docker --input /path/to/reads --outdir /path/to/output

# Show all options
./run-mag.sh --help
```

### Kitchen sink — all options with defaults

```bash
cd nextflow
./run-mag.sh --input /data/reads --outdir /data/output \
    --dedupe \
    --filtlong_size 40000000000 \
    --min_overlap 1000 \
    --run_maxbin true \
    --run_lorbin true \
    --run_comebin true \
    --lorbin_min_length 80000 \
    --metabat_min_cls 50000 \
    --checkm2_db /path/to/checkm2_db \
    --genomad_db /path/to/genomad_db \
    --checkv_db /path/to/checkv_db \
    --kaiju_db /path/to/kaiju_db \
    --macsyfinder_models /path/to/macsyfinder_models \
    --defensefinder_models /path/to/defensefinder_models \
    --run_metabolism true \
    --kofam_db /path/to/kofam_db \
    --eggnog_db /path/to/eggnog_db \
    --dbcan_db /path/to/dbcan_db \
    --assembly_cpus 24 \
    --assembly_memory '64 GB'
```

## Pipeline Overview

```
Sample FASTQs (N files)
         |  collect()
   ASSEMBLY_FLYE           All reads -> 1 co-assembly (Flye --meta)
         |
   MAP_READS (xN)          Per-sample alignment (minimap2, -F 0x904)
         |  collect()
   CALCULATE_DEPTHS         Coverage depth table (CoverM)
         |
    +---------+---------+---------+---------+
    |         |         |         |         |
 SemiBin2  MetaBAT2  MaxBin2  LorBin   COMEBin   Five binners
    |         |         |         |         |
    +---------+---------+---------+---------+
         |
   DASTOOL_CONSENSUS        Best bin per contig
         |
   CHECKM2                  Quality assessment (optional, needs --checkm2_db)

  --- parallel annotation & MGE branches from assembly ---

   PROKKA_ANNOTATE          Gene annotation (CDS, rRNA, tRNA)
      |
      +-- KAIJU_CLASSIFY    Protein-level taxonomy (RefSeq)
      +-- ISLANDPATH_DIMOB  Genomic island detection (dinucleotide bias)
      +-- MACSYFINDER       Secretion + conjugation systems
      +-- DEFENSEFINDER     Anti-phage defense systems (CRISPR, R-M, BREX)
      +-- KOFAMSCAN         KEGG Orthology (adaptive HMM thresholds)
      +-- EMAPPER           COG/GO/EC/KEGG/Pfam (eggNOG-mapper)
      +-- DBCAN             CAZyme annotation (3-method consensus)
      |   |   |
      +---+---+
      MERGE_ANNOTATIONS     Unified per-protein annotation table
      |
      MAP_TO_BINS           Per-MAG annotation tables (via contig2bin)
      |
      KEGG_MODULES          Module completeness scoring + heatmap

   GENOMAD_CLASSIFY          Virus + plasmid + provirus detection
      |
      +-- CHECKV_QUALITY     Viral genome quality assessment

   INTEGRONFINDER           Integron + gene cassette detection
   CALCULATE_TNF            Tetranucleotide frequency profiles
```

## Output

```
results/
├── assembly/
│   ├── assembly.fasta            Co-assembly
│   └── tnf.tsv                   Tetranucleotide frequencies (136 features)
├── mapping/
│   ├── *.sorted.bam              Per-sample alignments
│   ├── *.sorted.bam.bai          BAM indices
│   └── depths.txt                CoverM depth table
├── binning/
│   ├── semibin/contig_bins.tsv   SemiBin2 assignments
│   ├── metabat/contig_bins.tsv   MetaBAT2 assignments
│   ├── maxbin/contig_bins.tsv    MaxBin2 assignments
│   ├── lorbin/contig_bins.tsv    LorBin assignments
│   ├── comebin/contig_bins.tsv   COMEBin assignments
│   ├── dastool/
│   │   ├── bins/*.fa             Final consensus MAG FASTAs
│   │   ├── contig2bin.tsv        Contig-to-bin assignments
│   │   ├── allbins.fa            All bins concatenated
│   │   ├── bin_quality.tsv       DAS Tool SCG-based quality scores
│   │   └── summary.tsv           Consensus winners with scores
│   └── checkm2/
│       └── quality_report.tsv    CheckM2 completeness/contamination (if --checkm2_db)
├── annotation/
│   └── prokka/                   Prokka gene annotation (if --run_prokka)
│       ├── *.gff                 GFF3 annotations
│       ├── *.gbk                 GenBank format
│       ├── *.faa                 Protein sequences
│       ├── *.ffn                 Nucleotide CDS sequences
│       └── *.tsv                 Tab-separated feature table
├── taxonomy/
│   └── kaiju/                    Contig-level taxonomy (if --kaiju_db)
│       ├── kaiju_summary.tsv     Per-contig taxonomic assignments
│       └── kaiju_names.tsv       Assignments with full taxon names
├── mge/
│   ├── genomad/                  Virus + plasmid detection (if --genomad_db)
│   │   ├── virus_summary.tsv     Virus contigs with scores + taxonomy
│   │   ├── plasmid_summary.tsv   Plasmid contigs with scores
│   │   ├── virus.fna             Viral contig sequences
│   │   ├── plasmid.fna           Plasmid contig sequences
│   │   ├── provirus.tsv          Provirus boundaries + integrase calls
│   │   └── taxonomy.tsv          Per-contig taxonomy assignments
│   ├── checkv/                   Viral QA (if --checkv_db)
│   │   ├── quality_summary.tsv   Completeness + contamination
│   │   ├── viruses.fna           Host-trimmed viral sequences
│   │   └── proviruses.fna        Extracted provirus sequences
│   ├── integrons/                Integron detection (if --run_integronfinder)
│   │   ├── integrons.tsv         Per-element annotations
│   │   └── summary.tsv           Counts per contig
│   ├── genomic_islands/          Genomic island detection (if --run_islandpath)
│   │   └── genomic_islands.tsv   Island coordinates
│   ├── macsyfinder/              Secretion + conjugation (if --macsyfinder_models)
│   │   ├── all_systems.tsv       Detected systems with component hits
│   │   └── all_systems.txt       Human-readable system descriptions
│   └── defensefinder/            Anti-phage defense (if --run_defensefinder)
│       ├── systems.tsv           Detected defense systems (CRISPR, R-M, BREX, etc.)
│       ├── genes.tsv             Per-gene assignments within systems
│       └── hmmer.tsv             Raw HMM hits across all models
├── metabolism/                    Metabolic profiling (if --run_metabolism)
│   ├── kofamscan/
│   │   └── kofamscan_results.tsv  Per-protein KO assignments (adaptive threshold)
│   ├── emapper/
│   │   └── emapper_results.emapper.annotations  COG/GO/EC/KEGG/Pfam
│   ├── dbcan/
│   │   └── overview.txt           CAZyme consensus (>=2/3 methods agree)
│   ├── merged/
│   │   └── merged_annotations.tsv Unified per-protein annotation table
│   ├── per_mag/
│   │   └── *.tsv                  Per-MAG annotation tables
│   ├── modules/
│   │   ├── module_completeness.tsv  MAG x module completeness matrix
│   │   └── module_heatmap.svg       Clustered heatmap visualization
│   └── community/
│       └── community_annotations.tsv  All proteins with bin_id column
└── pipeline_info/
    ├── run_command.sh            Exact re-runnable command (for -resume)
    ├── timeline.html
    ├── report.html
    └── trace.txt
```

## Parameters

### Assembly

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--input` | (required) | Directory containing `*.fastq.gz` files |
| `--outdir` | `results` | Output directory |
| `--min_overlap` | `1000` | Flye `--min-overlap` |
| `--polish` | `true` | Flye polishing iterations |
| `--dedupe` | `false` | BBDuk deduplication before assembly |
| `--filtlong_size` | (skip) | Filtlong target bases (e.g. `40000000000`) |

### Binning

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--run_maxbin` | `true` | Include MaxBin2 in consensus |
| `--run_lorbin` | `true` | Include LorBin in consensus (deep learning, long-read) |
| `--run_comebin` | `true` | Include COMEBin in consensus (contrastive learning) |
| `--lorbin_min_length` | `80000` | LorBin `--bin_length` minimum (bp) |
| `--metabat_min_cls` | `50000` | MetaBAT2 minimum cluster size |

### Annotation & Taxonomy

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--run_prokka` | `true` | Run Prokka gene annotation on co-assembly |
| `--run_kaiju` | `true` | Run Kaiju protein-level taxonomy (requires Prokka + `kaiju_db`) |
| `--kaiju_db` | (skip) | Path to Kaiju database (`*.fmi` + `nodes.dmp` + `names.dmp`) |

### Mobile Genetic Elements

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--run_genomad` | `true` | Run geNomad virus + plasmid detection |
| `--genomad_db` | (skip) | Path to geNomad database |
| `--run_checkv` | `true` | Run CheckV viral quality assessment (requires geNomad) |
| `--checkv_db` | (skip) | Path to CheckV database |
| `--run_integronfinder` | `true` | Run IntegronFinder integron detection |
| `--run_islandpath` | `true` | Run IslandPath-DIMOB genomic island detection (requires Prokka) |
| `--run_macsyfinder` | `true` | Run MacSyFinder secretion/conjugation detection (requires Prokka) |
| `--macsyfinder_models` | (skip) | Path to MacSyFinder models (TXSScan + CONJScan) |
| `--run_defensefinder` | `true` | Run DefenseFinder anti-phage defense detection (requires Prokka) |
| `--defensefinder_models` | (skip) | Path to DefenseFinder models; null = auto-download |

### Metabolic Profiling

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--run_metabolism` | `false` | Run metabolic profiling (KofamScan + eggNOG + dbCAN) |
| `--kofam_db` | (skip) | Path to KOfam profiles dir (contains `profiles/` + `ko_list`) |
| `--eggnog_db` | (skip) | Path to eggNOG-mapper database dir |
| `--dbcan_db` | (skip) | Path to dbCAN database dir |

### Quality & Resources

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--checkm2_db` | (skip) | Path to CheckM2 DIAMOND database |
| `--assembly_cpus` | `24` | CPUs for assembly |
| `--assembly_memory` | `64 GB` | Memory for assembly |

### Profiles

| Profile | Use case |
|---------|----------|
| `standard` | Local execution (default) |
| `test` | Small test data, reduced resources |
| `shipboard` | Production: 32 CPUs, 256 GB RAM |

## Databases

Large reference databases are not included with the pipeline. Download them before first use:

```bash
cd nextflow

# Interactive menu
./download-databases.sh

# Or download specific databases
./download-databases.sh --genomad          # ~3.5 GB (virus + plasmid detection)
./download-databases.sh --checkv           # ~1.4 GB (viral quality assessment)
./download-databases.sh --checkm2          # ~3.5 GB (MAG quality assessment)
./download-databases.sh --kaiju            # ~47 GB  (protein-level taxonomy)
./download-databases.sh --macsyfinder      # ~50 MB  (secretion + conjugation models)
./download-databases.sh --defensefinder    # ~100 MB (anti-phage defense models)
./download-databases.sh --kofam            # ~4 GB   (KEGG Orthology HMM profiles)
./download-databases.sh --eggnog           # ~12 GB  (eggNOG-mapper DIAMOND db)
./download-databases.sh --dbcan            # ~2 GB   (dbCAN HMM + DIAMOND db)
./download-databases.sh --all              # All databases
```

## Design Notes

**CoverM for depth calculation.** Replaces `jgi_summarize_bam_contig_depths`, which has an integer overflow bug in MetaBAT2 <=2.17 when processing long-read BAMs. CoverM handles supplementary alignments correctly.

**Supplementary alignment filtering.** The mapping step uses `samtools view -F 0x904` to drop unmapped, secondary, and supplementary alignments. Chimeric long reads produce supplementary records that cause massive depth overcounting.

**Dynamic binner architecture.** Each binner emits `[label, file]` tuples that are mixed and collected for DAS Tool. Adding a new binner requires only a process definition and one line in `main.nf`.

**GPU-accelerated ML binners.** The local conda envs include `pytorch-gpu` for GPU-accelerated SemiBin2, LorBin, and COMEBin. The Docker image uses CPU-only PyTorch to keep the image small.

**Graceful failure handling.** All processes handle edge cases (empty input, tool crashes, 0 bins) without crashing the pipeline. DefenseFinder and MacSyFinder produce empty TSVs with headers on failure.

**Sixteen isolated conda environments.** Each tool or group of compatible tools gets its own environment to avoid dependency conflicts. See `install.sh --check` for status.

## Conda Environments

| Environment | Tools |
|-------------|-------|
| `dana-mag-flye` | Flye, Filtlong, Nextflow, OpenJDK |
| `dana-mag-mapping` | minimap2, samtools, CoverM |
| `dana-mag-semibin` | SemiBin2, LorBin, PyTorch GPU |
| `dana-mag-comebin` | COMEBin (rec3141 fork), PyTorch GPU |
| `dana-mag-binning` | MetaBAT2, MaxBin2, DAS_Tool |
| `dana-mag-genomad` | geNomad |
| `dana-mag-checkv` | CheckV |
| `dana-mag-integron` | IntegronFinder |
| `dana-mag-islandpath` | Python + HMMER (genomic island detection) |
| `dana-mag-macsyfinder` | MacSyFinder v2 |
| `dana-mag-defensefinder` | DefenseFinder |
| `dana-mag-kofamscan` | KofamScan, HMMER (KEGG Orthology) |
| `dana-mag-emapper` | eggNOG-mapper, DIAMOND (COG/GO/EC/Pfam) |
| `dana-mag-dbcan` | run_dbcan, HMMER, DIAMOND (CAZyme) |
| `dana-mag-checkm2` | CheckM2 |
| `dana-mag-kaiju` | Kaiju |
| `dana-bbmap` | BBMap (optional dedupe) |

## MAG Quality Standards (MIMAG)

| Tier | Completeness | Contamination | Additional |
|------|-------------|---------------|------------|
| High quality | >90% | <5% | 23S, 16S, 5S rRNA + tRNAs |
| Medium quality | >50% | <10% | -- |
| Low quality | <50% | <10% | -- |

## References

- Flye: Kolmogorov et al., *Nature Biotechnology* 2019
- SemiBin2: Pan et al., *Nature Communications* 2023
- MetaBAT2: Kang et al., *PeerJ* 2019
- LorBin: Gao et al., *Briefings in Bioinformatics* 2024
- COMEBin: Xie et al., *Nature Communications* 2024
- DAS Tool: Sieber et al., *Nature Microbiology* 2018
- CheckM2: Chklovski et al., *Nature Methods* 2023
- geNomad: Camargo et al., *Nature Biotechnology* 2024
- CheckV: Nayfach et al., *Nature Biotechnology* 2021
- IntegronFinder: Cury et al., *Nucleic Acids Research* 2016
- IslandPath-DIMOB: Bertelli & Brinkman, *Bioinformatics* 2018
- MacSyFinder: Abby et al., *PLoS ONE* 2014
- DefenseFinder: Tesson et al., *Nucleic Acids Research* 2022
- Kaiju: Menzel et al., *Nature Communications* 2016
- Prokka: Seemann, *Bioinformatics* 2014
- KofamScan: Aramaki et al., *Bioinformatics* 2020
- eggNOG-mapper: Cantalapiedra et al., *Molecular Biology and Evolution* 2021
- dbCAN3: Zheng et al., *Nucleic Acids Research* 2023
- CoverM: [github.com/wwood/CoverM](https://github.com/wwood/CoverM)
- MIMAG: Bowers et al., *Nature Biotechnology* 2017
