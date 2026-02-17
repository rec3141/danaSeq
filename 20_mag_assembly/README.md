# MAG Assembly Pipeline

Metagenome-assembled genome (MAG) reconstruction from Oxford Nanopore long reads. Co-assembles reads with Flye, maps back with minimap2, runs five binning algorithms (SemiBin2, MetaBAT2, MaxBin2, LorBin, COMEBin), and integrates results with DAS Tool consensus. Includes gene annotation (Prokka/Bakta), four taxonomy classifiers (Kaiju, Kraken2, sendsketch, rRNA/SILVA), mobile genetic element detection (geNomad, CheckV, IntegronFinder, IslandPath, MacSyFinder, DefenseFinder), metabolic profiling with pathway analysis (KofamScan, eggNOG-mapper, dbCAN3, MinPath, KEGG-Decoder), eukaryotic analysis (Tiara, Whokaryote, MetaEuk), and optional LLM-guided bin refinement (NCLB).

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
    --annotator prokka \
    --run_maxbin true \
    --run_lorbin true \
    --run_comebin true \
    --lorbin_min_length 80000 \
    --metabat_min_cls 50000 \
    --checkm2_db /path/to/checkm2_db \
    --genomad_db /path/to/genomad_db \
    --checkv_db /path/to/checkv_db \
    --kaiju_db /path/to/kaiju_db \
    --run_kraken2 true --kraken2_db /path/to/kraken2_db \
    --run_sendsketch true --sendsketch_address http://host:3068/sketch \
    --run_rrna true --silva_ssu_db /path/to/SILVA_SSU.fasta \
    --macsyfinder_models /path/to/macsyfinder_models \
    --defensefinder_models /path/to/defensefinder_models \
    --run_metabolism true \
    --kofam_db /path/to/kofam_db \
    --eggnog_db /path/to/eggnog_db \
    --dbcan_db /path/to/dbcan_db \
    --run_eukaryotic true \
    --run_metaeuk true --metaeuk_db /path/to/metaeuk_db \
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
      +-- KEGG_MODULES      Module completeness scoring + heatmap
      +-- MINPATH            Parsimony pathway reconstruction (Ye & Doak 2009)
      +-- KEGG_DECODER       Biogeochemical function scoring + heatmap (Graham et al. 2018)

   KRAKEN2_CLASSIFY          k-mer contig-level taxonomy (no annotation needed)
   SENDSKETCH_CLASSIFY       GTDB MinHash taxonomy (BBSketch TaxServer)
   RRNA_CLASSIFY             rRNA gene detection (barrnap) + SILVA classification (vsearch)

   TIARA_CLASSIFY            Deep learning eukaryotic classification
   WHOKARYOTE_CLASSIFY       Gene structure-based eukaryotic classification
      +-- METAEUK_PREDICT    Multi-exon eukaryotic gene prediction

   GENOMAD_CLASSIFY          Virus + plasmid + provirus detection
      |
      +-- CHECKV_QUALITY     Viral genome quality assessment

   INTEGRONFINDER           Integron + gene cassette detection
   CALCULATE_TNF            Tetranucleotide frequency profiles

  --- optional bin refinement (after DAS Tool) ---

   NCLB_GATHER              Build contig identity cards
   NCLB_CONVERSE            LLM-guided bin placement conversations
   NCLB_ELDERS              SCG redundancy investigation (ecotype vs contamination)
   NCLB_INTEGRATE           Apply proposals, extract refined community FASTAs
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
│   └── prokka/ or bakta/         Gene annotation (depending on --annotator)
│       ├── *.gff                 GFF3 annotations
│       ├── *.faa                 Protein sequences
│       ├── *.ffn                 Nucleotide CDS sequences
│       └── *.tsv                 Tab-separated feature table
├── taxonomy/
│   ├── kaiju/                    Protein-level taxonomy (if --kaiju_db)
│   │   ├── kaiju_genes.tsv       Per-gene Kaiju classifications
│   │   └── kaiju_contigs.tsv     Per-contig taxonomy (majority vote)
│   ├── kraken2/                  k-mer taxonomy (if --kraken2_db)
│   │   ├── kraken2_contigs.tsv   Per-contig classifications + lineage
│   │   └── kraken2_report.txt    Standard Kraken2 report (for Krona/Pavian)
│   ├── sendsketch/               GTDB MinHash taxonomy (if --sendsketch_address)
│   │   └── sendsketch_contigs.tsv  Per-contig GTDB taxonomy + ANI
│   └── rrna/                     rRNA classification (if --silva_ssu_db)
│       ├── rrna_genes.tsv        Per-gene rRNA classifications (barrnap + vsearch)
│       ├── rrna_contigs.tsv      Per-contig rRNA summary (best SSU/LSU taxonomy)
│       └── rrna_sequences.fasta  Extracted rRNA gene sequences
├── eukaryotic/                   Eukaryotic analysis (if --run_eukaryotic)
│   ├── tiara/
│   │   └── tiara_output.tsv      Per-contig classification + probabilities
│   ├── whokaryote/
│   │   └── whokaryote_classifications.tsv  Per-contig classification
│   └── metaeuk/                  Eukaryotic gene prediction (if --metaeuk_db)
│       ├── metaeuk_proteins.fas  Multi-exon protein predictions
│       ├── metaeuk_codon.fas     Nucleotide coding sequences
│       ├── metaeuk.gff           Gene structures (exon boundaries)
│       └── metaeuk_headers.tsv   Internal ID mapping
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
│   │   └── overview.tsv           CAZyme consensus (>=2/3 methods agree)
│   ├── merged/
│   │   └── merged_annotations.tsv Unified per-protein annotation table
│   ├── per_mag/
│   │   └── *.tsv                  Per-MAG annotation tables
│   ├── modules/
│   │   ├── module_completeness.tsv  MAG x module completeness matrix
│   │   └── module_heatmap.svg       Clustered heatmap visualization
│   ├── minpath/
│   │   ├── minpath_pathways.tsv     MAG x pathway (naive vs parsimony counts)
│   │   └── details/                 Per-MAG MinPath reports
│   ├── kegg_decoder/
│   │   ├── kegg_decoder_output.tsv  MAG x function completeness (~80 functions)
│   │   └── function_heatmap.svg     Publication-quality biogeochemical heatmap
│   └── community/
│       └── community_annotations.tsv  All proteins with bin_id column
├── binning/nclb/                 NCLB bin refinement (if --run_nclb + --nclb_dir)
│   ├── communities/*.fa          Refined community FASTAs
│   ├── gathering.json            Identity cards + resonance data
│   ├── proposals.json            LLM conversation proposals
│   ├── elder_reports.json        SCG redundancy investigations
│   ├── chronicle.json            Machine-readable decision log
│   ├── chronicle.md              Human-readable narrative
│   ├── contig2community.tsv      Contig membership assignments
│   ├── quality_report.tsv        Community quality metrics
│   └── valence_report.tsv        Per-contig valence scores
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

### Annotation

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--annotator` | (auto) | Annotator: `prokka`, `bakta`, or `none`. Overrides legacy flags |
| `--run_prokka` | `true` | (deprecated) Run Prokka — use `--annotator` instead |
| `--run_bakta` | `false` | (deprecated) Run Bakta — use `--annotator` instead |
| `--bakta_db` | (skip) | Path to Bakta database (required when using Bakta) |
| `--bakta_full` | `false` | Also run full Bakta annotation (ncRNA/tRNA/CRISPR — slow) |

### Taxonomy

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--run_kaiju` | `true` | Run Kaiju protein-level taxonomy (requires annotation + `kaiju_db`) |
| `--kaiju_db` | (skip) | Path to Kaiju database (`*.fmi` + `nodes.dmp` + `names.dmp`) |
| `--run_kraken2` | `false` | Run Kraken2 k-mer taxonomy on contigs (no annotation needed) |
| `--kraken2_db` | (skip) | Path to Kraken2 database (`hash.k2d` + `nodes.dmp` + `names.dmp`) |
| `--kraken2_confidence` | `0.0` | Kraken2 confidence threshold (0.0 = any match) |
| `--run_sendsketch` | `false` | Run BBSketch/sendsketch GTDB taxonomy (requires TaxServer) |
| `--sendsketch_address` | (skip) | BBSketch TaxServer URL (e.g. `http://host:3068/sketch`) |
| `--run_rrna` | `false` | Run barrnap + vsearch rRNA gene classification |
| `--silva_ssu_db` | (skip) | Path to SILVA SSU NR99 FASTA (DNA, U→T converted) |
| `--silva_lsu_db` | (skip) | Path to SILVA LSU NR99 FASTA (optional) |
| `--rrna_min_identity` | `0.80` | Minimum vsearch identity for rRNA classification |

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

### Eukaryotic Analysis

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--run_eukaryotic` | `false` | Enable eukaryotic contig classification (Tiara + Whokaryote) |
| `--tiara_min_len` | `3000` | Minimum contig length for Tiara (bp) |
| `--whokaryote_min_len` | `5000` | Minimum contig length for Whokaryote (bp) |
| `--run_metaeuk` | `false` | Run MetaEuk eukaryotic gene prediction (requires `--run_eukaryotic`) |
| `--metaeuk_db` | (skip) | Path to MetaEuk protein reference database (MMseqs2 format) |
| `--metaeuk_mem_limit` | `50G` | MetaEuk `--split-memory-limit` |
| `--metaeuk_max_intron` | `10000` | Maximum intron length in bp |

### Bin Refinement (NCLB)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--run_nclb` | `false` | Run NCLB LLM-guided bin refinement after DAS Tool |
| `--nclb_dir` | (skip) | Path to NCLB repository (contains `bin/`, `lib/`, `envs/`) |
| `--nclb_base_url` | (auto) | LLM server URL (default: `http://localhost:1234/v1`) |
| `--nclb_model` | (auto) | LLM model name (default: auto-detect from server) |
| `--nclb_with_ani` | `false` | Run minimap2 ANI during Elder investigations |

### Quality & Resources

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--checkm2_db` | (skip) | Path to CheckM2 DIAMOND database |
| `--assembly_cpus` | `16` | CPUs for assembly |
| `--assembly_memory` | `60 GB` | Memory for assembly |

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
./download-databases.sh --genomad          # ~3.5 GB  (virus + plasmid detection)
./download-databases.sh --checkv           # ~1.4 GB  (viral quality assessment)
./download-databases.sh --checkm2          # ~3.5 GB  (MAG quality assessment)
./download-databases.sh --kaiju            # ~47 GB   (protein-level taxonomy)
./download-databases.sh --kraken2          # ~8-70 GB (k-mer taxonomy, varies by DB)
./download-databases.sh --silva            # ~1 GB    (rRNA classification - SSU + LSU)
./download-databases.sh --macsyfinder      # ~50 MB   (secretion + conjugation models)
./download-databases.sh --defensefinder    # ~100 MB  (anti-phage defense models)
./download-databases.sh --kofam            # ~4 GB    (KEGG Orthology HMM profiles)
./download-databases.sh --eggnog           # ~12 GB   (eggNOG-mapper DIAMOND db)
./download-databases.sh --dbcan            # ~2 GB    (dbCAN HMM + DIAMOND db)
./download-databases.sh --bakta            # ~1.5 GB  (light) or ~30 GB (full)
./download-databases.sh --metaeuk          # ~18 GB   (OrthoDB eukaryotic proteins)
./download-databases.sh --all              # All databases
```

## Design Notes

**CoverM for depth calculation.** Replaces `jgi_summarize_bam_contig_depths`, which has an integer overflow bug in MetaBAT2 <=2.17 when processing long-read BAMs. CoverM handles supplementary alignments correctly.

**Supplementary alignment filtering.** The mapping step uses `samtools view -F 0x904` to drop unmapped, secondary, and supplementary alignments. Chimeric long reads produce supplementary records that cause massive depth overcounting.

**Dynamic binner architecture.** Each binner emits `[label, file]` tuples that are mixed and collected for DAS Tool. Adding a new binner requires only a process definition and one line in `main.nf`.

**GPU-accelerated ML binners.** The local conda envs include `pytorch-gpu` for GPU-accelerated SemiBin2, LorBin, and COMEBin. The Docker image uses CPU-only PyTorch to keep the image small.

**Graceful failure handling.** All processes handle edge cases (empty input, tool crashes, 0 bins) without crashing the pipeline. DefenseFinder and MacSyFinder produce empty TSVs with headers on failure.

**Annotate once, map to bins.** Metabolic profiling tools (KofamScan, eggNOG-mapper, dbCAN) run on the full protein FASTA from annotation. Results are merged into a unified per-protein table, then partitioned to individual MAGs via DAS Tool contig-to-bin assignments. This avoids redundant computation and ensures unbinned contigs are also annotated.

**Three pathway analysis approaches.** KEGG module completeness provides step-by-step pathway evaluation. MinPath (Ye & Doak 2009) uses integer programming to find the minimum set of pathways explaining observed KOs, preventing inflation in draft MAGs. KEGG-Decoder (Graham et al. 2018) scores ~80 environmentally-curated biogeochemical functions.

**Four independent taxonomy classifiers.** Kaiju (protein-level, RefSeq), Kraken2 (k-mer, configurable DB), sendsketch (GTDB MinHash via TaxServer), and rRNA (barrnap + vsearch SILVA). Each provides orthogonal signal; users can run any combination.

**Eukaryotic contig filtering.** Tiara (deep learning k-mer NN) and Whokaryote (gene structure random forest) independently classify contigs. MetaEuk runs only on the union of non-prokaryotic contigs, avoiding wasted computation on bacterial sequences.

**NCLB bin refinement.** Optional LLM-guided iterative refinement of DAS Tool bins. Four-step process: gather contig identity cards, LLM conversations for placement proposals, Elder investigations for SCG redundancy, and integration to produce refined community FASTAs.

**Twenty-five isolated conda environments.** Each tool or group of compatible tools gets its own environment to avoid dependency conflicts. See `install.sh --check` for status.

## Conda Environments

Twenty-five isolated environments avoid dependency conflicts (27 YAML specs including CPU variants):

| Environment | Tools |
|-------------|-------|
| `dana-mag-flye` | Flye, Filtlong, Nextflow, OpenJDK |
| `dana-mag-mapping` | minimap2, samtools, CoverM |
| `dana-mag-semibin` | SemiBin2, LorBin, PyTorch GPU |
| `dana-mag-comebin` | COMEBin (rec3141 fork), PyTorch GPU |
| `dana-mag-binning` | MetaBAT2, MaxBin2, DAS_Tool |
| `dana-mag-prokka` | Prokka |
| `dana-mag-bakta` | Bakta (modern Prokka alternative) |
| `dana-mag-kaiju` | Kaiju |
| `dana-mag-kraken2` | Kraken2 (k-mer contig-level taxonomy) |
| `dana-mag-rrna` | barrnap, vsearch (rRNA gene classification) |
| `dana-mag-genomad` | geNomad |
| `dana-mag-checkv` | CheckV |
| `dana-mag-integron` | IntegronFinder |
| `dana-mag-islandpath` | Python + HMMER (genomic island detection) |
| `dana-mag-macsyfinder` | MacSyFinder v2 |
| `dana-mag-defensefinder` | DefenseFinder |
| `dana-mag-kofamscan` | KofamScan, HMMER (KEGG Orthology) |
| `dana-mag-emapper` | eggNOG-mapper, DIAMOND (COG/GO/EC/Pfam) |
| `dana-mag-dbcan` | run_dbcan, HMMER, DIAMOND (CAZyme) |
| `dana-mag-pathway` | MinPath, KEGG-Decoder (pathway analysis) |
| `dana-mag-tiara` | Tiara (deep learning eukaryotic classification) |
| `dana-mag-whokaryote` | Whokaryote, Prodigal (gene structure eukaryotic classification) |
| `dana-mag-metaeuk` | MetaEuk (multi-exon eukaryotic gene prediction) |
| `dana-mag-checkm2` | CheckM2 |
| `dana-bbmap` | BBMap (optional dedupe) |

## MAG Quality Standards (MIMAG)

| Tier | Completeness | Contamination | Additional |
|------|-------------|---------------|------------|
| High quality | >90% | <5% | 23S, 16S, 5S rRNA + tRNAs |
| Medium quality | >50% | <10% | -- |
| Low quality | <50% | <10% | -- |

## References

**Assembly & Binning:**
- Flye: Kolmogorov et al., *Nature Biotechnology* 2019
- SemiBin2: Pan et al., *Nature Communications* 2023
- MetaBAT2: Kang et al., *PeerJ* 2019
- MaxBin2: Wu et al., *Bioinformatics* 2016
- LorBin: Gao et al., *Briefings in Bioinformatics* 2024
- COMEBin: Xie et al., *Nature Communications* 2024
- DAS Tool: Sieber et al., *Nature Microbiology* 2018
- CheckM2: Chklovski et al., *Nature Methods* 2023
- CoverM: [github.com/wwood/CoverM](https://github.com/wwood/CoverM)

**Annotation:**
- Prokka: Seemann, *Bioinformatics* 2014
- Bakta: Schwengers et al., *Microbial Genomics* 2021

**Taxonomy:**
- Kaiju: Menzel et al., *Nature Communications* 2016
- Kraken2: Wood et al., *Genome Biology* 2019
- BBSketch/sendsketch: Bushnell, [sourceforge.net/projects/bbmap](https://sourceforge.net/projects/bbmap/)
- barrnap: Seemann, [github.com/tseemann/barrnap](https://github.com/tseemann/barrnap)
- vsearch: Rognes et al., *PeerJ* 2016
- SILVA: Quast et al., *Nucleic Acids Research* 2013

**Mobile Genetic Elements:**
- geNomad: Camargo et al., *Nature Biotechnology* 2024
- CheckV: Nayfach et al., *Nature Biotechnology* 2021
- IntegronFinder: Cury et al., *Nucleic Acids Research* 2016
- IslandPath-DIMOB: Bertelli & Brinkman, *Bioinformatics* 2018
- MacSyFinder: Abby et al., *PLoS ONE* 2014
- DefenseFinder: Tesson et al., *Nature Communications* 2022

**Metabolic Profiling:**
- KofamScan: Aramaki et al., *Bioinformatics* 2020
- eggNOG-mapper: Cantalapiedra et al., *Molecular Biology and Evolution* 2021
- dbCAN3: Zheng et al., *Nucleic Acids Research* 2023
- MinPath: Ye & Doak, *PLoS Computational Biology* 2009
- KEGG-Decoder: Graham et al., *bioRxiv* 2018

**Eukaryotic Analysis:**
- Tiara: Karlicki et al., *Bioinformatics* 2022
- Whokaryote: Pronk et al., *Microbial Genomics* 2022
- MetaEuk: Levy Karin et al., *Microbiome* 2020

**Standards:**
- MIMAG: Bowers et al., *Nature Biotechnology* 2017
