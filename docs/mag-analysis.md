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
| 4c | `BAKTA_EXTRA` | Bakta | Full-database annotation (CDS, sORFs, expert systems; rRNA/tRNA/ncRNA and CRISPR are skipped) |
| 4d | `BAKTA_EXTRA_IDS` | Python | Renames `BAKTA_EXTRA` features to the `BAKTA_BASIC` locus tag of the CDS at the same coordinates, since all protein tools run on the basic proteins; sORFs get new tags after the last basic one. Writes `locus_tag_map.tsv` |

### Taxonomy

| Step | Process | Tool | Description |
|------|---------|------|-------------|
| 5a | `KAIJU_CONTIG_CLASSIFY` | Kaiju | Protein-level contig classification |
| 5b | `KAIJU_CLASSIFY` | Kaiju | Protein-level MAG classification |
| 5c | `KRAKEN2_CLASSIFY` | Kraken2 | k-mer-based classification |
| 5d | `SENDSKETCH_CLASSIFY` | BBTools sendsketch | MinHash taxonomy (GTDB) |
| 5e | `RNA_CLASSIFY` | barrnap + VSEARCH | rRNA extraction + SILVA classification |
| 5f | `GTDBTK_CLASSIFY` | GTDB-Tk | Genome-based taxonomy (160 GB RAM declared) |

### Metabolism

| Step | Process | Tool | Description |
|------|---------|------|-------------|
| 6a | `KOFAMSCAN` | KofamScan | KEGG ortholog assignment via HMM profiles |
| 6b | `EMAPPER` | eggNOG-mapper | Ortholog annotation (COG, GO, KEGG, CAZy) |
| 6c | `DBCAN` | dbCAN | CAZyme annotation via HMM + diamond + HMMER |
| 6d | `PROTEIN_DB_HMM`, `PROTEIN_DB_DIAMOND` | hmmsearch, DIAMOND | Protein reference DBs under `--protein_db_dir` (CANT-HYD, HADEG, PGAP/TIGRFAM, ...), one task per DB |
| 6e | `MERGE_ANNOTATIONS` | Python | Merge KofamScan + eggNOG + dbCAN + protein-DB hits per gene |
| 6f | `MAP_TO_BINS` | Python | Assign annotations to MAGs |
| 6g | `KEGG_MODULES` | Python | KEGG module completeness per MAG (55 bundled modules) |
| 6h | `MINPATH` | MinPath | Pathway parsimony analysis |
| 6i | `KEGG_DECODER` | KEGG-Decoder | Metabolic pathway heatmaps |
| 6j | `ANTISMASH` | antiSMASH | Biosynthetic gene cluster detection |
| 6k | `ECOSSDB_MAP` | Python (ecossdb) | Ecosystem services mapping (CICES 5.2) |
| 6l | `ECOSSDB_SCORE` | Python (ecossdb) | Ecosystem service scoring |
| 6m | `ECOSSDB_SDG` | Python (ecossdb) | UN Sustainable Development Goal mapping |
| 6n | `ECOSSDB_VIZ` | Python (ecossdb) | Ecosystem services visualization |
| 6o | `METACYC_SEARCH` | DIAMOND | Proteins against Rhea-annotated UniProtKB proteins (`--metacyc_refs`) |
| 6p | `METACYC_PATHWAYS` | Python | MetaCyc reactions and pathway calls per DAS Tool MAG and community, with confidence (`--metacyc_db`) |

#### Protein reference databases

`--protein_db_dir DIR` (with `--run_metabolism`; `--db_dir` resolves
`<db_dir>/protein_dbs`) searches the annotation proteins against every database
in `DIR`, one subdirectory per database. The layout is the one nanopore_live's
`--pollutant_dbs` reads, so one directory serves both pipelines:

| File | Meaning |
|------|---------|
| `<name>/<name>.hmm` | HMM profiles, searched with `hmmsearch` |
| `<name>/<name>.dmnd` | DIAMOND database, searched with `diamond blastp` (best hit per protein) |
| `<name>/<name>.genes.tsv` | Optional catalogue: `target gene category pathway subpathway description [threshold score_type motif]` |
| `<name>/<name>.opts` | Optional cutoff: `tc`, `ga`, `nc` or an E-value for HMMs, optionally followed by `resolve`; `id=,qcov=,scov=,evalue=,len=,bits=` for DIAMOND |

Without an `.opts` file, HMMs use `--hmm_cut` (default `tc`), or the
catalogue's per-target thresholds when `<name>.genes.tsv` has a `threshold`
column; DIAMOND uses `--protein_db_opts`. Each database's hits are in
`metabolism/protein_dbs/<name>/<name>.protein_db.tsv`, and
`merged_annotations.tsv` gains a `protein_db_hits` column (`<db>:<gene>`,
comma-separated). The results are stored per database, so adding one later
searches only that one; after changing a database or its `.opts`, delete its
`metabolism/protein_dbs/<name>/` to search it again.

An HMM search over more than `--protein_db_chunk` proteins (default 50,000) is
split into chunks searched in parallel, 4 `hmmsearch` threads per chunk, since
one `hmmsearch` gains little from more threads. E-values are computed against
the whole protein set, so the hits are those of a single search. Give the step
the node with `--protein_db_cpus`: on 72,231 proteins PGAP took ~2.2 CPU-hours
and the full CATH FunFam library ~8.7, and both scale with total protein length.

With `resolve` (e.g. `1e-3 resolve`), the per-domain hits go through
`cath-resolve-hits`, and only each protein's best non-overlapping set of models
is kept. This is for libraries of overlapping families, such as CATH FunFam,
where one protein otherwise matches dozens of related models.

`./download-databases.sh --funfam` writes `<db_dir>/protein_dbs/CATH-FunFam/`:
CATH FunFam v4.3.0 (212,872 models), searched at `1e-3 resolve`;
`--funfam-superfamilies 3.40.50.300,1.10.8` keeps only some superfamilies.

`./download-databases.sh --pgap` writes `<db_dir>/protein_dbs/PGAP/`: NCBI's
PGAP HMMs, which hold the maintained TIGRFAMs (accessions `TIGR*`) plus
NCBIfam and PRK models, searched at their noise cutoffs (`nc`; the trusted
cutoffs NCBI uses for genomes are often too strict for metagenomes).

#### MetaCyc pathways

With `--metacyc_db` (and `--metacyc_refs`), mag_analysis predicts which MetaCyc pathways
each DAS Tool MAG, and the community as a whole, carries. It does not run Pathway Tools on
your data. Proteins get MetaCyc reactions through Rhea: DIAMOND hits to Rhea-annotated
UniProtKB proteins, plus Bakta/eggNOG EC numbers and KofamScan/eggNOG KOs. Pathways are then
called with PathoLogic's published rules (Karp et al. 2011,
[doi:10.4056/sigs.1794338](https://doi.org/10.4056/sigs.1794338)). Those rules are
permissive, so every call carries a confidence value built from its evidence:

- each reaction is weighted by its strongest evidence: a DIAMOND hit scores 0.5 at 35%
  identity rising to 1.0 at 80%; EC-only and KO-only evidence score 0.45; an EC mapped
  straight to a MetaCyc reaction without a Rhea link scores 0.30
- `confidence` is the mean weight over the pathway's reactions (0-1)
- tier **high** means confidence >= 0.6, at least half of the present reactions backed by a
  sequence hit, and every key reaction sequence-backed; **medium** means confidence >= 0.35;
  everything else is **low**

The weights are first estimates. Treat low-tier calls as possibilities rather than findings.
The step needs Bakta annotation (`--annotator bakta`). It is skipped, with one log line, when
`--metacyc_db` is not given.

**Getting MetaCyc.** MetaCyc is licensed by SRI International and is never downloaded
automatically. You need one of:

| You have | Where it comes from |
|---|---|
| MetaCyc flat files (a tarball or directory containing `pathways.dat`, `reactions.dat`, `classes.dat`) | a BioCyc subscription (flat-file download for MetaCyc) |
| the Pathway Tools Linux installer (`pathway-tools-XX.X-linux-64-tier1-install`) | a Pathway Tools licence from SRI (free for academic use) |
| an existing Pathway Tools Apptainer image (`.sif`), e.g. from MetaPathways `build_pt` | your site |

**Setup (once per site):**

```bash
./download-databases.sh --db_dir /path/to/databases --metacyc-refs          # open references
./download-databases.sh --db_dir /path/to/databases --metacyc <tarball|dir|installer|.sif>
# installer on an Apptainer-only host (e.g. an HPC cluster): give the MetaPathways image as the runtime
apptainer pull mp.sif docker://quay.io/hallamlab/metapathways:4.0.3
./download-databases.sh --db_dir /path/to/databases --metacyc pathway-tools-30.0-linux-64-tier1-install \
    --metacyc-runtime mp.sif
```

| Step | Writes | Time | Disk |
|---|---|---|---|
| `--metacyc-refs` | `<db_dir>/metacyc_refs/`: Rhea, KEGG ko->reaction, GTDB r226->NCBI taxids, and the DIAMOND reference (Swiss-Prot plus one Rhea-annotated protein per UniRef90 cluster, 9.4M sequences) | ~1 h, mostly downloads (UniProt idmapping 7 GB; TrEMBL ~40 GB is streamed through a filter, never stored) | ~15 GB scratch at peak, 5.5 GB kept |
| `--metacyc` (flat files) | `<db_dir>/metacyc/` | seconds | ~140 MB |
| `--metacyc` (installer or .sif) | `<db_dir>/metacyc/`, exported in batch mode from a temporary Pathway Tools install | ~3 min | 3.5 GB temporary |

The installer runs as your own user, with no root, inside the Pathway Tools runtime that ships
in the MetaPathways image, and the install is deleted afterwards. `<db_dir>/metacyc/` is
created owner-only and carries a `VERSION.json` (MetaCyc release, source, date).

**Licence.** Keep `<db_dir>/metacyc/` private. Don't put it in container images, git
repositories or shared or public locations. The pipeline outputs carry MetaCyc pathway and
reaction identifiers and names, so share them under the same terms as MetaCyc itself.

`--db_dir` finds `<db_dir>/metacyc` and `<db_dir>/metacyc_refs` automatically. Taxonomic-range
rules use GTDB-Tk's classification (`--run_gtdbtk`). Without GTDB-Tk,
`--metacyc_taxonomy sendsketch` uses SendSketch contig calls instead, which are less reliable;
`none` turns the range rules off.

Outputs in `metabolism/metacyc/`:

- `metacyc_pathways.tsv`: one row per bin (plus `_community`) and pathway with any
  evidence: the call and the rule that decided it, `confidence` and `confidence_tier`,
  reaction counts by evidence (`n_seq`, `n_ec_only`, `n_ko_only`, `n_ec_direct_only`),
  sequence identity, key and unique reaction support, up to 12 supporting hits
  (`reaction=protein:UniProt:identity`), and a MinPath-style parsimony flag (`minpath`)
- `metacyc_reactions.tsv.gz`: one row per protein and MetaCyc reaction, with its evidence
  classes, DIAMOND identity, coverage and hit, and the ECs and KOs behind it
- `metacyc_bin_taxonomy.tsv`: each bin's GTDB lineage and the NCBI taxids used for the
  taxonomic ranges
- `metacyc_hits.tsv.gz`: the DIAMOND hits

On 72,231 proteins (7 MAGs), the DIAMOND search took 21 min on 8 threads with a 17 GB peak
(`sensitive`, `-b 2`), and the scorer took 14 s and 0.6 GB. For a co-assembly of ~13M
proteins, expect the search to take ~1.5-2.5 h on a 192-CPU node (`--metacyc_cpus 192`;
raise `--metacyc_block_size` to 4-6 and `--metacyc_memory` to ~96 GB for fewer passes over
the reference). The scorer then takes ~20 min, ~5 GB, and as much scratch in its work dir as
the input tables. Results do not depend on `--metacyc_chunk`; a smaller value only lowers memory.

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
| 9 | `CALCULATE_GENE_DEPTHS` | samtools + Python | Per-gene coverage from BAMs |
| 10 | `VIZ_PREPROCESS` | Python + Svelte | Interactive dashboard generation |

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

Note: the bin-size cut-offs `--binette_min_bin_bp`, `--binette_max_bin_bp`,
`--checkm2_max_bin_bp` and `--gtdbtk_max_bin_bp` are compared against the bin
FASTA's size in bytes, which approximates its length in bp (headers and line
breaks included). `--derep_min_size_bp` counts sequence bases.

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
| `--bakta_light_db` | (required if bakta) | Path to the Bakta light database (`BAKTA_BASIC`) |
| `--bakta_db` | (required if `--bakta_extra`) | Path to the Bakta full database (`BAKTA_EXTRA`) |
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
| `--rrna_min_identity` | `0.80` | Minimum identity for the rRNA VSEARCH search |
| `--run_gtdbtk` | `false` | Run GTDB-Tk classification |
| `--gtdbtk_db` | (required if gtdbtk) | GTDB-Tk database path |

### Metabolism

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--run_metabolism` | `false` | Enable metabolic profiling |
| `--kofam_db` | (required if metabolism) | KofamScan database path |
| `--eggnog_db` | (required if metabolism) | eggNOG database path |
| `--dbcan_db` | (required if metabolism) | dbCAN database path |
| `--protein_db_dir` | `null` | Directory of protein reference DBs (see Protein reference databases) |
| `--hmm_cut` | `tc` | Default hmmsearch cutoff: `tc`, `ga`, `nc` or an E-value |
| `--protein_db_opts` | `id=50,qcov=50,scov=50,evalue=1e-5` | Default DIAMOND cutoffs |
| `--protein_db_chunk` | `50000` | Proteins per parallel hmmsearch job (`0`: one search per database) |
| `--protein_db_cpus` | `--medium_cpus` | CPUs for each protein DB search |
| `--emapper_batch_size` | `50000` | eggNOG-mapper batch size (proteins) |
| `--run_antismash` | `false` | Run antiSMASH BGC detection |
| `--antismash_db` | (required if antismash) | antiSMASH database path |
| `--run_ecossdb` | `true` | Enable ecosystem services mapping |
| `--metacyc_db` | `null` | Licensed MetaCyc flat files (see MetaCyc pathways); step skipped without it |
| `--metacyc_refs` | `null` | Open references for MetaCyc pathways (`download-databases.sh --metacyc-refs`) |
| `--metacyc_taxonomy` | `gtdbtk` | Taxonomic ranges from `gtdbtk`, `sendsketch` (contig calls) or `none` |
| `--metacyc_sensitivity` | `sensitive` | DIAMOND mode (`fast`, `sensitive`, `more-sensitive`) |
| `--metacyc_block_size` | `2.0` | DIAMOND `-b`; memory grows with it |
| `--metacyc_chunk` | `500000` | Proteins per scorer chunk (~1.7 GB of `METACYC_PATHWAYS` memory each) |
| `--metacyc_cpus` / `--metacyc_memory` | `--medium_cpus` / `32 GB` | `METACYC_SEARCH` resources |

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
| `--all` | `false` | `run-mag-analysis.sh` only: expands to every `--run_*` flag (not a Nextflow param) |
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
│   ├── protein_dbs/<name>/     Protein reference DB hits (--protein_db_dir)
│   ├── merged/                 Merged annotation table
│   ├── per_mag/                Per-MAG annotation summaries
│   ├── modules/                KEGG module completeness
│   ├── minpath/                MinPath pathway analysis
│   ├── kegg_decoder/           Pathway heatmaps
│   ├── metacyc/                MetaCyc pathway calls + confidence (--metacyc_db)
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

Defaults from `nextflow.config`; the `--medium_cpus/--medium_memory`,
`--binning_cpus/--binning_memory`, `--bakta_cpus/--bakta_memory` and
`--assembly_cpus/--assembly_memory` parameters resize whole groups of processes.

| Component | CPUs | RAM | Notes |
|-----------|------|-----|-------|
| MetaBAT2 | 8 | 16 GB | `process_medium` label (`--medium_*`) |
| MaxBin2 | 8 | 16 GB | Sized with the binners (`--binning_*`) |
| SemiBin2 / LorBin / VAMB | 8 | 16 GB | `process_gpu` label (`--binning_*`); GPU optional |
| COMEBin | 8 | 16 GB | `process_gpu`; skipped on CPU unless `--comebin_on_cpu` |
| Prokka | 8 | 16 GB | `process_medium` |
| Bakta (basic and extra) | 8 | 24 GB | `--bakta_*` |
| eggNOG-mapper | 8 | 48 GB | Batched protein annotation |
| MetaEuk | 8 | 32 GB | Eukaryotic gene prediction |
| GTDB-Tk | 16 | 160 GB | Most memory-intensive step (R232 placement needs >=140 GB) |
| Kraken2 / Kaiju | 8 | 24 GB | `process_kraken` label; depends on database |
| CheckM2 | 8 | 16 GB | `process_medium` |
| MetaCyc search (DIAMOND) | 8 | 32 GB | `--metacyc_cpus` / `--metacyc_memory` |
| MetaCyc pathways | 2 | 8 GB | Single-threaded; memory set by `--metacyc_chunk` |

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
./download-databases.sh --pgap      # NCBI PGAP HMMs incl. TIGRFAM (~2.7 GB)
./download-databases.sh --metacyc-refs        # MetaCyc open references (~5.5 GB)
./download-databases.sh --metacyc <source>    # licensed MetaCyc you supply (see MetaCyc pathways)
```

Use `--db_dir` to auto-resolve paths from a standard layout:

```bash
./run-mag-analysis.sh --db_dir /path/to/databases ...
```

## Design Notes

- **Seven-binner consensus.** SemiBin2, MetaBAT2, MaxBin2, LorBin, COMEBin, VAMB, VAMB-tax feed into DAS Tool + Binette + MAGScoT for consensus refinement.
- **Modular activation.** MetaBAT2, DAS Tool, Bakta annotation (`--annotator bakta`) and ECOSSDB (`--run_ecossdb`, when metabolism runs) are on by default; every other module is off until enabled with its `--run_*` flag, or with the launcher's `--all`.
- **Technology-agnostic.** Accepts assemblies from Nanopore, Illumina, HiFi, or any external assembler. Only requires FASTA + MetaBAT2-format depth table.
- **Persistent caching.** Use `--store_dir` to skip completed processes across runs, even after `work/` cleanup.
- **Quality tiers.** The viz preprocessing tiers MAGs from CheckM2 completeness/contamination as HQ (>=90%, <5%), MQ (>=50%, <10%) or LQ; the MIMAG rRNA/tRNA criteria are not applied.
