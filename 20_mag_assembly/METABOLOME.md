# Metabolic Profiling Tools for MAG Analysis

Reference guide for metabolic annotation, functional profiling, and pathway reconstruction tools evaluated for integration with the danaSeq MAG assembly pipeline. Focused on **environmental/aquatic metagenomes** assembled from **Oxford Nanopore long reads** with Flye, binned via DAS Tool consensus.

## Context: Current Pipeline State

The danaSeq pipeline produces quality-filtered MAGs with gene predictions and taxonomic assignments, but **no functional annotation beyond Prokka product descriptions**:

```
CURRENT PIPELINE OUTPUTS (available for metabolic profiling)
├── Assembly
│   ├── assembly.fasta                    Co-assembly (Flye)
│   └── tnf.tsv                           136 tetranucleotide frequencies per contig
│
├── Gene Predictions (Prokka)
│   ├── PROKKA_*.faa                      Predicted protein sequences (all genes)
│   ├── PROKKA_*.gff                      Gene features with coordinates + products
│   └── PROKKA_*.tsv                      Feature summary table
│
├── Taxonomy (Kaiju)
│   ├── kaiju_genes.tsv                   Per-gene protein-level classification
│   └── kaiju_contigs.tsv                 Per-contig majority-vote taxonomy
│
├── MAG Bins (DAS Tool consensus)
│   ├── bins/*.fa                         Final MAG sequences
│   ├── contig2bin.tsv                    Contig → MAG assignments
│   ├── bin_quality.tsv                   SCG completeness/redundancy
│   └── checkm2/quality_report.tsv       CheckM2 completeness + contamination
│
├── Mobile Genetic Elements
│   ├── genomad/                          Virus + plasmid proteins + gene annotations
│   ├── defensefinder/                    Anti-phage defense systems
│   ├── macsyfinder/                      Secretion + conjugation systems
│   ├── integrons/                        Integron + gene cassette annotations
│   └── genomic_islands/                  Island coordinates
│
└── WHAT IS MISSING (this document addresses)
    ├── KEGG Ortholog (KO) assignments
    ├── COG functional categories
    ├── GO term annotations
    ├── Pfam domain annotations (beyond Prokka's limited set)
    ├── CAZyme (carbohydrate-active enzyme) profiles
    ├── Secondary metabolite biosynthetic gene clusters
    ├── Peptidase/protease classification
    ├── Pathway completeness scoring (KEGG Modules, MetaCyc)
    ├── Biogeochemical cycling trait profiles
    └── Metabolic network reconstruction
```

### What Prokka Already Provides (and Its Limitations)

Prokka annotates genes against a hierarchy of databases: ISfinder, NCBI AMR, UniProtKB (organism-specific), then Pfam/TIGRFAM HMMs. However:

- **~49% of predicted proteins** in typical environmental metagenomes are annotated as "hypothetical protein" by Prokka
- Prokka's UniProtKB search uses a curated subset, not the full database
- No KEGG KO assignments, no COG categories, no GO terms in standard output
- No pathway-level analysis (individual gene annotations only)
- Pfam/TIGRFAM HMMs provide domain-level annotation but not metabolic context

The gap between "gene X encodes a dehydrogenase" (Prokka) and "this MAG can perform complete denitrification" (metabolic profiling) is what the tools in this document bridge.

### Goals for Metabolic Profiling Integration

1. **Assign standardized functional identifiers** (KO, COG, EC numbers) to predicted proteins
2. **Score pathway completeness** per MAG (which KEGG Modules are present/absent)
3. **Profile biogeochemical capabilities** (carbon fixation, nitrogen cycling, sulfur metabolism)
4. **Detect specialized metabolism** (CAZymes, secondary metabolites, peptidases)
5. **Enable cross-sample comparison** via standardized functional profiles
6. **Maintain shipboard deployment feasibility** (manageable database sizes, reasonable compute)

---

## Tool Rankings

Each category lists tools in order of preference. Tiers:

- **Recommended**: Install and integrate into the Nextflow pipeline
- **Alternative**: Good fallback; install if primary tool fails or for consensus
- **Investigate Later**: Promising but insufficient evidence or integration complexity
- **Not Recommended**: Poor fit for environmental metagenomes, unmaintained, or superseded

---

## 1. Gene Prediction Foundation

Before functional annotation, gene prediction quality directly determines annotation completeness. The pipeline currently uses Prokka, which calls Prodigal internally.

### Current: Prokka

- **Gene caller**: Prodigal (v2.6.3, single mode adapted for metagenomes with `--metagenome`)
- **Annotation**: Hierarchical BLAST/HMM search against curated databases
- **Strengths**: Fast, widely used, produces standard output formats (GFF3, GBK, FAA, FFN, FNA, TSV)
- **Limitations**: Prodigal's gene models are optimized for Illumina assemblies; 49% hypothetical proteins in environmental metagenomes; limited database breadth
- **Conda**: Already installed (`dana-mag-prokka`)
- **Output used downstream**: `.faa` (proteins), `.gff` (features), `.gbk` (GenBank format for IslandPath)

### Recommended Upgrade: Bakta

- **Purpose**: Next-generation genome annotation, designed as Prokka replacement
- **Gene caller**: Pyrodigal (Cython reimplementation of Prodigal, 3-4x faster)
- **Annotation databases**: AMRFinderPlus, COG, IS, Pfam-A, UniRef90 clusters, expert protein sequences
- **Key advantage**: Reduces hypothetical proteins from ~49% to ~24% by using comprehensive UniRef90 clusters with pre-computed HMM profiles
- **Environmental suitability**: Good — database is not taxon-biased; designed for any bacterial/archaeal genome
- **Additional features**:
  - Assigns consistent, stable locus tags (important for cross-study comparisons)
  - Detects sORFs (small open reading frames) often missed by Prodigal
  - Annotates ncRNAs (tRNA, rRNA, tmRNA, ncRNA, CRISPR arrays, riboswitches)
  - Outputs INSDC-compliant GenBank/EMBL files for direct submission
  - Produces JSON output with structured annotations (machine-parseable)
- **Database**: ~70 GB full database; ~2 GB "light" database (sufficient for most uses)
- **Speed**: ~7 min per genome (full DB) vs ~4 min for Prokka; Pyrodigal gene calling itself is 3-4x faster than Prodigal
- **Conda**: `bioconda::bakta`
- **Input**: FASTA (genome or MAG)
- **Output**: GFF3, GBK, FAA, FFN, TSV, JSON, EMBL, hypotheticals.faa, hypotheticals.tsv
- **Maintenance**: Active (v1.10+, 2024-2025)
- **Verdict**: **Strong upgrade candidate.** The reduction in hypothetical proteins directly improves all downstream functional annotation. The structured JSON output simplifies pipeline integration. However, the 70 GB full database is a consideration for shipboard deployment; the 2 GB light database may suffice.
- **Citation**: Schwengers et al. (2021) *Microbial Genomics* 7:000685

### Alternative Gene Callers

#### Pyrodigal (standalone)

- **Purpose**: Direct gene calling without the annotation layer
- **Method**: Cython reimplementation of Prodigal; identical results, 3-4x faster
- **Use case**: When you want gene predictions only (e.g., feeding directly into KofamScan or eggNOG-mapper without Prokka's annotation overhead)
- **Conda**: `bioconda::pyrodigal`
- **Metagenome mode**: `pyrodigal -p meta` for metagenome gene calling
- **Note**: Bakta uses Pyrodigal internally; standalone use only makes sense if skipping Bakta/Prokka entirely
- **Citation**: Larralde (2022) *JOSS* 7:4296

#### Prodigal-GV

- **Purpose**: Gene calling for giant viruses and other divergent organisms
- **Method**: Modified Prodigal with expanded genetic code support (code 15, 16, etc.)
- **Use case**: If metagenome contains giant viruses (Nucleocytoviricota) where standard Prodigal underpredicts
- **Conda**: `bioconda::prodigal-gv`
- **Note**: geNomad already uses prodigal-gv internally for viral contigs; no need to run separately unless analyzing viral contigs independently
- **Citation**: Camargo et al. (2023) *Nature Biotechnology*

### Recommendation

**Short term**: Continue using Prokka; all downstream tools accept Prokka `.faa` output. **Medium term**: Evaluate Bakta on a representative MAG set and compare hypothetical protein rates. If Bakta's light database reduces hypotheticals by >15 percentage points, switch. The `.faa` output format is identical, so downstream tools need no changes.

---

## 2. Core Annotation Tools

These tools assign standardized functional identifiers (KO, COG, EC, Pfam, GO) to predicted proteins. They form the annotation layer between gene prediction and pathway reconstruction.

### Recommended

#### KofamScan

- **Purpose**: KEGG Ortholog (KO) assignment via profile HMMs
- **Method**: Searches predicted proteins against KOfam — a library of ~25,000 KO-specific HMM profiles. Each profile has an adaptive score threshold calibrated to distinguish true orthologs from paralogs.
- **Key advantage**: Adaptive thresholds per KO family. Unlike BLAST-based methods that use a single E-value cutoff, each KO profile has its own calibrated threshold, dramatically reducing false positives for promiscuous protein families (e.g., ABC transporters, kinases).
- **Environmental suitability**: Excellent — KO profiles are built from diverse organisms, not taxon-biased
- **Speed**: ~69x faster than BlastKOALA for equivalent results; processes ~1M proteins in 2-4 hours on 16 cores
- **Accuracy**: Comparable to BlastKOALA/GhostKOALA but freely available and local (no web submission)
- **Conda**: `bioconda::kofamscan` (includes HMMER dependency)
- **Database**: KOfam profiles (~4 GB download; `ftp://ftp.genome.jp/pub/db/kofam/`)
- **Input**: Protein FASTA (`.faa` from Prokka/Bakta)
- **Output**: TSV with protein_id, KO, threshold, score, E-value; asterisk marks hits above adaptive threshold
- **Parallelization**: Built-in `--cpu` flag; embarrassingly parallel across proteins
- **Licensing**: **Free for academic use** — KOfam profiles are freely downloadable; no KEGG license required
- **Limitations**: Only assigns KOs; does not provide COG, GO, or Pfam annotations; requires separate pathway reconstruction step (MinPath, KEGG-Decoder) to interpret results
- **Maintenance**: Active (profiles updated with each KEGG release)
- **Verdict**: **Top pick for KO assignment.** Fastest, most accurate freely available method. The adaptive thresholds are a significant quality advantage over BLAST-based approaches. Essential foundation for KEGG Module completeness scoring.
- **Citation**: Aramaki et al. (2020) *Bioinformatics* 36:2251-2252

#### eggNOG-mapper v2

- **Purpose**: Orthology-based functional annotation — assigns COG categories, KEGG KOs, GO terms, EC numbers, Pfam domains, CAZy families, and KEGG pathways/modules in a single pass
- **Method**: Maps query proteins to eggNOG v5 orthologous groups via fast diamond/MMseqs2 search, then transfers pre-computed annotations from the best ortholog group
- **Key advantage**: One tool, many annotations. A single run provides COG, KEGG, GO, EC, Pfam, CAZy, and KEGG pathway/module assignments. Coverage spans >5000 species from all domains of life.
- **Environmental suitability**: Excellent — eggNOG orthologous groups include environmental and uncultured organisms
- **Speed**: ~2-4 hours for 1M proteins on 16 cores (diamond mode); ~10-20 min with MMseqs2 mode
- **Accuracy**: High precision for well-characterized orthologs; lower sensitivity than HMM-based methods (KofamScan, InterProScan) for divergent proteins
- **Annotation types provided**:
  - COG functional categories (A-Z)
  - KEGG KOs, pathways, modules, reactions
  - GO terms (molecular function, biological process, cellular component)
  - EC numbers (enzyme classification)
  - Pfam domains
  - CAZy families
  - KEGG BRITE hierarchies
  - Free-text functional description
- **Conda**: `bioconda::eggnog-mapper`
- **Database**: eggNOG v5 diamond/MMseqs2 database (~45 GB for full bacteria/archaea; ~12 GB for diamond DB only)
- **Input**: Protein FASTA (`.faa`) or CDS nucleotide FASTA
- **Output**: TSV with 22 annotation columns per protein
- **Modes**:
  - `diamond` — default, BLAST-like sensitivity
  - `mmseqs` — faster, slightly lower sensitivity
  - `hmmer` — highest sensitivity for divergent proteins, slower
- **Licensing**: **Fully free** — no KEGG license required; annotations are pre-computed from public databases
- **Limitations**: Orthology-based transfer can miss novel functional variants; less accurate for highly divergent environmental proteins; database is large (45 GB)
- **Maintenance**: Active (v2.1.12+, 2024-2025)
- **Verdict**: **Essential complement to KofamScan.** Provides the broadest annotation coverage in a single tool. COG categories are particularly valuable for comparing community functional profiles. Use alongside KofamScan — KofamScan has better KO precision, eggNOG-mapper provides breadth (COG, GO, EC, Pfam).
- **Citation**: Cantalapiedra et al. (2021) *Molecular Biology and Evolution* 38:5825-5829

### Alternative

#### InterProScan

- **Purpose**: Multi-database protein domain scanning — integrates 14+ member databases into a single search
- **Method**: Runs query proteins against all InterPro member databases simultaneously: Pfam, TIGRFAM, PANTHER, SMART, CDD, SUPERFAMILY, Gene3D, Hamap, PIRSF, PRINTS, ProSitePatterns, ProSiteProfiles, SFLD, NCBIfam
- **Key advantage**: Most comprehensive domain annotation available; integrates results from all major protein family databases with consistent InterPro accessions and GO term mappings
- **Environmental suitability**: Excellent — domain-based, no taxon bias
- **Speed**: **Slow** — 8-24 hours for 1M proteins on 16 cores (bottleneck: PANTHER, Gene3D); can be sped up by disabling slow analyses (`-appl Pfam,TIGRFAM,CDD`)
- **Accuracy**: Highest sensitivity for detecting protein domains, especially in divergent/novel proteins
- **Annotation types provided**:
  - InterPro accessions (IPR) with descriptions
  - GO terms mapped from domain hits
  - Pfam, TIGRFAM, PANTHER, SMART, CDD identifiers
  - Pathway annotations via MetaCyc/Reactome mappings
- **Conda**: Partial — `bioconda::interproscan` exists but is complex; Java-based, requires separate database download
- **Database**: ~90 GB total (all member databases); ~15 GB for Pfam+TIGRFAM+CDD subset
- **Input**: Protein FASTA (`.faa`)
- **Output**: TSV, GFF3, XML, or JSON
- **Parallelization**: Built-in `-cpu` flag; can split input for cluster submission
- **Licensing**: **Free** — all member databases are freely available
- **Limitations**: Very large database; slow runtime; complex installation (Java + multiple binaries); diminishing returns over eggNOG-mapper + KofamScan for most use cases
- **Maintenance**: Active (v5.72+, EBI)
- **Verdict**: Provides the most comprehensive domain annotation but at significant computational cost. Run only if Pfam/TIGRFAM domain-level detail is specifically needed beyond what eggNOG-mapper provides. Consider the Pfam+TIGRFAM+CDD subset to reduce runtime.
- **Citation**: Paysan-Lafosse et al. (2023) *Nucleic Acids Research* 51:D483-D489

#### dbCAN3 (run_dbcan)

- **Purpose**: Carbohydrate-active enzyme (CAZyme) annotation with substrate specificity prediction
- **Method**: Triple-tool consensus — HMMER (dbCAN HMMs), DIAMOND (CAZy DB), and eCAMI (substrate specificity). Reports hits found by >=2 of 3 methods.
- **Key advantage**: Not just "this is a glycoside hydrolase" but "this GH5 targets cellulose" — substrate specificity prediction is unique to dbCAN3
- **Environmental suitability**: Excellent — CAZymes are central to aquatic carbon cycling (polysaccharide degradation, biofilm formation, chitin utilization)
- **Specialized features**:
  - PUL (polysaccharide utilization loci) detection for Bacteroidetes
  - CGC (CAZyme gene cluster) detection for co-localized degradation operons
  - Substrate mapping via dbCAN-sub (subfamily-level substrate prediction)
  - Signal peptide detection (secreted vs intracellular CAZymes)
- **CAZyme classes annotated**:
  - GH (glycoside hydrolases) — polysaccharide degradation
  - GT (glycosyltransferases) — polysaccharide biosynthesis
  - PL (polysaccharide lyases) — polysaccharide cleavage
  - CE (carbohydrate esterases) — deacetylation
  - AA (auxiliary activities) — lignin/oxidative degradation
  - CBM (carbohydrate-binding modules) — substrate binding
- **Conda**: `bioconda::dbcan` (v4+) or pip install
- **Database**: ~2 GB (CAZy DB + dbCAN HMMs + dbCAN-sub + eCAMI models)
- **Input**: Protein FASTA (`.faa`) and/or genome FASTA (`.fna`)
- **Output**: TSV per method (HMMER, DIAMOND, eCAMI), consensus overview.txt, CGC clusters, substrate predictions
- **Speed**: ~30 min for 50k proteins on 8 cores
- **Licensing**: **Free**
- **Maintenance**: Active (v4.1+, 2024-2025)
- **Verdict**: Essential for aquatic/environmental metagenomics where carbon cycling is a primary research question. Substrate specificity predictions add ecological interpretability beyond generic CAZy family assignments. Lightweight and fast.
- **Citation**: Zheng et al. (2023) *Nucleic Acids Research* 51:W115-W121

#### MEROPS (batch BLAST/HMMER)

- **Purpose**: Peptidase and protease classification
- **Method**: BLAST or HMMER search against the MEROPS database — classifies peptidases into clans, families, and subfamilies with mechanistic class assignment (serine, cysteine, metallo, aspartic, threonine, mixed, unknown)
- **Environmental suitability**: Good — peptidases are important in marine nutrient cycling (protein remineralization, dissolved organic nitrogen release)
- **Database**: ~500 MB (MEROPS sequence library + HMM profiles)
- **Input**: Protein FASTA (`.faa`)
- **Output**: BLAST/HMMER tabular output with MEROPS family assignments
- **Access**: Database download from https://www.ebi.ac.uk/merops/download_list.shtml (free registration required)
- **Integration**: Can be searched via HMMER directly; DRAM includes MEROPS in its integrated pipeline
- **Maintenance**: Active (EBI-maintained, regular updates)
- **Note**: Not a standalone tool — requires wrapping HMMER/DIAMOND search against the MEROPS database. DRAM automates this.
- **Verdict**: Important for marine/aquatic metagenomics (protein turnover is a major nutrient flux). DRAM integrates MEROPS automatically; standalone search only needed if running DRAM is not planned.
- **Citation**: Rawlings et al. (2018) *Nucleic Acids Research* 46:D624-D632

#### antiSMASH v8

- **Purpose**: Biosynthetic gene cluster (BGC) detection — identifies gene clusters encoding secondary metabolite biosynthesis pathways
- **Method**: Profile HMM scanning for core biosynthetic enzymes + rule-based cluster boundary detection + comparative genomics against MIBiG reference clusters
- **Cluster types**: 81+ recognized types including:
  - NRPS (nonribosomal peptide synthetases)
  - PKS (polyketide synthases, types I-III)
  - RiPPs (ribosomally synthesized and post-translationally modified peptides)
  - Terpenes, siderophores, bacteriocins, aryl polyenes
  - Hybrid clusters (NRPS-PKS, etc.)
- **Environmental suitability**: Good — detects BGCs regardless of taxonomic origin; marine bacteria are prolific secondary metabolite producers
- **Key features**:
  - MIBiG comparison identifies known compound families
  - ClusterBlast finds similar clusters in GenBank
  - KnownClusterBlast matches experimentally characterized clusters
  - Generates interactive HTML visualization
  - Supports `--taxon bacteria` and `--taxon fungi` modes
- **Conda**: `bioconda::antismash` (v8.0+; complex install with many dependencies)
- **Database**: ~15 GB (Pfam, MIBiG, ClusterBlast, Resfam, TIGRFAM, etc.)
- **Input**: GenBank-format annotated genome (`.gbk` from Prokka/Bakta) or nucleotide FASTA
- **Output**: GenBank + JSON per region, HTML visualization, TSV summary of detected clusters
- **Speed**: 5-30 min per genome depending on BGC density
- **Per-MAG usage**: Run on individual MAGs (not raw assembly) for biologically meaningful results — antiSMASH expects single-organism input
- **Licensing**: **Free** (academic + commercial)
- **Limitations**: Designed for single genomes, not metagenomes directly; run per-MAG after binning; complex installation; large database; incomplete MAGs will have fragmented BGCs
- **Maintenance**: Active (v8.0, 2024; Blin lab, most-cited BGC tool)
- **Verdict**: Industry standard for secondary metabolite detection. Important for marine metagenomics (siderophores, toxins, antimicrobials). Run per-MAG. The complex installation and large database make it a lower priority than KofamScan/eggNOG-mapper for initial integration, but high value for specialized analyses.
- **Citation**: Blin et al. (2023) *Nucleic Acids Research* 51:W46-W50

### Investigate Later

#### PGAP (NCBI Prokaryotic Genome Annotation Pipeline)

- **Purpose**: NCBI's official genome annotation pipeline
- **Method**: GeneMarkS-2 gene prediction + comprehensive BLAST/HMM annotation against NCBI databases
- **Note**: Gold standard for NCBI submissions; full annotation pipeline but designed for isolate genomes, not MAGs. Docker-based, complex setup. Consider only for publication-quality annotation of high-completeness MAGs.
- **Citation**: Tatusova et al. (2016) *Nucleic Acids Research* 44:6614-6624

#### METABOLIC-G annotators

- **Purpose**: The annotation components within METABOLIC (see Section 3) can be run independently
- **Note**: Uses HMMER against custom curated HMM profiles for biogeochemical cycling genes; specialized for environmental metabolic trait detection. Better accessed through the full METABOLIC pipeline.

---

## 3. Integrated Metabolic Profilers

These tools combine multiple databases and annotation steps into unified pipelines, providing both individual gene annotations and higher-level metabolic summaries.

### Recommended

#### DRAM (Distilled and Refined Annotation of Metabolism)

- **Purpose**: Comprehensive metabolic annotation integrating 8+ databases with distillation into interpretable summaries
- **Method**: Searches predicted proteins against multiple databases simultaneously, then distills raw annotations into metabolic module completeness scores and visual summaries
- **Databases integrated**:
  - KEGG (KOfam HMMs) — requires separate KEGG license or use `--use_uniref`
  - UniRef90 (free alternative to KEGG)
  - Pfam (protein families)
  - dbCAN (CAZymes with substrate prediction)
  - MEROPS (peptidases)
  - VOGDB (virus orthologous groups)
  - Methyl (custom methylotrophy HMMs)
  - Custom HMM sets (extensible)
- **Key advantage**: The "distillation" step — transforms raw annotations into:
  - **Product TSV**: Per-gene annotation with best hits from all databases
  - **Genome statistics**: Gene counts, coding density, annotation completeness
  - **Module completeness**: KEGG Module step coverage fractions
  - **Metabolism summary**: Heatmap-ready matrix of metabolic capabilities per genome
  - **ETC (electron transport chain) summary**: Complex I-V completeness
- **Environmental suitability**: Excellent — designed specifically for MAG annotation; handles draft-quality genomes gracefully
- **DRAM-v**: Viral variant for annotating viral contigs (identifies auxiliary metabolic genes — AMGs)
- **Speed**: 2-8 hours for 100 MAGs on 16 cores (database-dependent)
- **Conda**: `bioconda::dram`
- **Database setup**: `DRAM-setup.py prepare_databases` (downloads all databases; ~100 GB with KEGG, ~60 GB without)
- **Input**: FASTA files (MAG genomes); calls Prodigal internally for gene prediction
- **Output**: Per-genome and community-level annotation TSVs, distillation summaries, metabolism heatmaps
- **KEGG licensing considerations**:
  - `DRAM --use_uniref` mode bypasses KEGG database requirement entirely
  - UniRef90 provides comparable functional coverage for most metabolic pathways
  - If KEGG license available, use it for best KO assignments; if not, UniRef90 mode is fully functional
- **Limitations**: Complex database setup; large disk footprint; KEGG license ambiguity (KOfam profiles are free, but DRAM's KEGG integration historically expected full KEGG subscription); slower than individual tools
- **Maintenance**: Active (v1.5+; Shaffer et al., maintained by community)
- **Verdict**: **Best integrated solution for MAG metabolic profiling.** The distillation step is uniquely valuable — it transforms thousands of raw annotations into interpretable metabolic summaries. Particularly well-suited for our pipeline since it accepts MAG FASTAs directly and handles draft-quality genomes. Run on DAS_Tool consensus bins.
- **Citation**: Shaffer et al. (2020) *Nucleic Acids Research* 48:8883-8900

#### METABOLIC

- **Purpose**: Biogeochemical trait profiling — systematically profiles metabolic and biogeochemical functional traits from genomes and metagenomes
- **Method**: HMMER search against 146 curated HMM profiles covering biogeochemical cycling pathways, followed by pathway module completeness scoring
- **Two modes**:
  - **METABOLIC-G** (genome-centric): Annotates individual MAGs; reports metabolic capability per genome
  - **METABOLIC-C** (community-centric): Combines MAG annotations with read coverage data to estimate community-level metabolic potential and element cycling rates
- **Biogeochemical cycles profiled**:
  - Carbon: fixation (CBB, rTCA, Wood-Ljungdahl, 3-HP, DC/HB), degradation (chitin, cellulose, starch), fermentation, methanogenesis/methanotrophy
  - Nitrogen: fixation, nitrification, denitrification, DNRA, anammox, assimilation
  - Sulfur: oxidation (Sox), reduction (Dsr, Apr), thiosulfate disproportionation, DMSP cycling
  - Iron: Fe(II) oxidation, Fe(III) reduction
  - Hydrogen: hydrogenases (NiFe, FeFe)
  - Other: arsenic, selenium, halogenation
- **Key advantage**: Purpose-built for environmental microbial ecology — directly answers "what biogeochemical transformations can this community perform?" rather than generic functional annotation
- **Community-level outputs** (METABOLIC-C):
  - Biogeochemical cycling diagrams with gene abundances
  - Community metabolic network visualizations
  - Coverage-weighted metabolic potential per pathway
- **Conda**: Not standard — GitHub install with Perl + R dependencies
- **Database**: Curated HMM profiles bundled (~200 MB); also runs KEGG/Pfam annotation internally
- **Input**: MAG FASTAs + BAM files (for community mode); calls Prodigal internally
- **Output**: Per-genome functional trait table, community metabolism diagrams, R-generated visualizations
- **Speed**: 4-12 hours for 100 MAGs on 16 cores
- **Limitations**: Perl-based, complex installation; R dependencies for visualization; less actively maintained than DRAM; narrower scope (biogeochemical traits only, not general-purpose annotation)
- **Maintenance**: Moderate (v4.0; last major update 2021-2022)
- **Verdict**: **Best tool specifically for biogeochemical cycling analysis.** If the primary research question is "what element cycles does this community drive?", METABOLIC-C provides the most ecologically interpretable output. Complements DRAM (which is broader but less specialized for biogeochemistry). The community-centric mode using read coverage is particularly powerful for quantitative ecology.
- **Citation**: Zhou et al. (2022) *Microbiome* 10:33

### Alternative

#### MicrobeAnnotator

- **Purpose**: Multi-database consensus annotation with KEGG module completeness
- **Method**: Cascading search strategy — queries SwissProt first (highest quality), then TrEMBL, then RefSeq KO, then KOfam HMMs. Takes the best hit from the highest-quality database that returns a result.
- **Key advantage**: Annotation confidence hierarchy — SwissProt annotations are trusted over TrEMBL over RefSeq, providing built-in quality stratification
- **Annotation types**:
  - KO assignments (from multiple sources with confidence ranking)
  - EC numbers
  - GO terms
  - Pfam domains
  - KEGG Module completeness heatmaps (built-in visualization)
- **Conda**: `pip install microbeannotator` (conda: partial)
- **Database**: ~30 GB (SwissProt + TrEMBL + RefSeq + KOfam)
- **Input**: Protein FASTA (`.faa`)
- **Output**: Annotation TSV, KEGG Module completeness matrix, pathway heatmaps (matplotlib/seaborn)
- **Speed**: 3-6 hours for 100 MAGs on 16 cores
- **Licensing**: **Free** (uses only public databases)
- **Limitations**: Less actively maintained than DRAM; does not integrate CAZy/MEROPS/VOGDB; Python-based but complex dependencies
- **Maintenance**: Moderate (last update 2022-2023)
- **Verdict**: Good alternative if DRAM's database setup proves too complex. The cascading confidence strategy is elegant, but DRAM provides broader database coverage and better distillation.
- **Citation**: Ruiz-Perez et al. (2021) *BMC Bioinformatics* 22:11

#### Mantis

- **Purpose**: Consensus-based protein annotation with highest reported F1 score
- **Method**: Annotates proteins against multiple reference databases (NCBI COG, Pfam, KOfam, TIGRfam, NPFAM, eggNOG) and applies consensus integration with weighted scoring
- **Key advantage**: Highest F1 score reported in benchmark comparisons — 0.131 above eggNOG-mapper in the CAMI II challenge
- **Consensus approach**: When databases disagree, Mantis applies a weighted-voting algorithm considering database reliability, coverage, and annotation specificity
- **Conda**: `bioconda::mantis_pfa` (or pip)
- **Database**: Modular — download only the databases you need (~5-50 GB depending on selection)
- **Input**: Protein FASTA (`.faa`)
- **Output**: Consensus annotation TSV with identifiers from all queried databases
- **Speed**: Comparable to eggNOG-mapper
- **Licensing**: **Free**
- **Limitations**: Less widely adopted than DRAM/eggNOG-mapper; smaller user community; output format less standardized
- **Maintenance**: Active (2023-2024)
- **Verdict**: Intriguing for its consensus approach and benchmark performance. Consider if annotation accuracy is more important than ecosystem/tool maturity. The consensus strategy could reduce false annotations on novel environmental proteins.
- **Citation**: Queirós et al. (2021) *GigaScience* 10:giab042

### Investigate Later

#### Prokka2 / PGAP-lite

- **Status**: In development; will integrate NCBI's PGAP annotation logic into a Prokka-like interface
- **Note**: Could unify gene prediction + functional annotation in a single tool. Monitor for release.

#### MetaErg 2

- **Purpose**: Metagenome annotation pipeline
- **Method**: Integrates gene calling + functional annotation + metabolic profiling in one workflow
- **Note**: Uses diamond against KEGG, COG, TIGRFAM; produces KEGG module completeness. Less mature than DRAM.
- **Citation**: Dong & Strous (2019) *Frontiers in Genetics* 10:999

---

## 4. Pathway Reconstruction

These tools convert gene-level annotations (KO numbers, EC numbers) into pathway-level completeness scores, answering "which metabolic pathways are complete/incomplete in this MAG?"

### Recommended

#### MinPath

- **Purpose**: Parsimony-based pathway reconstruction — infers the minimum set of pathways needed to explain observed gene annotations
- **Method**: Given a set of KO or EC numbers, MinPath solves an integer programming problem to find the smallest set of KEGG or MetaCyc pathways that accounts for all observed annotations
- **Key advantage**: Avoids pathway inflation — tools that simply check "any gene in pathway X is present" massively overcount active pathways. MinPath asks "what is the minimum set of pathways consistent with these genes?" which is biologically more realistic for incomplete MAGs.
- **Environmental suitability**: Excellent — particularly important for draft-quality MAGs where gene absence could be real (pathway not present) or artifactual (assembly gap)
- **Conda**: `bioconda::minpath`
- **Input**: List of KO numbers or EC numbers (from KofamScan/eggNOG-mapper output)
- **Output**: Inferred pathways with member genes, pathway completeness status
- **Speed**: Seconds to minutes per genome
- **Databases**: KEGG pathways or MetaCyc pathways (bundled with MinPath)
- **Limitations**: Binary (pathway present/absent), not quantitative; does not score partial completeness; requires pre-computed KO/EC assignments
- **Maintenance**: Moderate (stable, infrequently updated)
- **Verdict**: **Essential post-processing step** for any KO-based annotation. Without parsimony filtering, pathway enumeration from environmental metagenomes is severely inflated. Fast, lightweight, easy to integrate.
- **Citation**: Ye & Doak (2009) *PLoS Computational Biology* 5:e1000465

#### KEGG Module Completeness Scoring

- **Purpose**: Score completeness of KEGG Modules (functional units smaller than pathways)
- **Method**: KEGG Modules are defined as boolean expressions of KOs (e.g., "K00844 + K00845 + (K00886,K00918)"). Scoring evaluates what fraction of required KOs are present in each module.
- **Multiple implementations**:
  - **KEGGDecoder** — Python script; reads KofamScan output, scores ~150 curated modules, outputs TSV + heatmap
  - **DRAM distillation** — built into DRAM (see Section 3); scores modules as part of the distillation step
  - **anvi'o `anvi-estimate-metabolism`** — part of the anvi'o ecosystem (see below)
  - **MicrobeAnnotator** — scores modules as part of its annotation pipeline
  - **Custom scripting** — KEGG module definitions can be parsed programmatically (available via KEGG REST API)
- **Key advantage**: Modules are smaller, better-defined functional units than pathways. "Complete module for dissimilatory nitrate reduction" is more informative and less error-prone than "denitrification pathway present."
- **Verdict**: Every metabolic profiling pipeline should include module completeness scoring. DRAM and KEGGDecoder both provide this; no separate tool needed if using either.

### Alternative

#### anvi'o Metabolism Module

- **Purpose**: KEGG module completeness estimation within the anvi'o ecosystem
- **Method**: `anvi-run-kegg-kofams` assigns KOs via KofamScan; `anvi-estimate-metabolism` scores KEGG modules per genome/bin; produces stepwise completeness scores with customizable thresholds
- **Key advantage**: Deep integration with anvi'o's pangenomic, phylogenomic, and visualization capabilities. If already using anvi'o for other analyses, this is the natural choice.
- **Features**:
  - User-defined module definitions (extend beyond KEGG)
  - Stepwise completeness (not just binary)
  - Enrichment analysis across genome groups
  - Interactive visualization in anvi'o interface
- **Conda**: `bioconda::anvio` (large install, ~2 GB with all dependencies)
- **Database**: KEGG KOfam profiles + KEGG module definitions (downloaded via `anvi-setup-kegg-data`)
- **Input**: anvi'o contigs database (requires format conversion from FASTA)
- **Output**: TSV module completeness scores, enrichment results
- **Limitations**: Requires buying into the anvi'o ecosystem (contigs database format, specific workflows); large installation; steeper learning curve
- **Maintenance**: Active (Meren lab, well-maintained)
- **Verdict**: Excellent if anvi'o is part of the workflow. For a Nextflow-based pipeline like danaSeq, the format conversion overhead and ecosystem lock-in make standalone tools (KofamScan + MinPath/KEGGDecoder) more practical.
- **Citation**: Eren et al. (2021) *Nature Microbiology* 6:3-6

#### Pathway Tools / MetaCyc

- **Purpose**: Metabolic pathway reconstruction using the MetaCyc database — the gold standard for curated pathway definitions
- **Method**: PathoLogic algorithm predicts metabolic pathways from annotated genomes using MetaCyc's 2969 pathways (7x KEGG's ~400 reference pathways), then constructs organism-specific Pathway/Genome Databases (PGDBs) with metabolic flux models
- **Key advantages**:
  - 2969 pathways vs KEGG's ~400 — far greater pathway resolution
  - Experimentally curated with literature evidence codes
  - Metabolic flux balance analysis capability
  - Transport reaction modeling
  - Metabolic network gap-filling algorithms
- **Conda**: Not available — standalone Java application with GUI + CLI
- **Database**: MetaCyc (free academic; ~2 GB)
- **Input**: GenBank/GFF annotated genomes with EC number assignments
- **Output**: Pathway/Genome Database (PGDB) with predicted pathways, reactions, metabolites, transport
- **Speed**: 10-30 min per genome
- **Licensing**: **Free for academic research** (requires registration); commercial license available
- **Limitations**: Java GUI-oriented (CLI exists but less documented); not designed for batch processing of hundreds of MAGs; requires EC number annotations (from eggNOG-mapper or similar); PGDB format is proprietary
- **Maintenance**: Active (SRI International, regular MetaCyc updates)
- **Verdict**: Gold standard for pathway curation quality, but integration into a Nextflow pipeline is impractical. Best used for targeted analysis of high-quality MAGs of particular interest, not bulk processing. Consider for publication-quality pathway analysis of key organisms.
- **Citation**: Karp et al. (2024) *Nucleic Acids Research* 52:D649-D656

### Investigate Later

#### HUMAnN3

- **Purpose**: Read-level functional profiling — quantifies pathway abundance and coverage directly from metagenomic reads without assembly
- **Method**: Maps reads to pangenomes (ChocoPhlAn) and protein databases (UniRef90) to quantify gene families and pathways per sample
- **Key advantage**: No assembly required — works on raw reads; provides quantitative abundance estimates (RPK) per pathway per sample
- **Output**: Gene family abundances (RPK), pathway abundances, pathway coverage, stratified by contributing organisms
- **Conda**: `bioconda::humann`
- **Database**: ChocoPhlAn pangenomes (~15 GB) + UniRef90 diamond DB (~30 GB)
- **Speed**: 2-8 hours per sample on 16 cores
- **Limitations**: Reference-dependent (novel organisms have low mapping rates in environmental samples — typically 10-40% for marine metagenomes); designed for short reads (long-read compatibility not officially supported); large databases
- **Note**: Complementary to MAG-based annotation — provides sample-level pathway quantification that MAG-based approaches cannot (since not all reads bin into MAGs). Consider for comparative analysis across sampling sites/timepoints.
- **Verdict**: Not a priority for the MAG pipeline but valuable as a complementary analysis for sample-level functional comparisons. Low mapping rates in aquatic environmental metagenomes limit utility.
- **Citation**: Beghini et al. (2021) *eLife* 10:e65088

#### GhostKOALA / BlastKOALA

- **Purpose**: KEGG's official KO assignment web services
- **Method**: GHOSTX (GhostKOALA) or BLAST (BlastKOALA) against KEGG GENES database
- **Limitations**: Web-only (no local install); 500K sequence limit per submission; slow turnaround (hours to days); cannot be integrated into Nextflow
- **Note**: Useful for validation — compare KofamScan results against KEGG's own assignment for a subset of proteins
- **Verdict**: Not suitable for pipeline integration; use for spot-checking KofamScan accuracy.
- **Citation**: Kanehisa et al. (2016) *Journal of Molecular Biology* 428:726-731

---

## 5. Existing Pipeline Integration Patterns

How other metagenome pipelines integrate metabolic profiling — patterns to learn from.

### nf-core/mag (v3.2+)

The nf-core community metagenome pipeline includes:

- **Gene calling**: Prodigal (metagenome mode) or Prokka (optional)
- **Functional annotation**: None built-in (by design — delegates to nf-core/funcscan)
- **Taxonomy**: GTDB-Tk for MAG classification
- **Integration pattern**: Outputs MAG FASTAs + protein FASTAs as pipeline outputs; downstream functional annotation is a separate pipeline

**Lesson**: nf-core separates assembly/binning from functional annotation. This is a clean architectural choice — the MAG pipeline produces genomes, a separate pipeline annotates them.

### nf-core/funcscan (v2.0+)

Dedicated functional annotation pipeline in the nf-core ecosystem:

- **Annotation tools**: Bakta, Prokka, Prodigal, Pyrodigal (gene prediction); AMRFinderPlus, DeepARG, fARGene (AMR); antiSMASH, DeepBGC, GECCO (BGCs); ABRicate (general screening)
- **No metabolic profiling**: Focused on AMR + BGC detection, not general metabolic pathway analysis
- **Integration pattern**: Takes genome FASTAs as input; runs tools in parallel; produces per-genome annotation outputs

**Lesson**: Even nf-core does not have a general-purpose metabolic profiling pipeline. This is a genuine gap in the Nextflow ecosystem.

### ATLAS (v2.18+)

Metagenome-atlas, a Snakemake-based pipeline:

- **Annotation**: eggNOG-mapper (default) or DRAM (optional)
- **Gene calling**: Prodigal (metagenome mode)
- **Genomes**: Outputs MAG FASTAs with GTDB-Tk taxonomy
- **Integration pattern**: eggNOG-mapper runs on co-assembly gene catalog; per-MAG annotation is done via gene-to-bin mapping

**Lesson**: ATLAS annotates the co-assembly gene catalog once, then maps annotations to bins. This is efficient — avoid re-annotating the same genes per-MAG.

### MetaWRAP

- **Annotation**: Prokka + custom BRITE/KO scripts
- **Note**: Includes a `reassemble_bins` module and `classify_bins` but minimal functional annotation
- **Lesson**: MetaWRAP focuses on binning quality, not functional annotation

### anvi'o Metagenomics Workflow

- **Annotation**: KofamScan (via `anvi-run-kegg-kofams`), COGs (via `anvi-run-ncbi-cogs`), Pfam (via `anvi-run-pfams`), HMMs (via `anvi-run-hmms`)
- **Metabolism**: `anvi-estimate-metabolism` for KEGG module completeness
- **Integration pattern**: Everything goes into an anvi'o contigs database; annotations are layers that can be queried interactively

**Lesson**: anvi'o's layered annotation approach is powerful for exploration but requires format conversion. For a batch-processing Nextflow pipeline, standalone tools are more practical.

### Common Pattern: Annotate Once, Map to Bins

Most mature pipelines annotate the co-assembly gene catalog once (all predicted proteins from the full assembly), then map annotations to individual MAGs via contig-to-bin assignments. This avoids redundant annotation of the same genes and is computationally efficient.

**Our approach should follow this pattern:**

```
Assembly proteins (Prokka .faa)     ← Already produced
          │
    Functional annotation            ← KofamScan + eggNOG-mapper (one-time, full assembly)
          │
    contig2bin mapping               ← Already produced (DAS_Tool)
          │
    Per-MAG annotation tables        ← Filter by bin membership
          │
    Pathway scoring per MAG          ← MinPath / KEGG module scoring
```

---

## 6. Database Requirements

### Size and Licensing Summary

| Database | Size | License | Required By | Update Frequency |
|----------|------|---------|-------------|-----------------|
| KOfam profiles | ~4 GB | Free (academic) | KofamScan | Quarterly (KEGG releases) |
| eggNOG v5 (diamond) | ~12 GB | Free | eggNOG-mapper | ~Annual |
| eggNOG v5 (full) | ~45 GB | Free | eggNOG-mapper (hmmer mode) | ~Annual |
| InterPro (full) | ~90 GB | Free | InterProScan | Monthly |
| InterPro (Pfam+TIGRFAM+CDD) | ~15 GB | Free | InterProScan (subset) | Monthly |
| dbCAN | ~2 GB | Free | run_dbcan | ~Biannual |
| MEROPS | ~500 MB | Free (registration) | HMMER/DIAMOND search | ~Annual |
| antiSMASH | ~15 GB | Free | antiSMASH | With releases |
| UniRef90 (diamond) | ~30 GB | Free | DRAM, HUMAnN3 | Monthly |
| DRAM databases (all) | ~100 GB | Mixed (see below) | DRAM | Variable |
| DRAM databases (no KEGG) | ~60 GB | Free | DRAM `--use_uniref` | Variable |
| MetaCyc | ~2 GB | Free (academic) | Pathway Tools | ~Annual |
| Bakta (full) | ~70 GB | Free | Bakta | With releases |
| Bakta (light) | ~2 GB | Free | Bakta | With releases |

### KEGG Licensing — Critical Consideration

**The KEGG database itself requires a paid license** for bulk download and local installation (~$5000/year academic, ~$25,000/year commercial). However:

- **KOfam HMM profiles**: **Free to download** from KEGG FTP — this is what KofamScan uses
- **KEGG pathway definitions**: **Free via KEGG REST API** — pathway structures, module definitions, and maps are freely accessible
- **KEGG GENES database**: **Requires license** — the sequence database used by BlastKOALA/GhostKOALA
- **DRAM with KEGG**: Historically required the KEGG GENES database for full functionality; **DRAM `--use_uniref` bypasses this entirely** using free UniRef90 clusters

**Practical resolution for danaSeq:**

1. **KofamScan** — fully functional with free KOfam profiles (**no license needed**)
2. **eggNOG-mapper** — pre-computed annotations from free eggNOG database; provides KO assignments derived from orthology (**no license needed**)
3. **DRAM** — use `--use_uniref` mode for free UniRef90-based annotation (**no license needed**)
4. **KEGG module scoring** — module definitions freely available via REST API or bundled with KEGGDecoder (**no license needed**)

**Bottom line: Full metabolic profiling is achievable without any KEGG license.**

### Database Download Procedures

#### KOfam Profiles

```bash
# Download KOfam profiles (~4 GB)
mkdir -p databases/kofam
wget -P databases/kofam ftp://ftp.genome.jp/pub/db/kofam/ko_list.gz
wget -P databases/kofam ftp://ftp.genome.jp/pub/db/kofam/profiles.tar.gz
gunzip databases/kofam/ko_list.gz
tar xzf databases/kofam/profiles.tar.gz -C databases/kofam/

# Verify
ls databases/kofam/profiles/ | wc -l    # Should be ~25,000 HMM files
head databases/kofam/ko_list             # KO → threshold mapping
```

#### eggNOG-mapper Databases

```bash
# Using eggNOG-mapper's built-in downloader (~12-45 GB depending on mode)
mamba run -p conda-envs/dana-mag-emapper \
    download_eggnog_data.py --data_dir databases/eggnog -y

# Diamond-only mode (smaller, faster)
download_eggnog_data.py --data_dir databases/eggnog -y -D

# Verify
ls -lh databases/eggnog/eggnog.db        # Should be ~45 GB (full) or ~12 GB (diamond)
ls -lh databases/eggnog/eggnog_proteins.dmnd
```

#### dbCAN

```bash
# Using run_dbcan's built-in downloader (~2 GB)
mamba run -p conda-envs/dana-mag-dbcan \
    dbcan_build --db_dir databases/dbcan --clean

# Or manual download
mkdir -p databases/dbcan
cd databases/dbcan
wget https://bcb.unl.edu/dbCAN2/download/dbCAN-HMMdb-V12.txt
wget https://bcb.unl.edu/dbCAN2/download/CAZyDB.fa
wget https://bcb.unl.edu/dbCAN2/download/tcdb.fa
wget https://bcb.unl.edu/dbCAN2/download/stp.hmm
hmmpress dbCAN-HMMdb-V12.txt
diamond makedb --in CAZyDB.fa -d CAZy
```

#### DRAM Databases

```bash
# Full setup (with UniRef90, no KEGG license needed)
mamba run -p conda-envs/dana-mag-dram \
    DRAM-setup.py prepare_databases \
    --output_dir databases/dram \
    --skip_uniref false \
    --threads 16

# Verify
DRAM-setup.py print_config
```

### Total Disk Requirements

**Minimum viable set** (KofamScan + eggNOG-mapper diamond + dbCAN):

```
KOfam profiles:           ~4 GB
eggNOG diamond DB:       ~12 GB
dbCAN:                    ~2 GB
─────────────────────────────
Total:                   ~18 GB
```

**Recommended set** (add DRAM with UniRef90):

```
Minimum viable:          ~18 GB
DRAM (UniRef90 mode):    ~60 GB
─────────────────────────────
Total:                   ~78 GB
```

**Full set** (all tools):

```
Recommended:             ~78 GB
antiSMASH:               ~15 GB
InterProScan subset:     ~15 GB
Bakta full:              ~70 GB
─────────────────────────────
Total:                  ~178 GB
```

**For shipboard deployment**: The minimum viable set (18 GB) or recommended set (78 GB) is practical. The full set (178 GB) requires significant storage but is feasible with modern NVMe drives.

---

## 7. Computational Requirements

### Per-Tool Resource Estimates

Estimates based on a typical metagenome co-assembly with ~500k predicted proteins and ~50-100 MAGs.

| Tool | CPUs | RAM | Time (16 cores) | Parallelization | GPU |
|------|------|-----|-----------------|-----------------|-----|
| KofamScan | 1-32 | 4-8 GB | 2-4 hours | Per-protein (built-in) | No |
| eggNOG-mapper (diamond) | 1-32 | 8-16 GB | 2-4 hours | Per-protein (built-in) | No |
| eggNOG-mapper (mmseqs) | 1-32 | 16-32 GB | 20-40 min | Per-protein (built-in) | No |
| InterProScan (full) | 1-32 | 16-32 GB | 8-24 hours | Per-protein (built-in) | No |
| InterProScan (Pfam+TF+CDD) | 1-32 | 8-16 GB | 2-6 hours | Per-protein (built-in) | No |
| dbCAN3 | 1-16 | 4-8 GB | 30-60 min | Per-protein (built-in) | No |
| antiSMASH (per MAG) | 1-8 | 4-8 GB | 5-30 min/MAG | Per-MAG (external) | No |
| DRAM (annotate) | 1-32 | 16-64 GB | 4-8 hours | Per-protein (built-in) | No |
| DRAM (distill) | 1-4 | 4-8 GB | 5-15 min | Sequential | No |
| METABOLIC-G | 1-32 | 8-16 GB | 4-8 hours | Per-MAG (built-in) | No |
| METABOLIC-C | 1-32 | 16-32 GB | 6-12 hours | Per-MAG (built-in) | No |
| MinPath | 1 | <1 GB | Seconds/MAG | Per-MAG (external) | No |
| Bakta (per MAG) | 1-16 | 8-16 GB | 5-15 min/MAG | Per-protein (built-in) | No |
| Mantis | 1-32 | 8-16 GB | 2-4 hours | Per-protein (built-in) | No |
| MicrobeAnnotator | 1-32 | 16-32 GB | 3-6 hours | Per-protein (built-in) | No |

### Bottleneck Analysis

For a pipeline processing 500k proteins from a co-assembly:

```
CRITICAL PATH (sequential dependencies)
Assembly → Prokka (gene prediction) → Functional annotation → Pathway scoring → Visualization
                                      ↑
                                  BOTTLENECK
                                  (2-8 hours)
```

**The annotation step is the bottleneck.** Gene prediction (Prokka) is already cached. Pathway scoring (MinPath) is fast. The annotation step (KofamScan, eggNOG-mapper) takes 2-8 hours but can run in parallel.

### Parallelization Strategy

```
                  Prokka .faa (all assembly proteins)
                          │
              ┌───────────┼───────────┐
              │           │           │
         KofamScan   eggNOG-mapper  dbCAN3        ← Run in parallel (independent)
              │           │           │
              └───────────┼───────────┘
                          │
                  Merge annotations
                          │
                  Map to bins (contig2bin.tsv)
                          │
              ┌───────────┼───────────┐
              │           │           │
         MinPath    KEGG modules   Community       ← Run in parallel (per-MAG)
                                  profiles
```

All three core annotation tools (KofamScan, eggNOG-mapper, dbCAN3) can run simultaneously on the same input `.faa` file. They are independent and can occupy different CPU/RAM resources.

### HPC/Cluster Considerations

- **KofamScan**: Memory-efficient; can run on any node with 8 GB RAM
- **eggNOG-mapper**: Needs ~16 GB RAM for diamond DB loading; ~32 GB for MMseqs2 mode
- **DRAM**: Needs ~64 GB RAM when using multiple databases simultaneously; benefits from fast SSD for database I/O
- **antiSMASH**: Embarrassingly parallel per-MAG; ideal for array jobs on HPC
- **No GPU requirements**: None of the recommended tools benefit from GPU acceleration

### Shipboard Deployment Considerations

The danaSeq pipeline is designed for shipboard deployment with constrained resources:

- **CPU**: 32 cores available (current pipeline uses 32x1 parallelization for real-time processing)
- **RAM**: 256 GB available (shipboard profile)
- **Storage**: NVMe drives with reasonable capacity
- **Network**: Intermittent or no internet connectivity during cruises

**Implications for metabolic profiling:**

1. **All databases must be pre-downloaded** before departure — no on-demand database access
2. **KofamScan + eggNOG-mapper** fit within 18 GB database storage and 16 GB RAM
3. **DRAM with UniRef90** adds 60 GB storage but stays within 64 GB RAM
4. **InterProScan full** (90 GB) may be impractical; the Pfam+TIGRFAM+CDD subset (15 GB) is viable
5. **Annotation runs can proceed overnight** — not time-critical like real-time read processing
6. **No internet dependency** — all tools run locally with pre-downloaded databases

---

## 8. Visualization Tools

### Recommended

#### KEGG-Decoder

- **Purpose**: Heatmap visualization of KEGG Module completeness across genomes
- **Method**: Reads KO assignments, scores ~150 curated biogeochemically-relevant modules, generates heatmap
- **Key advantage**: Pre-curated module set focused on environmentally relevant pathways (carbon, nitrogen, sulfur, energy metabolism, vitamin synthesis, organic nitrogen, misc)
- **Output**: Static heatmap (matplotlib) + TSV of module completeness scores
- **Conda**: pip install (`pip install kegg-decoder`)
- **Input**: TSV of genome × KO presence/absence (from KofamScan output)
- **Speed**: Seconds (visualization only)
- **Verdict**: **Quick, informative visualization** that answers "what can each MAG do?" at a glance. Essential output for any metabolic profiling analysis.
- **Citation**: Graham et al. (2018) *ISME Communications* 8:16

#### ggkegg (R/Bioconductor)

- **Purpose**: KEGG pathway map visualization using ggplot2 grammar
- **Method**: Downloads KEGG pathway maps, overlays gene presence/expression data using customizable ggplot2 aesthetics
- **Key advantage**: Full ggplot2 integration — combine with other tidyverse data; publication-quality static figures; flexible aesthetics (color by completeness, abundance, taxonomy)
- **Output**: KEGG pathway diagrams with genes colored by presence/abundance
- **Conda**: R package (`BiocManager::install("ggkegg")`)
- **Input**: KO assignments + KEGG pathway definitions
- **Verdict**: Best option for publication-quality KEGG pathway figures. Requires R familiarity.
- **Citation**: Ishigami & Kume (2024) *Bioinformatics* 40:btae060

### Alternative

#### KEGGCharter

- **Purpose**: Automated KEGG pathway map coloring
- **Method**: Takes KO + taxonomy assignments, colors KEGG reference pathway maps by taxonomic origin, generates multi-genome pathway comparison maps
- **Key advantage**: Automated — takes annotations in, produces colored pathway maps out; multi-genome comparison built-in
- **Output**: PDF/PNG KEGG pathway maps with genes colored by genome/taxon
- **Conda**: pip install (`pip install keggcharter`)
- **Input**: TSV with protein_id, KO, taxon columns
- **Limitations**: Depends on KEGG REST API for pathway maps (needs internet access for initial map download; maps can be cached)
- **Verdict**: Good for automated pathway visualization. Cache maps before shipboard deployment.
- **Citation**: Sequeira et al. (2022) *Computational and Structural Biotechnology Journal* 20:1798-1810

#### iPath3

- **Purpose**: Interactive metabolic pathway explorer — web-based visualization of metabolic networks
- **Method**: Maps KO/EC numbers onto a global metabolic network, highlights complete/incomplete pathways
- **Output**: Interactive SVG metabolic maps
- **Access**: Web-based (https://pathways.embl.de) + downloadable SVG templates
- **Limitations**: **Web-based** — requires internet. SVG templates can be downloaded for offline use but lose interactivity.
- **Verdict**: Excellent for exploration and presentations. Download SVG templates before deployment for offline use.
- **Citation**: Darzi et al. (2018) *Nucleic Acids Research* 46:W510-W514

### Investigate Later

#### anvi'o Interactive Interface

- **Purpose**: Interactive exploration of metabolic annotations in genomic context
- **Note**: Powerful but requires anvi'o ecosystem buy-in. Best for detailed exploration of individual MAGs, not batch visualization.

#### METABOLIC Community Diagrams

- **Purpose**: Biogeochemical cycling diagrams with gene abundance overlays
- **Note**: Auto-generated by METABOLIC-C; R-based visualization. Publication-ready but requires METABOLIC pipeline.

#### MicrobiomeAnalyst

- **Purpose**: Web-based multi-omics visualization platform
- **Note**: Accepts functional profiles for interactive exploration. Web-based — not suitable for shipboard use.

---

## 9. Recommendations for danaSeq Integration

### Evaluation Criteria

Tools are ranked based on:

1. **Environmental suitability** — performance on aquatic/marine metagenomes
2. **ONT long-read compatibility** — validated on long-read assemblies
3. **Draft MAG tolerance** — handles 50-90% complete genomes gracefully
4. **Shipboard feasibility** — reasonable database/compute requirements for deployment
5. **Nextflow integrability** — conda-installable, CLI-based, standard I/O formats
6. **Free licensing** — no commercial database subscriptions required
7. **Community adoption** — widely used, well-maintained, citable

### Tier 1: Core Annotation (Implement First)

These three tools form the minimum viable metabolic profiling layer:

| Priority | Tool | Purpose | Database | Why |
|----------|------|---------|----------|-----|
| 1 | **KofamScan** | KEGG KO assignment | 4 GB (free) | Best KO precision; adaptive thresholds; fast; foundation for pathway scoring |
| 2 | **eggNOG-mapper v2** | COG + GO + EC + Pfam | 12 GB (free) | Broadest annotation in one pass; COG categories for community comparison |
| 3 | **dbCAN3** | CAZyme profiling | 2 GB (free) | Essential for carbon cycling; substrate specificity unique to this tool |

**Total database requirement**: ~18 GB
**Total compute**: ~4-6 hours on 16 cores (running in parallel)
**What you get**: Per-protein KO, COG, GO, EC, Pfam, CAZy annotations; per-MAG KEGG module completeness; CAZyme substrate profiles

### Tier 2: Integrated Profiling (Implement Second)

| Priority | Tool | Purpose | Database | Why |
|----------|------|---------|----------|-----|
| 4 | **DRAM** | Integrated annotation + distillation | 60 GB (free, UniRef90 mode) | Best summary output; MEROPS + VOGDB + distillation step |
| 5 | **MinPath** | Parsimony pathway reconstruction | Bundled (~50 MB) | Prevents pathway inflation; essential for draft MAGs |
| 6 | **KEGG-Decoder** | Module completeness heatmaps | None (uses KofamScan output) | Quick visualization of metabolic capabilities |

**Additional database**: ~60 GB (DRAM with UniRef90)
**Additional compute**: ~4-8 hours on 16 cores
**What you get**: Integrated multi-database annotations; distilled metabolic summaries; parsimony-corrected pathways; publication-ready heatmaps

### Tier 3: Specialized Analysis (Implement as Needed)

| Priority | Tool | Purpose | When to add |
|----------|------|---------|-------------|
| 7 | **antiSMASH** | Secondary metabolite BGCs | When studying antimicrobial production, siderophores |
| 8 | **METABOLIC-C** | Community biogeochemical cycling | When quantifying element cycling rates |
| 9 | **Bakta** | Improved gene prediction | When hypothetical protein rate limits analysis |
| 10 | **InterProScan** (subset) | Comprehensive domain annotation | When domain-level detail needed beyond Pfam |

### Recommended Implementation Order

```
Phase 1: Core Annotation Layer
├── KofamScan (KO assignment)
├── eggNOG-mapper v2 (COG/GO/EC/Pfam)
├── dbCAN3 (CAZymes)
├── MinPath (parsimony pathways)
└── KEGG-Decoder (visualization)

Phase 2: Integrated Profiling
├── DRAM (multi-database + distillation)
└── Per-MAG annotation tables

Phase 3: Specialized Analysis
├── antiSMASH (BGCs, per-MAG)
├── METABOLIC-C (community biogeochemistry)
└── Bakta (improved gene prediction, if needed)
```

---

## 10. Proposed Integration Architecture

### Nextflow Module Design

Following the existing pipeline architecture (see `modules/mge.nf` for reference):

```
modules/
├── annotation.nf          ← EXISTING: PROKKA_ANNOTATE
├── metabolism.nf           ← NEW: Metabolic profiling module
│   ├── KOFAMSCAN           KO assignment from Prokka .faa
│   ├── EMAPPER             eggNOG-mapper from Prokka .faa
│   ├── DBCAN               CAZyme annotation from Prokka .faa
│   ├── MERGE_ANNOTATIONS   Combine KO + COG + CAZy per protein
│   ├── MAP_TO_BINS         Filter annotations by contig2bin
│   ├── MINPATH             Parsimony pathways per MAG
│   ├── KEGG_MODULES        Module completeness scoring per MAG
│   └── KEGG_DECODER        Heatmap visualization
└── mge.nf                 ← EXISTING: MGE detection
```

### Process Definitions (Sketch)

```groovy
// modules/metabolism.nf

process KOFAMSCAN {
    conda "${projectDir}/conda-envs/dana-mag-kofamscan"
    label 'process_high'
    publishDir "${params.outdir}/metabolism/kofamscan", mode: 'copy'

    input:
    path proteins     // Prokka .faa (all assembly proteins)
    path kofam_db     // KOfam profiles directory

    output:
    path "kofamscan_results.tsv", emit: ko_assignments

    script:
    """
    exec_annotation \\
        --profile ${kofam_db}/profiles \\
        --ko-list ${kofam_db}/ko_list \\
        --cpu ${task.cpus} \\
        --format detail-tsv \\
        -o kofamscan_results.tsv \\
        ${proteins}
    """
}

process EMAPPER {
    conda "${projectDir}/conda-envs/dana-mag-emapper"
    label 'process_high'
    publishDir "${params.outdir}/metabolism/emapper", mode: 'copy'

    input:
    path proteins     // Prokka .faa
    path eggnog_db    // eggNOG database directory

    output:
    path "emapper_results.emapper.annotations", emit: annotations

    script:
    """
    emapper.py \\
        -i ${proteins} \\
        --data_dir ${eggnog_db} \\
        -m diamond \\
        --cpu ${task.cpus} \\
        --output emapper_results \\
        --override
    """
}

process DBCAN {
    conda "${projectDir}/conda-envs/dana-mag-dbcan"
    label 'process_medium'
    publishDir "${params.outdir}/metabolism/dbcan", mode: 'copy'

    input:
    path proteins     // Prokka .faa
    path dbcan_db     // dbCAN database directory

    output:
    path "overview.txt",    emit: cazyme_overview
    path "*.tsv",           emit: cazyme_details

    script:
    """
    run_dbcan \\
        ${proteins} \\
        protein \\
        --db_dir ${dbcan_db} \\
        --out_dir . \\
        --dia_cpu ${task.cpus} \\
        --hmm_cpu ${task.cpus}
    """
}

process MERGE_ANNOTATIONS {
    label 'process_low'
    publishDir "${params.outdir}/metabolism/merged", mode: 'copy'

    input:
    path ko_tsv       // KofamScan output
    path emapper_tsv  // eggNOG-mapper output
    path cazyme_tsv   // dbCAN output

    output:
    path "merged_annotations.tsv", emit: merged

    script:
    """
    merge_annotations.py \\
        --kofamscan ${ko_tsv} \\
        --emapper ${emapper_tsv} \\
        --dbcan ${cazyme_tsv} \\
        -o merged_annotations.tsv
    """
}

process MAP_TO_BINS {
    label 'process_low'
    publishDir "${params.outdir}/metabolism/per_mag", mode: 'copy'

    input:
    path merged_annotations
    path contig2bin    // DAS_Tool contig2bin.tsv
    path gff           // Prokka .gff (for gene→contig mapping)

    output:
    path "per_mag_annotations/*.tsv", emit: mag_annotations
    path "community_annotations.tsv", emit: community

    script:
    """
    map_annotations_to_bins.py \\
        --annotations ${merged_annotations} \\
        --contig2bin ${contig2bin} \\
        --gff ${gff} \\
        --outdir per_mag_annotations \\
        --community community_annotations.tsv
    """
}

process KEGG_MODULES {
    label 'process_low'
    publishDir "${params.outdir}/metabolism/modules", mode: 'copy'

    input:
    path mag_annotations  // Per-MAG annotation TSVs

    output:
    path "module_completeness.tsv", emit: modules
    path "module_heatmap.svg",      emit: heatmap

    script:
    """
    kegg_module_completeness.py \\
        --input per_mag_annotations/ \\
        --output module_completeness.tsv \\
        --heatmap module_heatmap.svg
    """
}
```

### Workflow Integration (main.nf sketch)

```groovy
// In main.nf, after PROKKA_ANNOTATE and DASTOOL_CONSENSUS

if (params.run_metabolism) {
    // Core annotation (parallel on same .faa input)
    KOFAMSCAN(
        PROKKA_ANNOTATE.out.proteins,
        params.kofam_db
    )
    EMAPPER(
        PROKKA_ANNOTATE.out.proteins,
        params.eggnog_db
    )
    DBCAN(
        PROKKA_ANNOTATE.out.proteins,
        params.dbcan_db
    )

    // Merge and map to bins
    MERGE_ANNOTATIONS(
        KOFAMSCAN.out.ko_assignments,
        EMAPPER.out.annotations,
        DBCAN.out.cazyme_overview
    )
    MAP_TO_BINS(
        MERGE_ANNOTATIONS.out.merged,
        DASTOOL_CONSENSUS.out.contig2bin,
        PROKKA_ANNOTATE.out.gff
    )

    // Pathway scoring
    KEGG_MODULES(
        MAP_TO_BINS.out.mag_annotations
    )
}
```

### Conda Environments

New environments needed:

| Environment | Tools | Estimated Size |
|-------------|-------|----------------|
| `dana-mag-kofamscan` | KofamScan, HMMER | ~200 MB |
| `dana-mag-emapper` | eggNOG-mapper, diamond | ~500 MB |
| `dana-mag-dbcan` | run_dbcan, HMMER, diamond | ~300 MB |
| `dana-mag-metabolism` | Python scripts for merge/map/scoring | ~100 MB |

```yaml
# envs/kofamscan.yml
name: dana-mag-kofamscan
channels:
  - conda-forge
  - bioconda
dependencies:
  - kofamscan
  - hmmer

# envs/emapper.yml
name: dana-mag-emapper
channels:
  - conda-forge
  - bioconda
dependencies:
  - eggnog-mapper

# envs/dbcan.yml
name: dana-mag-dbcan
channels:
  - conda-forge
  - bioconda
dependencies:
  - dbcan
  - hmmer
  - diamond
```

### Pipeline Parameters

Add to `nextflow.config`:

```groovy
params {
    // Metabolism annotation
    run_metabolism      = false
    kofam_db            = null    // Path to KOfam profiles directory
    eggnog_db           = null    // Path to eggNOG database directory
    dbcan_db            = null    // Path to dbCAN database directory

    // Optional: DRAM integration
    run_dram            = false
    dram_db             = null    // Path to DRAM database directory

    // Optional: antiSMASH per-MAG
    run_antismash       = false
    antismash_db        = null    // Path to antiSMASH database directory
}
```

### Database Management

Extend `download-databases.sh` with metabolic profiling databases:

```bash
# Add to download-databases.sh menu:
# 5) KOfam profiles
# 6) eggNOG-mapper databases
# 7) dbCAN databases
# 8) DRAM databases (UniRef90 mode)
# 9) antiSMASH databases
```

### Expected Output Structure

```
results/
├── ... (existing outputs) ...
└── metabolism/                              ← NEW
    ├── kofamscan/
    │   └── kofamscan_results.tsv            Per-protein KO assignments
    ├── emapper/
    │   └── emapper_results.emapper.annotations   Per-protein COG/GO/EC/Pfam
    ├── dbcan/
    │   ├── overview.txt                     CAZyme consensus (≥2/3 methods)
    │   ├── hmmer.out                        HMMER hits
    │   ├── diamond.out                      DIAMOND hits
    │   └── substrate.out                    Substrate predictions
    ├── merged/
    │   └── merged_annotations.tsv           Combined KO+COG+CAZy per protein
    ├── per_mag/
    │   ├── MAG_001_annotations.tsv          Per-MAG annotation tables
    │   ├── MAG_002_annotations.tsv
    │   └── ...
    ├── modules/
    │   ├── module_completeness.tsv          KEGG module completeness per MAG
    │   └── module_heatmap.svg               Visualization
    └── community/
        └── community_annotations.tsv        All annotations, all MAGs
```

---

## 11. Tool Comparison Matrix

### Annotation Coverage

| Feature | KofamScan | eggNOG-mapper | InterProScan | DRAM | METABOLIC | Mantis |
|---------|-----------|---------------|-------------|------|-----------|--------|
| KEGG KOs | **Best** (adaptive thresholds) | Good (orthology transfer) | Via mapping | Good (KOfam or UniRef) | Good (HMMER) | Good |
| COG categories | No | **Yes** | No | No | No | **Yes** |
| GO terms | No | **Yes** | **Yes** | No | No | Partial |
| EC numbers | No | **Yes** | **Yes** | Yes | Yes | **Yes** |
| Pfam domains | No | **Yes** | **Best** (source DB) | **Yes** | No | **Yes** |
| CAZymes | No | Basic | No | **Yes** (dbCAN) | No | No |
| MEROPS | No | No | No | **Yes** | No | No |
| VOGDB | No | No | No | **Yes** | No | No |
| Pathway modules | No (needs post-processing) | Yes (KEGG modules) | No | **Yes** (distillation) | **Yes** (curated) | No |
| Distillation/summary | No | No | No | **Best** | **Good** (biogeochem) | No |

### Speed and Practicality

| Feature | KofamScan | eggNOG-mapper | InterProScan | DRAM | METABOLIC |
|---------|-----------|---------------|-------------|------|-----------|
| Speed (500k proteins) | 2-4 hrs | 2-4 hrs | 8-24 hrs | 4-8 hrs | 4-12 hrs |
| Database size | 4 GB | 12-45 GB | 15-90 GB | 60-100 GB | ~5 GB |
| RAM requirement | 4-8 GB | 8-32 GB | 16-32 GB | 16-64 GB | 8-16 GB |
| Conda install | Easy | Easy | Complex | Moderate | Complex |
| KEGG license | **Free** (KOfam) | **Free** | **Free** | **Free** (UniRef mode) | **Free** |
| Maintenance | Active | Active | Active | Active | Moderate |

### Recommendation Summary

**For danaSeq, the optimal combination is:**

```
KofamScan      →  Best KO assignments (precision)
eggNOG-mapper  →  Broadest annotation (COG, GO, EC, Pfam, KO)
dbCAN3         →  CAZyme profiling (carbon cycling)
MinPath        →  Parsimony pathway reconstruction
KEGG-Decoder   →  Visualization
```

**Combined, these five tools provide:**
- KEGG KO assignments (high-precision, adaptive threshold)
- COG functional categories (community comparison)
- GO terms (standardized function ontology)
- EC numbers (enzyme classification)
- Pfam domains (protein family membership)
- CAZyme families + substrates (carbohydrate metabolism)
- KEGG Module completeness (pathway-level analysis)
- Parsimony-corrected pathway presence (biologically realistic)
- Publication-ready heatmaps

**With the addition of DRAM (Phase 2):**
- MEROPS peptidase classification
- VOGDB viral orthologous groups
- UniRef90 deep homology search
- Distilled metabolic summaries
- Electron transport chain completeness

---

## 12. ONT Long-Read Considerations

### Advantages for Metabolic Profiling

Oxford Nanopore long reads provide several advantages for metabolic annotation:

1. **Operon integrity**: Long reads preserve multi-gene operons in single contigs, enabling detection of complete biosynthetic gene clusters (antiSMASH), polysaccharide utilization loci (dbCAN3 CGC mode), and secretion system gene clusters (MacSyFinder)

2. **Reduced gene fragmentation**: Fewer interrupted coding sequences compared to short-read assemblies, leading to more complete protein predictions and better HMM/BLAST hits

3. **Longer contigs**: Flye assemblies produce longer contigs (N50 typically 50-200 kb for metagenomes), which:
   - Improve binning accuracy (better tetranucleotide signals)
   - Preserve genomic context for MGE detection (genomic islands, integrons)
   - Enable full-length BGC detection (antiSMASH clusters can span 10-100+ kb)

4. **Circular genome assembly**: Flye can produce circular contigs for complete genomes, plasmids, and phages, enabling pathway completeness assessment without assembly gaps

### Challenges

1. **Higher error rate**: ONT reads have ~5-10% error rate (R10 chemistry) which can:
   - Introduce frameshifts in Prodigal/Pyrodigal gene predictions → truncated proteins → missed annotations
   - Create false stop codons → gene fragmentation
   - Affect HMM scores (KofamScan, dbCAN) due to amino acid substitutions

2. **Mitigations already in pipeline**:
   - Flye performs self-correction during assembly (consensus polishing)
   - Prokka/Prodigal handle moderate error rates in assembled contigs
   - HMM-based tools (KofamScan, dbCAN) are more tolerant of substitutions than exact-match methods

3. **Draft MAG completeness**: 50-90% complete MAGs will have missing genes:
   - MinPath's parsimony approach handles this correctly (does not over-infer pathways)
   - KEGG module completeness scoring naturally accounts for missing genes
   - Report completeness alongside metabolic profiles for proper interpretation

---

## 13. Quality Considerations for Metabolic Profiling of MAGs

### MIMAG Standards and Metabolic Inference

The Minimum Information about a Metagenome-Assembled Genome (MIMAG) standard defines quality tiers:

| Quality | Completeness | Contamination | Metabolic Profiling Applicability |
|---------|-------------|---------------|-----------------------------------|
| High | ≥90% | <5% | Full metabolic reconstruction; pathway completeness scores are reliable |
| Medium | ≥50% | <10% | Partial metabolic profiling; pathway absence may be artifactual (assembly gap) |
| Low | <50% | <10% | Gene-level annotation only; pathway completeness unreliable |
| Contaminated | Any | ≥10% | **Do not profile** — metabolic annotations may come from contaminant organism |

### Reporting Best Practices

When reporting metabolic profiles from draft MAGs:

1. **Always report MAG quality alongside metabolic findings** — a "complete denitrification pathway" in a 55% complete MAG is less certain than in a 95% complete MAG
2. **Distinguish pathway presence vs absence**: Presence of pathway genes is informative; absence could be real or artifactual
3. **Use MinPath for conservative reconstruction** — avoids inflating pathway counts in incomplete genomes
4. **Report module completeness fractions**, not binary presence/absence — "3/4 steps of nitrate reduction present" is more informative than "nitrate reduction detected"
5. **Flag high-contamination bins** — metabolic mixing from multiple organisms creates false pathway completeness

---

## 14. Key References

### Gene Prediction

- Schwengers O et al. (2021) Bakta: rapid and standardized annotation of bacterial genomes via a comprehensive and consistent sequence database. *Microbial Genomics* 7:000685.
- Seemann T (2014) Prokka: rapid prokaryotic genome annotation. *Bioinformatics* 30:2068-2069.
- Hyatt D et al. (2010) Prodigal: prokaryotic gene recognition and translation initiation site identification. *BMC Bioinformatics* 11:119.
- Larralde M (2022) Pyrodigal: Python bindings and CLI for Prodigal. *Journal of Open Source Software* 7:4296.

### Core Annotation

- Aramaki T et al. (2020) KofamKOALA: KEGG Ortholog assignment based on profile HMM and adaptive score threshold. *Bioinformatics* 36:2251-2252.
- Cantalapiedra CP et al. (2021) eggNOG-mapper v2: Functional annotation, orthology assignments, and domain prediction at the metagenomic scale. *Molecular Biology and Evolution* 38:5825-5829.
- Paysan-Lafosse T et al. (2023) InterPro in 2022. *Nucleic Acids Research* 51:D483-D489.
- Zheng J et al. (2023) dbCAN3: automated carbohydrate-active enzyme and substrate annotation. *Nucleic Acids Research* 51:W115-W121.
- Rawlings ND et al. (2018) The MEROPS database of proteolytic enzymes, their substrates and inhibitors in 2017 and a comparison with peptidases in the PANTHER database. *Nucleic Acids Research* 46:D624-D632.
- Blin K et al. (2023) antiSMASH 7.0: new and improved predictions for detection, regulation, chemical structures and visualisation. *Nucleic Acids Research* 51:W46-W50.

### Integrated Profilers

- Shaffer M et al. (2020) DRAM for distilling microbial metabolism to automate the curation of microbiome function. *Nucleic Acids Research* 48:8883-8900.
- Zhou Z et al. (2022) METABOLIC: high-throughput profiling of microbial genomes for functional traits, metabolism, biogeochemistry, and community-scale functional networks. *Microbiome* 10:33.
- Ruiz-Perez CA et al. (2021) MicrobeAnnotator: a user-friendly, comprehensive functional annotation pipeline for microbial genomes. *BMC Bioinformatics* 22:11.
- Queirós P et al. (2021) Mantis: flexible and consensus-driven genome annotation. *GigaScience* 10:giab042.

### Pathway Reconstruction

- Ye Y & Doak TG (2009) A parsimony approach to biological pathway reconstruction/inference for genomes and metagenomes. *PLoS Computational Biology* 5:e1000465.
- Eren AM et al. (2021) Community-led, integrated, reproducible multi-omics with anvi'o. *Nature Microbiology* 6:3-6.
- Karp PD et al. (2024) The BioCyc collection of microbial genomes and metabolic pathways. *Nucleic Acids Research* 52:D649-D656.
- Beghini F et al. (2021) Integrating taxonomic, functional, and strain-level profiling of diverse microbial communities with bioBakery 3. *eLife* 10:e65088.
- Kanehisa M et al. (2016) BlastKOALA and GhostKOALA: KEGG tools for functional characterization of genome and metagenome sequences. *Journal of Molecular Biology* 428:726-731.

### Visualization

- Graham ED et al. (2018) Potential for primary productivity in a globally-distributed bacterial phototroph. *ISME Communications* 8:16.
- Sequeira JC et al. (2022) KEGGCharter: a tool for representing genomic potential in KEGG pathways. *Computational and Structural Biotechnology Journal* 20:1798-1810.
- Darzi Y et al. (2018) iPath3.0: interactive pathways explorer v3. *Nucleic Acids Research* 46:W510-W514.
- Ishigami K & Kume K (2024) ggkegg: analysis and visualization of KEGG pathway data. *Bioinformatics* 40:btae060.

### Pipeline Integration References

- nf-core/mag: Krakau S et al. (2022) nf-core/mag: a best-practice pipeline for metagenome hybrid assembly and binning. *NAR Genomics and Bioinformatics* 4:lqac007.
- ATLAS: Kieser S et al. (2020) ATLAS: a Snakemake workflow for assembly, annotation, and genomic binning of metagenome sequence data. *BMC Bioinformatics* 21:257.
- MetaWRAP: Uritskiy GV et al. (2018) MetaWRAP — a flexible pipeline for genome-resolved metagenomic data analysis. *Microbiome* 6:158.

### Quality Standards

- Bowers RM et al. (2017) Minimum information about a single amplified genome (MISAG) and a metagenome-assembled genome (MIMAG) of bacteria and archaea. *Nature Biotechnology* 35:725-731.

---

*Last updated: 2026-02-14*
