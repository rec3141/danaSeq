# Mobile Genetic Element Detection Tools

Reference guide for virus, plasmid, and mobile genetic element detection tools evaluated for integration with the danaSeq MAG assembly pipeline. Focused on **environmental/aquatic metagenomes** assembled from **Oxford Nanopore long reads** with Flye.

## Context: Why MGE Detection Matters for MAGs

A metagenome co-assembly contains three classes of contigs:

```
ASSEMBLY CONTIGS
├── Chromosomal contigs (existing pipeline: MAG binning → DAS_Tool consensus)
│   ├── Genomic islands within chromosomes
│   │   ├── Pathogenicity islands (PAIs)
│   │   ├── Resistance islands (ARGs)
│   │   ├── Metabolic islands
│   │   └── Prophages (integrated viruses)
│   ├── Integrons (gene cassette capture)
│   ├── ICEs (integrative conjugative elements)
│   └── Defense systems (CRISPR, restriction-modification)
│
├── Plasmid contigs (→ plasmid detection)
│   ├── May carry: ARGs, virulence genes, transposons
│   └── Mobility: conjugative, mobilizable, or non-mobilizable
│
└── Viral contigs (→ virus detection)
    ├── Free phages (lytic)
    └── Temperate phages (may also appear as prophages above)
```

geNomad handles the top-level classification (virus vs plasmid vs chromosome). Genomic island tools then characterize what is *inside* chromosomal contigs. CheckV assesses viral contig quality.

---

## Tool Rankings

Each category lists tools in order of preference. Tiers:

- **Recommended**: Install and integrate into the Nextflow pipeline
- **Alternative**: Good fallback; install if primary tool fails or for consensus
- **Investigate Later**: Promising but insufficient evidence or environmental validation
- **Not Recommended**: Poor fit for environmental metagenomes or unmaintained

---

## 1. Virus Detection

### Recommended

#### geNomad
- **Method**: Neural network trained on >200k marker gene profiles from IMG/VR + IMG/PR databases
- **Detects**: Viruses, plasmids, and proviruses in a single run
- **Environmental suitability**: Excellent — trained on diverse environmental data (marine, freshwater, soil, sediment), not clinical-biased
- **Accuracy**: 95.3% MCC in 2024 benchmark (best in class)
- **Min contig**: ~1 kb
- **Conda**: `bioconda::genomad`
- **Database**: ~3.5 GB download (includes marker profiles, MMseqs2 databases)
- **Input**: FASTA (contigs or assembly)
- **Output**: TSV classifications (virus/plasmid/chromosome) with scores, provirus coordinates, taxonomy
- **Maintenance**: Active (2024+)
- **Verdict**: **Top pick** — handles virus + plasmid + provirus detection in one pass. Best accuracy, best environmental generalization. Primary tool for the pipeline.
- **Citation**: Camargo et al. (2024) *Nature Biotechnology*

### Alternative

#### VirSorter2
- **Method**: Multi-classifier ensemble for different viral groups (dsDNA, ssDNA, RNA, NCLDV)
- **Environmental suitability**: Good — broad viral group coverage
- **Accuracy**: Catches ~20% of sequences geNomad misses (complementary)
- **Conda**: `bioconda::virsorter=2`
- **Input**: FASTA contigs
- **Output**: TSV with viral scores, boundary coordinates, viral group classification
- **Maintenance**: Active
- **Note**: Higher false positive rate than geNomad; slower (3-4 hrs); originally designed for Illumina but works on long-read assemblies
- **Verdict**: Good second opinion for consensus approach. Consider running alongside geNomad if viral detection is a primary research question.
- **Citation**: Guo et al. (2021) *Microbiome*

#### VIBRANT
- **Method**: Iterative annotation with neural network scoring
- **Environmental suitability**: Good
- **Conda**: `bioconda::vibrant`
- **Maintenance**: Active (Feb 2025)
- **Note**: Best functional annotation of viral metabolic genes (AMGs). Lower detection sensitivity than geNomad.
- **Verdict**: Useful for functional characterization of detected viruses, not primary detection.
- **Citation**: Kieft et al. (2020) *Microbiome*

### Investigate Later

#### DeepVirFinder
- **Method**: CNN on k-mer frequencies
- **Environmental suitability**: Moderate — trained on diverse data but limited output
- **Conda**: `bioconda::deepvirfinder`
- **Maintenance**: Active (Mar 2025)
- **Note**: Fast with GPU; strong on short fragments; black-box CNN with no plasmid detection. Complementary to geNomad.

### Not Recommended

#### viralVerify
- **Method**: Gene-content classification (part of Flye/SPAdes ecosystem)
- **Environmental suitability**: Limited
- **Conda**: Part of SPAdes package
- **Verdict**: Verification/post-filter tool only, not suitable as primary detector.

---

## 2. Plasmid Detection (Binary Classification)

### Recommended

#### geNomad (plasmid mode)
- **Method**: Neural network (same run as virus detection)
- **Environmental suitability**: Excellent — no taxon bias, works on novel environmental plasmids
- **Accuracy**: 77.8% MCC (best general-purpose in benchmarks)
- **Min contig**: ~1 kb; excellent on short fragments (<6 kb)
- **Note**: No conjugation/mobility typing (use MOBFinder post-detection if needed)
- **Verdict**: **Top pick** — already running for virus detection, plasmid output comes free.

### Alternative

#### Plasmer
- **Method**: Random forest on shared k-mers + genomic features
- **Environmental suitability**: Good (generic, not taxon-biased)
- **Min contig**: 500 bp (best performance on short contigs in Plasmer benchmark, Figure 4)
- **Conda**: `bioconda::plasmer`
- **Maintenance**: Active (2023)
- **Benchmark**: Best balanced performance across all contig lengths in direct comparison with 6 other tools (Deeplasmid, mlplasmids, PlasFlow, Platon, RFPlasmid-generic, RFPlasmid-specific)
- **Verdict**: Strong second opinion for consensus. Worth running if plasmid detection is a primary research question.

#### RFPlasmid
- **Method**: Random forest on k-mers + marker genes; 17 species-specific models + generic model
- **Environmental suitability**: Good (use generic model for environmental samples)
- **Min contig**: ~3 kb
- **Conda**: Not standard; pip install
- **Maintenance**: Active
- **Note**: Excellent when species is known (species-specific models), performance drops for unknown taxa with generic model
- **Verdict**: Useful for characterized taxa; less useful for novel environmental organisms.

#### Platon
- **Method**: Replicon distribution scores + marker gene identification
- **Environmental suitability**: Good (has `--meta` mode for metagenome contigs)
- **Min contig**: 1 kb
- **Conda**: `bioconda::platon`
- **Maintenance**: Active (2023)
- **Note**: Best characterization output among plasmid classifiers — reports replicon groups, oriT sites, ARGs, and mobilization genes alongside classification.
- **Verdict**: Good complement to geNomad for plasmid characterization.

#### PlasmidHunter
- **Method**: Gene content profile + logistic regression
- **Environmental suitability**: Both clinical and environmental
- **Min contig**: 1 kb
- **Conda**: pip + conda dependencies
- **Maintenance**: Active (2024)
- **Note**: 97.6% accuracy on 100 kb contigs; fastest tool in benchmarks. Well-suited for long-read assemblies producing large contigs.
- **Verdict**: Worth evaluating for long-read pipelines specifically.

### Investigate Later

#### PlasClass
- **Method**: Logistic regression, length-specific models
- **Environmental suitability**: Both
- **Min contig**: 1 kb
- **Conda**: Yes
- **Maintenance**: Moderate
- **Note**: Good performance, optimized for long contigs.

#### PLASMe
- **Method**: Transformer model
- **Environmental suitability**: Unknown — insufficient environmental benchmarks
- **Note**: Designed for short-read assemblies; limited nanopore validation.

#### Deeplasmid
- **Method**: CNN + gene annotations
- **Min contig**: ~300 bp
- **Conda**: Docker only
- **Maintenance**: Active (2024)
- **Note**: Validated on long reads but Docker-only installation is a barrier. Competitive at longer contigs, weaker on short ones.

#### PPR-Meta
- **Method**: BiPathCNN (3-class: phage/plasmid/chromosome)
- **Environmental suitability**: Environmental
- **Min contig**: ~500 bp
- **Note**: Interesting 3-class classification approach but requires Docker/MATLAB. Hard to install and integrate.

### Not Recommended

#### PlasFlow
- **Method**: Neural network on k-mers
- **Environmental suitability**: Environmental (in theory)
- **Min contig**: 1 kb
- **Maintenance**: **Inactive** — installation broken, unmaintained
- **Verdict**: **Avoid** — superseded by geNomad and Plasmer; broken conda package.

#### mlplasmids
- **Method**: SVM on pentamers
- **Environmental suitability**: **Clinical only** — trained on 3 species (*E. coli*, *K. pneumoniae*, *E. faecium*)
- **Maintenance**: Old (2018)
- **Conda**: R only
- **Verdict**: **Avoid** — not applicable to environmental/aquatic samples.

#### PlaScope
- **Method**: Centrifuge against species-specific database
- **Environmental suitability**: **Clinical only**
- **Maintenance**: Limited
- **Verdict**: **Avoid** — species-specific, no environmental application.

#### PlasmidFinder
- **Method**: BLAST against clinical replicon database
- **Environmental suitability**: **Very poor** — clinical database, ~50% sensitivity even on known taxa
- **Verdict**: **Avoid** — clinical isolate tool with poor environmental coverage.

#### cBar
- **Method**: Pentamer frequencies
- **Maintenance**: **Obsolete**
- **Verdict**: **Avoid** — superseded by modern ML approaches.

---

## 3. Plasmid Post-Detection & Characterization

These tools work on contigs already classified as plasmids.

### Investigate Later

#### MOBFinder
- **Purpose**: Mobility typing — classifies plasmids into 11 MOB classes
- **Environmental suitability**: Environmental (sequence-based, not database-limited)
- **Conda**: GitHub
- **Note**: Post-detection only. Classifies already-identified plasmids by conjugation/mobilization type. Useful for understanding horizontal gene transfer potential.

#### PlasMAAG
- **Purpose**: Metagenomic plasmid binning via assembly-alignment graphs
- **Environmental suitability**: Excellent (designed specifically for metagenomes)
- **Conda**: GitHub
- **Note**: Needs multi-sample assembly graphs — may not fit Flye co-assembly workflow. Interesting for future multi-sample approaches.

#### gplas2
- **Purpose**: Graph-based plasmid reconstruction/binning
- **Environmental suitability**: Both
- **Conda**: pip
- **Note**: Needs assembly graph (GFA format). Could work with Flye GFA output if available. Useful for reconstructing complete plasmid sequences.

### Not Recommended

#### MOB-suite
- **Purpose**: Replicon + relaxase + conjugation typing; plasmid reconstruction from isolates
- **Environmental suitability**: **Poor** — databases heavily biased toward Enterobacteriaceae (E. coli, Klebsiella, Salmonella)
- **Conda**: `bioconda::mob_suite`
- **Maintenance**: Active
- **Note**: Designed for isolate genomes, not metagenome co-assemblies. MOB-recon's "winner-take-all" clustering assumes all contigs come from one organism — will mis-cluster cross-species contigs on a Flye co-assembly. Database has minimal coverage of aquatic taxa (Pseudomonas, Vibrio, Alteromonas, SAR11, etc.). Must run per-MAG after DAS_Tool binning, not on raw assembly. Conjugation typing is its unique strength, but alternatives like MOBFinder can provide mobility typing without the clinical database bias.
- **Verdict**: **Not recommended** for environmental/aquatic metagenomes. Database gaps for aquatic taxa will produce low sensitivity and potentially misleading results.

#### metaplasmidSPAdes
- **Purpose**: Circular contig extraction from SPAdes assembly graphs
- **Conda**: Yes (part of SPAdes)
- **Note**: Depends on SPAdes-format assembly graphs; we use Flye, so not applicable.
- **Verdict**: **Not compatible** with Flye-based pipeline.

---

## 4. Genomic Island & Mobile Element Detection

These detect elements *within* chromosomal contigs, complementing the virus/plasmid detection layer.

### Recommended

#### IntegronFinder
- **Purpose**: Integron detection — finds attI recombination sites and gene cassettes
- **Environmental suitability**: Good — sequence-based, reference-free
- **Conda**: `bioconda::integron_finder` + pip
- **Maintenance**: Active, well-maintained
- **Input**: FASTA contigs (needs gene calls, runs Prodigal internally)
- **Output**: TSV with integron coordinates, gene cassette annotations
- **Verdict**: Fast, reliable, easy to integrate. Integrons are key drivers of ARG spread in aquatic environments.

#### IslandPath-DIMOB
- **Purpose**: Genomic island detection via dinucleotide bias + mobility gene co-occurrence
- **Environmental suitability**: Good — standalone, reference-free
- **Conda**: `bioconda::islandpath`
- **Input**: GenBank-format annotated contigs (needs gene calls from Prokka/Prodigal)
- **Output**: GFF3 with island coordinates
- **Note**: Requires annotated input (GenBank format). Our pipeline already runs Prokka, so this input is available.
- **Verdict**: Good reference-free island detector. Complements other tools.

### Alternative

#### MacSyFinder v2
- **Purpose**: Detection of secretion systems (T1SS-T9SS), conjugation systems, and other macromolecular systems
- **Environmental suitability**: Good — model-based, works on any taxa
- **Conda**: `bioconda::macsyfinder`
- **Maintenance**: Active
- **Input**: Protein FASTA (from Prodigal/Prokka)
- **Output**: TSV with system type, component genes, genomic coordinates
- **Note**: Flexible model framework; community-contributed models for various systems.
- **Verdict**: Useful for characterizing conjugation and secretion machinery in environmental organisms.

#### DefenseFinder
- **Purpose**: Anti-phage defense system detection (CRISPR arrays, restriction-modification, BREX, etc.)
- **Environmental suitability**: Good — auto-updating models cover diverse defense systems
- **Conda**: pip
- **Maintenance**: Active
- **Input**: Protein FASTA
- **Output**: TSV with defense system type, genes, coordinates
- **Verdict**: Fast, auto-updating. Understanding defense systems helps interpret phage-host dynamics.

### Investigate Later

#### ICEfinder / ICEberg 3.0
- **Purpose**: Integrative and conjugative element (ICE) detection
- **Environmental suitability**: Good
- **Conda**: Not standard
- **Input**: GenBank-format annotated contigs
- **Note**: Web + standalone modes. Useful but harder to install than the recommended tools.

#### Alien_Hunter
- **Purpose**: Horizontal gene transfer detection via interpolated variable order motifs (IVOMs)
- **Environmental suitability**: Good — reference-free, compositional approach
- **Conda**: Debian package
- **Note**: Older tool but reference-free approach is valuable for novel environmental organisms.

#### mobileOG-db
- **Purpose**: MGE protein family annotation database
- **Input**: Query via HMMER against database
- **Note**: Annotation resource, not a standalone detector. Use for characterizing detected MGEs.

### Not Recommended

#### PHASTEST
- **Purpose**: Prophage detection
- **Note**: **Web only** — cannot integrate into Nextflow pipeline. No standalone CLI.
- **Verdict**: Not usable in automated pipelines.

#### GIPSy
- **Purpose**: Lifestyle-specific genomic islands (PAI, resistance, metabolic)
- **Note**: Java GUI + CLI, manual install. No conda package.
- **Verdict**: Too difficult to automate reliably.

#### VRprofile2
- **Purpose**: Prophages + PAIs + virulence/resistance
- **Note**: Web-primary with API; difficult to pipeline reliably.
- **Verdict**: Not suitable for Nextflow integration.

---

## 5. Quality Assessment

### Recommended

#### CheckV
- **Purpose**: Viral genome completeness and contamination estimation
- **Categories**: Complete, high-quality (>90%), medium-quality (50-90%), low-quality (<50%)
- **Database**: Includes environmental viruses from IMG/VR
- **Conda**: `bioconda::checkv`
- **Input**: FASTA (viral contigs from geNomad/VirSorter2)
- **Output**: TSV with completeness, contamination, gene counts, host contamination boundaries
- **Maintenance**: Active
- **Verdict**: **Standard tool** — essential complement to geNomad. Tells you if detected viral contigs are real and how complete they are.
- **Citation**: Nayfach et al. (2021) *Nature Biotechnology*

#### CheckM2 (already in pipeline)
- **Purpose**: MAG completeness for chromosomal bins
- **Note**: **Not for plasmids or viruses** — designed for cellular genomes only. Already integrated in the existing pipeline.

### Investigate Later

#### ViralQC
- **Purpose**: Alternative viral quality assessment; 38% better contamination detection than CheckV
- **Status**: 2025 preprint — promising but too new for production use
- **Verdict**: Monitor for future adoption once peer-reviewed and stable.

### No Plasmid QA Tool Exists
- There is no equivalent of CheckV/CheckM2 for assessing plasmid completeness
- Use geNomad confidence scores + contig circularization status (DTR detection) as proxies for plasmid quality

---

## Pipeline Integration Plan

### Phase 1: Core MGE Detection (geNomad + CheckV) -- COMPLETED

Integrated in `modules/mge.nf`. geNomad detects viruses, plasmids, and proviruses from assembly contigs. CheckV assesses viral genome quality.

**Nextflow params:** `--run_genomad`, `--genomad_db`, `--run_checkv`, `--checkv_db`
**Conda environments:** `dana-mag-genomad`, `dana-mag-checkv`

### Phase 2: Genomic Island & Defense System Detection -- COMPLETED

All four tools integrated in `modules/mge.nf`, running on Prokka protein/annotation output:

| Tool | Purpose | Input | Param |
|------|---------|-------|-------|
| IntegronFinder | Integron + gene cassette detection | Assembly FASTA | `--run_integronfinder` |
| IslandPath-DIMOB | Genomic island detection (dinucleotide bias + mobility genes) | Prokka GBK | `--run_islandpath` |
| MacSyFinder v2 | Secretion systems (TXSScan) + conjugation (CONJScan) | Prokka FAA | `--macsyfinder_models` |
| DefenseFinder | Anti-phage defense systems (CRISPR, R-M, BREX, Abi, etc.) | Prokka FAA | `--run_defensefinder` |

**Conda environments:** `dana-mag-integron`, `dana-mag-islandpath`, `dana-mag-macsyfinder`, `dana-mag-defensefinder`

**Known issue:** DefenseFinder 2.0.1 bundles macsyfinder 2.1.4 which is incompatible with CasFinder model definition version 2.1. The `download-databases.sh` script automatically downgrades CasFinder 3.1.1 to 3.1.0 as a workaround (see [defense-finder#91](https://github.com/mdmparis/defense-finder/issues/91)).

### Phase 3: Consensus Plasmid Detection (future)

If plasmid detection is a primary research question, add Plasmer and/or Platon alongside geNomad for a consensus approach similar to the multi-binner strategy used for MAG binning.

---

## Installation Commands

### geNomad

```bash
# Conda environment
mamba create -p conda-envs/dana-mag-genomad -c conda-forge -c bioconda genomad

# Download database (~3.5 GB)
genomad download-database /path/to/databases
```

### CheckV

```bash
# Conda environment
mamba create -p conda-envs/dana-mag-checkv -c conda-forge -c bioconda checkv

# Download database (~1.4 GB)
checkv download_database /path/to/databases
```

### VirSorter2 (alternative)

```bash
mamba create -p conda-envs/dana-mag-virsorter2 -c conda-forge -c bioconda virsorter=2

# Download database
virsorter2 setup -d /path/to/vs2_db -j 4
```

### Plasmer (alternative)

```bash
mamba create -p conda-envs/dana-mag-plasmer -c conda-forge -c bioconda plasmer
```

### IntegronFinder (Phase 2)

```bash
mamba create -p conda-envs/dana-mag-integron -c conda-forge -c bioconda integron_finder
```

---

## Input/Output Format Notes

### geNomad Output
- `*_summary.tsv` — per-contig classification (seq_name, length, topology, coordinates, genetic_code, virus_score, plasmid_score, chromosome_score, n_hallmarks, marker_classification, taxonomy)
- `*_virus_summary.tsv` — virus-specific details (host taxonomy, provirus boundaries)
- `*_plasmid_summary.tsv` — plasmid-specific details
- `*_virus.fna` — viral contig sequences
- `*_plasmid.fna` — plasmid contig sequences
- Scores range 0-1; default threshold 0.7 for classification

### CheckV Output
- `quality_summary.tsv` — per-contig quality (contig_id, contig_length, provirus, gene_count, viral_genes, host_genes, checkv_quality [Complete/High/Medium/Low/Not-determined], miuvig_quality, completeness, completeness_method, contamination, kmer_freq, warnings)
- `viruses.fna` — trimmed viral sequences (host contamination removed)
- `proviruses.fna` — extracted provirus sequences

### DAS_Tool Integration
- geNomad contigs flagged as virus/plasmid should be **excluded** from DAS_Tool consensus binning (they are not chromosomal MAG contigs)
- Filter assembly FASTA before binning, or filter binner TSV outputs post-hoc

---

## Key References

- Camargo AP et al. (2024) Identification of mobile genetic elements with geNomad. *Nature Biotechnology* 42:1303-1312.
- Nayfach S et al. (2021) CheckV assesses the quality and completeness of metagenome-assembled viral genomes. *Nature Biotechnology* 39:578-585.
- Guo J et al. (2021) VirSorter2: a multi-classifier, expert-guided approach to detect diverse DNA and RNA viruses. *Microbiome* 9:37.
- Kieft K et al. (2020) VIBRANT: automated recovery, annotation and curation of microbial viruses, and evaluation of viral community function from genomic sequences. *Microbiome* 8:90.
- Pradier L et al. (2023) PlasX: a curated and comprehensive resource for plasmid identification from metagenomes. *mSystems* 8:e00561-23.
- Dong C et al. (2023) Plasmer: an accurate and sensitive bacterial plasmid prediction tool based on machine learning of shared k-mers and genomic features. *Microbiology Spectrum* 11:e04645-22.
- Royer G et al. (2021) PlaScope: a targeted approach to assess the plasmidome from genome assemblies at the species level. *Microbial Genomics* 4.
- Arredondo-Alonso S et al. (2018) mlplasmids: a user-friendly tool to predict plasmid- and chromosome-derived sequences for single species. *Microbial Genomics* 4.
- Robertson J & Nash JHE (2018) MOB-suite: software tools for clustering, reconstruction and typing of plasmids from draft assemblies. *Microbial Genomics* 4.

---

*Last updated: 2026-02-14*
