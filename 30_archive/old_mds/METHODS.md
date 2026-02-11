# Methods

## dānaSeq: Real-Time Metagenomic Analysis Pipeline

### Overview

dānaSeq implements a two-stage pipeline for shipboard metagenomic analysis: (1) real-time processing during sequencing and (2) post-expedition MAG assembly and characterization.

---

## Real-Time Processing

### Sample Collection and Sequencing

Water samples are collected, filtered (typically 0.22 μm pore size), and DNA extracted using standard protocols. Libraries are prepared according to Oxford Nanopore Technologies protocols and sequenced on MinION, GridION, or PromethION devices.

### Data Validation

Raw FASTQ files from MinKNOW output are validated for gzip integrity using `gzip -t`. Corrupted files are repaired using BBMap's `reformat.sh` with error correction enabled. Validated files are cached to improve resume performance on re-analysis.

### Quality Control

#### Adapter Removal
BBDuk (BBTools suite) removes sequencing adapters and PhiX control sequences:
```
bbduk.sh in=<input> out=<output> \
  ref=adapters.fa ktrim=r k=23 mink=11 hdist=1 tpe tbo
```

#### Length and Quality Filtering
Filtlong filters reads based on length (≥1000bp) and mean quality score (≥Q7):
```
filtlong --min_length 1000 --min_mean_q 7 <input> > <output>
```

### Taxonomic Classification

Kraken2 (Wood et al. 2019) performs k-mer based taxonomic assignment against custom marine-focused databases. The pipeline serializes Kraken2 calls using GNU `sem` semaphores to prevent memory exhaustion:

```bash
sem --id kraken_db_lock kraken2 \
  --db <database> --threads 4 \
  --report <report> --output <output> <input>
```

### Gene Annotation

Prokka (Seemann 2014) rapidly annotates predicted ORFs:
```
prokka --outdir <outdir> --prefix <prefix> \
  --kingdom Bacteria --metagenome <input>
```

### Functional Gene Profiling

HMMER3 (Eddy 2011) searches translated ORFs against curated HMM databases:
```
hmmsearch --cut_tc --tblout <output> \
  --domtblout <domout> <hmm_db> <protein_fasta>
```

Multiple HMM databases are supported:
- **FOAM**: Functional ontology assignments (Prestat et al. 2014)
- **CANT-HYD**: Hydrocarbon degradation genes (Khot et al. 2022)
- **NCycDB**: Nitrogen cycling genes (Tu et al. 2019)
- **HADEG**: Hydrocarbon and aromatic degradation (Rojas-Vargas et al. 2023)
- **HMDB**: Heme-binding motifs (Wang et al. 2023)
- **TASmania**: Toluene/aromatic degradation (Akarsu 2019)
- **IDOPS**: Iron-dependent oxygenases (Díaz-Valerio et al. 2021)

### Data Integration

All results are stored in DuckDB, an embedded analytical database optimized for SQL queries on large datasets. This enables real-time querying during expeditions without requiring a database server.

---

## MAG Assembly

### Read Assembly

Flye assembler (Kolmogorov et al. 2019) performs overlap-layout-consensus assembly in metagenomic mode:

```
flye --nano-raw <input> --out-dir <output> \
  --genome-size <estimate> --meta --threads <cpus> \
  --min-overlap 1000
```

### Read Mapping

minimap2 (Li 2018) aligns reads back to assembled contigs:
```
minimap2 -ax map-ont -t <threads> <assembly> <reads> | \
samtools sort -@ <threads> -o <output.bam>
```

Coverage profiles are calculated using JGI's method (Kang et al. 2019):
```
jgi_summarize_bam_contig_depths --outputDepth <output> <bam_files>
```

### Genome Binning

Three complementary binning algorithms generate initial bins:

#### SemiBin2
Deep learning-based binner trained on thousands of reference genomes (Pan et al. 2023):
```
SemiBin single_easy_bin -i <contigs> -b <bam> \
  -o <output> --sequencing-type=long_read
```

#### MetaBAT2
Tetranucleotide frequency and differential coverage-based binning (Kang et al. 2019):
```
metabat2 -i <contigs> -a <depths> -o <output> \
  -m 1500 -t <threads>
```

#### MaxBin2
Marker gene-based binning using 107 bacterial single-copy genes (Wu et al. 2016):
```
run_MaxBin.pl -contig <contigs> -abund <coverage> \
  -out <output> -thread <threads>
```

### Consensus Binning

DAS Tool (Sieber et al. 2018) integrates results from all binners, selecting the best bin for each contig based on single-copy gene analysis:

```
DAS_Tool -i <bin_list> -c <contigs> -o <output> \
  --search_engine diamond --threads <threads>
```

### Genome Refinement

Selected bins undergo iterative polishing:

1. **Racon** (2 rounds): Consensus correction using read alignments (Vaser et al. 2017)
```
racon -t <threads> <reads> <alignments> <genome> > polished.fasta
```

2. **Medaka**: Neural network-based polishing trained on Nanopore errors
```
medaka_consensus -i <reads> -d <draft> -o <output> \
  -m r941_min_high_g360
```

### Quality Assessment

CheckM2 (Chklovski et al. 2023) evaluates completeness and contamination using machine learning on a set of phylogenetically-informed marker genes:

```
checkm2 predict --threads <cpus> --input <bins> \
  --output-directory <output>
```

Bins are classified according to MIMAG standards (Bowers et al. 2017):
- **High-quality**: >90% complete, <5% contaminated, presence of 23S, 16S, 5S rRNA and ≥18 tRNAs
- **Medium-quality**: >50% complete, <10% contaminated
- **Low-quality**: <50% complete, <10% contaminated

### Taxonomic Classification

Kaiju (Menzel et al. 2016) assigns taxonomy using translated search against NCBI nr:
```
kaiju -t <nodes.dmp> -f <database> -i <genome> -o <output>
```

For high-quality MAGs, GTDB-Tk (Chaumeil et al. 2020) provides standardized taxonomy based on the Genome Taxonomy Database.

---

## Statistical Analysis

### Ordination

t-SNE (van der Maaten & Hinton 2008) reduces dimensionality of compositional profiles for visualization:
```R
library(Rtsne)
tsne_result <- Rtsne(distance_matrix, perplexity=30,
                      dims=2, check_duplicates=FALSE)
```

UMAP (McInnes et al. 2018) provides alternative ordination with preserved global structure:
```R
library(umap)
umap_result <- umap(data_matrix, n_neighbors=15, min_dist=0.1)
```

### Clustering

Graph-based clustering identifies MAG groups with similar compositional profiles:
```R
library(igraph)
graph <- graph_from_adjacency_matrix(similarity_matrix)
clusters <- cluster_louvain(graph)
```

---

## Computational Requirements

### Real-Time Processing

- **CPU**: 32+ cores recommended for parallel processing
- **RAM**: 128GB minimum (Kraken2 database loads 50-100GB)
- **Storage**: 1TB+ for expedition data
- **Network**: Not required (fully offline capable)

### MAG Assembly

- **CPU**: 64+ cores for optimal performance
- **RAM**: 256GB recommended for large assemblies
- **Storage**: 2-5TB depending on dataset size
- **Time**: Hours to days depending on coverage and diversity

---

## Software Versions

Key software versions used in development:

- **BBTools**: 38.90
- **Filtlong**: 0.2.1
- **Kraken2**: 2.1.2
- **Prokka**: 1.14.6
- **HMMER**: 3.3.2
- **Flye**: 2.9.1
- **minimap2**: 2.24
- **SemiBin2**: 1.5.0
- **MetaBAT2**: 2.15
- **MaxBin2**: 2.2.7
- **DAS Tool**: 1.1.5
- **Racon**: 1.5.0
- **Medaka**: 1.7.2
- **CheckM2**: 1.0.1
- **Kaiju**: 1.9.0
- **GTDB-Tk**: 2.1.1
- **R**: 4.2+
- **Python**: 3.9+

---

## References

Bowers RM, et al. (2017) Minimum information about a single amplified genome (MISAG) and a metagenome-assembled genome (MIMAG). *Nat Biotechnol* 35:725-731.

Chaumeil PA, et al. (2020) GTDB-Tk: a toolkit to classify genomes with the Genome Taxonomy Database. *Bioinformatics* 36:1925-1927.

Chklovski A, et al. (2023) CheckM2: a rapid, scalable and accurate tool for assessing microbial genome quality using machine learning. *Nat Methods* 20:1203-1212.

Eddy SR. (2011) Accelerated profile HMM searches. *PLoS Comput Biol* 7:e1002195.

Kang DD, et al. (2019) MetaBAT 2: an adaptive binning algorithm for robust and efficient genome reconstruction from metagenome assemblies. *PeerJ* 7:e7359.

Kolmogorov M, et al. (2019) Assembly of long, error-prone reads using repeat graphs. *Nat Biotechnol* 37:540-546.

Li H. (2018) Minimap2: pairwise alignment for nucleotide sequences. *Bioinformatics* 34:3094-3100.

McInnes L, Healy J, Melville J. (2018) UMAP: Uniform Manifold Approximation and Projection for Dimension Reduction. *arXiv* 1802.03426.

Menzel P, et al. (2016) Fast and sensitive taxonomic classification for metagenomics with Kaiju. *Nat Commun* 7:11257.

Pan S, et al. (2023) SemiBin2: self-supervised contrastive learning leads to better MAGs for short- and long-read sequencing. *Nat Commun* 14:5632.

Prestat E, et al. (2014) FOAM (Functional Ontology Assignments for Metagenomes): a Hidden Markov Model (HMM) database with environmental focus. *Nucleic Acids Res* 42:10702-10713.

Seemann T. (2014) Prokka: rapid prokaryotic genome annotation. *Bioinformatics* 30:2068-2069.

Sieber CMK, et al. (2018) Recovery of genomes from metagenomes via a dereplication, aggregation and scoring strategy. *Nat Microbiol* 3:836-843.

van der Maaten L, Hinton G. (2008) Visualizing data using t-SNE. *J Mach Learn Res* 9:2579-2605.

Vaser R, et al. (2017) Fast and accurate de novo genome assembly from long uncorrected reads. *Genome Res* 27:737-746.

Wood DE, et al. (2019) Improved metagenomic analysis with Kraken 2. *Genome Biol* 20:257.

Wu YW, et al. (2016) MaxBin 2.0: an automated binning algorithm to recover genomes from multiple metagenomic datasets. *Bioinformatics* 32:605-607.
