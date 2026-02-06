```
╔═══════════════════════════════════════════════════════════════════════════╗
║                                                                           ║
║          🧬  MAG ASSEMBLY PIPELINE  🧬                                    ║
║                                                                           ║
║              METAGENOME-ASSEMBLED GENOMES                                 ║
║                                                                           ║
║     ╔═══╗ ╔═══╗ ╔═══╗ ╔═══╗ ╔═══╗ ╔═══╗ ╔═══╗ ╔═══╗                     ║
║     ║ A ║═║ T ║═║ C ║═║ G ║═║ A ║═║ T ║═║ C ║═║ G ║                     ║
║     ╚═══╝ ╚═══╝ ╚═══╝ ╚═══╝ ╚═══╝ ╚═══╝ ╚═══╝ ╚═══╝                     ║
║     ║   ║ ║   ║ ║   ║ ║   ║ ║   ║ ║   ║ ║   ║ ║   ║                     ║
║     ╔═══╗ ╔═══╗ ╔═══╗ ╔═══╗ ╔═══╗ ╔═══╗ ╔═══╗ ╔═══╗                     ║
║     ║ T ║═║ A ║═║ G ║═║ C ║═║ T ║═║ A ║═║ G ║═║ C ║                     ║
║     ╚═══╝ ╚═══╝ ╚═══╝ ╚═══╝ ╚═══╝ ╚═══╝ ╚═══╝ ╚═══╝                     ║
║                                                                           ║
║        Reconstruct complete microbial genomes from metagenomic soup!     ║
║                                                                           ║
╚═══════════════════════════════════════════════════════════════════════════╝
```

## 🎯 The Grand Challenge

You've got **millions of DNA reads** from **hundreds of species** all mixed together. 🌪️

Your mission: **Separate them into individual genomes**. 🎯

```
    Mixed Soup of Reads              Individual Genomes

    🧬🧬🧬🧬🧬🧬🧬              →      Genome A: 🧬🧬🧬 (Prochlorococcus)
    🧬🧬🧬🧬🧬🧬🧬              →      Genome B: 🧬🧬🧬 (Pelagibacter)
    🧬🧬🧬🧬🧬🧬🧬              →      Genome C: 🧬🧬🧬 (Synechococcus)
    🧬🧬🧬🧬🧬🧬🧬              →      Genome D: 🧬🧬🧬 (SAR86)
    🧬🧬🧬🧬🧬🧬🧬              →      ...and more!

    CHAOS                             ORDER
```

---

## 🌊 The Assembly Flow

```
┌─────────────────────────────────────────────────────────────────────┐
│                                                                     │
│  📦 Reads → 🏗️ Assemble → 📍 Map → 🗂️ Bin → ✨ Polish → 🏆 MAGs   │
│                                                                     │
│  Short     Long         Coverage  Separate  Refine    High-        │
│  pieces    contigs      profiles  genomes   quality   quality      │
│                                                                     │
└─────────────────────────────────────────────────────────────────────┘
```

---

## 📚 The Complete MAG Assembly Cookbook

### 🏗️ **10s: Assembly** - Build longer contigs!

```
╔══════════════════════════════════════════════════════════════╗
║                                                              ║
║   🧩 READS  +  🧩 READS  +  🧩 READS                         ║
║         ↓          ↓          ↓                              ║
║   ═══════════════════════════════════                        ║
║         FLYE ASSEMBLER (overlap-based)                       ║
║   ═══════════════════════════════════                        ║
║         ↓                                                    ║
║   🏗️ LONG CONTIGS (100kb - 10Mb!)                           ║
║                                                              ║
╚══════════════════════════════════════════════════════════════╝
```

**Assembly Scripts:**

| Script | Strategy | Use Case |
|--------|----------|----------|
| `10_assembly_flye.sh` | 🌍 Co-assembly (all samples) | **RECOMMENDED** - Maximum contiguity |
| `11_assembly_flye29.sh` | 🔧 Alternative parameters | Different overlap settings |
| `12_assembly_flye_each.sh` | 📦 Per-sample assembly | Strain variation analysis |
| `13_assembly_flye_bins.sh` | 🎯 Individual bin assembly | Polish specific MAGs |

**The Assembly Process:**
```
1. 🧹 Deduplicate reads (BBMap dedupe)
2. 🔍 Quality filter (Filtlong: Q7+, length >1kb)
3. 🏗️ Metagenomic assembly (Flye --meta, --min-overlap 1000)
4. 📊 Generate contigs (typically 1,000-50,000 contigs)
```

**Co-assembly vs Per-Sample:**
```
Co-Assembly (RECOMMENDED)              Per-Sample Assembly
├─ ✅ Better for low-abundance          ├─ ✅ Better for high diversity
├─ ✅ Longer contigs                    ├─ ✅ Preserves strain variants
├─ ✅ More complete genomes             ├─ ✅ Less chimeric
└─ ⚠️ Can create chimeras               └─ ⚠️ Shorter, more fragmented
```

---

### 📍 **20s: Read Mapping** - Who's covering what?

```
╔══════════════════════════════════════════════════════════════╗
║                                                              ║
║   Assembly: ══════════════════════════════════════           ║
║                                                              ║
║   Sample 1:   🧬🧬🧬        🧬🧬     🧬🧬🧬🧬                ║
║   Sample 2:      🧬🧬🧬🧬🧬    🧬     🧬🧬                  ║
║   Sample 3:   🧬    🧬🧬🧬🧬🧬🧬           🧬🧬             ║
║                                                              ║
║   Coverage:   ████████████████████████████████               ║
║   Profile:    High    Medium    High    Low                 ║
║                                                              ║
╚══════════════════════════════════════════════════════════════╝
```

**Mapping Scripts:**
- `20_mapping.sh` - 🎯 Map all samples to assembly (minimap2)
- `21_mapping_bins.sh` - 🔄 Remap to refined bins
- `22_calculate_coverage.sh` - 📊 Depth calculation (JGI method)
- `23_coverage.sh` - 📈 Alternative coverage metrics

**Why Map?**
- 📊 **Coverage profiles** - Who's abundant where?
- 🗂️ **Binning signal** - Contigs with similar coverage patterns = same genome!
- ✅ **Quality metrics** - Even coverage = good bin

---

### 🗂️ **30s: Binning** - Separate the genomes!

```
╔══════════════════════════════════════════════════════════════════════╗
║                                                                      ║
║                  🧬 BINNING - THE MAGIC STEP 🧬                      ║
║                                                                      ║
║   We group contigs into MAGs using TWO signals:                      ║
║                                                                      ║
║   1️⃣ COMPOSITION (Tetranucleotide Frequency)                        ║
║      Each species has a unique "DNA fingerprint"                    ║
║      ATCG GGTA CCAT TTAG ... (256 possible tetramers)               ║
║                                                                      ║
║   2️⃣ COVERAGE (Abundance across samples)                            ║
║      Species with similar abundance profiles = likely same organism ║
║                                                                      ║
║   ┌─────────────┐     ┌─────────────┐     ┌─────────────┐         ║
║   │  SemiBin2   │     │  MetaBAT2   │     │  MaxBin2    │         ║
║   │ (Deep Learn)│     │  (TNF+COV)  │     │ (Markers)   │         ║
║   └──────┬──────┘     └──────┬──────┘     └──────┬──────┘         ║
║          │                   │                    │                 ║
║          └───────────────────┴────────────────────┘                 ║
║                              │                                      ║
║                    ┌─────────▼──────────┐                          ║
║                    │     DAS TOOL        │                          ║
║                    │  (Consensus Best)   │                          ║
║                    └────────────────────┘                          ║
║                              │                                      ║
║                         🏆 BEST MAGS! 🏆                            ║
║                                                                      ║
╚══════════════════════════════════════════════════════════════════════╝
```

**The Binning Trinity:**

```
30_binning_semibin.sh    →  🤖 Deep learning approach
                            ├─ Trained on thousands of genomes
                            ├─ Best for complex communities
                            └─ Long-read optimized mode

31_binning_metabat.sh    →  📊 Classic TNF + coverage
                            ├─ Tetranucleotide frequency
                            ├─ Differential coverage
                            └─ Fast and reliable

32_binning_maxbin.sh     →  🧬 Marker gene approach
                            ├─ 107 bacterial marker genes
                            ├─ Phylogenetic signal
                            └─ Good for known taxa
```

**Why Use All Three?**
```
Each tool has strengths and weaknesses:
├─ SemiBin2: Great overall, but needs training
├─ MetaBAT2: Robust, but can split genomes
└─ MaxBin2: Precise, but misses novel taxa

DAS Tool picks the BEST bins from each! 🏆
```

---

### ✨ **40s: Polishing** - Make it shine!

```
╔═══════════════════════════════════════════════════════════════╗
║                                                               ║
║   Draft Genome:  ATCGATCGXXTCGATCGA  ← Errors! 😱            ║
║                         ↓                                     ║
║           ╔═══════════════════════╗                          ║
║           ║  RACON (2 rounds)     ║  ← Consensus correction  ║
║           ╚═══════════════════════╝                          ║
║                         ↓                                     ║
║   Better Genome: ATCGATCGATCGATCGA                           ║
║                         ↓                                     ║
║           ╔═══════════════════════╗                          ║
║           ║  MEDAKA               ║  ← Neural network polish ║
║           ╚═══════════════════════╝                          ║
║                         ↓                                     ║
║   Perfect!       ATCGATCGATCGATCGA  ← Publication ready! 🏆  ║
║                                                               ║
╚═══════════════════════════════════════════════════════════════╝
```

- `40_polish_assemblies.sh` - 💎 Racon (2x) + Medaka pipeline
- `41_medaka_split.sh` - ⚡ Parallelized for speed

**Polishing = Accuracy:**
- 🎯 Fixes sequencing errors
- 📈 Improves gene calling
- ✅ Essential for publication
- 🔬 Reduces false SNPs

---

### 🏷️ **50s: Bin Characterization** - Who are they?

```
╔═══════════════════════════════════════════════════════════════╗
║                                                               ║
║   Mystery Genome                                              ║
║   ============                                                ║
║   ATCGATCGATCGATCG...                                         ║
║          ↓                                                    ║
║   ╔═════════════════╗                                         ║
║   ║  KRAKEN2        ║  → "Prochlorococcus marinus"           ║
║   ╚═════════════════╝                                         ║
║          ↓                                                    ║
║   ╔═════════════════╗                                         ║
║   ║  SENDSKETCH     ║  → "98.5% match to MIT9313"            ║
║   ╚═════════════════╝                                         ║
║          ↓                                                    ║
║   ╔═════════════════╗                                         ║
║   ║  CHECKM2        ║  → "95% complete, 2% contamination"    ║
║   ╚═════════════════╝                                         ║
║                                                               ║
║   Result: HIGH-QUALITY PROCHLOROCOCCUS MAG! 🏆               ║
║                                                               ║
╚═══════════════════════════════════════════════════════════════╝
```

- `50_run_kraken_all.sh` - 🏷️ Taxonomic classification
- `51_sketch_bins.sh` - 🎯 Fast species-level ID
- `52_tetra_frequency.sh` - 🧮 Composition profiles

---

### 🎼 **60s: Complete Pipelines** - One-stop shop!

```
╔══════════════════════════════════════════════════════════════════╗
║                                                                  ║
║             🎼 THE FULL ORCHESTRA 🎼                             ║
║                                                                  ║
║   60_map_and_bin_complete.sh                                     ║
║   61_map_and_bin_optimized.sh  ← 🏆 RECOMMENDED                  ║
║                                                                  ║
║   Does EVERYTHING:                                               ║
║   ├─ 🏗️  Assemble (if needed)                                    ║
║   ├─ 📍 Map all samples                                          ║
║   ├─ 📊 Calculate coverage                                       ║
║   ├─ 🗂️  Bin with 3 tools                                        ║
║   ├─ 🤝 DAS Tool consensus                                       ║
║   ├─ ✅ CheckM2 quality                                          ║
║   ├─ ✨ Polish bins                                              ║
║   └─ 🏷️  Classify taxonomy                                       ║
║                                                                  ║
║   Input:  Raw reads                                              ║
║   Output: Publication-ready MAGs! 🎉                             ║
║                                                                  ║
╚══════════════════════════════════════════════════════════════════╝
```

**When to use:**
```bash
# Just starting? Use this!
./61_map_and_bin_optimized.sh

# It handles everything from reads → MAGs
# ✅ Smart checkpointing (resume if crashed)
# ✅ Quality filtering at each step
# ✅ Optimized parameters
# ✅ Detailed logging
```

---

### 🔄 **70s: Format Conversion** - Play nice with others

```
    ╔════════════════╗
    ║  Kraken Output ║
    ╚═══════╤════════╝
            │
     71_kraken2_to_bandage.sh
            │
            ▼
    ╔════════════════╗       ╔══════════════════╗
    ║  GFA Format    ║  →    ║   BANDAGE        ║
    ╚════════════════╝       ║   Visualization  ║
                             ╚══════════════════╝

                             Beautiful assembly graphs! 📊
```

- `70_kraken2_to_bandage.sh` - 🔄 Convert for Bandage visualization
- `71_circular_gfa.awk` - ⭕ Handle circular contigs (plasmids!)

---

### 📊 **80s: Visualization** - Make it beautiful!

```
╔════════════════════════════════════════════════════════════════════╗
║                                                                    ║
║                    📊 VISUALIZATION SUITE 📊                       ║
║                                                                    ║
║   80_plot_bins.R                                                   ║
║   ├─ 🎨 PCA of MAG composition                                     ║
║   ├─ 🗺️  t-SNE clustering                                          ║
║   ├─ 🔥 Heatmaps of abundance                                      ║
║   └─ 📈 Quality score distributions                                ║
║                                                                    ║
║   81_plot_mapping.R                                                ║
║   ├─ 📊 Coverage plots per sample                                  ║
║   ├─ 📉 Mapping statistics                                         ║
║   └─ 🎯 Read depth distributions                                   ║
║                                                                    ║
║   82_inter_binning_analysis.r  ← 🌟 ADVANCED                       ║
║   ├─ 🌀 t-SNE in 2D/3D                                             ║
║   ├─ 🗺️  UMAP projections                                          ║
║   ├─ 🕸️  Graph-based clustering                                    ║
║   └─ 🎨 Interactive visualizations                                 ║
║                                                                    ║
╚════════════════════════════════════════════════════════════════════╝
```

**Example outputs:**

```
    PCA of MAG Composition              Abundance Heatmap

      PC2
       │                                Sample1  Sample2  Sample3
       │   🔴🔴                          MAG_01  ████     ░░░░     ██
       │     🔵🔵🔵                       MAG_02  ██       ████     ██
    ───┼────────── PC1                  MAG_03  ░░       ██       ████
       │ 🟢🟢                            MAG_04  ████     ████     ░░░░
       │    🟡                           MAG_05  ██       ░░       ██
       │

    Species cluster beautifully!        Abundance varies across samples!
```

---

### 🌐 **90s: Integration** - Connect to bigger picture

```
╔══════════════════════════════════════════════════════════════╗
║                                                              ║
║   Your MAGs  →  90_edna_schema.r  →  Database                ║
║                                                              ║
║              →  91_foam_ecoserv.r →  Ecosystem Services     ║
║                                                              ║
║   Link MAGs to:                                              ║
║   ├─ 🗺️  Geographic locations                               ║
║   ├─ 🌡️  Environmental parameters                           ║
║   ├─ 📅 Time series data                                     ║
║   ├─ 🧬 Metabolic pathways                                   ║
║   └─ 🌍 Ecosystem function predictions                       ║
║                                                              ║
╚══════════════════════════════════════════════════════════════╝
```

- `90_edna_schema.r` - 💾 Database schema for eDNA integration
- `91_foam_ecoserv.r` - 🌍 FOAM ecosystem services analysis

---

## 🏆 What Makes a Good MAG?

```
╔═══════════════════════════════════════════════════════════════╗
║                                                               ║
║   MAG QUALITY STANDARDS (MIMAG - Minimum Information)        ║
║                                                               ║
║   🥇 HIGH QUALITY                                             ║
║      ├─ >90% completeness                                    ║
║      ├─ <5% contamination                                    ║
║      └─ Presence of 23S, 16S, 5S rRNA & tRNAs               ║
║                                                               ║
║   🥈 MEDIUM QUALITY                                           ║
║      ├─ >50% completeness                                    ║
║      └─ <10% contamination                                   ║
║                                                               ║
║   🥉 LOW QUALITY (still useful!)                             ║
║      ├─ <50% completeness                                    ║
║      └─ <10% contamination                                   ║
║                                                               ║
║   ❌ POOR QUALITY (discard)                                   ║
║      └─ >10% contamination                                   ║
║                                                               ║
╚═══════════════════════════════════════════════════════════════╝
```

**CheckM2 evaluates:**
- ✅ Presence of single-copy marker genes
- ❌ Duplicated markers (= contamination)
- 📊 Overall genome completeness
- 🧬 Taxonomic consistency

---

## 🚀 Quick Start

### Complete Pipeline (Easiest!)

```bash
cd 20_mag_assembly

# 🎼 Run the full orchestra
./61_map_and_bin_optimized.sh

# Go grab coffee ☕ (this takes hours-days depending on data size)

# 📊 Visualize when done
Rscript 80_plot_bins.R
Rscript 82_inter_binning_analysis.r
```

### Step-by-Step (More Control)

```bash
# 1. 🏗️ Assemble
./10_assembly_flye.sh

# 2. 📍 Map
./20_mapping.sh

# 3. 📊 Coverage
./22_calculate_coverage.sh

# 4. 🗂️ Bin (run all three!)
./30_binning_semibin.sh
./31_binning_metabat.sh
./32_binning_maxbin.sh

# DAS Tool consensus happens automatically in binning scripts

# 5. ✨ Polish
./40_polish_assemblies.sh

# 6. 🏷️ Characterize
./50_run_kraken_all.sh
./51_sketch_bins.sh

# 7. 📊 Visualize
Rscript 80_plot_bins.R
```

---

## 💡 Pro Tips

### 🎯 Co-Assembly Best Practices

```
✅ DO co-assemble when:
   ├─ Same environment type across samples
   ├─ Looking for core community members
   └─ Want maximum contiguity

❌ DON'T co-assemble when:
   ├─ Vastly different environments
   ├─ Studying strain-level variation
   └─ Samples separated by time/space
```

### 🗂️ Binning Optimization

```
For BEST results:
├─ 📊 Use samples with differential abundance
│   (don't bin samples that are too similar!)
├─ 🎯 Aim for 10-50x coverage
│   (too low = poor binning, too high = wasted $)
├─ 🧬 More samples = better signal
│   (5+ samples ideal)
└─ 🔍 Remove obvious contaminants first
    (human, plant, PhiX)
```

### ✨ Polishing Wisdom

```
Racon + Medaka is the gold standard for Nanopore:
├─ Racon: Fast, improves ~1-2% errors → 0.5%
├─ Medaka: Slower, neural net trained on Nanopore
└─ Result: ~99.9%+ accuracy

For publication:
✅ Always polish MAGs
✅ Report polishing methods
✅ Check gene calling improves
```

---

## 📊 Expected Output

### Directory Structure

```
mag_assembly_output/
├── 📁 assembly/
│   ├── assembly.fasta              ← Your contigs!
│   └── assembly_graph.gfa
│
├── 📁 mapping/
│   ├── sample1.sorted.bam
│   ├── sample2.sorted.bam
│   └── coverage_table.txt
│
├── 📁 bins/
│   ├── 📁 semibin/
│   │   ├── bin.1.fa
│   │   ├── bin.2.fa
│   │   └── ...
│   ├── 📁 metabat/
│   │   └── ...
│   ├── 📁 maxbin/
│   │   └── ...
│   └── 📁 das_tool_consensus/      ← 🏆 USE THESE!
│       ├── MAG_00001.fa
│       ├── MAG_00002.fa
│       └── ...
│
├── 📁 polished/
│   ├── MAG_00001.polished.fa
│   └── ...
│
├── 📁 checkm2/
│   └── quality_report.tsv          ← Completeness/contamination
│
└── 📁 taxonomy/
    ├── kraken_results.txt
    └── taxonomy_assignments.txt
```

### Statistics to Report

```
Our expedition recovered:
├─ 127 high-quality MAGs (>90% complete, <5% contamination)
├─ 213 medium-quality MAGs (>50% complete, <10% contamination)
├─ Representing 15 phyla
├─ Including 34 novel species (no close relatives)
└─ Total: 450 Mbp of novel genomic content

Dominant taxa:
1. 🦠 Alphaproteobacteria (45 MAGs) - SAR11, Rhodobacterales
2. 🦠 Cyanobacteria (23 MAGs) - Prochlorococcus, Synechococcus
3. 🦠 Bacteroidota (18 MAGs) - Flavobacteriaceae
4. 🦠 Verrucomicrobiota (12 MAGs) - Puniceicoccaceae
```

---

## 🐛 Troubleshooting

### Problem: Assembly too fragmented (N50 < 10kb)
```
🔍 Possible causes:
   ├─ Low coverage (<10x)
   ├─ High diversity (many strains)
   └─ Poor read quality

✅ Solutions:
   ├─ Sequence deeper
   ├─ Try per-sample assembly
   └─ Filter reads more strictly (Q10+)
```

### Problem: Bins have high contamination (>10%)
```
🔍 Possible causes:
   ├─ Insufficient coverage variation
   ├─ Closely related strains
   └─ Chimeric contigs

✅ Solutions:
   ├─ Add more samples with different abundances
   ├─ Use strain-resolution binning
   └─ Manual refinement in tool like Anvi'o
```

### Problem: Few bins recovered
```
🔍 Possible causes:
   ├─ Community too diverse (>1000 species)
   ├─ Coverage too low
   └─ All samples too similar

✅ Solutions:
   ├─ Focus on abundant taxa (>1% abundance)
   ├─ Sequence deeper (aim for 10Gb+ per sample)
   └─ Collect samples across gradients
```

---

## 🎓 Further Reading

### Essential Papers
- 📄 **Flye**: Kolmogorov et al., Nature Biotechnology 2019
- 📄 **SemiBin2**: Pan et al., Nature Communications 2023
- 📄 **MetaBAT2**: Kang et al., PeerJ 2019
- 📄 **DAS Tool**: Sieber et al., Nature Microbiology 2018
- 📄 **CheckM2**: Chklovski et al., Nature Methods 2023

### Standards
- 📖 **MIMAG**: Minimum Information about MAGs (Bowers et al., 2017)
- 📖 **GTDB**: Genome Taxonomy Database

### Tools Documentation
- 🔗 Flye: [github.com/fenderglass/Flye](https://github.com/fenderglass/Flye)
- 🔗 minimap2: [github.com/lh3/minimap2](https://github.com/lh3/minimap2)
- 🔗 SemiBin: [github.com/BigDataBiology/SemiBin](https://github.com/BigDataBiology/SemiBin)
- 🔗 CheckM2: [github.com/chklovski/CheckM2](https://github.com/chklovski/CheckM2)

---

```
╔═══════════════════════════════════════════════════════════════════╗
║                                                                   ║
║           🧬 FROM CHAOS TO GENOMES 🧬                             ║
║                                                                   ║
║        Assemble • Map • Bin • Polish • Characterize              ║
║                                                                   ║
║              Now go build some MAGs! 🏗️                          ║
║                                                                   ║
╚═══════════════════════════════════════════════════════════════════╝
```
