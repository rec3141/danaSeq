# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## ⚠️ CRITICAL SAFETY RULES ⚠️

**File Deletion Policy:**
1. **ALWAYS use `trash` or `trash-put` instead of `rm`** - trash-cli is installed on this system
2. Files can be recovered from trash with `trash-restore` or `trash-list`
3. NEVER use `rm -rf` without explicit user confirmation
4. Before any deletion, create a manifest: `ls -laR path/ > manifest_$(date +%Y%m%d_%H%M%S).txt`
5. "Clean up" means organize/tidy, NOT delete - ask for clarification
6. "Skip" means ignore/don't commit, NOT delete - ask for clarification

**Examples:**
```bash
# Good - recoverable
trash old_file.txt
trash-put deprecated_directory/

# Bad - permanent
rm -rf directory/  # DON'T DO THIS
```

## Project Overview

**Dana Pipeline** - A real-time metagenomic analysis system for Oxford Nanopore sequencing on oceanographic expeditions. Processes DNA reads as they stream from the sequencer, providing live taxonomic classification, gene annotation, MAG (Metagenome-Assembled Genome) assembly, and interactive geographic visualization.

This is a production bioinformatics pipeline used on research vessels and field stations for marine eDNA analysis.

## Repository Structure

```
danav2/
├── 10_realtime_processing/   # Real-time analysis at sea (26 scripts)
├── 20_mag_assembly/           # MAG reconstruction pipeline (26 scripts)
├── 30_archive/                # Deprecated/experimental code
├── agents/                    # Expert advisor scripts
├── python_pipeline/           # Modern Python rewrite (in progress)
├── banner.sh                  # Welcome screen utility
├── status.sh                  # Dependency checker
└── agents.sh                  # Agent launcher
```

### Numbered Directory System

Scripts use a **numbered prefix with intentional gaps** (10, 20, 30...) allowing future expansion:
- **10s-20s:** Core processing steps
- **30s-40s:** Secondary analysis
- **50s-60s:** Integration & complete pipelines
- **70s+:** Visualization & reporting

Example: `24_process_reads_optimized.sh` indicates step 24 in the processing workflow.

## Common Commands

### Real-Time Processing (Shipboard/Field)

```bash
cd 10_realtime_processing

# RECOMMENDED: Full pipeline with AI-enhanced processing
./24_process_reads_optimized.sh -i /path/to/nanopore/data -K -P

# Fast mode for urgent analysis (screening, initial assessment)
./22_process_reads_fast.sh -i /path/to/barcode_dir -P -S

# With HMM search for functional genes
./24_process_reads_optimized.sh -i /path/to/data -P --hmm /path/to/CANT-HYD.hmm

# Force re-run optional stages
./24_process_reads_optimized.sh -i /path/to/data -P --hmm HMM.hmm --force

# Launch interactive dashboard (separate terminal)
Rscript 60_edna_mapping_viz.r
```

**Critical Flags:**
- `-K`: Kraken2 taxonomic classification (ONLY use with `24_process_reads_optimized.sh`)
- `-P`: Prokka gene annotation
- `-S`: Sendsketch taxonomic profiling
- `-T`: Tetranucleotide frequency calculation
- `--hmm`: HMM search on predicted proteins
- `--force`: Re-run optional stages even if output exists

### MAG Assembly (Post-Expedition)

```bash
cd 20_mag_assembly

# RECOMMENDED: Complete pipeline (assembly → binning → polishing)
./61_map_and_bin_optimized.sh

# Step-by-step approach (for more control):
./10_assembly_flye.sh              # Assemble contigs
./20_mapping.sh                    # Map reads to assembly
./22_calculate_coverage.sh         # Calculate coverage
./30_binning_semibin.sh            # Bin with SemiBin2
./31_binning_metabat.sh            # Bin with MetaBAT2
./32_binning_maxbin.sh             # Bin with MaxBin2
# DAS Tool consensus runs automatically in binning scripts
./40_polish_assemblies.sh          # Polish with Racon + Medaka
./50_run_kraken_all.sh             # Classify bins taxonomically

# Visualize results
Rscript 80_plot_bins.R
Rscript 82_inter_binning_analysis.r
```

### Utilities

```bash
# Check installed dependencies and versions
./status.sh

# Display welcome banner
./banner.sh

# Consult expert advisors
./agents.sh
```

### Python Pipeline (Modern Rewrite - In Progress)

```bash
cd python_pipeline

# Install in development mode
pip install -e .

# With visualization extras
pip install -e ".[viz]"

# With development tools
pip install -e ".[dev]"

# Run tests (when implemented)
pytest tests/
```

## Architecture & Key Concepts

### 1. Two-Pipeline Design

**10_realtime_processing/** - Streaming analysis during sequencing
- Input: Raw FASTQ from MinKNOW sequencer
- Output: QC'd reads, taxonomic classifications, gene annotations, DuckDB database
- Use case: Shipboard monitoring, HAB detection, real-time decisions

**20_mag_assembly/** - Genome reconstruction post-expedition
- Input: QC'd reads from realtime processing
- Output: High-quality MAGs (bins), polished genomes, taxonomic assignments
- Use case: Publication-quality genomes, comparative genomics

### 2. Critical Kraken2 Memory Issue

**IMPORTANT:** Kraken2 loads 50-100GB database into RAM. The `24_process_reads_optimized.sh` script uses semaphores to serialize Kraken2 calls while keeping other steps parallel. Other scripts will spawn multiple Kraken2 instances and crash the system.

**Rule:** When using `-K` flag, ONLY use `24_process_reads_optimized.sh`.

See `10_realtime_processing/CRITICAL_KRAKEN_BUG.md` for details.

### 3. Parallelization Strategy (32×1 Model)

Uses GNU parallel with 32 workers processing 1 file each:
- Most steps: Fully parallel (BBDuk, Filtlong, FASTA conversion)
- Kraken2: Serialized with `sem --id kraken_db_lock`
- DuckDB writes: Serialized with `sem --id duckdb_lock`

### 4. Resume Logic & Stage-Aware Processing

All scripts have intelligent resume capability:
- Basic QC: Skips if `.fa` file exists
- Optional stages (Prokka, HMM, Kraken): Each checks independently
- Atomic operations: Uses `.tmp` and `.partial` files to prevent corruption on interrupt

Example workflow:
```bash
# Initial run with Prokka
./24_process_reads_optimized.sh -i data -P

# Later add HMM - only HMM runs, rest skipped
./24_process_reads_optimized.sh -i data -P --hmm CANT-HYD.hmm

# Force re-run of all optional stages
./24_process_reads_optimized.sh -i data -P --hmm CANT-HYD.hmm --force
```

See `10_realtime_processing/RESUME_LOGIC.md` for detailed scenarios.

### 5. MAG Assembly: The Binning Trinity

Three binning tools provide consensus for best results:
- **SemiBin2** (deep learning) - Best for complex communities
- **MetaBAT2** (TNF + coverage) - Fast, reliable classic approach
- **MaxBin2** (marker genes) - Good for known taxa

**DAS Tool** combines all three and selects the best bins. Always use all three binners for publication-quality results.

### 6. FASTQ Validation & Caching

Nanopore sequencers can produce corrupted gzip files. The pipeline:
1. Validates with `gzip -t`
2. Repairs with BBMap `reformat.sh` if needed
3. Caches validated files in `/data/.fastq_pass` (configurable via `CACHE_FASTQ`)

### 7. DuckDB Integration

R scripts (`4X_*_db.r`) integrate results into DuckDB for real-time SQL queries:
- `40_kraken_db.r` - Kraken classifications
- `42_prokka_db.r` - Prokka annotations
- `43_sketch_db.r` - Sendsketch profiles
- `44_tetra_db.r` - Tetranucleotide frequencies

DuckDB enables SQL queries on growing datasets during expeditions without a database server.

## Expected Input/Output Formats

### Input Structure (Oxford Nanopore)

```
input_dir/
└── fastq_pass/
    ├── barcode01/
    │   └── FLOWCELL_pass_barcode01_xxx.fastq.gz
    ├── barcode02/
    │   └── FLOWCELL_pass_barcode02_xxx.fastq.gz
    └── ...
```

### Output Structure

```
output/
├── FLOWCELL/
│   ├── barcode01/
│   │   ├── fa/           # Final QC'd FASTA sequences
│   │   ├── fq/           # Intermediate FASTQ (BBDuk)
│   │   ├── kraken/       # Kraken2 classifications (.tsv, .report)
│   │   ├── prokka/       # Prokka annotations (per-sample dirs)
│   │   ├── hmm/          # HMM search results (.tsv, .tbl)
│   │   ├── sketch/       # Sendsketch profiles
│   │   ├── tetra/        # Tetranucleotide frequencies
│   │   └── log.txt       # Processing log
│   └── barcode02/
│       └── ...
└── failed_files.txt      # Failure diagnostics
```

### MAG Assembly Output

```
mag_assembly_output/
├── assembly/
│   └── assembly.fasta              # Assembled contigs
├── mapping/
│   ├── sample*.sorted.bam          # Read alignments
│   └── coverage_table.txt          # Coverage per contig
├── bins/
│   ├── semibin/
│   ├── metabat/
│   ├── maxbin/
│   └── das_tool_consensus/         # ← USE THESE! Best consensus bins
│       ├── MAG_00001.fa
│       └── ...
├── polished/
│   └── MAG_*.polished.fa           # Racon + Medaka polished
├── checkm2/
│   └── quality_report.tsv          # Completeness/contamination
└── taxonomy/
    └── taxonomy_assignments.txt     # Taxonomic classifications
```

## Environment & Dependencies

### Required Tools

**Real-Time Processing:**
- BBMap (`BBMAP=/work/apps/bbmap`)
- Filtlong (`FILTLONG=/work/apps/Filtlong/bin/filtlong`)
- GNU parallel (must be GNU version, not moreutils)

**Optional (flag-dependent):**
- Kraken2 (`KRAKEN2=/usr/bin/kraken2`, `KRAKEN_DB=/path/to/db`) - for `-K`
- Prokka (`PROKKA_BIN=/work/apps/prokka/bin/prokka`) - for `-P`
- HMMER (`HMMSEARCH=/usr/bin/hmmsearch`) - for `--hmm`
- Perl with tetramer scripts (`APPS=/work/apps`) - for `-T`

**MAG Assembly:**
- Flye (metagenomic assembler)
- minimap2 (read mapping)
- SemiBin2, MetaBAT2, MaxBin2 (binning)
- DAS Tool (consensus binning)
- Racon, Medaka (polishing)
- CheckM2 (quality assessment)
- Kaiju (taxonomic classification)

**Analysis & Visualization:**
- R with tidyverse, DuckDB, leaflet packages
- Python 3.9+ (for Python pipeline)

### Hardcoded Paths

Scripts contain hardcoded paths for active expeditions:
- CMO2025: California to Mexico Oceanographic Survey
- QEI2025: Queen Elizabeth Islands Arctic Expedition

**Update paths before running on new projects!** See `10_realtime_processing/DEPLOYMENT_ISSUES.md`.

## Common Development Tasks

### Adding a New Processing Step

1. Choose appropriate number in sequence (gaps allow insertion)
2. Create script with numbered prefix: `25_new_step.sh`
3. Add epic header following existing style
4. Implement resume logic with file existence checks
5. Update README.md in relevant directory
6. Test with debug mode: `./script.sh -d`

### Debugging Processing Issues

```bash
# Enable debug mode
./24_process_reads_optimized.sh -i data -d

# Check status of all tools
./status.sh

# View processing logs
tail -f output/FLOWCELL/barcode01/log.txt

# Check failures
cat output/failed_files.txt

# Test single file
parallel --dry-run ... # Shows commands without running
```

### Working with HMM Databases

HMM databases must be:
- HMMER3 format (`.hmm` files)
- Uncompressed (not `.hmm.gz`)
- Contain trusted cutoffs (TC scores)
- Specified with **absolute paths** (not relative)

See `10_realtime_processing/HMM_SEARCH_GUIDE.md` for details.

### MAG Quality Standards (MIMAG)

- **High Quality:** >90% complete, <5% contamination, rRNA + tRNAs present
- **Medium Quality:** >50% complete, <10% contamination
- **Low Quality:** <50% complete, <10% contamination
- **Discard:** >10% contamination

CheckM2 reports completeness/contamination based on single-copy marker genes.

## Testing & Validation

### Before Expeditions

```bash
# Check all dependencies
./status.sh

# Test with small dataset
./24_process_reads_optimized.sh -i test_data/ -P

# Verify dashboard loads
Rscript 60_edna_mapping_viz.r
```

### During Expeditions

Monitor:
- Processing logs in real-time
- Dashboard updates
- Disk space (`df -h`)
- Memory usage (`free -h`)
- Failed files log

### After Expeditions

```bash
# Run complete MAG pipeline
cd 20_mag_assembly
./61_map_and_bin_optimized.sh

# Verify MAG quality
cat checkm2/quality_report.tsv

# Generate visualizations
Rscript 80_plot_bins.R
Rscript 82_inter_binning_analysis.r
```

## Important Documentation Files

Each directory contains detailed documentation:

**10_realtime_processing/:**
- `README.md` - User-facing workflow guide
- `CLAUDE.md` - Detailed architecture notes (for AI assistants)
- `RESUME_LOGIC.md` - Stage-aware resume behavior
- `HMM_SEARCH_GUIDE.md` - HMM integration details
- `OUTPUT_FORMATS.md` - Progress line formats
- `DEPLOYMENT_ISSUES.md` - Hardcoded paths, portability
- `CRASH_SAFETY.md` - Atomic operations, resume guarantees
- `CRITICAL_KRAKEN_BUG.md` - Memory management issue

**20_mag_assembly/:**
- `README.md` - Complete MAG assembly cookbook

**Root level:**
- `README.md` - Project overview with ASCII art
- `EPIC_TRANSFORMATION.md` - Transformation history
- `CHANGELOG.md` - Version history

When helping users, **always check relevant documentation first** - files contain detailed implementation notes and known issues.

## Common Pitfalls

1. **Using Kraken2 with wrong script** - Only use `-K` with `24_process_reads_optimized.sh`
2. **Relative paths for HMM files** - Must use absolute paths
3. **Forgetting to check dependencies** - Run `./status.sh` first
4. **Not using all three binners** - Use SemiBin2 + MetaBAT2 + MaxBin2 for best results
5. **Skipping polishing step** - Always polish MAGs before publication
6. **Insufficient coverage for binning** - Need 10-50x coverage and multiple samples with differential abundance
7. **Co-assembly of incompatible samples** - Don't co-assemble vastly different environments

## Script Selection Quick Reference

**Real-Time Processing:**
- Use `24_process_reads_optimized.sh` for most work (recommended)
- Use `22_process_reads_fast.sh` for urgent screening
- All other `2X_*.sh` scripts are legacy/experimental

**MAG Assembly:**
- Use `61_map_and_bin_optimized.sh` for complete pipeline
- Use `60_map_and_bin_complete.sh` as alternative
- Use individual scripts (10s, 20s, 30s...) for step-by-step control

## Getting Help

- Consult expert advisors: `./agents.sh`
- Check documentation in directory READMEs
- Review BUGFIX_*.md files for historical issue fixes
- Run `./status.sh` to verify tool installation
