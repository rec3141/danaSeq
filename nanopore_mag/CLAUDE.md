# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

**Important**: If user poses a question, agent responds with an answer, not a codebase change.

## Project Overview

**dānaSeq MAG Assembly** is a metagenome-assembled genome (MAG) reconstruction pipeline that runs alongside the real-time processing pipeline. It co-assembles nanopore reads (Flye, metaMDBG, or myloasm — selectable via `--assembler`), maps reads back, runs up to seven binning algorithms (SemiBin2, MetaBAT2, MaxBin2, LorBin, COMEBin, VAMB, VAMB taxvamb), and integrates results with three parallel consensus methods (DAS Tool, Binette, MAGScoT).

The pipeline is implemented in **Nextflow DSL2** in `nextflow/`. Legacy bash scripts are preserved in the parent directory for reference but are not actively maintained.

## Repository Structure

```
nanopore_mag/
├── nextflow/                    Primary pipeline (Nextflow DSL2)
│   ├── main.nf                 Pipeline entry point
│   ├── nextflow.config         Params, profiles, resources
│   ├── modules/
│   │   ├── preprocess.nf      CONCAT_READS (per-barcode concat + optional dedupe)
│   │   ├── assembly.nf         ASSEMBLY_FLYE, ASSEMBLY_METAMDBG, ASSEMBLY_MYLOASM, CALCULATE_TNF
│   │   ├── mapping.nf          MAP_READS, CALCULATE_DEPTHS, CALCULATE_GENE_DEPTHS
│   │   ├── binning.nf          BIN_SEMIBIN2, BIN_METABAT2, BIN_MAXBIN2,
│   │   │                       BIN_LORBIN, BIN_COMEBIN, BIN_VAMB, BIN_VAMB_TAX,
│   │   │                       DASTOOL_CONSENSUS, BINETTE_CONSENSUS,
│   │   │                       MAGSCOT_CONSENSUS, CHECKM2
│   │   ├── annotation.nf       PROKKA_ANNOTATE, BAKTA_BASIC, BAKTA_EXTRA
│   │   ├── taxonomy.nf         KAIJU_CONTIG_CLASSIFY, KAIJU_CLASSIFY,
│   │   │                       KRAKEN2_CLASSIFY, SENDSKETCH_CLASSIFY
│   │   ├── eukaryotic.nf      TIARA_CLASSIFY, WHOKARYOTE_CLASSIFY, METAEUK_PREDICT,
│   │   │                       MARFERRET_CLASSIFY
│   │   ├── rrna.nf             RNA_CLASSIFY (barrnap + vsearch + Aragorn)
│   │   ├── mge.nf              GENOMAD_CLASSIFY, CHECKV_QUALITY, INTEGRONFINDER,
│   │   │                       ISLANDPATH_DIMOB, MACSYFINDER, DEFENSEFINDER
│   │   ├── metabolism.nf       KOFAMSCAN, EMAPPER, DBCAN, MERGE_ANNOTATIONS,
│   │   │                       MAP_TO_BINS, KEGG_MODULES, MINPATH, KEGG_DECODER,
│   │   │                       ANTISMASH
│   │   ├── phylogeny.nf       GTDBTK_CLASSIFY (GTDB-Tk phylogenetic classification)
│   │   └── viz.nf              VIZ_PREPROCESS (dashboard JSON + static site build)
│   ├── envs/                   Conda YAML specs
│   │   ├── flye.yml            Flye, Filtlong, Nextflow, OpenJDK
│   │   ├── mapping.yml         minimap2, samtools, CoverM
│   │   ├── semibin.yml         SemiBin2, LorBin, PyTorch GPU
│   │   ├── semibin-cpu.yml     SemiBin2, LorBin, PyTorch CPU (for Docker)
│   │   ├── comebin.yml         COMEBin (rec3141 fork, PyTorch GPU)
│   │   ├── comebin-cpu.yml     COMEBin (CPU-only, for Docker)
│   │   ├── binning.yml         MetaBAT2, MaxBin2, DAS_Tool
│   │   ├── genomad.yml          geNomad (virus + plasmid detection)
│   │   ├── checkv.yml           CheckV (viral quality assessment)
│   │   ├── integron.yml        IntegronFinder (integron detection)
│   │   ├── islandpath.yml     IslandPath-DIMOB (genomic island detection)
│   │   ├── macsyfinder.yml    MacSyFinder (secretion systems + conjugation)
│   │   ├── defensefinder.yml  DefenseFinder (anti-phage defense systems)
│   │   ├── bakta.yml          Bakta (modern alternative to Prokka)
│   │   ├── kofamscan.yml      KofamScan (KEGG Orthology)
│   │   ├── emapper.yml        eggNOG-mapper (COG/GO/EC/KEGG/Pfam)
│   │   ├── dbcan.yml          dbCAN3 (CAZyme annotation)
│   │   ├── antismash.yml     antiSMASH (biosynthetic gene cluster detection)
│   │   ├── kaiju.yml            Kaiju (protein-level taxonomy)
│   │   ├── kraken2.yml         Kraken2 (k-mer contig-level taxonomy)
│   │   ├── prokka.yml          Prokka (gene annotation)
│   │   ├── checkm2.yml         CheckM2
│   │   ├── gtdbtk.yml          GTDB-Tk (phylogenetic MAG classification)
│   │   ├── tiara.yml           Tiara (eukaryotic contig classification, deep learning)
│   │   ├── whokaryote.yml     Whokaryote (eukaryotic classification, gene structure RF)
│   │   ├── metaeuk.yml        MetaEuk (eukaryotic gene prediction, multi-exon)
│   │   ├── marferret.yml      MarFERReT (DIAMOND + Python/pandas)
│   │   ├── rrna.yml            barrnap + vsearch (rRNA gene classification)
│   │   ├── bbmap.yml           BBMap (optional dedupe)
│   │   ├── vamb.yml            VAMB (variational autoencoder binner)
│   │   ├── binette.yml         Binette (consensus bin refinement)
│   │   ├── magscot.yml         MAGScoT (R + HMMER + prodigal, marker gene scoring)
│   │   └── viz.yml             Node.js + Python/pandas/scipy (viz dashboard)
│   ├── bin/                    Pipeline scripts (tetramer_freqs.py, islandpath_dimob.py,
│   │                           merge_annotations.py, map_annotations_to_bins.py,
│   │                           kegg_module_completeness.py, parse_marferret_results.py,
│   │                           parallel_defensefinder.py, parse_aragorn_results.py,
│   │                           parse_rrna_results.py, prepare_keggdecoder_input.py,
│   │                           run_minpath_per_mag.py)
│   ├── conda-envs/             Pre-built envs (created by install.sh)
│   ├── install.sh              Conda environment builder
│   ├── Dockerfile              Pipeline image (thin layer on base, rebuilt on push)
│   ├── Dockerfile.base         Base image (all conda envs + wrappers, rebuilt on-demand)
│   ├── entrypoint.sh           Docker entrypoint
│   ├── run-mag.sh              Pipeline launcher (local/Docker)
│   ├── seed-store-dir.sh       Seed storeDir from existing results (hardlink/symlink/copy)
│   ├── download-databases.sh   Database downloader (geNomad, CheckV, CheckM, CheckM2, etc.)
│   ├── viz/                    Interactive dashboard (Svelte + Vite)
│   │   ├── preprocess/         preprocess.py (TSV → JSON) + run_preprocess.sh
│   │   ├── public/              Static assets
│   │   │   ├── data/            JSON data files (generated by preprocess)
│   │   │   └── phylocanvas-bundle.js  Patched Phylocanvas.gl UMD bundle (see below)
│   │   ├── scripts/
│   │   │   └── patch-phylocanvas.sh   Patches Phylocanvas bundle for per-node label colors
│   │   └── src/                Svelte components + pages
│   ├── .dockerignore           Excludes conda-envs/, work/ from image
│   └── .gitignore              Excludes runtime artifacts
├── CLAUDE.md                   This file
├── 61_map_and_bin_optimized.sh Legacy orchestrator (reference only)
└── ...                         Other legacy scripts
```

## Running the Pipeline

```bash
cd nextflow

# Install conda environments
./install.sh
./install.sh --check

# Download databases (interactive menu or specify --genomad, --checkv, --checkm, --checkm2, --all)
./download-databases.sh --dir /data/scratch/refdbs

# Basic run (work dir defaults to /tmp/nanopore_mag_work)
./run-mag.sh --input /path/to/reads --outdir /path/to/output

# With Filtlong pre-filtering, no MaxBin2
./run-mag.sh --input /path/to/reads --outdir /path/to/output \
    --filtlong_size 40000000000 --run_maxbin false

# With persistent caching (storeDir) — survives work/ cleanup
./run-mag.sh --input /path/to/reads --outdir /path/to/output \
    --store_dir /data/scratch/mag_store

# Kitchen sink — compact form (--all enables all optional modules, --db_dir auto-resolves paths)
./run-mag.sh --input /data/minknow/QEI2025 \
    --outdir /data/minknow/QEI2025/nanopore_mag/tmpdir \
    --store_dir /data/minknow/QEI2025/nanopore_mag/final \
    --workdir /data/scratch/work \
    --all \
    --db_dir /data/scratch/refdbs \
    --sendsketch_address http://10.151.50.41:3068/sketch \
    --filtlong_size 40000000000 \
    --annotator bakta \
    --assembly_cpus 24 \
    --assembly_memory '120 GB'

# Kitchen sink — explicit form (same command, all flags spelled out for reference)
./run-mag.sh --input /data/minknow/QEI2025 \
    --outdir /data/minknow/QEI2025/nanopore_mag/tmpdir \
    --store_dir /data/minknow/QEI2025/nanopore_mag/final \
    --workdir /data/scratch/work \
    --dedupe \
    --filtlong_size 40000000000 \
    --annotator bakta \
    --bakta_light_db /data/scratch/refdbs/bakta/db-light \
    --bakta_db /data/scratch/refdbs/bakta/db \
    --bakta_extra \
    --genomad_db /data/scratch/refdbs/genomad_db \
    --checkv_db /data/scratch/refdbs/checkv_db \
    --checkm2_db /data/scratch/refdbs/checkm2 \
    --kaiju_db /data/scratch/refdbs/kaiju/refseq_ref \
    --run_kraken2 true \
    --kraken2_db /data/scratch/refdbs/krakendb/pluspfp_08gb \
    --run_sendsketch true \
    --sendsketch_address http://10.151.50.41:3068/sketch \
    --run_rrna true \
    --silva_ssu_db /data/scratch/refdbs/silva_db/SILVA_138.2_SSURef_NR99.fasta \
    --silva_lsu_db /data/scratch/refdbs/silva_db/SILVA_138.2_LSURef_NR99.fasta \
    --run_metabolism true \
    --kofam_db /data/scratch/refdbs/kofam_db \
    --eggnog_db /data/scratch/refdbs/eggnog_db \
    --dbcan_db /data/scratch/refdbs/dbcan_db \
    --macsyfinder_models /data/scratch/refdbs/macsyfinder_models \
    --defensefinder_models /data/scratch/refdbs/defensefinder_models \
    --run_antismash true \
    --antismash_db /data/scratch/refdbs/antismash_db \
    --run_eukaryotic true \
    --run_metaeuk true \
    --metaeuk_db /data/scratch/refdbs/metaeuk_db/metaeuk_db \
    --run_marferret true \
    --marferret_db /data/scratch/refdbs/marferret_db \
    --assembly_cpus 24 \
    --assembly_memory '120 GB'

# Docker mode (pre-built image from GHCR, no local build needed)
docker pull ghcr.io/rec3141/danaseq-mag:latest
./run-mag.sh --docker --input /path/to/reads --outdir /path/to/output

# HPC (Apptainer/Singularity)
apptainer pull danaseq-mag.sif docker://ghcr.io/rec3141/danaseq-mag:latest

# Local build (if modifying envs)
docker build -f Dockerfile.base -t danaseq-mag-base .
docker build -t danaseq-mag .

# Quick test with bundled data (direct nextflow)
mamba run -p conda-envs/dana-mag-flye \
    nextflow run main.nf --input test-data -profile test -resume

# Show all options
./run-mag.sh --help
```

### Re-running a Pipeline with -resume

`run-mag.sh` automatically records the Nextflow session ID in each run command saved to:

```
<outdir>/pipeline_info/run_command.txt
```

To re-run (e.g. after adding a new process or fixing a bug), use the **last line** of that file:

```bash
cd nextflow
bash -c "$(tail -1 <outdir>/pipeline_info/run_command.txt)"
```

The saved command is a `run-mag.sh` invocation (not raw mamba/nextflow) with
`--session <uuid>`, so cached tasks are always found — even after code changes.
To add new flags, just append them to the saved command.

**Session handling:** `run-mag.sh` handles sessions three ways:
1. **Auto-detect (default):** reads the last session ID from `run_command.txt` if it exists
2. **Explicit:** `--session <uuid>` to resume from a specific session
3. **Post-run capture:** after each run, extracts the actual session UUID from `.nextflow.log`
   and writes it into `run_command.txt` so future runs can resume from it

**Common pitfall:** changing `--input` path, `--assembly_cpus`, `--assembly_memory`, or other
params will invalidate the task hash and force a full re-run. Assembly takes **hours** on real
data (even "test" data with ~3 GB of reads). Always use the saved command from
`run_command.txt` and only append new flags — never change existing ones.

Only the processes whose script block changed will re-run; all others will be cached.

### Persistent Caching with storeDir

`-resume` caches are tied to the work directory and session ID — if `work/` is cleaned up,
all cached results are lost. `storeDir` provides persistent caching that survives across
runs, sessions, and work directory cleanup.

```bash
# Enable storeDir (opt-in, off by default)
./run-mag.sh --input /data/reads --outdir /data/output --store_dir /data/scratch/mag_store
```

**How it works:** When `--store_dir` is set, each process stores its outputs directly in
the store directory. On subsequent runs, if all declared output files already exist in
storeDir, the process is skipped entirely (Nextflow shows "Stored" status). No work
directory task is created.

**Key behavior:** When storeDir is active, publishDir is NOT ignored by Nextflow — both
run. To avoid doubling disk usage, all processes use `publishDir mode: 'link'` (hardlinks).
On the same filesystem, hardlinks cost zero extra disk space. Do not change publishDir
mode to 'copy' when using --store_dir.

**Seeding from existing results:** To populate a storeDir from a previous run without
re-executing the pipeline:

```bash
# Same filesystem → hardlinks (zero extra disk space)
./seed-store-dir.sh full_test_20260216 /data/mag_store

# Cross-filesystem → auto-detects and uses symlinks
./seed-store-dir.sh full_test_20260216 /data/scratch/mag_store

# Force a specific mode
./seed-store-dir.sh --mode copy full_test_20260216 /data/scratch/mag_store
```

The seed script handles legacy naming conventions (e.g. `contig_bins.tsv` →
`{name}_bins.tsv`, `PROKKA_*` → `annotation.*`).

### Verifying Process Output (Testing)

When testing a new or modified Nextflow process, **always check the raw output in the work
directory** before trusting the published results. Published files go through `publishDir`
copy/link logic which can silently produce empty or incomplete outputs.

```bash
# Find the work directory from the Nextflow log
grep 'PROCESS_NAME' .nextflow.log | tail -1   # shows [ab/cd1234] hash

# Check the actual outputs in the work directory
ls -la work/ab/cd1234*/                         # list all files
cat work/ab/cd1234*/.command.log                # stdout
cat work/ab/cd1234*/.command.err                # stderr
wc -l work/ab/cd1234*/*.tsv                     # line counts of output files
```

If published results are empty but work directory has data, the bug is in the `publishDir`
directive or in the output file copy logic (e.g. wrong filename pattern in the script block).

## Pipeline Architecture

### Processing DAG (fan-in → fan-out → fan-in)

```
Sample FASTQs (N files)
         │ CONCAT_READS (per barcode)
         │ collect()
   ASSEMBLY_{FLYE|METAMDBG|MYLOASM}   Fan-in: all reads → 1 co-assembly (--assembler)
         ├──────────────────────┬──────────────────────────┬──────────────────────────────────┐
   MAP_READS (×N)     CALCULATE_TNF   GENOMAD_CLASSIFY     INTEGRONFINDER                    │
         │ collect()        │                │              KAIJU_CONTIG_CLASSIFY              │
   CALCULATE_DEPTHS        │         CHECKV_QUALITY        KRAKEN2_CLASSIFY                   │
   CALCULATE_GENE_DEPTHS   │                               SENDSKETCH_CLASSIFY                │
                            │                               RNA_CLASSIFY                      │
                      PROKKA|BAKTA ─────────────────────────────────────────────────────┐     │
                            │                                                           │     │
                  ┌─────────┼──────────┬─────────────────────┐                          │     │
            ISLANDPATH   KOFAMSCAN  EMAPPER  DBCAN   TIARA + WHOKARYOTE   METAEUK → MARFERRET │
            MACSYFINDER  DEFENSEFINDER │              ANTISMASH (assembly + GFF)              │
            KAIJU_CLASSIFY   └─────┬───┘                                                      │
                          MERGE_ANNOTATIONS                                                   │
                                  │                                                           │
                           MAP_TO_BINS ──────────── (needs DAS_Tool contig2bin)                │
                           ┌──────┼──────┐                                                    │
                     KEGG_MODULES MINPATH KEGG_DECODER                                        │
         │                                                                                    │
    ┌────┼────┬────┬────┬────┬────┐                                                           │
 SemiBin2 MetaBAT2 MaxBin2 LorBin COMEBin VAMB VAMB_TAX   Binning (parallel, all optional)   │
    └────┼────┴────┴────┴────┴────┘                                                           │
   ┌─────┼──────────────┐                                                                     │
 DAS_Tool  Binette  MAGScoT   3 parallel consensus methods                                   │
   └─────┼──────────────┘                                                                     │
   CHECKM2                Quality assessment (optional, needs --checkm2_db)                   │
   GTDBTK_CLASSIFY        GTDB-Tk phylogenetic classification (optional, needs --gtdbtk_db)   │
         │                                                                                    │
   VIZ_PREPROCESS (×4)    Incremental dashboard (barrier: TNF → annotation → binning → all)   │
```

### Key Design Decisions

**CoverM for depth calculation:** Replaces `jgi_summarize_bam_contig_depths` which has an integer overflow bug in MetaBAT2 <=2.17 on long reads. CoverM handles supplementary alignments correctly and outputs MetaBAT2-compatible depth tables.

**Supplementary alignment filtering:** MAP_READS uses `samtools view -F 0x904` to drop unmapped, secondary, and supplementary alignments before sorting. Long reads produce chimeric alignments that cause massive depth overcounting.

**Dynamic binner architecture:** Binners emit `[label, file]` tuples that are mixed into a single channel (`ch_binner_results`) and collected for all three consensus methods (DAS Tool, Binette, MAGScoT). New binners can be added by appending to `ch_binner_results` in `main.nf` — no changes needed in the consensus processes. The three consensus methods run in parallel after all binners complete, and all outputs (raw binner + consensus) go to CheckM2 and GTDB-Tk for quality assessment.

**Graceful failure handling:** SemiBin2 catches crashes on small datasets (0 bins → empty ORFs → hmmsearch fail) and produces an empty output file. DAS_Tool filters out empty binner inputs and handles the "no bins above score threshold" case. The pipeline completes successfully even when individual binners fail.

**GPU vs CPU PyTorch:** The local conda envs (`semibin.yml`, `comebin.yml`) include `pytorch-gpu` for GPU-accelerated SemiBin2, LorBin, and COMEBin. The Docker image uses the `-cpu.yml` variants since most deployments won't have `--gpus`.

**Split Docker build:** The Docker image is split into a heavy base image (`Dockerfile.base`, all conda envs, rebuilt on-demand via `workflow_dispatch`) and a thin pipeline image (`Dockerfile`, just code, rebuilt on every push to main in <5 min). Both are published to `ghcr.io/rec3141/`. The base image uses a `CHECKM_DATA_DIR` sentinel to skip the CheckM v1 280 MB post-link download at build time; the real data should be provided via `download-databases.sh --checkm`.

**Resume:** Nextflow's built-in `-resume` uses task hashing. No manual checkpoint logic needed.

**Persistent caching:** All processes support `storeDir` for caching that survives work directory cleanup. Opt-in via `--store_dir`. When active, publishDir is ignored and outputs go directly to the store path.

### Conda Environments

Isolated conda environments avoid dependency conflicts:

| Environment | Tools | Rationale |
|-------------|-------|-----------|
| `dana-mag-flye` | Flye, Filtlong, Nextflow, OpenJDK | Python version conflicts; also hosts Nextflow runtime |
| `dana-mag-mapping` | minimap2, samtools, CoverM | Universal mapping tools |
| `dana-mag-semibin` | SemiBin2, LorBin, PyTorch GPU | ML dependencies isolated; LorBin shares PyTorch |
| `dana-mag-comebin` | COMEBin (rec3141 fork), PyTorch GPU | Contrastive learning binner; cloned from fork at install |
| `dana-mag-binning` | MetaBAT2, MaxBin2, DAS_Tool | Binning suite |
| `dana-mag-genomad` | geNomad | Virus + plasmid + provirus detection (neural network) |
| `dana-mag-checkv` | CheckV | Viral genome quality assessment |
| `dana-mag-integron` | IntegronFinder | Integron + gene cassette detection (attC/attI + HMM) |
| `dana-mag-islandpath` | Python + HMMER | Genomic island detection via dinucleotide bias (Python reimplementation) |
| `dana-mag-macsyfinder` | MacSyFinder v2 | Secretion systems (TXSScan) + conjugation (CONJScan) |
| `dana-mag-defensefinder` | DefenseFinder | Anti-phage defense systems (CRISPR, R-M, BREX, Abi, etc.) |
| `dana-mag-kofamscan` | KofamScan, HMMER | KEGG Orthology assignment via adaptive HMM thresholds |
| `dana-mag-emapper` | eggNOG-mapper, DIAMOND | COG/GO/EC/KEGG/Pfam functional annotation |
| `dana-mag-dbcan` | run_dbcan, HMMER, DIAMOND | CAZyme annotation (3-method consensus) |
| `dana-mag-kaiju` | Kaiju | Protein-level taxonomy (six-frame translation or pre-annotated .faa) |
| `dana-mag-kraken2` | Kraken2 | k-mer contig-level taxonomy (no annotation needed; maxForks 1) |
| `dana-mag-prokka` | Prokka | Gene annotation (alternative to Bakta) |
| `dana-mag-checkm2` | CheckM2 | Quality assessment (optional, needs `--checkm2_db`) |
| `dana-mag-gtdbtk` | GTDB-Tk 2.4.0 | Phylogenetic MAG classification (optional, needs `--gtdbtk_db`) |
| `dana-mag-tiara` | Tiara | Deep learning k-mer eukaryotic classification (98%+ accuracy) |
| `dana-mag-whokaryote` | Whokaryote, Prodigal | Gene structure-based eukaryotic classification (random forest) |
| `dana-mag-metaeuk` | MetaEuk | Eukaryotic gene prediction (multi-exon, intron-aware, homology-based) |
| `dana-mag-marferret` | DIAMOND, Python, pandas | Marine eukaryotic taxonomy + Pfam via MarFERReT |
| `dana-mag-rrna` | barrnap, vsearch | rRNA gene detection (barrnap) + SILVA classification (vsearch) |
| `dana-mag-viz` | Node.js, Python, pandas, scipy | Dashboard preprocessing (TSV→JSON) + static site build |
| `dana-mag-metamdbg` | metaMDBG | Long-read metagenome assembler (alternative to Flye) |
| `dana-mag-myloasm` | myloasm | Long-read assembler |
| `dana-mag-derep` | galah, skani, sourmash | Fast MAG dereplication (ANI-based) |
| `dana-mag-drep` | dRep | Quality-aware MAG dereplication (heavy: CheckM v1, mummer4) |
| `dana-mag-instrain` | inStrain | Strain-level population genomics |
| `dana-mag-strainy` | strainy | Strain-aware assembly phasing |
| `dana-mag-floria` | floria | Strain-resolved metagenome analysis |
| `dana-mag-skder` | skder | Genome database dereplication |
| `dana-mag-antismash` | antiSMASH 8.0.4 | Biosynthetic gene cluster detection (NRPS, PKS, terpenes, etc.) |
| `dana-mag-pathway` | MinPath, KEGG-Decoder | Pathway analysis (parsimony reconstruction + biogeochemical scoring) |
| `dana-bbmap` | BBMap | Optional dedupe (only if `params.dedupe`) |
| `dana-mag-vamb` | VAMB | Variational autoencoder binner (depth-based) |
| `dana-mag-binette` | Binette 1.2.1 | Consensus bin refinement via set operations + internal CheckM2 |
| `dana-mag-magscot` | R, HMMER, prodigal | MAGScoT consensus via marker gene scoring (R script from GitHub) |

### Nextflow Config Profiles

| Profile | Use case |
|---------|----------|
| `standard` | Local execution |
| `test` | Small test data, reduced resources (4 CPUs, 8GB) |

## Expected Input/Output

### Input

```
input_dir/
├── sample1.fastq.gz
├── sample2.fastq.gz
└── ...
```

`--input` must point to a directory containing `*.fastq.gz` files. All files are co-assembled.

### Output

```
results/
├── assembly/
│   ├── assembly.fasta         Co-assembly
│   ├── tnf.tsv                Tetranucleotide frequencies (136 features)
│   └── gc.tsv                 Per-contig GC content
├── mapping/
│   ├── *.sorted.bam           Per-sample alignments
│   ├── *.sorted.bam.bai       BAM indices
│   ├── depths.txt             CoverM depth table
│   └── gene_depths.tsv        Per-gene mean depths (from annotation BED + BAMs)
├── binning/
│   ├── semibin/semibin_bins.tsv
│   ├── metabat/metabat_bins.tsv
│   ├── maxbin/maxbin_bins.tsv
│   ├── lorbin/lorbin_bins.tsv
│   ├── comebin/comebin_bins.tsv
│   ├── vamb/vamb_bins.tsv              (if --run_vamb)
│   ├── vamb_tax/vamb_tax_bins.tsv      (if --run_vamb_tax + --sendsketch_address)
│   ├── dastool/
│   │   ├── bins/*.fa           Final consensus MAG FASTAs
│   │   ├── contig2bin.tsv      Contig-to-bin assignments
│   │   ├── allbins.fa          All bins concatenated
│   │   ├── bin_quality.tsv     Per-bin SCG completeness/redundancy (all binners)
│   │   ├── summary.tsv         DAS_Tool consensus winners with scores
│   │   ├── bacteria.scg        Bacterial single-copy gene assignments (protein_id\tSCG_name)
│   │   └── archaea.scg         Archaeal single-copy gene assignments (protein_id\tSCG_name)
│   ├── binette/                         (if --run_binette + --checkm2_db)
│   │   ├── binette_bins.tsv    Binette consensus contig-to-bin
│   │   └── bins/*.fa           Binette refined MAG FASTAs
│   ├── magscot/                         (if --run_magscot)
│   │   ├── magscot_bins.tsv    MAGScoT consensus contig-to-bin
│   │   └── bins/*.fa           MAGScoT refined MAG FASTAs
│   └── checkm2/
│       └── quality_report.tsv  CheckM2 quality for all binners + consensus (if --checkm2_db)
├── mge/                           MGE detection (if --genomad_db set)
│   ├── genomad/
│   │   ├── virus_summary.tsv      Virus contigs with scores + taxonomy
│   │   ├── plasmid_summary.tsv    Plasmid contigs with scores
│   │   ├── virus.fna              Viral contig sequences
│   │   ├── plasmid.fna            Plasmid contig sequences
│   │   ├── virus_proteins.faa     Virus protein sequences
│   │   ├── plasmid_proteins.faa   Plasmid protein sequences
│   │   ├── virus_genes.tsv        Per-gene annotations (markers, AMR, taxonomy)
│   │   ├── plasmid_genes.tsv      Per-gene annotations (markers, conjugation, AMR)
│   │   ├── provirus.tsv           Provirus boundaries + integrase calls
│   │   ├── provirus.fna           Excised provirus sequences
│   │   ├── taxonomy.tsv           Per-contig taxonomy assignments
│   │   └── genomad_summary.tsv    Per-contig classification scores
│   ├── checkv/                    Viral QA (if --checkv_db set)
│   │   ├── quality_summary.tsv    Completeness + contamination
│   │   ├── viruses.fna            Host-trimmed viral sequences
│   │   └── proviruses.fna         Extracted provirus sequences
│   ├── integrons/                 Integron detection (if --run_integronfinder)
│   │   ├── integrons.tsv          Per-element annotations (integrase, attC, attI, cassettes)
│   │   └── summary.tsv            Counts of complete/In0/CALIN integrons per contig
│   ├── islandpath/                Genomic island detection (if --run_islandpath)
│   │   └── genomic_islands.tsv    Island coordinates (id, contig, start, end)
│   ├── macsyfinder/               Secretion + conjugation (if --macsyfinder_models)
│   │   ├── all_systems.tsv        Detected systems with component hits
│   │   └── all_systems.txt        Human-readable system descriptions
│   └── defensefinder/             Anti-phage defense (if --run_defensefinder)
│       ├── systems.tsv            Detected defense systems (CRISPR, R-M, BREX, etc.)
│       ├── genes.tsv              Per-gene assignments within systems
│       └── hmmer.tsv              Raw HMM hits across all models
├── metabolism/                    Metabolic profiling (if --run_metabolism)
│   ├── kofamscan/
│   │   └── kofamscan_results.tsv  Per-protein KO assignments (adaptive threshold)
│   ├── emapper/
│   │   └── emapper_results.emapper.annotations  COG/GO/EC/KEGG/Pfam
│   ├── dbcan/
│   │   └── overview.tsv           CAZyme consensus (≥2/3 methods)
│   ├── merged/
│   │   └── merged_annotations.tsv Unified per-protein annotation table
│   ├── per_mag/
│   │   └── *.tsv                  Per-MAG annotation tables
│   ├── modules/
│   │   ├── module_completeness.tsv  MAG × module completeness matrix
│   │   └── module_heatmap.svg       Clustered heatmap
│   ├── minpath/
│   │   ├── minpath_pathways.tsv     MAG × pathway (naive vs parsimony)
│   │   └── details/                 Per-MAG MinPath reports
│   ├── kegg_decoder/
│   │   ├── kegg_decoder_output.tsv  MAG × function completeness (~80 functions)
│   │   └── function_heatmap.svg     Publication-quality heatmap
│   ├── antismash/                         (if --run_antismash)
│   │   ├── antismash_summary.tsv    Per-region BGC summary (type, known cluster, similarity)
│   │   ├── antismash_geneclusters/  Region GenBank files
│   │   └── antismash_json/          Full antiSMASH JSON output
│   └── community/
│       └── community_annotations.tsv  All proteins with bin_id column
├── taxonomy/                      Taxonomy classification
│   ├── kaiju/                     Protein-level taxonomy (if --kaiju_db set)
│   │   ├── kaiju_genes.tsv        Per-gene Kaiju classifications
│   │   └── kaiju_contigs.tsv      Per-contig taxonomy (majority vote)
│   ├── kraken2/                   k-mer taxonomy (if --kraken2_db set)
│   │   ├── kraken2_contigs.tsv    Per-contig Kraken2 classifications + lineage
│   │   └── kraken2_report.txt     Standard Kraken2 report (for Krona/Pavian)
│   ├── sendsketch/                GTDB MinHash taxonomy (if --sendsketch_address set)
│   │   └── sendsketch_contigs.tsv Per-contig GTDB taxonomy + ANI
│   ├── rrna/                      RNA gene classification (if --silva_ssu_db set)
│   │   ├── rrna_genes.tsv         Per-gene rRNA classifications (barrnap + vsearch)
│   │   ├── rrna_contigs.tsv       Per-contig rRNA summary (best SSU/LSU taxonomy)
│   │   ├── rrna_sequences.fasta   Extracted rRNA gene sequences
│   │   └── trna_genes.tsv         tRNA/tmRNA genes (Aragorn)
│   └── gtdbtk/                    GTDB-Tk classification (if --gtdbtk_db set)
│       ├── gtdbtk_taxonomy.tsv    Per-bin GTDB classification + placement info
│       └── gtdbtk_trees/          Newick placement trees (bac120/ar53)
├── eukaryotic/                    Eukaryotic analysis (if --run_eukaryotic or --run_metaeuk)
│   ├── tiara/
│   │   └── tiara_output.tsv       Per-contig Tiara classification + probabilities
│   ├── whokaryote/
│   │   └── whokaryote_classifications.tsv  Per-contig Whokaryote classification
│   ├── metaeuk/                   Eukaryotic gene prediction (if --metaeuk_db set)
│   │   ├── metaeuk_proteins.fas    Multi-exon eukaryotic protein predictions
│   │   ├── metaeuk_codon.fas       Nucleotide coding sequences
│   │   ├── metaeuk.gff             Gene structures (exon boundaries)
│   │   └── metaeuk_headers.tsv     Internal ID mapping
│   └── marferret/                 Marine eukaryotic taxonomy (if --marferret_db set)
│       ├── marferret_proteins.tsv  Per-protein MarFERReT taxonomy + Pfam
│       └── marferret_contigs.tsv   Per-contig aggregated taxonomy + Pfam domains
├── viz/                           Interactive dashboard (if --run_viz)
│   ├── data/                      12 JSON files for dashboard
│   │   ├── overview.json          Assembly + binning summary stats
│   │   ├── mags.json              Per-MAG metrics (completeness, contamination, taxonomy)
│   │   ├── checkm2_all.json       CheckM2 quality for all binners
│   │   ├── kegg_heatmap.json      KEGG module completeness matrix
│   │   ├── coverage.json          Per-sample contig coverage
│   │   ├── taxonomy_sunburst.json Taxonomy hierarchy for sunburst plot
│   │   ├── phylotree.json         GTDB-Tk placement trees + per-bin classification
│   │   └── ...                    mge_summary, mge_per_bin, eukaryotic, contig_explorer, contig_lengths
│   └── site/                      Static site (index.html + assets/)
└── pipeline_info/
    ├── run_command.txt         Exact re-runnable command (for -resume)
    ├── timeline.html
    ├── report.html
    └── trace.txt
```

## Development Notes

### Running the Pipeline (IMPORTANT)

**Always use `./run-mag.sh` to launch the pipeline, never `nextflow run main.nf` directly.**
`run-mag.sh` records every invocation to `<outdir>/pipeline_info/run_command.txt`, which is
essential for resuming runs with the exact same parameters. Running Nextflow directly skips
this logging, making it impossible to reliably resume later.

### Adding a New Binner

1. Add a process to `modules/binning.nf` that outputs `path("LABEL_bins.tsv")` in DAS_Tool format (contig\tbin, tab-separated)
2. In `main.nf`, invoke the process and mix its output into `ch_binner_results`:
   ```groovy
   BIN_NEWTOOL(ch_assembly, CALCULATE_DEPTHS.out.jgi_depth)
   ch_binner_results = ch_binner_results.mix(
       BIN_NEWTOOL.out.bins.map { ['newtool', it] }
   )
   ```
3. Add the binner's `out.fastas` to `ch_all_bins` (CheckM2) and `ch_gtdbtk_bins` (GTDB-Tk)
4. Add a conda env YAML if needed, or add the tool to an existing env
5. Add a `--run_newtool` param to `nextflow.config` if it should be optional
6. Add `'newtool'` to the binner lists in `viz/preprocess/preprocess.py`

Binners added to `ch_binner_results` automatically feed into all three consensus methods
(DAS Tool, Binette, MAGScoT). Consensus methods run in parallel after all binners complete.

### Modifying an Existing Module

Module files are in `nextflow/modules/*.nf`. Each process has:
- `conda` directive pointing to the pre-built env in `conda-envs/`
- `publishDir` for output routing
- `storeDir` for persistent caching (active when `params.store_dir` is set)
- Resource labels (`process_low`, `process_medium`, `process_high`, `process_kraken`)
- Graceful error handling for edge cases (empty output, tool crashes)

### Known Limitations

- **Small test data:** The bundled test data (~17MB) is too small for meaningful binning. SemiBin2 produces 0 bins, DAS_Tool finds no bins above the 0.5 score threshold. All handled gracefully.
- **DAS_Tool hanging:** On very small assemblies, DAS_Tool's ruby SCG annotation script can hang. The `set +e` error handling catches this on retry.
- **Docker GPU:** The Docker image uses CPU-only PyTorch. For GPU SemiBin2 in Docker, build with `semibin.yml` instead of `semibin-cpu.yml` and use `--gpus all`.

## Development Rules

### NEVER pollute the source tree with pipeline results or data files
- **NEVER** write preprocess output, pipeline results, or data files into `viz/public/data/` or any directory under `nextflow/viz/`
- Pipeline results go to `/data/scratch/<pipeline_name>/<unique_run_name>/`
- For vite preview/dev testing, serve data from `/tmp/viz/<unique_run_name>/`
- The `viz/public/data/` directory should only contain committed placeholder/example files, NOT real run data
- Mixing data from different runs in the dev directory causes subtle ID mismatch bugs

## Viz Dashboard: Phylocanvas.gl

The Phylogeny tab renders GTDB-Tk Newick placement trees (9 trees, 13K–44K leaves each)
using [Phylocanvas.gl](https://www.phylocanvas.gl/) (MIT, WebGL, handles 1M+ nodes).

### Architecture

- **Phylocanvas.gl** is loaded as a standalone UMD bundle (`public/phylocanvas-bundle.js`)
  via a `<script>` tag in `index.html`, NOT as an ESM import. Vite's CJS→ESM conversion
  adds strict mode which triggers a deck.gl bug (`TypeError: setting getter-only property
  "componentName"`). The UMD bundle bypasses Vite entirely and sets `window.phylocanvas.PhylocanvasGL`.
- **PhylocanvasTree.svelte** wraps the WebGL tree. It creates/destroys tree instances via
  `$effect`, handles click/hover via `tree.pickNodeFromLayer()`, and uses `ResizeObserver`
  to resize the canvas when the info panel opens/closes.
- **PhyloTreeView.svelte** provides the toolbar (tree selector, layout/color cycle, filters),
  computes per-node `styles` (colored bins, gold reference species), and shows a side panel
  with Quality-tab-style bin details on click.

### Per-node label color patch

Phylocanvas.gl renders all leaf labels with a single global `fontColour`. We patch the
bundle to make the leaf-label TextLayer use each node's `fillColour` (set via the `styles`
prop) when available, falling back to `fontColour`. This enables:
- **Gold labels** for validly-named reference species (e.g. "Brevinema andersonii")
- **Bright palette-colored labels** for user MAG bins
- **Dim labels** for unnamed placeholder references (e.g. "SKYB106 sp025061835")

The patch is applied automatically by `scripts/patch-phylocanvas.sh`, which runs via the
`postinstall` npm hook. It makes two sed replacements on the minified UMD bundle:
1. Changes `getColor` from static `fontColour` to `function(d){return d.fillColour||fc}`
2. Adds `getColor` to `updateTriggers` so colors update when styles change

If Phylocanvas releases a new version that changes the bundle format, the patch script
will warn but not fail. Check the sed patterns in `scripts/patch-phylocanvas.sh`.

### preprocess.py: phylotree.json

`build_phylotree()` produces `phylotree.json` containing:
- `newick` — dict of 9 Newick strings (backbone + 7 bac120 class-level + 1 ar53), with
  reference genome accessions relabeled to GTDB species names (requires `--gtdbtk-db`)
- `bins` — per-bin GTDB-Tk classification, lineage, method, closest genome reference/ANI/AF
- `leaf_info` — display label → {accession, lineage} for reference taxa (powers click info)
- `tree_metadata` — per-tree user bin counts (powers tree selector labels)
- `hierarchy` — d3-compatible taxonomy tree (legacy, unused by Phylocanvas)

## Common Issues

### No FASTQ files found
`--input` must point to a directory containing `*.fastq.gz` files directly (not subdirectories).

### Conda environment build fails
Run `./install.sh --check` to diagnose. Ensure mamba is on PATH. Build requires internet access.

**Channel order matters under strict priority:** Bioconda's official setup requires
`conda-forge` listed **before** `bioconda` in the channel list. With `channel_priority: strict`
(the modern default), the solver refuses to pull a dependency from a lower-priority channel
even if no compatible version exists in the higher one. If a YAML has `bioconda` first and
the package pins a dependency only available on conda-forge (e.g., antismash pins
`biopython==1.81`), the solve will fail with a misleading "no viable options" error. The fix
is to list `conda-forge` before `bioconda` in the YAML's `channels:` block. This matches
bioconda's own recommended configuration (`conda config --add channels bioconda` then
`conda config --add channels conda-forge`, which puts conda-forge highest).

### Depth values are wrong (overflow/negative)
This pipeline uses CoverM instead of `jgi_summarize_bam_contig_depths` to avoid the MetaBAT2 integer overflow bug. If you see overflow values, ensure the CALCULATE_DEPTHS process is using CoverM (check `modules/mapping.nf`).

### eggNOG-mapper / DIAMOND temp file management
DIAMOND writes intermediate alignments to unlinked temp files in the work directory
(visible via `/proc/<pid>/fd/` as `(deleted)` symlinks, invisible to `du`/`ls`). These
temp files grow at ~5 MB/s during search and can reach **20-30+ GB** for 133K proteins
against the 8.7 GB eggnog database.

**Fix (built-in batching):** The EMAPPER process uses two-stage batched mode:
- **Stage 1:** Splits proteins into chunks of `params.emapper_batch_size` (default 50K),
  runs DIAMOND search per chunk (`--no_annot`), bounds temp files to ~3-4 GB per chunk
- **Stage 2:** Annotates all merged seed orthologs at once (`-m no_search --dbmem`)

Performance flags:
- `--dmnd_iterate no` — single sensitivity pass instead of 4 rounds (eliminates temp bloat)
- `--dmnd_algo ctg` — contiguous seed algorithm, faster for protein search
- `--dbmem` — loads eggnog.db SQLite (~44 GB) into RAM for fast annotation (needs 48 GB)
- `--temp_dir .` — keeps temp files in the Nextflow work dir (explicit, not /tmp default)

**Result:** 61x speedup (7,950s → 129s for 100 proteins). EMAPPER memory is set to 48 GB
to accommodate `--dbmem`.

**Tuning batch size:** Adjust `--emapper_batch_size` to trade temp disk usage vs overhead.
Rule of thumb: 50K proteins → ~3-4 GB temp per chunk. Default 50K works for most systems.
For very constrained `/tmp`, lower to 25K. For systems with ample disk, raise to 100K+.

**Monitoring:** Check DIAMOND's hidden temp usage with:
```bash
# Total bytes in unlinked temp files held open by DIAMOND
for fd in $(ls /proc/<diamond_pid>/fd/); do
    grep pos /proc/<diamond_pid>/fdinfo/$fd 2>/dev/null
done | awk '{s+=$2} END {printf "%.1f GB\n", s/1024/1024/1024}'
```

## Database Locations (this system)

Reference databases live in `/data/scratch/refdbs/`, NOT in the pipeline's `databases/` directory. When downloading new databases, use:

```bash
./download-databases.sh --dir /data/scratch/refdbs
```

Current database paths for pipeline flags:
- `--bakta_db /data/scratch/refdbs/bakta/db`
- `--bakta_light_db /data/scratch/refdbs/bakta/db-light`
- `--genomad_db /data/scratch/refdbs/genomad_db`
- `--checkv_db /data/scratch/refdbs/checkv_db`
- `--checkm2_db /data/scratch/refdbs/checkm2`
- `--kaiju_db /data/scratch/refdbs/kaiju`
- `--kofam_db /data/scratch/refdbs/kofam_db`
- `--eggnog_db /data/scratch/refdbs/eggnog_db`
- `--dbcan_db /data/scratch/refdbs/dbcan_db`
- `--macsyfinder_models /data/scratch/refdbs/macsyfinder_models`
- `--defensefinder_models /data/scratch/refdbs/defensefinder_models`
- `--metaeuk_db /data/scratch/refdbs/metaeuk/orthodb_v11_euk/metaeuk_db`
- `--kraken2_db /data/scratch/refdbs/krakendb/pluspfp_08gb`
- `--silva_ssu_db /data/scratch/refdbs/silva_db/SILVA_138.2_SSURef_NR99.fasta`
- `--silva_lsu_db /data/scratch/refdbs/silva_db/SILVA_138.2_LSURef_NR99.fasta`
- `--marferret_db /data/scratch/refdbs/marferret_db`
- `--gtdbtk_db /data/scratch/refdbs/gtdbtk_db`
- `--antismash_db /data/scratch/refdbs/antismash_db`
- `--sendsketch_address http://10.151.50.41:3068/sketch`  (Ratnakara GTDB TaxServer)

---

**Repository:** https://github.com/rec3141/danaSeq
**License:** MIT (see LICENSE file)
