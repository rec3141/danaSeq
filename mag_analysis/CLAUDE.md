# CLAUDE.md

## Project Overview

**dānaSeq MAG Analysis** is a technology-agnostic downstream analysis pipeline that accepts pre-computed assembly outputs (FASTA + depth table + optional BAMs) from any assembler and runs binning, annotation, taxonomy, metabolic profiling, mobile genetic element detection, eukaryotic analysis, ecosystem services mapping, phylogenetics, and interactive visualization.

Implemented in **Nextflow DSL2** in `nextflow/`.

## Repository Structure

```
mag_analysis/
└── nextflow/
    ├── main.nf              Pipeline entry point (40+ processes)
    ├── nextflow.config      Params, profiles, resources
    ├── run-mag-analysis.sh  Pipeline launcher (local/Docker/Apptainer)
    ├── install.sh           Conda env builder
    ├── modules/
    │   ├── binning.nf       BIN_SEMIBIN2, BIN_METABAT2, BIN_MAXBIN2,
    │   │                    BIN_LORBIN, BIN_COMEBIN, BIN_VAMB, BIN_VAMB_TAX,
    │   │                    DASTOOL_CONSENSUS, BINETTE_CONSENSUS,
    │   │                    MAGSCOT_CONSENSUS, CHECKM2
    │   ├── annotation.nf    PROKKA_ANNOTATE, BAKTA_BASIC, BAKTA_EXTRA
    │   ├── taxonomy.nf      KAIJU_CONTIG/CLASSIFY, KRAKEN2, SENDSKETCH
    │   ├── mge.nf           GENOMAD, CHECKV, INTEGRONFINDER, ISLANDPATH,
    │   │                    MACSYFINDER, DEFENSEFINDER
    │   ├── eukaryotic.nf    TIARA, WHOKARYOTE, METAEUK, MARFERRET
    │   ├── metabolism.nf    KOFAMSCAN, EMAPPER, DBCAN, MERGE_ANNOTATIONS,
    │   │                    MAP_TO_BINS, KEGG_MODULES, MINPATH, KEGG_DECODER,
    │   │                    ANTISMASH, ECOSSDB_MAP/SCORE/SDG/VIZ
    │   ├── rrna.nf          RNA_CLASSIFY (barrnap + vsearch + Aragorn)
    │   ├── phylogeny.nf     GTDBTK_CLASSIFY
    │   ├── gene_depths.nf   CALCULATE_GENE_DEPTHS (samtools bedcov from GFF)
    │   └── viz.nf           VIZ_PREPROCESS (4-stage incremental dashboard build)
    ├── bin/                  Pipeline scripts (merge, map, MinPath, etc.)
    ├── viz/                  Interactive Svelte + Vite dashboard
    ├── ecossdb/              Ecosystem services database (git submodule)
    ├── envs/                 Conda YAML specs (35+ environments)
    ├── Dockerfile            Thin layer on danaseq-mag-base-slim
    └── Dockerfile.base-slim  Heavy base image (9 merged conda envs)
```

## Running

```bash
cd nextflow
./install.sh && ./install.sh --check

# Basic run
./run-mag-analysis.sh \
    --assembly /path/to/assembly.fasta \
    --depths /path/to/depths.txt \
    --outdir /path/to/output \
    --annotator bakta --db_dir /path/to/databases

# Kitchen sink
./run-mag-analysis.sh \
    --assembly /path/to/assembly.fasta \
    --depths /path/to/depths.txt \
    --bam_dir /path/to/mapping/ \
    --outdir /path/to/output \
    --all --db_dir /path/to/databases \
    --sendsketch_address http://host:3068/sketch
```

## Required Inputs

| Flag | Required? | Format |
|------|-----------|--------|
| `--assembly` | YES | FASTA |
| `--depths` | YES | MetaBAT2-format TSV (from CoverM or jgi_summarize_bam_contig_depths) |
| `--bam_dir` | Optional | Directory with `*.sorted.bam` + `.bai` (for SemiBin2, LorBin, COMEBin) |

## Key Design Decisions

- **STAGE_INPUTS** process generates assembly_info.txt + gc.tsv from the input FASTA for viz compatibility
- **CALCULATE_GENE_DEPTHS** parses GFF (not TSV) for Prokka/Bakta compatibility
- **VIZ race condition fix**: only VIZ_STAGE4 (final) builds the static site; stages 1-3 write JSON only
- **Graceful degradation**: binners emit empty TSV on failure; pipeline continues
- **Dynamic binner architecture**: new binners add to `ch_binner_results` and automatically feed all three consensus methods
- **ECOSSDB** bundled as git submodule for ecosystem services profiling

## Adding a New Binner

1. Add process to `modules/binning.nf` outputting `LABEL_bins.tsv` (contig\tbin)
2. In `main.nf`, mix into `ch_binner_results`: `ch_binner_results = ch_binner_results.mix(BIN_NEW.out.bins.map { ['new', it] })`
3. Add `out.fastas` to CheckM2 and GTDB-Tk bin collection blocks
4. Add `--run_new` param to `nextflow.config`

## Database Locations (this system)

Databases live in `/data/scratch/refdbs/`. Use `--db_dir /data/scratch/refdbs` for auto-resolution, or see `download-databases.sh` for individual paths.
