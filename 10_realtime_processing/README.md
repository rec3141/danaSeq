# Real-Time Processing Pipeline

Shipboard metagenomic analysis for Oxford Nanopore sequencing data. Processes reads as they stream from MinKNOW, providing live taxonomic classification, gene annotation, and functional profiling.

## Quick Start

### Nextflow (recommended)

```bash
cd nextflow

# Install conda environments (~15 min first time)
./install.sh
./install.sh --check

# Activate and run
conda activate conda-envs/dana-tools
nextflow run main.nf --input /path/to/nanopore/run \
    --run_kraken --kraken_db /path/to/db \
    --run_prokka --run_sketch --run_tetra \
    -resume

# Show all options
nextflow run main.nf --help
```

### Launcher script

```bash
cd nextflow

# Local conda (default)
./run-realtime.sh --input /path/to/data --outdir /path/to/output \
    --run_kraken --kraken_db /path/to/krakendb \
    --run_prokka

# Docker mode
docker build -t danaseq-realtime .
./run-realtime.sh --docker --input /path/to/data --outdir /path/to/output \
    --run_kraken --kraken_db /path/to/krakendb \
    --run_prokka

# Show all options
./run-realtime.sh --help
```

### Kitchen sink — all modules, all options with defaults

```bash
cd nextflow

# All analysis modules enabled, QC defaults shown explicitly
./run-realtime.sh --input /data/run1 --outdir /data/output \
    --run_kraken --kraken_db /path/to/krakendb \
    --run_prokka \
    --hmm_databases /path/to/CANT-HYD.hmm,/path/to/FOAM.hmm \
    --run_sketch \
    --run_tetra \
    --run_db_integration \
    --cleanup \
    --min_readlen 1500 \
    --keep_percent 80 \
    --min_file_size 1000000

# Watch mode — all modules, live sequencing
./run-realtime.sh --input /data/runs --outdir /data/output \
    --watch --db_sync_minutes 10 \
    --run_kraken --kraken_db /path/to/krakendb \
    --run_prokka \
    --hmm_databases /path/to/CANT-HYD.hmm \
    --run_sketch \
    --run_tetra \
    --run_db_integration
```

<details>
<summary>Manual docker run (without helper script)</summary>

```bash
docker run --user $(id -u):$(id -g) \
    -v /path/to/data:/data/input:ro \
    -v /path/to/output:/data/output \
    -v /path/to/krakendb:/kraken_db:ro \
    danaseq-realtime \
    run /pipeline/main.nf \
        --input /data/input --outdir /data/output \
        --kraken_db /kraken_db \
        --run_kraken --run_prokka -resume
```

Always use `--user $(id -u):$(id -g)` so output files are owned by your user.

</details>

## Pipeline Overview

```
MinKNOW FASTQ
    │
    ▼
VALIDATE_FASTQ          Gzip integrity + BBMap repair
    │
    ▼
QC_BBDUK                Adapter/quality trimming
    │
    ▼
QC_FILTLONG             Length filtering (>=1500bp), quality (Q7+, top 80%)
    │
    ▼
CONVERT_TO_FASTA        Header cleanup
    │
    ├──► KRAKEN2_CLASSIFY     Taxonomic classification (maxForks=1)
    ├──► PROKKA_ANNOTATE      ORF prediction + annotation
    ├──► HMM_SEARCH           HMMER3 against user-supplied databases
    ├──► SENDSKETCH           Rapid taxonomic sketching
    └──► TETRAMER_FREQ        Tetranucleotide composition
              │
              ▼
         DB_INTEGRATION        Load results into DuckDB
              │
              ▼
         CLEANUP               Compress/delete source files (optional)
```

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--input` | (required) | Directory **containing** `fastq_pass/` |
| `--outdir` | `results` | Output directory |
| `--run_kraken` | `false` | Enable Kraken2 classification |
| `--kraken_db` | (required if kraken) | Path to Kraken2 database |
| `--run_prokka` | `false` | Enable Prokka annotation |
| `--hmm_databases` | (skip) | Path to HMM file(s) for functional profiling |
| `--run_sketch` | `false` | Enable Sendsketch profiling |
| `--run_tetra` | `false` | Enable tetranucleotide frequency |
| `--watch` | `false` | Monitor for new files during live sequencing |
| `--watch_glob` | `*/fastq_pass/barcode*/*.fastq.gz` | Glob pattern for watch mode |
| `--db_sync_minutes` | `10` | DuckDB sync interval in watch mode |
| `--run_db_integration` | `false` | Load results into DuckDB |
| `--cleanup` | `false` | Compress/delete files after DB import |

### Profiles

| Profile | Use case |
|---------|----------|
| `standard` | Local execution, auto-detect resources |
| `test` | Small test files, reduced resources |
| `shipboard` | Production: 32 CPUs, 256 GB RAM, 100 GB Kraken |

## Watch Mode

Monitor a directory for new FASTQ files during active sequencing:

```bash
nextflow run main.nf --input /path/to/runs \
    --watch --db_sync_minutes 10 \
    --run_kraken --kraken_db /path/to/db \
    --run_prokka --run_db_integration
```

`DB_SYNC` runs as a long-lived process that periodically loads new results into DuckDB. R scripts are idempotent and track imports via `import_log`.

## Output

```
results/
├── FLOWCELL/
│   ├── barcode01/
│   │   ├── fa/           FASTA (quality-filtered sequences)
│   │   ├── fq/           Intermediate FASTQ (BBDuk output)
│   │   ├── kraken/       *.tsv (per-read), *.report (summary)
│   │   ├── prokka/       SAMPLE/PROKKA_* (GFF, FAA, FFN, TSV)
│   │   ├── hmm/          *.DBNAME.tsv, *.DBNAME.tbl
│   │   ├── sketch/       Sendsketch profiles
│   │   ├── tetra/        Tetranucleotide frequencies
│   │   └── log.txt
│   └── ...
├── dana.duckdb           Integrated database
└── pipeline_info/        Nextflow reports (timeline, trace, DAG)
```

## Post-DB Cleanup

The `--cleanup` flag compresses or deletes source files after DuckDB import:

| Directory/Files | Action |
|----------------|--------|
| `fa/*.fa` | Gzip in place (kept as compressed backup) |
| `kraken/`, `sketch/`, `tetra/` | Delete (data lives in DuckDB) |
| `prokka/*/PROKKA_*.tsv` | Delete (loaded into DuckDB) |
| `prokka/*/PROKKA_*.gff`, `.faa`, `.ffn` | Gzip in place |
| `hmm/`, `dana.duckdb`, `log.txt` | Kept (not cleaned) |

Safe for watch mode -- operates per-file and checks `import_log` before deleting.

## DuckDB Integration

R scripts in `nextflow/bin/` load results into DuckDB:

| Script | Data |
|--------|------|
| `40_kraken_db.r` | Kraken2 classifications |
| `41_krakenreport_db.r` | Kraken2 summary reports |
| `42_prokka_db.r` | Prokka annotations |
| `43_sketch_db.r` | Sendsketch profiles |
| `44_tetra_db.r` | Tetranucleotide frequencies |
| `45_stats_db.r` | Assembly statistics |
| `46_log_db.r` | Processing logs |
| `47_merge_db.r` | Merge per-run databases |

Scripts are idempotent and track imports via `import_log`.

## Design Notes

**Kraken2 serialization.** Uses `maxForks = 1` in the module definition. Kraken2 loads 50-100 GB into RAM, so only one instance runs at a time.

**Watch mode.** Uses `Channel.watchPath()` with a flat glob pattern (Java WatchService limitation -- no recursive `**`).

**Conda environments.** Three isolated environments avoid dependency conflicts:
- **dana-bbmap** -- BBMap (samtools conflicts with R)
- **dana-prokka** -- Prokka (BioPerl pins perl 5.26)
- **dana-tools** -- Filtlong, Kraken2, HMMER, R/DuckDB, Nextflow, OpenJDK

## Common Issues

**No FASTQ files found.** `--input` must point at the directory **containing** `fastq_pass/`. The pipeline gives targeted advice when files aren't found.

**Kraken2 out of memory.** Nextflow serializes Kraken2 via `maxForks = 1`. If running outside Nextflow, never launch multiple Kraken2 instances.

**Conda environment build fails.** Run `./install.sh --check` to diagnose. Ensure mamba/conda is on PATH.

**Watch mode not picking up files.** Java WatchService requires a flat glob pattern (no `**`). Adjust `--watch_glob` to match your directory depth.

## References

- BBTools: Bushnell B, [sourceforge.net/projects/bbmap](https://sourceforge.net/projects/bbmap/)
- Filtlong: Wick R, [github.com/rrwick/Filtlong](https://github.com/rrwick/Filtlong)
- Kraken2: Wood DE, Lu J, Langmead B, Genome Biology 2019
- Prokka: Seemann T, Bioinformatics 2014
- HMMER3: Eddy SR, PLoS Computational Biology 2011
- DuckDB: [duckdb.org](https://duckdb.org)
