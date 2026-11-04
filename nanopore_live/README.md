# Real-Time Processing Pipeline

Real-time metagenomic analysis for Oxford Nanopore sequencing data. Processes reads as they stream from MinKNOW, providing live taxonomic classification, gene annotation, and functional profiling.

## Storage conventions

Data for each sequencing run lives across three locations, by data tier:

| Tier | Location | Contents | Drive | Lifetime |
|---|---|---|---|
| **Raw reads** | `/data/<PROJECT>/<EXPERIMENT>/<RUN>/fastq_pass/` | Sequencer-output fastq files as they come off MinKNOW. Basecalled + demultiplexed per-barcode fastqs.  The working source for downstream pipelines (mapping, sketching, etc.). | SSD | Permanent |
| **Archived reads** | `/matika/seqs/metagenome/nanopore/reads/` | <RUN> folders are rsync'd here; run_logger.sh adds new folders to the run_log.tsv | NAS | Permanent. 
| **Pipeline outputs (active)** | `/data/scratch/nanopore_live/<RUN>/<FC>/<BC>` | Per-run workspace for `run-realtime.sh` / nextflow — sketches, viz JSONs, intermediate files. | SSD-backed for speed during active analysis. | SSD | Until run is completed |
| **Pipeline outputs (archived)** | `/vistara/dana/out_dana_bc/<RUN>/<FC>/<BC>` | Same per-run dirs after the pipeline is done — moved off the scratch SSD to free working space. Same layout, longer retention. | USB3 | Long-term (manual cleanup) |
| **Mapping references** | `/data/scratch/refdbs/mapping_refs/<ref>/` | `.fna` + `.idx` + sidecar files for each mapping reference (AIS / SAR / virus / etc.). See [docs/mapping-references.md](../docs/mapping-references.md). | SSD | Permanent |

### Archiving a completed run

**Automatic (watch mode):** add `--archive` to `run-realtime.sh`. The
launcher spawns `bin/watch_for_completion.sh`, which polls the input dir
for MinKNOW's `final_summary_*.txt` (one per flow cell). Once every
observed FC has one and nextflow has drained (no in-flight tasks, trace
stable), it SIGTERMs nextflow and invokes `bin/archive_run.sh`. After the
archive succeeds, it runs `nextflow clean <session> -f` to reclaim this run's
work dirs from the shared `-w` tree (scoped to the run's session UUID, so other
runs sharing the work dir are untouched).

```bash
./run-realtime.sh --input /data/<RUN> --outdir /data/scratch/nanopore_live/<RUN> \
    --watch --run_db_integration ... \
    --archive \
    --archive_raw_dest /matika/seqs/metagenome/nanopore/reads \
    --archive_out_dest /vistara/dana/out_dana_bc
```

The destination roots are configurable via flags or
`$ARCHIVE_RAW_DEST` / `$ARCHIVE_OUT_DEST` so the mechanism is portable
across sites.

**Manual** (for runs that finished before the flag existed):

```bash
RUN=2026_05_25_Winnipeg_River
# Verify nothing is writing
lsof +D /data/scratch/nanopore_live/$RUN | head
# Call archive_run.sh directly:
bin/archive_run.sh \
    --input  /data/$RUN \
    --outdir /data/scratch/nanopore_live/$RUN \
    --raw-dest /matika/seqs/metagenome/nanopore/reads \
    --out-dest /vistara/dana/out_dana_bc
```

Raw reads move to `/matika/`, processed reads and pipeline outputs move to `/vistara`,

### Per-barcode concatenation

Before the raw copy goes out, `archive_run.sh` calls `bin/concat_run_fastqs.sh`
on the run while it is still on local disk, which replaces MinKNOW's chunked `fastq_pass/<barcode>/*_0.fastq.gz …
*_571.fastq.gz` with a single `.fastq.gz` per barcode. MinKNOW writes plain
gzip, so the concatenation is an ordinary multi-member gzip that zlib, pigz,
zcat, BBTools and `fastq_filter` all read.

Nothing is deleted on trust:

1. every chunk is checked with `gzip -t` first — a chunk that is cut off, the
   common real failure, makes any decompressor stop there and report a short
   count, and the concatenated file would report the same short count, so the
   two would otherwise agree on a truncated total;
2. a cut-off chunk has its complete records salvaged into a repaired copy,
   which is used in its place (a chunk with nothing readable in it fails the
   barcode instead, since the archive holds the only copy);
3. the result must be exactly as long as the chunks that went into it (raw
   concatenation is byte-exact, so this costs nothing and catches a short
   write);
4. BBTools `reformat.sh` reads the finished file and must parse every record
   and report a non-zero count.

Each chunk is read once, for the append and the `gzip -t` in the same pass,
and the result is read once by `reformat.sh`. Doing this before the copy
means all of that happens on the SSD and only one file per barcode crosses
to the NAS. For a run that is already archived, `--work DIR` stages each
barcode to local disk, builds and verifies there, and copies just the result
back, comparing checksums afterwards since the chunks are about to go.

Each collapsed directory gets a `.concat_manifest.tsv` recording the chunks,
their sizes, the totals and any salvage. Its presence also marks the directory
as done, so re-running is a no-op. Run it by hand over an already-archived run
with:

```bash
bin/concat_run_fastqs.sh --run /matika/seqs/metagenome/nanopore/reads/<RUN>
```


## Quick Start

```bash
# Install conda environments (~15 min first time)
./install.sh
./install.sh --check

# Run (local conda, handles activation automatically)
./run-realtime.sh --input /path/to/nanopore/run --outdir /path/to/output \
    --run_kraken --kraken_db /path/to/krakendb \
    --annotator bakta --run_sketch --run_tetra

# Or with Docker
docker build -t danaseq-realtime .
./run-realtime.sh --docker --input /path/to/nanopore/run --outdir /path/to/output \
    --run_kraken --kraken_db /path/to/krakendb \
    --annotator bakta

# Show all options
./run-realtime.sh --help
```

### Kitchen sink — all modules, all options with defaults

```bash
# All analysis modules enabled, QC defaults shown explicitly
./run-realtime.sh --input /data/run1 --outdir /data/output \
    --run_kraken --kraken_db /path/to/krakendb \
    --annotator bakta \
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
    --annotator bakta \
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
        --run_kraken --annotator bakta -resume
```

Always use `--user $(id -u):$(id -g)` so output files are owned by your user.

</details>

## Pipeline Overview

```
MinKNOW FASTQ
    │
    ▼
VALIDATE_FASTQ          Gzip integrity check; damaged files salvaged (bin/fastq_gz_salvage.py)
    │
    ▼
QC_FASTQ_FILTER         Length + quality filtering (C, streaming)
    │
    ▼
CONVERT_TO_FASTA        Header cleanup
    │
    ├──► KRAKEN2_CLASSIFY     Taxonomic classification (batched per sample)
    ├──► BAKTA_CDS / PROKKA   Gene annotation (Bakta default)
    │       ├──► HMM_SEARCH         HMMER3 against user-supplied / pollutant HMM databases
    │       └──► DIAMOND_DB_SEARCH  DIAMOND blastp against .dmnd protein databases
    ├──► SENDSKETCH           Per-read GTDB classification (local sendsketch server)
    ├──► TETRAMER_FREQ        Tetranucleotide composition (C)
    └──► MAP_REFERENCE        minimap2 against --mapping_refs (from the filtered FASTQ)
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
| `--annotator` | `bakta` | Annotator: `bakta`, `prokka`, or `none` |
| `--hmm_databases` | (skip) | Path to HMM file(s) for functional profiling |
| `--hmm_cut` | `tc` | hmmsearch cutoff for HMM DBs (`tc`/`ga`/`nc`/E-value); per-DB `<name>.opts` / catalogue thresholds override |
| `--protein_dbs` | (skip) | DIAMOND `.dmnd` file(s) searched with blastp against annotation proteins |
| `--protein_db_opts` | `id=50,qcov=50,scov=50,evalue=1e-5` | DIAMOND cutoffs (HADEG-recommended); per-DB `<name>.opts` overrides |
| `--pollutant_dbs` | (skip) | Directory of pollutant DBs (`<name>/<name>.hmm` or `.dmnd` + sidecars), e.g. `/data/scratch/refdbs/pollutants` |
| `--run_sketch` | `false` | Enable per-read GTDB sendsketch classification |
| `--sendsketch_address` | `$DANA_SENDSKETCH_ADDRESS`, else `http://localhost:3068/sketch` | Local GTDB sendsketch server URL |
| `--run_tetra` | `false` | Enable tetranucleotide frequency |
| `--watch` | `false` | Monitor for new files during live sequencing (watch globs are auto-detected from the depth of `fastq_pass/` below `--input`) |
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
    --annotator bakta --run_db_integration
```

`DB_SYNC` runs as a long-lived process that periodically loads new results into DuckDB. The Python loaders in `bin/` are idempotent and track imports via `import_log`. Watch vs batch mode is locked per `--outdir` (recorded in `<outdir>/.pipeline_mode`) since switching modes changes per-file vs per-sample grouping and would corrupt the database.

## Output

```
results/
├── FLOWCELL/
│   ├── barcode01/
│   │   ├── fa/           FASTA (quality-filtered sequences)
│   │   ├── kraken/       *.tsv (per-read), *.report (summary)
│   │   ├── bakta/        Bakta annotations (GFF3, FAA, FFN, TSV)
│   │   ├── hmm/          *.DBNAME.tsv, *.DBNAME.tbl
│   │   ├── protdb/       *.DBNAME.tsv, *.DBNAME.tbl (DIAMOND)
│   │   ├── sketch/       *.sendsketch_reads.tsv (per-read GTDB classification)
│   │   ├── tetra/        *.lrn, *.lengths (tetranucleotide frequencies)
│   │   ├── map/<file>/   <ref>.txt (minimap2 hits per input file and reference)
│   │   └── dana.duckdb   Integrated database (one per barcode)
│   └── ...
├── viz/                  Dashboard JSON (watch mode; DB_SYNC regenerates each tick)
├── run_manifest.json     Pipeline version, parameters, containers
└── pipeline_info/        Nextflow reports (timeline, trace, DAG)
```

## Post-DB Cleanup

The `--cleanup` flag compresses or deletes source files after DuckDB import:

| Directory/Files | Action |
|----------------|--------|
| `kraken/`, `sketch/`, `tetra/` (except `tetra/*.lengths`) | Delete (data lives in DuckDB) |
| `prokka/*/PROKKA_*.tsv`, `bakta/*/*.tsv` | Delete (loaded into DuckDB) |
| `fa/*.fa`, `.gff`/`.gff3`, `.faa`, `.ffn` | Kept uncompressed (Nextflow `publishDir` would re-copy them on `-resume`) |
| `hmm/`, `protdb/`, `map/`, `dana.duckdb` | Kept (not cleaned) |

Cleanup runs only for barcode directories whose `import_log` is non-empty (i.e. after at least one successful import).

## DuckDB Integration

`bin/import_all.py <barcode_dir>` loads everything in one process and one connection into the per-barcode `dana.duckdb`:

| Step (in `import_all.py`) | Data |
|--------|------|
| `import_kraken`, `import_krakenreport` | Kraken2 per-read classifications and summary reports with full taxonomy |
| `import_annotations` | Bakta/Prokka gene annotations, stats, read-contig mapping |
| `import_sketch` | Sendsketch per-read GTDB classifications |
| `import_tetra` | Tetranucleotide frequencies + sequence index |
| `import_mapping.py` | minimap2 hits (`map/<file>/<ref>.txt`) into `mapping` |
| `import_protein_hits.py` | HMM / DIAMOND protein-DB hits (see below) |
| `db_schema.py` | Shared schema definitions |

The loader is idempotent and tracks imports via `import_log`. Contig names are prefixed with the batch file ID to prevent collisions across batches. (The single-table scripts `kraken_db.py`, `krakenreport_db.py`, `annotation_db.py`, `sketch_db.py`, `tetra_db.py` are legacy and not called by the pipeline.)

`import_protein_hits.py` (called by `import_all.py`, or standalone with a barcode dir) loads the
protein-DB searches (`hmm/<chunk>.<db>.tsv`, `protdb/<chunk>.<db>.tsv`, written by
`bin/protein_db_search.sh`) into `protein_hits` (one row per protein x target, with a `db`
column), `protein_db_annot` (target -> gene/category/pathway) and `protein_db_searches`
(proteins searched per db/chunk, the normalisation denominator).

## Pollutants tab

`--pollutant_dbs /data/scratch/refdbs/pollutants` searches every annotation protein against
the pollutant marker DBs there (CANT-HYD, HADEG, Hg-MATE, KEGG-xeno, BacMet, PAZy, RDase,
CARD — layout, cutoffs and sources in that directory's `README.md`). `viz/preprocess/preprocess.py`
turns the tables into `pollutants.json` (gene occurrences per read: best target per protein
per DB, one count per read+gene; hgcAB same-read pairs; hgcA/hgcB carrier taxa; mapping
references with a `Pollutant-*` meta.json category such as crAssphage), which the SPA's
**Pollutants** tab groups as database > category > pathway > gene with per-M-protein,
per-M-read or raw abundance on the map and tables.

## Deploy to microscape.app

The viz SPA + preprocessed JSONs are published to the portal via the ingest API. In watch mode, DB_SYNC fires a per-run deploy hook each tick (after preprocess regenerates `<outdir>/viz/*.json*`).

**Batch runs never deploy on their own** — there is no DB_SYNC tick, so
`<outdir>/viz` is never generated. Preprocess and push manually once the
run is done:

```bash
# 1. Generate <outdir>/viz JSONs (sketch distances, DuckDB->JSON, read t-SNE)
viz/preprocess/run_preprocess.sh \
    --input <outdir> --output <outdir>/viz \
    --mapping-refs /data/scratch/refdbs/mapping_refs \
    --env conda-envs/dana-bbmap

# 2. Push to the portal
viz/deploy.sh --preprocess-dir <outdir>/viz \
              --slug grace_20260826 --name "Grace 2026-08-26" --visibility public
```

**One-off manual push** (viz JSONs already present):

```bash
viz/deploy.sh --preprocess-dir <outdir>/viz \
              --slug genice_ci --name "GenIce CI"
```

`deploy.sh` runs `npm install` (first time only), `vite build`, stages `dist/data/`, tarballs the flat tree, and POSTs to `https://microscape.app/api/v1/deploy`.

**Auto-deploy per sync tick:** pass `--deploy_slug` and `--deploy_name` to `run-realtime.sh` and it writes the hook for you:

```bash
./run-realtime.sh --input /data/run1 --outdir /data/output \
    --watch --run_db_integration ... \
    --deploy_slug genice_ci --deploy_name "GenIce CI"
```

Or drop the hook by hand at `<sync dir>/deploy.sh`, where the sync dir is `--store_dir` when one is set (DB_SYNC scans the store dir and builds `viz/` there), otherwise `--outdir`:

```bash
#!/usr/bin/env bash
exec /data/dev/danaSeq/nanopore_live/viz/deploy.sh \
    --preprocess-dir "$1/viz" \
    --slug genice_ci --name "GenIce CI"
```

DB_SYNC detects and invokes any executable `<sync dir>/deploy.sh` at the end of each sync tick, so every live-sequencing refresh rebuilds + re-pushes.

**API key:** read from `$MICROSCAPE_API_KEY` if set, else `~/.config/microscape/api-key`.

## Design Notes

**Kraken2 batching.** In batch mode, FASTAs are grouped per sample (`groupTuple`) so the DB loads once per barcode instead of once per file. In watch mode, files stream individually for live results. Uses `maxForks = 1` since kraken2 loads 50-100 GB into RAM.

**Compiled C tools.** `fastq_filter` (QC) and `tetramer_freqs` (TNF) are compiled C binaries in `bin/`, replacing filtlong and a Python script respectively. Both run without conda environments.

**Watch mode.** Uses `Channel.watchPath()` with a flat glob pattern (Java WatchService limitation -- no recursive `**`).

**Conda environments.** Four isolated environments avoid dependency conflicts:
- **dana-bbmap** -- BBMap (samtools conflicts with R)
- **dana-prokka** -- Prokka (BioPerl pins perl 5.26)
- **dana-bakta** -- Bakta (modern alternative to Prokka)
- **dana-tools** -- Kraken2, HMMER, Python+pandas+DuckDB, Nextflow, OpenJDK

## Common Issues

**No FASTQ files found.** `--input` must point at the directory **containing** `fastq_pass/`. The pipeline gives targeted advice when files aren't found.

**Kraken2 out of memory.** Nextflow serializes Kraken2 via `maxForks = 1`. If running outside Nextflow, never launch multiple Kraken2 instances.

**Conda environment build fails.** Run `./install.sh --check` to diagnose. Ensure mamba/conda is on PATH.

**Watch mode not picking up files.** Java WatchService requires a flat glob pattern (no `**`). The pipeline derives the depth of `fastq_pass/` from the files present at start-up (or watches depths 0-3 if none exist yet); check the `Watch glob:` lines in the log and point `--input` at the run directory or its parent.

## References

- BBTools: Bushnell B, [sourceforge.net/projects/bbmap](https://sourceforge.net/projects/bbmap/)
- Kraken2: Wood DE, Lu J, Langmead B, Genome Biology 2019
- Prokka: Seemann T, Bioinformatics 2014
- HMMER3: Eddy SR, PLoS Computational Biology 2011
- DuckDB: [duckdb.org](https://duckdb.org)
