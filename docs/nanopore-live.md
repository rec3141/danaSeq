# Real-Time Nanopore Pipeline

Real-time metagenomic analysis for Oxford Nanopore sequencing data. Processes reads as they stream from MinKNOW, providing live taxonomic classification, gene annotation, and functional profiling.

## Quick Start

```bash
cd nanopore_live
./install.sh && ./install.sh --check

./run-realtime.sh --input /path/to/nanopore/run --outdir /path/to/output \
    --run_kraken --kraken_db /path/to/krakendb \
    --annotator bakta --bakta_db /path/to/bakta/db-light --run_sketch --run_tetra
```

## Pipeline Stages

| Step | Process | Tool | Description |
|------|---------|------|-------------|
| 1 | `VALIDATE_FASTQ` | gzip + `fastq_gz_salvage.py` | Gzip integrity check; damaged files salvaged record by record |
| 2 | `QC_FASTQ_FILTER` | fastq_filter (C++) | Length + quality filtering (streaming, Filtlong-style scoring) |
| 3 | `CONVERT_TO_FASTA` | awk | Header cleanup and format conversion |
| 4a | `KRAKEN2_CLASSIFY` | Kraken2 | Taxonomic classification (batched per sample) |
| 4b | `PROKKA_ANNOTATE` | Prokka | ORF prediction + functional annotation |
| 4c | `BAKTA_CDS` / `BAKTA_FULL` | Bakta | CDS-only or full annotation (default annotator) |
| 4d | `HMM_SEARCH` | HMMER3 | Profile HMM search against user-supplied / pollutant databases |
| 4e | `DIAMOND_DB_SEARCH` | DIAMOND | blastp of annotation proteins against `.dmnd` databases |
| 4f | `SENDSKETCH` | BBTools | Per-read GTDB classification via a sendsketch server |
| 4g | `TETRAMER_FREQ` | tetramer_freqs (C) | Tetranucleotide composition profiles |
| 4h | `MAP_REFERENCE` | minimap2 (ONT) | Per-reference alignment for AIS / custom genomes — see [mapping-references.md](mapping-references.md) |
| 5 | `DB_INTEGRATION` | Python + DuckDB | Load all results into DuckDB (batch mode) |
| 6 | `DB_SYNC` | Python + DuckDB | Periodic sync, dashboard preprocess and deploy hook during watch mode |
| 7 | `CLEANUP` | bash | Delete imported source files after DB import |

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--input` | (required) | Directory **containing** `fastq_pass/` |
| `--outdir` | `results` | Output directory |
| `--run_kraken` | `false` | Enable Kraken2 classification |
| `--kraken_db` | (required if kraken) | Path to Kraken2 database |
| `--annotator` | `bakta` | Gene annotator: `prokka`, `bakta`, or `none` |
| `--bakta_db` | `/data/scratch/refdbs/bakta/db-light` | Path to Bakta database (site-specific default) |
| `--bakta_full` | `false` | Run full Bakta annotation (ncRNA/tRNA/CRISPR) |
| `--hmm_databases` | (skip) | Comma-delimited paths to HMM files |
| `--hmm_cut` | `tc` | hmmsearch cutoff (`tc`/`ga`/`nc`/E-value); per-DB `<name>.opts` overrides |
| `--protein_dbs` | (skip) | Comma-delimited DIAMOND `.dmnd` paths |
| `--protein_db_opts` | `id=50,qcov=50,scov=50,evalue=1e-5` | DIAMOND cutoffs |
| `--pollutant_dbs` | (skip) | Directory of `<name>/<name>.hmm|.dmnd` pollutant databases |
| `--mapping_refs` | (skip) | Directory of reference subdirs for minimap2 mapping (AIS / custom) — see [mapping-references.md](mapping-references.md) |
| `--run_sketch` | `false` | Enable per-read GTDB sendsketch classification |
| `--sendsketch_address` | `$DANA_SENDSKETCH_ADDRESS`, else `http://localhost:3068/sketch` | Sendsketch server URL |
| `--run_tetra` | `false` | Enable tetranucleotide frequency |
| `--watch` | `false` | Monitor for new files during live sequencing (watch globs are auto-detected from the depth of `fastq_pass/` below `--input`) |
| `--db_sync_minutes` | `10` | DuckDB sync interval in watch mode |
| `--run_db_integration` | `false` | Load results into DuckDB |
| `--cleanup` | `false` | Delete imported files after DB import |
| `--metadata` | (none) | Sample metadata TSV passed to the dashboard preprocess |
| `--min_readlen` | `1500` | Minimum read length after filtering |
| `--keep_percent` | `80` | Percent of reads to keep (by quality) |
| `--min_file_size` | `1000000` | Minimum FASTQ file size in bytes (1 MB) |
| `--store_dir` | (none) | Persistent cache directory (storeDir) |

## Outputs

```
results/
├── FLOWCELL/
│   ├── barcode01/
│   │   ├── fa/           FASTA (quality-filtered sequences)
│   │   ├── kraken/       *.tsv (per-read), *.report (summary)
│   │   ├── bakta/ or prokka/   <id>/<id>.gff3|.faa|.tsv (Bakta) or <id>/PROKKA_* (Prokka)
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

### Read map (t-SNE of every read)

In watch mode, DB_SYNC runs `viz/preprocess/compute_read_tsne.py` whenever
`tetra_data` grows. With [fastsne](https://github.com/rec3141/fastsne) installed
(the container ships the CPU build; `install.sh` builds it into `dana-tools`),
every read is mapped, not a sample:

- each sample's `tetra_data` is exported once to `viz/read_map/fsb/`
  (re-exported only when its `dana.duckdb` changes)
- a reference t-SNE is built from a uniform random sample of up to `--ref-max`
  reads (10M on CPU, 20M on a GPU) and rebuilt when the read count has doubled
  (`--rebuild-factor`)
- the other reads are projected onto it, only for samples whose data changed
- `viz/read_explorer.json` keeps a `--max-reads` sample for the dashboard, with
  coordinates from the full map
- `viz/read_map/viewer.html?data=reads.json` shows every read (WebGL, paged),
  coloured by kraken genus, sample or distance to the reference; projected reads
  can be hidden

fastsne_cpu uses half the available cores by default (`--threads`), leaving room
for basecalling and classification. Without fastsne, or with `--engine sklearn`,
the sampled sklearn t-SNE runs as before.

## Profiles

| Profile | Use case |
|---------|----------|
| `standard` | Local execution |
| `test` | Lower size/length thresholds (no test data are bundled) |
| `shipboard` | Production workstation: 32 CPUs, 256 GB RAM, 100 GB Kraken2 |
| `local_tools` | Legacy: conda disabled, tools expected on `PATH` |

## Resource Requirements

Defaults from `nextflow.config` (batch / watch mode):

| Component | CPUs | RAM | Notes |
|-----------|------|-----|-------|
| Kraken2 | 4 / 1 | 10 GB (100 GB in `shipboard`) | Depends on database size; serialized via `maxForks=1` |
| Bakta (`BAKTA_CDS`) | all cores / 1 | 16 GB / 4 GB | Batch: one task at a time; watch: one task per core |
| minimap2 (`MAP_REFERENCE`) | 1 | 4 GB | At most 16 / 4 concurrent tasks |
| Other `process_medium` tasks (Prokka, HMMER3, sendsketch, DB sync) | 1 | 4 GB | |
| Other `process_low` tasks (QC, DIAMOND, tetramers) | 1 | 2 GB | |

## Watch Mode

Monitor a directory for new FASTQ files during active sequencing:

```bash
./run-realtime.sh --input /path/to/runs --outdir /path/to/output \
    --watch --db_sync_minutes 10 \
    --run_kraken --kraken_db /path/to/db \
    --annotator bakta --run_db_integration
```

`DB_SYNC` runs as a long-lived process that periodically loads new results into DuckDB. The Python loader (`bin/import_all.py`) is idempotent and tracks imports via `import_log`. Watch vs batch mode is locked per `--outdir` (recorded in `<outdir>/.pipeline_mode`).

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

`bin/import_all.py <barcode_dir>` loads all results into the per-barcode `dana.duckdb` in one process:

| Step | Data |
|--------|------|
| `import_kraken`, `import_krakenreport` | Kraken2 per-read classifications and summary reports |
| `import_annotations` | Bakta/Prokka gene annotations, stats, read-contig mapping |
| `import_sketch` | Sendsketch per-read GTDB classifications |
| `import_tetra` | Tetranucleotide frequencies + sequence index |
| `import_mapping.py` | minimap2 hits into the `mapping` table |
| `import_protein_hits.py` | HMM / DIAMOND protein-DB hits |
| `db_schema.py` | Shared schema definitions |

The loader is idempotent and tracks imports via `import_log`. The legacy R scripts in `bin/` (`47_merge_db.r` etc.) are not called by the pipeline.

## Input

Oxford Nanopore directory structure with multiplexed barcodes:

```
input_dir/fastq_pass/
├── barcode01/*.fastq.gz
├── barcode02/*.fastq.gz
└── ...
```

## Design Notes

- **Kraken2 batching.** In batch mode, FASTAs are grouped per sample (`groupTuple`) so the DB loads once per barcode. In watch mode, files stream individually. Uses `maxForks = 1` since Kraken2 loads 50-100 GB into RAM.
- **Compiled C tools.** `fastq_filter` (QC) and `tetramer_freqs` (TNF) are compiled C binaries in `bin/`, replacing filtlong and a Python script.
- **Watch mode.** Uses `Channel.watchPath()` with a flat glob pattern (Java WatchService limitation -- no recursive `**`).
- **Conda environments.** Four isolated environments avoid dependency conflicts: dana-bbmap, dana-prokka, dana-bakta, dana-tools.
