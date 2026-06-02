# Reference Mapping

Map nanopore reads against arbitrary reference genomes via minimap2 (ONT
preset). Used by the [`AIS` page in the SPA](https://microscape.app/) to
surface invasive species detection, but the module is generic — any
reference works.

The pipeline takes a **directory of references** (`--mapping_refs`) and
fans out across every reference it finds. Each per-barcode dana.duckdb
gains a `mapping` table keyed on `reference` (the reference name), so
downstream tools (the SPA's `preprocess_ais.py`, ad-hoc duckdb queries,
etc.) can pull per-sample stats with a single `WHERE reference = '...'`.

## Quick start

```bash
nextflow run nanopore_live/main.nf \
    --input        /path/to/run \
    --outdir       /data/scratch/nanopore_live/<run-id> \
    --mapping_refs /matika/projects/project_zebra/refs \
    --run_db_integration --danadir nanopore_live/bin \
    -resume
```

Drop a new `refs/<name>/` directory with `<name>.idx` + `meta.json` and
the next pipeline run (or watch-mode tick) picks it up automatically. No
code changes needed.

## Reference directory layout

```
<refs_dir>/
└── <reference_name>/
    ├── <reference_name>.fna          # source fasta (or multifasta)
    ├── <reference_name>.idx          # minimap2 ONT preset index — required
    ├── <reference_name>.contigs.json # linear contig offsets (SPA preprocess)
    ├── meta.json                     # human-readable labels (SPA preprocess)
    └── source.txt                    # free-text provenance (optional)
```

`<reference_name>` is the canonical handle. **Keep it lowercase, alnum,
no spaces** — it becomes:

- the `mapping.reference` value in dana.duckdb;
- the output filename in each barcode dir (`map/<name>.txt`);
- the SPA JSON filename (`ais_<name>.json`);
- the species `id` in the SPA's `AIS_SPECIES` list.

Examples: `zebramussel`, `quaggamussel`, `sealamprey`, `roundgoby`.

## Adding a new reference

The end-to-end workflow has five short steps. Examples use sea lamprey;
substitute `REFNAME` and `ACC` for your target.

### 1. Acquire the fasta

```bash
REFNAME=sealamprey
ACC=GCF_048934315.1     # RefSeq / GenBank accession
REFDIR=/matika/projects/project_zebra/refs/$REFNAME
mkdir -p "$REFDIR" && cd "$REFDIR"

# NCBI datasets API → zip → extract the genomic fasta
curl -sS --fail -o ncbi.zip \
  "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/$ACC/download?include_annotation_type=GENOME_FASTA"
FNA_PATH=$(unzip -l ncbi.zip | awk '/_genomic\.fna$/ {print $4}' | head -1)
unzip -p ncbi.zip "$FNA_PATH" > "$REFNAME.fna"
rm ncbi.zip
```

For species with no whole-genome assembly (e.g. some AIS invertebrates),
fall back to a multifasta of all GenBank records under the taxon:

```bash
TAXID=77747
ES=$(curl -sS "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=nuccore&term=txid${TAXID}[Organism]&usehistory=y&retmax=0")
WENV=$(echo "$ES" | grep -oP '<WebEnv>\K[^<]+')
QK=$(echo  "$ES" | grep -oP '<QueryKey>\K[^<]+')
curl -sS "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&WebEnv=$WENV&query_key=$QK&rettype=fasta&retmode=text" \
  > "$REFNAME.fna"
```

Detection sensitivity is then limited to whatever regions are covered by
those markers — typically COI / 18S / mitochondrion. Note this in
`meta.json.notes` so future-you remembers.

### 2. Build the minimap2 index

The pipeline requires a pre-built `.idx` (memory-mapping a 4 GB index per
task is much faster than rebuilding it 600× per pipeline run):

```bash
minimap2 -d "$REFDIR/$REFNAME.idx" "$REFDIR/$REFNAME.fna"
```

If you want to use a non-ONT preset, build the index with that preset's
flags here — the pipeline only invokes `minimap2 -ax map-ont` against the
.idx, but minimap2 picks up most index-time parameters from the .idx
itself.

### 3. Build the contig offset map

The SPA's per-sample position histogram needs to know how to linearise
hits across multiple contigs. `build_contig_index.py` writes a small
JSON sidecar:

```bash
python3 /matika/projects/project_zebra/build_contig_index.py "$REFDIR/$REFNAME.fna"
# → writes $REFDIR/$REFNAME.contigs.json
```

Resulting shape:
```json
{
  "total_length": 1808348452,
  "contigs": [
    {"name": "NC_133667.1", "length": 37671560, "offset": 0},
    {"name": "NC_133668.1", "length": 29380143, "offset": 37671560},
    ...
  ]
}
```

### 4. Write the helper file `meta.json`

The pipeline itself doesn't need `meta.json` (it uses the directory name
as the canonical ID). The SPA preprocess reads it to label markers /
species cards.

```json
{
  "id": "sealamprey",
  "species": "Petromyzon marinus",
  "common_name": "Sea lamprey",
  "category": "AIS",
  "accession": "GCF_048934315.1",
  "assembly_name": "UKy_Petmar_22M1.pri1.0",
  "assembly_level": "Chromosome",
  "source": "NCBI Datasets API",
  "added": "2026-05-27",
  "notes": "Free-form. Close native sister taxa, recommended HQ cutoffs, etc."
}
```

| Field | Required | Used by | Notes |
|---|---|---|---|
| `id` | yes | SPA | Must match directory + `.idx` basename |
| `species` | yes | SPA | Latin binomial; rendered italicised |
| `common_name` | yes | SPA | UI label |
| `category` | no | SPA | Future grouping (e.g. `AIS`, `pathogen`, `host`) |
| `accession` | recommended | docs | NCBI/RefSeq accession |
| `assembly_name`, `assembly_level` | recommended | docs | Provenance |
| `source` | recommended | docs | Where the fasta came from |
| `added` | recommended | docs | ISO date |
| `notes` | recommended | docs | Cross-mapping caveats, recommended HQ cutoff, etc. |

### 5. Run the pipeline

Once `<refs_dir>/<name>/<name>.idx` exists, point the pipeline at
`<refs_dir>` and re-run (or it'll pick up automatically on the next
watch-mode tick):

```bash
nextflow run nanopore_live/main.nf \
    --input        /path/to/run \
    --outdir       /data/scratch/nanopore_live/<run-id> \
    --mapping_refs <refs_dir> \
    --run_db_integration --danadir nanopore_live/bin \
    -resume
```

Outputs land at:
- **Per-barcode** : `<outdir>/<FC>/<BC>/map/<refname>.txt` (filtered SAM-ish: mapq≥1, aligned_len≥10)
- **DuckDB** : `<outdir>/<FC>/<BC>/dana.duckdb` → table `mapping` (loaded by `import_all.py → import_mapping`)

### 6. Surface in the SPA

Edit `nanopore_live/viz/src/views/AISView.svelte`:

```js
const AIS_SPECIES = [
  { id: 'zebramussel', label: 'Zebra mussel', latin: 'Dreissena polymorpha' },
  { id: 'sealamprey',  label: 'Sea lamprey',  latin: 'Petromyzon marinus' },
  // ... add new id here. Must match meta.json `id` exactly.
];
```

Then run `preprocess_ais.py` to aggregate the per-barcode duckdb rows
into the SPA-shaped JSON:

```bash
DPY=nanopore_live/conda-envs/dana-tools/bin/python
$DPY preprocess_ais.py \
    --reference $REFNAME \
    --species   "$(jq -r .species     refs/$REFNAME/meta.json)" \
    --common    "$(jq -r .common_name  refs/$REFNAME/meta.json)" \
    --contigs   refs/$REFNAME/$REFNAME.contigs.json \
    --output    <preprocess_dir>/ais_$REFNAME.json
```

Then `viz/deploy.sh` to push to microscape.app.

## HAB compound references

The same `mapping_refs` plumbing powers the SPA's `/habs` view for
cyanobacterial HAB compound BGC detection. Differences vs. AIS:

- Each ref is a **biosynthetic gene cluster** (full BGC fasta with
  intergenic regions, not just CDS) named after the chemical compound:
  `microcystin`, `cylindrospermopsin`, `saxitoxin`, `anatoxin_a`,
  `nodularin`.
- The fasta concatenates the BGC from **multiple producer taxa** per
  compound (e.g. `microcystin.fna` contains the mcy cluster from
  *Microcystis*, *Planktothrix*, *Anabaena/Dolichospermum*, …) so a hit
  reflects bloom potential regardless of which producer is present. The
  per-producer offsets are surfaced as ticks on the position histogram.
- `meta.json` uses `category: "HAB"` and adds `compound` (chemical name)
  + `class` (chemistry / syndrome class, e.g. `cyclic heptapeptide`) in
  place of `species` / `common_name`.
- Preprocess: run `preprocess_habs.py` (sister to `preprocess_ais.py`)
  to emit `hab_<id>.json`:

  ```bash
  $DPY preprocess_habs.py \
      --reference $REFNAME \
      --compound  "$(jq -r .compound refs/$REFNAME/meta.json)" \
      --class     "$(jq -r .class    refs/$REFNAME/meta.json)" \
      --contigs   refs/$REFNAME/$REFNAME.contigs.json \
      --output    <preprocess_dir>/hab_$REFNAME.json
  ```

- SPA: `nanopore_live/viz/src/views/HABsView.svelte` mirrors AISView
  (same Leaflet + identity + position histograms, amber accent instead
  of emerald). Add new compounds to `HAB_COMPOUNDS` there; `id` must
  match the ref dirname / `.idx` basename / `hab_<id>.json` filename.

## Standalone scripts (no Nextflow)

For one-off batches outside the pipeline:

| Script | What it does |
|---|---|
| `/matika/projects/project_zebra/run-nanopore-mapping.sh` | Walks `out_dana_bc/<FC>/<BC>/fq/`, runs minimap2, writes `map/<refname>.txt` to scratch. Args: `REFERENCE_INDEX_PATH`. |
| `/matika/projects/project_zebra/load_mapping_to_duckdb.py` | Walks the scratch tree, appends to each barcode's `dana.duckdb`. Same `mapping` table, same `import_log` key. Idempotent. |
| `/matika/projects/project_zebra/build_contig_index.py` | Writes `<name>.contigs.json` from `<name>.fna`. |
| `/matika/projects/project_zebra/preprocess_ais.py` | Aggregates duckdb → `ais_<name>.json` for the SPA. |

The pipeline integration uses the same `mapping` schema and `import_log`
key (`map/<refname>.txt`), so files imported by either path round-trip
cleanly.

## DuckDB schema

```sql
CREATE TABLE mapping (
    reference    TEXT,      -- e.g. 'zebramussel'
    qname        TEXT,      -- read id
    flag         INTEGER,
    rname        TEXT,      -- target sequence
    pos          INTEGER,
    mapq         INTEGER,
    cigar        TEXT,
    nm           INTEGER,   -- NM tag (edit distance)
    as_score     INTEGER,   -- AS tag
    de           REAL,      -- minimap2 gap-compressed divergence (de:f:)
    identity_pct REAL,      -- (1 - de) * 100
    aligned_len  INTEGER
);
```

Useful aggregations:

```sql
-- Per-reference summary across all barcodes (run inside one dana.duckdb)
SELECT reference,
       COUNT(*)                 AS hits,
       ROUND(AVG(identity_pct), 2) AS mean_id,
       SUM(CASE WHEN identity_pct > 90 THEN 1 ELSE 0 END) AS hq_hits,
       ROUND(AVG(mapq), 1)      AS mean_mapq
FROM mapping
GROUP BY reference;
```
