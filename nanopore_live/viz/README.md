# nanopore_live viz

Single-page dashboard for a `nanopore_live` run — sample overview, taxonomy,
function, reads, t-SNE. Svelte + Vite. Consumes a directory of preprocessed
JSON files produced by `preprocess/run_preprocess.sh`, and ships as a flat
`dist/` tree that can be dropped under any path on a static host.

## Data shape the SPA expects

All fetches are relative (`./data/*.json[.gz]`), so the SPA works at any
deploy path (`/`, `/chesterfield/`, `/genice_ci/`, …) without rebuilding.
At runtime it probes `.json.gz` first and falls back to plain `.json` if
absent:

| File | Source |
|---|---|
| `data/overview.json` | aggregate counts (flowcells, samples, reads) |
| `data/samples.json` | per-barcode rows |
| `data/metadata.json` | optional external sample-metadata TSV, parsed |
| `data/sample_taxonomy.json` | per-sample taxonomy counts |
| `data/sample_function.json` | per-sample KO/function counts |
| `data/sample_tsne.json` | sample t-SNE coordinates |
| `data/taxonomy_sunburst.json` | hierarchical taxonomy counts |
| `data/read_explorer.json` | per-read detail (large; lazy-loaded) |

`preprocess/run_preprocess.sh` writes these to an arbitrary output
directory (pass `--output`). Its outputs are **flat** (no `data/`
wrapper). The `data/` subdir convention is a deploy-time reshape — see
below.

## Dev loop

Serve the built SPA against a live preprocess output via Vite's
`VIZ_DATA_DIR` middleware (maps `/data/*` → files in that dir):

```bash
cd nanopore_live/viz
npm install

# one-shot dev
VIZ_DATA_DIR=/data/scratch/nanopore_live/genice_ci npm run dev

# preview a production build
npm run build
VIZ_DATA_DIR=/data/scratch/nanopore_live/genice_ci npm run serve
```

## Production deploy → microscape.app

`deploy.sh` in this directory builds the SPA, stages the preprocess output
under `dist/data/`, tars it, and pushes to
[microscape.app](https://microscape.app)'s bearer-authenticated ingest
endpoint (`POST /api/v1/deploy`). One tarball = one run.

### One-time: API key

1. Sign in to https://microscape.app as a lab-admin.
2. Go to **Admin → API keys**, click **Mint key**, name it after the
   pipeline (`"danaseq nanopore_live"` etc.).
3. Copy the plaintext (shown once) and save it as
   `~/.config/microscape/api-key` with mode 600:

   ```bash
   mkdir -p ~/.config/microscape
   install -m 600 /dev/null ~/.config/microscape/api-key
   cat > ~/.config/microscape/api-key     # paste key, Ctrl-D
   ```

   Or put it in `$MICROSCAPE_API_KEY` via your shell profile / CI secrets
   store / systemd `EnvironmentFile`.

### Per-run deploy

```bash
./deploy.sh \
    --preprocess-dir /data/scratch/nanopore_live/genice_ci \
    --slug genice_ci \
    --name "GenIce CI"
```

Slug rules: lowercase alnum with `-` or `_`, up to 64 chars. It becomes the
run's URL — `https://microscape.app/<slug>/` — and must not collide with
app-reserved paths (the server will 400 a collision). The slug should
generally match the preprocess output's directory name so operators can
find runs.

See `./deploy.sh --help` for flags. Notables:

- `--skip-build` reuses an existing `dist/` — fast iteration during CI
  plumbing.
- `--dry-run` builds + tarballs but skips the POST, and leaves the tarball
  on disk for inspection.
- `--endpoint` / `--pipeline` override the defaults
  (`https://microscape.app/api/v1/deploy`, `danaseq-nanopore-live`).

Re-deploying the same slug is idempotent — the microscape-app server
upserts the run row and `rsync -a --delete`s the new tree into place, so
stale files from a previous deploy are cleaned up.

### Wiring into the pipeline

After each preprocess tick, have the pipeline execute a per-run hook at
`<outdir>/deploy.sh`. Minimal template:

```bash
#!/usr/bin/env bash
# <outdir>/deploy.sh — invoked by DB_SYNC after each preprocess pass
exec /path/to/nanopore_live/viz/deploy.sh \
    --preprocess-dir "$1/viz" \
    --slug "$(basename "$1")" \
    --name "$(basename "$1")"
```

If the preprocess output lands as mixed `.json`/`.json.gz` (the default),
the SPA's `.gz`-first probe will 404 on files that ship only as plain
`.json` and fall back to the plain fetch. Functional but noisy in the
console. To silence: gzip everything in the output dir before deploy,

```bash
find <preprocess-dir> -type f \( -name '*.json' -o -name '*.tsv' \) \
    ! -name '*.gz' -exec gzip -f {} +
```

or bake that into `run_preprocess.sh` as a final step.

## What deploy.sh does end-to-end

1. `npm install` (if `node_modules/` missing) and `npm run build` →
   `dist/index.html` + `dist/assets/*` (relative-URL bundle thanks to
   `base: './'` in `vite.config.js`).
2. `rm -rf dist/data && mkdir dist/data` and copy every
   `*.json*` from `--preprocess-dir` into `dist/data/`.
3. `tar -czf $TMP_TGZ -C dist .` (flat root: `index.html`, `assets/`,
   `data/` directly under it — no wrapper dir).
4. `curl --fail -X POST` the tarball to
   `https://microscape.app/api/v1/deploy` with
   `Authorization: Bearer $API_KEY` and the metadata headers
   (`X-Microscape-Slug`, `X-Microscape-Pipeline`, `X-Microscape-Name`).
5. The microscape-app server authenticates the key, extracts the tarball
   into an isolated staging dir, `rsync --delete`s it into
   `RUNS_ROOT/<slug>/`, upserts the run row, and returns the canonical
   URL.

See the
[microscape-app README](https://github.com/rec3141/microscape-app#pipeline-ingest-api)
for the server side of the ingest contract (request/response shape,
visibility rules, `.gz` handling).
