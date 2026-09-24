"""Chunked, columnar writers for the per-contig viz data.

A co-assembly can hold millions of contigs. A single JSON file for the contig
explorer, per-sample depths or gene features then decompresses to more than a
browser can hold as one string (~512 MB), and the dashboard hangs. These
writers split the data into files that the SPA fetches only as needed:

  contig_explorer_manifest.json         columns, vocabularies, extents, chunk list
  contig_explorer.part-NNN.json.gz      one chunk of contigs, columnar
  contig_sample_depths_manifest.json    samples, maxDepth, per-chunk counts
  contig_sample_depths.part-NNN.json.gz sparse per-sample depths, aligned to
                                        the explorer chunk with the same NNN
  genes_manifest.json                   shard list with first/last sort keys
  genes.part-NNN.json.gz                {contig_id: [feature, ...]} for a
                                        contiguous run of contigs

Every file sits directly in the output directory: the deploy tool publishes
only files at the top level of viz/data.

Contigs are ordered by (length descending, id ascending) in all three sets, so
the contigs at or above any length threshold are a prefix of the chunk list,
and the SPA fetches only the chunks its minimum-length filter needs. The same
key locates a contig's gene shard.

t-SNE and UMAP coordinates are columns of the explorer chunks rather than
separate files. They describe the same rows in the same order, the SPA needs
them together with the other columns, and keeping them beside the contig ids
lets a --skip-tsne / --skip-umap rerun read them back by id however the
chunking changes.

Explorer chunk layout:
  {"index": k, "offset": first global row, "count": n,
   "ids": [delta-encoded integers] | ["id", ...],
   "columns": {name: [value per row], ...}}
Numeric columns hold numbers (null where absent, e.g. a contig with no
embedding). Coded columns hold integers indexing manifest["vocab"][name]; code
0 is always "". When every id is <prefix><zero-padded integer> with one prefix
and width, "ids" holds the integers delta-encoded (first value relative to 0)
and manifest["id_encoding"] = {"prefix", "width"}; otherwise ids are listed.

Sample-depth chunk layout (only nonzero depths are stored):
  {"index": k, "offset": ..., "count": n,
   "idx": [[delta-encoded row-in-chunk, ...] per sample],
   "val": [[depth, ...] per sample]}
"""

import glob
import gzip
import json
import os
import re

import numpy as np

EXPLORER = 'contig_explorer'
DEPTHS = 'contig_sample_depths'
GENES = 'genes'
EMBEDDINGS = ('tsne', 'umap')

DEFAULT_CONTIG_CHUNK_SIZE = 100_000
DEFAULT_GENE_SHARD_SIZE = 250_000

FORMAT_VERSION = 1
SORT_ORDER = 'length_desc_id_asc'

_ID_RE = re.compile(r'^(.*?)(\d+)$')


# ---------------------------------------------------------------------------
# File helpers
# ---------------------------------------------------------------------------

def _dumps(obj):
    return json.dumps(obj, separators=(',', ':'), allow_nan=False)


def write_manifest(path, obj):
    """Manifests are small: write .json and .json.gz like the other viz files."""
    text = _dumps(obj)
    with open(path, 'w') as f:
        f.write(text)
    with gzip.open(path + '.gz', 'wt', compresslevel=6) as f:
        f.write(text)


def write_part(path, obj):
    """Chunk files are written gzip-only (<path>.gz); the SPA fetches .gz first."""
    tmp = path + '.gz.tmp'
    with gzip.open(tmp, 'wt', compresslevel=6) as f:
        f.write(_dumps(obj))
    os.replace(tmp, path + '.gz')


def read_json_any(path):
    """Read <path>.gz if present, else <path>; None if neither exists."""
    if os.path.exists(path + '.gz'):
        with gzip.open(path + '.gz', 'rt') as f:
            return json.load(f)
    if os.path.exists(path):
        with open(path) as f:
            return json.load(f)
    return None


def part_name(prefix, index):
    return f'{prefix}.part-{index:03d}.json'


def _remove(path):
    for p in (path, path + '.gz'):
        if os.path.exists(p):
            os.remove(p)


def remove_stale_parts(output_dir, prefix, keep):
    """Delete <prefix>.part-*.json[.gz] files whose index is not below `keep`."""
    pat = re.compile(re.escape(prefix) + r'\.part-(\d+)\.json(\.gz)?$')
    for p in glob.glob(os.path.join(output_dir, prefix + '.part-*')):
        m = pat.search(os.path.basename(p))
        if m and int(m.group(1)) >= keep:
            os.remove(p)


def remove_chunked(output_dir, prefix):
    """Remove a chunked data set (manifest + parts)."""
    _remove(os.path.join(output_dir, f'{prefix}_manifest.json'))
    remove_stale_parts(output_dir, prefix, 0)


def remove_legacy(output_dir, *names):
    """Remove single-file outputs superseded by the chunked layout.

    The SPA prefers the manifest, so a leftover legacy file is never read but
    would still be staged by deploy.sh and can be larger than a browser can load.
    """
    for name in names:
        _remove(os.path.join(output_dir, name))


# ---------------------------------------------------------------------------
# Ordering and id encoding
# ---------------------------------------------------------------------------

def sort_key(length, cid):
    return (-int(length), cid)


def id_encoding(ids):
    """Return {"prefix", "width"} when every id is prefix + fixed-width digits."""
    if not ids:
        return None
    m = _ID_RE.match(ids[0])
    if not m:
        return None
    prefix, width = m.group(1), len(m.group(2))
    for cid in ids:
        m = _ID_RE.match(cid)
        if not m or m.group(1) != prefix or len(m.group(2)) != width:
            return None
    return {'prefix': prefix, 'width': width}


def encode_ids(ids, enc):
    if enc is None:
        return list(ids)
    plen = len(enc['prefix'])
    out, prev = [], 0
    for cid in ids:
        n = int(cid[plen:])
        out.append(n - prev)
        prev = n
    return out


def decode_ids(encoded, enc):
    if enc is None:
        return list(encoded)
    out, n = [], 0
    for d in encoded:
        n += d
        out.append(f"{enc['prefix']}{n:0{enc['width']}d}")
    return out


def _is_number(v):
    return isinstance(v, (int, float)) and not isinstance(v, bool)


def _num(v):
    if v is None:
        return None
    if isinstance(v, float) and v != v:
        return None
    return v


# ---------------------------------------------------------------------------
# Contig explorer (+ embeddings) and per-sample depths
# ---------------------------------------------------------------------------

def write_contig_explorer(output_dir, explorer, embeddings=None, sample_depths=None,
                          chunk_size=DEFAULT_CONTIG_CHUNK_SIZE, min_length=0):
    """Write the chunked contig explorer and, if given, per-sample depths.

    explorer:      {"contigs": [record, ...], other metadata...} as built by
                   build_contig_explorer (records: id, length, numeric and
                   string fields).
    embeddings:    {"tsne": {id: [x, y]} or None, "umap": ...}
    sample_depths: {"samples": [...], "depths": {id: [...]}, "maxDepth": x} or None
    min_length:    drop contigs shorter than this (0 keeps every contig).

    Returns the explorer manifest.
    """
    chunk_size = max(1, int(chunk_size))
    embeddings = embeddings or {}
    records = explorer.get('contigs', [])
    if min_length > 0:
        records = [r for r in records if r['length'] >= min_length]
    records = sorted(records, key=lambda r: sort_key(r['length'], r['id']))
    n = len(records)

    # Column classification: the first record fixes the column order; a column
    # is numeric when every value is a number (None allowed), else coded.
    base_cols = [k for k in (records[0].keys() if records else []) if k != 'id']
    numeric, coded = [], []
    for col in base_cols:
        if all(_is_number(r.get(col)) or r.get(col) is None for r in records):
            numeric.append(col)
        else:
            coded.append(col)

    has_emb = {}
    for method in EMBEDDINGS:
        emb = embeddings.get(method)
        has_emb[method] = bool(emb) and any(r['id'] in emb for r in records)
        if has_emb[method]:
            numeric += [f'{method}_x', f'{method}_y']

    # Materialise columns once in sorted order.
    ids = [r['id'] for r in records]
    cols = {}
    for col in base_cols:
        cols[col] = [r.get(col) for r in records]
    for method in EMBEDDINGS:
        if has_emb[method]:
            emb = embeddings[method]
            xs, ys = [], []
            for cid in ids:
                xy = emb.get(cid)
                xs.append(_num(xy[0]) if xy else None)
                ys.append(_num(xy[1]) if xy else None)
            cols[f'{method}_x'], cols[f'{method}_y'] = xs, ys

    vocab, codes = {}, {}
    for col in coded:
        vals = ['' if v is None else str(v) for v in cols[col]]
        voc = [''] + sorted(set(vals) - {''})
        lookup = {v: i for i, v in enumerate(voc)}
        vocab[col] = voc
        codes[col] = [lookup[v] for v in vals]

    # Minimum contig length per bin, so a view can load just enough chunks to
    # see every contig of one bin.
    bin_min_length = {}
    lengths = cols.get('length', [0] * n)
    for col in coded:
        if col == 'bin' or col.endswith('_bin'):
            mins = {}
            voc = vocab[col]
            for code, ln in zip(codes[col], lengths):
                if code:
                    name = voc[code]
                    if name not in mins or ln < mins[name]:
                        mins[name] = ln
            if mins:
                bin_min_length[col] = mins

    extents = {}
    for col in numeric:
        vals = [v for v in map(_num, cols[col]) if v is not None]
        if vals:
            extents[col] = [min(vals), max(vals)]

    enc = id_encoding(ids)
    chunks = []
    n_chunks = (n + chunk_size - 1) // chunk_size
    for k in range(n_chunks):
        lo, hi = k * chunk_size, min(n, (k + 1) * chunk_size)
        body = {
            'index': k, 'offset': lo, 'count': hi - lo,
            'ids': encode_ids(ids[lo:hi], enc),
            'columns': {},
        }
        for col in base_cols + [c for c in numeric if c not in base_cols]:
            body['columns'][col] = codes[col][lo:hi] if col in codes else [_num(v) for v in cols[col][lo:hi]]
        fname = part_name(EXPLORER, k)
        write_part(os.path.join(output_dir, fname), body)
        chunks.append({
            'index': k, 'file': fname, 'offset': lo, 'count': hi - lo,
            'max_length': lengths[lo], 'min_length': lengths[hi - 1],
            'first': [lengths[lo], ids[lo]], 'last': [lengths[hi - 1], ids[hi - 1]],
        })
    remove_stale_parts(output_dir, EXPLORER, n_chunks)

    manifest = {
        'format': 'danaseq-contig-explorer', 'version': FORMAT_VERSION,
        'sort': SORT_ORDER,
        'n_contigs': n,
        'chunk_size': chunk_size,
        'min_length_cutoff': int(min_length),
        'id_encoding': enc,
        'column_order': base_cols + [c for c in numeric if c not in base_cols],
        'numeric_columns': numeric,
        'coded_columns': coded,
        'vocab': vocab,
        'bin_min_length': bin_min_length,
        'extents': extents,
        'has_tsne': has_emb['tsne'],
        'has_umap': has_emb['umap'],
        'chunks': chunks,
    }
    for key, val in explorer.items():
        if key not in ('contigs', 'has_tsne', 'has_umap'):
            manifest[key] = val

    if sample_depths is not None:
        write_sample_depths(output_dir, ids, chunks, sample_depths)
    else:
        remove_chunked(output_dir, DEPTHS)

    # Manifest last: a reader never sees a manifest naming parts not yet written.
    write_manifest(os.path.join(output_dir, f'{EXPLORER}_manifest.json'), manifest)
    remove_legacy(output_dir, f'{EXPLORER}.json', 'contig_tsne.json', 'contig_umap.json')
    return manifest


def write_sample_depths(output_dir, ids, chunks, sample_depths):
    """Write sparse per-sample depth chunks aligned to the explorer chunks."""
    samples = list(sample_depths['samples'])
    depths = sample_depths['depths']
    n_s = len(samples)
    zero = [0.0] * n_s
    out_chunks = []
    max_depth = 0.0
    for ch in chunks:
        lo, hi = ch['offset'], ch['offset'] + ch['count']
        mat = np.array([depths.get(cid) or zero for cid in ids[lo:hi]], dtype=np.float64)
        if mat.size:
            max_depth = max(max_depth, float(mat.max()))
        idx_lists, val_lists, counts = [], [], []
        for s in range(n_s):
            col = mat[:, s] if mat.size else np.zeros(0)
            nz = np.flatnonzero(col)
            deltas = np.diff(nz, prepend=0) if nz.size else nz
            idx_lists.append(deltas.tolist())
            val_lists.append(col[nz].tolist())
            counts.append(int(nz.size))
        fname = part_name(DEPTHS, ch['index'])
        write_part(os.path.join(output_dir, fname), {
            'index': ch['index'], 'offset': lo, 'count': hi - lo,
            'idx': idx_lists, 'val': val_lists,
        })
        out_chunks.append({'index': ch['index'], 'file': fname, 'offset': lo,
                           'count': hi - lo, 'nonzero': counts})
    remove_stale_parts(output_dir, DEPTHS, len(chunks))
    manifest = {
        'format': 'danaseq-contig-sample-depths', 'version': FORMAT_VERSION,
        'samples': samples,
        'maxDepth': sample_depths.get('maxDepth', round(max_depth, 4)),
        'aligned_to': f'{EXPLORER}_manifest.json',
        'chunks': out_chunks,
    }
    write_manifest(os.path.join(output_dir, f'{DEPTHS}_manifest.json'), manifest)
    remove_legacy(output_dir, f'{DEPTHS}.json')
    return manifest


def load_existing_embeddings(output_dir):
    """Read previously written t-SNE / UMAP coordinates as {method: {id: [x, y]}}.

    Looks in the chunked explorer first, then in the legacy contig_<method>.json
    files, so a --skip-tsne / --skip-umap rerun keeps embeddings computed by an
    earlier stage whichever layout that stage wrote. Call before writing the
    explorer, which replaces both.
    """
    found = {m: None for m in EMBEDDINGS}
    manifest = read_json_any(os.path.join(output_dir, f'{EXPLORER}_manifest.json'))
    if manifest and any(manifest.get(f'has_{m}') for m in EMBEDDINGS):
        enc = manifest.get('id_encoding')
        for m in EMBEDDINGS:
            if manifest.get(f'has_{m}'):
                found[m] = {}
        for ch in manifest.get('chunks', []):
            body = read_json_any(os.path.join(output_dir, ch['file']))
            if body is None:
                continue
            ids = decode_ids(body['ids'], enc)
            for m in EMBEDDINGS:
                if found[m] is None:
                    continue
                xs = body['columns'].get(f'{m}_x')
                ys = body['columns'].get(f'{m}_y')
                if xs is None or ys is None:
                    continue
                for cid, x, y in zip(ids, xs, ys):
                    if x is not None and y is not None:
                        found[m][cid] = [x, y]
    for m in EMBEDDINGS:
        if not found[m]:
            legacy = read_json_any(os.path.join(output_dir, f'contig_{m}.json'))
            found[m] = legacy or None
    return found


# ---------------------------------------------------------------------------
# Gene features
# ---------------------------------------------------------------------------

def write_gene_shards(output_dir, genes, len_map, shard_size=DEFAULT_GENE_SHARD_SIZE,
                      min_length=0):
    """Write genes as shards capped by feature count.

    genes:      {contig_id: [feature, ...]}
    len_map:    {contig_id: length}; must match the explorer's lengths, since the
                SPA finds a contig's shard by its (length, id) sort key.
    shard_size: maximum features per shard. A contig is never split, so a
                contig with more features than the cap gets a shard of its own.
    """
    shard_size = max(1, int(shard_size))
    contigs = [c for c in genes if genes[c]]
    if min_length > 0:
        contigs = [c for c in contigs if int(len_map.get(c, 0)) >= min_length]
    contigs.sort(key=lambda c: sort_key(len_map.get(c, 0), c))

    shards = []
    cur, cur_n = [], 0

    def flush():
        k = len(shards)
        fname = part_name(GENES, k)
        write_part(os.path.join(output_dir, fname), {c: genes[c] for c in cur})
        first, last = cur[0], cur[-1]
        shards.append({
            'index': k, 'file': fname, 'n_contigs': len(cur), 'n_genes': cur_n,
            'first': [int(len_map.get(first, 0)), first],
            'last': [int(len_map.get(last, 0)), last],
        })

    for c in contigs:
        g = len(genes[c])
        if cur and cur_n + g > shard_size:
            flush()
            cur, cur_n = [], 0
        cur.append(c)
        cur_n += g
    if cur:
        flush()
    remove_stale_parts(output_dir, GENES, len(shards))

    manifest = {
        'format': 'danaseq-genes', 'version': FORMAT_VERSION,
        'sort': SORT_ORDER,
        'n_contigs': len(contigs),
        'n_genes': sum(s['n_genes'] for s in shards),
        'shard_gene_cap': shard_size,
        'min_length_cutoff': int(min_length),
        'shards': shards,
    }
    write_manifest(os.path.join(output_dir, f'{GENES}_manifest.json'), manifest)
    remove_legacy(output_dir, f'{GENES}.json')
    return manifest
