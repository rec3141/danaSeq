#!/usr/bin/env python3
"""Import protein-DB search results (HMM_SEARCH / DIAMOND_DB_SEARCH) into dana.duckdb.

Called from import_all.py with the per-barcode directory as cwd. Reads
`hmm/<chunk>.<db>.tsv` and `protdb/<chunk>.<db>.tsv` as written by
bin/protein_db_search.sh:

    #db=<db>  method=...  <cutoffs>  n_proteins=<N>
    gene_id db target gene category pathway subpathway evalue score pident qcov scov
    ...

Rows go to `protein_hits`; distinct (db, target) annotations are upserted into
`protein_db_annot`; one `protein_db_searches` row per file records how many
proteins were searched (the normalisation denominator). Idempotent per file via
`import_log`, like every other loader.

Legacy hmm/*.tsv files from the pre-catalogue HMM_SEARCH (header
`gene_id hmm_name evalue score`, no provenance line) are still accepted: the DB
name comes from the filename and n_proteins is counted from the chunk's
annotation .faa when it is still on disk.
"""

import glob
import gzip
import os
import sys


def _num(v, cast=float):
    try:
        return cast(v) if v not in ('', None, '-') else None
    except ValueError:
        return None


def _count_faa(chunk):
    """Proteins in the annotation .faa this chunk was searched from (legacy files)."""
    cands = [f'bakta/{chunk}/{chunk}.faa', f'bakta/{chunk}/{chunk}.faa.gz']
    cands += [f for f in glob.glob(f'prokka/{chunk}/*.faa*') if '.tmp.' not in f]
    for f in cands:
        if os.path.exists(f):
            op = gzip.open if f.endswith('.gz') else open
            with op(f, 'rt') as fh:
                return sum(1 for line in fh if line.startswith('>'))
    return None


def _parse(path):
    base = os.path.basename(path)[:-len('.tsv')]
    meta = {}
    rows, annot = [], {}
    with open(path) as fh:
        first = fh.readline().rstrip('\n')
        if first.startswith('#'):
            for kv in first[1:].split('\t'):
                k, _, v = kv.partition('=')
                meta[k.strip()] = v.strip()
            header = fh.readline().rstrip('\n').split('\t')
        else:
            header = first.split('\t')
        legacy = 'hmm_name' in header
        db = meta.get('db') or base.rsplit('.', 1)[-1]
        chunk = base[:-(len(db) + 1)] if base.endswith('.' + db) else base.rsplit('.', 1)[0]
        col = {h: i for i, h in enumerate(header)}
        for line in fh:
            f = line.rstrip('\n').split('\t')
            if len(f) < len(header) or not f[0]:
                continue
            if legacy:
                target = f[col['hmm_name']]
                rows.append((db, chunk, f[0], target, _num(f[col['evalue']]),
                             _num(f[col['score']]), None, None, None))
                continue
            target = f[col['target']]
            rows.append((db, chunk, f[0], target, _num(f[col['evalue']]), _num(f[col['score']]),
                         _num(f[col['pident']]), _num(f[col['qcov']]), _num(f[col['scov']])))
            if f[col['gene']] or f[col['category']]:
                annot[(db, target)] = (db, target, f[col['gene']], f[col['category']],
                                       f[col['pathway']], f[col['subpathway']])
    n_prot = _num(meta.get('n_proteins'), int)
    if n_prot is None:
        n_prot = _count_faa(chunk)
    params = '\t'.join(f'{k}={v}' for k, v in meta.items() if k not in ('db', 'n_proteins'))
    return db, chunk, rows, annot, n_prot, params


def import_protein_hits(con, imported):
    files = sorted(glob.glob('hmm/*.tsv') + glob.glob('protdb/*.tsv'))
    pending = [f for f in files if f not in imported]
    if not pending:
        return

    n_rows = 0
    for path in pending:
        try:
            db, chunk, rows, annot, n_prot, params = _parse(path)
        except Exception as e:
            print(f'[WARNING] protein hits read failed {path}: {e}', file=sys.stderr)
            continue
        # A chunk re-searched with other cutoffs (same filename) replaces its
        # previous rows instead of adding to them.
        con.execute('DELETE FROM protein_hits WHERE db = ? AND chunk = ?', [db, chunk])
        if rows:
            con.executemany('INSERT INTO protein_hits VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)', rows)
        if annot:
            con.executemany(
                'INSERT INTO protein_db_annot VALUES (?, ?, ?, ?, ?, ?) '
                'ON CONFLICT (db, target) DO UPDATE SET gene = excluded.gene, '
                'category = excluded.category, pathway = excluded.pathway, '
                'subpathway = excluded.subpathway',
                list(annot.values()))
        con.execute(
            'INSERT INTO protein_db_searches VALUES (?, ?, ?, ?, ?) '
            'ON CONFLICT (db, chunk) DO UPDATE SET n_proteins = excluded.n_proteins, '
            'n_hits = excluded.n_hits, params = excluded.params',
            [db, chunk, n_prot, len(rows), params])
        con.execute('INSERT INTO import_log (filename) VALUES (?) ON CONFLICT DO NOTHING', [path])
        n_rows += len(rows)
    print(f'[import_protein_hits] {len(pending)} files, {n_rows} hits', file=sys.stderr)


def _stale(con, imported):
    """Already-imported files whose provenance line no longer matches what was
    loaded (e.g. CANT-HYD re-searched with --cut_nc instead of --cut_tc)."""
    loaded = {(db, chunk): params for db, chunk, params in
              con.execute('SELECT db, chunk, params FROM protein_db_searches').fetchall()}
    stale = set()
    for path in sorted(glob.glob('hmm/*.tsv') + glob.glob('protdb/*.tsv')):
        if path not in imported:
            continue
        with open(path) as fh:
            first = fh.readline().rstrip('\n')
        if not first.startswith('#'):
            continue  # legacy file: nothing to compare
        meta = dict(kv.partition('=')[::2] for kv in first[1:].split('\t'))
        db = meta.get('db', '')
        base = os.path.basename(path)[:-len('.tsv')]
        chunk = base[:-(len(db) + 1)] if base.endswith('.' + db) else base
        params = '\t'.join(f'{k}={v}' for k, v in meta.items() if k not in ('db', 'n_proteins'))
        if loaded.get((db, chunk)) != params:
            stale.add(path)
    return stale


def main():
    """Standalone: import ONLY protein-DB hits for one barcode dir (used by
    backfills that must not touch the other tables, e.g. mapping).
    --recheck also re-imports files whose cutoffs changed since loading."""
    import duckdb
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from db_schema import ensure_schema
    args = [a for a in sys.argv[1:] if not a.startswith('--')]
    os.chdir(args[0])
    con = duckdb.connect('dana.duckdb')
    ensure_schema(con)
    imported = {r[0] for r in con.execute(
        "SELECT filename FROM import_log WHERE filename LIKE 'hmm/%' OR filename LIKE 'protdb/%'"
    ).fetchall()}
    if '--recheck' in sys.argv:
        imported -= _stale(con, imported)
    import_protein_hits(con, imported)
    con.close()


if __name__ == '__main__':
    main()
