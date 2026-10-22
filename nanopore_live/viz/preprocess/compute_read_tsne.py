#!/usr/bin/env python3
"""Compute read-level t-SNE from tetramer frequencies across all samples.

Reads tetra_data from per-sample DuckDB files, joins with kraken taxonomy
and stats, computes t-SNE on tetramer vectors, and outputs read_explorer.json
for the WebGL scatter plot.

With fastsne installed (https://github.com/rec3141/fastsne; `fastsne` on a GPU,
else `fastsne_cpu`) every read is mapped, not just a sample: each sample's
tetra_data is exported once to <output>/read_map/fsb/, a reference t-SNE is
built from a random sample of up to --ref-max reads (rebuilt when the read
count has grown --rebuild-factor times since), and the other reads are
projected onto it, only for samples whose data changed. read_explorer.json
keeps a --max-reads sample with coordinates from that map, and
<output>/read_map/viewer.html?data=reads.json shows all reads (fastsne viewer,
paged; projected reads can be hidden). Without fastsne, or with
--engine sklearn, the sampled sklearn t-SNE below runs as before.
"""

import argparse
import json
import os
import re
import shutil
import struct
import subprocess
import sys
import time
from pathlib import Path

try:
    import duckdb
except ImportError:
    print("[ERROR] duckdb required: pip install duckdb", file=sys.stderr)
    sys.exit(1)

try:
    import numpy as np
    from sklearn.decomposition import PCA
    from sklearn.manifold import TSNE
except ImportError:
    print("[ERROR] numpy + scikit-learn required: pip install numpy scikit-learn", file=sys.stderr)
    sys.exit(1)


def discover_samples(input_dir):
    """Find all sample directories with dana.duckdb."""
    samples = []
    input_path = Path(input_dir)
    for flowcell_dir in sorted(input_path.iterdir()):
        if not flowcell_dir.is_dir():
            continue
        flowcell = flowcell_dir.name
        for barcode_dir in sorted(flowcell_dir.iterdir()):
            if not barcode_dir.is_dir() or not barcode_dir.name.startswith('barcode'):
                continue
            db_path = barcode_dir / 'dana.duckdb'
            if db_path.exists():
                samples.append({
                    'id': f"{flowcell}:{barcode_dir.name}",
                    'db_path': str(db_path),
                    'dir': str(barcode_dir),
                })
    return samples


def load_lengths(sample_dir):
    """Load seqid → length from .lengths files in tetra/ directory."""
    lengths = {}
    tetra_dir = Path(sample_dir) / 'tetra'
    if not tetra_dir.is_dir():
        return lengths
    for f in tetra_dir.glob('*.lengths'):
        with open(f) as fh:
            for line in fh:
                parts = line.rstrip('\n').split('\t')
                if len(parts) == 2:
                    lengths[parts[0]] = int(parts[1])
    return lengths


# GTDB lineage prefix → field name on each read record
_GTDB_RANK_FIELD = {'d': 'gtdb_domain', 'p': 'gtdb_phylum', 'c': 'gtdb_class',
                    'o': 'gtdb_order', 'f': 'gtdb_family', 'g': 'gtdb_genus',
                    's': 'gtdb_species'}


def parse_gtdb_lineage_fields(lineage):
    """Parse a GTDB 'd__X;p__Y;...;s__Z' string into a flat dict keyed by
    gtdb_domain / gtdb_phylum / ... Empty rank tokens are skipped so the SPA
    shows '?' rather than a phantom 'g__' blank.
    """
    out = {}
    if not lineage:
        return out
    for tok in lineage.split(';'):
        tok = tok.strip()
        if len(tok) < 4 or tok[1:3] != '__':
            continue
        field = _GTDB_RANK_FIELD.get(tok[0])
        if not field:
            continue
        name = tok[3:].strip()
        if name:
            out[field] = name
    return out


def build_gtdb_read_lookup(con):
    """read_id → {gtdb_domain, gtdb_phylum, ...} from the sendsketch table.

    Only classified rows (status='C') contribute; unclassified reads simply
    lack gtdb_* fields on their scatter record and show as 'Unclassified'.
    """
    lookup = {}
    try:
        rows = con.execute(
            "SELECT read_id, lineage FROM sendsketch WHERE status = 'C' AND lineage IS NOT NULL"
        ).fetchall()
        for read_id, lineage in rows:
            fields = parse_gtdb_lineage_fields(lineage)
            if fields:
                lookup[read_id] = fields
    except Exception:
        pass  # sendsketch may not exist / be empty
    return lookup


def build_lineage_lookup(con):
    """Build taxid → {domain, phylum, class, order, family, genus} from krakenreport.

    Uses the taxonomy (semicolon-separated lineage) and rank_full (space-separated
    rank codes) columns which are parallel arrays, so zipping them maps each
    ancestor name to its rank."""
    lookup = {}
    try:
        rows = con.execute(
            "SELECT taxid, taxonomy, rank_full FROM krakenreport "
            "WHERE taxonomy IS NOT NULL AND rank_full IS NOT NULL"
        ).fetchall()
        for taxid, taxonomy, rank_full in rows:
            tax_parts = [x.strip() for x in taxonomy.split(';')]
            rank_parts = rank_full.split()
            if len(tax_parts) != len(rank_parts):
                continue
            lineage = dict(zip(rank_parts, tax_parts))
            # R2 = domain/superkingdom (Bacteria/Archaea/Eukaryota)
            lookup[taxid] = {
                'domain':  lineage.get('R2'),
                'phylum':  lineage.get('P'),
                'class':   lineage.get('C'),
                'order':   lineage.get('O'),
                'family':  lineage.get('F'),
                'genus':   lineage.get('G'),
                'species': lineage.get('S'),
            }
    except Exception:
        pass  # krakenreport may not exist
    return lookup


def extract_reads(db_path, sample_id, max_per_sample, seq_lengths=None, want_tetra=True):
    """Extract reads with tetramer, stats, and taxonomy from a single DuckDB.
    want_tetra=False skips the tetramer values (the fastsne path only needs the metadata)."""
    reads = []
    tetra_rows = []
    try:
        con = duckdb.connect(db_path, read_only=True)

        # Check what tables exist
        tables = [r[0] for r in con.execute("SHOW TABLES").fetchall()]

        if 'tetra_data' not in tables:
            con.close()
            return [], []

        # Get tetramer column names (exclude seqid)
        cols = con.execute("SELECT column_name FROM information_schema.columns WHERE table_name='tetra_data' AND column_name != 'seqid'").fetchall()
        tetra_cols = [c[0] for c in cols]

        if not tetra_cols:
            con.close()
            return [], []

        # Build taxid → lineage lookup from krakenreport
        has_kraken = 'kraken' in tables
        has_krakenreport = 'krakenreport' in tables
        lineage_lookup = build_lineage_lookup(con) if has_krakenreport else {}

        # Build read_id → GTDB lineage fields from sendsketch (optional)
        has_sendsketch = 'sendsketch' in tables
        gtdb_lookup = build_gtdb_read_lookup(con) if has_sendsketch else {}

        # Build read_id → gene annotations lookup via read_contig_map
        has_rcm = 'read_contig_map' in tables
        has_annotations = 'prokka_annotations' in tables
        gene_lookup = {}  # read_id → "gene1; gene2; ..."
        product_lookup = {}  # read_id → "product1; product2; ..."
        if has_rcm and has_annotations:
            try:
                ann_rows = con.execute("""
                    SELECT m.read_id,
                           GROUP_CONCAT(DISTINCT CASE WHEN p.gene != '' THEN p.gene END, '; '),
                           GROUP_CONCAT(DISTINCT CASE WHEN p.product != '' AND p.product != 'hypothetical protein' THEN p.product END, '; ')
                    FROM read_contig_map m
                    JOIN locus_index l ON m.contig_id = l.seqid
                    JOIN prokka_annotations p ON l.locus_tag = p.locus_tag
                    GROUP BY m.read_id
                """).fetchall()
                for read_id, genes, products in ann_rows:
                    if genes:
                        gene_lookup[read_id] = genes
                    if products:
                        product_lookup[read_id] = products
            except Exception:
                pass  # read_contig_map may be empty

        # Check which optional columns exist in stats
        has_stats = 'stats' in tables
        stats_cols = set()
        if has_stats:
            stats_cols = {r[0] for r in con.execute(
                "SELECT column_name FROM information_schema.columns WHERE table_name='stats'"
            ).fetchall()}

        has_gc = 'gc' in stats_cols
        has_length = 'length' in stats_cols

        # Build SELECT clause dynamically
        tetra_col_str = ', '.join(f't."{c}"' for c in tetra_cols)

        select_parts = ['t.seqid']
        if has_stats and has_length:
            select_parts.append('s.length')
        else:
            select_parts.append('NULL as length')
        if has_stats and has_gc:
            select_parts.append('s.gc')
        else:
            select_parts.append('NULL as gc')
        if has_kraken:
            select_parts.append('k.taxa_name')
        else:
            select_parts.append('NULL as taxa_name')
        if want_tetra:
            select_parts.append(tetra_col_str)
        else:
            # GC estimate from the tetramer frequencies, computed in the query
            gc_terms = ' + '.join(f't."{c}" * {(c.count("G") + c.count("C")) / 4.0}' for c in tetra_cols
                                  if c.count("G") + c.count("C"))
            select_parts.append(f'({gc_terms}) * 100 AS gc_est')

        select_str = ', '.join(select_parts)

        # Build FROM/JOIN clause
        from_str = 'tetra_data t'
        if has_stats:
            from_str += ' LEFT JOIN stats s ON t.seqid = s.seqid'
        if has_kraken:
            from_str += ' LEFT JOIN kraken k ON t.seqid = k.seqid'

        query = f"SELECT {select_str} FROM {from_str} USING SAMPLE {max_per_sample}"
        rows = con.execute(query).fetchall()

        for row in rows:
            seqid = row[0]
            length = row[1]
            gc = row[2]
            taxa_name = row[3]
            tetra_values = [float(v) if v is not None else 0.0 for v in row[4:]] if want_tetra else []

            # Use .lengths file if stats join didn't provide length
            if length is None and seq_lengths:
                length = seq_lengths.get(seqid)

            # Estimate GC from tetramer frequencies if not available
            if gc is None and not want_tetra and row[4] is not None:
                gc = round(float(row[4]), 1)
            elif gc is None and tetra_values:
                gc_est = sum(
                    tetra_values[i] * (tetra_cols[i].count('G') + tetra_cols[i].count('C')) / 4.0
                    for i in range(len(tetra_cols))
                ) * 100
                gc = round(gc_est, 1)
            elif gc is not None:
                gc = round(gc, 1)

            # Resolve full lineage via taxid → krakenreport lookup
            lineage = {}
            if taxa_name:
                m = re.search(r'\(taxid[_ ](\d+)\)', str(taxa_name))
                if m:
                    lineage = lineage_lookup.get(int(m.group(1)), {})

            read = {
                'id': seqid,
                'sample': sample_id,
                'length': length,
                'gc': gc,
                'kraken_domain':  lineage.get('domain'),
                'kraken_phylum':  lineage.get('phylum'),
                'kraken_class':   lineage.get('class'),
                'kraken_order':   lineage.get('order'),
                'kraken_family':  lineage.get('family'),
                'kraken_genus':   lineage.get('genus'),
                'kraken_species': lineage.get('species'),
            }
            # Overlay GTDB fields keyed by read_id (sendsketch was per-read)
            g = gtdb_lookup.get(seqid)
            if g:
                read.update(g)
            genes = gene_lookup.get(seqid)
            products = product_lookup.get(seqid)
            if genes:
                read['genes'] = genes
            if products:
                read['products'] = products
            reads.append(read)
            if want_tetra:
                tetra_rows.append(tetra_values)

        con.close()
    except Exception as e:
        print(f"  [WARNING] Failed to extract reads from {db_path}: {e}", file=sys.stderr)

    return reads, tetra_rows


# ---------------------------------------------------------------------------
# fastsne: map every read (reference t-SNE + projection)
# ---------------------------------------------------------------------------

def find_fastsne(engine):
    """(binary, 'gpu'|'cpu') for the requested engine, or None. FASTSNE_BIN overrides."""
    if engine == 'sklearn':
        return None
    if os.environ.get('FASTSNE_BIN'):
        b = os.environ['FASTSNE_BIN']
        return b, 'gpu' if os.path.basename(b) == 'fastsne' else 'cpu'
    gpu = shutil.which('fastsne')
    if gpu and engine in ('auto', 'gpu'):
        try:
            if 'GPU' in subprocess.run(['nvidia-smi', '-L'], capture_output=True, text=True, timeout=20).stdout:
                return gpu, 'gpu'
        except (OSError, subprocess.SubprocessError):
            pass
    cpu = shutil.which('fastsne_cpu')
    if cpu and engine in ('auto', 'cpu'):
        return cpu, 'cpu'
    if engine in ('gpu', 'cpu'):
        print(f"[ERROR] --engine {engine}: no usable fastsne binary", file=sys.stderr)
        sys.exit(1)
    return None


def fastsne_share(binary):
    """<prefix>/share/fastsne (viewer.html, tools/) for an installed binary."""
    return Path(os.path.realpath(binary)).parent.parent / 'share' / 'fastsne'


def export_fsb(sample, fsb_dir):
    """tetra_data (+ kraken taxon) of one sample -> <flowcell>__<barcode>.fsb and .meta.tsv.
    .fsb: b"FSB1", int32 D, int64 N, N*D float32, then the N ids joined by newlines. Skipped when
    both files are newer than the database. Returns (fsb path, rows) or None."""
    name = sample['id'].replace(':', '__')
    fsb, meta = fsb_dir / f'{name}.fsb', fsb_dir / f'{name}.meta.tsv'
    db_mtime = os.path.getmtime(sample['db_path'])
    if fsb.exists() and meta.exists() and fsb.stat().st_mtime >= db_mtime:
        with open(fsb, 'rb') as f:
            n = struct.unpack('<q', f.read(16)[8:16])[0]
        return (fsb, n) if n > 0 else None
    try:
        con = duckdb.connect(sample['db_path'], read_only=True)
        tables = {r[0] for r in con.execute("SHOW TABLES").fetchall()}
        if 'tetra_data' not in tables:
            con.close()
            return None
        cols = [r[0] for r in con.execute(
            "SELECT column_name FROM information_schema.columns "
            "WHERE table_name='tetra_data' AND column_name != 'seqid' ORDER BY ordinal_position").fetchall()]
        sel = ', '.join(f'coalesce(t."{c}", 0)::FLOAT AS "c{i}"' for i, c in enumerate(cols))
        taxon = 'k.taxa_name' if 'kraken' in tables else 'NULL'
        join = 'LEFT JOIN kraken k ON t.seqid = k.seqid' if 'kraken' in tables else ''
        res = con.execute(f"SELECT t.seqid AS id, {taxon} AS taxon, {sel} FROM tetra_data t {join}").fetchnumpy()
        con.close()
    except Exception as e:
        print(f"  [WARNING] export failed for {sample['db_path']}: {e}", file=sys.stderr)
        return None
    n = len(res['id'])
    if n == 0:
        return None
    x = np.nan_to_num(np.column_stack([np.asarray(res[f'c{i}'], dtype=np.float32) for i in range(len(cols))]))
    ids = [str(v) for v in res['id']]
    tx = res['taxon']   # a masked array where taxa are NULL
    taxa = ['' if m or v is None else str(v) for v, m in zip(np.ma.getdata(tx), np.ma.getmaskarray(tx))]
    tmp, tmpm = fsb.with_suffix('.fsb.tmp'), meta.with_suffix('.tsv.tmp')
    with open(tmp, 'wb') as f:
        f.write(b'FSB1' + struct.pack('<iq', x.shape[1], x.shape[0]))
        x.tofile(f)
        f.write('\n'.join(ids).encode())
    with open(tmpm, 'w') as f:
        f.write('id\ttaxon\n')
        f.writelines(f"{i}\t{t}\n" for i, t in zip(ids, taxa))
    os.replace(tmpm, meta)
    os.replace(tmp, fsb)
    return fsb, n


def run_fastsne(cmd, log):
    print(f"[INFO] {' '.join(str(c) for c in cmd)}", file=sys.stderr)
    with open(log, 'a') as lf:
        r = subprocess.run([str(c) for c in cmd], stdout=lf, stderr=subprocess.STDOUT)
    if r.returncode != 0:
        print(f"[ERROR] fastsne failed (exit {r.returncode}), see {log}", file=sys.stderr)
        sys.exit(1)


_HIGHER = {"Bacteria", "Archaea", "Eukaryota", "Viruses", "cellular", "root", "Bacilli", "Clostridia"}
_HIGHER_SUFFIX = ("aceae", "ales", "ota", "ati", "idae", "ineae", "proteobacteria", "mycetes", "bacteriia", "phyceae")


def kraken_genus(name, cache={}):
    """Kraken taxon name -> genus ('Ca. X' for Candidatus; '(above genus)' for higher ranks; '' unclassified)."""
    g = cache.get(name)
    if g is None:
        t = re.sub(r"\s*\(taxid[ _]\d+\)\s*$", "", str(name).replace("_", " ")).split()
        while t and t[0] == "unclassified":
            t = t[1:]
        if not t:
            g = ""
        elif t[0] == "Candidatus" and len(t) > 1:
            g = "Ca. " + t[1]
        elif t[0] in ("cellular", "root") or (len(t) == 1 and (t[0] in _HIGHER or t[0].endswith(_HIGHER_SUFFIX))):
            g = "(above genus)"
        else:
            g = t[0]
        cache[name] = g
    return g


def build_read_map(samples, out_dir, fs, args):
    """Export, (re)build the reference when needed, project changed samples. Returns the map table
    (id, x, y, nn_dist, kind, sample, genus) for every read, or None if there are too few reads."""
    import pandas as pd
    binary, kind = fs
    work = Path(out_dir) / 'read_map'
    fsb_dir, proj_dir = work / 'fsb', work / 'proj'
    fsb_dir.mkdir(parents=True, exist_ok=True)
    proj_dir.mkdir(exist_ok=True)
    log = work / 'fastsne.log'
    threads = ['--threads', args.threads] if kind == 'cpu' else []

    t0 = time.time()
    files = []   # (sample id, fsb path, rows)
    for s in samples:
        r = export_fsb(s, fsb_dir)
        if r:
            files.append((s['id'], r[0], r[1]))
    total = sum(n for _, _, n in files)
    print(f"[INFO] fastsne ({kind}): {total:,} reads in {len(files)} samples (export {time.time() - t0:.0f}s)",
          file=sys.stderr)
    if total < 1000:
        return None

    # reference: a uniform random sample (fastsne --sample hashes each row, so file order does not matter)
    state_path, model, ref_tsv = work / 'state.json', work / 'ref.model', work / 'ref.tsv'
    state = json.loads(state_path.read_text()) if state_path.exists() else {}
    if not model.exists() or total >= args.rebuild_factor * state.get('total_at_build', 0):
        ref_max = args.ref_max or (20_000_000 if kind == 'gpu' else 10_000_000)
        frac = min(1.0, ref_max / total)
        lst = work / 'all_fsb.txt'
        lst.write_text(''.join(f'{f}\n' for _, f, _ in files))
        run_fastsne([binary, '-i', f'@{lst}', '--sample', f'{frac:.6g}', '--seed', 7, '--dedup', 2,
                     '-o', ref_tsv, '--save-model', model] + threads, log)
        state = {'total_at_build': total, 'sample': frac, 'built': time.strftime('%Y-%m-%d %H:%M:%S')}
        state_path.write_text(json.dumps(state))
        print(f"[INFO] reference built from {frac:.1%} of {total:,} reads", file=sys.stderr)

    # project samples whose export or the model is newer than their projection, in one fastsne call
    stale = [(sid, f, n) for sid, f, n in files
             if not (proj_dir / (f.stem + '.tsv')).exists()
             or (proj_dir / (f.stem + '.tsv')).stat().st_mtime < max(f.stat().st_mtime, model.stat().st_mtime)]
    if stale:
        lst, tmp = work / 'project_fsb.txt', work / 'project.tsv'
        lst.write_text(''.join(f'{f}\n' for _, f, _ in stale))
        run_fastsne([binary, '--project', model, '-i', f'@{lst}', '-o', tmp] + threads, log)
        with open(tmp) as src:
            header = src.readline()
            for _, f, n in stale:
                with open(proj_dir / (f.stem + '.tsv'), 'w') as dst:
                    dst.write(header)
                    for _ in range(n):
                        dst.write(src.readline())
        tmp.unlink()
        print(f"[INFO] projected {sum(n for _, _, n in stale):,} reads of {len(stale)} changed samples", file=sys.stderr)

    # one table for all reads; reference reads keep their embedded positions
    parts = []
    for sid, f, n in files:
        p = pd.read_csv(proj_dir / (f.stem + '.tsv'), sep='\t', dtype={'id': str})
        m = pd.read_csv(f.with_suffix('.meta.tsv'), sep='\t', dtype=str, keep_default_na=False)
        p['sample'] = sid
        p['genus'] = [kraken_genus(t) if t else '' for t in m['taxon']]
        parts.append(p)
    tab = pd.concat(parts, ignore_index=True)
    ref = pd.read_csv(ref_tsv, sep='\t', dtype={'id': str}).drop_duplicates('id').set_index('id')
    isref = tab['id'].isin(ref.index).to_numpy()
    tab.loc[isref, ['x', 'y']] = ref.loc[tab.loc[isref, 'id'], ['x', 'y']].to_numpy()
    tab['kind'] = np.where(isref, 'reference', 'projected')
    return tab


def write_read_viewer(tab, out_dir, fs):
    """Paged fastsne viewer of every read: <out_dir>/read_map/viewer.html?data=reads.json"""
    share = fastsne_share(fs[0])
    pack, viewer = share / 'tools' / 'pack_viewer.py', share / 'viewer.html'
    if not pack.exists() or not viewer.exists():
        print(f"[WARNING] {share} has no pack_viewer.py / viewer.html; skipping the all-reads viewer", file=sys.stderr)
        return
    vdir = Path(out_dir) / 'read_map'
    tsv = vdir / 'map_table.tsv'
    tab[['id', 'x', 'y', 'nn_dist', 'kind', 'sample', 'genus']].to_csv(tsv, sep='\t', index=False)
    for old in vdir.glob('reads.p*.fsv'):
        old.unlink()
    cmd = [sys.executable, pack, vdir / 'reads.json', '--emb', f't-SNE={tsv}', '--kind', f'{tsv}:kind',
           '--layer', f'genus={tsv}:genus', '--layer', f'sample={tsv}:sample',
           '--layer', f'distance to reference={tsv}:nn_dist~log',
           '--title', f'{len(tab):,} reads ({(tab.kind == "reference").sum():,} reference + projected)']
    r = subprocess.run([str(c) for c in cmd], capture_output=True, text=True)
    tsv.unlink()
    if r.returncode != 0:
        print(f"[WARNING] pack_viewer failed: {r.stderr[-500:]}", file=sys.stderr)
        return
    shutil.copy(viewer, vdir / 'viewer.html')
    print(f"[SUCCESS] all-reads viewer: {vdir / 'viewer.html'}?data=reads.json", file=sys.stderr)


def write_read_explorer_fastsne(samples, tab, args):
    """read_explorer.json: a --max-reads sample (metadata as for the sklearn path) placed on the full map."""
    xy = tab.set_index('id')[['x', 'y']]
    xy = xy[~xy.index.duplicated()]
    max_per_sample = max(100, args.max_reads // max(1, len(samples)))
    all_reads = []
    for s in samples:
        reads, _ = extract_reads(s['db_path'], s['id'], max_per_sample, load_lengths(s['dir']), want_tetra=False)
        all_reads.extend(reads)
    all_reads = [r for r in all_reads if r['id'] in xy.index][:args.max_reads]
    if all_reads:
        c = xy.loc[[r['id'] for r in all_reads]].to_numpy()
        for r, (x, y) in zip(all_reads, c):
            r['tsne_x'] = round(float(x), 3)
            r['tsne_y'] = round(float(y), 3)
    output_path = os.path.join(args.output, 'read_explorer.json')
    tmp = output_path + '.tmp'
    with open(tmp, 'w') as f:
        json.dump({'reads': all_reads}, f, separators=(',', ':'))
    os.replace(tmp, output_path)
    size = os.path.getsize(output_path) / 1024 / 1024
    print(f"[SUCCESS] read_explorer.json: {size:.1f} MB ({len(all_reads)} of {len(tab):,} mapped reads)", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(description='Compute read-level t-SNE')
    parser.add_argument('--input', required=True, help='Directory with nanopore_live outputs')
    parser.add_argument('--output', required=True, help='Output directory')
    parser.add_argument('--max-reads', type=int, default=1000000,
                        help='Maximum reads in read_explorer.json (and in the t-SNE itself without fastsne)')
    parser.add_argument('--engine', choices=['auto', 'gpu', 'cpu', 'sklearn'], default='auto',
                        help='auto: fastsne on a GPU, else fastsne_cpu, else sklearn')
    parser.add_argument('--ref-max', type=int, default=0,
                        help='fastsne reference size (default 20M on a GPU, 10M on CPU); the rest is projected')
    parser.add_argument('--rebuild-factor', type=float, default=2.0,
                        help='rebuild the fastsne reference when reads grow this many times since it was built')
    parser.add_argument('--no-viewer', action='store_true', help='skip the all-reads fastsne viewer pages')
    parser.add_argument('--threads', type=int, default=max(1, len(os.sched_getaffinity(0)) // 2),
                        help='fastsne_cpu threads (default: half the available cores, leaving room for '
                             'basecalling and classification running alongside)')
    args = parser.parse_args()

    samples = discover_samples(args.input)
    print(f"[INFO] Found {len(samples)} samples", file=sys.stderr)

    fs = find_fastsne(args.engine)
    if fs:
        tab = build_read_map(samples, args.output, fs, args)
        if tab is not None:
            if not args.no_viewer:
                write_read_viewer(tab, args.output, fs)
            write_read_explorer_fastsne(samples, tab, args)
            return
        print("[INFO] too few reads for fastsne; using sklearn", file=sys.stderr)

    # Calculate per-sample budget
    max_per_sample = max(100, args.max_reads // max(1, len(samples)))
    print(f"[INFO] Max {max_per_sample} reads per sample, {args.max_reads} total", file=sys.stderr)

    # Extract reads from all samples
    all_reads = []
    all_tetra = []

    for i, s in enumerate(samples):
        print(f"  [{i+1}/{len(samples)}] {s['id']}...", file=sys.stderr, end='', flush=True)
        seq_lengths = load_lengths(s['dir'])
        reads, tetra = extract_reads(s['db_path'], s['id'], max_per_sample, seq_lengths)
        all_reads.extend(reads)
        all_tetra.extend(tetra)
        print(f" {len(reads)} reads", file=sys.stderr)

        # Stop if we hit the total limit
        if len(all_reads) >= args.max_reads:
            all_reads = all_reads[:args.max_reads]
            all_tetra = all_tetra[:args.max_reads]
            print(f"[INFO] Reached {args.max_reads} read limit", file=sys.stderr)
            break

    if len(all_reads) < 10:
        print("[WARNING] Too few reads for t-SNE. Writing empty read_explorer.json", file=sys.stderr)
        output_path = os.path.join(args.output, 'read_explorer.json')
        with open(output_path, 'w') as f:
            json.dump({'reads': []}, f)
        return

    print(f"[INFO] Computing t-SNE on {len(all_reads)} reads...", file=sys.stderr)

    # Build tetramer matrix
    tetra_matrix = np.array(all_tetra, dtype=np.float32)

    # Handle any NaN/inf, then fourth-root transform to reduce dominance of
    # high-frequency tetramers and better resolve compositional differences
    tetra_matrix = np.nan_to_num(tetra_matrix, nan=0.0, posinf=0.0, neginf=0.0)
    tetra_matrix = np.power(tetra_matrix, 0.25)

    # PCA to 50 dims first (matches Rtsne default, massive speedup)
    n_pca = min(50, tetra_matrix.shape[1], tetra_matrix.shape[0])
    print(f"[INFO] PCA {tetra_matrix.shape[1]} → {n_pca} dims", file=sys.stderr)
    tetra_matrix = PCA(n_components=n_pca, svd_solver='randomized', random_state=42).fit_transform(tetra_matrix)

    # Compute t-SNE (Barnes-Hut approximation for O(n log n) scaling)
    perplexity = min(30, max(5, len(all_reads) // 10))
    tsne = TSNE(n_components=2, perplexity=perplexity, random_state=42,
                max_iter=1000, learning_rate='auto', init='pca',
                method='barnes_hut', n_jobs=-1)
    coords = tsne.fit_transform(tetra_matrix)

    # Add coordinates to reads
    for i, read in enumerate(all_reads):
        read['tsne_x'] = round(float(coords[i, 0]), 3)
        read['tsne_y'] = round(float(coords[i, 1]), 3)

    # Write output
    output_path = os.path.join(args.output, 'read_explorer.json')
    with open(output_path, 'w') as f:
        json.dump({'reads': all_reads}, f, separators=(',', ':'))

    size = os.path.getsize(output_path) / 1024 / 1024
    print(f"[SUCCESS] read_explorer.json: {size:.1f} MB ({len(all_reads)} reads)", file=sys.stderr)


if __name__ == '__main__':
    main()
