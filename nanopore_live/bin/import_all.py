#!/usr/bin/env python3
"""Import all pipeline results into DuckDB in a single process.

Usage: python3 import_all.py <barcode_dir>

Combines kraken, krakenreport, annotation, sketch, and tetra imports
into one process with one DB connection, avoiding repeated Python
startup and connection overhead.
"""

import sys
import os
import re
import glob

import duckdb
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from db_schema import ensure_schema

TAXID_RE = re.compile(r'\(taxid (\d+)\)$')
EC_RE = re.compile(r'EC:(\S+)')
COG_RE = re.compile(r'COG:(\S+)')
LOCUS_TAG_RE = re.compile(r'locus_tag=([^;]+)')
GFF_ID_RE = re.compile(r'ID=([^;]+)')


# ---- Kraken ----

def import_kraken(con, imported):
    files = sorted(glob.glob('kraken/*.tsv'))
    pending = [f for f in files if f not in imported]
    if not pending:
        return

    all_rows = []
    log = []
    for file in pending:
        try:
            with open(file) as fh:
                for line in fh:
                    if line.startswith('#'):
                        continue
                    parts = line.rstrip('\n').split('\t')
                    if len(parts) < 2:
                        continue
                    seqid = parts[0]
                    taxa_w_id = parts[1]
                    m = TAXID_RE.search(taxa_w_id)
                    taxid = int(m.group(1)) if m else None
                    taxa_name = TAXID_RE.sub('', taxa_w_id).rstrip(' ')
                    all_rows.append((seqid, taxid, taxa_name))
        except Exception as e:
            print(f"[WARNING] kraken read failed {file}: {e}", file=sys.stderr)
            continue
        log.append(file)

    if all_rows:
        # Deduplicate by seqid
        seen = set()
        deduped = []
        for r in all_rows:
            if r[0] not in seen:
                seen.add(r[0])
                deduped.append(r)
        df = pd.DataFrame(deduped, columns=['seqid', 'taxid', 'taxa_name'])
        con.execute("CREATE TEMP TABLE _kr AS SELECT * FROM kraken WHERE false")
        con.append('_kr', df)
        con.execute("INSERT INTO kraken SELECT * FROM _kr ON CONFLICT (seqid) DO NOTHING")
        con.execute("DROP TABLE _kr")

    for f in log:
        con.execute("INSERT INTO import_log (filename) VALUES (?) ON CONFLICT DO NOTHING", [f])


# ---- Krakenreport ----

def import_krakenreport(con, imported):
    files = sorted(glob.glob('kraken/*.report'))
    pending = [f for f in files if f not in imported]
    if not pending:
        return

    all_rows = []
    log = []
    for filepath in pending:
        try:
            basename = os.path.basename(filepath)
            parts = basename.split('_')
            flowcell = parts[0] if len(parts) > 0 else ''
            barcode = parts[2] if len(parts) > 2 else ''

            taxonomy_stack = []
            rank_stack = []

            with open(filepath) as fh:
                for line in fh:
                    line = line.rstrip('\n')
                    fields = line.split('\t')
                    if len(fields) < 6:
                        continue
                    try:
                        pct = float(fields[0])
                    except ValueError:
                        pct = 0.0
                    try:
                        reads = int(fields[1])
                    except ValueError:
                        reads = 0
                    try:
                        direct = int(fields[2])
                    except ValueError:
                        direct = 0
                    rank = fields[3]
                    try:
                        taxid = int(fields[4])
                    except ValueError:
                        taxid = 0
                    raw_name = fields[5]
                    stripped = raw_name.lstrip(' ')
                    n_spaces = len(raw_name) - len(stripped)
                    level = n_spaces // 2 + 1
                    name = stripped.strip()

                    if level <= len(taxonomy_stack):
                        taxonomy_stack = taxonomy_stack[:level - 1]
                        rank_stack = rank_stack[:level - 1]
                    taxonomy_stack.append(name)
                    rank_stack.append(rank)

                    all_rows.append((
                        pct, reads, direct, rank, taxid, name,
                        level, flowcell, barcode,
                        '; '.join(taxonomy_stack),
                        ' '.join(rank_stack)
                    ))
        except Exception as e:
            print(f"[WARNING] krakenreport read failed {filepath}: {e}", file=sys.stderr)
            continue
        log.append(filepath)

    if all_rows:
        df = pd.DataFrame(all_rows, columns=[
            'percent', 'reads', 'direct', 'rank', 'taxid', 'name',
            'level', 'flowcell', 'barcode', 'taxonomy', 'rank_full'])
        con.append('krakenreport', df)

    for f in log:
        con.execute("INSERT INTO import_log (filename) VALUES (?) ON CONFLICT DO NOTHING", [f])


# ---- Sketch ----

def import_sketch(con, imported):
    files = sorted(glob.glob('sketch/*.txt'))
    pending = [f for f in files if f not in imported]
    if not pending:
        return

    all_rows = []
    log = []
    for file in pending:
        try:
            with open(file) as fh:
                for line in fh:
                    if line.startswith('#'):
                        continue
                    parts = line.rstrip('\n').split('\t')
                    if len(parts) < 3:
                        continue
                    fileid = parts[0].replace('.fa', '')
                    try:
                        ani = float(parts[2])
                    except ValueError:
                        continue
                    all_rows.append((fileid, parts[1], ani))
        except Exception as e:
            print(f"[WARNING] sketch read failed {file}: {e}", file=sys.stderr)
            continue
        log.append(file)

    if all_rows:
        df = pd.DataFrame(all_rows, columns=['fileid', 'ref_name', 'ani'])
        con.append('sendsketch', df)

    for f in log:
        con.execute("INSERT INTO import_log (filename) VALUES (?) ON CONFLICT DO NOTHING", [f])


# ---- Tetra ----

def import_tetra(con, imported):
    if not os.path.exists('tnfs.txt'):
        return

    with open('tnfs.txt') as f:
        tnfs = f.readline().rstrip('\n').split('\t')

    # Create tetra_data table dynamically
    col_defs = ', '.join(
        f'"{c}" TEXT' if i == 0 else f'"{c}" REAL'
        for i, c in enumerate(tnfs)
    )
    con.execute(f"CREATE TABLE IF NOT EXISTS tetra_data ({col_defs})")

    files = sorted(glob.glob('tetra/*.lrn'))
    pending = [f for f in files if f not in imported]
    if not pending:
        return

    all_tetra = []
    all_seq = []
    log = []
    for file in pending:
        fileid = os.path.splitext(os.path.basename(file))[0]
        try:
            with open(file) as fh:
                for line in fh:
                    if line.startswith('%'):
                        continue
                    parts = line.rstrip('\n').split('\t')
                    if len(parts) != len(tnfs):
                        continue
                    seqid = parts[0]
                    values = [seqid]
                    for v in parts[1:]:
                        try:
                            values.append(float(v))
                        except ValueError:
                            values.append(0.0)
                    all_tetra.append(values)
                    all_seq.append((seqid, fileid))
        except Exception as e:
            print(f"[WARNING] tetra read failed {file}: {e}", file=sys.stderr)
            continue
        log.append(file)

    if all_seq:
        df_seq = pd.DataFrame(all_seq, columns=['seqid', 'fileid'])
        con.execute("CREATE TEMP TABLE _seq AS SELECT * FROM sequence_index WHERE false")
        con.append('_seq', df_seq)
        con.execute("INSERT INTO sequence_index SELECT * FROM _seq ON CONFLICT DO NOTHING")
        con.execute("DROP TABLE _seq")

    if all_tetra:
        df = pd.DataFrame(all_tetra, columns=tnfs)
        # First column is TEXT, rest REAL — ensure types
        for c in tnfs[1:]:
            df[c] = pd.to_numeric(df[c], errors='coerce').fillna(0.0)
        con.append('tetra_data', df)

    for f in log:
        con.execute("INSERT INTO import_log (filename) VALUES (?) ON CONFLICT DO NOTHING", [f])


# ---- Annotations (bakta/prokka) ----

def import_annotations(con, imported):
    """Import bakta and prokka annotations with contig collision fix."""

    all_ann = []
    all_stats = []
    all_locus = []
    all_map = []
    log_entries = []

    for pattern, gff_ext, is_bakta in [
        ('prokka/**/*.tsv', '.gff', False),
        ('bakta/**/*.tsv', '.gff3', True),
    ]:
        tsvs = sorted(glob.glob(pattern, recursive=True))
        if is_bakta:
            tsvs = [f for f in tsvs
                    if not f.endswith('.hypotheticals.tsv')
                    and not f.endswith('.inference.tsv')]

        for tsv_path in tsvs:
            if tsv_path in imported:
                continue
            gff_path = tsv_path.replace('.tsv', gff_ext)
            if not os.path.exists(gff_path):
                continue

            fileid = os.path.basename(os.path.dirname(tsv_path))
            try:
                # Parse TSV
                ann_rows = _read_annotation_tsv(tsv_path, is_bakta)
                # Parse GFF
                stats_rows, locus_rows, contig_order = _read_gff(gff_path, fileid, is_bakta)
            except Exception as e:
                print(f"[WARNING] annotation import failed {tsv_path}: {e}", file=sys.stderr)
                continue

            all_ann.extend(ann_rows)
            all_stats.extend(stats_rows)
            all_locus.extend(locus_rows)

            # Build read_contig_map
            fa_path = os.path.join('fa', f"{fileid}.fa")
            if os.path.exists(fa_path) and contig_order:
                uuids = []
                with open(fa_path) as f:
                    for line in f:
                        if line.startswith('>'):
                            uuids.append(line[1:].strip().split()[0])
                for i in range(min(len(uuids), len(contig_order))):
                    all_map.append((uuids[i], f"{fileid}:{contig_order[i]}", fileid))

            log_entries.append(tsv_path)
            log_entries.append(gff_path)

    if not log_entries:
        return

    n = len(log_entries) // 2
    print(f"[INFO] annotations: {n} batches, {len(all_ann)} genes, {len(all_map)} read mappings",
          file=sys.stderr)

    if all_ann:
        df = pd.DataFrame(all_ann, columns=[
            'locus_tag', 'ftype', 'length_bp', 'gene', 'ec_number', 'cog', 'product'])
        con.execute("CREATE TEMP TABLE _ann AS SELECT * FROM prokka_annotations WHERE false")
        con.append('_ann', df)
        con.execute("INSERT INTO prokka_annotations SELECT * FROM _ann ON CONFLICT (locus_tag) DO NOTHING")
        con.execute("DROP TABLE _ann")

    if all_stats:
        df = pd.DataFrame(all_stats, columns=['seqid', 'length'])
        con.append('stats', df)

    if all_locus:
        df = pd.DataFrame(all_locus, columns=['seqid', 'locus_tag'])
        con.execute("CREATE TEMP TABLE _loc AS SELECT * FROM locus_index WHERE false")
        con.append('_loc', df)
        con.execute("INSERT INTO locus_index SELECT * FROM _loc ON CONFLICT (locus_tag) DO NOTHING")
        con.execute("DROP TABLE _loc")

    if all_map:
        df = pd.DataFrame(all_map, columns=['read_id', 'contig_id', 'fileid'])
        con.append('read_contig_map', df)

    for entry in log_entries:
        con.execute("INSERT INTO import_log (filename) VALUES (?) ON CONFLICT DO NOTHING", [entry])


def _read_annotation_tsv(path, is_bakta):
    rows = []
    if is_bakta:
        header = None
        with open(path) as fh:
            for line in fh:
                if line.startswith('#Sequence') or line.startswith('#sequence'):
                    header = line.lstrip('#').rstrip('\n').lower().replace(' ', '_').split('\t')
                    continue
                if line.startswith('#') or header is None:
                    continue
                fields = line.rstrip('\n').split('\t')
                row = dict(zip(header, fields))
                ftype = row.get('type', 'CDS').upper()
                if ftype not in ('CDS', 'TRNA', 'RRNA'):
                    continue
                gene = row.get('gene', '') or ''
                product = row.get('product', '') or ''
                try:
                    length_bp = abs(int(row.get('stop', 0)) - int(row.get('start', 0))) + 1
                except (ValueError, TypeError):
                    length_bp = 0
                ec, cog = '', ''
                dbxrefs = row.get('dbxrefs', '') or ''
                if dbxrefs:
                    m = EC_RE.search(dbxrefs)
                    if m: ec = m.group(1)
                    m = COG_RE.search(dbxrefs)
                    if m: cog = m.group(1)
                rows.append((row.get('locus_tag', ''), ftype, length_bp, gene, ec, cog, product))
    else:
        header = None
        with open(path) as fh:
            for line in fh:
                if header is None:
                    header = line.rstrip('\n').lower().split('\t')
                    continue
                fields = line.rstrip('\n').split('\t')
                row = dict(zip(header, fields))
                ftype = row.get('ftype', '')
                if ftype not in ('CDS', 'rRNA', 'tRNA'):
                    continue
                rows.append((
                    row.get('locus_tag', ''), ftype,
                    int(row.get('length_bp', 0) or 0),
                    row.get('gene', '') or '', row.get('ec_number', '') or '',
                    row.get('cog', '') or '', row.get('product', '') or '',
                ))
    return rows


def _read_gff(path, fileid, is_bakta):
    stats_rows, locus_rows, contig_order = [], [], []
    with open(path) as fh:
        in_fasta = False
        for line in fh:
            if line.startswith('##FASTA'):
                in_fasta = True
                continue
            if in_fasta:
                continue
            if line.startswith('##sequence-region'):
                parts = line.strip().split()
                if len(parts) >= 4:
                    raw = parts[1]
                    prefixed = f"{fileid}:{raw}"
                    try:
                        length = int(parts[3])
                    except ValueError:
                        length = 0
                    stats_rows.append((prefixed, length))
                    contig_order.append(raw)
                continue
            if line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 9:
                continue
            ftype = fields[2]
            valid = ('CDS', 'tRNA', 'rRNA', 'cds') if is_bakta else ('CDS', 'tRNA', 'rRNA')
            if ftype not in valid:
                continue
            prefixed = f"{fileid}:{fields[0]}"
            m = LOCUS_TAG_RE.search(fields[8])
            if not m and is_bakta:
                m = GFF_ID_RE.search(fields[8])
            if m:
                locus_rows.append((prefixed, m.group(1)))
    return stats_rows, locus_rows, contig_order


# ---- Main ----

def main():
    barcode_dir = sys.argv[1]
    os.chdir(barcode_dir)

    con = duckdb.connect('dana.duckdb')
    ensure_schema(con)

    imported = {r[0] for r in con.execute("SELECT filename FROM import_log").fetchall()}

    import_kraken(con, imported)
    import_krakenreport(con, imported)
    import_annotations(con, imported)
    import_sketch(con, imported)
    import_tetra(con, imported)

    con.close()


if __name__ == '__main__':
    main()
