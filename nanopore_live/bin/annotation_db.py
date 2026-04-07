#!/usr/bin/env python3
"""Load Bakta/Prokka annotations into DuckDB.

Usage: python3 annotation_db.py <barcode_dir>

Handles both bakta/ and prokka/ output directories.
Fixes contig name collision by prefixing contig_N with the batch file ID.
Builds read_contig_map linking original read UUIDs to prefixed contig names.

Tables written: prokka_annotations, locus_index, stats, read_contig_map
"""

import sys
import os
import re
import glob

import duckdb

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from db_schema import ensure_schema

EC_RE = re.compile(r'EC:(\S+)')
COG_RE = re.compile(r'COG:(\S+)')
LOCUS_TAG_RE = re.compile(r'locus_tag=([^;]+)')
GFF_ID_RE = re.compile(r'ID=([^;]+)')


def read_bakta_tsv(path):
    """Parse bakta annotation TSV into prokka_annotations rows."""
    rows = []
    header = None
    with open(path) as fh:
        for line in fh:
            if line.startswith('#Sequence') or line.startswith('#sequence'):
                header = line.lstrip('#').rstrip('\n').lower().replace(' ', '_').split('\t')
                continue
            if line.startswith('#'):
                continue
            if header is None:
                continue
            fields = line.rstrip('\n').split('\t')
            row = dict(zip(header, fields))

            ftype = row.get('type', 'CDS').upper()
            if ftype not in ('CDS', 'TRNA', 'RRNA'):
                continue

            locus_tag = row.get('locus_tag', '')
            gene = row.get('gene', '') or ''
            product = row.get('product', '') or ''
            try:
                length_bp = abs(int(row.get('stop', 0)) - int(row.get('start', 0))) + 1
            except (ValueError, TypeError):
                length_bp = 0

            ec_number = ''
            cog = ''
            dbxrefs = row.get('dbxrefs', '') or ''
            if dbxrefs:
                m = EC_RE.search(dbxrefs)
                if m:
                    ec_number = m.group(1)
                m = COG_RE.search(dbxrefs)
                if m:
                    cog = m.group(1)

            rows.append((locus_tag, ftype, length_bp, gene, ec_number, cog, product))

    return rows


def read_bakta_gff(path, fileid):
    """Parse bakta GFF3 for sequence lengths and locus_tag→contig mapping.

    Returns (stats_rows, locus_rows, contig_names_in_order).
    Contig names are prefixed with fileid to prevent collision.
    """
    stats_rows = []
    locus_rows = []
    contig_order = []  # contig names in order of appearance

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
                    raw_contig = parts[1]
                    prefixed = f"{fileid}:{raw_contig}"
                    try:
                        length = int(parts[3])
                    except ValueError:
                        length = 0
                    stats_rows.append((prefixed, length))
                    contig_order.append(raw_contig)
                continue

            if line.startswith('#'):
                continue

            fields = line.rstrip('\n').split('\t')
            if len(fields) < 9:
                continue

            ftype = fields[2]
            if ftype not in ('CDS', 'tRNA', 'rRNA', 'cds'):
                continue

            raw_contig = fields[0]
            prefixed = f"{fileid}:{raw_contig}"
            attrs = fields[8]

            m = LOCUS_TAG_RE.search(attrs)
            if not m:
                m = GFF_ID_RE.search(attrs)
            if m:
                locus_rows.append((prefixed, m.group(1)))

    return stats_rows, locus_rows, contig_order


def read_prokka_tsv(path):
    """Parse prokka annotation TSV into prokka_annotations rows."""
    rows = []
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
                row.get('locus_tag', ''),
                ftype,
                int(row.get('length_bp', 0) or 0),
                row.get('gene', '') or '',
                row.get('ec_number', '') or '',
                row.get('cog', '') or '',
                row.get('product', '') or '',
            ))

    return rows


def read_prokka_gff(path, fileid):
    """Parse prokka GFF for sequence lengths and locus_tag→contig mapping.

    Returns (stats_rows, locus_rows, contig_names_in_order).
    """
    stats_rows = []
    locus_rows = []
    contig_order = []

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
                    raw_contig = parts[1]
                    prefixed = f"{fileid}:{raw_contig}"
                    try:
                        length = int(parts[3])
                    except ValueError:
                        length = 0
                    stats_rows.append((prefixed, length))
                    contig_order.append(raw_contig)
                continue

            if line.startswith('#'):
                continue

            fields = line.rstrip('\n').split('\t')
            if len(fields) < 9:
                continue

            ftype = fields[2]
            if ftype not in ('CDS', 'tRNA', 'rRNA'):
                continue

            raw_contig = fields[0]
            prefixed = f"{fileid}:{raw_contig}"
            attrs = fields[8]

            m = LOCUS_TAG_RE.search(attrs)
            if m:
                locus_rows.append((prefixed, m.group(1)))

    return stats_rows, locus_rows, contig_order


def read_fa_ids(fa_path):
    """Extract sequence IDs from a FASTA file, in order."""
    ids = []
    with open(fa_path) as f:
        for line in f:
            if line.startswith('>'):
                ids.append(line[1:].strip().split()[0])
    return ids


def build_read_contig_map(fa_path, fileid, contig_order):
    """Build UUID→prefixed_contig mapping from FA file and contig order.

    The Nth sequence in the FA corresponds to contig_N in bakta/prokka output.
    """
    read_ids = read_fa_ids(fa_path)
    rows = []
    for i, read_id in enumerate(read_ids):
        if i < len(contig_order):
            contig_id = f"{fileid}:{contig_order[i]}"
        else:
            contig_id = f"{fileid}:contig_{i + 1}"
        rows.append((read_id, contig_id, fileid))
    return rows


def import_batch(con, tsv_path, gff_path, fa_dir, is_bakta=True):
    """Import one annotation batch (TSV + GFF + optional FA mapping)."""
    # Derive fileid from directory name
    fileid = os.path.basename(os.path.dirname(tsv_path))

    # Parse TSV
    if is_bakta:
        ann_rows = read_bakta_tsv(tsv_path)
    else:
        ann_rows = read_prokka_tsv(tsv_path)

    if not ann_rows:
        return

    # Parse GFF
    if is_bakta:
        stats_rows, locus_rows, contig_order = read_bakta_gff(gff_path, fileid)
    else:
        stats_rows, locus_rows, contig_order = read_prokka_gff(gff_path, fileid)

    # Insert annotations (locus_tags are unique across batches)
    con.executemany(
        "INSERT INTO prokka_annotations VALUES (?, ?, ?, ?, ?, ?, ?) "
        "ON CONFLICT (locus_tag) DO NOTHING",
        ann_rows
    )

    # Insert stats with prefixed contig IDs
    if stats_rows:
        con.executemany("INSERT INTO stats VALUES (?, ?)", stats_rows)

    # Insert locus_index with prefixed contig IDs
    if locus_rows:
        con.executemany(
            "INSERT INTO locus_index VALUES (?, ?) ON CONFLICT (locus_tag) DO NOTHING",
            locus_rows
        )

    # Build read_contig_map from FA file
    fa_path = os.path.join(fa_dir, f"{fileid}.fa")
    if os.path.exists(fa_path) and contig_order:
        map_rows = build_read_contig_map(fa_path, fileid, contig_order)
        if map_rows:
            con.executemany(
                "INSERT INTO read_contig_map VALUES (?, ?, ?)",
                map_rows
            )
    elif not os.path.exists(fa_path):
        print(f"[WARNING] FA file not found for read_contig_map: {fa_path}", file=sys.stderr)


def collect_pending(pattern, gff_ext, imported, is_bakta):
    """Find pending TSV+GFF pairs, parse them, return accumulated rows."""
    tsvs = sorted(glob.glob(pattern, recursive=True))
    if is_bakta:
        tsvs = [f for f in tsvs
                if not f.endswith('.hypotheticals.tsv')
                and not f.endswith('.inference.tsv')]

    all_ann = []
    all_stats = []
    all_locus = []
    all_map = []
    log_entries = []

    for tsv_path in tsvs:
        if tsv_path in imported:
            continue
        gff_path = tsv_path.replace('.tsv', gff_ext)
        if not os.path.exists(gff_path):
            continue

        fileid = os.path.basename(os.path.dirname(tsv_path))
        try:
            if is_bakta:
                ann_rows = read_bakta_tsv(tsv_path)
                stats_rows, locus_rows, contig_order = read_bakta_gff(gff_path, fileid)
            else:
                ann_rows = read_prokka_tsv(tsv_path)
                stats_rows, locus_rows, contig_order = read_prokka_gff(gff_path, fileid)
        except Exception as e:
            print(f"[WARNING] Import failed for {tsv_path}: {e}", file=sys.stderr)
            continue

        if ann_rows:
            all_ann.extend(ann_rows)
        all_stats.extend(stats_rows)
        all_locus.extend(locus_rows)

        fa_path = os.path.join('fa', f"{fileid}.fa")
        if os.path.exists(fa_path) and contig_order:
            map_rows = build_read_contig_map(fa_path, fileid, contig_order)
            all_map.extend(map_rows)

        log_entries.append(tsv_path)
        log_entries.append(gff_path)

    return all_ann, all_stats, all_locus, all_map, log_entries


def main():
    barcode_dir = sys.argv[1]
    os.chdir(barcode_dir)

    con = duckdb.connect('dana.duckdb')
    ensure_schema(con)

    imported = {r[0] for r in con.execute(
        "SELECT filename FROM import_log"
    ).fetchall()}

    # Collect all pending rows from prokka + bakta
    all_ann, all_stats, all_locus, all_map, log_entries = [], [], [], [], []

    for pattern, gff_ext, is_bakta in [
        ('prokka/**/*.tsv', '.gff', False),
        ('bakta/**/*.tsv', '.gff3', True),
    ]:
        ann, stats, locus, rcm, logs = collect_pending(
            pattern, gff_ext, imported, is_bakta
        )
        all_ann.extend(ann)
        all_stats.extend(stats)
        all_locus.extend(locus)
        all_map.extend(rcm)
        log_entries.extend(logs)

    if not log_entries:
        con.close()
        return

    n = len(log_entries) // 2
    print(f"[INFO] annotation_db: importing {n} batches "
          f"({len(all_ann)} annotations, {len(all_map)} read mappings)",
          file=sys.stderr)

    # Bulk insert via appender (much faster than executemany for large batches)
    def bulk_append(table, rows):
        appender = con.appender(table)
        for row in rows:
            appender.append_row(row)
        appender.close()

    def bulk_upsert(table, rows, cols, pk):
        """Insert rows, skipping primary key conflicts via temp table."""
        tmp = f"_tmp_{table}"
        col_defs = ', '.join(f"{c} {t}" for c, t in cols)
        con.execute(f"CREATE TEMP TABLE {tmp} ({col_defs})")
        appender = con.appender(tmp)
        for row in rows:
            appender.append_row(row)
        appender.close()
        con.execute(f"INSERT INTO {table} SELECT * FROM {tmp} ON CONFLICT ({pk}) DO NOTHING")
        con.execute(f"DROP TABLE {tmp}")

    if all_ann:
        bulk_upsert('prokka_annotations', all_ann, [
            ('locus_tag', 'TEXT'), ('ftype', 'TEXT'), ('length_bp', 'INTEGER'),
            ('gene', 'TEXT'), ('ec_number', 'TEXT'), ('cog', 'TEXT'), ('product', 'TEXT'),
        ], 'locus_tag')

    if all_stats:
        bulk_append('stats', all_stats)

    if all_locus:
        bulk_upsert('locus_index', all_locus, [
            ('seqid', 'TEXT'), ('locus_tag', 'TEXT'),
        ], 'locus_tag')

    if all_map:
        bulk_append('read_contig_map', all_map)

    for entry in log_entries:
        con.execute(
            "INSERT INTO import_log (filename) VALUES (?) ON CONFLICT DO NOTHING",
            [entry])

    con.close()


if __name__ == '__main__':
    main()
