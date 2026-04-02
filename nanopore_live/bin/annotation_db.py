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


def main():
    barcode_dir = sys.argv[1]
    os.chdir(barcode_dir)

    con = duckdb.connect('dana.duckdb')
    ensure_schema(con)

    imported = {r[0] for r in con.execute(
        "SELECT filename FROM import_log"
    ).fetchall()}

    # --- Prokka ---
    prokka_tsvs = sorted(glob.glob('prokka/**/*.tsv', recursive=True))
    for tsv_path in prokka_tsvs:
        if tsv_path in imported:
            continue
        gff_path = tsv_path.replace('.tsv', '.gff')
        if not os.path.exists(gff_path):
            continue
        try:
            import_batch(con, tsv_path, gff_path, 'fa', is_bakta=False)
        except Exception as e:
            print(f"[WARNING] Prokka import failed for {tsv_path}: {e}", file=sys.stderr)
            continue
        con.execute("INSERT INTO import_log (filename) VALUES (?) ON CONFLICT DO NOTHING", [tsv_path])
        con.execute("INSERT INTO import_log (filename) VALUES (?) ON CONFLICT DO NOTHING", [gff_path])

    # --- Bakta ---
    bakta_tsvs = sorted(glob.glob('bakta/**/*.tsv', recursive=True))
    # Filter out .hypotheticals.tsv and .inference.tsv
    bakta_tsvs = [f for f in bakta_tsvs
                  if not f.endswith('.hypotheticals.tsv')
                  and not f.endswith('.inference.tsv')]
    for tsv_path in bakta_tsvs:
        if tsv_path in imported:
            continue
        gff_path = tsv_path.replace('.tsv', '.gff3')
        if not os.path.exists(gff_path):
            continue
        try:
            import_batch(con, tsv_path, gff_path, 'fa', is_bakta=True)
        except Exception as e:
            print(f"[WARNING] Bakta import failed for {tsv_path}: {e}", file=sys.stderr)
            continue
        con.execute("INSERT INTO import_log (filename) VALUES (?) ON CONFLICT DO NOTHING", [tsv_path])
        con.execute("INSERT INTO import_log (filename) VALUES (?) ON CONFLICT DO NOTHING", [gff_path])

    con.close()


if __name__ == '__main__':
    main()
