#!/usr/bin/env python3
"""Load Kraken2 classification TSVs into DuckDB.

Usage: python3 kraken_db.py <barcode_dir>

Reads kraken/*.tsv, parses seqid + taxid + taxa_name, upserts into kraken table.
"""

import sys
import os
import re
import glob

import duckdb

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from db_schema import ensure_schema

TAXID_RE = re.compile(r'\(taxid (\d+)\)$')


def main():
    barcode_dir = sys.argv[1]
    os.chdir(barcode_dir)

    con = duckdb.connect('dana.duckdb')
    ensure_schema(con)

    kraken_files = sorted(glob.glob('kraken/*.tsv'))
    imported = {r[0] for r in con.execute(
        "SELECT filename FROM import_log WHERE filename LIKE 'kraken%.tsv'"
    ).fetchall()}
    pending = [f for f in kraken_files if f not in imported]

    for file in pending:
        rows = []
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
                    rows.append((seqid, taxid, taxa_name))
        except Exception as e:
            print(f"[WARNING] Failed to read {file}: {e}", file=sys.stderr)
            continue

        if not rows:
            con.execute(
                "INSERT INTO import_log (filename) VALUES (?) ON CONFLICT DO NOTHING",
                [file]
            )
            continue

        # Deduplicate by seqid within file
        seen = set()
        deduped = []
        for r in rows:
            if r[0] not in seen:
                seen.add(r[0])
                deduped.append(r)

        try:
            con.execute("BEGIN TRANSACTION")
            con.executemany(
                "INSERT INTO kraken (seqid, taxid, taxa_name) VALUES (?, ?, ?) "
                "ON CONFLICT (seqid) DO NOTHING",
                deduped
            )
            con.execute(
                "INSERT INTO import_log (filename) VALUES (?) ON CONFLICT DO NOTHING",
                [file]
            )
            con.execute("COMMIT")
        except Exception as e:
            con.execute("ROLLBACK")
            print(f"[WARNING] Import failed for {file}: {e}", file=sys.stderr)

    con.close()


if __name__ == '__main__':
    main()
