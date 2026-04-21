#!/usr/bin/env python3
"""Load per-read sendsketch GTDB classifications into DuckDB.

Usage: python3 sketch_db.py <barcode_dir>

Reads sketch/*.sendsketch_reads.tsv (produced by modules/sketch.nf) and
inserts rows into the sendsketch table.
"""

import sys
import os
import glob

import duckdb

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from db_schema import ensure_schema


def main():
    barcode_dir = sys.argv[1]
    os.chdir(barcode_dir)

    con = duckdb.connect('dana.duckdb')
    ensure_schema(con)

    tsvs = sorted(glob.glob('sketch/*.sendsketch_reads.tsv'))
    imported = {r[0] for r in con.execute(
        "SELECT filename FROM import_log WHERE filename LIKE 'sketch%'"
    ).fetchall()}
    pending = [f for f in tsvs if f not in imported]

    for file in pending:
        rows = []
        try:
            with open(file) as fh:
                header = fh.readline()  # skip header
                for line in fh:
                    parts = line.rstrip('\n').split('\t')
                    if len(parts) < 5:
                        continue
                    read_id, status, ani_s, ref_name, lineage = parts[:5]
                    try:
                        ani = float(ani_s)
                    except ValueError:
                        continue
                    rows.append((read_id, status, ani, ref_name, lineage))
        except Exception as e:
            print(f"[WARNING] Failed to read {file}: {e}", file=sys.stderr)
            continue

        if rows:
            con.executemany(
                "INSERT INTO sendsketch (read_id, status, ani, ref_name, lineage) VALUES (?, ?, ?, ?, ?)",
                rows
            )
        con.execute(
            "INSERT INTO import_log (filename) VALUES (?) ON CONFLICT DO NOTHING",
            [file]
        )

    con.close()


if __name__ == '__main__':
    main()
