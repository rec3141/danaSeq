#!/usr/bin/env python3
"""Load Sendsketch similarity results into DuckDB.

Usage: python3 sketch_db.py <barcode_dir>

Reads sketch/*.txt, inserts into sendsketch table.
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

    sketch_files = sorted(glob.glob('sketch/*.txt'))
    imported = {r[0] for r in con.execute(
        "SELECT filename FROM import_log WHERE filename LIKE 'sketch%'"
    ).fetchall()}
    pending = [f for f in sketch_files if f not in imported]

    for file in pending:
        rows = []
        try:
            with open(file) as fh:
                for line in fh:
                    if line.startswith('#'):
                        continue
                    parts = line.rstrip('\n').split('\t')
                    if len(parts) < 3:
                        continue
                    fileid = parts[0].replace('.fa', '')
                    ref_name = parts[1]
                    try:
                        ani = float(parts[2])
                    except ValueError:
                        continue
                    rows.append((fileid, ref_name, ani))
        except Exception as e:
            print(f"[WARNING] Failed to read {file}: {e}", file=sys.stderr)
            continue

        if rows:
            con.executemany(
                "INSERT INTO sendsketch (fileid, ref_name, ani) VALUES (?, ?, ?)",
                rows
            )
        con.execute(
            "INSERT INTO import_log (filename) VALUES (?) ON CONFLICT DO NOTHING",
            [file]
        )

    con.close()


if __name__ == '__main__':
    main()
