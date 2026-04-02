#!/usr/bin/env python3
"""Load tetranucleotide frequency data into DuckDB.

Usage: python3 tetra_db.py <barcode_dir>

Reads tnfs.txt for column names, creates tetra_data table dynamically,
loads tetra/*.lrn files, populates tetra_data and sequence_index tables.
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

    # Read column names from tnfs.txt
    with open('tnfs.txt') as f:
        tnfs = f.readline().rstrip('\n').split('\t')

    # Create tetra_data table dynamically
    col_defs = ', '.join(
        f'"{c}" TEXT' if i == 0 else f'"{c}" REAL'
        for i, c in enumerate(tnfs)
    )
    con.execute(f"CREATE TABLE IF NOT EXISTS tetra_data ({col_defs})")

    # Find pending LRN files
    lrn_files = sorted(glob.glob('tetra/*.lrn'))
    imported = {r[0] for r in con.execute(
        "SELECT filename FROM import_log WHERE filename LIKE 'tetra%'"
    ).fetchall()}
    pending = [f for f in lrn_files if f not in imported]

    for file in pending:
        fileid = os.path.splitext(os.path.basename(file))[0]

        # Parse LRN file (skip % comment lines)
        rows = []
        seq_rows = []
        try:
            with open(file) as fh:
                for line in fh:
                    if line.startswith('%'):
                        continue
                    parts = line.rstrip('\n').split('\t')
                    if len(parts) != len(tnfs):
                        continue
                    seqid = parts[0]
                    values = []
                    for v in parts[1:]:
                        try:
                            values.append(float(v))
                        except ValueError:
                            values.append(0.0)
                    rows.append([seqid] + values)
                    seq_rows.append((seqid, fileid))
        except Exception as e:
            print(f"[WARNING] Failed to read {file}: {e}", file=sys.stderr)
            continue

        if rows:
            # Insert into sequence_index
            con.executemany(
                "INSERT INTO sequence_index (seqid, fileid) VALUES (?, ?) "
                "ON CONFLICT DO NOTHING",
                seq_rows
            )
            # Insert into tetra_data
            placeholders = ', '.join(['?'] * len(tnfs))
            con.executemany(
                f"INSERT INTO tetra_data VALUES ({placeholders})",
                rows
            )

        con.execute(
            "INSERT INTO import_log (filename) VALUES (?) ON CONFLICT DO NOTHING",
            [file]
        )

    con.close()


if __name__ == '__main__':
    main()
