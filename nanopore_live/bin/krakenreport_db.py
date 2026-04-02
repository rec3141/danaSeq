#!/usr/bin/env python3
"""Load Kraken2 report files into DuckDB with reconstructed taxonomy paths.

Usage: python3 41_krakenreport_db.py <barcode_dir>

Reads kraken/*.report, builds full taxonomy lineage from indentation,
inserts into krakenreport table.
"""

import sys
import os
import glob

import duckdb

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from db_schema import ensure_schema


def read_kraken_report(filepath):
    """Parse a kraken2 report file, reconstructing full taxonomy paths."""
    basename = os.path.basename(filepath)
    parts = basename.split('_')
    flowcell = parts[0] if len(parts) > 0 else ''
    barcode = parts[2] if len(parts) > 2 else ''

    taxonomy_stack = []
    rank_stack = []
    rows = []

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

            # Compute indent level from leading spaces
            stripped = raw_name.lstrip(' ')
            n_spaces = len(raw_name) - len(stripped)
            level = n_spaces // 2 + 1
            name = stripped.strip()

            # Maintain taxonomy stack
            if level <= len(taxonomy_stack):
                taxonomy_stack = taxonomy_stack[:level - 1]
                rank_stack = rank_stack[:level - 1]
            taxonomy_stack.append(name)
            rank_stack.append(rank)

            taxonomy = '; '.join(taxonomy_stack)
            rank_full = ' '.join(rank_stack)

            rows.append((
                pct, reads, direct, rank, taxid, name,
                level, flowcell, barcode, taxonomy, rank_full
            ))

    return rows


def main():
    barcode_dir = sys.argv[1]
    os.chdir(barcode_dir)

    con = duckdb.connect('dana.duckdb')
    ensure_schema(con)

    report_files = sorted(glob.glob('kraken/*.report'))
    imported = {r[0] for r in con.execute(
        "SELECT filename FROM import_log WHERE filename LIKE 'kraken%report'"
    ).fetchall()}
    pending = [f for f in report_files if f not in imported]

    for file in pending:
        try:
            rows = read_kraken_report(file)
        except Exception as e:
            print(f"[WARNING] Failed to read {file}: {e}", file=sys.stderr)
            continue

        if rows:
            con.executemany(
                "INSERT INTO krakenreport "
                "(percent, reads, direct, rank, taxid, name, level, "
                "flowcell, barcode, taxonomy, rank_full) "
                "VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
                rows
            )
        con.execute(
            "INSERT INTO import_log (filename) VALUES (?) ON CONFLICT DO NOTHING",
            [file]
        )

    con.close()


if __name__ == '__main__':
    main()
