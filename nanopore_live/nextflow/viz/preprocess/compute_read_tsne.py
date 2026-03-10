#!/usr/bin/env python3
"""Compute read-level t-SNE from tetramer frequencies across all samples.

Reads tetra_data from per-sample DuckDB files, joins with kraken taxonomy
and stats, computes t-SNE on tetramer vectors, and outputs read_explorer.json
for the WebGL scatter plot.
"""

import argparse
import json
import os
import re
import sys
from pathlib import Path

try:
    import duckdb
except ImportError:
    print("[ERROR] duckdb required: pip install duckdb", file=sys.stderr)
    sys.exit(1)

try:
    import numpy as np
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
                    'id': f"{flowcell}_{barcode_dir.name}",
                    'db_path': str(db_path),
                })
    return samples


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
            }
    except Exception:
        pass  # krakenreport may not exist
    return lookup


def extract_reads(db_path, sample_id, max_per_sample):
    """Extract reads with tetramer, stats, and taxonomy from a single DuckDB."""
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
        select_parts.append(tetra_col_str)

        select_str = ', '.join(select_parts)

        # Build FROM/JOIN clause
        from_str = 'tetra_data t'
        if has_stats:
            from_str += ' LEFT JOIN stats s ON t.seqid = s.seqid'
        if has_kraken:
            from_str += ' LEFT JOIN kraken k ON t.seqid = k.seqid'

        query = f"SELECT {select_str} FROM {from_str} LIMIT {max_per_sample}"
        rows = con.execute(query).fetchall()

        for row in rows:
            seqid = row[0]
            length = row[1]
            gc = row[2]
            taxa_name = row[3]
            tetra_values = [float(v) if v is not None else 0.0 for v in row[4:]]

            # Estimate GC from tetramer frequencies if not available
            if gc is None and tetra_values:
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

            reads.append({
                'id': seqid,
                'sample': sample_id,
                'length': length,
                'gc': gc,
                'kraken_domain': lineage.get('domain'),
                'kraken_phylum': lineage.get('phylum'),
                'kraken_class':  lineage.get('class'),
                'kraken_order':  lineage.get('order'),
                'kraken_family': lineage.get('family'),
                'kraken_genus':  lineage.get('genus'),
            })
            tetra_rows.append(tetra_values)

        con.close()
    except Exception as e:
        print(f"  [WARNING] Failed to extract reads from {db_path}: {e}", file=sys.stderr)

    return reads, tetra_rows


def main():
    parser = argparse.ArgumentParser(description='Compute read-level t-SNE')
    parser.add_argument('--input', required=True, help='Directory with nanopore_live outputs')
    parser.add_argument('--output', required=True, help='Output directory')
    parser.add_argument('--max-reads', type=int, default=200000, help='Maximum total reads for t-SNE')
    args = parser.parse_args()

    samples = discover_samples(args.input)
    print(f"[INFO] Found {len(samples)} samples", file=sys.stderr)

    # Calculate per-sample budget
    max_per_sample = max(100, args.max_reads // max(1, len(samples)))
    print(f"[INFO] Max {max_per_sample} reads per sample, {args.max_reads} total", file=sys.stderr)

    # Extract reads from all samples
    all_reads = []
    all_tetra = []

    for i, s in enumerate(samples):
        print(f"  [{i+1}/{len(samples)}] {s['id']}...", file=sys.stderr, end='', flush=True)
        reads, tetra = extract_reads(s['db_path'], s['id'], max_per_sample)
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

    # Compute t-SNE
    perplexity = min(30, max(5, len(all_reads) // 10))
    tsne = TSNE(n_components=2, perplexity=perplexity, random_state=42,
                max_iter=1000, learning_rate='auto', init='pca')
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
