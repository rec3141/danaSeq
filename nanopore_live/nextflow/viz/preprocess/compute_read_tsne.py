#!/usr/bin/env python3
"""Compute read-level t-SNE from tetramer frequencies across all samples.

Reads tetra_data from per-sample DuckDB files, joins with kraken taxonomy
and stats, computes t-SNE on tetramer vectors, and outputs read_explorer.json
for the WebGL scatter plot.
"""

import argparse
import json
import os
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
                    'id': f"{flowcell}/{barcode_dir.name}",
                    'db_path': str(db_path),
                })
    return samples


def extract_reads(db_path, sample_id, max_per_sample):
    """Extract reads with tetramer, stats, and taxonomy from a single DuckDB."""
    reads = []
    tetra_rows = []
    try:
        con = duckdb.connect(db_path, read_only=True)

        # Check what tables exist
        tables = [r[0] for r in con.execute("SHOW TABLES").fetchall()]

        if 'tetra_data' not in tables or 'stats' not in tables:
            con.close()
            return [], []

        # Get tetramer column names (exclude seqid)
        cols = con.execute("SELECT column_name FROM information_schema.columns WHERE table_name='tetra_data' AND column_name != 'seqid'").fetchall()
        tetra_cols = [c[0] for c in cols]

        if not tetra_cols:
            con.close()
            return [], []

        # Join tetra + stats + kraken
        tetra_col_str = ', '.join(f't."{c}"' for c in tetra_cols)
        has_kraken = 'kraken' in tables

        if has_kraken:
            query = f"""
                SELECT s.seqid, s.length, s.gc,
                       k.taxa_name,
                       {tetra_col_str}
                FROM stats s
                JOIN tetra_data t ON s.seqid = t.seqid
                LEFT JOIN kraken k ON s.seqid = k.seqid
                LIMIT {max_per_sample}
            """
        else:
            query = f"""
                SELECT s.seqid, s.length, s.gc,
                       NULL as taxa_name,
                       {tetra_col_str}
                FROM stats s
                JOIN tetra_data t ON s.seqid = t.seqid
                LIMIT {max_per_sample}
            """

        rows = con.execute(query).fetchall()

        for row in rows:
            seqid = row[0]
            length = row[1]
            gc = row[2]
            taxa_name = row[3]

            # Parse taxonomy
            kraken_phylum = None
            kraken_class = None
            kraken_domain = None
            if taxa_name:
                # Extract taxon name (remove taxid suffix)
                import re
                match = re.match(r'(.+?)(?:\s*\(taxid \d+\))?$', str(taxa_name))
                name = match.group(1).strip() if match else str(taxa_name).strip()
                kraken_phylum = name  # simplified: use full name as phylum

            # Tetramer values
            tetra_values = [float(v) if v is not None else 0.0 for v in row[4:]]

            reads.append({
                'id': seqid,
                'sample': sample_id,
                'length': length,
                'gc': round(gc, 1) if gc else None,
                'kraken_phylum': kraken_phylum,
                'kraken_domain': kraken_domain,
                'kraken_class': kraken_class,
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

    # Handle any NaN/inf
    tetra_matrix = np.nan_to_num(tetra_matrix, nan=0.0, posinf=0.0, neginf=0.0)

    # Compute t-SNE
    perplexity = min(30, max(5, len(all_reads) // 10))
    tsne = TSNE(n_components=2, perplexity=perplexity, random_state=42,
                n_iter=1000, learning_rate='auto', init='pca')
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
