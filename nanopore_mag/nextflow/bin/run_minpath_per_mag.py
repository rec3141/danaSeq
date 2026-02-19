#!/usr/bin/env python3
"""Run MinPath parsimony pathway reconstruction on per-MAG annotation TSVs.

For each MAG TSV in the input directory, extracts KO assignments, runs MinPath
to find the minimum set of KEGG pathways consistent with observed KOs, and
collects results into a summary table.

MinPath (Ye & Doak, 2009) prevents pathway inflation in draft MAGs by using
integer programming (GLPK) to find the parsimonious pathway set.

Usage:
    run_minpath_per_mag.py \\
        --input per_mag/ \\
        --minpath_dir /path/to/MinPath \\
        --output minpath_pathways.tsv \\
        --details minpath_details/
"""

import argparse
import csv
import os
import subprocess
import sys
import tempfile


def extract_kos_from_mag(tsv_path):
    """Extract KO assignments from a per-MAG annotation TSV.

    Returns a list of (protein_id, ko) tuples for proteins with KO assignments.
    The KO column may contain 'ko:' prefixed or comma-separated entries.
    """
    results = []
    with open(tsv_path) as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        for row in reader:
            ko_raw = row.get('KO', '').strip()
            if not ko_raw:
                continue
            protein_id = row.get('protein_id', '').strip()
            if not protein_id:
                continue
            # Handle comma-separated or ko:-prefixed KOs
            for k in ko_raw.replace('ko:', '').split(','):
                k = k.strip()
                if k.startswith('K') and len(k) == 6:
                    results.append((protein_id, k))
    return results


def write_minpath_input(ko_pairs, out_path):
    """Write MinPath input file: one line per protein-KO pair.

    Format: protein_id<tab>KO_number
    """
    with open(out_path, 'w') as fh:
        for protein_id, ko in ko_pairs:
            fh.write(f"{protein_id}\t{ko}\n")


def run_minpath(input_path, report_path, details_path, minpath_dir):
    """Run MinPath.py on a KO input file.

    Returns True on success, False on failure.
    """
    minpath_script = os.path.join(minpath_dir, 'MinPath.py')
    if not os.path.isfile(minpath_script):
        print(f"[ERROR] MinPath.py not found at: {minpath_script}", file=sys.stderr)
        return False

    cmd = [
        sys.executable, minpath_script,
        '-ko', input_path,
        '-report', report_path,
        '-details', details_path,
    ]

    try:
        result = subprocess.run(
            cmd,
            capture_output=True, text=True, timeout=600
        )
        if result.returncode != 0:
            print(f"[WARNING] MinPath failed (exit {result.returncode}): "
                  f"{result.stderr[:500]}", file=sys.stderr)
            return False
        return True
    except subprocess.TimeoutExpired:
        print("[WARNING] MinPath timed out after 600s", file=sys.stderr)
        return False
    except Exception as e:
        print(f"[WARNING] MinPath execution error: {e}", file=sys.stderr)
        return False


def parse_minpath_report(report_path):
    """Parse MinPath report file.

    Actual report format (key-value pairs, double-space separated):
      path 00010 kegg n/a  naive 1  minpath 0  fam0  105  fam-found  2  name  Glycolysis...

    Returns list of dicts with pathway info.
    """
    pathways = []
    if not os.path.isfile(report_path):
        return pathways

    with open(report_path) as fh:
        for line in fh:
            line = line.strip()
            if not line or not line.startswith('path'):
                continue

            # Parse key-value pairs from the line
            # Format: path <id> kegg <src>  naive <0|1>  minpath <0|1>  fam0  <N>  fam-found  <N>  name  <name>
            tokens = line.split()
            if len(tokens) < 2:
                continue

            pathway_id = tokens[1]

            # Extract values by finding keywords in the token list
            kv = {}
            i = 0
            while i < len(tokens):
                if tokens[i] == 'naive' and i + 1 < len(tokens):
                    kv['naive'] = tokens[i + 1]
                    i += 2
                elif tokens[i] == 'minpath' and i + 1 < len(tokens):
                    kv['minpath'] = tokens[i + 1]
                    i += 2
                elif tokens[i] == 'fam0' and i + 1 < len(tokens):
                    kv['fam0'] = tokens[i + 1]
                    i += 2
                elif tokens[i] == 'fam-found' and i + 1 < len(tokens):
                    kv['fam-found'] = tokens[i + 1]
                    i += 2
                elif tokens[i] == 'name' and i + 1 < len(tokens):
                    # Everything after 'name' is the pathway name
                    kv['name'] = ' '.join(tokens[i + 1:])
                    break
                else:
                    i += 1

            is_naive = int(kv.get('naive', '0'))
            is_minpath = int(kv.get('minpath', '0'))

            pathways.append({
                'pathway_id': pathway_id,
                'pathway_name': kv.get('name', ''),
                'naive': is_naive,
                'minpath': is_minpath,
                'total_families': kv.get('fam0', '0'),
                'found_families': kv.get('fam-found', '0'),
            })

    return pathways


def main():
    parser = argparse.ArgumentParser(
        description='Run MinPath parsimony pathway reconstruction per MAG')
    parser.add_argument('--input', required=True,
                        help='Directory of per-MAG annotation TSVs')
    parser.add_argument('--minpath_dir', required=True,
                        help='Path to cloned MinPath repository')
    parser.add_argument('--output', required=True,
                        help='Output summary TSV (all MAGs)')
    parser.add_argument('--details', required=True,
                        help='Output directory for per-MAG detail files')
    args = parser.parse_args()

    os.makedirs(args.details, exist_ok=True)

    # Find per-MAG TSVs
    tsv_files = sorted(f for f in os.listdir(args.input) if f.endswith('.tsv'))
    if not tsv_files:
        print("[WARNING] No per-MAG TSV files found", file=sys.stderr)
        with open(args.output, 'w') as fh:
            fh.write('mag_id\tpathway_id\tpathway_name\tnaive\tminpath\t'
                     'total_families\tfound_families\n')
        return

    print(f"[INFO] Processing {len(tsv_files)} MAGs with MinPath", file=sys.stderr)

    all_results = []

    for tsv_file in tsv_files:
        mag_id = tsv_file.replace('.tsv', '')
        tsv_path = os.path.join(args.input, tsv_file)

        # Extract KOs
        ko_pairs = extract_kos_from_mag(tsv_path)
        if not ko_pairs:
            print(f"[INFO] {mag_id}: no KO assignments — skipping", file=sys.stderr)
            continue

        print(f"[INFO] {mag_id}: {len(ko_pairs)} KO assignments", file=sys.stderr)

        # Write MinPath input to temp file
        with tempfile.NamedTemporaryFile(mode='w', suffix='.tsv',
                                         delete=False) as tmp:
            tmp_input = tmp.name
            for protein_id, ko in ko_pairs:
                tmp.write(f"{protein_id}\t{ko}\n")

        report_path = os.path.join(args.details, f"{mag_id}.report.txt")
        details_path = os.path.join(args.details, f"{mag_id}.details.txt")

        try:
            success = run_minpath(tmp_input, report_path, details_path,
                                  args.minpath_dir)
            if success:
                pathways = parse_minpath_report(report_path)
                for pw in pathways:
                    pw['mag_id'] = mag_id
                    all_results.append(pw)
                print(f"[INFO] {mag_id}: {sum(1 for p in pathways if p['minpath'])} "
                      f"minpath / {sum(1 for p in pathways if p['naive'])} naive pathways",
                      file=sys.stderr)
            else:
                print(f"[WARNING] {mag_id}: MinPath failed", file=sys.stderr)
        finally:
            os.unlink(tmp_input)

    # Write summary TSV
    with open(args.output, 'w', newline='') as fh:
        writer = csv.DictWriter(fh, delimiter='\t',
                                fieldnames=['mag_id', 'pathway_id', 'pathway_name',
                                            'naive', 'minpath', 'total_families',
                                            'found_families'])
        writer.writeheader()
        writer.writerows(all_results)

    n_mags = len(set(r['mag_id'] for r in all_results))
    print(f"[INFO] Wrote {len(all_results)} pathway entries for {n_mags} MAGs",
          file=sys.stderr)


if __name__ == '__main__':
    main()
