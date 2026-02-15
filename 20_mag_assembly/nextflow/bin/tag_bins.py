#!/usr/bin/env python3
"""Tag bins as prokaryotic, eukaryotic, mixed, or organellar based on contig classifications.

For each binner's output, counts prokaryotic vs eukaryotic vs organellar contigs per bin,
calculates domain fractions by both contig count and base pair length, and assigns a tag.

Tagging rules:
  - "prokaryotic": >90% prokaryotic (by bp)
  - "eukaryotic":  >50% eukaryotic (by bp)
  - "organellar":  >50% organellar (by bp)
  - "mixed":       everything else

Input:
  - contig_classifications.tsv: consensus classification from classify_consensus.py
  - One or more contig2bin TSV files (contig_id<TAB>bin_id format, DAS_Tool compatible)
  - Contigs FASTA: for calculating sequence lengths

Output:
  - Per-binner: {binner}_bin_domain_tags.tsv
  - Merged: all_bin_domain_tags.tsv

Usage:
    tag_bins.py --classifications contig_classifications.tsv \
                --contigs assembly.fasta \
                --bins metabat:metabat_bins.tsv,semibin:semibin_bins.tsv,dastool:dastool_contig2bin.tsv \
                --outdir bin_tags/
"""

import argparse
import csv
import os
import sys


def parse_fasta_lengths(path):
    """Parse FASTA file and return dict of {seq_id: length_bp}."""
    lengths = {}
    current_id = None
    current_len = 0

    try:
        with open(path) as fh:
            for line in fh:
                line = line.strip()
                if line.startswith('>'):
                    if current_id is not None:
                        lengths[current_id] = current_len
                    current_id = line[1:].split()[0]
                    current_len = 0
                else:
                    current_len += len(line)
            if current_id is not None:
                lengths[current_id] = current_len
    except FileNotFoundError:
        print(f"[WARNING] FASTA file not found: {path}", file=sys.stderr)
    except Exception as e:
        print(f"[WARNING] Error parsing FASTA: {e}", file=sys.stderr)

    return lengths


def parse_classifications(path):
    """Parse contig_classifications.tsv → {contig_id: consensus_class}."""
    results = {}
    try:
        with open(path) as fh:
            reader = csv.DictReader(fh, delimiter='\t')
            for row in reader:
                contig_id = row.get('contig_id', '').strip()
                consensus = row.get('consensus_class', 'unknown').strip()
                if contig_id:
                    results[contig_id] = consensus
    except FileNotFoundError:
        print(f"[WARNING] Classifications file not found: {path}", file=sys.stderr)
    except Exception as e:
        print(f"[WARNING] Error parsing classifications: {e}", file=sys.stderr)

    return results


def parse_contig2bin(path):
    """Parse DAS_Tool-format contig2bin TSV → {contig_id: bin_id}."""
    results = {}
    try:
        with open(path) as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                parts = line.split('\t')
                if len(parts) >= 2:
                    results[parts[0]] = parts[1]
    except FileNotFoundError:
        print(f"[WARNING] Contig2bin file not found: {path}", file=sys.stderr)
    except Exception as e:
        print(f"[WARNING] Error parsing contig2bin: {e}", file=sys.stderr)

    return results


def tag_bin(prok_bp, euk_bp, org_bp, total_bp):
    """Assign domain tag based on base pair fractions."""
    if total_bp == 0:
        return 'unknown'

    prok_frac = prok_bp / total_bp
    euk_frac = euk_bp / total_bp
    org_frac = org_bp / total_bp

    if prok_frac > 0.9:
        return 'prokaryotic'
    elif euk_frac > 0.5:
        return 'eukaryotic'
    elif org_frac > 0.5:
        return 'organellar'
    else:
        return 'mixed'


def process_binner(binner_name, contig2bin, classifications, contig_lengths):
    """Process one binner's results and return per-bin domain tags."""
    # Group contigs by bin
    bin_contigs = {}
    for contig_id, bin_id in contig2bin.items():
        bin_contigs.setdefault(bin_id, []).append(contig_id)

    rows = []
    for bin_id in sorted(bin_contigs.keys()):
        contigs = bin_contigs[bin_id]
        n_contigs = len(contigs)

        prok_bp = euk_bp = org_bp = unknown_bp = 0
        prok_n = euk_n = org_n = unknown_n = 0

        for contig_id in contigs:
            cls = classifications.get(contig_id, 'unknown')
            bp = contig_lengths.get(contig_id, 0)

            if cls == 'prokaryotic':
                prok_bp += bp
                prok_n += 1
            elif cls == 'eukaryotic':
                euk_bp += bp
                euk_n += 1
            elif cls == 'organellar':
                org_bp += bp
                org_n += 1
            else:
                unknown_bp += bp
                unknown_n += 1

        total_bp = prok_bp + euk_bp + org_bp + unknown_bp
        domain_tag = tag_bin(prok_bp, euk_bp, org_bp, total_bp)

        rows.append({
            'binner': binner_name,
            'bin_id': bin_id,
            'n_contigs': n_contigs,
            'total_bp': total_bp,
            'prok_bp': prok_bp,
            'euk_bp': euk_bp,
            'org_bp': org_bp,
            'unknown_bp': unknown_bp,
            'prok_frac': f"{prok_bp / total_bp:.3f}" if total_bp > 0 else '0.000',
            'euk_frac': f"{euk_bp / total_bp:.3f}" if total_bp > 0 else '0.000',
            'org_frac': f"{org_bp / total_bp:.3f}" if total_bp > 0 else '0.000',
            'tag': domain_tag,
        })

    return rows


FIELDNAMES = [
    'binner', 'bin_id', 'n_contigs', 'total_bp',
    'prok_bp', 'euk_bp', 'org_bp', 'unknown_bp',
    'prok_frac', 'euk_frac', 'org_frac', 'tag',
]


def main():
    parser = argparse.ArgumentParser(
        description='Tag bins as prokaryotic/eukaryotic/mixed/organellar')
    parser.add_argument('--classifications', required=True,
                        help='Consensus contig_classifications.tsv')
    parser.add_argument('--contigs', required=True,
                        help='Assembly contigs FASTA (for sequence lengths)')
    parser.add_argument('--bins', required=True,
                        help='Comma-separated binner:path pairs (e.g. metabat:file.tsv,semibin:file.tsv)')
    parser.add_argument('--outdir', default='.',
                        help='Output directory (default: current directory)')
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    # Parse inputs
    contig_lengths = parse_fasta_lengths(args.contigs)
    classifications = parse_classifications(args.classifications)
    print(f"[INFO] Loaded {len(contig_lengths)} contig lengths, "
          f"{len(classifications)} classifications", file=sys.stderr)

    # Process each binner
    all_rows = []
    for binner_spec in args.bins.split(','):
        binner_spec = binner_spec.strip()
        if not binner_spec:
            continue
        if ':' not in binner_spec:
            print(f"[WARNING] Invalid binner spec (expected name:path): {binner_spec}",
                  file=sys.stderr)
            continue

        binner_name, bin_path = binner_spec.split(':', 1)
        if not os.path.exists(bin_path):
            print(f"[WARNING] Binner file not found: {bin_path}", file=sys.stderr)
            continue

        contig2bin = parse_contig2bin(bin_path)
        if not contig2bin:
            print(f"[INFO] {binner_name}: no bins found (empty contig2bin)", file=sys.stderr)
            continue

        rows = process_binner(binner_name, contig2bin, classifications, contig_lengths)
        all_rows.extend(rows)

        # Write per-binner file
        per_binner_path = os.path.join(args.outdir, f"{binner_name}_bin_domain_tags.tsv")
        with open(per_binner_path, 'w', newline='') as fh:
            writer = csv.DictWriter(fh, fieldnames=FIELDNAMES, delimiter='\t')
            writer.writeheader()
            writer.writerows(rows)

        n_prok = sum(1 for r in rows if r['tag'] == 'prokaryotic')
        n_euk = sum(1 for r in rows if r['tag'] == 'eukaryotic')
        n_mix = sum(1 for r in rows if r['tag'] == 'mixed')
        n_org = sum(1 for r in rows if r['tag'] == 'organellar')
        print(f"[INFO] {binner_name}: {len(rows)} bins — "
              f"{n_prok} prok, {n_euk} euk, {n_org} org, {n_mix} mixed",
              file=sys.stderr)

    # Write merged file
    merged_path = os.path.join(args.outdir, 'all_bin_domain_tags.tsv')
    with open(merged_path, 'w', newline='') as fh:
        writer = csv.DictWriter(fh, fieldnames=FIELDNAMES, delimiter='\t')
        writer.writeheader()
        writer.writerows(all_rows)

    print(f"[INFO] Total: {len(all_rows)} bins tagged across all binners", file=sys.stderr)


if __name__ == '__main__':
    main()
