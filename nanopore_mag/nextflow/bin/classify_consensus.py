#!/usr/bin/env python3
"""DEPRECATED: No longer used by the Nextflow pipeline.

Merge Tiara and Whokaryote contig classifications into a consensus.

Consensus logic:
  - If both agree → high confidence
  - If they disagree → use Tiara (higher accuracy, 98%+ vs 95%)
  - Organellar calls from Tiara always kept (Whokaryote has no organellar detection)

Tiara classes:  archaea, bacteria, prokarya, eukarya, organelle, mitochondrion, plastid, unknown
Whokaryote classes: prokaryote, eukaryote

Output files:
  - contig_classifications.tsv: full per-contig classification table
  - prokaryotic_contigs.txt: list of prokaryotic contig IDs
  - eukaryotic_contigs.txt: list of eukaryotic contig IDs
  - organellar_contigs.txt: list of organellar contig IDs (with subdivision)

Usage:
    classify_consensus.py --tiara tiara_output.tsv --whokaryote whokaryote_classifications.tsv
                          --outdir consensus_out/
"""

import argparse
import csv
import os
import sys


# Map Tiara classes to domain
TIARA_TO_DOMAIN = {
    'archaea': 'prokaryotic',
    'bacteria': 'prokaryotic',
    'prokarya': 'prokaryotic',
    'eukarya': 'eukaryotic',
    'organelle': 'organellar',
    'mitochondrion': 'organellar',
    'plastid': 'organellar',
    'unknown': 'unknown',
}

# Map Whokaryote classes to domain
WHOKARYOTE_TO_DOMAIN = {
    'prokaryote': 'prokaryotic',
    'eukaryote': 'eukaryotic',
}


def parse_tiara(path):
    """Parse Tiara output TSV.

    Tiara output format (tab-separated):
      sequence_id  class_fst_stage  class_snd_stage
    First stage: prokarya, eukarya, organelle, unknown
    Second stage: archaea, bacteria, eukarya, mitochondrion, plastid, unknown
    """
    results = {}
    try:
        with open(path) as fh:
            reader = csv.DictReader(fh, delimiter='\t')
            for row in reader:
                seq_id = row.get('sequence_id', '').strip()
                if not seq_id:
                    continue

                # Prefer second-stage classification (more specific)
                cls = row.get('class_snd_stage', '').strip().lower()
                if not cls or cls == 'n/a':
                    cls = row.get('class_fst_stage', '').strip().lower()

                results[seq_id] = cls
    except FileNotFoundError:
        print(f"[WARNING] Tiara file not found: {path}", file=sys.stderr)
    except Exception as e:
        print(f"[WARNING] Error parsing Tiara output: {e}", file=sys.stderr)

    return results


def parse_whokaryote(path):
    """Parse Whokaryote prediction output.

    Whokaryote classifications format (tab-separated):
      contig  prediction  (+ feature columns)
    Predictions: prokaryote, eukaryote
    """
    results = {}
    try:
        with open(path) as fh:
            reader = csv.DictReader(fh, delimiter='\t')
            for row in reader:
                contig = row.get('contig', '').strip()
                if not contig:
                    continue
                pred = row.get('prediction', '').strip().lower()
                results[contig] = pred
    except FileNotFoundError:
        print(f"[WARNING] Whokaryote file not found: {path}", file=sys.stderr)
    except Exception as e:
        print(f"[WARNING] Error parsing Whokaryote output: {e}", file=sys.stderr)

    return results


def consensus_classify(tiara_class, whokaryote_class):
    """Determine consensus classification and confidence.

    Returns (consensus_domain, confidence) tuple.
    """
    tiara_domain = TIARA_TO_DOMAIN.get(tiara_class, 'unknown')
    whokaryote_domain = WHOKARYOTE_TO_DOMAIN.get(whokaryote_class, 'unknown')

    # Organellar calls from Tiara are always kept (Whokaryote can't detect them)
    if tiara_domain == 'organellar':
        return 'organellar', 'high'

    # Both classifiers agree
    if tiara_domain == whokaryote_domain and tiara_domain != 'unknown':
        return tiara_domain, 'high'

    # One is unknown, use the other
    if tiara_domain == 'unknown' and whokaryote_domain != 'unknown':
        return whokaryote_domain, 'low'
    if whokaryote_domain == 'unknown' and tiara_domain != 'unknown':
        return tiara_domain, 'medium'

    # Disagreement: use Tiara (higher accuracy)
    if tiara_domain != 'unknown':
        return tiara_domain, 'medium'

    # Both unknown
    return 'unknown', 'low'


def organellar_subtype(tiara_class):
    """Return organellar subtype from Tiara classification."""
    if tiara_class == 'mitochondrion':
        return 'mitochondrial'
    elif tiara_class == 'plastid':
        return 'plastid'
    elif tiara_class == 'organelle':
        return 'unspecified'
    return ''


def main():
    parser = argparse.ArgumentParser(
        description='Merge Tiara and Whokaryote classifications into consensus')
    parser.add_argument('--tiara', required=True,
                        help='Tiara output TSV')
    parser.add_argument('--whokaryote', required=True,
                        help='Whokaryote classifications TSV')
    parser.add_argument('--outdir', default='.',
                        help='Output directory (default: current directory)')
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    # Parse both classifiers
    tiara_results = parse_tiara(args.tiara)
    whokaryote_results = parse_whokaryote(args.whokaryote)

    print(f"[INFO] Tiara: {len(tiara_results)} contigs classified", file=sys.stderr)
    print(f"[INFO] Whokaryote: {len(whokaryote_results)} contigs classified", file=sys.stderr)

    # Union of all contig IDs
    all_contigs = sorted(set(tiara_results.keys()) | set(whokaryote_results.keys()))
    print(f"[INFO] Total unique contigs: {len(all_contigs)}", file=sys.stderr)

    # Classify each contig
    classifications = []
    prok_ids = []
    euk_ids = []
    org_ids = []

    for contig_id in all_contigs:
        tiara_class = tiara_results.get(contig_id, 'unknown')
        whokaryote_class = whokaryote_results.get(contig_id, 'unknown')
        consensus, confidence = consensus_classify(tiara_class, whokaryote_class)
        subtype = organellar_subtype(tiara_class) if consensus == 'organellar' else ''

        classifications.append({
            'contig_id': contig_id,
            'tiara_class': tiara_class,
            'whokaryote_class': whokaryote_class,
            'consensus_class': consensus,
            'confidence': confidence,
            'organellar_subtype': subtype,
        })

        if consensus == 'prokaryotic':
            prok_ids.append(contig_id)
        elif consensus == 'eukaryotic':
            euk_ids.append(contig_id)
        elif consensus == 'organellar':
            org_ids.append(contig_id)

    # Write main classification table
    tsv_path = os.path.join(args.outdir, 'contig_classifications.tsv')
    with open(tsv_path, 'w', newline='') as fh:
        writer = csv.DictWriter(fh, fieldnames=[
            'contig_id', 'tiara_class', 'whokaryote_class',
            'consensus_class', 'confidence', 'organellar_subtype'
        ], delimiter='\t')
        writer.writeheader()
        writer.writerows(classifications)

    # Write contig ID lists
    for filename, ids in [
        ('prokaryotic_contigs.txt', prok_ids),
        ('eukaryotic_contigs.txt', euk_ids),
        ('organellar_contigs.txt', org_ids),
    ]:
        with open(os.path.join(args.outdir, filename), 'w') as fh:
            for cid in ids:
                fh.write(cid + '\n')

    # Summary stats
    n_total = len(all_contigs)
    print(f"[INFO] Consensus: {len(prok_ids)} prokaryotic, "
          f"{len(euk_ids)} eukaryotic, {len(org_ids)} organellar, "
          f"{n_total - len(prok_ids) - len(euk_ids) - len(org_ids)} unknown "
          f"(of {n_total} total)", file=sys.stderr)


if __name__ == '__main__':
    main()
