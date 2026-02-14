#!/usr/bin/env python3
"""Map merged protein annotations to individual MAGs via DAS_Tool contig2bin.tsv.

Uses GFF to build gene->contig mapping, then contig->bin from DAS_Tool.
Outputs per-MAG annotation TSVs and a community-wide table.

Usage:
    map_annotations_to_bins.py \
        --annotations merged_annotations.tsv \
        --contig2bin contig2bin.tsv \
        --gff annotation.gff \
        --outdir per_mag/ \
        --community community_annotations.tsv
"""

import argparse
import csv
import os
import re
import sys


def parse_contig2bin(path):
    """Parse DAS_Tool contig2bin.tsv: contig_id -> bin_id."""
    mapping = {}
    with open(path) as fh:
        for line in fh:
            line = line.rstrip('\n')
            if not line or line.startswith('#'):
                continue
            parts = line.split('\t')
            if len(parts) >= 2:
                mapping[parts[0]] = parts[1]
    return mapping


def parse_gff_gene_contigs(path):
    """Parse GFF to build protein_id -> contig_id mapping from CDS features.

    Extracts the ID attribute from CDS features. Prokka uses ID=<locus_tag>_N,
    Bakta uses ID=<contig>_<N>. The parent contig is column 0 (seqid).
    """
    mapping = {}
    if not path:
        return mapping

    try:
        with open(path) as fh:
            for line in fh:
                if line.startswith('#'):
                    continue
                line = line.rstrip('\n')
                if not line:
                    continue
                parts = line.split('\t')
                if len(parts) < 9:
                    continue

                feature_type = parts[2]
                if feature_type != 'CDS':
                    continue

                contig_id = parts[0]
                attributes = parts[8]

                # Extract ID from attributes
                id_match = re.search(r'ID=([^;]+)', attributes)
                if id_match:
                    gene_id = id_match.group(1)
                    mapping[gene_id] = contig_id

                # Also try locus_tag (Prokka uses this as protein name in .faa)
                lt_match = re.search(r'locus_tag=([^;]+)', attributes)
                if lt_match:
                    locus_tag = lt_match.group(1)
                    mapping[locus_tag] = contig_id
    except FileNotFoundError:
        print(f"[WARNING] GFF file not found: {path}", file=sys.stderr)
    except Exception as e:
        print(f"[WARNING] Error parsing GFF: {e}", file=sys.stderr)

    return mapping


def main():
    parser = argparse.ArgumentParser(
        description='Map merged annotations to MAG bins')
    parser.add_argument('--annotations', required=True,
                        help='Merged annotations TSV (from merge_annotations.py)')
    parser.add_argument('--contig2bin', required=True,
                        help='DAS_Tool contig2bin.tsv')
    parser.add_argument('--gff', help='Prokka/Bakta GFF file (for gene->contig mapping)')
    parser.add_argument('--outdir', required=True,
                        help='Output directory for per-MAG TSVs')
    parser.add_argument('--community', required=True,
                        help='Output path for community-wide TSV')
    args = parser.parse_args()

    # Parse contig-to-bin assignments
    contig2bin = parse_contig2bin(args.contig2bin)
    print(f"[INFO] Loaded {len(contig2bin)} contig-to-bin assignments", file=sys.stderr)

    # Parse GFF for gene->contig mapping (supplements contig_id column in merged TSV)
    gff_gene_map = parse_gff_gene_contigs(args.gff)
    print(f"[INFO] GFF gene-to-contig map: {len(gff_gene_map)} entries", file=sys.stderr)

    os.makedirs(args.outdir, exist_ok=True)

    # Read merged annotations and assign to bins
    bin_rows = {}  # bin_id -> list of rows
    unbinned = 0
    total = 0
    header = None

    with open(args.annotations) as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        header = reader.fieldnames

        for row in reader:
            total += 1
            protein_id = row['protein_id']
            contig_id = row.get('contig_id', '')

            # Try GFF mapping first (more accurate), fall back to annotation column
            if protein_id in gff_gene_map:
                contig_id = gff_gene_map[protein_id]
                row['contig_id'] = contig_id

            bin_id = contig2bin.get(contig_id, '')
            row['bin_id'] = bin_id

            if bin_id:
                bin_rows.setdefault(bin_id, []).append(row)
            else:
                unbinned += 1

    out_header = header + ['bin_id']

    # Write per-MAG TSVs
    for bin_id, rows in sorted(bin_rows.items()):
        out_path = os.path.join(args.outdir, f"{bin_id}.tsv")
        with open(out_path, 'w', newline='') as fh:
            writer = csv.DictWriter(fh, fieldnames=out_header, delimiter='\t')
            writer.writeheader()
            writer.writerows(rows)

    # Write community-wide TSV (all proteins with bin_id column)
    with open(args.community, 'w', newline='') as fh:
        # Re-read to include unbinned proteins too
        pass

    # Re-read and write community file with all proteins
    with open(args.annotations) as in_fh, \
         open(args.community, 'w', newline='') as out_fh:
        reader = csv.DictReader(in_fh, delimiter='\t')
        writer = csv.DictWriter(out_fh, fieldnames=out_header, delimiter='\t')
        writer.writeheader()

        for row in reader:
            protein_id = row['protein_id']
            contig_id = row.get('contig_id', '')

            if protein_id in gff_gene_map:
                contig_id = gff_gene_map[protein_id]
                row['contig_id'] = contig_id

            row['bin_id'] = contig2bin.get(contig_id, '')
            writer.writerow(row)

    binned_proteins = sum(len(rows) for rows in bin_rows.values())
    print(f"[INFO] Total proteins: {total}", file=sys.stderr)
    print(f"[INFO] Binned: {binned_proteins} across {len(bin_rows)} MAGs", file=sys.stderr)
    print(f"[INFO] Unbinned: {unbinned}", file=sys.stderr)
    print(f"[INFO] Per-MAG TSVs written to: {args.outdir}", file=sys.stderr)
    print(f"[INFO] Community TSV written to: {args.community}", file=sys.stderr)


if __name__ == '__main__':
    main()
