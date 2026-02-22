#!/usr/bin/env python3
"""Convert Bakta annotation TSV + barrnap rRNA + Aragorn tRNA to compact JSON for viz.

Output format: { "contig_1": [{ s, e, d, t, g, p }, ...], ... }
  s = start, e = end, d = strand (+1/-1), t = type, g = gene name, p = product

Usage:
  python genes_to_json.py <bakta_annotation.tsv> <output.json> [rrna_genes.tsv] [trna_genes.tsv]
"""
import json
import os
import sys


def load_bakta(tsv_path):
    """Parse Bakta annotation TSV into feature dict."""
    genes = {}
    n = 0

    with open(tsv_path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 8:
                continue

            contig = parts[0]
            ftype = parts[1].lower()
            if ftype == 'region':
                continue

            try:
                start = int(parts[2])
                end = int(parts[3])
            except ValueError:
                continue

            strand = 1 if parts[4] == '+' else -1
            locus_tag = parts[5] if len(parts) > 5 and parts[5] else ''
            gene = parts[6] if len(parts) > 6 and parts[6] else ''
            product = parts[7] if len(parts) > 7 and parts[7] else ''

            type_map = {
                'cds': 'cds', 'trna': 'trna', 'rrna': 'rrna',
                'tmrna': 'tmrna', 'ncrna': 'ncrna', 'ncrna-region': 'ncrna',
                'crispr': 'crispr', 'crispr-repeat': 'crispr',
                'crispr-spacer': 'crispr',
                'regulatory_region': 'reg', 'oric': 'ori', 'orit': 'ori',
            }
            ftype = type_map.get(ftype, ftype)

            feat = {'s': start, 'e': end, 'd': strand, 't': ftype}
            if locus_tag:
                feat['id'] = locus_tag
            if gene:
                feat['g'] = gene
            if product and product != 'hypothetical protein':
                feat['p'] = product

            genes.setdefault(contig, []).append(feat)
            n += 1

    return genes, n


def load_rrna(tsv_path):
    """Parse barrnap rrna_genes.tsv into feature dict.

    Only includes 'best' hits (one per overlapping region) if the 'best'
    column is present (added by parse_rrna_results.py deduplication).
    Falls back to local dedup when 'best' column is absent.
    """
    raw_genes = {}  # contig -> list of (feat_dict, completeness, identity)
    n = 0
    has_best_col = False

    with open(tsv_path) as f:
        header = f.readline().rstrip('\n').split('\t')
        best_col = header.index('best') if 'best' in header else -1
        has_best_col = best_col >= 0
        compl_col = header.index('gene_completeness') if 'gene_completeness' in header else 8
        ident_col = header.index('vsearch_identity') if 'vsearch_identity' in header else 11
        tax_col = header.index('taxonomy') if 'taxonomy' in header else 13

        for line in f:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 7:
                continue

            # Skip non-best hits if the column exists
            if has_best_col and len(parts) > best_col and parts[best_col] != 'True':
                continue

            contig = parts[1]
            try:
                start = int(parts[2])
                end = int(parts[3])
            except ValueError:
                continue

            strand = 1 if parts[4] == '+' else -1
            rrna_type = parts[5]  # e.g. "16S_rRNA", "23S_rRNA", "5S_rRNA"
            kingdom = parts[6] if len(parts) > 6 else ''

            # Completeness and identity for dedup ranking
            try:
                completeness = float(parts[compl_col]) if len(parts) > compl_col else 0.0
            except ValueError:
                completeness = 0.0
            try:
                identity = float(parts[ident_col]) if len(parts) > ident_col else 0.0
            except ValueError:
                identity = 0.0

            # Map type to compact name
            gene_name = rrna_type.replace('_rRNA', 'S').replace('_', ' ')  # "16S", "23S", "5S"
            product = f"{rrna_type} ({kingdom})" if kingdom else rrna_type

            # Add taxonomy if available
            if len(parts) > tax_col and parts[tax_col]:
                tax = parts[tax_col].split(';')
                deepest = [t for t in tax if t]
                if deepest:
                    product += f" - {deepest[-1]}"

            feat = {'s': start, 'e': end, 'd': strand, 't': 'rrna', 'g': gene_name}
            if product:
                feat['p'] = product

            raw_genes.setdefault(contig, []).append((feat, completeness, identity))
            n += 1

    # Deduplicate overlapping hits if 'best' column was absent
    genes = {}
    n_kept = 0
    if has_best_col:
        # Already filtered by 'best' column above
        for contig, entries in raw_genes.items():
            genes[contig] = [e[0] for e in entries]
            n_kept += len(entries)
    else:
        # Local dedup: for overlapping rRNA on same contig+strand, keep best
        for contig, entries in raw_genes.items():
            # Group by strand
            by_strand = {}
            for feat, compl, ident in entries:
                by_strand.setdefault(feat['d'], []).append((feat, compl, ident))

            kept = []
            for strand_entries in by_strand.values():
                strand_entries.sort(key=lambda x: x[0]['s'])
                clusters = []
                for entry in strand_entries:
                    feat = entry[0]
                    merged = False
                    for cl in clusters:
                        last = cl[-1][0]
                        ov = max(0, min(feat['e'], last['e']) - max(feat['s'], last['s']) + 1)
                        shorter = min(feat['e'] - feat['s'] + 1, last['e'] - last['s'] + 1)
                        if shorter > 0 and ov / shorter > 0.5:
                            cl.append(entry)
                            merged = True
                            break
                    if not merged:
                        clusters.append([entry])
                for cl in clusters:
                    best = max(cl, key=lambda x: (x[1], x[2]))
                    kept.append(best[0])

            genes[contig] = kept
            n_kept += len(kept)

    if not has_best_col and n != n_kept:
        print(f"  rRNA dedup: {n} → {n_kept} (removed {n - n_kept} overlapping hits)")

    return genes, n_kept


def load_trna(tsv_path):
    """Parse Aragorn trna_genes.tsv into feature dict."""
    genes = {}
    n = 0

    with open(tsv_path) as f:
        header = f.readline()
        for line in f:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 7:
                continue

            contig = parts[1]
            try:
                start = int(parts[2])
                end = int(parts[3])
            except ValueError:
                continue

            strand = 1 if parts[4] == '+' else -1
            gene_type = parts[5].lower()  # tRNA or tmRNA
            amino_acid = parts[6] if len(parts) > 6 else ''
            anticodon = parts[7] if len(parts) > 7 else ''

            ftype = 'trna' if gene_type == 'trna' else 'tmrna'
            gene_name = f"tRNA-{amino_acid}" if ftype == 'trna' else 'tmRNA'
            product = f"tRNA-{amino_acid}({anticodon})" if ftype == 'trna' and anticodon else ''

            feat = {'s': start, 'e': end, 'd': strand, 't': ftype, 'g': gene_name}
            if product:
                feat['p'] = product

            genes.setdefault(contig, []).append(feat)
            n += 1

    return genes, n


def main():
    if len(sys.argv) < 3:
        print(f"Usage: {sys.argv[0]} <bakta.tsv> <output.json> [rrna_genes.tsv] [trna_genes.tsv]", file=sys.stderr)
        sys.exit(1)

    tsv_path = sys.argv[1]
    out_path = sys.argv[2]
    rrna_path = sys.argv[3] if len(sys.argv) > 3 else None
    trna_path = sys.argv[4] if len(sys.argv) > 4 else None

    genes, n_bakta = load_bakta(tsv_path)
    print(f"  Bakta: {n_bakta} features from {len(genes)} contigs")

    n_rrna = 0
    if rrna_path and os.path.isfile(rrna_path):
        rrna_genes, n_rrna = load_rrna(rrna_path)
        for contig, feats in rrna_genes.items():
            genes.setdefault(contig, []).extend(feats)
        print(f"  rRNA: {n_rrna} features from {len(rrna_genes)} contigs (best hits only)")

    n_trna = 0
    if trna_path and os.path.isfile(trna_path):
        trna_genes, n_trna = load_trna(trna_path)
        for contig, feats in trna_genes.items():
            genes.setdefault(contig, []).extend(feats)
        print(f"  tRNA/tmRNA: {n_trna} features from {len(trna_genes)} contigs")

    # Sort features by start position within each contig
    for contig in genes:
        genes[contig].sort(key=lambda f: f['s'])

    with open(out_path, 'w') as f:
        json.dump(genes, f, separators=(',', ':'))

    n_contigs = len(genes)
    n_features = n_bakta + n_rrna + n_trna
    size_mb = os.path.getsize(out_path) / 1e6
    print(f"  Wrote {out_path}: {n_contigs} contigs, {n_features} features, {size_mb:.1f} MB")


if __name__ == '__main__':
    main()
