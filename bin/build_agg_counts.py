#!/usr/bin/env python3
"""Generate aggregated_counts.json.gz from existing viz data files.

Usage: build_agg_counts.py <viz_dir>
"""
import sys, json, gzip, os
from collections import defaultdict

def load_json(path):
    if path.endswith('.gz') or os.path.exists(path + '.gz'):
        p = path if path.endswith('.gz') else path + '.gz'
        with gzip.open(p, 'rt') as f:
            return json.load(f)
    with open(path) as f:
        return json.load(f)

def main():
    viz_dir = sys.argv[1]

    print("Loading data...")
    counts = load_json(os.path.join(viz_dir, 'counts.json.gz'))
    taxonomy = load_json(os.path.join(viz_dir, 'taxonomy.json'))
    asvs = load_json(os.path.join(viz_dir, 'asvs.json.gz'))

    count_data = counts['data']       # [[si, ai, count, prop], ...]
    count_samples = counts['samples'] # sample IDs by index
    count_asvs = counts['asvs']       # ASV IDs by index

    print(f"  {len(count_data):,} count entries, {len(count_samples):,} samples, {len(count_asvs):,} ASVs")

    # Get taxonomy info
    db = next(iter(taxonomy), None)
    if not db:
        print("No taxonomy data found")
        return

    tax_data = taxonomy[db]
    levels = tax_data.get('levels', [])
    assignments = tax_data.get('assignments', {})
    print(f"  Taxonomy: {db}, levels={levels}")

    # Build ASV group lookup
    asv_group = {}
    for a in asvs:
        asv_group[a['id']] = a.get('group', 'prokaryote')

    # Pre-compute sample totals
    sample_totals = defaultdict(int)
    for si, ai, count, prop in count_data:
        sample_totals[si] += count

    result = {}

    for level_idx, level_name in enumerate(levels):
        print(f"  Aggregating {level_name}...")

        # Aggregate: (sample_idx, taxon) -> total count
        agg = defaultdict(int)
        for si, ai, count, prop in count_data:
            asv_id = count_asvs[ai]
            tax_arr = assignments.get(asv_id)
            taxon = (tax_arr[level_idx] if tax_arr and level_idx < len(tax_arr) and tax_arr[level_idx] else 'unclassified')
            agg[(si, taxon)] += count

        # Build output
        taxa = sorted(set(t for _, t in agg.keys()))
        taxon_idx = {t: i for i, t in enumerate(taxa)}

        data = []
        for (si, taxon), cnt in agg.items():
            if cnt == 0:
                continue
            ti = taxon_idx[taxon]
            total = sample_totals[si]
            prop = round(cnt / total, 6) if total > 0 else 0.0
            data.append([si, ti, cnt, prop])

        result[level_name] = {
            'data': data,
            'samples': count_samples,
            'taxa': taxa,
        }
        print(f"    {len(taxa)} taxa, {len(data):,} entries")

    # Add 'group' level
    print("  Aggregating group...")
    agg = defaultdict(int)
    for si, ai, count, prop in count_data:
        asv_id = count_asvs[ai]
        group = asv_group.get(asv_id, 'prokaryote')
        agg[(si, group)] += count

    groups = sorted(set(g for _, g in agg.keys()))
    group_idx = {g: i for i, g in enumerate(groups)}
    data = []
    for (si, group), cnt in agg.items():
        if cnt == 0:
            continue
        ti = group_idx[group]
        total = sample_totals[si]
        prop = round(cnt / total, 6) if total > 0 else 0.0
        data.append([si, ti, cnt, prop])

    result['group'] = {
        'data': data,
        'samples': count_samples,
        'taxa': groups,
    }
    print(f"    {len(groups)} groups, {len(data):,} entries")

    # Write output
    out_path = os.path.join(viz_dir, 'aggregated_counts.json.gz')
    print(f"Writing {out_path}...")
    with gzip.open(out_path, 'wt') as f:
        json.dump(result, f)

    size_mb = os.path.getsize(out_path) / 1e6
    print(f"Done! {size_mb:.1f} MB")

if __name__ == '__main__':
    main()
