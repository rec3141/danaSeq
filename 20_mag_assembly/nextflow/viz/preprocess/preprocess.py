#!/usr/bin/env python3
"""
Preprocess MAG pipeline TSV outputs into JSON for the web dashboard.

Generates 12 JSON files from Nextflow pipeline results:
  overview.json, mags.json, checkm2_all.json, taxonomy_sunburst.json,
  kegg_heatmap.json, coverage.json, mge_summary.json, mge_per_bin.json,
  eukaryotic.json, contig_explorer.json, contig_lengths.json

Usage:
  python3 preprocess.py --results <path> --output <path> [--skip-tsne]
"""

import argparse
import json
import os
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import pdist


def load_tsv(path, **kwargs):
    """Load a TSV file, returning None if it doesn't exist."""
    if not os.path.exists(path):
        print(f"  [WARNING] Missing: {path}", file=sys.stderr)
        return None
    return pd.read_csv(path, sep='\t', **kwargs)


def parse_args():
    p = argparse.ArgumentParser(description="Preprocess MAG pipeline results to JSON")
    p.add_argument('--results', '-r', required=True, help='Pipeline results directory')
    p.add_argument('--output', '-o', required=True, help='Output directory for JSON files')
    p.add_argument('--skip-tsne', action='store_true', help='Skip t-SNE (saves ~5 min)')
    return p.parse_args()


def build_overview(results_dir, assembly_info, depths_df, dastool_summary, checkm2_df,
                   virus_df, plasmid_df, contig2bin):
    """Build overview.json with assembly stats and summary counts."""
    print("Building overview.json ...")

    # Assembly stats from assembly_info
    lengths = assembly_info['length'].values
    total_size = int(lengths.sum())
    n_contigs = len(lengths)

    # Compute N50
    sorted_lens = np.sort(lengths)[::-1]
    cumsum = np.cumsum(sorted_lens)
    n50 = int(sorted_lens[np.searchsorted(cumsum, total_size / 2)])

    # MAG counts by quality tier: prefer CheckM2, fall back to SCG
    n_mags = len(dastool_summary)
    hq = mq = lq = 0
    checkm2_lookup = {}
    if checkm2_df is not None:
        for _, row in checkm2_df.iterrows():
            checkm2_lookup[row['Name']] = row

    for _, ds_row in dastool_summary.iterrows():
        mag_name = ds_row['bin']
        short_name = mag_name.replace('dastool-', '')
        cm_row = checkm2_lookup.get(short_name)
        if cm_row is not None:
            comp = cm_row['Completeness']
            cont = cm_row['Contamination']
        else:
            comp = ds_row.get('SCG_completeness', 0)
            cont = ds_row.get('SCG_redundancy', 0)
        if comp >= 90 and cont < 5:
            hq += 1
        elif comp >= 50 and cont < 10:
            mq += 1
        else:
            lq += 1

    # Binner origin counts
    binner_counts = Counter()
    for b in dastool_summary['bin'].values:
        # dastool-{binner}_{num} or dastool-{binner}_{num}_sub
        m = re.match(r'dastool-(\w+?)_\d+', b)
        if m:
            binner_counts[m.group(1)] += 1
        else:
            binner_counts['other'] += 1

    overview = {
        'assembly_size': total_size,
        'n50': n50,
        'n_contigs': n_contigs,
        'n_mags': n_mags,
        'hq': hq,
        'mq': mq,
        'lq': lq,
        'n_virus': len(virus_df) if virus_df is not None else 0,
        'n_plasmid': len(plasmid_df) if plasmid_df is not None else 0,
        'binner_counts': dict(binner_counts),
    }
    return overview


def build_contig_lengths(assembly_info):
    """Build contig_lengths.json with pre-bucketed histogram data."""
    print("Building contig_lengths.json ...")
    lengths = assembly_info['length'].values

    # Log-spaced bins from min to max
    min_len = max(int(lengths.min()), 1)
    max_len = int(lengths.max())
    bins = np.logspace(np.log10(min_len), np.log10(max_len + 1), num=50)
    counts, edges = np.histogram(lengths, bins=bins)

    return {
        'bin_edges': [int(e) for e in edges],
        'counts': [int(c) for c in counts],
        'n50': int(np.sort(lengths)[::-1][np.searchsorted(
            np.cumsum(np.sort(lengths)[::-1]), lengths.sum() / 2)]),
    }


def build_mags(dastool_summary, checkm2_df, contig2bin, kaiju_df, depths_df,
               virus_df, plasmid_df, defense_df, integron_df):
    """Build mags.json with per-MAG joined records."""
    print("Building mags.json ...")

    mags = []
    dastool_names = set(dastool_summary['bin'].values)

    # Build contig-to-MAG lookup
    c2b = {}
    if contig2bin is not None:
        for _, row in contig2bin.iterrows():
            c2b[row.iloc[0]] = row.iloc[1]

    # Build per-MAG contig lists
    mag_contigs = defaultdict(list)
    for contig, mag in c2b.items():
        mag_contigs[mag].append(contig)

    # Parse sample names from depths header
    depth_cols = [c for c in depths_df.columns if c.endswith('.sorted.bam') and not c.endswith('-var')]
    sample_names = [re.sub(r'\.sorted\.bam$', '', c) for c in depth_cols]

    # Build virus/plasmid/defense/integron sets by contig
    virus_contigs = set(virus_df['seq_name'].values) if virus_df is not None else set()
    plasmid_contigs = set(plasmid_df['seq_name'].values) if plasmid_df is not None else set()

    defense_contig_map = defaultdict(list)
    if defense_df is not None:
        for _, row in defense_df.iterrows():
            # sys_id format: contig_NNNNN_TYPE_N
            sid = row['sys_id']
            parts = sid.split('_')
            # Reconstruct contig name: everything before the type/number suffix
            # Pattern: contig_12345_RM_Type_I_1 → contig_12345
            contig = '_'.join(parts[:2])  # contig_NNNNN
            defense_contig_map[contig].append(row['type'])

    integron_contigs = set()
    if integron_df is not None:
        for _, row in integron_df.iterrows():
            integron_contigs.add(row['ID_replicon'])

    # Build per-MAG taxonomy by majority vote from kaiju
    mag_taxonomy = {}
    if kaiju_df is not None:
        for mag_name, contigs in mag_contigs.items():
            lineages = []
            for c in contigs:
                match = kaiju_df[kaiju_df['contig_id'] == c]
                if len(match) > 0 and pd.notna(match.iloc[0].get('lineage', None)):
                    lineages.append(match.iloc[0]['lineage'])
            if lineages:
                # Majority vote at each rank
                ranks = ['domain', 'phylum', 'class', 'order', 'family', 'genus', 'species']
                tax = {}
                for i, rank in enumerate(ranks):
                    rank_values = []
                    for lin in lineages:
                        parts = [p.strip() for p in lin.split(';')]
                        if i < len(parts) and parts[i]:
                            rank_values.append(parts[i])
                    if rank_values:
                        tax[rank] = Counter(rank_values).most_common(1)[0][0]
                    else:
                        tax[rank] = ''
                mag_taxonomy[mag_name] = tax

    for _, row in dastool_summary.iterrows():
        mag_name = row['bin']
        contigs = mag_contigs.get(mag_name, [])

        # CheckM2 quality
        short_name = mag_name.replace('dastool-', '')
        checkm2_row = None
        if checkm2_df is not None:
            matches = checkm2_df[checkm2_df['Name'] == short_name]
            if len(matches) > 0:
                checkm2_row = matches.iloc[0]

        completeness = float(checkm2_row['Completeness']) if checkm2_row is not None else float(row.get('SCG_completeness', 0))
        contamination = float(checkm2_row['Contamination']) if checkm2_row is not None else float(row.get('SCG_redundancy', 0))
        gc = float(checkm2_row['GC_Content']) if checkm2_row is not None else 0

        # Quality tier
        if completeness >= 90 and contamination < 5:
            quality = 'HQ'
        elif completeness >= 50 and contamination < 10:
            quality = 'MQ'
        else:
            quality = 'LQ'

        # Binner origin
        m = re.match(r'dastool-(\w+?)_\d+', mag_name)
        binner = m.group(1) if m else 'unknown'

        # Per-sample coverage
        coverage = {}
        if contigs and len(depth_cols) > 0:
            contig_depths = depths_df[depths_df['contigName'].isin(contigs)]
            for col, sname in zip(depth_cols, sample_names):
                coverage[sname] = round(float(contig_depths[col].mean()), 4) if len(contig_depths) > 0 else 0

        # MGE counts
        n_virus = sum(1 for c in contigs if c in virus_contigs)
        n_plasmid = sum(1 for c in contigs if c in plasmid_contigs)
        n_defense = sum(len(defense_contig_map.get(c, [])) for c in contigs)
        n_integron = sum(1 for c in contigs if c in integron_contigs)

        tax = mag_taxonomy.get(mag_name, {})

        mags.append({
            'name': mag_name,
            'binner': binner,
            'quality': quality,
            'completeness': round(completeness, 2),
            'contamination': round(contamination, 2),
            'gc': round(gc, 4),
            'size': int(row['size']),
            'n_contigs': int(row['contigs']),
            'n50': int(row['N50']),
            'bin_score': round(float(row['bin_score']), 4),
            'taxonomy': tax,
            'coverage': coverage,
            'n_virus': n_virus,
            'n_plasmid': n_plasmid,
            'n_defense': n_defense,
            'n_integron': n_integron,
        })

    return mags


def build_checkm2_all(checkm2_df, dastool_summary):
    """Build checkm2_all.json for the scatter plot of all 573 bins."""
    print("Building checkm2_all.json ...")
    if checkm2_df is None:
        return []

    dastool_names = set(dastool_summary['bin'].values)
    records = []
    for _, row in checkm2_df.iterrows():
        name = row['Name']
        full_name = f"dastool-{name}"
        is_dastool = full_name in dastool_names

        # Determine binner
        m = re.match(r'(\w+?)_\d+', name)
        binner = m.group(1) if m else 'unknown'

        records.append({
            'name': name,
            'completeness': round(float(row['Completeness']), 2),
            'contamination': round(float(row['Contamination']), 2),
            'genome_size': int(row['Genome_Size']),
            'gc': round(float(row['GC_Content']), 4),
            'n50': int(row['Contig_N50']),
            'binner': binner,
            'is_dastool': is_dastool,
            'dastool_name': full_name if is_dastool else None,
        })
    return records


def build_taxonomy_sunburst(kaiju_df, assembly_info):
    """Build taxonomy_sunburst.json as a D3 hierarchy tree."""
    print("Building taxonomy_sunburst.json ...")
    if kaiju_df is None:
        return {'name': 'root', 'children': []}

    # Build contig length lookup
    len_map = dict(zip(assembly_info['#seq_name'], assembly_info['length']))

    # Count size-weighted taxonomy at each rank
    tree = {}
    ranks = ['domain', 'phylum', 'class', 'order']

    for _, row in kaiju_df.iterrows():
        lineage = row.get('lineage', '')
        if pd.isna(lineage) or not lineage:
            continue
        parts = [p.strip() for p in lineage.split(';')]
        contig_len = len_map.get(row['contig_id'], 1000)

        node = tree
        for i, rank_name in enumerate(ranks):
            if i < len(parts) and parts[i]:
                taxon = parts[i]
            else:
                taxon = 'Unclassified'
            if taxon not in node:
                node[taxon] = {'_size': 0, '_children': {}}
            node[taxon]['_size'] += contig_len
            node = node[taxon]['_children']

    # Convert to D3 hierarchy, collapsing <0.5% nodes
    total_size = sum(v['_size'] for v in tree.values())
    threshold = total_size * 0.005

    def to_hierarchy(node_dict, parent_name='root'):
        children = []
        other_size = 0
        for name, data in sorted(node_dict.items(), key=lambda x: -x[1]['_size']):
            if data['_size'] < threshold:
                other_size += data['_size']
            else:
                child = {'name': name, 'value': data['_size']}
                if data['_children']:
                    sub = to_hierarchy(data['_children'], name)
                    if sub:
                        child['children'] = sub
                children.append(child)
        if other_size > 0:
            children.append({'name': 'Other', 'value': other_size})
        return children

    return {
        'name': 'Life',
        'children': to_hierarchy(tree),
    }


def build_kegg_heatmap(results_dir):
    """Build kegg_heatmap.json with Ward-clustered module completeness matrix."""
    print("Building kegg_heatmap.json ...")
    path = os.path.join(results_dir, 'metabolism', 'modules', 'module_completeness.tsv')
    df = load_tsv(path)
    if df is None:
        return {'mag_ids': [], 'module_ids': [], 'module_names': [], 'matrix': [],
                'row_order': [], 'col_order': []}

    # Filter out _community and _unbinned rows
    df = df[~df['mag_id'].str.startswith('_')].copy()
    df = df.set_index('mag_id')

    # Parse module IDs and names from column headers like "M00001:Glycolysis ..."
    module_ids = []
    module_names = []
    for col in df.columns:
        parts = col.split(':', 1)
        module_ids.append(parts[0])
        module_names.append(parts[1] if len(parts) > 1 else parts[0])

    matrix = df.values.astype(float)

    # Ward hierarchical clustering on rows (MAGs) and columns (modules)
    # Only cluster if we have enough data
    if matrix.shape[0] > 2:
        row_dist = pdist(matrix, metric='euclidean')
        row_linkage = linkage(row_dist, method='ward')
        row_order = leaves_list(row_linkage).tolist()
    else:
        row_order = list(range(matrix.shape[0]))

    if matrix.shape[1] > 2:
        col_dist = pdist(matrix.T, metric='euclidean')
        col_linkage = linkage(col_dist, method='ward')
        col_order = leaves_list(col_linkage).tolist()
    else:
        col_order = list(range(matrix.shape[1]))

    # Reorder
    reordered = matrix[row_order][:, col_order].tolist()
    reordered_mags = [df.index[i] for i in row_order]
    reordered_modules = [module_ids[i] for i in col_order]
    reordered_names = [module_names[i] for i in col_order]

    return {
        'mag_ids': reordered_mags,
        'module_ids': reordered_modules,
        'module_names': reordered_names,
        'matrix': [[round(v, 3) for v in row] for row in reordered],
        'row_order': row_order,
        'col_order': col_order,
    }


def build_coverage(dastool_summary, contig2bin, depths_df):
    """Build coverage.json with per-MAG per-sample depth data."""
    print("Building coverage.json ...")

    depth_cols = [c for c in depths_df.columns if c.endswith('.sorted.bam') and not c.endswith('-var')]
    sample_names = [re.sub(r'\.sorted\.bam$', '', c) for c in depth_cols]

    # Build contig-to-MAG lookup
    c2b = {}
    if contig2bin is not None:
        for _, row in contig2bin.iterrows():
            c2b[row.iloc[0]] = row.iloc[1]

    mag_names = list(dastool_summary['bin'].values)
    # Per-MAG per-sample mean depth
    matrix = []
    for mag in mag_names:
        contigs = [c for c, b in c2b.items() if b == mag]
        contig_depths = depths_df[depths_df['contigName'].isin(contigs)]
        row = []
        for col in depth_cols:
            if len(contig_depths) > 0:
                row.append(round(float(contig_depths[col].mean()), 4))
            else:
                row.append(0)
        matrix.append(row)

    return {
        'mag_ids': mag_names,
        'sample_names': sample_names,
        'matrix': matrix,
    }


def build_mge_summary(virus_df, plasmid_df, defense_df, integron_df):
    """Build mge_summary.json with aggregate MGE stats."""
    print("Building mge_summary.json ...")

    # Defense system type counts
    defense_types = Counter()
    if defense_df is not None:
        for _, row in defense_df.iterrows():
            defense_types[row['type']] += 1

    # Viral taxonomy donut
    viral_taxonomy = Counter()
    if virus_df is not None:
        for _, row in virus_df.iterrows():
            tax = row.get('taxonomy', '')
            if pd.isna(tax) or not tax:
                viral_taxonomy['Unclassified'] += 1
            else:
                parts = [p.strip() for p in tax.split(';') if p.strip()]
                # Use the most specific available rank
                if len(parts) >= 4:
                    viral_taxonomy[parts[3]] += 1  # Order level
                elif len(parts) >= 2:
                    viral_taxonomy[parts[1]] += 1
                else:
                    viral_taxonomy[parts[0] if parts else 'Unclassified'] += 1

    return {
        'n_virus': len(virus_df) if virus_df is not None else 0,
        'n_plasmid': len(plasmid_df) if plasmid_df is not None else 0,
        'n_defense': len(defense_df) if defense_df is not None else 0,
        'n_integron_elements': len(integron_df) if integron_df is not None else 0,
        'defense_types': dict(defense_types),
        'viral_taxonomy': dict(viral_taxonomy),
    }


def build_mge_per_bin(contig2bin, virus_df, plasmid_df, defense_df, integron_df):
    """Build mge_per_bin.json with per-MAG MGE element lists."""
    print("Building mge_per_bin.json ...")

    c2b = {}
    if contig2bin is not None:
        for _, row in contig2bin.iterrows():
            c2b[row.iloc[0]] = row.iloc[1]

    per_bin = defaultdict(lambda: {'viruses': [], 'plasmids': [], 'defense': [], 'integrons': []})

    if virus_df is not None:
        for _, row in virus_df.iterrows():
            contig = row['seq_name']
            mag = c2b.get(contig, '_unbinned')
            per_bin[mag]['viruses'].append({
                'contig': contig,
                'length': int(row['length']),
                'score': round(float(row['virus_score']), 4),
                'taxonomy': str(row.get('taxonomy', '')),
                'n_hallmarks': int(row.get('n_hallmarks', 0)),
            })

    if plasmid_df is not None:
        for _, row in plasmid_df.iterrows():
            contig = row['seq_name']
            mag = c2b.get(contig, '_unbinned')
            per_bin[mag]['plasmids'].append({
                'contig': contig,
                'length': int(row['length']),
                'score': round(float(row['plasmid_score']), 4),
            })

    if defense_df is not None:
        for _, row in defense_df.iterrows():
            sid = row['sys_id']
            parts = sid.split('_')
            contig = '_'.join(parts[:2])
            mag = c2b.get(contig, '_unbinned')
            per_bin[mag]['defense'].append({
                'sys_id': sid,
                'type': row['type'],
                'subtype': row['subtype'],
                'contig': contig,
            })

    if integron_df is not None:
        for _, row in integron_df.iterrows():
            contig = row['ID_replicon']
            mag = c2b.get(contig, '_unbinned')
            if row.get('type_elt') == 'protein' or row.get('element', '').startswith('attc'):
                continue  # Skip individual elements, just track integrons
            per_bin[mag]['integrons'].append({
                'contig': contig,
                'integron_id': str(row.get('ID_integron', '')),
                'type': str(row.get('type', '')),
            })

    return dict(per_bin)


def build_eukaryotic(results_dir, assembly_info):
    """Build eukaryotic.json with Tiara/Whokaryote/MarFERReT summaries."""
    print("Building eukaryotic.json ...")

    len_map = dict(zip(assembly_info['#seq_name'], assembly_info['length']))

    # Tiara
    tiara_path = os.path.join(results_dir, 'eukaryotic', 'tiara', 'tiara_output.tsv')
    tiara_df = load_tsv(tiara_path)
    tiara_summary = Counter()
    tiara_size_weighted = Counter()
    if tiara_df is not None:
        for _, row in tiara_df.iterrows():
            cls = row['class_fst_stage']
            length = len_map.get(row['sequence_id'], 1000)
            tiara_summary[cls] += 1
            tiara_size_weighted[cls] += length

    # Whokaryote
    who_path = os.path.join(results_dir, 'eukaryotic', 'whokaryote', 'whokaryote_classifications.tsv')
    who_df = load_tsv(who_path)
    who_summary = Counter()
    if who_df is not None:
        for _, row in who_df.iterrows():
            cls = row.get('classification', row.get('predicted_class', 'unknown'))
            who_summary[cls] += 1

    # MarFERReT
    mar_path = os.path.join(results_dir, 'eukaryotic', 'marferret', 'marferret_contigs.tsv')
    mar_df = load_tsv(mar_path)
    mar_taxonomy = Counter()
    mar_contigs = []
    if mar_df is not None:
        for _, row in mar_df.iterrows():
            tax = row.get('top_taxonomy', '')
            if pd.notna(tax) and tax:
                # Parse phylum from taxonomy
                parts = [p.strip() for p in str(tax).split(';')]
                phylum = parts[1] if len(parts) > 1 else parts[0]
                if phylum:
                    mar_taxonomy[phylum] += 1
            mar_contigs.append({
                'contig': str(row['contig_id']),
                'n_proteins': int(row.get('n_proteins', 0)),
                'n_classified': int(row.get('n_classified', 0)),
                'taxonomy': str(tax) if pd.notna(tax) else '',
                'pfam': str(row.get('pfam_domains', '')) if pd.notna(row.get('pfam_domains', '')) else '',
            })

    return {
        'tiara_counts': dict(tiara_summary),
        'tiara_size_weighted': {k: int(v) for k, v in tiara_size_weighted.items()},
        'whokaryote_counts': dict(who_summary),
        'marferret_taxonomy': dict(mar_taxonomy),
        'marferret_contigs': mar_contigs[:500],  # Limit for performance
    }


def build_contig_explorer(results_dir, assembly_info, depths_df, contig2bin, kaiju_df,
                          skip_tsne=False):
    """Build contig_explorer.json with PCA + t-SNE + UMAP on fourth-root transformed TNF."""
    print("Building contig_explorer.json ...")

    from sklearn.decomposition import PCA
    from sklearn.preprocessing import StandardScaler

    # Load TNF (raw tetranucleotide frequencies — already relative abundance)
    tnf_path = os.path.join(results_dir, 'assembly', 'tnf.tsv')
    tnf_df = load_tsv(tnf_path, header=None)
    if tnf_df is None:
        return {'contigs': [], 'has_tsne': False, 'has_umap': False}

    contig_ids = tnf_df.iloc[:, 0].values
    tnf_matrix = tnf_df.iloc[:, 1:].values.astype(float)

    # Fourth-root transformation of relative abundances
    # Stabilizes variance for compositional data (Hellinger-like),
    # down-weights dominant k-mers, amplifies rare signal
    print("  Applying fourth-root transformation ...")
    tnf_transformed = np.power(np.abs(tnf_matrix), 0.25) * np.sign(tnf_matrix)

    # Build metadata lookups
    len_map = dict(zip(assembly_info['#seq_name'], assembly_info['length']))
    c2b = {}
    if contig2bin is not None:
        for _, row in contig2bin.iterrows():
            c2b[row.iloc[0]] = row.iloc[1]

    # Per-source taxonomy maps: each classifier gets its own domain/phylum per contig
    from collections import Counter

    def strip_gtdb_prefix(s):
        """Remove GTDB-style rank prefix like 'd__', 'p__', etc."""
        if s and len(s) > 3 and s[1:3] == '__':
            return s[3:]
        return s

    def parse_lineage(lineage):
        """Extract (domain, phylum) from semicolon-delimited lineage string."""
        parts = [p.strip() for p in lineage.split(';')]
        domain = parts[0] if len(parts) > 0 else ''
        phylum = parts[1] if len(parts) > 1 else ''
        return domain, phylum

    # Kaiju (six-frame or protein-level)
    kaiju_tax = {}  # contig_id -> {domain, phylum}
    if kaiju_df is not None:
        for _, row in kaiju_df.iterrows():
            lineage = row.get('lineage', '')
            if pd.notna(lineage) and lineage and lineage != 'Unclassified':
                domain, phylum = parse_lineage(lineage)
                if domain or phylum:
                    kaiju_tax[row['contig_id']] = {'domain': domain, 'phylum': phylum}
    print(f"  Kaiju taxonomy: {len(kaiju_tax)} contigs")

    # Kraken2 (GTDB k-mer)
    kraken2_tax = {}
    kraken2_path = os.path.join(results_dir, 'taxonomy', 'kraken2', 'kraken2_contigs.tsv')
    kraken2_df = load_tsv(kraken2_path)
    if kraken2_df is not None:
        for _, row in kraken2_df.iterrows():
            if row.get('status') != 'C':
                continue
            lineage = row.get('lineage', '')
            if pd.notna(lineage) and lineage and lineage != 'Unclassified':
                domain, phylum = parse_lineage(lineage)
                domain = strip_gtdb_prefix(domain)
                phylum = strip_gtdb_prefix(phylum)
                if domain or phylum:
                    kraken2_tax[row['contig_id']] = {'domain': domain, 'phylum': phylum}
    print(f"  Kraken2 taxonomy: {len(kraken2_tax)} contigs")

    # rRNA (SILVA SSU/LSU)
    rrna_tax = {}
    rrna_path = os.path.join(results_dir, 'taxonomy', 'rrna', 'rrna_contigs.tsv')
    rrna_df = load_tsv(rrna_path)
    if rrna_df is not None:
        for _, row in rrna_df.iterrows():
            lineage = row.get('best_ssu_taxonomy', '')
            if not lineage or pd.isna(lineage):
                lineage = row.get('best_lsu_taxonomy', '')
            if pd.notna(lineage) and lineage:
                domain, phylum = parse_lineage(lineage)
                if domain or phylum:
                    rrna_tax[row['contig_id']] = {'domain': domain, 'phylum': phylum}
    print(f"  rRNA taxonomy: {len(rrna_tax)} contigs")

    # Merged (Kraken2 > rRNA > Kaiju priority)
    tax_map = {}
    for cid, tax in kaiju_tax.items():
        tax_map[cid] = {**tax, 'source': 'kaiju'}
    for cid, tax in rrna_tax.items():
        tax_map[cid] = {**tax, 'source': 'rrna'}
    for cid, tax in kraken2_tax.items():
        tax_map[cid] = {**tax, 'source': 'kraken2'}
    sources = Counter(v['source'] for v in tax_map.values())
    print(f"  Merged taxonomy: {dict(sources)} (total {len(tax_map)} contigs)")

    depth_map = {}
    depth_cols = [c for c in depths_df.columns if c.endswith('.sorted.bam') and not c.endswith('-var')]
    if depth_cols:
        for _, row in depths_df.iterrows():
            depth_map[row['contigName']] = float(row['totalAvgDepth'])

    # StandardScaler + PCA on fourth-root transformed data
    print("  Running PCA on fourth-root transformed TNF ...")
    scaler = StandardScaler()
    tnf_scaled = scaler.fit_transform(tnf_transformed)

    n_components = min(50, tnf_scaled.shape[1], tnf_scaled.shape[0])
    pca = PCA(n_components=n_components)
    pca_coords = pca.fit_transform(tnf_scaled)
    pca_input = pca_coords[:, :min(50, n_components)]

    # t-SNE on PCA output of fourth-root transformed data
    tsne_coords = None
    if not skip_tsne:
        try:
            from sklearn.manifold import TSNE
            print("  Running t-SNE on fourth-root transformed TNF ...")
            tsne = TSNE(n_components=2, perplexity=30, random_state=42, max_iter=1000,
                        learning_rate='auto', init='pca')
            tsne_coords = tsne.fit_transform(pca_input)
            print("  t-SNE complete.")
        except Exception as e:
            print(f"  [WARNING] t-SNE failed: {e}", file=sys.stderr)

    # UMAP on PCA output of fourth-root transformed data
    umap_coords = None
    if not skip_tsne:  # UMAP gated by same flag (both are slow)
        try:
            import umap
            print("  Running UMAP on fourth-root transformed TNF ...")
            reducer = umap.UMAP(n_components=2, n_neighbors=15, min_dist=0.1,
                                metric='euclidean', random_state=42)
            umap_coords = reducer.fit_transform(pca_input)
            print("  UMAP complete.")
        except Exception as e:
            print(f"  [WARNING] UMAP failed: {e}", file=sys.stderr)

    # Build output records
    contigs = []
    for i, cid in enumerate(contig_ids):
        rec = {
            'id': cid,
            'length': int(len_map.get(cid, 0)),
            'bin': c2b.get(cid, ''),
            'depth': round(depth_map.get(cid, 0), 4),
            'pca_x': round(float(pca_coords[i, 0]), 4),
            'pca_y': round(float(pca_coords[i, 1]), 4),
        }
        # Merged taxonomy (best available)
        tax = tax_map.get(cid, {})
        rec['domain'] = tax.get('domain', '')
        rec['phylum'] = tax.get('phylum', '')
        rec['tax_source'] = tax.get('source', '')
        # Per-source taxonomy for individual classifier views
        kt = kaiju_tax.get(cid, {})
        rec['kaiju_domain'] = kt.get('domain', '')
        rec['kaiju_phylum'] = kt.get('phylum', '')
        k2t = kraken2_tax.get(cid, {})
        rec['kraken2_domain'] = k2t.get('domain', '')
        rec['kraken2_phylum'] = k2t.get('phylum', '')
        rt = rrna_tax.get(cid, {})
        rec['rrna_domain'] = rt.get('domain', '')
        rec['rrna_phylum'] = rt.get('phylum', '')

        if tsne_coords is not None:
            rec['tsne_x'] = round(float(tsne_coords[i, 0]), 4)
            rec['tsne_y'] = round(float(tsne_coords[i, 1]), 4)

        if umap_coords is not None:
            rec['umap_x'] = round(float(umap_coords[i, 0]), 4)
            rec['umap_y'] = round(float(umap_coords[i, 1]), 4)

        contigs.append(rec)

    return {
        'contigs': contigs,
        'has_tsne': tsne_coords is not None,
        'has_umap': umap_coords is not None,
        'pca_variance_explained': [round(float(v), 4) for v in pca.explained_variance_ratio_[:5]],
        'transform': 'fourth_root',
    }


def main():
    args = parse_args()
    results_dir = args.results
    output_dir = args.output
    os.makedirs(output_dir, exist_ok=True)

    print(f"Loading data from {results_dir} ...")

    # Load core data files
    assembly_info = load_tsv(os.path.join(results_dir, 'assembly', 'assembly_info.txt'),
                             comment='#', header=None,
                             names=['#seq_name', 'length', 'cov.', 'circ.', 'repeat',
                                    'mult.', 'alt_group', 'graph_path'])
    # Assembly info has header line starting with #, reload properly
    assembly_info = pd.read_csv(os.path.join(results_dir, 'assembly', 'assembly_info.txt'),
                                sep='\t')

    depths_df = load_tsv(os.path.join(results_dir, 'mapping', 'depths.txt'))
    contig2bin = load_tsv(os.path.join(results_dir, 'binning', 'dastool', 'contig2bin.tsv'),
                          header=None, names=['contig', 'bin'])
    dastool_summary = load_tsv(os.path.join(results_dir, 'binning', 'dastool', 'summary.tsv'))
    checkm2_df = load_tsv(os.path.join(results_dir, 'binning', 'checkm2', 'quality_report.tsv'))
    kaiju_df = load_tsv(os.path.join(results_dir, 'taxonomy', 'kaiju', 'kaiju_contigs.tsv'))
    virus_df = load_tsv(os.path.join(results_dir, 'mge', 'genomad', 'virus_summary.tsv'))
    plasmid_df = load_tsv(os.path.join(results_dir, 'mge', 'genomad', 'plasmid_summary.tsv'))
    defense_df = load_tsv(os.path.join(results_dir, 'mge', 'defensefinder', 'systems.tsv'))
    integron_df = load_tsv(os.path.join(results_dir, 'mge', 'integrons', 'integrons.tsv'),
                           comment='#')

    if assembly_info is None or depths_df is None or dastool_summary is None:
        print("[ERROR] Missing critical input files (assembly_info, depths, dastool summary)")
        sys.exit(1)

    # 1. overview.json
    overview = build_overview(results_dir, assembly_info, depths_df, dastool_summary,
                              checkm2_df, virus_df, plasmid_df, contig2bin)
    with open(os.path.join(output_dir, 'overview.json'), 'w') as f:
        json.dump(overview, f)
    print(f"  Wrote overview.json")

    # 2. contig_lengths.json
    contig_lengths = build_contig_lengths(assembly_info)
    with open(os.path.join(output_dir, 'contig_lengths.json'), 'w') as f:
        json.dump(contig_lengths, f)
    print(f"  Wrote contig_lengths.json")

    # 3. mags.json
    mags = build_mags(dastool_summary, checkm2_df, contig2bin, kaiju_df, depths_df,
                      virus_df, plasmid_df, defense_df, integron_df)
    with open(os.path.join(output_dir, 'mags.json'), 'w') as f:
        json.dump(mags, f)
    print(f"  Wrote mags.json ({len(mags)} MAGs)")

    # 4. checkm2_all.json
    checkm2_all = build_checkm2_all(checkm2_df, dastool_summary)
    with open(os.path.join(output_dir, 'checkm2_all.json'), 'w') as f:
        json.dump(checkm2_all, f)
    print(f"  Wrote checkm2_all.json ({len(checkm2_all)} bins)")

    # 5. taxonomy_sunburst.json
    sunburst = build_taxonomy_sunburst(kaiju_df, assembly_info)
    with open(os.path.join(output_dir, 'taxonomy_sunburst.json'), 'w') as f:
        json.dump(sunburst, f)
    print(f"  Wrote taxonomy_sunburst.json")

    # 6. kegg_heatmap.json
    kegg = build_kegg_heatmap(results_dir)
    with open(os.path.join(output_dir, 'kegg_heatmap.json'), 'w') as f:
        json.dump(kegg, f)
    print(f"  Wrote kegg_heatmap.json ({len(kegg['mag_ids'])} MAGs × {len(kegg['module_ids'])} modules)")

    # 7. coverage.json
    coverage = build_coverage(dastool_summary, contig2bin, depths_df)
    with open(os.path.join(output_dir, 'coverage.json'), 'w') as f:
        json.dump(coverage, f)
    print(f"  Wrote coverage.json")

    # 8. mge_summary.json
    mge_summary = build_mge_summary(virus_df, plasmid_df, defense_df, integron_df)
    with open(os.path.join(output_dir, 'mge_summary.json'), 'w') as f:
        json.dump(mge_summary, f)
    print(f"  Wrote mge_summary.json")

    # 9. mge_per_bin.json
    mge_per_bin = build_mge_per_bin(contig2bin, virus_df, plasmid_df, defense_df, integron_df)
    with open(os.path.join(output_dir, 'mge_per_bin.json'), 'w') as f:
        json.dump(mge_per_bin, f)
    print(f"  Wrote mge_per_bin.json ({len(mge_per_bin)} bins)")

    # 10. eukaryotic.json
    eukaryotic = build_eukaryotic(results_dir, assembly_info)
    with open(os.path.join(output_dir, 'eukaryotic.json'), 'w') as f:
        json.dump(eukaryotic, f)
    print(f"  Wrote eukaryotic.json")

    # 11. contig_explorer.json (largest, do last)
    contig_explorer = build_contig_explorer(results_dir, assembly_info, depths_df,
                                            contig2bin, kaiju_df,
                                            skip_tsne=args.skip_tsne)
    with open(os.path.join(output_dir, 'contig_explorer.json'), 'w') as f:
        json.dump(contig_explorer, f)
    print(f"  Wrote contig_explorer.json ({len(contig_explorer['contigs'])} contigs, "
          f"t-SNE={'yes' if contig_explorer['has_tsne'] else 'skipped'})")

    print("\nDone! All JSON files written to", output_dir)


if __name__ == '__main__':
    main()
