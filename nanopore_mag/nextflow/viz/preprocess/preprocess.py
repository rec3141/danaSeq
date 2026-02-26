#!/usr/bin/env python3
"""
Preprocess MAG pipeline TSV outputs into JSON for the web dashboard.

Generates 12 JSON files from Nextflow pipeline results:
  overview.json, mags.json, checkm2_all.json, taxonomy_sunburst.json,
  kegg_heatmap.json, coverage.json, mge_summary.json, mge_per_bin.json,
  eukaryotic.json, contig_explorer.json, contig_lengths.json

Usage:
  python3 preprocess.py --results <path> --output <path> [--store-dir <path>] [--skip-tsne] [--skip-umap]
"""

import argparse
import gzip
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


def write_json_gz(path, data, **dump_kwargs):
    """Write JSON to path and a gzip-compressed sidecar at path + '.gz'."""
    text = json.dumps(data, **dump_kwargs)
    with open(path, 'w') as f:
        f.write(text)
    with gzip.open(path + '.gz', 'wt', compresslevel=6) as f:
        f.write(text)


# ---------------------------------------------------------------------------
# Path resolution: storeDir (persistent cache) overlaid on results dir
# ---------------------------------------------------------------------------
_STORE_DIR = None   # set by main() when --store-dir is provided


def resolve_path(results_dir, *parts):
    """Return the best path for a results sub-path.

    If --store-dir was given and the file exists there, use it.
    Otherwise fall back to results_dir.  This lets preprocess.py
    transparently pick up outputs that haven't been published yet
    (e.g. BAKTA_EXTRA, metabolism sub-dirs).
    """
    rel = os.path.join(*parts)
    if _STORE_DIR:
        store_path = os.path.join(_STORE_DIR, rel)
        if os.path.exists(store_path):
            return store_path
    return os.path.join(results_dir, rel)


def load_annotation_tsv(results_dir):
    """Load best available annotation TSV: BAKTA_EXTRA > BAKTA_BASIC > PROKKA."""
    candidates = [
        resolve_path(results_dir, 'annotation', 'bakta', 'extra', 'annotation.tsv'),
        resolve_path(results_dir, 'annotation', 'bakta', 'basic', 'annotation.tsv'),
        resolve_path(results_dir, 'annotation', 'prokka', 'annotation.tsv'),
    ]
    for path in candidates:
        if os.path.exists(path):
            print(f"  Loading annotation from: {path}")
            return pd.read_csv(path, sep='\t', comment='#', header=None,
                               names=['contig', 'type', 'start', 'stop', 'strand',
                                      'locus_tag', 'gene', 'product', 'dbxrefs'])
    print("  [WARNING] No annotation TSV found (BAKTA_EXTRA/BAKTA_BASIC/PROKKA)", file=sys.stderr)
    return None


RANKS = ['domain', 'phylum', 'class', 'order', 'family', 'genus']


def strip_gtdb_prefix(s):
    """Remove GTDB-style rank prefix like 'd__', 'p__', etc."""
    if s and len(s) >= 3 and s[1:3] == '__':
        return s[3:]  # empty string for bare prefixes like 'f__'
    return s


def parse_lineage(lineage):
    """Extract taxonomy dict from semicolon-delimited lineage string."""
    parts = [p.strip() for p in lineage.split(';')]
    tax = {}
    for i, rank in enumerate(RANKS):
        tax[rank] = parts[i] if i < len(parts) and parts[i] else ''
    return tax


def build_taxonomy_maps(results_dir, kaiju_df, sendsketch_df=None):
    """Build per-contig taxonomy dicts for Kaiju, Kraken2, rRNA, and SendSketch."""
    kaiju_tax = {}
    if kaiju_df is not None:
        for _, row in kaiju_df.iterrows():
            lineage = row.get('lineage', '')
            if pd.notna(lineage) and lineage and lineage != 'Unclassified':
                tax = parse_lineage(lineage)
                if any(tax.values()):
                    kaiju_tax[row['contig_id']] = tax

    kraken2_tax = {}
    kraken2_path = resolve_path(results_dir, 'taxonomy', 'kraken2', 'kraken2_contigs.tsv')
    kraken2_df = load_tsv(kraken2_path)
    if kraken2_df is not None:
        for _, row in kraken2_df.iterrows():
            if row.get('status') != 'C':
                continue
            lineage = row.get('lineage', '')
            if pd.notna(lineage) and lineage and lineage != 'Unclassified':
                tax = parse_lineage(lineage)
                tax = {k: strip_gtdb_prefix(v) for k, v in tax.items()}
                if any(tax.values()):
                    kraken2_tax[row['contig_id']] = tax

    rrna_tax = {}
    rrna_path = resolve_path(results_dir, 'taxonomy', 'rrna', 'rrna_contigs.tsv')
    rrna_df = load_tsv(rrna_path)
    if rrna_df is not None:
        for _, row in rrna_df.iterrows():
            lineage = row.get('best_ssu_taxonomy', '')
            if not lineage or pd.isna(lineage):
                lineage = row.get('best_lsu_taxonomy', '')
            if pd.notna(lineage) and lineage:
                tax = parse_lineage(lineage)
                if any(tax.values()):
                    rrna_tax[row['contig_id']] = tax

    sendsketch_tax = {}
    if sendsketch_df is not None:
        for _, row in sendsketch_df.iterrows():
            if row.get('status') != 'C':
                continue
            lineage = row.get('lineage', '')
            if pd.notna(lineage) and lineage and lineage != 'Unclassified':
                tax = parse_lineage(lineage)
                tax = {k: strip_gtdb_prefix(v) for k, v in tax.items()}
                if any(tax.values()):
                    sendsketch_tax[row['contig_id']] = tax

    return {'kaiju': kaiju_tax, 'kraken2': kraken2_tax, 'rrna': rrna_tax, 'sendsketch': sendsketch_tax}


def majority_vote_taxonomy(contigs, tax_map, len_map=None):
    """Compute majority-vote taxonomy from a contig->tax dict.

    Votes on the full lineage (length-weighted if len_map provided) to avoid
    chimeric assignments where each rank comes from a different organism.
    """
    entries = [(c, tax_map[c]) for c in contigs if c in tax_map]
    if not entries:
        return {}
    # Build a lineage string per contig and accumulate weight
    lineage_weights = Counter()
    for cid, lin in entries:
        weight = len_map.get(cid, 1) if len_map else 1
        # Use the deepest non-empty rank as the lineage key
        key = tuple(lin.get(r, '') for r in RANKS)
        lineage_weights[key] += weight
    # Winner is the full lineage with the most weight
    winner = lineage_weights.most_common(1)[0][0]
    tax = {}
    for i, rank in enumerate(RANKS):
        if winner[i]:
            tax[rank] = winner[i]
    return tax


def parse_args():
    p = argparse.ArgumentParser(description="Preprocess MAG pipeline results to JSON")
    p.add_argument('--results', '-r', required=True, help='Pipeline results directory')
    p.add_argument('--output', '-o', required=True, help='Output directory for JSON files')
    p.add_argument('--store-dir', '-s', default=None,
                   help='Persistent storeDir overlay (checked first for each file)')
    p.add_argument('--skip-tsne', action='store_true', help='Skip t-SNE (saves ~5 min)')
    p.add_argument('--skip-umap', action='store_true', help='Skip UMAP (saves ~5-15 min)')
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
    if dastool_summary is None:
        n_mags = hq = mq = lq = 0
        binner_counts = Counter()
        overview = {
            'assembly_size': total_size,
            'n50': n50,
            'n_contigs': n_contigs,
            'n_mags': 0,
            'hq': 0,
            'mq': 0,
            'lq': 0,
            'n_virus': len(virus_df) if virus_df is not None else 0,
            'n_plasmid': len(plasmid_df) if plasmid_df is not None else 0,
            'binner_counts': {},
        }
        return overview

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


def build_contig_lengths(assembly_info, depths_df=None, results_dir=None):
    """Build contig_lengths.json with pre-bucketed histogram data, coverage histogram, and scatter."""
    print("Building contig_lengths.json ...")
    lengths = assembly_info['length'].values

    # Log-spaced bins from min to max
    min_len = max(int(lengths.min()), 1)
    max_len = int(lengths.max())
    bins = np.logspace(np.log10(min_len), np.log10(max_len + 1), num=50)
    counts, edges = np.histogram(lengths, bins=bins)

    result = {
        'bin_edges': [int(e) for e in edges],
        'counts': [int(c) for c in counts],
        'n50': int(np.sort(lengths)[::-1][np.searchsorted(
            np.cumsum(np.sort(lengths)[::-1]), lengths.sum() / 2)]),
    }

    # Coverage histogram and length-coverage scatter (if depths available)
    if depths_df is not None:
        depth_map = {}
        for _, row in depths_df.iterrows():
            depth_map[row['contigName']] = float(row['totalAvgDepth'])

        contig_names = assembly_info['#seq_name'].values
        depths = np.array([depth_map.get(c, 0.0) for c in contig_names])
        positive = depths > 0

        # Coverage histogram: log10-space bins, equal-width bars
        if positive.any():
            log_depths = np.log10(depths[positive])
            cov_min = float(np.floor(log_depths.min()))
            cov_max = float(np.ceil(log_depths.max()))
            cov_bins = np.linspace(cov_min, cov_max, num=50)
            cov_counts, cov_edges = np.histogram(log_depths, bins=cov_bins)
            result['cov_log_edges'] = [round(float(e), 4) for e in cov_edges]
            result['cov_counts'] = [int(c) for c in cov_counts]
            print(f"  Coverage histogram: {int(positive.sum())} contigs with depth > 0")

        # Length-coverage scatter: all contigs assigned to any bin
        if positive.any():
            binned_contigs = set()
            if results_dir:
                for binner in ['semibin', 'metabat', 'maxbin', 'lorbin', 'comebin']:
                    btsv = resolve_path(results_dir, 'binning', binner, f'{binner}_bins.tsv')
                    bdf = load_tsv(btsv, header=None, names=['contig', 'bin'])
                    if bdf is not None:
                        binned_contigs.update(bdf['contig'].values)
                c2b_path = resolve_path(results_dir, 'binning', 'dastool', 'contig2bin.tsv')
                c2b_df = load_tsv(c2b_path, header=None, names=['contig', 'bin'])
                if c2b_df is not None:
                    binned_contigs.update(c2b_df['contig'].values)

            pos_idx = np.where(positive)[0]
            if binned_contigs:
                pos_idx = np.array([i for i in pos_idx if contig_names[i] in binned_contigs])
            pos_idx.sort()
            result['scatter_log_length'] = [round(float(np.log10(lengths[i])), 2) for i in pos_idx]
            result['scatter_log_depth'] = [round(float(np.log10(depths[i])), 2) for i in pos_idx]
            # GC% for coloring
            if results_dir:
                gc_path = resolve_path(results_dir, 'assembly', 'gc.tsv')
                gc_df = load_tsv(gc_path)
                if gc_df is not None:
                    gc_map = dict(zip(gc_df['contig_id'], gc_df['gc_pct'].astype(float)))
                    result['scatter_gc'] = [round(float(gc_map.get(contig_names[i], 0)), 1) for i in pos_idx]
            print(f"  Length-coverage scatter: {len(pos_idx)} points")

    return result


def build_mags(dastool_summary, checkm2_df, contig2bin, kaiju_df, depths_df,
               virus_df, plasmid_df, defense_df, integron_df, annotation_df=None):
    """Build mags.json with per-MAG joined records."""
    print("Building mags.json ...")

    if dastool_summary is None:
        return []

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
    if depths_df is not None:
        depth_cols = [c for c in depths_df.columns if c.endswith('.sorted.bam') and not c.endswith('-var')]
    else:
        depth_cols = []
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

    # Build per-MAG taxonomy by majority vote from kaiju (full lineage vote)
    mag_taxonomy = {}
    if kaiju_df is not None:
        kaiju_tax = {}
        for _, row in kaiju_df.iterrows():
            lineage = row.get('lineage', '')
            if pd.notna(lineage) and lineage and lineage != 'Unclassified':
                kaiju_tax[row['contig_id']] = parse_lineage(lineage)
        for mag_name, contigs in mag_contigs.items():
            tax = majority_vote_taxonomy(contigs, kaiju_tax)
            if tax:
                mag_taxonomy[mag_name] = tax

    # Per-MAG gene counts from annotation TSV
    mag_gene_stats = {}
    if annotation_df is not None:
        contig_genes = defaultdict(lambda: {'n_cds': 0, 'n_hypo': 0, 'n_trna': 0, 'n_ncrna': 0})
        for _, row in annotation_df.iterrows():
            contig = row['contig']
            gene_type = str(row['type']).lower()
            product = str(row.get('product', '')).lower()
            if gene_type == 'cds':
                contig_genes[contig]['n_cds'] += 1
                if 'hypothetical' in product:
                    contig_genes[contig]['n_hypo'] += 1
            elif gene_type == 'trna':
                contig_genes[contig]['n_trna'] += 1
            elif gene_type in ('ncrna', 'tmrna', 'crispr'):
                contig_genes[contig]['n_ncrna'] += 1
        for mag_name_key, contigs_list in mag_contigs.items():
            n_cds = sum(contig_genes[c]['n_cds'] for c in contigs_list)
            n_hypo = sum(contig_genes[c]['n_hypo'] for c in contigs_list)
            n_trna = sum(contig_genes[c]['n_trna'] for c in contigs_list)
            n_ncrna = sum(contig_genes[c]['n_ncrna'] for c in contigs_list)
            mag_gene_stats[mag_name_key] = {
                'n_genes': n_cds,
                'hypo_pct': round(n_hypo / n_cds * 100, 1) if n_cds > 0 else 0,
                'n_trna': n_trna,
                'n_ncrna': n_ncrna,
            }

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

        gene_stats = mag_gene_stats.get(mag_name, {'n_genes': 0, 'hypo_pct': 0, 'n_trna': 0, 'n_ncrna': 0})
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
            'n_genes': gene_stats['n_genes'],
            'hypo_pct': gene_stats['hypo_pct'],
            'n_trna': gene_stats['n_trna'],
            'n_ncrna': gene_stats['n_ncrna'],
        })

    return mags


def build_checkm2_all(results_dir, checkm2_df, dastool_summary, contig2bin,
                      kaiju_df, virus_df, plasmid_df, defense_df, integron_df,
                      assembly_info=None, sendsketch_df=None):
    """Build checkm2_all.json for all bins with quality, taxonomy, and MGE."""
    print("Building checkm2_all.json ...")
    if checkm2_df is None:
        return []

    dastool_names = set(dastool_summary['bin'].values)

    # Separate contig maps for DAS Tool vs raw binner bins
    binners = ['semibin', 'metabat', 'maxbin', 'lorbin', 'comebin']
    dastool_bin_contigs = {}  # dastool full name -> [contig_ids]
    binner_bin_contigs = {}   # raw bin name -> [contig_ids]

    # DAS Tool contig2bin (keyed by full dastool- name)
    if contig2bin is not None:
        for _, row in contig2bin.iterrows():
            mag = row.iloc[1]  # e.g., dastool-semibin_022
            dastool_bin_contigs.setdefault(mag, []).append(row.iloc[0])

    # Per-binner contig2bin (keyed by raw bin name)
    for binner in binners:
        bins_path = resolve_path(results_dir, 'binning', binner, f'{binner}_bins.tsv')
        bdf = load_tsv(bins_path, header=None, names=['contig', 'bin'])
        if bdf is not None:
            for _, row in bdf.iterrows():
                binner_bin_contigs.setdefault(row['bin'], []).append(row['contig'])

    # MGE contig sets
    virus_contigs = set(virus_df['seq_name'].values) if virus_df is not None else set()
    plasmid_contigs = set(plasmid_df['seq_name'].values) if plasmid_df is not None else set()

    defense_contig_map = defaultdict(list)
    if defense_df is not None:
        for _, row in defense_df.iterrows():
            sid = row['sys_id']
            parts = sid.split('_')
            contig = '_'.join(parts[:2])
            defense_contig_map[contig].append(row['type'])

    integron_contigs = set()
    if integron_df is not None:
        for _, row in integron_df.iterrows():
            integron_contigs.add(row['ID_replicon'])

    # Per-source taxonomy maps
    tax_maps = build_taxonomy_maps(results_dir, kaiju_df, sendsketch_df)

    # Contig length map for length-weighted composition
    len_map = {}
    if assembly_info is not None:
        len_map = dict(zip(assembly_info['#seq_name'], assembly_info['length']))

    records = []
    for _, row in checkm2_df.iterrows():
        name = row['Name']
        # Only dastool- prefixed CheckM2 entries are DAS Tool consensus bins
        if name.startswith('dastool-'):
            is_dastool = name in dastool_names
            short_name = name.replace('dastool-', '')
            contigs = dastool_bin_contigs.get(name, [])
        else:
            is_dastool = False
            short_name = name
            contigs = binner_bin_contigs.get(name, [])

        m = re.match(r'(\w+?)_\d+', short_name)
        binner = m.group(1) if m else 'unknown'

        # Per-source majority vote taxonomy (length-weighted full lineage)
        taxonomy = {}
        for source, tmap in tax_maps.items():
            tax = majority_vote_taxonomy(contigs, tmap, len_map)
            if tax:
                taxonomy[source] = tax

        # Per-source per-rank length-weighted composition
        composition = {}
        for source, tmap in tax_maps.items():
            src_comp = {}
            for cid in contigs:
                clen = len_map.get(cid, 0)
                if clen == 0:
                    continue
                ctax = tmap.get(cid)
                for rk in RANKS:
                    taxon = (ctax.get(rk, '') if ctax else '') or 'Unclassified'
                    src_comp.setdefault(rk, {})
                    src_comp[rk][taxon] = src_comp[rk].get(taxon, 0) + clen
            if src_comp:
                composition[source] = src_comp

        n_virus = sum(1 for c in contigs if c in virus_contigs)
        n_plasmid = sum(1 for c in contigs if c in plasmid_contigs)
        n_defense = sum(len(defense_contig_map.get(c, [])) for c in contigs)
        n_integron = sum(1 for c in contigs if c in integron_contigs)

        records.append({
            'name': name,
            'completeness': round(float(row['Completeness']), 2),
            'contamination': round(float(row['Contamination']), 2),
            'genome_size': int(row['Genome_Size']),
            'gc': round(float(row['GC_Content']), 4),
            'n50': int(row['Contig_N50']),
            'binner': binner,
            'is_dastool': is_dastool,
            'dastool_name': name if is_dastool else None,
            'taxonomy': taxonomy,
            'composition': composition,
            'n_contigs': len(contigs),
            'n_virus': n_virus,
            'n_plasmid': n_plasmid,
            'n_defense': n_defense,
            'n_integron': n_integron,
        })
    return records


def _build_sunburst_tree(tax_map, len_map):
    """Build a D3 hierarchy tree from a contig->taxonomy dict."""
    tree = {}
    ranks = ['domain', 'phylum', 'class', 'order']

    for contig_id, tax in tax_map.items():
        contig_len = len_map.get(contig_id, 1000)
        node = tree
        for rank in ranks:
            taxon = tax.get(rank, '') or 'Unclassified'
            if taxon not in node:
                node[taxon] = {'_size': 0, '_children': {}}
            node[taxon]['_size'] += contig_len
            node = node[taxon]['_children']

    if not tree:
        return {'name': 'Life', 'children': []}

    total_size = sum(v['_size'] for v in tree.values())
    threshold = total_size * 0.005

    def to_hierarchy(node_dict):
        children = []
        other_size = 0
        for name, data in sorted(node_dict.items(), key=lambda x: -x[1]['_size']):
            if data['_size'] < threshold:
                other_size += data['_size']
            else:
                child = {'name': name, 'value': data['_size']}
                if data['_children']:
                    sub = to_hierarchy(data['_children'])
                    if sub:
                        child['children'] = sub
                children.append(child)
        if other_size > 0:
            children.append({'name': 'Other', 'value': other_size})
        return children

    return {'name': 'Life', 'children': to_hierarchy(tree)}


def build_taxonomy_sunburst(results_dir, kaiju_df, assembly_info, contig2bin, sendsketch_df=None):
    """Build taxonomy_sunburst.json with per-source sunburst trees + composition."""
    print("Building taxonomy_sunburst.json ...")
    len_map = dict(zip(assembly_info['#seq_name'], assembly_info['length']))
    tax_maps = build_taxonomy_maps(results_dir, kaiju_df, sendsketch_df)

    result = {'sunbursts': {}, 'composition': {}}
    for source, tmap in tax_maps.items():
        result['sunbursts'][source] = _build_sunburst_tree(tmap, len_map)
        print(f"  {source}: {len(tmap)} contigs")

    # Per-binner length-weighted composition at each rank
    # Load per-binner contig sets
    binners = ['semibin', 'metabat', 'maxbin', 'lorbin', 'comebin']
    binner_contigs = {}  # binner -> set of contig ids

    # DAS Tool
    if contig2bin is not None:
        dt_contigs = set()
        for _, row in contig2bin.iterrows():
            dt_contigs.add(row.iloc[0])
        binner_contigs['dastool'] = dt_contigs

    for binner in binners:
        bins_path = resolve_path(results_dir, 'binning', binner, f'{binner}_bins.tsv')
        bdf = load_tsv(bins_path, header=None, names=['contig', 'bin'])
        if bdf is not None:
            binner_contigs[binner] = set(bdf['contig'].values)

    # Compute: composition[binner][source][rank] = {taxon: total_bp, ...}
    comp = {}
    for binner, contigs in binner_contigs.items():
        comp[binner] = {}
        for source, tmap in tax_maps.items():
            rank_data = {}
            for rank in RANKS:
                counts = {}
                for cid in contigs:
                    tax = tmap.get(cid)
                    if tax:
                        val = tax.get(rank, '') or 'Unclassified'
                    else:
                        val = 'Unclassified'
                    length = len_map.get(cid, 0)
                    counts[val] = counts.get(val, 0) + length
                rank_data[rank] = counts
            comp[binner][source] = rank_data

    result['composition'] = comp
    return result


# Merge fine-grained KEGG categories into 5 display groups
_CATEGORY_GROUPS = {
    'Central carbon':     'Carbon',
    'Carbon fixation':    'Carbon',
    'Methane metabolism':  'Carbon',
    'Fermentation':       'Carbon',
    'Nitrogen metabolism': 'Nitrogen',
    'Sulfur metabolism':   'Sulfur',
    'Energy metabolism':   'Energy',
    'Photosynthesis':     'Energy',
    'Hydrogen metabolism': 'Energy',
    'Vitamins':           'Biosynthesis',
    'Amino acids':        'Biosynthesis',
    'Xenobiotics':        'Other',
    'Transporters':       'Other',
}


def _load_kegg_categories(results_dir):
    """Load module->category and module->group mappings from kegg_module_completeness.py."""
    import importlib.util
    script = os.path.join(results_dir, '..', 'bin', 'kegg_module_completeness.py')
    if not os.path.exists(script):
        # Try relative to this file
        script = os.path.join(os.path.dirname(__file__), '..', '..', 'bin', 'kegg_module_completeness.py')
    if os.path.exists(script):
        spec = importlib.util.spec_from_file_location('kegg_mod', script)
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        cats = {mid: cat for mid, (name, cat, steps) in mod.KEGG_MODULES.items()}
        groups = {mid: _CATEGORY_GROUPS.get(cat, 'Other') for mid, cat in cats.items()}
        return cats, groups
    return {}, {}


def build_kegg_heatmap(results_dir):
    """Build kegg_heatmap.json with Ward-clustered module completeness matrix."""
    print("Building kegg_heatmap.json ...")
    path = resolve_path(results_dir, 'metabolism', 'modules', 'module_completeness.tsv')
    df = load_tsv(path)
    if df is None:
        return {'mag_ids': [], 'module_ids': [], 'module_names': [], 'matrix': [],
                'row_order': [], 'col_order': [], 'module_categories': []}

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

    # Load module categories from pipeline script
    cat_map, group_map = _load_kegg_categories(results_dir)
    module_categories = [cat_map.get(mid, 'Other') for mid in module_ids]
    module_groups = [group_map.get(mid, 'Other') for mid in module_ids]

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

    # Reorder data by clustering; send identity order so D3Heatmap doesn't double-reorder
    reordered = matrix[row_order][:, col_order].tolist()
    reordered_mags = [df.index[i] for i in row_order]
    reordered_modules = [module_ids[i] for i in col_order]
    reordered_names = [module_names[i] for i in col_order]
    reordered_categories = [module_categories[i] for i in col_order]
    reordered_groups = [module_groups[i] for i in col_order]

    n_rows = len(reordered_mags)
    n_cols = len(reordered_modules)

    return {
        'mag_ids': reordered_mags,
        'module_ids': reordered_modules,
        'module_names': reordered_names,
        'module_categories': reordered_categories,
        'module_groups': reordered_groups,
        'matrix': [[round(v, 3) for v in row] for row in reordered],
        'row_order': list(range(n_rows)),
        'col_order': list(range(n_cols)),
    }


def build_scg_heatmap(results_dir):
    """Build scg_heatmap.json with Ward-clustered marker gene presence/absence matrix."""
    print("Building scg_heatmap.json ...")

    # Load contig2bin maps for all binners
    binners = ['semibin', 'metabat', 'maxbin', 'lorbin', 'comebin']
    contig_to_bins = defaultdict(set)  # contig -> set of bin names

    # DAS Tool consensus (prefixed names like dastool-semibin_001)
    c2b_path = resolve_path(results_dir, 'binning', 'dastool', 'contig2bin.tsv')
    c2b_df = load_tsv(c2b_path, header=None, names=['contig', 'bin'])
    if c2b_df is not None:
        for _, row in c2b_df.iterrows():
            contig_to_bins[row['contig']].add(row['bin'])

    # Per-binner maps
    for binner in binners:
        bins_path = resolve_path(results_dir, 'binning', binner, f'{binner}_bins.tsv')
        bdf = load_tsv(bins_path, header=None, names=['contig', 'bin'])
        if bdf is not None:
            for _, row in bdf.iterrows():
                contig_to_bins[row['contig']].add(row['bin'])

    result = {}
    for scg_set in ['bacteria', 'archaea']:
        scg_path = resolve_path(results_dir, 'binning', 'dastool', f'{scg_set}.scg')
        scg_df = load_tsv(scg_path, header=None, names=['protein_id', 'scg_name'])
        if scg_df is None or scg_df.empty:
            result[scg_set] = {
                'mag_ids': [], 'module_ids': [], 'module_names': [],
                'matrix': [], 'row_order': [], 'col_order': [],
            }
            continue

        # Build bin x marker count matrix
        bin_marker_counts = defaultdict(Counter)  # bin -> {marker: count}
        all_markers = set()

        for _, row in scg_df.iterrows():
            protein_id = str(row['protein_id'])
            marker = str(row['scg_name'])
            all_markers.add(marker)

            # Extract contig from protein ID: "contig_000035_269" -> "contig_000035"
            parts = protein_id.rsplit('_', 1)
            contig = parts[0] if len(parts) > 1 else protein_id

            # Assign to all bins this contig belongs to
            for bin_name in contig_to_bins.get(contig, []):
                bin_marker_counts[bin_name][marker] += 1

        if not bin_marker_counts:
            result[scg_set] = {
                'mag_ids': [], 'module_ids': [], 'module_names': [],
                'matrix': [], 'row_order': [], 'col_order': [],
            }
            continue

        # Sort for deterministic order
        mag_ids = sorted(bin_marker_counts.keys())
        marker_ids = sorted(all_markers)

        # Build count matrix
        matrix = []
        for mag in mag_ids:
            row = [bin_marker_counts[mag].get(m, 0) for m in marker_ids]
            matrix.append(row)
        mat = np.array(matrix, dtype=float)

        # Ward hierarchical clustering on rows (bins) and columns (markers)
        if mat.shape[0] > 2:
            row_dist = pdist(mat, metric='euclidean')
            row_linkage = linkage(row_dist, method='ward')
            row_order = leaves_list(row_linkage).tolist()
        else:
            row_order = list(range(mat.shape[0]))

        if mat.shape[1] > 2:
            col_dist = pdist(mat.T, metric='euclidean')
            col_linkage = linkage(col_dist, method='ward')
            col_order = leaves_list(col_linkage).tolist()
        else:
            col_order = list(range(mat.shape[1]))

        # Reorder by clustering
        reordered = mat[row_order][:, col_order].tolist()
        reordered_mags = [mag_ids[i] for i in row_order]
        reordered_markers = [marker_ids[i] for i in col_order]

        n_rows = len(reordered_mags)
        n_cols = len(reordered_markers)

        result[scg_set] = {
            'mag_ids': reordered_mags,
            'module_ids': reordered_markers,
            'module_names': reordered_markers,
            'matrix': [[int(v) for v in row] for row in reordered],
            'row_order': list(range(n_rows)),
            'col_order': list(range(n_cols)),
        }
        print(f"  {scg_set}: {n_rows} bins × {n_cols} markers")

    return result


def build_coverage(results_dir, dastool_summary, contig2bin, depths_df):
    """Build coverage.json with per-bin per-sample depth for all binners, Bray-Curtis clustered."""
    print("Building coverage.json ...")

    if depths_df is None or dastool_summary is None:
        return {'bins': [], 'sample_names': [], 'matrix': []}

    depth_cols = [c for c in depths_df.columns if c.endswith('.sorted.bam') and not c.endswith('-var')]
    sample_names = [re.sub(r'\.sorted\.bam$', '', c) for c in depth_cols]

    def compute_bin_depths(c2b_map, bin_names):
        """Compute per-bin per-sample mean depth matrix."""
        rows = []
        for b in bin_names:
            contigs = [c for c, bn in c2b_map.items() if bn == b]
            cd = depths_df[depths_df['contigName'].isin(contigs)]
            row = []
            for col in depth_cols:
                if len(cd) > 0:
                    row.append(round(float(cd[col].mean()), 4))
                else:
                    row.append(0)
            rows.append(row)
        return rows

    # DAS Tool consensus bins
    dastool_c2b = {}
    if contig2bin is not None:
        for _, row in contig2bin.iterrows():
            dastool_c2b[row.iloc[0]] = row.iloc[1]
    dastool_bins = list(dastool_summary['bin'].values)

    # All binners
    binners = ['semibin', 'metabat', 'maxbin', 'lorbin', 'comebin']
    all_bins = []   # list of {id, binner, depths: []}
    all_matrix = []

    # Add DAS Tool bins
    dastool_depths = compute_bin_depths(dastool_c2b, dastool_bins)
    for name, depths in zip(dastool_bins, dastool_depths):
        all_bins.append({'id': name, 'binner': 'dastool'})
        all_matrix.append(depths)

    # Add per-binner bins
    for binner in binners:
        bins_path = resolve_path(results_dir, 'binning', binner, f'{binner}_bins.tsv')
        bdf = load_tsv(bins_path, header=None, names=['contig', 'bin'])
        if bdf is None:
            continue
        c2b_map = dict(zip(bdf['contig'], bdf['bin']))
        bin_names = sorted(set(c2b_map.values()))
        depths = compute_bin_depths(c2b_map, bin_names)
        for name, d in zip(bin_names, depths):
            all_bins.append({'id': name, 'binner': binner})
            all_matrix.append(d)
        print(f"  {binner}: {len(bin_names)} bins")

    mat = np.array(all_matrix)
    print(f"  Total bins: {len(all_bins)} x {len(sample_names)} samples")

    # Fourth-root transform of proportions for Bray-Curtis clustering
    row_sums = mat.sum(axis=1, keepdims=True)
    row_sums[row_sums == 0] = 1
    proportions = mat / row_sums
    transformed = np.power(proportions, 0.25)

    # Bray-Curtis hierarchical clustering on rows (all bins)
    row_order = list(range(len(all_bins)))
    if transformed.shape[0] > 2:
        row_dist = pdist(transformed, metric='braycurtis')
        row_dist = np.nan_to_num(row_dist, nan=0.0)
        row_link = linkage(row_dist, method='average')
        row_order = leaves_list(row_link).tolist()

    # Bray-Curtis clustering on columns (samples)
    col_order = list(range(len(sample_names)))
    if transformed.shape[1] > 2:
        col_dist = pdist(transformed.T, metric='braycurtis')
        col_dist = np.nan_to_num(col_dist, nan=0.0)
        col_link = linkage(col_dist, method='average')
        col_order = leaves_list(col_link).tolist()

    # Reorder
    reordered_bins = [all_bins[i] for i in row_order]
    reordered_samples = [sample_names[i] for i in col_order]
    reordered_matrix = mat[row_order][:, col_order].tolist()
    reordered_matrix = [[round(v, 4) for v in row] for row in reordered_matrix]

    return {
        'bins': reordered_bins,         # [{id, binner}, ...]
        'sample_names': reordered_samples,
        'matrix': reordered_matrix,
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
    tiara_path = resolve_path(results_dir, 'eukaryotic', 'tiara', 'tiara_output.tsv')
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
    who_path = resolve_path(results_dir, 'eukaryotic', 'whokaryote', 'whokaryote_classifications.tsv')
    who_df = load_tsv(who_path)
    who_summary = Counter()
    if who_df is not None:
        for _, row in who_df.iterrows():
            cls = row.get('classification', row.get('predicted_class', 'unknown'))
            who_summary[cls] += 1

    # MarFERReT
    mar_path = resolve_path(results_dir, 'eukaryotic', 'marferret', 'marferret_contigs.tsv')
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
                          skip_tsne=False, skip_umap=False, output_dir=None, sendsketch_df=None):
    """Build contig_explorer.json with PCA + t-SNE + UMAP on fourth-root transformed TNF."""
    print("Building contig_explorer.json ...")

    from sklearn.decomposition import PCA
    from sklearn.preprocessing import StandardScaler

    # Load TNF (raw tetranucleotide frequencies — already relative abundance)
    tnf_path = resolve_path(results_dir, 'assembly', 'tnf.tsv')
    tnf_df = load_tsv(tnf_path, header=None)
    if tnf_df is None:
        return {'contigs': [], 'has_tsne': False, 'has_umap': False}

    contig_ids = tnf_df.iloc[:, 0].values
    tnf_matrix = tnf_df.iloc[:, 1:].values.astype(float)

    # Load GC% from pipeline output (computed alongside TNF)
    gc_path = resolve_path(results_dir, 'assembly', 'gc.tsv')
    gc_df = load_tsv(gc_path)
    if gc_df is not None:
        gc_map = dict(zip(gc_df['contig_id'], gc_df['gc_pct'].astype(float)))
        gc_vals = np.array(list(gc_map.values()))
        print(f"  GC%: median={np.median(gc_vals):.1f}%, range=[{gc_vals.min():.1f}%, {gc_vals.max():.1f}%]")
    else:
        gc_map = {}
        print("  [WARNING] gc.tsv not found, GC% will be unavailable")

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

    # Load per-binner contig2bin maps
    binners = ['semibin', 'metabat', 'maxbin', 'lorbin', 'comebin']
    binner_maps = {}
    for binner in binners:
        bins_path = resolve_path(results_dir, 'binning', binner, f'{binner}_bins.tsv')
        bdf = load_tsv(bins_path, header=None, names=['contig', 'bin'])
        if bdf is not None:
            bmap = dict(zip(bdf['contig'], bdf['bin']))
            binner_maps[binner] = bmap
            print(f"  {binner}: {len(set(bmap.values()))} bins, {len(bmap)} contigs")

    # Per-source taxonomy maps (shared helper)
    tax_sources = build_taxonomy_maps(results_dir, kaiju_df, sendsketch_df)
    kaiju_tax = tax_sources['kaiju']
    kraken2_tax = tax_sources['kraken2']
    rrna_tax = tax_sources['rrna']
    sendsketch_tax = tax_sources['sendsketch']
    print(f"  Kaiju taxonomy: {len(kaiju_tax)} contigs")
    print(f"  Kraken2 taxonomy: {len(kraken2_tax)} contigs")
    print(f"  rRNA taxonomy: {len(rrna_tax)} contigs")
    print(f"  SendSketch taxonomy: {len(sendsketch_tax)} contigs")

    # Merged (Kraken2 > SendSketch > rRNA > Kaiju priority)
    tax_map = {}
    for cid, tax in kaiju_tax.items():
        tax_map[cid] = {**tax, 'source': 'kaiju'}
    for cid, tax in rrna_tax.items():
        tax_map[cid] = {**tax, 'source': 'rrna'}
    for cid, tax in sendsketch_tax.items():
        tax_map[cid] = {**tax, 'source': 'sendsketch'}
    for cid, tax in kraken2_tax.items():
        tax_map[cid] = {**tax, 'source': 'kraken2'}
    sources = Counter(v['source'] for v in tax_map.values())
    print(f"  Merged taxonomy: {dict(sources)} (total {len(tax_map)} contigs)")

    # geNomad replicon classification: only use FDR-filtered calls from summary files
    # Default all contigs to 'chromosome', then overlay confident virus/plasmid/provirus
    replicon_map = {cid: 'chromosome' for cid in assembly_info['#seq_name'].values}

    virus_summary_path = resolve_path(results_dir, 'mge', 'genomad', 'virus_summary.tsv')
    virus_sum_df = load_tsv(virus_summary_path)
    if virus_sum_df is not None:
        for cid in virus_sum_df['seq_name'].values:
            replicon_map[cid] = 'virus'

    plasmid_summary_path = resolve_path(results_dir, 'mge', 'genomad', 'plasmid_summary.tsv')
    plasmid_sum_df = load_tsv(plasmid_summary_path)
    if plasmid_sum_df is not None:
        for cid in plasmid_sum_df['seq_name'].values:
            replicon_map[cid] = 'plasmid'

    # Proviruses: contigs hosting integrated phage (overrides chromosome)
    provirus_path = resolve_path(results_dir, 'mge', 'genomad', 'provirus.tsv')
    provirus_df = load_tsv(provirus_path)
    if provirus_df is not None:
        for source_seq in provirus_df['source_seq'].values:
            if source_seq in replicon_map:
                replicon_map[source_seq] = 'provirus'

    counts = Counter(replicon_map.values())
    print(f"  Replicon: {dict(counts)} ({len(replicon_map)} contigs)")

    # geNomad viral taxonomy (semicolon-delimited lineage per virus contig)
    genomad_tax = {}
    genomad_tax_path = resolve_path(results_dir, 'mge', 'genomad', 'taxonomy.tsv')
    genomad_tax_df = load_tsv(genomad_tax_path)
    if genomad_tax_df is not None:
        for _, row in genomad_tax_df.iterrows():
            lineage = row.get('lineage', '')
            if pd.notna(lineage) and lineage:
                tax = parse_lineage(lineage)
                if any(tax.values()):
                    genomad_tax[row['seq_name']] = tax
        print(f"  geNomad taxonomy: {len(genomad_tax)} contigs")

    # Tiara eukaryotic classification
    tiara_map = {}
    tiara_path = resolve_path(results_dir, 'eukaryotic', 'tiara', 'tiara_output.tsv')
    tiara_df = load_tsv(tiara_path)
    if tiara_df is not None:
        for _, row in tiara_df.iterrows():
            tiara_map[row['sequence_id']] = row['class_fst_stage']
        print(f"  Tiara: {len(tiara_map)} contigs ({dict(Counter(tiara_map.values()))})")

    # Whokaryote classification
    whokaryote_map = {}
    whokaryote_path = resolve_path(results_dir, 'eukaryotic', 'whokaryote', 'whokaryote_classifications.tsv')
    whokaryote_df = load_tsv(whokaryote_path)
    if whokaryote_df is not None:
        for _, row in whokaryote_df.iterrows():
            whokaryote_map[row['contig']] = row.iloc[-1]  # 'predicted' column is last
        print(f"  Whokaryote: {len(whokaryote_map)} contigs ({dict(Counter(whokaryote_map.values()))})")

    depth_map = {}
    if depths_df is not None:
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

    # Build output records (without embeddings — those go in a separate file)
    contigs = []
    for i, cid in enumerate(contig_ids):
        contig_gc = gc_map.get(cid, 0)
        rec = {
            'id': cid,
            'length': int(len_map.get(cid, 0)),
            'bin': c2b.get(cid, ''),
            'depth': round(depth_map.get(cid, 0), 4),
            'gc': round(float(contig_gc), 2),
            'pca_x': round(float(pca_coords[i, 0]), 4),
            'pca_y': round(float(pca_coords[i, 1]), 4),
        }
        # Per-source taxonomy at all ranks
        for prefix, src_map in [('kaiju', kaiju_tax), ('kraken2', kraken2_tax), ('rrna', rrna_tax), ('sendsketch', sendsketch_tax)]:
            st = src_map.get(cid, {})
            for rank in RANKS:
                rec[f'{prefix}_{rank}'] = st.get(rank, '')
        # Per-binner bin assignments
        for binner, bmap in binner_maps.items():
            rec[f'{binner}_bin'] = bmap.get(cid, '')
        # geNomad replicon type and viral taxonomy
        rec['replicon'] = replicon_map.get(cid, '')
        for prefix, src_map in [('genomad', genomad_tax)]:
            st = src_map.get(cid, {})
            for rank in RANKS:
                rec[f'{prefix}_{rank}'] = st.get(rank, '')
        # Eukaryotic classifiers
        rec['tiara'] = tiara_map.get(cid, '')
        rec['whokaryote'] = whokaryote_map.get(cid, '')
        contigs.append(rec)

    explorer = {
        'contigs': contigs,
        'pca_variance_explained': [round(float(v), 4) for v in pca.explained_variance_ratio_[:5]],
        'transform': 'fourth_root',
    }

    # Embeddings written as separate files so --skip-tsne / --skip-umap are independent
    tsne_emb = None
    umap_emb = None

    if not skip_tsne:
        try:
            from sklearn.manifold import TSNE
            print("  Running t-SNE on fourth-root transformed TNF ...")
            tsne = TSNE(n_components=2, perplexity=30, random_state=42, max_iter=1000,
                        learning_rate='auto', init='pca')
            tsne_coords = tsne.fit_transform(pca_input)
            print("  t-SNE complete.")
            tsne_emb = {}
            for i, cid in enumerate(contig_ids):
                tsne_emb[cid] = [round(float(tsne_coords[i, 0]), 4),
                                 round(float(tsne_coords[i, 1]), 4)]
        except Exception as e:
            print(f"  [WARNING] t-SNE failed: {e}", file=sys.stderr)
    else:
        print("  --skip-tsne: skipping t-SNE (existing contig_tsne.json preserved)")

    if not skip_umap:
        try:
            import umap
            print("  Running UMAP on fourth-root transformed TNF ...")
            reducer = umap.UMAP(n_components=2, n_neighbors=15, min_dist=0.1,
                                metric='euclidean', random_state=42)
            umap_coords = reducer.fit_transform(pca_input)
            print("  UMAP complete.")
            umap_emb = {}
            for i, cid in enumerate(contig_ids):
                umap_emb[cid] = [round(float(umap_coords[i, 0]), 4),
                                 round(float(umap_coords[i, 1]), 4)]
        except Exception as e:
            print(f"  [WARNING] UMAP failed: {e}", file=sys.stderr)
    else:
        print("  --skip-umap: skipping UMAP (existing contig_umap.json preserved)")

    return explorer, tsne_emb, umap_emb


def main():
    args = parse_args()
    results_dir = args.results
    output_dir = args.output
    os.makedirs(output_dir, exist_ok=True)

    global _STORE_DIR
    _STORE_DIR = args.store_dir
    if _STORE_DIR:
        print(f"Store dir overlay: {_STORE_DIR}")

    print(f"Loading data from {results_dir} ...")

    # Load core data files
    _assembly_info_path = resolve_path(results_dir, 'assembly', 'assembly_info.txt')
    try:
        assembly_info = pd.read_csv(_assembly_info_path, sep='\t')
    except Exception:
        assembly_info = None
        print(f"  [WARNING] Missing: {_assembly_info_path}", file=sys.stderr)

    depths_df = load_tsv(resolve_path(results_dir, 'mapping', 'depths.txt'))
    contig2bin = load_tsv(resolve_path(results_dir, 'binning', 'dastool', 'contig2bin.tsv'),
                          header=None, names=['contig', 'bin'])
    dastool_summary = load_tsv(resolve_path(results_dir, 'binning', 'dastool', 'summary.tsv'))
    checkm2_df = load_tsv(resolve_path(results_dir, 'binning', 'checkm2', 'quality_report.tsv'))
    kaiju_df = load_tsv(resolve_path(results_dir, 'taxonomy', 'kaiju', 'kaiju_contigs.tsv'))
    sendsketch_df = load_tsv(resolve_path(results_dir, 'taxonomy', 'sendsketch', 'sendsketch_contigs.tsv'))
    virus_df = load_tsv(resolve_path(results_dir, 'mge', 'genomad', 'virus_summary.tsv'))
    plasmid_df = load_tsv(resolve_path(results_dir, 'mge', 'genomad', 'plasmid_summary.tsv'))
    defense_df = load_tsv(resolve_path(results_dir, 'mge', 'defensefinder', 'systems.tsv'))
    integron_df = load_tsv(resolve_path(results_dir, 'mge', 'integrons', 'integrons.tsv'),
                           comment='#')

    if assembly_info is None:
        print("[ERROR] Missing assembly_info.txt — cannot proceed without contig list")
        sys.exit(1)

    if depths_df is None:
        print("  [WARNING] depths.txt not found — depth/coverage fields will be empty", file=sys.stderr)
    if dastool_summary is None:
        print("  [WARNING] DAS_Tool summary not found — mags.json will be empty", file=sys.stderr)

    # 1. overview.json
    overview = build_overview(results_dir, assembly_info, depths_df, dastool_summary,
                              checkm2_df, virus_df, plasmid_df, contig2bin)
    write_json_gz(os.path.join(output_dir, 'overview.json'), overview)
    print(f"  Wrote overview.json")

    # 2. contig_lengths.json
    contig_lengths = build_contig_lengths(assembly_info, depths_df, results_dir)
    write_json_gz(os.path.join(output_dir, 'contig_lengths.json'), contig_lengths)
    print(f"  Wrote contig_lengths.json")

    # 3. mags.json
    annotation_df = load_annotation_tsv(results_dir)
    mags = build_mags(dastool_summary, checkm2_df, contig2bin, kaiju_df, depths_df,
                      virus_df, plasmid_df, defense_df, integron_df, annotation_df)
    write_json_gz(os.path.join(output_dir, 'mags.json'), mags)
    print(f"  Wrote mags.json ({len(mags)} MAGs)")

    # 3b. genes.json (compact gene features for contig detail view)
    print("Building genes.json ...")
    from genes_to_json import load_bakta, load_rrna, load_trna
    annot_candidates = [
        resolve_path(results_dir, 'annotation', 'bakta', 'extra', 'annotation.tsv'),
        resolve_path(results_dir, 'annotation', 'bakta', 'basic', 'annotation.tsv'),
        resolve_path(results_dir, 'annotation', 'prokka', 'annotation.tsv'),
    ]
    annot_tsv = next((p for p in annot_candidates if os.path.exists(p)), None)
    rrna_tsv = resolve_path(results_dir, 'taxonomy', 'rrna', 'rrna_genes.tsv')
    trna_tsv = resolve_path(results_dir, 'taxonomy', 'rrna', 'trna_genes.tsv')
    if annot_tsv:
        genes, n_bakta = load_bakta(annot_tsv)
        print(f"  Bakta: {n_bakta} features from {len(genes)} contigs")
        if os.path.isfile(rrna_tsv):
            rrna_genes, n_rrna = load_rrna(rrna_tsv)
            for contig, feats in rrna_genes.items():
                genes.setdefault(contig, []).extend(feats)
            print(f"  rRNA: {n_rrna} features from {len(rrna_genes)} contigs")
        if os.path.isfile(trna_tsv):
            trna_genes, n_trna = load_trna(trna_tsv)
            for contig, feats in trna_genes.items():
                genes.setdefault(contig, []).extend(feats)
            print(f"  tRNA/tmRNA: {n_trna} features from {len(trna_genes)} contigs")
        for contig in genes:
            genes[contig].sort(key=lambda f: f['s'])
        write_json_gz(os.path.join(output_dir, 'genes.json'), genes, separators=(',', ':'))
        print(f"  Wrote genes.json ({len(genes)} contigs)")
    else:
        write_json_gz(os.path.join(output_dir, 'genes.json'), {})
        print("  No annotation found — wrote empty genes.json")

    # 4. checkm2_all.json
    checkm2_all = build_checkm2_all(results_dir, checkm2_df, dastool_summary, contig2bin,
                                     kaiju_df, virus_df, plasmid_df, defense_df, integron_df,
                                     assembly_info=assembly_info, sendsketch_df=sendsketch_df)
    write_json_gz(os.path.join(output_dir, 'checkm2_all.json'), checkm2_all)
    print(f"  Wrote checkm2_all.json ({len(checkm2_all)} bins)")

    # 5. taxonomy_sunburst.json
    sunburst = build_taxonomy_sunburst(results_dir, kaiju_df, assembly_info, contig2bin, sendsketch_df)
    write_json_gz(os.path.join(output_dir, 'taxonomy_sunburst.json'), sunburst)
    print(f"  Wrote taxonomy_sunburst.json")

    # 6. kegg_heatmap.json
    kegg = build_kegg_heatmap(results_dir)
    write_json_gz(os.path.join(output_dir, 'kegg_heatmap.json'), kegg)
    print(f"  Wrote kegg_heatmap.json ({len(kegg['mag_ids'])} MAGs × {len(kegg['module_ids'])} modules)")

    # 6b. scg_heatmap.json
    scg = build_scg_heatmap(results_dir)
    write_json_gz(os.path.join(output_dir, 'scg_heatmap.json'), scg)
    bact_n = len(scg.get('bacteria', {}).get('mag_ids', []))
    bact_m = len(scg.get('bacteria', {}).get('module_ids', []))
    arch_n = len(scg.get('archaea', {}).get('mag_ids', []))
    arch_m = len(scg.get('archaea', {}).get('module_ids', []))
    print(f"  Wrote scg_heatmap.json (bacteria: {bact_n}×{bact_m}, archaea: {arch_n}×{arch_m})")

    # 7. coverage.json
    coverage = build_coverage(results_dir, dastool_summary, contig2bin, depths_df)
    write_json_gz(os.path.join(output_dir, 'coverage.json'), coverage)
    print(f"  Wrote coverage.json")

    # 8. mge_summary.json
    mge_summary = build_mge_summary(virus_df, plasmid_df, defense_df, integron_df)
    write_json_gz(os.path.join(output_dir, 'mge_summary.json'), mge_summary)
    print(f"  Wrote mge_summary.json")

    # 9. mge_per_bin.json
    mge_per_bin = build_mge_per_bin(contig2bin, virus_df, plasmid_df, defense_df, integron_df)
    write_json_gz(os.path.join(output_dir, 'mge_per_bin.json'), mge_per_bin)
    print(f"  Wrote mge_per_bin.json ({len(mge_per_bin)} bins)")

    # 10. eukaryotic.json
    eukaryotic = build_eukaryotic(results_dir, assembly_info)
    write_json_gz(os.path.join(output_dir, 'eukaryotic.json'), eukaryotic)
    print(f"  Wrote eukaryotic.json")

    # 11. contig_explorer.json (largest, do last)
    contig_explorer, tsne_emb, umap_emb = build_contig_explorer(
        results_dir, assembly_info, depths_df, contig2bin, kaiju_df,
        skip_tsne=args.skip_tsne, skip_umap=args.skip_umap, output_dir=output_dir,
        sendsketch_df=sendsketch_df)
    write_json_gz(os.path.join(output_dir, 'contig_explorer.json'), contig_explorer)
    print(f"  Wrote contig_explorer.json ({len(contig_explorer['contigs'])} contigs)")
    if tsne_emb is not None:
        write_json_gz(os.path.join(output_dir, 'contig_tsne.json'), tsne_emb,
                      separators=(',', ':'))
        print(f"  Wrote contig_tsne.json ({len(tsne_emb)} contigs)")
    if umap_emb is not None:
        write_json_gz(os.path.join(output_dir, 'contig_umap.json'), umap_emb,
                      separators=(',', ':'))
        print(f"  Wrote contig_umap.json ({len(umap_emb)} contigs)")

    print("\nDone! All JSON files written to", output_dir)


if __name__ == '__main__':
    main()
