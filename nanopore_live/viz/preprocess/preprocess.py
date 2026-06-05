#!/usr/bin/env python3
"""Preprocess nanopore_live outputs into JSON for the sample dashboard.

Scans all dana.duckdb files in the input directory, aggregates per-sample
statistics, taxonomy, and functional profiles, then outputs JSON files
consumed by the Svelte frontend.
"""

import argparse
import gzip
import json
import os
import re
import sys
import math
from collections import defaultdict
from pathlib import Path

try:
    import duckdb
except ImportError:
    print("[ERROR] duckdb Python package required. Install with: pip install duckdb", file=sys.stderr)
    sys.exit(1)

try:
    import numpy as np
    HAS_NUMPY = True
except ImportError:
    HAS_NUMPY = False


def discover_samples(input_dir):
    """Find all sample directories with dana.duckdb files."""
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
                    'flowcell': flowcell,
                    'barcode': barcode_dir.name,
                    'id': f"{flowcell}:{barcode_dir.name}",
                    'db_path': str(db_path),
                    'dir': str(barcode_dir),
                })
    return samples


def extract_start_time(sample_dir):
    """Extract earliest start_time from the first non-empty fastq.gz in fq/."""
    fq_dir = Path(sample_dir) / 'fq'
    if not fq_dir.is_dir():
        return None

    # Sort by batch number suffix (e.g. _0.fastq.gz, _1.fastq.gz, ...)
    fq_files = sorted(
        fq_dir.glob('*.fastq.gz'),
        key=lambda f: int(re.search(r'_(\d+)\.fastq\.gz$', f.name).group(1))
            if re.search(r'_(\d+)\.fastq\.gz$', f.name) else 0
    )

    for fq in fq_files:
        if fq.stat().st_size < 100:
            continue  # skip empty gzip stubs
        try:
            with gzip.open(fq, 'rt') as f:
                header = f.readline().strip()
            m = re.search(r'start_time=(\S+)', header)
            if m:
                # Parse ISO timestamp, strip sub-second precision for JSON
                ts = m.group(1)
                # Normalize: keep YYYY-MM-DDTHH:MM:SS, drop fractional seconds + tz
                ts_clean = re.sub(r'\.\d+([+-]\d{2}:\d{2})?$', '', ts)
                return ts_clean
        except Exception:
            continue
    return None


def query_sample_stats(db_path):
    """Query per-sample statistics from DuckDB."""
    stats = {
        'read_count': 0, 'total_bases': 0, 'avg_length': 0, 'gc': 0,
        'n_classified': 0, 'n_genes': 0,
    }
    try:
        con = duckdb.connect(db_path, read_only=True)

        # Check which columns stats table has
        try:
            col_info = con.execute(
                "SELECT column_name FROM information_schema.columns WHERE table_name='stats'"
            ).fetchall()
            stat_cols = {c[0] for c in col_info}

            has_gc = 'gc' in stat_cols
            if has_gc:
                result = con.execute("SELECT COUNT(*), SUM(length), AVG(length), AVG(gc) FROM stats").fetchone()
            else:
                result = con.execute("SELECT COUNT(*), SUM(length), AVG(length) FROM stats").fetchone()

            if result and result[0]:
                stats['read_count'] = result[0]
                stats['total_bases'] = result[1] or 0
                stats['avg_length'] = round(result[2] or 0, 1)
                if has_gc and result[3] is not None:
                    stats['gc'] = round(result[3], 1)
        except Exception:
            pass

        # Classified reads
        try:
            result = con.execute("SELECT COUNT(*) FROM kraken").fetchone()
            stats['n_classified'] = result[0] if result else 0
        except Exception:
            pass

        # Gene count
        try:
            result = con.execute("SELECT COUNT(*) FROM prokka_annotations").fetchone()
            stats['n_genes'] = result[0] if result else 0
        except Exception:
            pass

        con.close()
    except Exception as e:
        print(f"  [WARNING] Failed to query {db_path}: {e}", file=sys.stderr)

    return stats


def query_taxonomy(db_path):
    """Get per-sample taxonomy at phylum and class level from krakenreport.

    Uses clade 'reads' count (all reads in subtree) rather than 'direct'
    count (reads at exactly that rank).

    Returns (phylum_tax, class_tax, class_to_phylum) where class_to_phylum
    maps each class name to its parent phylum name via the taxonomy lineage.
    """
    phylum_tax = {}
    class_tax = {}
    class_to_phylum = {}
    try:
        con = duckdb.connect(db_path, read_only=True)
        try:
            # Phylum level (clade counts)
            rows = con.execute("""
                SELECT name, SUM(reads) as count
                FROM krakenreport
                WHERE rank = 'P'
                GROUP BY name
                ORDER BY count DESC
            """).fetchall()
            for name, count in rows:
                clean_name = name.strip()
                if clean_name and count > 0:
                    phylum_tax[clean_name] = count

            # Class level (clade counts) + taxonomy lineage for parent phylum
            rows = con.execute("""
                SELECT name, SUM(reads) as count, MAX(taxonomy) as taxonomy
                FROM krakenreport
                WHERE rank = 'C'
                GROUP BY name
                ORDER BY count DESC
            """).fetchall()
            for name, count, taxonomy in rows:
                clean_name = name.strip()
                if clean_name and count > 0:
                    class_tax[clean_name] = count
                    # Extract parent phylum from lineage
                    if taxonomy:
                        parts = [p.strip() for p in taxonomy.split(';')]
                        for p in parts:
                            if p in phylum_tax:
                                class_to_phylum[clean_name] = p
                                break
        except Exception:
            pass  # No krakenreport table — leave phylum/class empty
        con.close()
    except Exception as e:
        print(f"  [WARNING] Failed to query taxonomy from {db_path}: {e}", file=sys.stderr)

    return phylum_tax, class_tax, class_to_phylum


# GTDB lineage prefix → rank code. No 'K' (kingdom) — GTDB starts at Domain.
GTDB_PREFIX_RANK = {'d': 'D', 'p': 'P', 'c': 'C', 'o': 'O', 'f': 'F', 'g': 'G', 's': 'S'}


def parse_gtdb_lineage(lineage):
    """Parse 'd__Bacteria;p__Foo;...;s__Bar' → list of (rank_code, name).

    Empty tokens (e.g. 'g__' with no name) are skipped; the lineage still
    threads through so deeper ranks reach their real ancestors.
    """
    out = []
    for tok in lineage.split(';'):
        tok = tok.strip()
        if len(tok) < 4 or tok[1:3] != '__':
            continue
        rank = GTDB_PREFIX_RANK.get(tok[0])
        if not rank:
            continue
        name = tok[3:].strip()
        if name:
            out.append((rank, name))
    return out


def query_taxonomy_gtdb(db_path):
    """Per-sample phylum/class counts from sendsketch GTDB classifications.

    Mirrors query_taxonomy() shape: (phylum_tax, class_tax, class_to_phylum).
    Sources from the sendsketch table's classified rows; empty if the table
    is missing or holds no classified reads.
    """
    phylum_tax = {}
    class_tax = {}
    class_to_phylum = {}
    try:
        con = duckdb.connect(db_path, read_only=True)
        try:
            rows = con.execute("""
                SELECT lineage, COUNT(*) AS count
                FROM sendsketch
                WHERE status = 'C'
                GROUP BY lineage
            """).fetchall()
            for lineage, count in rows:
                if not lineage or count <= 0:
                    continue
                phylum = None
                for rank, name in parse_gtdb_lineage(lineage):
                    if rank == 'P':
                        phylum_tax[name] = phylum_tax.get(name, 0) + count
                        phylum = name
                    elif rank == 'C':
                        class_tax[name] = class_tax.get(name, 0) + count
                        if phylum:
                            class_to_phylum[name] = phylum
        except Exception:
            pass  # No sendsketch table — leave empty
        con.close()
    except Exception as e:
        print(f"  [WARNING] Failed to query GTDB taxonomy from {db_path}: {e}", file=sys.stderr)
    return phylum_tax, class_tax, class_to_phylum


def query_sunburst_data_gtdb(db_path):
    """Full-lineage counts from sendsketch, in (path_counts, name_rank) shape
    compatible with build_taxonomy_sunburst() / build_taxonomy_parent_map().

    Leaf = deepest non-empty rank in the lineage. Taxonomy string = remaining
    ancestor names joined by ';', matching the Kraken path format.
    """
    path_counts = {}
    name_rank = {}
    try:
        con = duckdb.connect(db_path, read_only=True)
        try:
            rows = con.execute("""
                SELECT lineage, COUNT(*) AS cnt
                FROM sendsketch
                WHERE status = 'C'
                GROUP BY lineage
            """).fetchall()
            for lineage, cnt in rows:
                if not lineage or cnt <= 0:
                    continue
                parsed = parse_gtdb_lineage(lineage)
                if not parsed:
                    continue
                for rank, name in parsed:
                    name_rank[name] = rank
                leaf_name = parsed[-1][1]
                # Match Kraken's convention: `taxonomy` includes the leaf name
                # at the end. build_taxonomy_sunburst puts counts on
                # filtered[-1] and build_taxonomy_parent_map derives parent
                # edges from adjacent pairs in the filtered chain — both
                # require the leaf in the string to work past the genus.
                taxonomy = ';'.join(n for _, n in parsed)
                path_counts[(taxonomy, leaf_name)] = \
                    path_counts.get((taxonomy, leaf_name), 0) + cnt
        except Exception:
            pass
        con.close()
    except Exception as e:
        print(f"  [WARNING] Failed to query GTDB sunburst from {db_path}: {e}", file=sys.stderr)
    return path_counts, name_rank


def query_function(db_path):
    """Get per-sample functional annotation summary."""
    func = {
        'cds_count': 0, 'rrna_count': 0, 'trna_count': 0,
        'hypothetical_pct': 0, 'n_genes': 0, 'ec_counts': {},
    }
    try:
        con = duckdb.connect(db_path, read_only=True)
        try:
            # Feature type counts
            rows = con.execute("""
                SELECT ftype, COUNT(*) as count
                FROM prokka_annotations
                GROUP BY ftype
            """).fetchall()
            for ftype, count in rows:
                if ftype == 'CDS':
                    func['cds_count'] = count
                elif ftype == 'rRNA':
                    func['rrna_count'] = count
                elif ftype == 'tRNA':
                    func['trna_count'] = count
            func['n_genes'] = func['cds_count'] + func['rrna_count'] + func['trna_count']

            # Hypothetical protein ratio
            if func['cds_count'] > 0:
                hyp = con.execute("""
                    SELECT COUNT(*) FROM prokka_annotations
                    WHERE ftype = 'CDS' AND product LIKE '%hypothetical%'
                """).fetchone()
                func['hypothetical_pct'] = round((hyp[0] / func['cds_count']) * 100, 1) if hyp else 0

            # EC number counts
            ec_rows = con.execute("""
                SELECT ec_number, COUNT(*) as count
                FROM prokka_annotations
                WHERE ec_number IS NOT NULL AND ec_number != ''
                GROUP BY ec_number
                ORDER BY count DESC
                LIMIT 100
            """).fetchall()
            func['ec_counts'] = {ec: count for ec, count in ec_rows}
        except Exception:
            pass
        con.close()
    except Exception as e:
        print(f"  [WARNING] Failed to query function from {db_path}: {e}", file=sys.stderr)

    return func


def compute_grid_layout(samples_list):
    """Arrange samples in a rectangular grid sorted by start_time, flowcell, barcode."""
    sorted_samples = sorted(samples_list, key=lambda s: (
        s.get('start_time', '9999'),  # samples without timestamp sort last
        s.get('flowcell', ''),
        s.get('barcode', ''),
    ))
    n = len(sorted_samples)
    if n == 0:
        return {}
    cols = math.ceil(math.sqrt(n))
    rows = math.ceil(n / cols)
    grid = {}
    for i, s in enumerate(sorted_samples):
        row = i // cols
        col = i % cols
        # Center the grid around origin, space points 2 units apart
        grid[s['id']] = [
            (col - cols / 2) * 2,
            -(row - rows / 2) * 2,  # negative so top-left is first
        ]
    return grid


STANDARD_RANKS = frozenset({'R2', 'D', 'K', 'P', 'C', 'O', 'F', 'G', 'S'})


def query_sunburst_data(db_path):
    """Query full taxonomy lineage data for sunburst tree construction.

    Returns a dict mapping filtered taxonomy path tuples to direct read counts,
    plus a name->rank mapping for standard ranks.
    """
    path_counts = {}  # tuple(filtered_path) -> direct count
    name_rank = {}    # taxon_name -> rank letter
    try:
        con = duckdb.connect(db_path, read_only=True)
        try:
            rows = con.execute("""
                SELECT name, rank, taxonomy, SUM(CAST(direct AS INTEGER)) as cnt
                FROM krakenreport
                WHERE rank IN ('R2', 'D', 'K', 'P', 'C', 'O', 'F', 'G', 'S')
                GROUP BY name, rank, taxonomy
            """).fetchall()
            for name, rank, taxonomy, cnt in rows:
                clean = name.strip()
                name_rank[clean] = rank
                if cnt and cnt > 0 and taxonomy:
                    path_counts[(taxonomy, clean)] = \
                        path_counts.get((taxonomy, clean), 0) + cnt
        except Exception:
            pass
        con.close()
    except Exception:
        pass
    return path_counts, name_rank


def build_taxonomy_sunburst(all_path_counts, all_name_rank):
    """Build hierarchical sunburst tree from aggregated taxonomy data.

    Uses direct counts at each node (not clade) so the sunburst correctly
    sums children upward without double-counting.
    """
    # Build tree from filtered taxonomy paths
    tree = {}  # nested dict: name -> {'_count': int, '_children': {name: ...}}

    for (taxonomy, leaf_name), count in all_path_counts.items():
        parts = [p.strip() for p in taxonomy.split(';')]
        # Filter to standard-rank names only (+ the leaf itself)
        filtered = [p for p in parts if p in all_name_rank]
        if not filtered:
            continue

        # Walk down tree, creating nodes as needed
        node = tree
        for part in filtered:
            if part not in node:
                node[part] = {'_count': 0, '_children': {}}
            if part == filtered[-1]:
                node[part]['_count'] += count
            node = node[part]['_children']

    # Convert nested dict to D3-compatible sunburst format
    def to_sunburst(children_dict):
        result = []
        for name, data in children_dict.items():
            node = {'name': name}
            child_list = to_sunburst(data['_children'])
            if child_list:
                node['children'] = child_list
            if data['_count'] > 0:
                node['value'] = data['_count']
            # Skip nodes with no value and no children
            if 'children' in node or 'value' in node:
                result.append(node)
        # Sort by subtree total (descending)
        def subtree_total(n):
            v = n.get('value', 0)
            for c in n.get('children', []):
                v += subtree_total(c)
            return v
        result.sort(key=subtree_total, reverse=True)
        return result

    children = to_sunburst(tree)
    return {
        'name': 'Community',
        'children': children,
    }


def build_taxonomy_parent_map(all_path_counts, all_name_rank):
    """Build a global parent map from taxonomy lineage paths.

    Returns {taxon_name: parent_taxon_name} for all standard-rank taxa.
    """
    parents = {}
    for (taxonomy, leaf_name), _count in all_path_counts.items():
        parts = [p.strip() for p in taxonomy.split(';')]
        filtered = [p for p in parts if p in all_name_rank]
        for i, name in enumerate(filtered):
            if i > 0 and name not in parents:
                parents[name] = filtered[i - 1]
    return parents


def filter_sample_sunburst(path_counts, min_frac=0.001):
    """Filter per-sample sunburst data to taxa with >min_frac of total reads.

    Collapses (taxonomy, leaf_name) to {leaf_name: count} and filters.
    """
    # Sum direct counts per taxon name
    taxon_counts = {}
    for (taxonomy, leaf_name), count in path_counts.items():
        taxon_counts[leaf_name] = taxon_counts.get(leaf_name, 0) + count
    total = sum(taxon_counts.values()) or 1
    return {k: v for k, v in taxon_counts.items() if v / total >= min_frac}


def compute_sample_tsne(sketch_distances_file, sample_ids):
    """Compute t-SNE from comparesketch distance matrix.

    Returns dict {sample_id: [x, y]} or None if sketch data is unavailable.
    """
    if not sketch_distances_file or not os.path.exists(sketch_distances_file):
        print("[INFO] No sketch distances file — skipping sample t-SNE", file=sys.stderr)
        return None

    if not HAS_NUMPY:
        print("[WARNING] numpy not available — skipping sample t-SNE", file=sys.stderr)
        return None

    try:
        from sklearn.manifold import TSNE
    except ImportError:
        print("[WARNING] scikit-learn not available — skipping sample t-SNE", file=sys.stderr)
        return None

    # Parse sketch distances (format=3: query\tref\tANI\tsizeRatio)
    distances = {}
    name_to_sample = {}

    # Build mapping from sketch filename to sample ID
    # Sketch filenames use underscore (FAZ82916_barcode01), sample IDs use colon
    for sid in sample_ids:
        sketch_name = sid.replace(':', '_')
        name_to_sample[sketch_name] = sid

    with open(sketch_distances_file) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split('\t')
            if len(parts) < 3:
                continue
            q_name = parts[0].replace('.fa', '')
            r_name = parts[1].replace('.fa', '')
            try:
                ani = float(parts[2])
            except ValueError:
                continue

            q_sid = name_to_sample.get(q_name)
            r_sid = name_to_sample.get(r_name)
            if q_sid and r_sid and q_sid != r_sid:
                key = tuple(sorted([q_sid, r_sid]))
                dist = max(0, 100 - ani)
                distances[key] = min(distances.get(key, 100), dist)

    if not distances:
        print("[WARNING] No valid distances parsed from sketch file — skipping sample t-SNE", file=sys.stderr)
        return None

    # Build distance matrix
    id_list = sorted(sample_ids)
    n = len(id_list)

    # sklearn t-SNE requires perplexity < n_samples and n_samples >= 3 in
    # practice (perplexity floor of 2 ⇒ need ≥3 points). Skip cleanly while
    # the run is still ramping up — preprocess.py downstream of this must
    # not crash, otherwise samples.json never gets rewritten.
    if n < 3:
        print(f"[INFO] Only {n} sample(s) with sketches — skipping sample t-SNE", file=sys.stderr)
        return None

    id_to_idx = {sid: i for i, sid in enumerate(id_list)}
    dist_matrix = np.full((n, n), 100.0)  # max distance for unrelated pairs (no ANI hit)
    np.fill_diagonal(dist_matrix, 0)

    for (s1, s2), dist in distances.items():
        i, j = id_to_idx.get(s1), id_to_idx.get(s2)
        if i is not None and j is not None:
            dist_matrix[i, j] = dist
            dist_matrix[j, i] = dist

    # t-SNE on distance matrix
    perplexity = min(30, max(2, n // 3))
    tsne = TSNE(n_components=2, metric='precomputed', perplexity=perplexity,
                random_state=42, max_iter=1000, init='random')
    coords = tsne.fit_transform(dist_matrix)

    # If taller than wide, swap axes so the layout is landscape
    x_span = coords[:, 0].max() - coords[:, 0].min()
    y_span = coords[:, 1].max() - coords[:, 1].min()
    if y_span > x_span:
        coords = coords[:, ::-1]  # swap x and y

    return {id_list[i]: [float(coords[i, 0]), float(coords[i, 1])] for i in range(n)}


def normalize_barcode(raw):
    """Normalize barcode to 'barcodeNN' format.

    Accepts: '1', '01', 'barcode1', 'barcode01' -> 'barcode01'."""
    if not raw:
        return ''
    s = str(raw).strip().lower()
    m = re.match(r'^(?:barcode)?(\d+)$', s)
    if not m:
        return ''
    return f"barcode{int(m.group(1)):02d}"


COLUMN_ALIASES = {
    'latitude': 'lat', 'longitude': 'lon', 'long': 'lon', 'lng': 'lon',
}


def normalize_column(name):
    """Normalize a column header to its canonical name (case-insensitive)."""
    lower = name.strip().lower()
    return COLUMN_ALIASES.get(lower, lower)


def load_metadata(metadata_file):
    """Load user-provided metadata TSV.

    Requires flowcell and barcode columns.
    Canonical sample key: FLOWCELL:barcodeNN.
    Column aliases (case-insensitive): latitude->lat, longitude/long/lng->lon."""
    if not metadata_file or not os.path.exists(metadata_file):
        return {}

    metadata = {}
    with open(metadata_file) as f:
        raw_header = f.readline().strip().split('\t')
        header = [normalize_column(h) for h in raw_header]
        has_flowcell = 'flowcell' in header
        has_barcode = 'barcode' in header
        if not has_flowcell or not has_barcode:
            print("[ERROR] Metadata TSV requires both flowcell (or sample_id) and barcode columns",
                  file=sys.stderr)
            return {}

        for line in f:
            fields = line.strip().split('\t')
            row = dict(zip(header, fields))

            flowcell = (row.get('flowcell', '')).strip()
            barcode = normalize_barcode(row.get('barcode', ''))
            if not flowcell or not barcode:
                continue
            key = f"{flowcell}:{barcode}"

            # Parse numeric fields
            parsed = {}
            for k, v in row.items():
                if k in ('flowcell', 'barcode'):
                    continue
                try:
                    parsed[k] = float(v)
                except (ValueError, TypeError):
                    if v:
                        parsed[k] = v

            metadata[key] = parsed

    print(f"[INFO] Loaded metadata for {len(metadata)} samples", file=sys.stderr)
    return metadata


def build_mapping_payloads(samples, mapping_refs_dir, output_dir):
    """Build per-reference AIS-style JSON payloads from the mapping table.

    Output filename = <category_prefix><refname>.json, where category is
    looked up from <mapping_refs_dir>/<refname>/meta.json and mapped to
    'ais_' / 'hab_' / 'sar_'. References without meta.json default to
    'ais_' so the AIS view still finds them (matches the existing UI's
    no-fallback fetch convention).

    Per-sample stats include identity_hist (40 bins, 60-100% in 1% steps)
    and pos_hist (200 bins across the reference's total length, contig
    offsets from <refname>.contigs.json if available). Without contigs
    metadata we drop pos_hist; the view falls back to identity-only.
    """
    IDENTITY_LO = 60
    IDENTITY_HI = 100
    IDENTITY_BINS = 40
    POS_BINS = 200
    HQ_DEFAULT_PCT = 90  # Precomputed hq_hits cutoff; the UI slider
                         # recomputes live from identity_hist so this is
                         # just an initial scalar for sorters.

    # Discover all references seen anywhere in this run.
    refs = set()
    for s in samples:
        try:
            con = duckdb.connect(s['db_path'], read_only=True)
            for (r,) in con.execute(
                "SELECT DISTINCT reference FROM mapping WHERE reference IS NOT NULL"
            ).fetchall():
                refs.add(r)
            con.close()
        except Exception:
            # No mapping table on this sample's DB — skip.
            continue

    if not refs:
        print("[INFO] No mapping data found across samples; skipping AIS payloads",
              file=sys.stderr)
        return

    # Resolve reference metadata (category + contigs).
    ref_meta = {}
    for ref in refs:
        meta = {'category': 'AIS', 'contigs': None, 'total_length': 0}
        if mapping_refs_dir:
            mdir = Path(mapping_refs_dir) / ref
            mj = mdir / 'meta.json'
            if mj.exists():
                try:
                    md = json.load(open(mj))
                    # 'SAR-proxy' / 'SAR' / 'AIS' / 'HAB' → uppercase, strip
                    # the proxy qualifier.
                    cat = (md.get('category') or 'AIS').split('-')[0].upper()
                    if cat in ('AIS', 'HAB', 'SAR'):
                        meta['category'] = cat
                except Exception as e:
                    print(f"  [WARNING] meta.json {ref}: {e}", file=sys.stderr)
            cj = mdir / f'{ref}.contigs.json'
            if cj.exists():
                try:
                    cd = json.load(open(cj))
                    meta['contigs'] = cd.get('contigs', [])
                    meta['total_length'] = cd.get('total_length',
                        sum(c.get('length', 0) for c in meta['contigs']))
                except Exception as e:
                    print(f"  [WARNING] contigs.json {ref}: {e}", file=sys.stderr)
        ref_meta[ref] = meta

    # Build per-ref payloads.
    for ref in sorted(refs):
        meta = ref_meta[ref]
        contigs = meta['contigs']
        total_length = meta['total_length']
        # Lookup contig name → offset for the per-hit position bin.
        contig_offset = {c['name']: c['offset'] for c in (contigs or [])}
        use_pos_hist = total_length > 0

        counts = {}
        stats = {}
        for s in samples:
            try:
                con = duckdb.connect(s['db_path'], read_only=True)
                rows = con.execute(
                    "SELECT identity_pct, mapq, rname, pos FROM mapping "
                    "WHERE reference = ?",
                    [ref]
                ).fetchall()
                con.close()
            except Exception:
                continue
            if not rows:
                continue
            sid = s['id']
            n = len(rows)
            counts[sid] = n
            sum_id = 0.0
            sum_mapq = 0.0
            id_hist = [0] * IDENTITY_BINS
            pos_hist = [0] * POS_BINS if use_pos_hist else None
            for ident, mapq, rname, pos in rows:
                sum_id += ident
                sum_mapq += mapq
                # Identity bin: clamp out-of-range to the edge bins so
                # the histogram totals match n exactly.
                bi = int(ident - IDENTITY_LO)
                if bi < 0:
                    bi = 0
                elif bi >= IDENTITY_BINS:
                    bi = IDENTITY_BINS - 1
                id_hist[bi] += 1
                if pos_hist is not None:
                    off = contig_offset.get(rname, 0)
                    gp = off + pos
                    pi = int(gp / total_length * POS_BINS)
                    if pi < 0:
                        pi = 0
                    elif pi >= POS_BINS:
                        pi = POS_BINS - 1
                    pos_hist[pi] += 1
            # Initial hq_hits at the 90% default cutoff; UI slider updates live.
            hq_bin_idx = max(0, min(IDENTITY_BINS, HQ_DEFAULT_PCT - IDENTITY_LO))
            hq_hits = sum(id_hist[hq_bin_idx:])
            stat = {
                'n': n,
                'hq_hits': hq_hits,
                'mean_identity': sum_id / n,
                'mean_mapq': sum_mapq / n,
                'identity_hist': id_hist,
            }
            if pos_hist is not None:
                stat['pos_hist'] = pos_hist
            stats[sid] = stat

        payload = {
            'counts': counts,
            'stats': stats,
            'identity_bins': {
                'lo_pct': IDENTITY_LO, 'hi_pct': IDENTITY_HI, 'n': IDENTITY_BINS
            },
        }
        if use_pos_hist:
            n_contigs = len(contigs or [])
            payload['pos_bins'] = {
                'n': POS_BINS,
                'total_length': total_length,
                'n_contigs': n_contigs,
            }
            # `contigs` array is only useful for chromosome-scale refs;
            # for fragmented refs the view renders a single linear axis.
            if contigs and n_contigs <= 64:
                payload['contigs'] = contigs

        prefix = meta['category'].lower() + '_'
        filename = f'{prefix}{ref}.json'
        path = os.path.join(output_dir, filename)
        with open(path, 'w') as f:
            json.dump(payload, f, separators=(',', ':'))
        size = os.path.getsize(path)
        total_hits = sum(counts.values())
        print(f"  {filename}: {size / 1024:.1f} KB "
              f"({len(counts)} samples, {total_hits} hits, "
              f"category={meta['category']})", file=sys.stderr)


def shannon_diversity(taxonomy):
    """Compute Shannon diversity index from taxonomy counts."""
    total = sum(taxonomy.values())
    if total == 0:
        return 0
    h = 0
    for count in taxonomy.values():
        if count > 0:
            p = count / total
            h -= p * math.log(p)
    return round(h, 3)


def main():
    parser = argparse.ArgumentParser(description='Preprocess nanopore_live outputs for viz')
    parser.add_argument('--input', required=True, help='Directory with nanopore_live outputs')
    parser.add_argument('--output', required=True, help='Output directory for JSON files')
    parser.add_argument('--metadata', help='Sample metadata TSV file')
    parser.add_argument('--sketch-distances', help='comparesketch distances TSV')
    parser.add_argument('--mapping-refs',
        help='Directory of reference subdirs (matches --mapping_refs in main.nf). '
             'Used to enrich the AIS/HAB/SAR payloads with contig offsets and '
             'category labels from per-ref meta.json.')
    args = parser.parse_args()

    # Fall back to a marker file the launcher dropped in the outdir at start —
    # avoids needing to thread --mapping-refs through every caller, and lets
    # us pick up the path for already-running DB_SYNC loops where the
    # Nextflow script body was templated before the flag existed.
    if not args.mapping_refs:
        marker = Path(args.input) / '.pipeline_state' / 'mapping_refs'
        if marker.exists():
            try:
                p = marker.read_text().strip()
                if p and Path(p).is_dir():
                    args.mapping_refs = p
                    print(f"[INFO] mapping refs auto-detected from {marker} -> {p}",
                          file=sys.stderr)
            except Exception:
                pass

    os.makedirs(args.output, exist_ok=True)

    # Discover samples
    print("[INFO] Discovering samples...", file=sys.stderr)
    samples = discover_samples(args.input)
    print(f"[INFO] Found {len(samples)} samples", file=sys.stderr)

    if not samples:
        print("[ERROR] No samples found. Check --input directory.", file=sys.stderr)
        sys.exit(1)

    # Load metadata
    metadata = load_metadata(args.metadata)

    # Query each sample — collect Kraken and GTDB taxonomy payloads in parallel.
    # Each source produces its own phylum/class/sunburst dicts; they're emitted
    # side-by-side in the final JSON under {kraken: ..., gtdb: ...} so the viz
    # can flip between taxonomies without a re-fetch.
    all_stats = {}
    all_function = {}
    # Kraken payloads
    all_phylum_tax = {}
    all_class_tax = {}
    all_sunburst_paths = {}   # (taxonomy, leaf_name) -> total proportional count
    all_sunburst_names = {}   # name -> rank
    all_sample_sunburst = {}  # sample_id -> {taxon: direct_count} (filtered >0.1%)
    # GTDB payloads (same shape, sourced from sendsketch)
    all_phylum_tax_gtdb = {}
    all_class_tax_gtdb = {}
    all_sunburst_paths_gtdb = {}
    all_sunburst_names_gtdb = {}
    all_sample_sunburst_gtdb = {}

    for i, s in enumerate(samples):
        sid = s['id']
        print(f"  [{i+1}/{len(samples)}] {sid}...", file=sys.stderr, end='', flush=True)

        stats = query_sample_stats(s['db_path'])
        phylum_tax, class_tax, class_to_phylum = query_taxonomy(s['db_path'])
        phylum_tax_g, class_tax_g, class_to_phylum_g = query_taxonomy_gtdb(s['db_path'])
        func = query_function(s['db_path'])

        # Collect full taxonomy lineage for sunburst (Kraken + GTDB)
        sb_paths, sb_names = query_sunburst_data(s['db_path'])
        all_sunburst_names.update(sb_names)
        # Normalize per-sample counts to proportions (equal weight per sample)
        sb_total = sum(sb_paths.values()) or 1
        for key, cnt in sb_paths.items():
            all_sunburst_paths[key] = all_sunburst_paths.get(key, 0) + cnt / sb_total
        all_sample_sunburst[sid] = filter_sample_sunburst(sb_paths)

        sb_paths_g, sb_names_g = query_sunburst_data_gtdb(s['db_path'])
        all_sunburst_names_gtdb.update(sb_names_g)
        sb_total_g = sum(sb_paths_g.values()) or 1
        for key, cnt in sb_paths_g.items():
            all_sunburst_paths_gtdb[key] = all_sunburst_paths_gtdb.get(key, 0) + cnt / sb_total_g
        all_sample_sunburst_gtdb[sid] = filter_sample_sunburst(sb_paths_g)

        # Extract start time from first fastq
        start_time = extract_start_time(s['dir'])
        if start_time:
            stats['start_time'] = start_time

        # Derived summary fields on samples.json — prefer Kraken when present
        # (back-compat), else fall back to GTDB so samples with only the new
        # sendsketch data still show a dominant phylum/class.
        primary_phylum = phylum_tax or phylum_tax_g
        primary_class = class_tax or class_tax_g
        primary_c2p = class_to_phylum if phylum_tax else class_to_phylum_g
        stats['diversity'] = shannon_diversity(primary_phylum)
        if primary_class:
            dominant_class = max(primary_class, key=primary_class.get)
            stats['dominant_class'] = dominant_class
            if dominant_class in primary_c2p:
                stats['dominant_phylum'] = primary_c2p[dominant_class]
            elif primary_phylum:
                stats['dominant_phylum'] = max(primary_phylum, key=primary_phylum.get)
        elif primary_phylum:
            stats['dominant_phylum'] = max(primary_phylum, key=primary_phylum.get)

        all_stats[sid] = stats
        all_phylum_tax[sid] = phylum_tax
        all_class_tax[sid] = class_tax
        all_phylum_tax_gtdb[sid] = phylum_tax_g
        all_class_tax_gtdb[sid] = class_tax_g
        all_function[sid] = func

        ts_str = f" [{start_time[:10]}]" if start_time else ""
        print(f" {stats['read_count']} reads{ts_str}", file=sys.stderr)

    # Build samples.json
    samples_list = []
    for s in samples:
        sid = s['id']
        entry = {'id': sid, 'flowcell': s['flowcell'], 'barcode': s['barcode']}
        entry.update(all_stats.get(sid, {}))
        samples_list.append(entry)

    # Filter out empty samples (no reads)
    samples_list = [s for s in samples_list if s.get('read_count', 0) > 0]
    sample_ids = [s['id'] for s in samples_list]

    print(f"[INFO] {len(samples_list)} samples with data", file=sys.stderr)

    # Compute sample embeddings
    sample_tsne = compute_sample_tsne(args.sketch_distances, sample_ids)
    print("[INFO] Computing grid layout...", file=sys.stderr)
    sample_grid = compute_grid_layout(samples_list)

    # Build overview
    n_flowcells = len(set(s.get('flowcell', '') for s in samples_list))
    overview = {
        'n_flowcells': n_flowcells,
        'n_samples': len(samples_list),
        'total_reads': sum(s.get('read_count', 0) for s in samples_list),
        'total_bases': sum(s.get('total_bases', 0) for s in samples_list),
    }

    # Build per-source taxonomy sunburst trees + parent maps
    print("[INFO] Building taxonomy sunbursts (kraken + gtdb)...", file=sys.stderr)
    sunburst_kraken = build_taxonomy_sunburst(all_sunburst_paths, all_sunburst_names)
    sunburst_gtdb   = build_taxonomy_sunburst(all_sunburst_paths_gtdb, all_sunburst_names_gtdb)
    parent_map_kraken = build_taxonomy_parent_map(all_sunburst_paths, all_sunburst_names)
    parent_map_gtdb   = build_taxonomy_parent_map(all_sunburst_paths_gtdb, all_sunburst_names_gtdb)

    # Write outputs
    def write_json(filename, data):
        path = os.path.join(args.output, filename)
        with open(path, 'w') as f:
            json.dump(data, f, separators=(',', ':'))
        size = os.path.getsize(path)
        print(f"  {filename}: {size / 1024:.1f} KB", file=sys.stderr)

    write_json('overview.json', overview)
    write_json('samples.json', samples_list)
    embeddings = {'grid': sample_grid}
    if sample_tsne:
        embeddings['tsne'] = sample_tsne
    write_json('sample_tsne.json', embeddings)

    def build_tax_payload(parents, names_rank, phylum, cls, sample_sunburst):
        """Shape one taxonomy source into the viz-consumed {parents,ranks,samples}."""
        return {
            'parents': parents,
            'ranks': dict(names_rank),
            'samples': {
                sid: {
                    'phylum': phylum.get(sid, {}),
                    'class':  cls.get(sid, {}),
                    'direct': sample_sunburst.get(sid, {}),
                }
                for sid in sample_ids
            },
        }

    # Taxonomy payloads — {kraken: ..., gtdb: ...}. Viz flips between the two
    # via a global source toggle; old flat-shape consumers should fall back to
    # the `kraken` branch.
    tax_data = {
        'kraken': build_tax_payload(parent_map_kraken, all_sunburst_names,
                                    all_phylum_tax, all_class_tax, all_sample_sunburst),
        'gtdb':   build_tax_payload(parent_map_gtdb, all_sunburst_names_gtdb,
                                    all_phylum_tax_gtdb, all_class_tax_gtdb, all_sample_sunburst_gtdb),
    }
    sunburst = {'kraken': sunburst_kraken, 'gtdb': sunburst_gtdb}
    write_json('sample_taxonomy.json', tax_data)
    write_json('taxonomy_sunburst.json', sunburst)
    write_json('sample_function.json', {k: v for k, v in all_function.items() if k in sample_ids})
    write_json('metadata.json', metadata)

    # Per-reference AIS/HAB/SAR payloads from the mapping table. Writes one
    # JSON per reference into the same output dir; the SPA's AIS/HAB/SAR
    # views fetch them on demand by name.
    print("[INFO] Building reference-mapping payloads (ais_/hab_/sar_)...", file=sys.stderr)
    build_mapping_payloads(samples, args.mapping_refs, args.output)

    print("[SUCCESS] Main preprocessing complete", file=sys.stderr)


if __name__ == '__main__':
    main()
