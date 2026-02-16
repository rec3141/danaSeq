#!/usr/bin/env python3
"""Merge KofamScan, eggNOG-mapper, and dbCAN annotations into a unified per-protein TSV.

Joins on protein_id across three annotation sources:
  - KofamScan:  KEGG Orthology assignments (adaptive threshold filtering)
  - eggNOG-mapper: COG, GO, EC, KEGG, Pfam, description
  - dbCAN:      CAZyme consensus (requires agreement from >= 2/3 methods)

Output columns:
  protein_id, contig_id, KO, COG_category, GOs, EC, KEGG_Pathway, PFAMs, CAZy, description

Usage:
    merge_annotations.py --kofamscan kofamscan_results.tsv \
                         --emapper emapper_results.emapper.annotations \
                         --dbcan overview.txt \
                         -o merged_annotations.tsv
"""

import argparse
import csv
import sys


def parse_kofamscan(path):
    """Parse KofamScan detail-tsv output.

    KofamScan detail-tsv format:
      Column 0: significance mark (* = above threshold)
      Column 1: gene name (protein_id)
      Column 2: KO number
      Column 3: thrshld (adaptive threshold)
      Column 4: score
      Column 5: E-value
      Column 6: KO definition

    Only rows marked with '*' (above adaptive threshold) are kept.
    """
    results = {}
    if not path:
        return results

    try:
        with open(path) as fh:
            for line in fh:
                line = line.rstrip('\n')
                if line.startswith('#') or not line.strip():
                    continue
                parts = line.split('\t')
                if len(parts) < 5:
                    # Try space-delimited (KofamScan default detail format)
                    parts = line.split(None, 6)
                if len(parts) < 5:
                    continue

                significance = parts[0].strip()
                if significance != '*':
                    continue

                protein_id = parts[1].strip()
                ko = parts[2].strip()

                # Keep first (best) hit per protein
                if protein_id not in results:
                    results[protein_id] = ko
    except FileNotFoundError:
        print(f"[WARNING] KofamScan file not found: {path}", file=sys.stderr)
    except Exception as e:
        print(f"[WARNING] Error parsing KofamScan: {e}", file=sys.stderr)

    return results


def parse_emapper(path):
    """Parse eggNOG-mapper .emapper.annotations file.

    eggNOG-mapper v2.1 output columns (21 columns, 0-indexed):
      0:  #query
      1:  seed_ortholog
      2:  evalue
      3:  score
      4:  eggNOG_OGs
      5:  max_annot_lvl
      6:  COG_category
      7:  Description
      8:  Preferred_name
      9:  GOs
     10:  EC
     11:  KEGG_ko
     12:  KEGG_Pathway
     13:  KEGG_Module
     14:  KEGG_Reaction
     15:  KEGG_rclass
     16:  BRITE
     17:  KEGG_TC
     18:  CAZy
     19:  BiGG_Reaction
     20:  PFAMs
    """
    results = {}
    if not path:
        return results

    try:
        with open(path) as fh:
            for line in fh:
                if line.startswith('#') or line.startswith('##'):
                    continue
                line = line.rstrip('\n')
                if not line:
                    continue

                parts = line.split('\t')
                if len(parts) < 21:
                    continue

                protein_id = parts[0].strip()
                results[protein_id] = {
                    'COG_category': parts[6] if parts[6] != '-' else '',
                    'description':  parts[7] if parts[7] != '-' else '',
                    'GOs':          parts[9] if parts[9] != '-' else '',
                    'EC':           parts[10] if parts[10] != '-' else '',
                    'KEGG_ko':      parts[11] if parts[11] != '-' else '',
                    'KEGG_Pathway': parts[12] if parts[12] != '-' else '',
                    'PFAMs':        parts[20] if parts[20] != '-' else '',
                }
    except FileNotFoundError:
        print(f"[WARNING] eggNOG-mapper file not found: {path}", file=sys.stderr)
    except Exception as e:
        print(f"[WARNING] Error parsing eggNOG-mapper: {e}", file=sys.stderr)

    return results


def parse_dbcan(path):
    """Parse dbCAN overview.txt — keep assignments with >= 2/3 method agreement.

    overview.txt columns:
      Gene ID, HMMER, dbCAN_sub, DIAMOND, #ofTools
    """
    results = {}
    if not path:
        return results

    try:
        with open(path) as fh:
            reader = csv.DictReader(fh, delimiter='\t')
            for row in reader:
                gene_id = row.get('Gene ID', '').strip()
                n_tools = row.get('#ofTools', '0').strip()

                if not gene_id:
                    continue

                try:
                    if int(n_tools) < 2:
                        continue
                except ValueError:
                    continue

                # Collect CAZy family assignments from the methods that hit
                # run_dbcan v4+ renamed HMMER → dbCAN_hmm
                cazymes = set()
                for col in ['HMMER', 'dbCAN_hmm', 'dbCAN_sub', 'DIAMOND']:
                    val = row.get(col, '-').strip()
                    if val and val != '-':
                        # May contain multiple families separated by +
                        for fam in val.split('+'):
                            fam = fam.strip()
                            # Strip sub-family detail (e.g. GH13_20 -> GH13_20)
                            if fam and fam != '-':
                                cazymes.add(fam)

                if cazymes:
                    results[gene_id] = ','.join(sorted(cazymes))
    except FileNotFoundError:
        print(f"[WARNING] dbCAN file not found: {path}", file=sys.stderr)
    except Exception as e:
        print(f"[WARNING] Error parsing dbCAN: {e}", file=sys.stderr)

    return results


def extract_contig_id(protein_id):
    """Extract contig ID from Prokka/Bakta protein IDs.

    Prokka:  contig_00001_1  -> contig_00001
    Bakta:   contig_00001_00001  -> contig_00001
    General: take everything up to the last underscore + digits
    """
    parts = protein_id.rsplit('_', 1)
    if len(parts) == 2 and parts[1].isdigit():
        return parts[0]
    return protein_id


def main():
    parser = argparse.ArgumentParser(
        description='Merge KofamScan, eggNOG-mapper, and dbCAN annotations')
    parser.add_argument('--kofamscan', help='KofamScan detail-tsv output')
    parser.add_argument('--emapper', help='eggNOG-mapper .emapper.annotations file')
    parser.add_argument('--dbcan', help='dbCAN overview.txt file')
    parser.add_argument('-o', '--output', default='-',
                        help='Output TSV (default: stdout)')
    args = parser.parse_args()

    if not any([args.kofamscan, args.emapper, args.dbcan]):
        print("[ERROR] At least one annotation source is required", file=sys.stderr)
        sys.exit(1)

    # Parse all sources
    ko_map = parse_kofamscan(args.kofamscan)
    emapper_map = parse_emapper(args.emapper)
    dbcan_map = parse_dbcan(args.dbcan)

    print(f"[INFO] KofamScan: {len(ko_map)} proteins with KO assignments", file=sys.stderr)
    print(f"[INFO] eggNOG-mapper: {len(emapper_map)} proteins annotated", file=sys.stderr)
    print(f"[INFO] dbCAN: {len(dbcan_map)} proteins with CAZyme assignments", file=sys.stderr)

    # Collect all protein IDs
    all_proteins = sorted(set(ko_map.keys()) | set(emapper_map.keys()) | set(dbcan_map.keys()))
    print(f"[INFO] Total unique proteins: {len(all_proteins)}", file=sys.stderr)

    # Write merged output
    out_fh = open(args.output, 'w') if args.output != '-' else sys.stdout
    header = ['protein_id', 'contig_id', 'KO', 'COG_category', 'GOs', 'EC',
              'KEGG_Pathway', 'PFAMs', 'CAZy', 'description']
    print('\t'.join(header), file=out_fh)

    for protein_id in all_proteins:
        contig_id = extract_contig_id(protein_id)
        ko = ko_map.get(protein_id, '')
        em = emapper_map.get(protein_id, {})
        cazy = dbcan_map.get(protein_id, '')

        # Use eggNOG KEGG_ko as fallback when KofamScan KO is absent
        if not ko:
            ko = em.get('KEGG_ko', '')

        row = [
            protein_id,
            contig_id,
            ko,
            em.get('COG_category', ''),
            em.get('GOs', ''),
            em.get('EC', ''),
            em.get('KEGG_Pathway', ''),
            em.get('PFAMs', ''),
            cazy,
            em.get('description', ''),
        ]
        print('\t'.join(row), file=out_fh)

    if args.output != '-':
        out_fh.close()

    print(f"[INFO] Wrote {len(all_proteins)} merged annotations", file=sys.stderr)


if __name__ == '__main__':
    main()
