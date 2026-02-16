#!/usr/bin/env python3
"""Convert per-MAG annotation TSVs to KEGG-Decoder input format.

KEGG-Decoder (Graham et al. 2018) expects a two-column tab-separated file:
    genome_proteinid<tab>KO

where 'genome' is everything before the first underscore in the protein ID.
Since our MAG IDs may contain underscores (e.g. DAS_Tool_bin_1), we replace
underscores with hyphens in the MAG ID so the first underscore correctly
delineates genome from protein: DAS-Tool-bin-1_contig00001_1

Note: KEGG-Decoder requires at least 3 genomes.

Usage:
    prepare_keggdecoder_input.py --input per_mag/ --output keggdecoder_input.tsv
"""

import argparse
import csv
import os
import sys


def sanitize_mag_id(mag_id):
    """Replace underscores with hyphens so KEGG-Decoder's split("_")[0] works.

    KEGG-Decoder uses info[0].split("_")[0] to extract the genome name,
    so the MAG ID must not contain underscores.
    """
    return mag_id.replace('_', '-')


def main():
    parser = argparse.ArgumentParser(
        description='Convert per-MAG annotation TSVs to KEGG-Decoder format')
    parser.add_argument('--input', required=True,
                        help='Directory of per-MAG annotation TSVs')
    parser.add_argument('--output', required=True,
                        help='Output TSV for KEGG-Decoder')
    args = parser.parse_args()

    tsv_files = sorted(f for f in os.listdir(args.input) if f.endswith('.tsv'))
    if not tsv_files:
        print("[WARNING] No per-MAG TSV files found", file=sys.stderr)
        with open(args.output, 'w') as fh:
            pass  # empty file
        return

    total_entries = 0
    mag_count = 0

    with open(args.output, 'w') as out_fh:
        for tsv_file in tsv_files:
            mag_id = tsv_file.replace('.tsv', '')
            safe_id = sanitize_mag_id(mag_id)
            tsv_path = os.path.join(args.input, tsv_file)
            mag_entries = 0

            with open(tsv_path) as in_fh:
                reader = csv.DictReader(in_fh, delimiter='\t')
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
                            # Format: sanitized-mag-id_protein_id<tab>KO
                            # KEGG-Decoder splits on first '_' to get genome name
                            out_fh.write(f"{safe_id}_{protein_id}\t{k}\n")
                            mag_entries += 1

            if mag_entries > 0:
                mag_count += 1
                total_entries += mag_entries

    if mag_count < 3:
        print(f"[WARNING] Only {mag_count} MAGs with KO assignments "
              f"(KEGG-Decoder requires >= 3)", file=sys.stderr)

    print(f"[INFO] Wrote {total_entries} KO entries for {mag_count} MAGs",
          file=sys.stderr)


if __name__ == '__main__':
    main()
