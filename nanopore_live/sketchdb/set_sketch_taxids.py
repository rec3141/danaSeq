#!/usr/bin/env python3
"""Stamp GTDB taxids into BBTools sketch headers by genome accession.

Usage: set_sketch_taxids.py accession2taxid.tsv IN.sketch OUT.sketch

Each sketch header (#SZ:... line) names its genome file in FN: (for example
FN:GCA_000008085.1_genomic.fna.gz). The accession from that file name is looked
up in accession2taxid.tsv (from build_gtdb_taxdump.py), and the header's ID:
field is set to that GTDB species taxid; NM0: is set to the accession so the
server reports which genome matched. A sketch whose accession is not in the map
is an error, since every GTDB representative has a taxonomy.
"""
import re
import sys

ACC = re.compile(r'(GC[AF]_\d{9}\.\d+)')


def main():
    map_path, in_path, out_path = sys.argv[1:]
    taxid = {}
    with open(map_path) as f:
        next(f)
        for line in f:
            acc, tid, _ = line.rstrip('\n').split('\t')
            taxid[acc] = tid
    n = 0
    missing = []
    with open(in_path) as fin, open(out_path, 'w') as fout:
        for line in fin:
            if line.startswith('#SZ:'):
                fields = line.rstrip('\n').split('\t')
                fn = next((x[3:] for x in fields if x.startswith('FN:')), '')
                m = ACC.search(fn)
                tid = taxid.get(m.group(1)) if m else None
                if tid is None:
                    missing.append(fn)
                else:
                    # A sketch made without a tree has no ID: field; it goes
                    # before FN:, where SketchMaker writes it.
                    out = []
                    for x in fields:
                        if x.startswith('ID:'):
                            continue
                        if x.startswith('FN:'):
                            out.append('ID:' + tid)
                        elif x.startswith('NM0:'):
                            x = 'NM0:' + m.group(1)
                        out.append(x)
                    line = '\t'.join(out) + '\n'
                    n += 1
            fout.write(line)
    print(f'{n} sketches stamped, {len(missing)} without a GTDB accession', file=sys.stderr)
    if missing:
        print('  e.g. ' + ', '.join(missing[:5]), file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()
