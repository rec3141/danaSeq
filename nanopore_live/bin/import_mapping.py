"""Import reference-mapping output into dana.duckdb.

Called from import_all.py with the per-barcode directory as cwd. Reads
`map/*.txt` — each file is a filtered SAM-ish stream produced by
modules/mapping.nf (filter_minimap2.awk: headers stripped, mapq>=1 and
aligned_len>=10). The file basename (sans .txt) is the canonical
reference name and lands in the `mapping.reference` column.

Idempotent via the shared `import_log` table — files keyed as
`map/<refname>.txt`. Re-importing the same key replaces that
reference's rows (DELETE + INSERT) so updates pick up cleanly.
"""

import glob
import os
import re
import sys


CIGAR_RE = re.compile(r'(\d+)([MIDNSHP=X])')
TAG_INT_RE = re.compile(r'\t(NM|AS):i:(-?\d+)')
TAG_FLT_RE = re.compile(r'\tde:f:([0-9.]+)')


def _aligned_length(cigar):
    n = 0
    for length, op in CIGAR_RE.findall(cigar):
        if op in 'M=X':
            n += int(length)
    return n


def _parse_txt(path, reference):
    """Yield row tuples for the `mapping` table."""
    with open(path, 'rb') as fh:
        for raw in fh:
            line = raw.decode('utf-8', 'replace').rstrip('\n')
            f = line.split('\t')
            if len(f) < 6:
                continue
            qname, flag, rname, pos, mapq, cigar = f[0], f[1], f[2], f[3], f[4], f[5]
            try:
                flag_i = int(flag); pos_i = int(pos); mapq_i = int(mapq)
            except ValueError:
                continue
            nm = as_score = None
            de = None
            for tag, val in TAG_INT_RE.findall(line):
                if tag == 'NM':
                    nm = int(val)
                elif tag == 'AS':
                    as_score = int(val)
            m = TAG_FLT_RE.search(line)
            if m:
                de = float(m.group(1))
            ident = (1.0 - de) * 100.0 if de is not None else None
            yield (reference, qname, flag_i, rname, pos_i, mapq_i, cigar,
                   nm, as_score, de, ident, _aligned_length(cigar))


def import_mapping(con, imported):
    files = sorted(glob.glob('map/*.txt'))
    pending = [f for f in files if f not in imported]
    if not pending:
        return

    for path in pending:
        ref = os.path.splitext(os.path.basename(path))[0]
        rows = list(_parse_txt(path, ref))
        # Replace any existing rows for this reference (idempotent updates).
        con.execute('DELETE FROM mapping WHERE reference = ?', [ref])
        if rows:
            con.executemany(
                'INSERT INTO mapping VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)',
                rows,
            )
        con.execute(
            'INSERT INTO import_log (filename) VALUES (?) ON CONFLICT DO NOTHING',
            [path],
        )
        print(f'[import_mapping] {path}: {len(rows)} rows', file=sys.stderr)
