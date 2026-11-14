#!/usr/bin/env python3
"""Regression test for mag_analysis/bin/parse_marferret_results.py contig IDs.

DIAMOND's qseqid is the MetaEuk protein header, T_acc|C_acc|strand|..., so the
contig_id columns must carry C_acc, never the reference target T_acc. The first
header below is verbatim MetaEuk easy-predict output (MetaEuk commit 1f7d7023)
for a contig named my_contig_42 and a reference protein TARGET_ACC_1.
"""
import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[1] / 'mag_analysis/bin/parse_marferret_results.py'

HITS = [
    # (qseqid, expected contig_id)
    ('TARGET_ACC_1|my_contig_42|+|440|1.409e-130|1|300|899|300[300]:899[899]:600[600]',
     'my_contig_42'),
    ('TARGET_ACC_1|my_contig_42|-|120|2e-30|2|10|500|500[500]:400[400]:101[101]|200[200]:10[10]:191[191]',
     'my_contig_42'),
    ('REF9|odd|name|+|50|1e-5|1|1|90|1[1]:90[90]:90[90]', 'odd|name'),
    # Control: a plain protein ID (no strand token) keeps its first field.
    ('contigB_7', 'contigB_7'),
]


class MarferretContigIds(unittest.TestCase):
    def test_contig_id_is_metaeuk_contig_not_target(self):
        with tempfile.TemporaryDirectory(prefix='marferret-test-') as tmp:
            tmp = Path(tmp)
            diamond = tmp / 'hits.tsv'
            diamond.write_text(''.join(
                f'{q}\tSUBJ{i}\t90.0\t100\t1e-20\t200\tsubject {i}\n'
                for i, (q, _) in enumerate(HITS)))
            prot, contig = tmp / 'proteins.tsv', tmp / 'contigs.tsv'
            subprocess.run([sys.executable, str(SCRIPT), '--diamond', str(diamond),
                            '--out-proteins', str(prot), '--out-contigs', str(contig)],
                           check=True, capture_output=True)
            with open(prot) as f:
                rows = {r['protein_id']: r['contig_id']
                        for r in csv.DictReader(f, delimiter='\t')}
            self.assertEqual(rows, {q: c for q, c in HITS})
            with open(contig) as f:
                contigs = {r['contig_id']: int(r['n_proteins'])
                           for r in csv.DictReader(f, delimiter='\t')}
            self.assertEqual(contigs, {'my_contig_42': 2, 'odd|name': 1, 'contigB_7': 1})


if __name__ == '__main__':
    unittest.main()
