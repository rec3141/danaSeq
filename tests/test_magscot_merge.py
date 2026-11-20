#!/usr/bin/env python3
"""Regression test for bin merging in mag_analysis/bin/magscot.py.

Iteration 1 merges original bins from different binners; each later iteration
merges the hybrid bins of the previous iteration (seeds) with original bins
again, and a hybrid bin is the seed's contigs plus the partner's contigs
(upstream MAGScoT, R/merge.R). Needs pandas.

Toy input, nine single-copy markers M1..M9 spread over five contigs:

    set1 b1 = c1 c2 c5     set2 b2 = c2 c3     set3 b3 = c1 c3 c4

Iteration 1 builds MAGScoT_1_1 = b1 + b2 = c1 c2 c3 c5 (6 markers). Iteration 2
seeds on MAGScoT_1_1, finds partner b3 through c1/c3, and must build
MAGScoT_2_1 = c1 c2 c3 c4 c5, the only bin carrying all nine markers.
"""
import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[1] / 'mag_analysis/bin/magscot.py'

C2B = ('b1\tc1\tset1\nb1\tc2\tset1\nb1\tc5\tset1\n'
       'b2\tc2\tset2\nb2\tc3\tset2\n'
       'b3\tc1\tset3\nb3\tc3\tset3\nb3\tc4\tset3\n')
HMM = ''.join(f'{c}_{i}\tM{m}\t1e-10\n' for c, i, m in [
    ('c1', 1, 1), ('c1', 2, 2), ('c2', 1, 3), ('c2', 2, 4), ('c3', 1, 5),
    ('c4', 1, 6), ('c4', 2, 7), ('c4', 3, 8), ('c5', 1, 9)])


def run(tmp, n_iterations):
    tmp = Path(tmp)
    (tmp / 'profiles').mkdir(exist_ok=True)
    (tmp / 'profiles/gtdb_rel207_default_markers.tsv').write_text(
        'marker\tset\n' + ''.join(f'M{i}\tA\n' for i in range(1, 10)))
    (tmp / 'c2b.tsv').write_text(C2B)
    (tmp / 'hmm.tsv').write_text(HMM)
    out = tmp / f'out_n{n_iterations}'
    subprocess.run([sys.executable, str(SCRIPT), '-i', str(tmp / 'c2b.tsv'),
                    '--hmm', str(tmp / 'hmm.tsv'), '--profile_dir', str(tmp / 'profiles'),
                    '-m', '2', '-s', '0.5', '-n', str(n_iterations), '-o', str(out)],
                   check=True, capture_output=True)
    with open(f'{out}.scores.out') as f:
        scores = {r['bin']: r for r in csv.DictReader(f, delimiter='\t')}
    with open(f'{out}.refined.contig_to_bin.out') as f:
        refined = {}
        for r in csv.DictReader(f, delimiter='\t'):
            refined.setdefault(r['binnew'], set()).add(r['contig'])
    return scores, refined


class MagscotMerge(unittest.TestCase):
    def test_second_iteration_merges_hybrid_seed_with_original_bin(self):
        with tempfile.TemporaryDirectory(prefix='magscot-test-') as tmp:
            scores, refined = run(tmp, 2)
        self.assertEqual(int(scores['MAGScoT_1_1']['n_contigs']), 4)
        self.assertIn('MAGScoT_2_1', scores)
        self.assertEqual(int(scores['MAGScoT_2_1']['n_contigs']), 5)
        self.assertEqual(float(scores['MAGScoT_2_1']['Completeness']), 1.0)
        self.assertEqual(list(refined.values()), [{'c1', 'c2', 'c3', 'c4', 'c5'}])

    def test_single_iteration_unchanged(self):
        # Control: with one merge iteration, the best bin is the original b3.
        with tempfile.TemporaryDirectory(prefix='magscot-test-') as tmp:
            scores, refined = run(tmp, 1)
        self.assertFalse([b for b in scores if b.startswith('MAGScoT_2_')])
        self.assertEqual(list(refined.values()), [{'c1', 'c3', 'c4'}])


if __name__ == '__main__':
    unittest.main()
