#!/usr/bin/env python3
"""Regression tests for single-contig bins in mag_analysis/bin/magscot.py.

Long-read assemblies have single-contig MAGs, and several binners may each
report the same contig as a bin. The input check compares bins with contigs
per binning tool, as upstream MAGScoT read_contig_to_bin() does (R/io.R), and
refinement continues while any contig row remains, so a last single-contig
bin can still be selected. Needs pandas.
"""
import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[1] / 'mag_analysis/bin/magscot.py'


def hits(contig, markers):
    return ''.join(f'{contig}_{m}\tM{m}\t1e-10\n' for m in markers)


def run_magscot(c2b, hmm):
    with tempfile.TemporaryDirectory(prefix='magscot-test-') as tmp:
        tmp = Path(tmp)
        (tmp / 'profiles').mkdir()
        (tmp / 'profiles/gtdb_rel207_default_markers.tsv').write_text(
            'marker\tset\n' + ''.join(f'M{i}\tA\n' for i in range(1, 11)))
        (tmp / 'c2b.tsv').write_text(c2b)
        (tmp / 'hmm.tsv').write_text(hmm)
        out = tmp / 'out'
        res = subprocess.run([sys.executable, str(SCRIPT), '-i', str(tmp / 'c2b.tsv'),
                              '--hmm', str(tmp / 'hmm.tsv'),
                              '--profile_dir', str(tmp / 'profiles'),
                              '--skip_merge_bins', '-o', str(out)],
                             capture_output=True, text=True)
        refined = None
        if res.returncode == 0:
            with open(f'{out}.refined.contig_to_bin.out') as f:
                refined = sorted(r['contig'] for r in csv.DictReader(f, delimiter='\t'))
        return res, refined


class MagscotSingleContigBins(unittest.TestCase):
    def test_same_contig_binned_by_two_tools(self):
        # One complete contig, reported as a bin by two binners.
        res, refined = run_magscot('m1\tc1\tmetabat\ns1\tc1\tsemibin\n',
                                   hits('c1', range(1, 11)))
        self.assertEqual(res.returncode, 0, res.stderr)
        self.assertEqual(refined, ['c1'])

    def test_swapped_columns_still_rejected(self):
        # Control: bin and contig columns swapped within one tool is an error.
        res, _ = run_magscot('c1\tm1\tmetabat\nc2\tm1\tmetabat\nc3\ts1\tsemibin\n',
                             hits('c1', range(1, 11)))
        self.assertEqual(res.returncode, 1)
        self.assertIn('More unique bins than contigs', res.stderr)

    def test_last_single_contig_bin_is_selected(self):
        # set1 A = a1+a2 (complete) wins first and removes a1 from set2's X;
        # set2's Y = d1 (7 of 10 markers) is then the only row left.
        c2b = ('A\ta1\tset1\nA\ta2\tset1\n'
               'X\ta1\tset2\nY\td1\tset2\n')
        hmm = hits('a1', range(1, 6)) + hits('a2', range(6, 11)) + hits('d1', range(1, 8))
        res, refined = run_magscot(c2b, hmm)
        self.assertEqual(res.returncode, 0, res.stderr)
        self.assertEqual(refined, ['a1', 'a2', 'd1'])


if __name__ == '__main__':
    unittest.main()
