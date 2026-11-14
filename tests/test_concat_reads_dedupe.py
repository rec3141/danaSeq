#!/usr/bin/env python3
"""Regression test for the CONCAT_READS dedupe step in nanopore_assembly.

The two awk programs (duplicate check and dedupe) are taken from
nanopore_assembly/modules/preprocess.nf, un-escaped the way Nextflow renders the
script block, and run on FASTQ whose headers carry SAM tags after a tab
(`samtools fastq -T` output). A read ID is the header up to the first space or
tab; header text, tabs included, must pass through unchanged.
"""
import re
import subprocess
import unittest
from pathlib import Path

NF = Path(__file__).resolve().parents[1] / 'nanopore_assembly/modules/preprocess.nf'


def concat_reads_block():
    text = NF.read_text()
    start = text.index('process CONCAT_READS')
    end = text.index('\nprocess ', start + 1)
    return text[start:end]


def render(groovy):
    # Nextflow script blocks: \\ -> \ and \$ -> $
    return groovy.replace('\\\\', '\\').replace('\\$', '$')


def dedupe_pipeline():
    """The shell commands between `zcat ..._raw.fastq.gz` and `pigz`, as rendered."""
    block = concat_reads_block()
    m = re.search(r'has_dupes=\\\$\(zcat \S+ \\\\\n\s*\| (awk .*?)\)\n', block, re.S)
    check = render(m.group(1))
    # The dedupe pipeline is the one that ends in pigz.
    body = re.findall(r'zcat \$\{meta\.id\}_raw\.fastq\.gz \\\\\n((?:\s*\|[^\n]*\n)+)', block)
    body = [b for b in body if 'pigz' in b][-1].split('| pigz')[0]
    stages = [render(s.strip().rstrip('\\').strip())
              for s in body.split('\n') if s.strip().startswith('|')]
    return check, ' '.join(stages).lstrip('| ').strip()


def sh(cmd, data):
    return subprocess.run(['bash', '-c', cmd], input=data, capture_output=True,
                          check=True).stdout


class ConcatReadsDedupe(unittest.TestCase):
    def setUp(self):
        self.check, self.dedupe = dedupe_pipeline()

    def test_tab_tagged_duplicates_are_found_and_removed(self):
        data = (b'@r1\tPU:Z:bc1\nACGT\n+\nIIII\n'
                b'@r2\tPU:Z:bc1\nGGGG\n+\nIIII\n'
                b'@r1\tPU:Z:bc2\nACGT\n+\nIIII\n')
        self.assertEqual(sh(self.check, data).strip(), b'yes')
        out = sh(self.dedupe, data)
        self.assertEqual(out, data.split(b'@r1\tPU:Z:bc2')[0])
        self.assertEqual(len(out.splitlines()) % 4, 0)

    def test_space_and_plain_headers(self):
        data = (b'@r1 runid=a\nACGT\n+\nIIII\n'
                b'@r1 runid=b\nACGT\n+\nIIII\n'
                b'@r3\nTTTT\n+\nIIII\n')
        self.assertEqual(sh(self.check, data).strip(), b'yes')
        self.assertEqual(sh(self.dedupe, data),
                         b'@r1 runid=a\nACGT\n+\nIIII\n@r3\nTTTT\n+\nIIII\n')

    def test_no_duplicates(self):
        data = b'@a\tX:Z:1\nAC\n+\nII\n@b\tX:Z:1\nAC\n+\nII\n'
        self.assertEqual(sh(self.check, data).strip(), b'')
        self.assertEqual(sh(self.dedupe, data), data)


if __name__ == '__main__':
    unittest.main()
