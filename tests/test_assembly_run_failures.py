#!/usr/bin/env python3
"""Regression tests for failures that nanopore_assembly must not report as success.

- MAP_READS (modules/mapping.nf): a minimap2 that dies after writing part of
  its output must fail the task. The script block is rendered from the .nf
  file and run under `bash -ue` (Nextflow's default task shell) with the real
  minimap2/samtools; a wrapper makes minimap2 exit 1 after it has written its
  SAM. Skipped when minimap2 or samtools is not installed.
- run-nanopore-assembly.sh: Nextflow history files are tab-separated (start
  time, duration, run name, status, revision, session, command) and the
  duration can contain a space ("5m 13s"). The session is field 6 by tab, and
  nf_run_succeeded() accepts a non-zero container exit only when the last run
  started after launch and is recorded as OK.
"""
import os
import random
import re
import shutil
import subprocess
import tempfile
import time
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
MAPPING_NF = REPO / 'nanopore_assembly/modules/mapping.nf'
RUN_SCRIPTS = [REPO / 'nanopore_assembly/run-nanopore-assembly.sh',
               REPO / 'mag_analysis/run-mag-analysis.sh',
               REPO / 'illumina_assembly/run-illumina-assembly.sh']

SESSION_A = '5ba0f450-1a34-45f5-a5d1-f4ce73bfd1ae'
SESSION_B = 'eb536729-73e3-4524-8089-6b1ec0544b40'


def history_line(start, duration, status, session):
    return f'{start}\t{duration}\tsoggy_roentgen\t{status}\tbac9347b\t{session}\tnextflow run main.nf\n'


def render_map_reads(meta_id, fastq, assembly):
    text = MAPPING_NF.read_text()
    block = text[text.index('process MAP_READS'):]
    block = block[block.index('script:'):]
    body = block.split('"""')[1]
    subs = {'meta.id': meta_id, 'fastq': fastq, 'assembly': assembly, 'task.cpus': '1'}
    body = re.sub(r'(?<!\\)\$\{([^}]+)\}', lambda m: subs[m.group(1)], body)
    return body.replace('\\\\', '\\').replace('\\$', '$')


@unittest.skipUnless(shutil.which('minimap2') and shutil.which('samtools'),
                     'minimap2 and samtools needed')
class MapReadsPipefail(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix='map-reads-test-')
        d = Path(self.tmp.name)
        rng = random.Random(7)
        ref = ''.join(rng.choice('ACGT') for _ in range(20000))
        (d / 'ref.fa').write_text('>c1\n' + ref + '\n')
        with open(d / 'reads.fq', 'w') as f:
            for i in range(40):
                s = rng.randrange(0, 19000)
                seq = ref[s:s + 1000]
                f.write(f'@r{i}\n{seq}\n+\n{"I" * len(seq)}\n')
        self.d = d
        self.script = d / 'task.sh'
        self.script.write_text(render_map_reads('s1', 'reads.fq', 'ref.fa'))

    def tearDown(self):
        self.tmp.cleanup()

    def run_task(self, crash):
        env = dict(os.environ)
        if crash:
            stub = self.d / 'bin'
            stub.mkdir()
            (stub / 'minimap2').write_text(
                f'#!/bin/bash\n{shutil.which("minimap2")} "$@"\nexit 1\n')
            (stub / 'minimap2').chmod(0o755)
            env['PATH'] = f'{stub}:{env["PATH"]}'
        return subprocess.run(['bash', '-ue', str(self.script)], cwd=self.d, env=env,
                              capture_output=True, text=True)

    def test_minimap2_failure_fails_the_task(self):
        res = self.run_task(crash=True)
        self.assertNotEqual(res.returncode, 0, res.stderr)

    def test_successful_mapping_passes(self):
        res = self.run_task(crash=False)
        self.assertEqual(res.returncode, 0, res.stderr)
        self.assertIn('alignments', res.stderr)


class NextflowHistoryParsing(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix='nf-history-test-')
        self.hist = Path(self.tmp.name) / 'history'

    def tearDown(self):
        self.tmp.cleanup()

    def test_session_field_in_every_run_script(self):
        self.hist.write_text(history_line('2026-10-20 05:39:32', '5m 13s', 'OK', SESSION_A) +
                             history_line('2026-10-20 08:21:54', '31.9s', 'OK', SESSION_B))
        for script in RUN_SCRIPTS:
            cmds = set(re.findall(r"awk -F'\\t' '\{print \$6\}'", script.read_text()))
            self.assertTrue(cmds, f'{script.name}: no tab-split history read')
            for cmd in cmds:
                out = subprocess.run(['bash', '-c', f'{cmd} "$1" | tail -1', '_', str(self.hist)],
                                     capture_output=True, text=True, check=True).stdout.strip()
                self.assertEqual(out, SESSION_B, script.name)
            self.assertNotRegex(script.read_text(), r"awk '\{print \$6\}'", script.name)

    def nf_run_succeeded(self, since):
        text = RUN_SCRIPTS[0].read_text()
        fn = re.search(r'^nf_run_succeeded\(\) \{\n.*?^\}\n', text, re.S | re.M).group(0)
        res = subprocess.run(['bash', '-c', fn + 'nf_run_succeeded "$1" "$2"', '_',
                              str(self.hist), str(since)], capture_output=True, text=True)
        return res.returncode == 0

    def test_nf_run_succeeded(self):
        now = int(time.time())
        stamp = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(now))
        self.assertFalse(self.nf_run_succeeded(now - 5), 'missing history')
        self.hist.write_text(history_line(stamp, '5m 13s', 'OK', SESSION_A))
        self.assertTrue(self.nf_run_succeeded(now - 5))
        self.hist.write_text(history_line(stamp, '1h 2m 3s', 'ERR', SESSION_A))
        self.assertFalse(self.nf_run_succeeded(now - 5), 'failed run')
        self.hist.write_text(history_line('2026-10-20 05:39:32', '5m 13s', 'OK', SESSION_A))
        self.assertFalse(self.nf_run_succeeded(now - 5), 'earlier successful run')


if __name__ == '__main__':
    unittest.main()
