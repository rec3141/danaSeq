#!/usr/bin/env python3
"""Small local regressions; compile source in a temporary directory, no big runs."""
import gzip
import itertools
import os
from pathlib import Path
import random
import resource
import signal
import subprocess
import tempfile
import unittest


SOURCE = Path(__file__).resolve().parents[1] / 'nanopore_assembly/bin/fastq_filter.cpp'


def record(name, length, base='A', quality='I'):
    return ('@' + name + '\n' + base * length + '\n+\n' + quality * length + '\n').encode()


class FastqFilterTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.build = tempfile.TemporaryDirectory(prefix='fastq-filter-tests-')
        if os.environ.get('FASTQ_FILTER_TEST_BINARY'):
            cls.binary = os.environ['FASTQ_FILTER_TEST_BINARY']
            return
        cls.binary = str(Path(cls.build.name) / 'fastq_filter')
        subprocess.run(['g++', '-O2', '-std=c++11', str(SOURCE), '-o', cls.binary,
                        '-lz', '-lpthread'], check=True)

    @classmethod
    def tearDownClass(cls):
        cls.build.cleanup()

    def run_filter(self, data, target=20000, buckets=512, extra=(), preexec=None):
        with tempfile.TemporaryDirectory(prefix='fastq-filter-spill-') as scratch:
            result = subprocess.run(
                [self.binary, '-t', str(target), '--buckets', str(buckets),
                 '--spill_dir', scratch, *extra], input=data, capture_output=True,
                preexec_fn=preexec)
            self.assertEqual(list(Path(scratch).iterdir()), [], 'spill files leaked')
            return result

    def successful(self, *args, **kwargs):
        result = self.run_filter(*args, **kwargs)
        self.assertEqual(result.returncode, 0, result.stderr.decode())
        return result.stdout

    def test_whole_buckets_have_canonical_order(self):
        reads = [record('a', 10000), record('b', 10000), record('c', 1000)]
        for budget in (20000, 50000):
            outputs = {self.successful(b''.join(order), target=budget)
                       for order in itertools.permutations(reads)}
            self.assertEqual(len(outputs), 1)

    def test_bucket_count_and_remainder(self):
        reads = [record('top', 20000), record('too_long', 12000), record('fits', 11900)]
        expected = reads[0] + reads[2]
        for buckets in (2, 64, 512, 4096):
            for order in (reads, reads[::-1]):
                self.assertEqual(self.successful(b''.join(order), target=31950,
                                                buckets=buckets), expected)

    def test_duplicate_ids_and_equal_scores(self):
        reads = [record('dup', 1000, 'A'), record('dup', 1000, 'C')]
        for extra in ((), ('--no_dedupe',)):
            for order in (reads, reads[::-1]):
                self.assertEqual(self.successful(b''.join(order), target=1000,
                                                extra=extra), reads[0])
        better = record('dup', 2000)
        for order in (reads + [better], [better] + reads):
            self.assertEqual(self.successful(b''.join(order), target=10000), better)
        self.assertEqual(self.successful(b''.join(reads), target=10000,
                                         extra=('--no_dedupe',)), b''.join(reads))

    def test_many_buckets_with_low_descriptor_limit(self):
        def limit_fds():
            resource.setrlimit(resource.RLIMIT_NOFILE, (32, 32))
        # More than 16 distinct buckets, revisited to exercise LRU reopen/append.
        reads = [record('r%d' % i, 100 + i * 31) for i in range(80)]
        data = b''.join(reads + reads[::-1])
        expected = self.successful(data, target=1000000, buckets=64)
        self.assertEqual(self.successful(data, target=1000000, buckets=4096,
                                        preexec=limit_fds), expected)

    def test_spill_failures_are_fatal_and_cleaned(self):
        def fail_writes():
            signal.signal(signal.SIGXFSZ, signal.SIG_IGN)
            resource.setrlimit(resource.RLIMIT_FSIZE, (0, 0))
        for length in (100, 100000):  # buffered close and direct write failures
            result = self.run_filter(record('a', length), target=200000,
                                     preexec=fail_writes)
            self.assertNotEqual(result.returncode, 0)
            self.assertNotIn(b'Accepted:', result.stderr)
            self.assertEqual(result.stdout, b'')

    def test_output_failures_are_fatal(self):
        for extra in (('-o', '/dev/full'),):
            self.assertNotEqual(self.run_filter(record('a', 100), extra=extra).returncode, 0)
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'output.gz'
            path.symlink_to('/dev/full')
            result = self.run_filter(record('a', 100), extra=('-o', str(path)))
            self.assertNotEqual(result.returncode, 0)
        with open('/dev/full', 'wb') as output, tempfile.TemporaryDirectory() as scratch:
            result = subprocess.run([self.binary, '-t', '1000', '--spill_dir', scratch],
                                    input=record('a', 100), stdout=output, stderr=subprocess.PIPE)
            self.assertNotEqual(result.returncode, 0)
            self.assertEqual(list(Path(scratch).iterdir()), [])

    def test_invalid_or_truncated_input_fails(self):
        for data in (b'@a\nAAA\n+\n', b'@a\nAAA\n+\nII\n',
                     gzip.compress(record('a', 100))[:-5]):
            self.assertNotEqual(self.run_filter(data).returncode, 0)
        self.assertNotEqual(self.run_filter(b'', extra=('/nonexistent/reads.fastq',)).returncode, 0)

    def test_broken_pipe_cleans_scratch(self):
        with tempfile.TemporaryDirectory() as scratch:
            read_fd, write_fd = os.pipe()
            os.close(read_fd)
            try:
                result = subprocess.run(
                    [self.binary, '-t', '20000', '--spill_dir', scratch],
                    input=record('a', 10000), stdout=write_fd, stderr=subprocess.PIPE)
            finally:
                os.close(write_fd)
            self.assertEqual(result.returncode, 1, result.stderr)
            self.assertEqual(list(Path(scratch).iterdir()), [])

    def test_permutations_and_bucket_counts_agree(self):
        rng = random.Random(7)
        reads = [record('r%03d' % i, rng.randrange(100, 3000),
                        quality=chr(rng.randrange(10, 41) + 33)) for i in range(100)]
        expected = self.successful(b''.join(reads), target=43000, buckets=2)
        for buckets in (64, 512, 4096):
            rng.shuffle(reads)
            self.assertEqual(self.successful(b''.join(reads), target=43000,
                                            buckets=buckets), expected)

    def test_empty_input_and_streaming_mode(self):
        self.assertEqual(self.successful(b''), b'')
        data = record('a', 100) + record('b', 200)
        self.assertEqual(self.successful(data, target=10000, extra=('--onepass',)), data)


if __name__ == '__main__':
    unittest.main()
