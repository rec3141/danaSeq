"""Test selection of CAMI ablation depth-table columns."""
import os
from pathlib import Path
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
HELPER = ROOT / "cami3/ablate-lib.sh"
JGI = ("contigName\tcontigLen\ttotalAvgDepth\tA.sorted\tA.sorted-var\tB.sorted\tB.sorted-var\n"
       "ctg\t500\t2.00543\t0.130435\t0.113669\t1.875\t5.87406\n"
       "empty\t500\t0\t0\t0\t0\t0\n")

def subset(table, kind, *samples):
    with tempfile.TemporaryDirectory() as folder:
        src, dst = (Path(folder) / name for name in ("input", "output"))
        src.write_text(table)
        dst.write_text("prior output\n")
        proc = subprocess.run(["bash", "-c", 'source "$1"; shift; subset_depth "$@"',
                               "test", HELPER, dst, *samples], text=True, capture_output=True,
                              env={**os.environ, "TYPE": kind, "DEP": str(src)})
        return proc, dst.read_text()

class DepthTest(unittest.TestCase):
    def test_real_jgi_table(self):
        result, text = subset(JGI, "short", "B", "A")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(text, JGI)
        result, text = subset(JGI, "short", "B")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("ctg\t500\t1.875\t1.875\t5.87406\n", text)

    def test_long_read_subset_averages_printed_means(self):
        table = ("contigName\tcontigLen\ttotalAvgDepth\tA.sorted.bam\tA.sorted.bam-var\tB.sorted.bam\tB.sorted.bam-var\n"
                 "ctg\t500\t1.1571\t0.1714\t0.1424\t2.1429\t6.14\n"
                 "empty\t500\t0\t0\t0\t0\t0\n")
        result, text = subset(table, "long", "B", "A")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("ctg\t500\t1.15715\t0.1714\t0.1424\t2.1429\t6.14\n", text)
        self.assertIn("empty\t500\t0\t0\t0\t0\t0\n", text)

    def test_invalid_subset(self):
        for table, names in [
            (JGI, []), (JGI, ["A","A"]), (JGI, ["A","C"]),
            (JGI.replace("B.sorted", "A.sorted.bam"), ["A","C"]),
            (JGI.replace("B.sorted-var", "C.sorted-var"), ["B"]),
            (JGI.replace("0.130435", "NaN"), ["A"]),
            (JGI.replace("\t5.87406", ""), ["B"]), ("", ["A"]),
        ]:
            with self.subTest(table=table, names=names):
                result, text = subset(table, "short", *names)
                self.assertNotEqual(result.returncode, 0)
                self.assertEqual(text, "prior output\n")

if __name__ == "__main__":
    unittest.main()
