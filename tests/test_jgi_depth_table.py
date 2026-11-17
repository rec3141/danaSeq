"""Test the dense CoverM to JGI converter."""
from pathlib import Path
import subprocess
import unittest

ROOT = Path(__file__).resolve().parents[1]
RESHAPE = ROOT / "nanopore_assembly/bin/jgi_depth_table.awk"
DENSE = ("Contig\tA.sorted Length\tA.sorted Mean\tA.sorted Variance"
         "\tB.sorted Length\tB.sorted Mean\tB.sorted Variance\n"
         "ctg\t500\t0.17142858\t0.14244781\t500\t2.142857\t6.1399918\n"
         "empty\t500\t0\t0\t500\t0\t0\n")
def reshape(table):
    return subprocess.run([RESHAPE], input=table, capture_output=True, text=True)


class CoverMDepthTest(unittest.TestCase):
    def test_coverm_rounding(self):
        proc = reshape(DENSE)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertEqual(proc.stdout,
                         "contigName\tcontigLen\ttotalAvgDepth\tA.sorted.bam\tA.sorted.bam-var\tB.sorted.bam\tB.sorted.bam-var\n"
                         "ctg\t500\t1.1571\t0.1714\t0.1424\t2.1429\t6.14\n"
                         "empty\t500\t0\t0\t0\t0\t0\n")

    def test_invalid_coverm(self):
        for table in [
            "", DENSE.replace("B.sorted Mean", "A.sorted Mean"),
            DENSE.replace("B.sorted", "A.sorted"), DENSE.replace("0.17142858", "NaN"),
            DENSE.replace("\t6.1399918", ""), DENSE.replace("500\t2.142857","499\t2.142857"),
        ]:
            with self.subTest(table=table):
                self.assertNotEqual(reshape(table).returncode, 0)

    def test_unmapped_is_a_contig_for_these_three_methods(self):
        result = reshape(DENSE.replace("ctg\t500", "unmapped\t500"))
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("unmapped\t500\t1.1571", result.stdout)

    def test_header_only_and_invalid_rows(self):
        header = DENSE.splitlines()[0] + "\n"
        for row in ["", "unmapped\t0", "\t500\t0\t0\t500\t0\t0",
                    "   \t500\t0\t0\t500\t0\t0", "c\t0\t0\t0\t0\t0\t0",
                    "c\t2147483648\t0\t0\t2147483648\t0\t0",
                    "c\t500\t0\t0\t500.0\t0\t0"]:
            with self.subTest(row=row):
                result = reshape(header + (row + "\n" if row else ""))
                self.assertNotEqual(result.returncode, 0)
                self.assertTrue(result.stderr.strip())

if __name__ == "__main__":
    unittest.main()
