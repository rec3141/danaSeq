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

def run_subset(src, dst, kind, *samples, umask=-1):
    return subprocess.run(["bash", "-c", 'source "$1"; shift; subset_depth "$@"',
                           "test", HELPER, dst, *samples], text=True, capture_output=True,
                          umask=umask, env={**os.environ, "TYPE": kind, "DEP": str(src)})


def subset(table, kind, *samples):
    with tempfile.TemporaryDirectory() as folder:
        src, dst = (Path(folder) / name for name in ("input", "output"))
        src.write_text(table)
        dst.write_text("prior output\n")
        proc = run_subset(src, dst, kind, *samples)
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

    def test_zero_contig_jgi_table_and_unmapped_reference(self):
        header = JGI.splitlines(keepends=True)[0]
        proc, output = subset(header, "short", "A", "B")
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertEqual(output, header)
        proc, output = subset(JGI.replace("ctg\t500", "unmapped\t500"), "short", "A")
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertIn("unmapped\t500\t", output)

    def test_sample_arguments_keep_their_boundaries(self):
        proc, output = subset(JGI.replace("A.sorted", "A one.sorted"), "short", "A one")
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertIn("A one.sorted-var", output)
        for names in [("", "A"), ("A", ""), ("A B",)]:
            with self.subTest(names=names):
                proc, output = subset(JGI, "short", *names)
                self.assertNotEqual(proc.returncode, 0)
                self.assertEqual(output, "prior output\n")

    def test_bad_ids_and_reference_lengths_preserve_destination(self):
        for replacement in ["\t500", "   \t500", "ctg\tabc", "ctg\t-1", "ctg\t0",
                            "ctg\t1.5", "ctg\t2147483648"]:
            with self.subTest(replacement=replacement):
                proc, output = subset(JGI.replace("ctg\t500", replacement), "short", "A")
                self.assertNotEqual(proc.returncode, 0)
                self.assertEqual(output, "prior output\n")

    def test_directory_destination_preserves_its_contents(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            src = root / "depth.tsv"
            src.write_text(JGI)
            dst = root / "existing"
            dst.mkdir()
            (dst / "table").write_text("prior output\n")
            alias = root / "alias"
            alias.symlink_to(dst, target_is_directory=True)
            for path in [dst, alias]:
                with self.subTest(path=path):
                    proc = run_subset(src, path, "short", "A")
                    self.assertNotEqual(proc.returncode, 0)
                    self.assertEqual((dst / "table").read_text(), "prior output\n")
            self.assertTrue(alias.is_symlink())

    def test_replacement_preserves_permissions(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            src, dst = root / "input", root / "output"
            src.write_text(JGI)
            dst.write_text("prior output\n")
            for mode in [0o600, 0o640]:
                with self.subTest(mode=mode):
                    dst.chmod(mode)
                    proc = run_subset(src, dst, "short", "A", umask=0o022)
                    self.assertEqual(proc.returncode, 0, proc.stderr)
                    self.assertEqual(dst.stat().st_mode & 0o777, mode)
                    self.assertIn("ctg\t500\t0.130435", dst.read_text())

    def test_selected_overflow_preserves_destination(self):
        for table in [JGI.replace("0.130435", "1e999"),
                      JGI.replace("0.113669", "1e999"),
                      JGI.replace("0.130435", "1e308").replace("1.875", "1e308")]:
            for kind in ["short", "long"]:
                with self.subTest(table=table, kind=kind):
                    proc, output = subset(table, kind, "A", "B")
                    self.assertNotEqual(proc.returncode, 0)
                    self.assertEqual(output, "prior output\n")

if __name__ == "__main__":
    unittest.main()
