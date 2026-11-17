"""binning/binner_status.tsv: every requested binner ends the run in one reported state.

The report helpers are taken verbatim from mag_analysis/main.nf (the block that starts
"// Binner outcome report." and ends at "workflow {") and run under real Nextflow with
small mock processes standing in for binners, so retry-then-ignore failures, -resume,
an aborted run and the atomic publication are exercised the way the pipeline uses them.
main.nf itself is also compiled with -preview.
"""
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from test_mag_storedir_failures import resolve_nextflow  # noqa: E402

MAIN = HERE.parent / "mag_analysis" / "main.nf"
START = "// Binner outcome report."
END = "workflow {\n\n    main:\n"


def helper_block() -> str:
    text = MAIN.read_text(encoding="utf-8")
    i, j = text.index(START), text.index(END)
    return text[i:j]


# Mock steps. ok_bin and empty_bin emit an assignments table; status_bin also emits a
# status.tsv; fail_bin exits non-zero (retried once, then ignored); late_bin waits so an
# aborting step can stop the run before it finishes.
MOCK = r'''
params.outdir = 'out'
params.run_semibin = true
params.bam_dir = null
params.run_maxbin = true
params.run_lorbin = false
params.run_comebin = false
params.run_vamb = true
params.run_vamb_tax = true
params.run_sendsketch = false
params.sendsketch_address = null
params.binners_only = false
params.run_binette = true
params.checkm2_db = null
params.run_magscot = true
params.abort = false

process OK_BIN {
    output: path("ok.tsv"), emit: bins
    script: "printf 'c1\\tb1\\n' > ok.tsv"
}
process EMPTY_BIN {
    output: path("empty.tsv"), emit: bins
    script: "touch empty.tsv"
}
process STATUS_BIN {
    output:
    path("s.tsv"), emit: bins
    path("status.tsv"), emit: status
    script: "touch s.tsv; printf 'status\\nno_input_bins\\n' > status.tsv"
}
process FAIL_BIN {
    errorStrategy { task.attempt <= 1 ? 'retry' : 'ignore' }
    maxRetries 1
    output: path("never.tsv"), emit: bins
    script: "exit 3"
}
process LATE_BIN {
    output: path("late.tsv"), emit: bins
    script: "sleep 20; printf 'c\\tb\\n' > late.tsv"
}
process CONSUME {
    input: path(x)
    output: path("consumed.tsv"), emit: bins
    script: "cp ${x} consumed.tsv"
}
process ABORT {
    errorStrategy 'terminate'
    script: "sleep 2; exit 5"
}
'''

WORKFLOW = r'''
workflow {
    binnerPlan.putAll(planBinners(false, false))
    removePreviousBinnerReport()
    // dastool is launched by the plan; here it consumes FAIL_BIN's output, which never
    // arrives, so it emits nothing although nothing about it failed: no_result.
    OK_BIN();     if (!params.abort) OK_BIN.out.bins.subscribe { binnerSeen['metabat'] = it }
    EMPTY_BIN();  EMPTY_BIN.out.bins.subscribe  { binnerSeen['maxbin'] = it }
    FAIL_BIN();   FAIL_BIN.out.bins.subscribe   { binnerSeen['vamb'] = it }
    CONSUME(FAIL_BIN.out.bins); CONSUME.out.bins.subscribe { binnerSeen['dastool'] = it }
    STATUS_BIN(); STATUS_BIN.out.bins.subscribe { binnerSeen['magscot'] = it }
                  STATUS_BIN.out.status.subscribe { binnerStatus['magscot'] = it }
    if (params.abort) {
        LATE_BIN(); LATE_BIN.out.bins.subscribe { binnerSeen['metabat'] = it }
        ABORT()
    }
    workflow.onComplete = {
        try {
            writeBinnerReport(workflow.success as boolean)
        } catch (Exception e) {
            log.warn "Could not write binning/binner_status.tsv: ${e.message}"
        }
    }
}
'''


class BinnerReportTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = Path(tempfile.mkdtemp(prefix="binner-report-"))
        cls.nf, extra = resolve_nextflow(cls.tmp)
        if not cls.nf:
            raise unittest.SkipTest("nextflow not available")
        cls.env = {**os.environ, "PATH": f"{extra}:{os.environ['PATH']}",
                   "NXF_OFFLINE": "true", "NXF_HOME": str(cls.tmp / ".nxf"),
                   "NXF_ANSI_LOG": "false"}

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp, ignore_errors=True)

    def run_nf(self, d: Path, *args):
        script = d / "main.nf"
        script.write_text(MOCK + helper_block() + WORKFLOW, encoding="utf-8")
        return subprocess.run([str(self.nf), "run", str(script), "-work-dir", str(d / "work"),
                               "--outdir", str(d / "out"), *args],
                              cwd=d, env=self.env, text=True, capture_output=True, timeout=300)

    def table(self, d: Path) -> dict:
        lines = (d / "out/binning/binner_status.tsv").read_text().splitlines()
        self.assertTrue(lines[0].startswith("# run_name="), lines[0])
        self.assertEqual(lines[1], "binner\tstate\tdetail")
        return {r.split("\t")[0]: r.split("\t")[1] for r in lines[2:]}

    def run_name(self, d: Path) -> str:
        first = (d / "out/binning/binner_status.tsv").read_text().splitlines()[0]
        return first.split("\t")[0].split("=", 1)[1]

    def assert_no_partial(self, d: Path):
        self.assertEqual([p.name for p in (d / "out/binning").iterdir() if p.name.startswith(".")], [])

    def test_states_fresh_and_resumed(self):
        d = self.tmp / "normal"; d.mkdir()
        expect = {
            "semibin": "skipped",       # requested, no --bam_dir
            "metabat": "ok",            # emitted a non-empty table
            "maxbin": "completed_empty",  # emitted an empty table, no status
            "lorbin": "disabled",
            "comebin": "disabled",
            "vamb": "no_result",        # task failed twice, ignored; run still succeeds
            "vamb_tax": "skipped",      # requested without SendSketch
            "dastool": "no_result",     # its only input came from the ignored failure
            "binette": "skipped",       # requested without --checkm2_db
            "magscot": "no_input_bins",  # value from the emitted status.tsv
        }
        r = self.run_nf(d)
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        self.assertEqual(self.table(d), expect)
        self.assert_no_partial(d)
        # -resume: cached tasks emit their stored outputs, so the table is unchanged.
        r = self.run_nf(d, "-resume")
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        self.assertIn("Cached", r.stdout + r.stderr)
        self.assertEqual(self.table(d), expect)
        self.assert_no_partial(d)

    def test_binners_only_marks_consensus_not_run(self):
        d = self.tmp / "binners_only"; d.mkdir()
        r = self.run_nf(d, "--binners_only", "true")
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        t = self.table(d)
        self.assertEqual((t["dastool"], t["magscot"]), ("not_run", "not_run"))

    def test_aborted_run_reports_not_completed(self):
        d = self.tmp / "abort"; d.mkdir()
        r = self.run_nf(d, "--abort", "true")
        self.assertNotEqual(r.returncode, 0)
        t = self.table(d)
        # In this mode metabat comes only from the slow step, which ABORT stops: nothing emitted,
        # and the run did not succeed, so it is not_completed, never ok or failed.
        self.assertEqual(t["metabat"], "not_completed")
        self.assertEqual(t["semibin"], "skipped")
        self.assert_no_partial(d)

    def test_previous_report_is_replaced_and_names_its_run(self):
        d = self.tmp / "previous"; (d / "out/binning").mkdir(parents=True)
        (d / "out/binning/binner_status.tsv").write_text("# run_name=old_run\nbinner\tstate\tdetail\nmetabat\tok\t\n")
        r = self.run_nf(d)
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        self.assertNotEqual(self.run_name(d), "old_run")
        self.assertIn(self.run_name(d), r.stdout)  # the run name Nextflow printed at launch

    def test_write_failure_leaves_no_table_that_passes_for_this_run(self):
        d = self.tmp / "unwritable"; (d / "out/binning").mkdir(parents=True)
        (d / "out/binning/binner_status.tsv").write_text("# run_name=old_run\nbinner\tstate\tdetail\n")
        os.chmod(d / "out/binning", 0o555)
        self.addCleanup(os.chmod, d / "out/binning", 0o755)
        r = self.run_nf(d)
        log = r.stdout + r.stderr
        self.assertIn("Could not remove the previous", log)
        self.assertIn("Could not write binning/binner_status.tsv", log)
        # The only table there is the old one, and it names a different run.
        self.assertEqual(self.run_name(d), "old_run")
        self.assertNotIn("[old_run]", r.stdout)

    def test_main_nf_compiles(self):
        r = subprocess.run([str(self.nf), "run", str(MAIN), "-preview", "--assembly", "x.fa",
                            "--depths", "d.txt", "--outdir", str(self.tmp / "preview")],
                           cwd=self.tmp, env=self.env, text=True, capture_output=True, timeout=300)
        self.assertNotIn("Script compilation error", r.stdout + r.stderr, r.stdout + r.stderr)
        self.assertNotIn("Unexpected input", r.stdout + r.stderr, r.stdout + r.stderr)


if __name__ == "__main__":
    unittest.main()
