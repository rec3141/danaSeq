#!/usr/bin/env python3
"""Regression test for nanopore_live/bin/watch_for_completion.sh, phase 4.

The watcher must SIGTERM only the nextflow process for its --outdir. It must
not signal itself or the launcher that started it, even though both carry the
same `--outdir DIR` on their command lines (run-realtime.sh starts the watcher
with `--outdir`, and is itself started with `--outdir`).

The test builds that process tree with stand-ins:

    driver.sh --input IN --outdir OUT        (stands in for run-realtime.sh)
      |- watch_for_completion.sh --input IN --outdir OUT --poll 0   (background)
      `- nextflow run main.nf --input IN --outdir OUT              (foreground)

The fake nextflow records SIGTERM and exits; the driver then reaps the watcher
and records its exit status, as run-realtime.sh does before archiving.

Only processes whose command line contains this test's private temp path can
match, so nothing else on the machine is signalled. Linux/macOS, needs pgrep.
"""
import os
import signal
import subprocess
import tempfile
import textwrap
import time
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
WATCHER = REPO / "nanopore_live" / "bin" / "watch_for_completion.sh"

FAKE_NEXTFLOW = textwrap.dedent("""\
    #!/usr/bin/env bash
    # Stand-in for nextflow in watch mode: runs until signalled.
    out=""
    while (( $# )); do [[ "$1" == --outdir ]] && out="$2"; shift; done
    trap 'echo TERM > "$out/nextflow_got_term"; exit 143' TERM
    while :; do sleep 0.1; done
    """)

DRIVER = textwrap.dedent("""\
    #!/usr/bin/env bash
    # Stand-in for run-realtime.sh: watcher in the background, nextflow in the
    # foreground, then reap the watcher.
    in="$2"; out="$4"
    "$WATCHER" --input "$in" --outdir "$out" --poll 0 &
    watch_pid=$!
    "$FAKEBIN/nextflow" run main.nf --input "$in" --outdir "$out"
    nf_rc=$?
    if wait "$watch_pid"; then watch_rc=0; else watch_rc=$?; fi
    echo "$watch_rc" > "$out/watcher_rc"
    echo "$nf_rc" > "$out/nextflow_rc"
    echo done > "$out/driver_finished"
    """)


class WatcherSignalsOnlyNextflow(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix="wfc-test-")
        root = Path(self.tmp.name)
        self.inp = root / "input"
        self.out = root / "output"
        fc = self.inp / "20261114_0000_X1_FAK00000_deadbeef"
        fc.mkdir(parents=True)
        (fc / "final_summary_FAK00000_deadbeef.txt").write_text("done\n")
        self.out.mkdir()
        self.fakebin = root / "fakebin"
        self.fakebin.mkdir()
        nf = self.fakebin / "nextflow"
        nf.write_text(FAKE_NEXTFLOW)
        nf.chmod(0o755)
        self.driver = root / "driver.sh"
        self.driver.write_text(DRIVER)
        self.driver.chmod(0o755)

    def tearDown(self):
        # Stop anything left over (the unfixed watcher leaves nextflow running).
        try:
            os.killpg(self.proc.pid, signal.SIGKILL)
        except (AttributeError, ProcessLookupError, PermissionError):
            pass
        subprocess.run(["pkill", "-KILL", "-f", str(self.out)],
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        self.tmp.cleanup()

    def test_phase4_signals_nextflow_not_watcher_or_launcher(self):
        env = dict(os.environ, WATCHER=str(WATCHER), FAKEBIN=str(self.fakebin))
        self.proc = subprocess.Popen(
            ["bash", str(self.driver), "--input", str(self.inp),
             "--outdir", str(self.out)],
            env=env, start_new_session=True,
            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        deadline = time.monotonic() + 60
        while time.monotonic() < deadline and self.proc.poll() is None:
            time.sleep(0.2)
        log = self.out / "pipeline_info" / "watch_for_completion.log"
        log_text = log.read_text() if log.exists() else "<no log>"

        self.assertTrue((self.out / "nextflow_got_term").exists(),
                        "nextflow was never signalled\n" + log_text)
        self.assertTrue((self.out / "driver_finished").exists(),
                        "the launcher did not survive phase 4\n" + log_text)
        self.assertEqual((self.out / "watcher_rc").read_text().strip(), "0",
                         "the watcher did not exit 0\n" + log_text)
        self.assertTrue(
            (self.out / "pipeline_info" / "final_summary_seen.txt").exists())


if __name__ == "__main__":
    unittest.main()
