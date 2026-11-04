#!/usr/bin/env python3
"""Place an already-published run's ASVs on the EarthDNA globe without re-running the pipeline.

A deployed viz bundle's data/asvs.json(.gz) carries every ASV's id and sequence, which is all PLACE_GLOBE
needs: this writes them as the seqtab place_globe.py reads, runs it, and builds globe.json.gz with
build_viz.build_globe, so the result matches what a fresh run publishes. Rerun it whenever the reference
map changes.

usage: backfill_globe.py ASVS_JSON EARTHDNA_DB WORKDIR THREADS [FASTSNE_BIN]
output (in WORKDIR): globe_placements.tsv, globe_summary.json, globe.json.gz
"""
import gzip, json, os, subprocess, sys
import pandas as pd

here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, here)
import build_viz  # noqa: E402

asvs_path, db, work, threads = sys.argv[1:5]
fastsne = sys.argv[5] if len(sys.argv) > 5 else "fastsne_cpu"
os.makedirs(work, exist_ok=True)

opener = gzip.open if asvs_path.endswith(".gz") else open
with opener(asvs_path, "rt") as fh:
    asvs = json.load(fh)
seq_to_id = {a["sequence"]: a["id"] for a in asvs if a.get("sequence")}
pd.DataFrame({"sequence": list(seq_to_id)}).to_pickle(os.path.join(work, "seqtab_final.pkl"))

subprocess.run([sys.executable, os.path.join(here, "place_globe.py"), "seqtab_final.pkl", os.path.abspath(db), threads,
                fastsne], cwd=work, check=True)
build_viz.write_json(build_viz.build_globe(os.path.join(work, "globe_placements.tsv"), seq_to_id),
                     os.path.join(work, "globe.json.gz"), compress=True)
