#!/usr/bin/env python3
"""Place a run's ASVs on the EarthDNA globe of full-length SSU references (https://github.com/rec3141/earthdna).

1. Each ASV (both strands) is mapped with vsearch to the database's panel of aligned references; the alignment
   columns it covers say which SSU region the amplicon is, so no primer list is needed and a run holding two
   assays (16S V4 + 18S V9) sorts itself. At least 80% of the ASV must align at 60% identity or more, which keeps
   out ITS and 28S LSU amplicons.
2. Each ASV goes to the prebuilt region whose columns it overlaps best (Jaccard >= 0.6); the others stay
   unplaced and are counted.
3. Per region, the ASVs (forward strand) are projected with fastsne onto that region's placement model: references
   cut to the same columns, carrying their full-length map positions.

usage: place_globe.py seqtab_final.pkl EARTHDNA_DB THREADS [FASTSNE_BIN]
output: globe_placements.tsv  sequence, lon, lat, nn_dist, region, identity
        globe_summary.json    per-region counts, unplaced count, database manifest
"""
import collections, json, os, pickle, struct, subprocess, sys
import numpy as np

seqtab_path, db, threads = sys.argv[1], sys.argv[2], sys.argv[3]
fastsne = sys.argv[4] if len(sys.argv) > 4 else "fastsne_cpu"
manifest = json.load(open(os.path.join(db, "manifest.json")))
regions = [tuple(int(x) for x in r.split("-")) for r in manifest["regions"]]
K = manifest.get("k", 5)

with open(seqtab_path, "rb") as fh:
    seqtab = pickle.load(fh)
seqs = list(seqtab["sequence"].unique())
qid = {f"q{i}": s for i, s in enumerate(seqs)}
with open("asvs.fa", "w") as f:
    for q, s in qid.items(): f.write(f">{q}\n{s}\n")
print(f"[INFO] {len(seqs)} ASVs", flush=True)

# 1. alignment columns covered by each ASV
names = [l[1:].split()[0] for l in open(os.path.join(db, "panel.fa")) if l[0] == ">"]
row = {n: i for i, n in enumerate(names)}
cols, offs = np.load(os.path.join(db, "panel_cols.npy")), np.load(os.path.join(db, "panel_offs.npy"))
subprocess.run(["vsearch", "--usearch_global", "asvs.fa", "--db", os.path.join(db, "panel.fa"), "--id", "0.6",
                "--query_cov", "0.8", "--strand", "both", "--maxaccepts", "4", "--maxrejects", "64", "--top_hits_only",
                "--threads", threads, "--userout", "panel_hits.tsv",
                "--userfields", "query+target+id+qstrand+tilo+tihi", "--quiet"], check=True)
best = {}
for line in open("panel_hits.tsv"):
    q, t, ident, st, lo, hi = line.rstrip("\n").split("\t")
    if q not in best or float(ident) > best[q][0]:
        c = cols[offs[row[t]]:offs[row[t] + 1]]
        best[q] = (float(ident), st, int(c[int(lo) - 1]), int(c[int(hi) - 1]))

# 2. region per ASV
def jacc(x, y): return max(0, min(x[1], y[1]) - max(x[0], y[0])) / (max(x[1], y[1]) - min(x[0], y[0]))
assign = collections.defaultdict(list)
for q, (ident, st, c0, c1) in best.items():
    j = [jacc((c0, c1), r) for r in regions]
    if j and max(j) >= 0.6: assign[int(np.argmax(j))].append(q)

# 3. project per region
lut = np.full(256, 255, np.uint8)
for i, ch in enumerate(b"ACGT"): lut[ch] = i; lut[ch + 32] = i
lut[ord("U")] = lut[ord("u")] = 3
pw = (4 ** np.arange(K - 1, -1, -1)).astype(np.int64)
def counts(s):
    c = lut[np.frombuffer(s.encode(), np.uint8)]
    bad = c > 3
    c = np.where(bad, 0, c).astype(np.int64)
    if len(c) < K: return np.zeros(4 ** K, np.float32)
    w = np.lib.stride_tricks.sliding_window_view(c, K) @ pw
    ok = ~np.lib.stride_tricks.sliding_window_view(bad, K).any(1)
    return np.bincount(w[ok], minlength=4 ** K).astype(np.float32)
rc = str.maketrans("ACGTURYKMBDHVN", "TGCAAYRMKVHDBN")

placed = {}
for k, qs in sorted(assign.items()):
    c0, c1 = regions[k]
    with open("region.fsb", "wb") as f:
        f.write(b"FSB1" + struct.pack("<iq", 4 ** K, len(qs)))
        for q in qs:
            s = qid[q].upper()
            if best[q][1] == "-": s = s.translate(rc)[::-1]
            f.write(counts(s).tobytes())
        f.write("\n".join(qs).encode())
    subprocess.run([fastsne, "--project", os.path.join(db, f"regions/{c0}-{c1}.model"), "-i", "region.fsb",
                    "-o", "region_placed.tsv", "--threads", threads], check=True)
    for line in list(open("region_placed.tsv"))[1:]:
        q, lo, la, nd = line.rstrip("\n").split("\t")
        placed[q] = (lo, la, nd, f"{c0}-{c1}")
    print(f"[INFO] region {c0}-{c1}: {len(qs)} ASVs", flush=True)

with open("globe_placements.tsv", "w") as f:
    f.write("sequence\tlon\tlat\tnn_dist\tregion\tidentity\n")
    for q, (lo, la, nd, reg) in placed.items():
        f.write(f"{qid[q]}\t{lo}\t{la}\t{nd}\t{reg}\t{best[q][0]}\n")
summary = dict(asvs=len(seqs), placed=len(placed), unmapped=len(seqs) - len(best),
               no_region=len(best) - sum(len(v) for v in assign.values()),
               regions={f"{regions[k][0]}-{regions[k][1]}": len(v) for k, v in sorted(assign.items())},
               reference=manifest)
json.dump(summary, open("globe_summary.json", "w"), indent=1)
print(f"[INFO] placed {len(placed)} of {len(seqs)} ASVs ({summary['unmapped']} matched no SSU reference, "
      f"{summary['no_region']} covered no prebuilt region)", flush=True)
