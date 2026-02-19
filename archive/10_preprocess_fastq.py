# preprocess.py
import gzip, csv, glob, re, sys
import pysam
import pandas as pd

# seqsum = sys.argv[1]        # sequencing_summary*.txt
# tsne   = sys.argv[2]        # tsne_coords.tsv
# out    = sys.argv[3]        # e.g. events.parquet

seqsum = "sequencing_summary_FAZ82285_bb2401c5_3f629e21.txt"
tsne = "tsne_xy.csv"
out = "events.parquet"

# 1) reads -> (start_time, duration, barcode)
cols = ["read_id", "start_time", "duration", "barcode_arrangement", "filename_fastq"]
keep = []
with open(seqsum, "r") as fh:
    reader = csv.DictReader(fh, delimiter="\t")
    for row in reader:
        rid   = row["read_id"]
        st    = float(row["start_time"])
        dur   = float(row["duration"])
        bc    = row.get("barcode_arrangement") or ""
        if not bc or bc == "n/a":
            # if barcode not filled, derive from filename if your runs are split per barcode
            fn = row.get("filename_fastq","")
            m = re.search(r'barcode(\d+)', fn)
            bc = f"barcode{m.group(1)}" if m else "unbarcoded"
        keep.append((rid, st, dur, bc))
reads_df = pd.DataFrame(keep, columns=["read_id","start_time","duration","barcode"])
reads_df["end_time"] = reads_df["start_time"] + reads_df["duration"]

# 2) best contig per read from BAMs (primary alignment only)
#    If each BAM is per-barcode file, great; otherwise pull bc from reads_df merge.
pairs = []
for bam in glob.glob("*.bam"):
    sam = pysam.AlignmentFile(bam, "rb")
    for aln in sam.fetch(until_eof=True):
        if aln.is_unmapped or aln.is_secondary or aln.is_supplementary:
            continue
        # primary alignment
        rid = aln.query_name
        contig = aln.reference_name
        pairs.append((rid, contig))
    sam.close()

map_df = pd.DataFrame(pairs, columns=["read_id","contig"]).drop_duplicates("read_id")

# 3) join
events = reads_df.merge(map_df, on="read_id", how="inner")

# 4) t-SNE coords
tsne_df = pd.read_csv(tsne, sep=None, engine="python")
tsne_df.columns = [c.strip().lower() for c in tsne_df.columns]
# expect columns: contig_id/contig, tsne_x, tsne_y
for cand in ["contig","contig_id","name"]:
    if cand in tsne_df.columns:
        tsne_df = tsne_df.rename(columns={cand:"contig"})
        break

events = events.merge(tsne_df[["contig","tsne_x","tsne_y"]], on="contig", how="inner")

# Optional: clip to a time window
# events = events.query("start_time >= 0 & end_time <= 1800")  # first 30min

events.to_parquet(out, index=False)
