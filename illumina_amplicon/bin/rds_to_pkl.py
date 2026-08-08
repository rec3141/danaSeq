#!/usr/bin/env python3
"""Convert a dada2 seqtab RDS (wide matrix) to a Python pickle (long-format DataFrame).

Usage: rds_to_pkl.py <input.rds> <output.pkl>

Reads the RDS via a subprocess call to R, converts to long format.
"""
import sys, subprocess, pickle, pandas as pd, tempfile, os

if len(sys.argv) < 3:
    print("Usage: rds_to_pkl.py <input.rds> <output.pkl>", file=sys.stderr)
    sys.exit(1)

rds_path = sys.argv[1]
pkl_path = sys.argv[2]

# Use R to dump the seqtab as a TSV (long format)
with tempfile.NamedTemporaryFile(suffix='.tsv', delete=False, mode='w') as tmp:
    tmp_path = tmp.name

r_code = f"""
st <- readRDS("{rds_path}")
if (is.matrix(st)) {{
    # Wide matrix: rows=samples, cols=sequences
    samples <- rownames(st)
    sequences <- colnames(st)
    records <- NULL
    for (i in seq_len(nrow(st))) {{
        nz <- which(st[i,] > 0)
        if (length(nz) > 0) {{
            records <- rbind(records, data.frame(
                sample = samples[i],
                sequence = sequences[nz],
                count = as.integer(st[i, nz]),
                stringsAsFactors = FALSE
            ))
        }}
    }}
    write.table(records, "{tmp_path}", sep="\\t", row.names=FALSE, quote=FALSE)
    cat("[INFO] Exported", nrow(records), "non-zero entries\\n")
}} else {{
    cat("[ERROR] Input is not a matrix\\n")
    q(status=1)
}}
"""

result = subprocess.run(["Rscript", "-e", r_code], capture_output=True, text=True)
if result.returncode != 0:
    print(result.stderr, file=sys.stderr)
    sys.exit(1)
print(result.stdout, end='')

# Read TSV and save as pickle
df = pd.read_csv(tmp_path, sep='\t')
os.unlink(tmp_path)

with open(pkl_path, 'wb') as f:
    pickle.dump(df, f)

print(f"[INFO] Saved {len(df)} rows to {pkl_path}")
