#!/usr/bin/env bash
set -euo pipefail

# need to check whether there is something about 2000 regions per chunk that blows chunks
# ============================
# Medaka split workflow (reuse existing BAM)
# Steps: inference (chunked) -> sequence
# ============================

# ---- user config (env or edit) ----
BAM="${BAM:-/path/to/aln.sorted.bam}"            # your precomputed BAM
ASSEMBLY="${ASSEMBLY:-/path/to/assembly.fasta}"  # same reference used for BAM
OUTDIR="${OUTDIR:-medaka_split_out}"
MODEL="${MODEL:-r1041_e82_400bps_sup_v5.2.0:consensus}"    # set to your R10.4.1 SUP model
# Parallelism controls:
JOBS="${JOBS:-8}"                # how many inference chunks to run at once
TPJ="${TPJ:-2}"                  # threads per inference job (Medaka suggests <=2)
# Chunking strategy (pick ONE):
REGIONS_PER_CHUNK="${REGIONS_PER_CHUNK:-50}"     # number of contigs per chunk
# or size-based chunking (approx bp target per chunk; set to 0 to disable)
TARGET_BP_PER_CHUNK="${TARGET_BP_PER_CHUNK:-0}"
# ------------------------------------

# ---- checks ----
need() { command -v "$1" >/dev/null 2>&1 || { echo "ERROR: '$1' not found"; exit 127; }; }
for t in samtools medaka; do need "$t"; done
[ -s "$BAM" ] || { echo "ERROR: BAM not found: $BAM"; exit 1; }
[ -s "${BAM}.bai" ] || { echo "ERROR: BAM index (.bai) missing: ${BAM}.bai"; exit 1; }
[ -s "$ASSEMBLY" ] || { echo "ERROR: ASSEMBLY not found: $ASSEMBLY"; exit 1; }

mkdir -p "$OUTDIR"/{regions,results,logs}

echo "=== Medaka split mode ==="
echo "BAM:      $BAM"
echo "Assembly: $ASSEMBLY"
echo "Model:    $MODEL"
echo "Outdir:   $OUTDIR"
echo "JOBS:     $JOBS   TPJ: $TPJ"
echo

# ---- derive contigs present (mapped >0) ----
REGTAB="$OUTDIR/regions/contigs.tsv"
if [ ! -s "$REGTAB" ]; then
  # idxstats: contig length mapped unmapped
  samtools idxstats "$BAM" \
    | awk '$1!="*" && $3>0 {print $1"\t"$2"\t"$3}' \
    | sort -k2,2nr > "$REGTAB"
fi

NCONTIGS=$(wc -l < "$REGTAB" | tr -d ' ')
[ "$NCONTIGS" -gt 0 ] || { echo "ERROR: no mapped contigs found in BAM."; exit 1; }
echo "Contigs with mapped reads: $NCONTIGS"

# ---- make chunk files (regions lists) ----
CHUNKDIR="$OUTDIR/regions"
rm -f "$CHUNKDIR"/chunk_*.regions

if [ "${TARGET_BP_PER_CHUNK}" -gt 0 ]; then
  # size-based packing
  awk -v tgt="$TARGET_BP_PER_CHUNK" '
    BEGIN{sum=0; c=1}
    { if(sum==0){fn=sprintf("chunk_%04d.regions", c)}
      print $1 >> fn; sum+=$2;
      if(sum>=tgt){close(fn); c++; sum=0}
    }
    END{ }
  ' "$REGTAB"
  mv chunk_*.regions "$CHUNKDIR/" 2>/dev/null || true
else
  # fixed number of contigs per chunk
  awk -v n="'$REGIONS_PER_CHUNK'" '
    BEGIN{c=1; k=0; fn=sprintf("chunk_%04d.regions", c)}
    { print $1 >> fn; k++;
      if(k>=n){ close(fn); c++; k=0; fn=sprintf("chunk_%04d.regions", c)}
    }
  ' "$REGTAB"
  mv chunk_*.regions "$CHUNKDIR/" 2>/dev/null || true
fi

CHUNKS=( "$CHUNKDIR"/chunk_*.regions )
[ ${#CHUNKS[@]} -gt 0 ] || { echo "ERROR: no region chunks created."; exit 1; }
echo "Region chunks: ${#CHUNKS[@]} (in $CHUNKDIR)"

# ---- run medaka inference per chunk (resumable) ----
# Each chunk produces results/chunk_XXXX.hdf
run_infer() {
  local regfile="$1"
  local tag="$(basename "$regfile" .regions)"
  local hdf="$OUTDIR/results/${tag}.hdf"
  if [ -s "$hdf" ]; then
    echo "[skip] $tag exists"
    return 0
  fi
#    --model "$MODEL" \
  echo "[infer] $tag -> $(basename "$hdf")"
  medaka inference \
    --debug \
    "$BAM" "$hdf" \
    --regions $(tr '\n' ' ' < "$regfile") \
    --threads "$TPJ" >"$OUTDIR/logs/${tag}.log" 2>&1
}
export -f run_infer
export OUTDIR MODEL TPJ BAM

# Use GNU parallel if available; otherwise xargs
if command -v parallel >/dev/null 2>&1; then
  parallel -j "$JOBS" run_infer ::: "${CHUNKS[@]}"
else
  # fall back to background jobs capped at $JOBS
  i=0
  for c in "${CHUNKS[@]}"; do
    run_infer "$c" &
    i=$((i+1))
    if [ $i -ge "$JOBS" ]; then wait; i=0; fi
  done
  wait
fi

# ---- stitch results -> consensus ----
CONS="$OUTDIR/consensus.fasta"
if [ -s "$CONS" ]; then
  echo "[seq] SKIP (exists): $CONS"
else
  echo "[seq] medaka sequence -> $(basename "$CONS")"
  medaka sequence "$OUTDIR"/results/*.hdf "$ASSEMBLY" "$CONS" \
    >"$OUTDIR/logs/sequence.log" 2>&1
fi

echo "=== Done ==="
echo "Consensus: $CONS"
