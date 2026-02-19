#!/usr/bin/env bash
set -euo pipefail

# ===========================================================
# ONT MAG polishing (resumable)
# metaFlye (polished) -> Racon (x1-2) -> Medaka
# ===========================================================

# -------- user inputs (env or edit) ----------
THREADS="${THREADS:-16}"
READS="${READS:-/path/to/reads.fastq.gz}"                # space-separated list ok
FLYE_POLISHED="${FLYE_POLISHED:-/path/to/metaflye/assembly.fasta}"
OUTDIR="${OUTDIR:-polish_out}"
RACON_ROUNDS="${RACON_ROUNDS:-1}"                        # 0,1,2
MEDAKA_MODEL="${MEDAKA_MODEL:-}"                         # e.g. r1041_e82_400bps_sup_v5.2.0
KEEP_SAM="${KEEP_SAM:-0}"                                # 1 keeps racon*.sam/paf
USE_PAF="${USE_PAF:-0}"                                  # 1 uses PAF (smaller) for racon input
FORCE="${FORCE:-0}"                                      # 1 recompute even if outputs exist
# ----------------------------------------------------------

req() { command -v "$1" >/dev/null 2>&1 || { echo "ERROR: '$1' not found in PATH"; exit 127; }; }
for tool in minimap2 samtools racon medaka_consensus; do
  command -v "$tool" >/dev/null 2>&1 || { echo "ERROR: '$tool' not found in PATH"; exit 127; }
done

[ -f "$FLYE_POLISHED" ] || { echo "ERROR: FLYE_POLISHED not found: $FLYE_POLISHED"; exit 1; }
for f in $READS; do
  [ -f "$f" ] || { echo "ERROR: READS file missing: $f"; exit 1; }
done

mkdir -p "$OUTDIR"
cp -f "$FLYE_POLISHED" "$OUTDIR/00_flye_polished.fasta"
CUR="$OUTDIR/00_flye_polished.fasta"

echo "=== Inputs ==="
echo "Threads:        $THREADS"
echo "Reads:          $READS"
echo "Flye polished:  $FLYE_POLISHED"
echo "Racon rounds:   $RACON_ROUNDS"
echo "Outdir:         $OUTDIR"
echo "Medaka model:   ${MEDAKA_MODEL:-<not set>}"
echo "KEEP_SAM:       $KEEP_SAM   USE_PAF: $USE_PAF   FORCE: $FORCE"
echo

# ---------- helpers ----------
file_ok() { [ -n "$1" ] && [ -s "$1" ]; }                 # file exists and non-empty
bam_ok()  { [ -n "$1" ] && [ -s "$1" ] && [ -s "${1}.bai" ]; }

map_to_both() {
  # usage: map_to_both draft.fa outprefix -> outprefix.sam/paf and outprefix.bam/.bai
  local draft="$1"; local prefix="$2"
  local ext="sam"; local sam_or_paf; local bam

  if [ "${USE_PAF:-0}" = "1" ]; then
    ext="paf"
  fi

  sam_or_paf="${prefix}.${ext}"
  bam="${prefix}.bam"

  if [ "${FORCE:-0}" != "1" ] && file_ok "$sam_or_paf" && bam_ok "$bam"; then
    echo "[map] SKIP (exists): $prefix.{${ext},bam}"
    return
  fi

  echo "[map] minimap2 -> $prefix.{${ext},bam}"

  if [ "$ext" = "sam" ]; then
    minimap2 -x map-ont --secondary=no -t "$THREADS" -a "$draft" $READS \
      | tee >(cat > "$sam_or_paf") \
      | samtools sort -@ "$THREADS" -o "$bam" -
  else
    minimap2 -x map-ont --secondary=no -t "$THREADS" "$draft" $READS > "$sam_or_paf"
    minimap2 -x map-ont --secondary=no -t "$THREADS" -a "$draft" $READS \
      | samtools sort -@ "$THREADS" -o "$bam" -
  fi
  samtools index -@ "$THREADS" "$bam"
}

map_to_bam() {
  # usage: map_to_bam draft.fa out.bam
  local draft="$1"; local bam_out="$2"
  if [ "${FORCE:-0}" != "1" ] && bam_ok "$bam_out"; then
    echo "[map] SKIP (exists): $bam_out"
    return
  fi
  echo "[map] minimap2 -> $bam_out"
  minimap2 -x map-ont --secondary=no -t "$THREADS" -a "$draft" $READS \
    | samtools sort -@ "$THREADS" -o "$bam_out" -
  samtools index -@ "$THREADS" "$bam_out"
}

racon_round() {
  local round="$1"
  local prefix="$OUTDIR/racon${round}"
  local racon_out="$OUTDIR/${round}.racon.fasta"
  local ext="sam"
  [ "${USE_PAF:-0}" = "1" ] && ext="paf"
  local aln="$prefix.$ext"

  if [ "${FORCE:-0}" != "1" ] && file_ok "$racon_out"; then
    echo "[racon] round $round: SKIP (exists) -> $racon_out"
    CUR="$racon_out"
    return
  fi

  map_to_both "$CUR" "$prefix"

  echo "[racon] round $round"
  racon -t "$THREADS" -m 8 -x -6 -g -8 -w 500 \
        "$READS" "$aln" "$CUR" > "$racon_out"

  if [ "${KEEP_SAM:-0}" != "1" ]; then
    rm -f "$aln"
  fi

  CUR="$racon_out"
}


# ---------- Racon loop (resumable) ----------
i=1
while [[ "$i" -le "$RACON_ROUNDS" ]]; do
  racon_round "$i"
  i=$((i+1))
done

# ---------- Medaka (resumable, reuse existing BAM if provided) ----------
echo
echo "=== Medaka ==="

# Require a model
if [ -z "${MEDAKA_MODEL:-}" ]; then
  echo "ERROR: MEDAKA_MODEL is not set."
  echo "Example: MEDAKA_MODEL=r1041_e82_400bps_sup_v5.2.0"
  exit 2
fi

MEDAKA_OUT="$OUTDIR/medaka_out"
FINAL="$MEDAKA_OUT/consensus.fasta"

# If caller passed a precomputed BAM, use it; else create one
if [ -n "${BAM_MEDAKA:-}" ] && [ -s "$BAM_MEDAKA" ]; then
  echo "[medaka] Using existing BAM_MEDAKA: $BAM_MEDAKA"
  [ -s "${BAM_MEDAKA}.bai" ] || samtools index -@ "${THREADS:-8}" "$BAM_MEDAKA"
else
  BAM_MEDAKA="$OUTDIR/medaka_input.bam"
  map_to_bam "$CUR" "$BAM_MEDAKA"
fi

# quick sanity: does BAM have any alignments?
#samtools view -c "$BAM_MEDAKA" | awk '{if($1==0){exit 1}}' \
#  || { echo "[medaka] ERROR: BAM has zero alignments to $CUR"; exit 1; }

# Skip Medaka if already done
if [ "${FORCE:-0}" != "1" ] && [ -s "$FINAL" ]; then
  echo "[medaka] SKIP (final exists): $FINAL"
  echo "=== Done ==="
  echo "Flye-polished input:      $FLYE_POLISHED"
  echo "After Racon rounds:       $CUR"
  echo "Medaka consensus (final): $FINAL"
  echo "All outputs under:        $OUTDIR"
  exit 0
fi

mkdir -p "$MEDAKA_OUT"
echo "[medaka_consensus] model=${MEDAKA_MODEL}"

set +e
medaka_consensus \
  -i $READS \
  -d "$CUR" \
  -o "$MEDAKA_OUT" \
  -t "$THREADS" \
  -m "$MEDAKA_MODEL"
rc=$?
set -e

if [ $rc -ne 0 ]; then
  echo "[medaka] ERROR (exit $rc). Check model name and inputs."
  exit $rc
fi

echo
echo "=== Done ==="
echo "Flye-polished input:      $FLYE_POLISHED"
echo "After Racon rounds:       $CUR"
echo "Medaka consensus (final): $FINAL"
echo "All outputs under:        $OUTDIR"


