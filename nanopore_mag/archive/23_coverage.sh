#!/usr/bin/env bash
set -euo pipefail

# ===========================================================
# ONT MAG post-polish helper
# 1) Build contig-depth tables for binning (CoverM; optional jgi)
# 2) (Optional) Run CheckM2 and GUNC on bins if provided
# ===========================================================

source activate checkm2

# ---------- user-configurable (env or edit) ----------------
THREADS="${THREADS:-32}"

# Polished contigs (from your first script; e.g., medaka_out/consensus.fasta)
ASSEMBLY="${ASSEMBLY:-/path/to/polished/consensus.fasta}"

# One or more ONT read files (space-separated). You can also provide BAMs instead (see below).
READS="${READS:-/path/to/reads1.fastq.gz /path/to/reads2.fastq.gz}"

# Output dir
OUTDIR="${OUTDIR:-coverage-out}"

# CoverM mapping & long-read friendly thresholds
COVERM_IDENTITY="${COVERM_IDENTITY:-0.85}"       # min read percent identity (ONT SUP non-duplex)
COVERM_ALIGNED_FRAC="${COVERM_ALIGNED_FRAC:-0.50}" # min fraction of read aligned to count
COVERM_METHODS="${COVERM_METHODS:-mean covered_bases variance tpm rpkm}"

# Optional: use precomputed BAMs instead of raw reads (space-separated)
# BAMs must be primary-alignments only, sorted & indexed
BAMS="${BAMS:-}"

# Optional: also emit a jgi-style depth table using MetaBAT’s jgi_summarize (tuned for ONT)
RUN_JGI="${RUN_JGI:-0}"
JGI_MIN_MAPQ="${JGI_MIN_MAPQ:-5}"
JGI_MIN_PID="${JGI_MIN_PID:-80}"

# Optional QC on bins (set BINS_DIR to directory of MAG FASTAs, e.g. *.fa)
BINS_DIR="${BINS_DIR:-}"         # e.g., /work/bins_dastool/nonredundant_bins
BINS_EXT="${BINS_EXT:-fa}"       # fa|fna|fasta

# Databases (required only if QC is enabled)
CHECKM2_DB="${CHECKM2_DB:-/path/to/CheckM2_database}"     # e.g., checkm2_db_2024_05
GUNC_DB="${GUNC_DB:-/path/to/gunc_db}"                    # e.g., gunc_db_progenomes2.1
# -----------------------------------------------------------

# ---------- tool checks ----------
need() { command -v "$1" >/dev/null 2>&1 || { echo "ERROR: '$1' not found"; exit 127; }; }
for t in coverm samtools minimap2; do need "$t"; done
[[ -f "$ASSEMBLY" ]] || { echo "ERROR: ASSEMBLY not found: $ASSEMBLY"; exit 1; }

mkdir -p "$OUTDIR"
LOG="$OUTDIR/run.log"

echo "=== Post-polish depth & QC ===" | tee "$LOG"
echo "Assembly: $ASSEMBLY" | tee -a "$LOG"
echo "Threads:  $THREADS" | tee -a "$LOG"
echo "Reads:    ${BAMS:+<using BAMs>} ${BAMS:+"$BAMS"} ${READS:+"$READS"}" | tee -a "$LOG"
echo "Outdir:   $OUTDIR" | tee -a "$LOG"
echo | tee -a "$LOG"

# ---------- 1) Contig coverage with CoverM ----------
COVERM_OUT="$OUTDIR/coverm_contig.tsv"
echo "[CoverM] computing per-contig depths -> $COVERM_OUT" | tee -a "$LOG"

if [[ -n "$BAMS" ]]; then
  # Use precomputed BAMs
  coverm contig \
    --bam-files $BAMS \
    -r "$ASSEMBLY" \
    --methods $COVERM_METHODS \
    --threads "$THREADS" \
    --min-read-percent-identity "$COVERM_IDENTITY" \
    --min-read-aligned-percent "$COVERM_ALIGNED_FRAC" \
    --output-file "$COVERM_OUT"
else
  # Map internally with minimap2 using ONT preset
  coverm contig \
    --reads $READS \
    -r "$ASSEMBLY" \
    --mapper minimap2 \
    --minimap2-args "-x map-ont --secondary=no" \
    --methods $COVERM_METHODS \
    --threads "$THREADS" \
    --min-read-percent-identity "$COVERM_IDENTITY" \
    --min-read-aligned-percent "$COVERM_ALIGNED_FRAC" \
    --output-file "$COVERM_OUT"
fi

# Helpful derivatives for binners
# A slim matrix many tools like (contig + per-sample mean coverage only):
awk 'BEGIN{FS=OFS="\t"} NR==1{print $1,$(NF-3)} NR>1{print $1,$(NF-3)}' "$COVERM_OUT" > "$OUTDIR/contig_mean_coverage.tsv"
# NOTE: for multi-sample, prefer using CoverM’s wide format directly; most binners accept the full table.

# ---------- (optional) 1b) jgi-style depth (ONT-tuned) ----------
if [[ "$RUN_JGI" == "1" ]]; then
  need jgi_summarize_bam_contig_depths
  echo "[jgi] generating ONT-tuned jgi depth (MAPQ>=$JGI_MIN_MAPQ, PID>=$JGI_MIN_PID%)" | tee -a "$LOG"

  # Build/collect BAMs against the same assembly (if not provided)
  if [[ -z "$BAMS" ]]; then
    BAMDIR="$OUTDIR/jgi_bams"
    mkdir -p "$BAMDIR"
    i=1
    for R in $READS; do
      BAM="$BAMDIR/sample${i}.bam"
      echo "  mapping $R -> $BAM" | tee -a "$LOG"
      minimap2 -t "$THREADS" -x map-ont --secondary=no "$ASSEMBLY" "$R" \
        | samtools sort -@ "$THREADS" -o "$BAM" -
      samtools index -@ "$THREADS" "$BAM"
      i=$((i+1))
    done
    BAMS="$(printf '%s ' $BAMDIR/*.bam)"
  fi

  # jgi summarize
  JGI_OUT="$OUTDIR/jgi_depth.txt"
  jgi_summarize_bam_contig_depths \
    --outputDepth "$JGI_OUT" \
    --minMapQual "$JGI_MIN_MAPQ" \
    --percentIdentity "$JGI_MIN_PID" \
    $BAMS

  echo "  jgi depth -> $JGI_OUT" | tee -a "$LOG"
fi

# ---------- 2) (optional) QC on bins ----------
if [[ -n "$BINS_DIR" ]]; then
  echo | tee -a "$LOG"
  echo "[QC] Bins directory provided: $BINS_DIR" | tee -a "$LOG"
  [[ -d "$BINS_DIR" ]] || { echo "ERROR: BINS_DIR is not a directory: $BINS_DIR"; exit 1; }

  # CheckM2
  if [[ -n "$CHECKM2_DB" ]]; then
    need checkm2
    echo "[CheckM2] running genome quality on bins (*.$BINS_EXT)" | tee -a "$LOG"
    CHECKM2_OUT="$OUTDIR/checkm2"
    mkdir -p "$CHECKM2_OUT"
    checkm2 predict \
      --threads "$THREADS" \
      --database_path "$CHECKM2_DB" \
      --extension "$BINS_EXT" \
      --input "$BINS_DIR" \
      --output_directory "$CHECKM2_OUT" \
      --remove_intermediates
    echo "  CheckM2 summary -> $CHECKM2_OUT/quality_report.tsv" | tee -a "$LOG"
  else
    echo "[CheckM2] skipped (set CHECKM2_DB to enable)" | tee -a "$LOG"
  fi

  # GUNC
  if [[ -n "$GUNC_DB" ]]; then
    need gunc
    echo "[GUNC] assessing chimerism/contamination" | tee -a "$LOG"
    GUNC_OUT="$OUTDIR/gunc"
    mkdir -p "$GUNC_OUT"
    gunc run \
      --input_dir "$BINS_DIR" \
      --input_pattern "*.$BINS_EXT" \
      --threads "$THREADS" \
      --db_file "$GUNC_DB" \
      --output_dir "$GUNC_OUT" \
      --detailed_output
    echo "  GUNC outputs in -> $GUNC_OUT" | tee -a "$LOG"
  else
    echo "[GUNC] skipped (set GUNC_DB to enable)" | tee -a "$LOG"
  fi
else
  echo "[QC] No BINS_DIR provided; skipping CheckM2/GUNC (run this after binning)." | tee -a "$LOG"
fi

echo
echo "=== Finished ==="
echo "CoverM table:           $COVERM_OUT"
[[ "$RUN_JGI" == "1" ]] && echo "jgi depth:              $JGI_OUT"
echo "Slim mean-coverage:    $OUTDIR/contig_mean_coverage.tsv"
[[ -n "$BINS_DIR" ]] && echo "QC under:               $OUTDIR/{checkm2,gunc}"
