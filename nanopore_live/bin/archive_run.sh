#!/usr/bin/env bash
# archive_run.sh
#
# Archive a completed nanopore run:
#   - Raw reads: every FC run dir under --input (one with a final_summary_*.txt)
#                rsync'd to <raw-dest>/<run_name>/
#   - Pipeline outputs: every FC subdir under --outdir (matching FA*/FB*/AT*/MC-*)
#                       rsync'd to <out-dest>/<FC>/
#
# Paths are configurable so this is portable to other sites. Defaults fall
# back to ARCHIVE_RAW_DEST / ARCHIVE_OUT_DEST environment variables. Refuses
# to clobber an existing destination — safe to re-run.
#
# Usage:
#   archive_run.sh --input DIR --outdir DIR \
#       [--raw-dest DIR] [--out-dest DIR] [--dry-run]
set -euo pipefail

INPUT=""
OUTDIR=""
RAW_DEST="${ARCHIVE_RAW_DEST:-}"
OUT_DEST="${ARCHIVE_OUT_DEST:-}"
DRY_RUN=false
SELF=$(basename "$0")

usage() {
    cat <<EOF >&2
Usage: $SELF --input DIR --outdir DIR [--raw-dest DIR] [--out-dest DIR] [--dry-run]

Required:
  --input DIR     MinKNOW input directory (contains FC run dirs).
  --outdir DIR    Pipeline output directory (contains FC/BC subdirs).

Paths (required via flag or environment):
  --raw-dest DIR  Destination root for raw reads.
                    Default: \$ARCHIVE_RAW_DEST
  --out-dest DIR  Destination root for pipeline FC outputs.
                    Default: \$ARCHIVE_OUT_DEST

Optional:
  --dry-run       Print what would happen, don't move files.
EOF
}

while (( $# )); do
    case "$1" in
        -h|--help)   usage; exit 0 ;;
        --input)     INPUT="${2:-}"; shift 2 ;;
        --outdir)    OUTDIR="${2:-}"; shift 2 ;;
        --raw-dest)  RAW_DEST="${2:-}"; shift 2 ;;
        --out-dest)  OUT_DEST="${2:-}"; shift 2 ;;
        --dry-run)   DRY_RUN=true; shift ;;
        *) echo "[$SELF] unknown arg: $1" >&2; usage; exit 1 ;;
    esac
done

[[ -d "$INPUT" ]]    || { echo "[$SELF] --input missing/not-a-dir: $INPUT"   >&2; exit 1; }
[[ -d "$OUTDIR" ]]   || { echo "[$SELF] --outdir missing/not-a-dir: $OUTDIR" >&2; exit 1; }
[[ -n "$RAW_DEST" ]] || { echo "[$SELF] --raw-dest or \$ARCHIVE_RAW_DEST required" >&2; exit 1; }
[[ -n "$OUT_DEST" ]] || { echo "[$SELF] --out-dest or \$ARCHIVE_OUT_DEST required" >&2; exit 1; }

now() { date "+%Y-%m-%dT%H:%M:%S%z"; }
log() { printf '%s [%s] %s\n' "$(now)" "$SELF" "$*" >&2; }

run() {
    if $DRY_RUN; then
        log "DRY-RUN: $*"
    else
        "$@"
    fi
}

if ! $DRY_RUN; then
    mkdir -p "$RAW_DEST" "$OUT_DEST"
fi

# Destinations we actually archived this run (for the post-archive integrity scan).
ARCHIVED_DESTS=()

archive_one() {
    local src="$1" dst="$2" kind="$3"
    if [[ -d "$dst" ]]; then
        log "SKIP $kind: destination already exists ($dst)"
        return 0
    fi
    log "$kind: rsync $src/ -> $dst/"
    run rsync -aH --remove-source-files "$src/" "$dst/"
    run find "$src" -depth -type d -empty -delete
    run rmdir "$src" 2>/dev/null || true
    ARCHIVED_DESTS+=("$dst")
}

# Verify gzip/BGZF integrity of every archived FASTQ. `gzip -t` catches
# truncation / hard corruption (non-zero exit); its "trailing garbage" warning
# catches the appended-junk corruption that silently truncated cat'd BGZF
# streams (see fastq_filter / mapping_refs/PROTOCOL.md). We verify AFTER the
# move (rsync is faithful), so this flags corruption that was already present
# in the pipeline output — letting us catch a bad run before the source dirs
# are gone for good. Warns rather than failing: the data is already archived,
# and the operator needs the report, not a dead exit.
verify_archives() {
    local report="$1"; shift
    local dirs=("$@")
    [[ ${#dirs[@]} -gt 0 ]] || return 0
    if $DRY_RUN; then log "DRY-RUN: would verify FASTQ integrity under ${#dirs[@]} dest(s)"; return 0; fi
    log "verifying FASTQ integrity under ${#dirs[@]} archived dest(s)..."
    : > "$report"
    find "${dirs[@]}" -name '*.fastq.gz' -print0 2>/dev/null \
      | xargs -0 -P 8 -n 1 bash -c '
            f="$1"; msg=$(gzip -t "$f" 2>&1)
            if [[ $? -ne 0 ]]; then printf "CORRUPT(error)\t%s\n" "$f"
            elif grep -q "trailing garbage" <<<"$msg"; then printf "CORRUPT(trailing-garbage)\t%s\n" "$f"; fi
        ' _ >> "$report"
    local nbad; nbad=$(wc -l < "$report")
    if [[ "$nbad" -gt 0 ]]; then
        log "WARNING: archive integrity — $nbad corrupt FASTQ file(s) detected; see $report"
    else
        log "archive integrity OK (no corrupt FASTQ)"
        rm -f "$report"
    fi
}

# ---------- Raw archive: each FC run dir (identified by final_summary) ----------
log "scanning $INPUT for FC run dirs with final_summary"
raw_count=0
while IFS= read -r marker; do
    [[ -n "$marker" ]] || continue
    fc_run=$(dirname "$marker")
    rn=$(basename "$fc_run")
    archive_one "$fc_run" "$RAW_DEST/$rn" "raw"
    raw_count=$((raw_count + 1))
done < <(find "$INPUT" -type f -name "final_summary_*.txt" 2>/dev/null)
log "raw archive: processed $raw_count FC run dir(s)"

# ---------- Pipeline outputs: each FC subdir under $OUTDIR ----------
log "scanning $OUTDIR for FC subdirs (FA*/FB*/AT*/MC-*)"
out_count=0
for fc in "$OUTDIR"/*/; do
    [[ -d "$fc" ]] || continue
    bn=$(basename "$fc")
    case "$bn" in
        FA*|FB*|AT*|MC-*) ;;
        *) continue ;;
    esac
    archive_one "$fc" "$OUT_DEST/$bn" "out"
    out_count=$((out_count + 1))
done
log "outputs archive: processed $out_count FC subdir(s)"

# ---------- Post-archive integrity scan ----------
verify_archives "${OUTDIR}/pipeline_info/archive_integrity_report.txt" "${ARCHIVED_DESTS[@]}"

log "archive_run done"
