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
# back to ARCHIVE_RAW_DEST / ARCHIVE_OUT_DEST environment variables. Never
# clobbers an existing destination — safe to re-run. A raw run dir that already
# exists is skipped. An output FC dir that already exists (flow cell re-used for
# a second run) is merged: files rsync'd alongside, each barcode's dana.duckdb
# kept and the new files appended to it with import_all.py (see merge_fc_out).
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
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# concat_run_fastqs.sh needs BBTools reformat.sh. run-realtime.sh calls this
# script outside the conda envs, so put the pipeline's dana-bbmap env on PATH
# when reformat.sh is not already there (otherwise barcodes stay chunked).
if ! command -v "${REFORMAT:-reformat.sh}" >/dev/null 2>&1 \
        && [[ -x "$SCRIPT_DIR/../conda-envs/dana-bbmap/bin/reformat.sh" ]]; then
    export PATH="$SCRIPT_DIR/../conda-envs/dana-bbmap/bin:$PATH"
fi

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

# Safety: never archive while a nextflow pipeline is still running against this
# outdir. The raw archive uses `rsync --remove-source-files`, which would delete
# input/intermediate files out from under live tasks (symlinked into --input),
# breaking VALIDATE_FASTQ / QC_FASTQ_FILTER etc. Match nextflow by its main.nf +
# --outdir args so this never matches archive_run.sh's own command line.
if pgrep -af -- 'main\.nf' 2>/dev/null | grep -q -- "--outdir ${OUTDIR}"; then
    echo "[$SELF] REFUSING: a nextflow pipeline is still running for --outdir $OUTDIR" >&2
    echo "[$SELF] Stop the pipeline (let it drain / SIGTERM) before archiving." >&2
    exit 3
fi

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

# Set when a merge could not be completed; the script exits non-zero at the end
# (after the integrity scan) so the caller knows the source was left in place.
ARCHIVE_FAILED=false

# import_all.py needs duckdb. run-realtime.sh calls this script outside the
# conda envs, so prefer the pipeline's dana-tools python when it is installed.
IMPORT_PY=python3
if [[ -x "$SCRIPT_DIR/../conda-envs/dana-tools/bin/python3" ]]; then
    IMPORT_PY="$SCRIPT_DIR/../conda-envs/dana-tools/bin/python3"
fi

# Merge a pipeline FC output dir into an archive FC dir that already exists.
#
# A flow cell can be washed and re-run, so a second run lands on the same
# <FC>/<barcode> paths. Per-run output files carry the run hash in their names
# (or, in batch mode, sit under per-sample names), so the trees union cleanly;
# the one file that genuinely overlaps is each barcode's dana.duckdb. Instead
# of overwriting it, rsync everything else and let import_all.py append the new
# files — it skips anything already recorded in import_log, so the existing
# rows are untouched.
#
# Safety:
#   - Preflight first: any file present on both sides with DIFFERENT content
#     (other than dana.duckdb) aborts the merge for this FC before anything is
#     touched. Identical files (e.g. the fixed tnfs.txt header) are fine.
#   - Each archived dana.duckdb is copied aside before import and restored if
#     the import fails; the run's own dana.duckdb is then kept in the source.
merge_fc_out() {
    local src="$1" dst="$2"
    log "out: $dst exists (flow cell re-used) — merging and importing into dana.duckdb"

    local conflicts=() rel
    while IFS= read -r -d '' f; do
        rel="${f#"$src"/}"
        [[ "$(basename "$rel")" == dana.duckdb ]] && continue
        if [[ -e "$dst/$rel" ]] && ! cmp -s "$f" "$dst/$rel"; then
            conflicts+=("$rel")
        fi
    done < <(find "$src" -type f -print0)
    if (( ${#conflicts[@]} )); then
        log "ERROR out: ${#conflicts[@]} file(s) differ between $src and $dst; not merging. First few:"
        printf '    %s\n' "${conflicts[@]:0:10}" >&2
        ARCHIVE_FAILED=true
        return 0
    fi

    # Files directly under the FC dir (not in a barcode dir), if any.
    run rsync -aH --remove-source-files --exclude='*/' "$src/" "$dst/"

    local bdir bn bak
    for bdir in "$src"/*/; do
        [[ -d "$bdir" ]] || continue
        bn=$(basename "$bdir")
        if [[ ! -d "$dst/$bn" ]]; then
            log "out: $bn new to this flow cell — copying as is"
            run rsync -aH --remove-source-files "$bdir" "$dst/$bn/"
            continue
        fi
        bak=""
        if [[ -f "$dst/$bn/dana.duckdb" ]]; then
            bak="$dst/$bn/dana.duckdb.pre-merge.$$"
            run cp -p "$dst/$bn/dana.duckdb" "$bak"
        fi
        log "out: $bn: rsync (except dana.duckdb), then import_all.py"
        run rsync -aH --remove-source-files --exclude=/dana.duckdb "$bdir" "$dst/$bn/"
        if run "$IMPORT_PY" "$SCRIPT_DIR/import_all.py" "$dst/$bn"; then
            [[ -n "$bak" ]] && run rm -f "$bak"
            run rm -f "$bdir/dana.duckdb"
        else
            log "ERROR out: import_all.py failed for $dst/$bn; restoring its dana.duckdb," \
                "keeping the run's own copy at $bdir/dana.duckdb"
            [[ -n "$bak" ]] && run mv -f "$bak" "$dst/$bn/dana.duckdb"
            ARCHIVE_FAILED=true
        fi
    done

    run find "$src" -depth -type d -empty -delete
    run rmdir "$src" 2>/dev/null || true
    ARCHIVED_DESTS+=("$dst")
}

archive_one() {
    local src="$1" dst="$2" kind="$3"
    if [[ -d "$dst" ]]; then
        if [[ "$kind" == out ]]; then
            merge_fc_out "${src%/}" "$dst"
        else
            log "SKIP $kind: destination already exists ($dst)"
        fi
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
# streams (see fastq_filter / docs/mapping-references.md). We verify AFTER the
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
      | xargs -0 -r -P 8 -n 1 bash -c '
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
    # Collapse MinKNOW's per-barcode chunks into one .fastq.gz per barcode
    # BEFORE the copy, while the run is still on local disk: the reading and
    # verifying all happen at SSD speed, and what crosses to the archive is
    # one file per barcode instead of hundreds. concat_run_fastqs.sh verifies
    # with BBTools before it removes anything, so a barcode it cannot prove is
    # simply copied as it stands.
    if ! $DRY_RUN; then
        if "$SCRIPT_DIR/concat_run_fastqs.sh" --run "$fc_run" --stage "${TMPDIR:-/tmp}"; then
            log "concat: collapsed per-barcode chunks in $fc_run"
        else
            log "WARNING concat: some barcodes in $fc_run stay chunked; see the log above"
        fi
    fi
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

if $ARCHIVE_FAILED; then
    log "archive_run finished WITH ERRORS — see ERROR lines above; affected sources left in place"
    exit 1
fi
log "archive_run done"
