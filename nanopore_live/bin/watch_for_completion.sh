#!/usr/bin/env bash
# watch_for_completion.sh
#
# Detect end-of-run for a live nanopore_live pipeline in watch mode.
#
# Polls the MinKNOW input dir for `final_summary_*.txt` (written by MinKNOW
# when acquisition + processing have stopped). When every observed flow-cell
# run dir has a final_summary, waits for nextflow to drain (no more in-flight
# tasks in the workdir + trace.txt row count stable), then SIGTERMs the
# nextflow process targeting --outdir so it exits cleanly via -resume state.
#
# Usage:
#   watch_for_completion.sh --input DIR --outdir DIR [--workdir DIR] [--poll N]
#
# Exit codes:
#   0  detected completion + signalled nextflow (or detected but no process to signal)
#   1  bad args / missing dirs
#   2  signal failed
set -euo pipefail

INPUT=""
OUTDIR=""
WORKDIR=""
POLL=30
SELF=$(basename "$0")

usage() {
    cat <<EOF >&2
Usage: $SELF --input DIR --outdir DIR [--workdir DIR] [--poll SECONDS]

Required:
  --input DIR       MinKNOW input directory (contains FC run dirs)
  --outdir DIR      Pipeline output directory (contains pipeline_info/)

Optional:
  --workdir DIR     Nextflow -w directory (for in-flight task detection).
                    Auto-detected from \$OUTDIR/.nextflow.log if omitted.
  --poll SECONDS    Poll interval (default: 30)
EOF
}

while (( $# )); do
    case "$1" in
        -h|--help)  usage; exit 0 ;;
        --input)    INPUT="${2:-}"; shift 2 ;;
        --outdir)   OUTDIR="${2:-}"; shift 2 ;;
        --workdir)  WORKDIR="${2:-}"; shift 2 ;;
        --poll)     POLL="${2:-}"; shift 2 ;;
        *) echo "[$SELF] unknown arg: $1" >&2; usage; exit 1 ;;
    esac
done

[[ -d "$INPUT" ]]  || { echo "[$SELF] --input missing/not-a-dir: $INPUT"  >&2; exit 1; }
[[ -d "$OUTDIR" ]] || { echo "[$SELF] --outdir missing/not-a-dir: $OUTDIR" >&2; exit 1; }
[[ "$POLL" =~ ^[0-9]+$ ]] || { echo "[$SELF] --poll must be a positive integer" >&2; exit 1; }

mkdir -p "$OUTDIR/pipeline_info"
TRACE="$OUTDIR/pipeline_info/trace.txt"
MARKER="$OUTDIR/pipeline_info/final_summary_seen.txt"
LOG="$OUTDIR/pipeline_info/watch_for_completion.log"

# portable timestamp (no GNU date -Iseconds dependency)
now() { date "+%Y-%m-%dT%H:%M:%S%z"; }
log() { printf '%s [%s] %s\n' "$(now)" "$SELF" "$*" | tee -a "$LOG" >&2; }

# Locate FC run dirs anywhere under $INPUT. MinKNOW names them
# <YYYYMMDD>_<HHMM>_<POSITION>_<FCID>_<HASH>/. Match conservatively on the
# trailing _<HASH> piece and the .../final_summary_*.txt sibling pattern;
# avoids reliance on GNU find -regex.
list_fc_run_dirs() {
    # Any dir that directly contains a final_summary_*.txt OR a fastq_pass/.
    # First find candidate markers, dirname them, dedup.
    {
        find "$INPUT" -type f -name "final_summary_*.txt" 2>/dev/null \
            | while IFS= read -r f; do dirname "$f"; done
        find "$INPUT" -mindepth 1 -maxdepth 4 -type d -name "fastq_pass" 2>/dev/null \
            | while IFS= read -r d; do dirname "$d"; done
    } | sort -u
}

fc_has_final_summary() {
    # POSIX-portable existence check for final_summary_*.txt in a single dir.
    local d="$1"
    set +e
    set -- "$d"/final_summary_*.txt
    [[ -f "$1" ]]
    local rc=$?
    set -e
    return $rc
}

# Auto-detect workdir from .nextflow.log if not given.
if [[ -z "$WORKDIR" ]]; then
    if [[ -f "$OUTDIR/.nextflow.log" ]]; then
        WORKDIR=$(grep -oE 'workDir[: ]+[^ ]+' "$OUTDIR/.nextflow.log" 2>/dev/null \
            | head -1 | awk '{print $NF}' || true)
    fi
fi
if [[ -n "$WORKDIR" && ! -d "$WORKDIR" ]]; then
    log "WARN: --workdir does not exist: $WORKDIR (drain check will be trace-only)"
    WORKDIR=""
fi

count_in_flight_tasks() {
    # A nextflow task is "in flight" when .command.begin exists but no
    # .exitcode (or .command.exitcode on some versions). Returns 0 if no
    # workdir available — drain check then falls back to trace-only.
    [[ -n "$WORKDIR" && -d "$WORKDIR" ]] || { echo 0; return; }
    # Use comm on sorted-unique lists of task dirs.
    local begins exits
    begins=$(mktemp); exits=$(mktemp)
    find "$WORKDIR" -type f -name ".command.begin" 2>/dev/null \
        | while IFS= read -r f; do dirname "$f"; done | sort -u > "$begins"
    find "$WORKDIR" -type f \( -name ".exitcode" -o -name ".command.exitcode" \) 2>/dev/null \
        | while IFS= read -r f; do dirname "$f"; done | sort -u > "$exits"
    comm -23 "$begins" "$exits" | wc -l | tr -d ' '
    rm -f "$begins" "$exits"
}

count_trace_rows() {
    if [[ -f "$TRACE" ]]; then
        # Subtract 1 for header row, clamp to 0.
        local n
        n=$(wc -l < "$TRACE" | tr -d ' ')
        if [[ "$n" -gt 0 ]]; then echo $((n - 1)); else echo 0; fi
    else
        echo 0
    fi
}

# ---------- Phase 1: wait until every observed FC dir has final_summary
log "phase 1: waiting for final_summary in every FC run dir under $INPUT"
while :; do
    fc_dirs=$(list_fc_run_dirs)
    if [[ -z "$fc_dirs" ]]; then
        sleep "$POLL"; continue
    fi
    pending=0
    while IFS= read -r fc; do
        fc_has_final_summary "$fc" || { pending=$((pending + 1)); }
    done <<<"$fc_dirs"
    if [[ "$pending" -eq 0 ]]; then
        log "all $(echo "$fc_dirs" | wc -l | tr -d ' ') FC run dirs have final_summary"
        now > "$MARKER"
        break
    fi
    log "still pending: $pending FC dir(s) without final_summary"
    sleep "$POLL"
done

# ---------- Phase 2: confirm input fastq count is stable (no late stragglers)
log "phase 2: confirming input fastq count is stable"
prev_count=""
stable=0
while :; do
    count=$(find "$INPUT" -type f -name "*.fastq.gz" 2>/dev/null | wc -l | tr -d ' ')
    if [[ "$count" = "$prev_count" ]]; then
        stable=$((stable + 1))
        if [[ "$stable" -ge 2 ]]; then
            log "input fastq count stable at $count for $stable polls"
            break
        fi
    else
        log "input fastq count: $count (was ${prev_count:-?})"
        stable=0
        prev_count="$count"
    fi
    sleep "$POLL"
done

# ---------- Phase 3: wait for nextflow drain (no in-flight + trace.txt static)
log "phase 3: waiting for nextflow drain (workdir=${WORKDIR:-<unavailable>})"
prev_trace=""
drain=0
while :; do
    in_flight=$(count_in_flight_tasks)
    trace_rows=$(count_trace_rows)
    if [[ "$trace_rows" = "$prev_trace" && "$in_flight" -eq 0 ]]; then
        drain=$((drain + 1))
        if [[ "$drain" -ge 2 ]]; then
            log "drained: trace_rows=$trace_rows in_flight=$in_flight (stable for $drain polls)"
            break
        fi
    else
        log "trace_rows=$trace_rows (was ${prev_trace:-?}), in_flight=$in_flight"
        drain=0
        prev_trace="$trace_rows"
    fi
    sleep "$POLL"
done

# ---------- Phase 4: signal nextflow targeting this outdir
log "phase 4: signalling nextflow for --outdir $OUTDIR"
# Match the nextflow java process by its --outdir argument.
nf_pids=$(pgrep -f -- "--outdir $OUTDIR( |$)" || true)
if [[ -z "$nf_pids" ]]; then
    log "WARN: no nextflow process found for outdir $OUTDIR — exiting OK"
    exit 0
fi
for pid in $nf_pids; do
    if kill -TERM "$pid" 2>/dev/null; then
        log "SIGTERM sent to PID $pid"
    else
        log "ERROR: kill -TERM $pid failed"
        exit 2
    fi
done
log "watch_for_completion exiting"
