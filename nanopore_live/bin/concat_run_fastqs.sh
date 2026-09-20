#!/bin/bash
# concat_run_fastqs.sh
#
# Collapse MinKNOW's chunked per-barcode FASTQ files into one .fastq.gz per
# barcode, in place, and remove the chunks — but only after BBTools has read
# both the chunks and the result and agreed on what they contain.
#
# MinKNOW writes plain gzip (not BGZF), so a raw `cat` of the chunks is a
# valid multi-member gzip. BBTools spans the members (verified: three 100-read
# members concatenated read back as 300 reads), as do zlib, pigz, zcat and
# fastq_filter.
#
# Per barcode:
#   1. gzip -t every chunk. A chunk that is cut off, or has junk appended,
#      makes any decompressor stop at that point and report a short count —
#      and the concatenated file would report the same short count, so the
#      two would agree on a truncated total. gzip -t is the only step that
#      tells "this member ended" apart from "this file ended", so it is never
#      skipped. A cut-off chunk (the common real failure) has its complete
#      records salvaged into a repaired copy, which is used in its place and
#      recorded in the manifest.
#   2. cat the chunks, writing the result while the identical byte stream is
#      piped through reformat.sh: that parses every record in the originals
#      and gives the expected read and base counts.
#   3. reformat.sh reads the finished file back from disk.
#   4. The counts must match and be non-zero. Only then are the chunks
#      removed. Anything unproven is left exactly as it was found.
#
# reformat.sh runs under `timeout` with stdin closed where possible: on a
# malformed .gz it can throw and then hang, and a hang must fail the barcode,
# not the script.
#
# Usage: concat_run_fastqs.sh --run RUN_DIR [options]
set -uo pipefail

SELF=$(basename "$0")
RUN=""; DRY_RUN=false; KEEP=false; KEEP_TRUNC=false; FORCE=false
STAGE="${TMPDIR:-/tmp}"; MIN_CHUNKS=2; BB_TO=21600
REFORMAT="${REFORMAT:-reformat.sh}"; BBMEM="${BBMEM:--Xmx2g}"

now() { date +%Y-%m-%dT%H:%M:%S%z; }
log() { echo "$(now) [$SELF] $*" >&2; }
usage() {
    cat <<EOF >&2
Usage: $SELF --run RUN_DIR [options]
  --run DIR           MinKNOW run directory (holds fastq_pass/ and/or fastq_fail/)
  --dry-run           report what would be done, touch nothing
  --keep-originals    concatenate and verify, but do not remove the chunks
  --keep-truncated    keep the original of any chunk that had to be salvaged
  --force             re-collapse a directory that already has a manifest
  --stage DIR         local scratch for repaired chunks (default: \$TMPDIR)
  --min-chunks N      only collapse barcodes with at least N chunks (default: 2)
  --timeout S         seconds reformat.sh may spend on one stream (default: 21600)
EOF
    exit 1
}
while [[ $# -gt 0 ]]; do
    case "$1" in
        --run) RUN="${2:-}"; shift 2 ;;
        --dry-run) DRY_RUN=true; shift ;;
        --keep-originals) KEEP=true; shift ;;
        --keep-truncated) KEEP_TRUNC=true; shift ;;
        --force) FORCE=true; shift ;;
        --stage) STAGE="${2:-}"; shift 2 ;;
        --min-chunks) MIN_CHUNKS="${2:-}"; shift 2 ;;
        --timeout) BB_TO="${2:-}"; shift 2 ;;
        -h|--help) usage ;;
        *) echo "[$SELF] unknown argument: $1" >&2; usage ;;
    esac
done
[[ -n "$RUN" && -d "$RUN" ]] || { echo "[$SELF] --run must name an existing directory" >&2; exit 1; }
command -v "$REFORMAT" >/dev/null 2>&1 || { echo "[$SELF] BBTools '$REFORMAT' not on PATH" >&2; exit 1; }
GZ=$(command -v pigz || command -v gzip)
mkdir -p "$STAGE" || exit 1

TOT_OK=0; TOT_SKIP=0; TOT_FAIL=0; TOT_CHUNKS=0; TOT_TRUNC=0

# "<reads> <bases>" from reformat.sh's Input: line, or non-zero exit.
bb_parse() { awk -F'\t' '/^Input:/ {gsub(/[^0-9]/,"",$2); gsub(/[^0-9]/,"",$3); print $2, $3; f=1}
                         END {exit !f}'; }
bb_file() {  # validate a file on disk
    local out rc
    out=$(timeout "$BB_TO" "$REFORMAT" in="$1" out=null ow=t int=f $BBMEM < /dev/null 2>&1); rc=$?
    (( rc != 0 )) && { printf '%s\n' "$out" | grep -E 'Exception|Error' | head -2 >&2; return "$rc"; }
    bb_parse <<<"$out"
}

# Salvage a cut-off chunk: inflate what there is, keep only whole 4-line
# records, recompress.
salvage_chunk() {
    local src="$1" dst="$2"
    $GZ -dc "$src" 2>/dev/null \
      | awk 'NR%4==1 {if (buf != "") printf "%s", buf; buf=""} {buf = buf $0 "\n"}
             END {n=split(buf, L, "\n"); if (n-1 >= 4) printf "%s", buf}' \
      | $GZ -c > "$dst" 2>/dev/null
    [[ -s "$dst" ]] && $GZ -t "$dst" 2>/dev/null
}

collapse_dir() {
    local d="$1"
    local -a chunks
    mapfile -t chunks < <(find "$d" -maxdepth 1 -name '*.fastq.gz' -type f -printf '%f\n' | sort -V)
    local n=${#chunks[@]}
    (( n == 0 )) && return 0
    local manifest="$d/.concat_manifest.tsv"
    if [[ -f "$manifest" ]] && ! $FORCE; then
        log "SKIP $d: already collapsed"; TOT_SKIP=$((TOT_SKIP + 1)); return 0
    fi
    if (( n < MIN_CHUNKS )); then
        log "SKIP $d: $n chunk(s), nothing to collapse"; TOT_SKIP=$((TOT_SKIP + 1)); return 0
    fi

    local first="${chunks[0]}" out c
    out=$(sed -E 's/_[0-9]+\.fastq\.gz$/.fastq.gz/' <<<"$first")
    [[ "$out" == "$first" ]] && out="$(basename "$d")_all.fastq.gz"
    for c in "${chunks[@]}"; do [[ "$c" == "$out" ]] && { out="$(basename "$d")_all.fastq.gz"; break; }; done

    if $DRY_RUN; then
        log "DRY-RUN $d: would collapse $n chunks -> $out"; TOT_SKIP=$((TOT_SKIP + 1)); return 0
    fi

    # Input list: chunk paths, with repaired copies substituted after a repair.
    local -a srcs=(); for c in "${chunks[@]}"; do srcs+=("$d/$c"); done
    local -A repaired=()
    local tmp_out="$d/.partial.$out" ntrunc=0
    local in_counts out_counts rc

    # Every chunk is tested on its own first. This is not optional: a chunk
    # that is cut off or has junk appended makes the decompressor stop there
    # and report a short count, and the same short count would come back from
    # the concatenated file — the two would agree on a truncated total and the
    # originals would be deleted. gzip -t is what distinguishes "this member
    # ended" from "this file ended".
    local i bad=0
    for i in "${!srcs[@]}"; do
        $GZ -t "${srcs[$i]}" 2>/dev/null && continue
        bad=$((bad + 1))
        local rep="$STAGE/.repair.$$.$i.fastq.gz" rep_counts=""
        if salvage_chunk "${srcs[$i]}" "$rep" && rep_counts=$(bb_file "$rep") \
           && [[ "${rep_counts%% *}" != "0" ]]; then
            log "TRUNCATED ${srcs[$i]}: salvaged ${rep_counts%% *} complete reads"
            repaired["${chunks[$i]}"]=1; srcs[$i]="$rep"; ntrunc=$((ntrunc + 1))
        else
            # Nothing readable in it at all. This is the only copy, so it is
            # left for a human to look at rather than quietly dropped.
            log "FAIL $d: ${srcs[$i]} holds no recoverable reads; nothing removed"
            rm -f "$STAGE/.repair.$$."* ; TOT_FAIL=$((TOT_FAIL + 1)); return 1
        fi
    done
    (( bad > 0 )) && log "$d: $bad of $n chunks needed repair"

    rm -f "$tmp_out"
    # One read of the chunks: write the result while BBTools reads the same
    # bytes straight from the originals.
    in_counts=$(cat "${srcs[@]}" | tee "$tmp_out" \
                | timeout "$BB_TO" "$REFORMAT" in=stdin.fq.gz out=null ow=t int=f $BBMEM 2>&1 | bb_parse)
    rc=$?
    if (( rc != 0 )) || [[ -z "$in_counts" ]]; then
        log "FAIL $d: reformat.sh could not read the chunk stream; nothing removed"
        rm -f "$tmp_out" "$STAGE/.repair.$$."*; TOT_FAIL=$((TOT_FAIL + 1)); return 1
    fi
    if ! out_counts=$(bb_file "$tmp_out"); then
        log "FAIL $d: reformat.sh could not read the concatenated file; nothing removed"
        rm -f "$tmp_out" "$STAGE/.repair.$$."*; TOT_FAIL=$((TOT_FAIL + 1)); return 1
    fi
    if [[ "$in_counts" != "$out_counts" || "${in_counts%% *}" == "0" ]]; then
        log "FAIL $d: chunks [$in_counts] vs result [$out_counts]; nothing removed"
        rm -f "$tmp_out" "$STAGE/.repair.$$."*; TOT_FAIL=$((TOT_FAIL + 1)); return 1
    fi

    mv -f "$tmp_out" "$d/$out" || { log "FAIL $d: could not install $out"; TOT_FAIL=$((TOT_FAIL+1)); return 1; }
    { printf '# collapsed\t%s\tchunks %s\treads %s\tbases %s\tsalvaged %s\tvalidator %s\n' \
             "$(now)" "$n" "${in_counts%% *}" "${in_counts##* }" "$ntrunc" \
             "$("$REFORMAT" --version 2>&1 | grep -iom1 'BBTools version [0-9.]*')"
      printf '#file\tbytes\tstatus\n'
      for c in "${chunks[@]}"; do
          printf '%s\t%s\t%s\n' "$c" "$(stat -c %s "$d/$c" 2>/dev/null || echo 0)" \
                 "${repaired[$c]:+TRUNCATED(salvaged)}${repaired[$c]:-OK}"
      done; } > "$manifest"
    rm -f "$STAGE/.repair.$$."*

    if $KEEP; then
        log "OK $d: $n chunks -> $out (${in_counts%% *} reads, $ntrunc salvaged) [originals kept]"
    else
        for c in "${chunks[@]}"; do
            [[ "$c" == "$out" ]] && continue
            $KEEP_TRUNC && [[ -n "${repaired[$c]:-}" ]] && continue
            rm -f "$d/$c"
        done
        log "OK $d: $n chunks -> $out (${in_counts%% *} reads, $ntrunc salvaged), chunks removed"
    fi
    TOT_OK=$((TOT_OK + 1)); TOT_CHUNKS=$((TOT_CHUNKS + n)); TOT_TRUNC=$((TOT_TRUNC + ntrunc))
    return 0
}

log "run $RUN"
shopt -s nullglob
for top in "$RUN"/fastq_pass "$RUN"/fastq_fail; do
    [[ -d "$top" ]] || continue
    for d in "$top"/*/; do [[ -d "$d" ]] && collapse_dir "${d%/}"; done
    collapse_dir "$top"
done
shopt -u nullglob
log "done $RUN: $TOT_OK collapsed ($TOT_CHUNKS chunks, $TOT_TRUNC salvaged), $TOT_SKIP skipped, $TOT_FAIL failed"
(( TOT_FAIL == 0 ))
