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
# Per barcode, the data is read twice and written once:
#   1. Each chunk is read once: its bytes are appended to the result and, in
#      the same pass, piped through gzip -t. A chunk that is cut off or has
#      junk appended makes any decompressor stop at that point and report a
#      short count — and the concatenated file would report the same short
#      count, so comparing only those two would accept a truncated total.
#      gzip -t is what tells "this member ended" apart from "this file
#      ended", so it is never skipped. A cut-off chunk (the common real
#      failure) has its complete records salvaged into a repaired copy, which
#      is used in its place and recorded in the manifest.
#   2. The result must be exactly as long as the chunks that went into it.
#      Raw concatenation is byte-exact, so this is a free check that the
#      write was faithful — no extra pass for it.
#   3. reformat.sh reads the result once: it parses every record and reports
#      the read and base counts that go into the manifest. A non-zero count
#      and a clean exit are required.
# Only then are the chunks removed. Anything unproven is left as found.
#
# --work DIR does the whole thing on local disk and copies only the finished
# file back, which is what you want when the run already lives on a NAS: one
# read and one write over the wire instead of three reads. When the run is
# still on the machine that produced it, concatenate before copying it
# anywhere — archive_run.sh does that.
#
# reformat.sh runs under `timeout` with stdin closed where possible: on a
# malformed .gz it can throw and then hang, and a hang must fail the barcode,
# not the script. It is given qin=33 because nanopore FASTQ is always
# Phred+33 and BBTools' autodetection guesses Phred+64 on a read whose
# quality line is all '@', then aborts on the first out-of-range value.
#
# Usage: concat_run_fastqs.sh --run RUN_DIR [options]
set -uo pipefail

SELF=$(basename "$0")
RUN=""; DRY_RUN=false; KEEP=false; KEEP_TRUNC=false; FORCE=false
STAGE="${TMPDIR:-/tmp}"; WORKDIR=""; MIN_CHUNKS=2; BB_TO=21600
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
  --work DIR          do the work on local disk: stage each barcode's chunks
                      there, build and verify the result there, copy only the
                      result back. Use for runs that already live on a NAS.
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
        --work) WORKDIR="${2:-}"; shift 2 ;;
        --min-chunks) MIN_CHUNKS="${2:-}"; shift 2 ;;
        --timeout) BB_TO="${2:-}"; shift 2 ;;
        -h|--help) usage ;;
        *) echo "[$SELF] unknown argument: $1" >&2; usage ;;
    esac
done
[[ -n "$RUN" && -d "$RUN" ]] || { echo "[$SELF] --run must name an existing directory" >&2; exit 1; }
command -v "$REFORMAT" >/dev/null 2>&1 || { echo "[$SELF] BBTools '$REFORMAT' not on PATH" >&2; exit 1; }
# Two decompressors on purpose. pigz parallelises compression well, but it
# cannot parallelise inflating a single gzip member, so it buys nothing on
# the read side -- and on a loaded machine `pigz -dc` into a pipe parks in
# futex_wait and never returns (observed: 4s of CPU in 3 minutes of wall
# clock on a file plain gzip inflates in 2.7s). Every read, test and salvage
# path therefore uses gzip; only the recompress in salvage_chunk uses pigz.
GUNZIP=$(command -v gzip || command -v pigz)
GZ=$(command -v pigz || command -v gzip)
mkdir -p "$STAGE" || exit 1

TOT_OK=0; TOT_SKIP=0; TOT_FAIL=0; TOT_CHUNKS=0; TOT_TRUNC=0

# "<reads> <bases>" from reformat.sh's Input: line, or non-zero exit.
bb_parse() { awk -F'\t' '/^Input:/ {gsub(/[^0-9]/,"",$2); gsub(/[^0-9]/,"",$3); print $2, $3; f=1}
                         END {exit !f}'; }
bb_file() {  # validate a file on disk
    # The stream is decompressed by gzip/pigz and handed to reformat.sh on
    # stdin rather than passed as in=<file>. BBTools 39.52's multithreaded
    # gzip reader (stream.bam.BgzfInputStreamMT2) deadlocks on a loaded
    # machine -- the reader thread blocks in JobQueue.take() while its
    # decompress worker parks waiting for work, and the run sits at 0% CPU
    # until the timeout fires. Same file, same counts, seconds instead of
    # hours when gzip does the inflating.
    # pipefail so a decompression error fails the check instead of looking
    # like a short-but-valid file.
    local out rc
    out=$( set -o pipefail
           $GUNZIP -dc "$1" 2>/dev/null \
           | timeout "$BB_TO" "$REFORMAT" in=stdin.fq out=null ow=t int=f qin=33 $BBMEM 2>&1 ); rc=$?
    (( rc != 0 )) && { printf '%s\n' "$out" | grep -E 'Exception|Error' | head -2 >&2; return "$rc"; }
    bb_parse <<<"$out"
}

# Salvage a cut-off chunk: inflate what there is, keep only records that are
# structurally whole, recompress.
#
# Counting lines is not enough. A file cut mid-record usually ends with all
# four lines present but the quality line short, and emitting that record
# gets the whole salvage rejected ("Mismatch between length of bases and
# qualities"). A record is kept only if the header starts '@', the separator
# starts '+', and the quality line is exactly as long as the bases. Anything
# out of frame ends the salvage there and keeps what came before; a trailing
# partial record is dropped.
salvage_chunk() {
    local src="$1" dst="$2"
    $GUNZIP -dc "$src" 2>/dev/null \
      | awk '{ L[++k] = $0
               if (k == 4) {
                   if (substr(L[1],1,1) == "@" && substr(L[3],1,1) == "+" \
                       && length(L[2]) == length(L[4]))
                       printf "%s\n%s\n%s\n%s\n", L[1], L[2], L[3], L[4]
                   else
                       exit
                   k = 0
               } }' \
      | $GZ -c > "$dst" 2>/dev/null
    [[ -s "$dst" ]] && $GUNZIP -t "$dst" 2>/dev/null
}

collapse_dir() {
    local d="$1"
    local -a chunks
    # -not -name '.*' matters. Unlike a shell glob, find's -name '*.fastq.gz'
    # matches dotfiles, so this script's own scratch file -- .partial.<out> left
    # behind by a write-back that failed -- came back as an input chunk on the
    # next pass. It sorts first (leading dot), so its salvaged reads landed as a
    # duplicate prefix and the run reported success. Observed on 20260430
    # barcode16: 810 "chunks", 1788010 reads duplicated. Never read a dotfile.
    mapfile -t chunks < <(find "$d" -maxdepth 1 -name '*.fastq.gz' -not -name '.*' -type f -printf '%f\n' | sort -V)
    # A leftover .partial.* means a previous write-back died partway. The chunks
    # are still here, so this is recoverable, but it needs a decision rather than
    # a silent retry: refuse the barcode and say so.
    local -a stale
    mapfile -t stale < <(find "$d" -maxdepth 1 -name '.partial.*' -type f -printf '%f\n')
    if (( ${#stale[@]} > 0 )); then
        log "FAIL $d: leftover ${stale[*]} from an interrupted write-back; remove it and re-run; nothing done"
        TOT_FAIL=$((TOT_FAIL + 1)); return 1
    fi
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
    local ntrunc=0 in_counts out_counts rc wd="" tmp_out="$d/.partial.$out"
    if [[ -n "$WORKDIR" ]]; then
        # Pull the barcode onto local disk once; everything below then reads
        # and writes at SSD speed and only the finished file goes back.
        wd="$WORKDIR/$(basename "$d").$$"
        rm -rf "$wd"; mkdir -p "$wd" || { log "FAIL $d: cannot create $wd"; TOT_FAIL=$((TOT_FAIL+1)); return 1; }
        local j
        for j in "${!srcs[@]}"; do
            cp -f "${srcs[$j]}" "$wd/${chunks[$j]}" || {
                log "FAIL $d: cannot stage ${chunks[$j]}"; rm -rf "$wd"
                TOT_FAIL=$((TOT_FAIL+1)); return 1; }
            srcs[$j]="$wd/${chunks[$j]}"
        done
        tmp_out="$wd/$out"
    fi

    # Build the result. Each chunk is read once: its bytes go to the result
    # and, in the same pass, through gzip -t. Testing is not optional — a
    # chunk that is cut off or has junk appended makes any decompressor stop
    # there and report a short count, and the concatenated file would report
    # the same short count, so the two would agree on a truncated total.
    # gzip -t is what tells "this member ended" from "this file ended".
    local i want got pass bad
    for pass in 1 2; do
        rm -f "$tmp_out"; : > "$tmp_out" || {
            log "FAIL $d: cannot write $tmp_out"; rm -rf "$wd"; TOT_FAIL=$((TOT_FAIL+1)); return 1; }
        want=0; bad=0; rc=0
        for i in "${!srcs[@]}"; do
            # gzip -t quits at the first bad byte, so the tail of the stage
            # has to keep draining or tee -- and then cat -- die of SIGPIPE
            # and a corrupt chunk looks like a read error instead of one the
            # repair pass below can salvage. Drain, then report gzip's status.
            cat "${srcs[$i]}" | tee -a "$tmp_out" \
              | { $GUNZIP -t 2>/dev/null; __g=$?; cat >/dev/null 2>&1; exit "$__g"; }
            local -a st=("${PIPESTATUS[@]}")
            if (( ${st[0]:-1} != 0 )); then rc=1; break; fi
            if (( ${st[2]:-1} != 0 )); then bad=$((bad + 1)); fi
            want=$((want + $(stat -c %s "${srcs[$i]}")))
        done
        (( rc != 0 )) && { log "FAIL $d: read error while concatenating; nothing removed"
                           rm -f "$tmp_out"; rm -rf "$wd"; TOT_FAIL=$((TOT_FAIL+1)); return 1; }
        got=$(stat -c %s "$tmp_out" 2>/dev/null || echo -1)
        if (( want != got )); then
            log "FAIL $d: short write ($got of $want bytes); nothing removed"
            rm -f "$tmp_out"; rm -rf "$wd"; TOT_FAIL=$((TOT_FAIL+1)); return 1
        fi
        (( bad == 0 )) && break
        (( pass == 2 )) && { log "FAIL $d: chunks still fail gzip -t after repair; nothing removed"
                             rm -f "$tmp_out"; rm -rf "$wd"; rm -f "$STAGE/.repair.$$."*
                             TOT_FAIL=$((TOT_FAIL+1)); return 1; }

        # Repair pass: name the bad chunks and salvage their whole records.
        log "$d: $bad of $n chunks failed gzip -t; repairing"
        for i in "${!srcs[@]}"; do
            $GUNZIP -t "${srcs[$i]}" 2>/dev/null && continue
            local rep="$STAGE/.repair.$$.$i.fastq.gz" rep_counts=""
            if salvage_chunk "${srcs[$i]}" "$rep" && rep_counts=$(bb_file "$rep") \
               && [[ "${rep_counts%% *}" != "0" ]]; then
                log "TRUNCATED ${chunks[$i]}: salvaged ${rep_counts%% *} complete reads"
                repaired["${chunks[$i]}"]=1; srcs[$i]="$rep"; ntrunc=$((ntrunc + 1))
            else
                # Nothing readable in it. The archive holds the only copy, so
                # it is left for a human rather than quietly dropped.
                log "FAIL $d: ${chunks[$i]} holds no recoverable reads; nothing removed"
                rm -f "$tmp_out" "$STAGE/.repair.$$."*; rm -rf "$wd"
                TOT_FAIL=$((TOT_FAIL + 1)); return 1
            fi
        done
    done
    if ! out_counts=$(bb_file "$tmp_out"); then
        log "FAIL $d: reformat.sh could not read the concatenated file; nothing removed"
        rm -f "$tmp_out" "$STAGE/.repair.$$."*; TOT_FAIL=$((TOT_FAIL + 1)); return 1
    fi
    if [[ "${out_counts%% *}" == "0" ]]; then
        log "FAIL $d: the concatenated file holds no reads; nothing removed"
        rm -f "$tmp_out" "$STAGE/.repair.$$."*; TOT_FAIL=$((TOT_FAIL + 1)); return 1
    fi
    in_counts="$out_counts"

    if [[ -n "$wd" ]]; then
        # The archive gets one write. Compare checksums afterwards, because
        # the chunks are about to go and this becomes the only copy.
        local local_md5 remote_md5
        local_md5=$(md5sum < "$tmp_out" | cut -d" " -f1)
        cp -f "$tmp_out" "$d/.partial.$out" || {
            # Drop the half-written file. Leaving it behind is what let a later
            # pass pick it up as an input chunk (see the find above).
            log "FAIL $d: could not write the result back"; rm -f "$d/.partial.$out"; rm -rf "$wd"
            TOT_FAIL=$((TOT_FAIL+1)); return 1; }
        remote_md5=$(md5sum < "$d/.partial.$out" | cut -d" " -f1)
        if [[ "$local_md5" != "$remote_md5" ]]; then
            log "FAIL $d: the copy back does not match ($local_md5 vs $remote_md5); nothing removed"
            rm -f "$d/.partial.$out"; rm -rf "$wd"; TOT_FAIL=$((TOT_FAIL+1)); return 1
        fi
        mv -f "$d/.partial.$out" "$d/$out" || { log "FAIL $d: could not install $out"
            rm -rf "$wd"; TOT_FAIL=$((TOT_FAIL+1)); return 1; }
        rm -rf "$wd"
    else
        mv -f "$tmp_out" "$d/$out" || { log "FAIL $d: could not install $out"; TOT_FAIL=$((TOT_FAIL+1)); return 1; }
    fi
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
