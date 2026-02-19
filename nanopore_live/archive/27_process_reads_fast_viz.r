#!/usr/bin/env bash
# quick-nano-barcode.fast.sh — fast, simple QC + annotation for Nanopore barcoded reads
# Goals vs original:
#  - 10x wall‑clock speed via parallelism, fewer passes over disk, and streaming I/O
#  - safer flags, clear args, and resumable idempotent outputs
#  - deterministic per-file ordering when needed (md5 over path)
#
# Requirements: GNU coreutils, awk, grep, sort, xargs, find, gzip; BBMap (bbduk.sh, reformat.sh),
#   Filtlong, Kraken2 (+ DB), Prokka, Rscript; optional: GNU parallel, pigz.
#
# Usage examples:
#   ./quick-nano-barcode.fast.sh -i /data/project_QEI2025 -o out_dana -S -K -P -T -t 8
#   ./quick-nano-barcode.fast.sh --help

set -euo pipefail
IFS=$'\n\t'

# -------------------------------
# Defaults & Env
# -------------------------------
: "${JAVA_TOOL_OPTIONS:= -Dhttps.protocols=TLSv1.2}"
export JAVA_TOOL_OPTIONS

THREADS=${THREADS:-$(nproc)}
BBMAP=${BBMAP:-/work/apps/bbmap}
DANADIR=${DANADIR:-/work/apps/dana}
PROKKA_BIN=${PROKKA_BIN:-/work/apps/prokka/bin/prokka}
FILTLONG=${FILTLONG:-/work/apps/Filtlong/bin/filtlong}
KRAKEN2=${KRAKEN2:-/usr/bin/kraken2}
KRAKEN_DB=${KRAKEN_DB:-/data/scratch/refdbs/krakendb/pluspfp_08gb}
APPS=${APPS:-/work/apps}

# run toggles
RUN_SKETCH=0
RUN_KRAKEN=0
RUN_TETRA=0
RUN_PROKKA=0
PROKKA_THREADS=${PROKKA_THREADS:-0}  # 0 = auto

MAX_DURATION=${MAX_DURATION:-3600}  # 1h safety cap per batch
MIN_READLEN=${MIN_READLEN:-1500}
KEEP_PCT=${KEEP_PCT:-80}

# -------------------------------
# CLI
# -------------------------------
print_help(){ cat <<EOF
  Fast QC + annotation for barcoded Nanopore FASTQ(.gz) in directories containing fastq_pass/.
  
  Required:
    -i, --input DIR       Project dir containing symlinks to run/barcode outputs
  Optional:
    -o, --output DIR      Output root (default: out_<basename(input)>_<rand>)
  -t, --threads N       Global threads (default: nproc)
  -S, --sketch          Run BBMap sendsketch + Dana DB update
  -K, --kraken          Run Kraken2 + Dana DB updates
  -P, --prokka          Run Prokka + Dana DB update (metagenome mode)
  -T, --tetra           Run tetramer ESOM pipeline + Dana DB
  --max-duration SEC    Stop after SEC for this batch (default: ${MAX_DURATION})
  --min-readlen N       Min read length for keep (default: ${MIN_READLEN})
  --keep-pct P          Filtlong keep percent (default: ${KEEP_PCT})
  --deterministic       Process files in md5(path) order (default: filesystem order)
  -h, --help            This help
  
  Env overrides:
    BBMAP, DANADIR, PROKKA_BIN, FILTLONG, KRAKEN2, KRAKEN_DB, THREADS, PROKKA_THREADS
  EOF
}

INPUT=""
OUTPUT=""
DETERMINISTIC=0

while (( $# )); do
         case "$1" in
         -i|--input)   INPUT="$2"; shift 2 ;;
       -o|--output)  OUTPUT="$2"; shift 2 ;;
-t|--threads) THREADS="$2"; shift 2 ;;
-S|--sketch)  RUN_SKETCH=1; shift ;;
-K|--kraken)  RUN_KRAKEN=1; shift ;;
-P|--prokka)  RUN_PROKKA=1; shift ;;
-T|--tetra)   RUN_TETRA=1; shift ;;
--max-duration) MAX_DURATION="$2"; shift 2 ;;
--min-readlen)  MIN_READLEN="$2"; shift 2 ;;
--keep-pct)     KEEP_PCT="$2"; shift 2 ;;
--deterministic) DETERMINISTIC=1; shift ;;
-h|--help) print_help; exit 0 ;;
*) echo "Unknown arg: $1" >&2; print_help; exit 2 ;;
esac
done

if [[ -z "${INPUT}" ]]; then echo "[ERR] --input is required" >&2; exit 2; fi
if [[ ! -d "${INPUT}" ]]; then echo "[ERR] input dir missing: ${INPUT}" >&2; exit 2; fi

if [[ -z "${OUTPUT}" ]]; then OUTPUT=$(mktemp -u "out_$(basename "${INPUT}")_XXXXXX"); fi
mkdir -p "${OUTPUT}"

# temp cache for validated FASTQs
CACHE_FASTQ=/data/.fastq_pass
mkdir -p "${CACHE_FASTQ}"

# -------------------------------
# Collect FASTQ files
# -------------------------------
mapfile -d '' FILES < <( \
                         find "${INPUT}" -follow -type f -name '*fastq.gz' -size +10k \
                         -regex '.*/fastq_pass/.*' -path '*barcode*' \
                         ! -path '*/data/minknow/*' ! -path "*${OUTPUT}*" ! -name '*.tmp.*' -print0 )

NUMFILES=${#FILES[@]}
  echo "[INFO] Found ${NUMFILES} FASTQ files"
  
  # deterministic ordering (optional)
  if (( DETERMINISTIC )); then
  # shellcheck disable=SC2016
  mapfile -t FILES_SORTED < <(printf '%s\n' "${FILES[@]}" | awk '{print | "md5sum"}' | sort | awk '{print $2}')
  FILES=("${FILES_SORTED[@]}")
  fi
  
  # Build the (flowcell, barcode) dir structure up front
  printf '%s\n' "${FILES[@]}" \
  | awk -F'/' '{n=split($NF,a,"_"); if(n>=3) print a[1]"\t"a[3];}' \
  | LC_ALL=C sort -u \
  | awk -v out="${OUTPUT}" -F"\t" '{base=out"/"$1"/"$2; print base"/fa\n"base"/fq\n"base"/sketch\n"base"/prokka\n"base"/tetra\n"base"/stats\n"base"/kraken"}' \
  | xargs -r -n 50 -P "$(( THREADS>16?16:THREADS ))" mkdir -p
  
  # -------------------------------
  # Helpers
  # -------------------------------
  exists_nonempty(){ [[ -s "$1" ]]; }
  log(){ echo "[$(date +%H:%M:%S)] $*"; }
  
  # Validate/compress cache (idempotent). Returns path in CACHE_FASTQ
  validate_fastq(){
    local src="$1" base="$(basename "$1")" cached="${CACHE_FASTQ}/${base}"
    if exists_nonempty "${cached}"; then printf '%s\n' "${cached}"; return 0; fi
    if gzip -t "${src}" 2>/dev/null; then
    ln -sf "${src}" "${cached}" && printf '%s\n' "${cached}" && return 0
    fi
    # try repair via reformat.sh (lossless rewrite)
    local tmp
    tmp=$(mktemp --suffix=.fastq.gz)
    "${BBMAP}/reformat.sh" in="${src}" out="${tmp}" ow 1>/dev/null 2>&1 || true
    if gzip -t "${tmp}" 2>/dev/null; then
    mv "${tmp}" "${cached}" && printf '%s\n' "${cached}" && return 0
    fi
    rm -f "${tmp}" || true
    return 1
  }
  
  # Core per-file pipeline (streaming; minimal disk I/O)
  process_one(){
    local fq="$1"
    local base="$(basename "${fq}" .fastq.gz)"
    local fc="$(awk -F'_' '{print $1}' <<<"${base}")"
    local barcode="$(awk -F'_' '{print $3}' <<<"${base}")"
    local bcdir="${OUTPUT}/${fc}/${barcode}"
    local fafile="${bcdir}/fa/${base}.fa"
    local fqout="${bcdir}/fq/${base}.fastq.gz"  # keep trimmed gz for provenance
    
    # Skip if final fasta already exists and non-empty
    if exists_nonempty "${fafile}"; then
    echo "SKIP: ${base} (exists)"; return 0
    fi
    
    echo "RUN : ${base}"
    
    # Stream: bbduk -> filtlong -> reformat (to fasta)
    # bbduk threads auto; set explicitly for consistency
    # Note reformat fastawrap=0 to avoid line wraps
    # Also tee gz of trimmed fastq without extra pass (optional if pigz present)
    if command -v pigz >/dev/null 2>&1; then
    "${BBMAP}/bbduk.sh" in="${fq}" out=stdout ref=adapters,artifacts,phix,lambda qtrim=rl trimq=15 entropy=0.75 qin=33 minlength="${MIN_READLEN}" t="${THREADS}" \
    2>>"${bcdir}/log.txt" \
    | "${FILTLONG}" --min_length "${MIN_READLEN}" --keep_percent "${KEEP_PCT}" \
    2>>"${bcdir}/log.txt" \
    | tee >(pigz -c > "${fqout}") \
    | "${BBMAP}/reformat.sh" in=stdin out="${fafile}.tmp.fa" fastawrap=0 ow 1>>"${bcdir}/log.txt" 2>&1
    else
      "${BBMAP}/bbduk.sh" in="${fq}" out=stdout ref=adapters,artifacts,phix,lambda qtrim=rl trimq=15 entropy=0.75 qin=33 minlength="${MIN_READLEN}" t="${THREADS}" \
    2>>"${bcdir}/log.txt" \
    | "${FILTLONG}" --min_length "${MIN_READLEN}" --keep_percent "${KEEP_PCT}" \
    2>>"${bcdir}/log.txt" \
    | tee >(gzip -c > "${fqout}") \
    | "${BBMAP}/reformat.sh" in=stdin out="${fafile}.tmp.fa" fastawrap=0 ow 1>>"${bcdir}/log.txt" 2>&1
    fi
    
    # Simplify headers (1 pass)
    awk '/^>/{sub(/ .*/,"",$0)}1' "${fafile}.tmp.fa" > "${fafile}"
    rm -f "${fafile}.tmp.fa"
    
    # ---- Optional sub-tools ----
    if (( RUN_SKETCH )); then
    local sk_out="${bcdir}/sketch/${base}.txt"
    if [[ ! -s "${sk_out}" ]]; then
    "${BBMAP}/sendsketch.sh" in="${fafile}" address=nt out="${sk_out}" format=3 \
    1>>"${bcdir}/log.txt" 2>&1 || true
    Rscript "${DANADIR}/sketch-db.r" "${bcdir}" 1>>"${bcdir}/log.txt" 2>&1 || true
    fi
    fi
    
    if (( RUN_KRAKEN )); then
    local k_dir="${bcdir}/kraken"; mkdir -p "${k_dir}"
    local k_tsv="${k_dir}/${base}.tsv"
    if [[ ! -e "${k_tsv}" ]]; then
    echo "[KRAKEN] DB: ${KRAKEN_DB}" >>"${bcdir}/log.txt"
    "${KRAKEN2}" --db "${KRAKEN_DB}" --use-names --threads 1 --report "${k_dir}/${base}.report" "${fafile}" 2>>"${bcdir}/log.txt" \
    | gawk -f "${DANADIR}/kraken_parse.awk" > "${k_tsv}" || true
    Rscript "${DANADIR}/kraken-db.r" "${bcdir}" 1>>"${bcdir}/log.txt" 2>&1 || true
    Rscript "${DANADIR}/krakenreport-db.r" "${bcdir}" 1>>"${bcdir}/log.txt" 2>&1 || true
    fi
    fi
    
    if (( RUN_PROKKA )); then
    local prokdir="${bcdir}/prokka/${base}"; mkdir -p "${prokdir}"
    shopt -s nullglob
    local have=("${prokdir}/PROKKA_"*.tsv)
    if (( ${#have[@]} == 0 )); then
      "${PROKKA_BIN}" --metagenome --fast --cpus "${PROKKA_THREADS}" --evalue 1e-20 \
      --outdir "${prokdir}" --force --quiet "${fafile}" 1>>"${bcdir}/log.txt" 2>&1 || true
      rm -f "${prokdir}"/*.err "${prokdir}"/*.fna "${prokdir}"/*.fsa "${prokdir}"/*.gbk \
      "${prokdir}"/*.log "${prokdir}"/*.sqn "${prokdir}"/*.txt || true
      Rscript "${DANADIR}/prokka-db.r" "${bcdir}" 1>>"${bcdir}/log.txt" 2>&1 || true
      fi
      fi
      
      if (( RUN_TETRA )); then
      local t_dir="${bcdir}/tetra"; mkdir -p "${t_dir}"
      local lrn="${t_dir}/${base}.lrn"
      if [[ ! -s "${lrn}" ]]; then
      grep '>' "${fafile}" | sed 's/>//' | paste - - - > "${t_dir}/annotation.${base}.txt"
      perl "${APPS}/tetramer_freqs_esom.pl" -f "${fafile}" -a "${t_dir}/annotation.${base}.txt" -min "${MIN_READLEN}" -max 10000000 \
      1>>"${bcdir}/log.txt" 2>&1 || true
      if [[ ! -s "${bcdir}/tnfs.txt" ]]; then
      paste <(echo "seqid") <(head -n4 Tetra_*.lrn | tail -n1 | cut -f2-) > "${bcdir}/tnfs.txt" || true
      fi
      paste <(awk '$1!~/^%/' Tetra_${base}*.names) <(awk '$1!~/^%/' Tetra_${base}*.lrn) | cut -f3,5- > "${lrn}" || true
      rm -f Tetra_${base}.* "${t_dir}/annotation.${base}.txt" || true
      Rscript "${DANADIR}/tetra-db.r" "${bcdir}" 1>>"${bcdir}/log.txt" 2>&1 || true
      fi
      fi
    }
    export -f process_one validate_fastq exists_nonempty log
    export OUTPUT BBMAP DANADIR PROKKA_BIN FILTLONG KRAKEN2 KRAKEN_DB APPS RUN_SKETCH RUN_KRAKEN RUN_TETRA RUN_PROKKA PROKKA_THREADS MIN_READLEN KEEP_PCT THREADS
    
    # -------------------------------
    # Validate/cache all FASTQs (parallel)
    # -------------------------------
    start_ts=$(date +%s)
    printf '%s\0' "${FILES[@]}" \
    | xargs -0 -I{} -P "${THREADS}" bash -c 'validate_fastq "$@" >/dev/null || echo "[BAD] $1" >&2' _ {}
    
    # After validation, map to cached paths
    mapfile -t CACHED < <(printf '%s\n' "${FILES[@]}" | awk -v c="${CACHE_FASTQ}" -F'/' '{print c"/"$NF}')
    
    echo "[INFO] Processing with ${THREADS} parallel workers"
    
    # -------------------------------
    # Process batch with time cap
    # -------------------------------
    run_batch(){
      local now elapsed
      now=$(date +%s); elapsed=$(( now - start_ts ))
      if (( elapsed > MAX_DURATION )); then echo "[WARN] Time cap hit: ${elapsed}s > ${MAX_DURATION}s"; return 0; fi
      if command -v parallel >/dev/null 2>&1; then
      parallel -j "${THREADS}" --halt now,fail=1 process_one ::: "${CACHED[@]}"
      else
        printf '%s\n' "${CACHED[@]}" | xargs -r -n 1 -P "${THREADS}" bash -c 'process_one "$@"' _
      fi
    }
    run_batch
    
    echo "[DONE] Output => ${OUTPUT}"
    