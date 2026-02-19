#!/usr/bin/env bash
# quick-nano-barcode.fast.sh — fast, simple QC + annotation for Nanopore barcoded reads
# Speed-ups: parallelism, streaming I/O, minimal temp files, idempotent outputs
# Requirements: GNU coreutils, awk, grep, sort, xargs, find, gzip; BBMap (bbduk.sh, reformat.sh),
#   Filtlong, Kraken2 (+ DB), Prokka, Rscript; optional: GNU parallel, pigz.
#
# Usage:
#   ./quick-nano-barcode.fast.sh -i /data/project_QEI2025 -o out_dana -S -K -P -T -t 8 --deterministic
#
# Deterministic mode note:
#   We sort by md5(PATH) to get a reproducible "shuffle" that interleaves flowcells in live runs.

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

VALIDATION=0
RUN_SKETCH=0
RUN_KRAKEN=0
RUN_TETRA=0
RUN_PROKKA=0
PROKKA_THREADS=${PROKKA_THREADS:-0}  # 0 = auto

MAX_DURATION=${MAX_DURATION:-3600}  # 1h safety cap per batch
MIN_READLEN=${MIN_READLEN:-1500}
KEEP_PCT=${KEEP_PCT:-80}
MIN_SIZE=${MIN_SIZE:-1000000}

print_help(){ cat <<EOF
Fast QC + annotation for barcoded Nanopore FASTQ(.gz) in directories containing fastq_pass/.

Required:
  -i, --input DIR       Project dir containing run/barcode outputs
Optional:
  -o, --output DIR      Output root (default: out_<basename(input)>_<rand>)
  -t, --threads N       Global threads (default: nproc)
  -S, --sketch          Run BBMap sendsketch + Dana DB update
  -K, --kraken          Run Kraken2 + Dana DB updates
  -P, --prokka          Run Prokka + Dana DB update (metagenome mode)
  -T, --tetra           Run tetramer ESOM pipeline + Dana DB
  -V, --validation      Validate input FASTQS
  --min-size            Minimum size of FASTQ file to process, e.g. 1k, 1M
  --max-duration SEC    Stop after SEC for this batch (default: ${MAX_DURATION})
  --min-readlen N       Min read length for keep (default: ${MIN_READLEN})
  --keep-pct P          Filtlong keep percent (default: ${KEEP_PCT})
  --deterministic       Process files in md5(path) order (stable "shuffle")
  -h, --help            This help

Env overrides:
  BBMAP, DANADIR, PROKKA_BIN, FILTLONG, KRAKEN2, KRAKEN_DB, THREADS, PROKKA_THREADS
EOF
}

INPUT=""
OUTPUT=""
DETERMINISTIC=1

while (( $# )); do
  case "$1" in
    -i|--input)   INPUT="$2"; shift 2 ;;
    -o|--output)  OUTPUT="$2"; shift 2 ;;
    -t|--threads) THREADS="$2"; shift 2 ;;
    -S|--sketch)  RUN_SKETCH=1; shift ;;
    -K|--kraken)  RUN_KRAKEN=1; shift ;;
    -P|--prokka)  RUN_PROKKA=1; shift ;;
    -T|--tetra)   RUN_TETRA=1; shift ;;
    -V|--validation) VALIDATION=1; shift ;;
    --min-size)     MIN_SIZE="$2"; shift 2 ;;
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

CACHE_FASTQ="/data/.fastq_pass"
mkdir -p "${CACHE_FASTQ}"

# -------------------------------
# Collect FASTQ files
# -------------------------------
echo "[INFO] $(date) Getting FASTQ files"
mapfile -d '' FILES < <(
  find "${INPUT}" -follow -type f -name '*fastq.gz' -size +"${MIN_SIZE}" \
    -regex '.*/fastq_pass/.*' -path '*barcode*' \
    ! -path '*/data/minknow/*' ! -path "*${OUTPUT}*" ! -name '*.tmp.*' -print0
)

NUMFILES=${#FILES[@]}
echo "[INFO] $(date) Found ${NUMFILES} FASTQ files"

echo "[INFO] $(date) Sorting FASTQ files"
# Deterministic pseudo-random shuffle by md5(PATH)
if (( DETERMINISTIC )); then
mapfile -t FILES < <(
  printf '%s\n' "${FILES[@]}" \
  | python3 -c 'import sys,hashlib;print("\n".join(p for _,p in sorted((hashlib.md5(s.encode()).hexdigest(),s) for s in sys.stdin.read().splitlines())))'
)
NUMFILES=${#FILES[@]}
echo "[INFO] $(date) Sorted ${NUMFILES} FASTQ files"

fi

#( printf '%s\n' "${FILES[@]}" | head ) || :


# -------------------------------
# Pre-create only needed dirs (flowcell/barcode)
# -------------------------------
echo "[INFO] $(date) Creating new directories"

# One-pass, deduped, safe; no sort/xargs needed
awk -v out="${OUTPUT}" -F'/' '
  BEGIN { OFS="/" }
  {
    f = $NF
    # flowcell = prefix before _pass_
    if (match(f, /^([^_]+)_pass_/, m) && match(f, /_(barcode[0-9]+)_/, b)) {
      base = out "/" m[1] "/" b[1]
      if (!seen[base]++) {
        print base "/fa"
        print base "/fq"
        print base "/sketch"
        print base "/prokka"
        print base "/tetra"
        print base "/stats"
        print base "/kraken"
      }
    }
  }
' < <(printf '%s\n' "${FILES[@]}") | while IFS= read -r d; do
  # Skip accidental empties just in case
  [[ -n "$d" && "$d" != "/" ]] && mkdir -p -- "$d"
done

# -------------------------------
# Helpers
# -------------------------------
exists_nonempty(){ [[ -s "$1" ]]; }
log(){ echo "[$(date +%H:%M:%S)] $*"; }

# Validate/compress cache (idempotent). Returns path in CACHE_FASTQ
validate_fastq(){
  local src="$1"
  local base="$(basename "$1")"
  local cached="${CACHE_FASTQ}/${base}"

  if exists_nonempty "${cached}"; then printf '%s\n' "${cached}"; return 0; fi
  if gzip -t "${src}" 2>/dev/null; then
    ln -sf "${src}" "${cached}" && printf '%s\n' "${cached}" && return 0
  fi
  local tmp
  tmp=$(mktemp --suffix=.fastq.gz)
  "${BBMAP}/reformat.sh" in="${src}" out="${tmp}" ow #1>/dev/null 2>&1 || true
  if gzip -t "${tmp}" 2>/dev/null; then
    mv "${tmp}" "${cached}" && printf '%s\n' "${cached}" && return 0
  fi
  rm -f "${tmp}" || true
  return 1
}

process_one(){
  local fq="$1"
  local base="$(basename "${fq}" .fastq.gz)"
  local fc="$(awk -F'_' '{print $1}' <<<"${base}")"
  local barcode="$(awk -F'_' '{print $3}' <<<"${base}")"
  local bcdir="${OUTPUT}/${fc}/${barcode}"
  local fafile="${bcdir}/fa/${base}.fa"
  local fqout="${bcdir}/fq/${base}.fastq.gz"

#echo $fq $base $fc $barcode $bcdir $fafile $fqout

if [[ -e "${fafile}" ]]; then
#  echo "SKIP: ${base} (exists)"; 
  return 0
fi

  echo -n "RUN : ${base}  : "

  # Stream: bbduk -> filtlong -> reformat -> write tmp fasta + trimmed fastq
    #2>>"${bcdir}/log.txt" \
    echo -n "BBDUK "
    

{
  "${FILTLONG}" --min_length "${MIN_READLEN}" --keep_percent "${KEEP_PCT}" \
    <( "${BBMAP}/bbduk.sh" \
          in="${fq}" out=stdout.fq \
          ref=adapters,artifacts,phix,lambda \
          qtrim=rl trimq=15 entropy=0.75 qin=33 minlength="${MIN_READLEN}" t="${THREADS}" ) \
  | tee >( command -v pigz >/dev/null && pigz -c > "${fqout}" || gzip -c > "${fqout}" ) \
  | "${BBMAP}/reformat.sh" in=stdin.fq int=f out="${fafile}.tmp.fa" fastawrap=0 ow

} >>"${bcdir}/log.txt" 2>&1


  # Simplify headers if tmp exists; else create empty marker
  if [[ -s "${fafile}.tmp.fa" ]]; then
    awk '/^>/{sub(/ .*/,"",$0)}1' "${fafile}.tmp.fa" > "${fafile}"
  else
    : > "${fafile}"
  fi
  rm -f "${fafile}.tmp.fa"

  if [[ ! -s "${fafile}" ]]; then
    #echo "SKIP: ${base} (empty)"; 
    return 0
  fi

  echo "here"

  # Optional tooling
  if (( RUN_SKETCH )); then
    echo -n "SKETCH "
    local sk_out="${bcdir}/sketch/${base}.txt"
    if [[ ! -s "${sk_out}" ]]; then
      "${BBMAP}/sendsketch.sh" in="${fafile}" address=nt out="${sk_out}" format=3 \
      #  1>>"${bcdir}/log.txt" 2>&1 || true
      Rscript "${DANADIR}/sketch-db.r" "${bcdir}" #1>>"${bcdir}/log.txt" 2>&1 || true
    fi
  fi

  if (( RUN_KRAKEN )); then
      #2>>"${bcdir}/log.txt" \
    echo -n "KRAKEN "
    local k_dir="${bcdir}/kraken"
    local k_tsv="${k_dir}/${base}.tsv"
    if [[ ! -e "${k_tsv}" ]]; then
      echo "[KRAKEN] DB: ${KRAKEN_DB}" #>>"${bcdir}/log.txt"
      "${KRAKEN2}" --db "${KRAKEN_DB}" --use-names --threads 1 --report "${k_dir}/${base}.report" "${fafile}" \
        | gawk -f "${DANADIR}/kraken_parse.awk" > "${k_tsv}" || true
      Rscript "${DANADIR}/kraken-db.r" "${bcdir}" #1>>"${bcdir}/log.txt" 2>&1 || true
      Rscript "${DANADIR}/krakenreport-db.r" "${bcdir}" #1>>"${bcdir}/log.txt" 2>&1 || true
    fi
  fi

  if (( RUN_PROKKA )); then
    echo -n "PROKKA "
    local prokdir="${bcdir}/prokka/${base}"
    shopt -s nullglob
    local have=("${prokdir}/PROKKA_"*.tsv)
    if (( ${#have[@]} == 0 )); then
      "${PROKKA_BIN}" --metagenome --fast --cpus "${PROKKA_THREADS}" --evalue 1e-20 \
        --outdir "${prokdir}" --force --quiet "${fafile}" #1>>"${bcdir}/log.txt" 2>&1 || true
      rm -f "${prokdir}"/*.err "${prokdir}"/*.fna "${prokdir}"/*.fsa "${prokdir}"/*.gbk \
            "${prokdir}"/*.log "${prokdir}"/*.sqn "${prokdir}"/*.txt || true
      Rscript "${DANADIR}/prokka-db.r" "${bcdir}" #1>>"${bcdir}/log.txt" 2>&1 || true
    fi
  fi

  if (( RUN_TETRA )); then
    echo -n "TETRA "
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
      Rscript "${DANADIR}/tetra-db.r" "${bcdir}" #1>>"${bcdir}/log.txt" 2>&1 || true
    fi
  fi
  echo ""
}

export -f process_one validate_fastq exists_nonempty log
export OUTPUT BBMAP DANADIR PROKKA_BIN FILTLONG KRAKEN2 KRAKEN_DB APPS RUN_SKETCH RUN_KRAKEN RUN_TETRA RUN_PROKKA PROKKA_THREADS MIN_READLEN KEEP_PCT THREADS CACHE_FASTQ VALIDATE MIN_SIZE

# -------------------------------
# Validate/cache all FASTQs (parallel)
# -------------------------------
start_ts=$(date +%s)

#skip validation
if (( VALIDATION )); then
echo "[INFO] $(date) Validating FASTQ files"
printf '%s\0' "${FILES[@]}" | parallel --null -j "${THREADS}" --bar 'validate_fastq {} >/dev/null || echo "[BAD] {}" >&2'
#printf '%s\0' "${FILES[@]}" | parallel --null -j "${THREADS}" --bar 'validate_fastq {} || echo "[BAD] {}"'
fi

# After validation, map to cached paths
mapfile -t CACHED < <(printf '%s\n' "${FILES[@]}" | awk -v c="${CACHE_FASTQ}" -F'/' '{print c"/"$NF}')

echo "[INFO] $(date) Processing with ${THREADS} parallel workers"

# -------------------------------
# Process batch with time cap
# -------------------------------
run_batch(){
  local now elapsed
  now=$(date +%s); elapsed=$(( now - start_ts ))
  (( elapsed > MAX_DURATION )) && { echo "[WARN] Time cap hit: ${elapsed}s > ${MAX_DURATION}s"; return 0; }

    # GNU parallel: read NUL-delimited args from stdin (no ::: expansion)
    printf '%s\0' "${CACHED[@]}" \
    | parallel --null -j "${THREADS}" --halt now,fail=1 --bar process_one {}

}
run_batch

echo "[DONE] Output => ${OUTPUT}"
