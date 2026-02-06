#!/usr/bin/env bash
# quick-nano-barcode.fast.sh — fast, simple QC + annotation for Nanopore barcoded reads

set -euo pipefail
IFS=$'\n\t'

# Defaults & Env
THREADS=${THREADS:-$(nproc)}
BBMAP=${BBMAP:-/work/apps/bbmap}
DANADIR=${DANADIR:-/work/apps/dana}
PROKKA_BIN=${PROKKA_BIN:-/work/apps/prokka/bin/prokka}
FILTLONG=${FILTLONG:-/work/apps/Filtlong/bin/filtlong}
KRAKEN2=${KRAKEN2:-/usr/bin/kraken2}
KRAKEN_DB=${KRAKEN_DB:-/data/scratch/refdbs/krakendb/pluspfp_08gb}
APPS=${APPS:-/work/apps}

MIN_SIZE="1M"
RUN_SKETCH=0
RUN_KRAKEN=0
RUN_TETRA=0
RUN_PROKKA=0
PROKKA_THREADS=${PROKKA_THREADS:-1}
MAX_DURATION=${MAX_DURATION:-3600}
MIN_READLEN=${MIN_READLEN:-1500}
KEEP_PCT=${KEEP_PCT:-80}

print_help(){ cat <<EOF
Fast QC + annotation for barcoded Nanopore FASTQ(.gz) in directories containing fastq_pass/.

Required:
  -i, --input DIR       Project dir containing run/barcode outputs
Optional:
  -o, --output DIR      Output root (default: out_<basename(input)>_<timestamp>)
  -t, --threads N       Global threads (default: nproc)
  -S, --sketch          Run BBMap sendsketch + Dana DB update
  -K, --kraken          Run Kraken2 + Dana DB updates  
  -P, --prokka          Run Prokka + Dana DB update (metagenome mode)
  -T, --tetra           Run tetramer ESOM pipeline + Dana DB
  --min-size            Minimum size of FASTQ file to process, e.g. 1k, 1M
  --max-duration SEC    Stop after SEC for this batch (default: ${MAX_DURATION})
  --min-readlen N       Min read length for keep (default: ${MIN_READLEN})
  --keep-pct P          Filtlong keep percent (default: ${KEEP_PCT})
  --deterministic       Process files in md5(path) order (stable "shuffle")
  -h, --help            This help
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
    --min-size)     MIN_SIZE="$2"; shift 2 ;;
    --max-duration) MAX_DURATION="$2"; shift 2 ;;
    --min-readlen)  MIN_READLEN="$2"; shift 2 ;;
    --keep-pct)     KEEP_PCT="$2"; shift 2 ;;
    --deterministic) DETERMINISTIC=1; shift ;;
    -h|--help) print_help; exit 0 ;;
    *) echo "Unknown arg: $1" >&2; print_help; exit 2 ;;
  esac
done

[[ -z "${INPUT}" ]] && { echo "[ERR] --input is required" >&2; exit 2; }
[[ ! -d "${INPUT}" ]] && { echo "[ERR] input dir missing: ${INPUT}" >&2; exit 2; }
[[ -z "${OUTPUT}" ]] && OUTPUT="out_$(basename "${INPUT}")_$(date +%Y%m%d_%H%M%S)"

mkdir -p "${OUTPUT}"

# Collect FASTQ files
echo "[INFO] $(date +%H:%M:%S) Getting FASTQ files"
mapfile -d '' FILES < <(
  find "${INPUT}" -follow -type f -name '*fastq.gz' -size +"${MIN_SIZE}" \
    -regex '.*/fastq_pass/.*' -path '*barcode*' \
    ! -path '*/data/minknow/*' ! -path "*${OUTPUT}*" ! -name '*.tmp.*' -print0
)

NUMFILES=${#FILES[@]}
echo "[INFO] $(date +%H:%M:%S) Found ${NUMFILES} FASTQ files"
(( NUMFILES == 0 )) && { echo "[WARN] No FASTQ files found"; exit 0; }

# Optional deterministic sorting
if (( DETERMINISTIC )); then
  echo "[INFO] $(date +%H:%M:%S) Sorting files deterministically"
  mapfile -t FILES < <(
    printf '%s\n' "${FILES[@]}" | python3 -c "
import sys, hashlib
paths = [line.strip() for line in sys.stdin if line.strip()]
for path in sorted(paths, key=lambda p: hashlib.md5(p.encode()).hexdigest()):
    print(path)
"
  )
fi

# Create output directories
echo "[INFO] $(date +%H:%M:%S) Creating directories"
printf '%s\n' "${FILES[@]}" \
| xargs -I {} basename {} .fastq.gz \
| cut -d'_' -f1,3 \
| awk -F'_' '$2 ~ /^barcode[0-9]+$/ {print $1"/"$2}' \
| sort -u \
| xargs -I {} mkdir -p "${OUTPUT}/{}"/{fa,fq,sketch,prokka,tetra,stats,kraken}

process_one(){
  local fq="$1"
  local base="$(basename "${fq}" .fastq.gz)"
  local fc="$(cut -d'_' -f1 <<<"${base}")"
  local barcode="$(cut -d'_' -f3 <<<"${base}")"
  
  [[ ! "$barcode" =~ ^barcode[0-9]+$ ]] && return 0
  
  local bcdir="${OUTPUT}/${fc}/${barcode}"
  local fafile="${bcdir}/fa/${base}.fa"
  local fqout="${bcdir}/fq/${base}.fastq.gz"
  local logfile="${bcdir}/log.txt"

  [[ -s "${fafile}" ]] && return 0

  echo -n "RUN: ${base} : "

  # Main pipeline: bbduk -> filtlong -> save fastq -> convert to fasta
  echo -n "BBDUK "
  local tmpfq="${fqout}.tmp.fq"
  local tmpfa="${fafile}.tmp.fa"
  
  # Step 1: bbduk + filtlong -> temp fastq
  if "${FILTLONG}" --min_length "${MIN_READLEN}" --keep_percent "${KEEP_PCT}" \
    <( "${BBMAP}/bbduk.sh" in="${fq}" out=stdout.fq \
          ref=adapters,artifacts,phix,lambda qtrim=rl trimq=15 entropy=0.75 \
          qin=33 minlength="${MIN_READLEN}" threads=1 2>>"${logfile}" ) 2>>"${logfile}" \
    > "${tmpfq}"; then
    
    # Step 2: convert to fasta  
    if "${BBMAP}/reformat.sh" in="${tmpfq}" out="${tmpfa}" fastawrap=0 ow >>"${logfile}" 2>&1; then
      # Simplify fasta headers and compress fastq
      awk '/^>/{sub(/ .*/,"")} 1' "${tmpfa}" > "${fafile}"
      gzip -c "${tmpfq}" > "${fqout}"
      rm -f "${tmpfa}" "${tmpfq}"
    else
      touch "${fafile}"
#      rm -f "${tmpfa}" "${tmpfq}"
    fi
  else
    touch "${fafile}" 
#    rm -f "${tmpfq}" "${tmpfa}"
  fi

  [[ ! -s "${fafile}" ]] && { echo "EMPTY"; return 0; }

  # Optional analyses
  if (( RUN_SKETCH )); then
    echo -n "SKETCH "
    local sk_out="${bcdir}/sketch/${base}.txt"
    if [[ ! -s "${sk_out}" ]]; then
      "${BBMAP}/sendsketch.sh" in="${fafile}" address=nt out="${sk_out}" format=3 \
        >>"${logfile}" 2>&1 && \
      Rscript "${DANADIR}/sketch-db.r" "${bcdir}" >>"${logfile}" 2>&1 || true
    fi
  fi

  if (( RUN_KRAKEN )); then
    echo -n "KRAKEN "
    local k_tsv="${bcdir}/kraken/${base}.tsv"
    if [[ ! -e "${k_tsv}" ]]; then
      "${KRAKEN2}" --db "${KRAKEN_DB}" --use-names --threads 1 \
        --report "${bcdir}/kraken/${base}.report" "${fafile}" 2>>"${logfile}" \
      | gawk -f "${DANADIR}/kraken_parse.awk" > "${k_tsv}" 2>>"${logfile}" && {
        Rscript "${DANADIR}/kraken-db.r" "${bcdir}" >>"${logfile}" 2>&1
        Rscript "${DANADIR}/krakenreport-db.r" "${bcdir}" >>"${logfile}" 2>&1
      } || true
    fi
  fi

  if (( RUN_PROKKA )); then
    echo -n "PROKKA "
    local prokdir="${bcdir}/prokka/${base}"
    mkdir -p "${prokdir}"
    if ! ls "${prokdir}/PROKKA_"*.tsv >/dev/null 2>&1; then
      "${PROKKA_BIN}" --metagenome --fast --cpus "${PROKKA_THREADS}" --evalue 1e-20 \
        --outdir "${prokdir}" --force --quiet "${fafile}" >>"${logfile}" 2>&1 && {
        rm -f "${prokdir}"/*.{err,fna,fsa,gbk,log,sqn,txt} 2>/dev/null
        Rscript "${DANADIR}/prokka-db.r" "${bcdir}" >>"${logfile}" 2>&1
      } || true
    fi
  fi

  if (( RUN_TETRA )); then
    echo -n "TETRA "
    local t_dir="${bcdir}/tetra"
    local lrn="${t_dir}/${base}.lrn"
    mkdir -p "${t_dir}"
    if [[ ! -s "${lrn}" ]]; then
      grep '>' "${fafile}" | sed 's/>//' > "${t_dir}/annotation.${base}.txt" && \
      perl "${APPS}/tetramer_freqs_esom.pl" -f "${fafile}" \
        -a "${t_dir}/annotation.${base}.txt" -min "${MIN_READLEN}" -max 10000000 \
        >>"${logfile}" 2>&1 && {
        
        [[ ! -s "${bcdir}/tnfs.txt" ]] && ls Tetra_*.lrn >/dev/null 2>&1 && \
          paste <(echo "seqid") <(head -n4 Tetra_*.lrn | tail -n1 | cut -f2-) > "${bcdir}/tnfs.txt"
        
        ls Tetra_${base}*.{names,lrn} >/dev/null 2>&1 && \
          paste <(awk '$1!~/^%/' Tetra_${base}*.names) <(awk '$1!~/^%/' Tetra_${base}*.lrn) \
          | cut -f3,5- > "${lrn}"
        
        rm -f Tetra_${base}.* "${t_dir}/annotation.${base}.txt"
        Rscript "${DANADIR}/tetra-db.r" "${bcdir}" >>"${logfile}" 2>&1
      } || true
    fi
  fi
  
  echo "DONE"
}

export -f process_one
export OUTPUT BBMAP DANADIR PROKKA_BIN FILTLONG KRAKEN2 KRAKEN_DB APPS 
export RUN_SKETCH RUN_KRAKEN RUN_TETRA RUN_PROKKA PROKKA_THREADS MIN_READLEN KEEP_PCT

echo "[INFO] $(date +%H:%M:%S) Processing with ${THREADS} parallel workers"
start_ts=$(date +%s)

# Process files with time cap
process_batch(){
  local elapsed=$(( $(date +%s) - start_ts ))
  if (( elapsed > MAX_DURATION )); then
    echo "[WARN] Time cap hit: ${elapsed}s > ${MAX_DURATION}s"
    return 0
  fi

  printf '%s\0' "${FILES[@]}" | parallel --null -j "${THREADS}" --halt now,fail=1 --bar process_one {}
}

process_batch
echo "[DONE] Output => ${OUTPUT}"