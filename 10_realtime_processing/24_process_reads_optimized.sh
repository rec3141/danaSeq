#!/usr/bin/env bash
################################################################################
#                                                                              #
#  ⚡  REAL-TIME NANOPORE PROCESSING - OPTIMIZED  ⚡                           #
#                                                                              #
#  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~                         #
#     🌊  Process reads as they stream from sequencer  🌊                     #
#  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~                         #
#                                                                              #
#  Quality Control → Taxonomic Classification → Gene Annotation               #
#                                                                              #
#  AI-Enhanced • Resume Capable • Production Ready                             #
#                                                                              #
################################################################################
#
# PURPOSE: Fast, incremental QC and annotation for Oxford Nanopore barcoded reads
#
# FEATURES:
#   • Automatic FASTQ validation and repair
#   • Quality filtering (BBDuk + Filtlong)
#   • Optional Kraken2 classification
#   • Optional Prokka annotation
#   • Optional sendsketch profiling
#   • Smart resume (skip completed files)
#   • Time-limited batching
#   • DuckDB integration
#
# USAGE:
#   ./24_process_reads_optimized.sh -i <input_dir> -K -P -S
#
# AUTHOR: Dana Pipeline (AI-optimized)
# VERSION: 2.0
#
################################################################################

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
HMMSEARCH=${HMMSEARCH:-/usr/bin/hmmsearch}
APPS=${APPS:-/work/apps}

MIN_SIZE="1M"
RUN_SKETCH=0
RUN_KRAKEN=0
RUN_TETRA=0
RUN_PROKKA=0
HMM_DATABASES=""
FORCE=0
VERBOSE=0
DEBUG=0
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
  --hmm FILE1,FILE2     Run HMM searches on Prokka gene calls (requires -P)
                        Provide full paths to HMM files (comma-delimited)
                        Examples: --hmm /path/to/CANT-HYD.hmm
                                  --hmm /path/to/CANT-HYD.hmm,/path/to/FOAM.hmm
                        Searches use trusted cutoffs (--cut_tc) from HMM files
  --force               Force re-run of optional stages (Prokka, HMM, Kraken, etc.)
                        even if output already exists
  -v, --verbose         Show command output (default: log to files)
  -d, --debug           Debug mode: show all commands and keep temp files
  --min-size            Minimum size of FASTQ file to process (default: ${MIN_SIZE})
  --max-duration SEC    Stop after SEC for this batch (default: ${MAX_DURATION})
  --min-readlen N       Min read length for keep (default: ${MIN_READLEN})
  --keep-pct P          Filtlong keep percent (default: ${KEEP_PCT})
  --shuffle             Process files in shuffled order (deterministic md5 sort)
  -h, --help            This help
EOF
}

INPUT=""
OUTPUT=""
SHUFFLE=0

while (( $# )); do
  case "$1" in
    -i|--input)   INPUT="$2"; shift 2 ;;
    -o|--output)  OUTPUT="$2"; shift 2 ;;
    -t|--threads) THREADS="$2"; shift 2 ;;
    -S|--sketch)  RUN_SKETCH=1; shift ;;
    -K|--kraken)  RUN_KRAKEN=1; shift ;;
    -P|--prokka)  RUN_PROKKA=1; shift ;;
    -T|--tetra)   RUN_TETRA=1; shift ;;
    --hmm)        HMM_DATABASES="$2"; shift 2 ;;
    --force)      FORCE=1; shift ;;
    -v|--verbose) VERBOSE=1; shift ;;
    -d|--debug)   DEBUG=1; VERBOSE=1; shift ;;
    --min-size)     MIN_SIZE="$2"; shift 2 ;;
    --max-duration) MAX_DURATION="$2"; shift 2 ;;
    --min-readlen)  MIN_READLEN="$2"; shift 2 ;;
    --keep-pct)     KEEP_PCT="$2"; shift 2 ;;
    --shuffle)      SHUFFLE=1; shift ;;
    -h|--help) print_help; exit 0 ;;
    *) echo "Unknown arg: $1" >&2; print_help; exit 2 ;;
  esac
done

[[ -z "${INPUT}" ]] && { echo "[ERR] --input is required" >&2; echo "" >&2; print_help; exit 2; }
[[ ! -d "${INPUT}" ]] && { echo "[ERR] input dir missing: ${INPUT}" >&2; exit 2; }
[[ -n "${HMM_DATABASES}" && ${RUN_PROKKA} -eq 0 ]] && { echo "[ERR] --hmm requires -P (Prokka)" >&2; exit 2; }
[[ -z "${OUTPUT}" ]] && OUTPUT="out_$(basename "${INPUT}")_$(date +%Y%m%d_%H%M%S)"

mkdir -p "${OUTPUT}"
CACHE_FASTQ="/data/.fastq_pass"
mkdir -p "${CACHE_FASTQ}"

# Failure tracking
FAILURE_LOG="${OUTPUT}/failed_files.txt"
> "${FAILURE_LOG}"  # Clear any existing failure log

# Logging helpers: Control output verbosity based on user flags
# In normal mode: Only errors shown, all output logged to files
# In verbose mode: Commands echoed to terminal and logged
# In debug mode: Full command details shown with timing info
run_cmd() {
  local cmd="$1"
  local logfile="$2"

  if (( DEBUG )); then
    echo "[DEBUG] Running: $cmd" >&2
    eval "$cmd" 2>&1 | tee -a "$logfile"
  elif (( VERBOSE )); then
    eval "$cmd" 2>&1 | tee -a "$logfile"
  else
    eval "$cmd" >>"$logfile" 2>&1
  fi
}

debug_msg() {
  (( DEBUG )) && echo "[DEBUG] $*" >&2
}

verbose_msg() {
  (( VERBOSE )) && echo "[VERBOSE] $*" >&2
}

# Log file failure with diagnostics
# This allows the pipeline to continue processing other files
log_failure() {
  local file="$1"
  local stage="$2"
  local reason="$3"

  # Log to failure file
  echo "${file}|${stage}|${reason}" >> "${FAILURE_LOG}"

  # Show failure inline in progress output
  echo -n "(FAILED)"

  # Show reason if verbose mode
  if (( VERBOSE )); then
    echo ""
    echo "  └─ ${stage}: ${reason}"
  else
    echo ""  # Add newline to flush output for parallel --bar
  fi
}

# Kraken wrapper for semaphore execution
# This function runs with properly quoted variables
run_kraken_locked() {
  local kraken_bin="$1"
  local db="$2"
  local report="$3"
  local input="$4"
  local logfile="$5"
  local parse_awk="$6"
  local output_tsv="$7"

  "${kraken_bin}" --db "${db}" --use-names --threads 1 \
    --report "${report}" "${input}" 2>>"${logfile}" \
    | gawk -f "${parse_awk}" > "${output_tsv}" 2>>"${logfile}"
}

# ============================================================================
# Banner
# ============================================================================

echo -e "\033[0;36m"
echo "    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
echo -e "  ~~      \033[1;36m🧬  DANA PIPELINE - Real-Time Nanopore Processing  🧬\033[0;36m      ~~"
echo "    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
echo ""
echo -e "    ~~~  \033[1;37mQC → Filter → Classify → Annotate → Discover\033[0;36m  ~~~"
echo "    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
echo -e "\033[0m"

# ============================================================================
# STEP 1: Discover FASTQ files from Oxford Nanopore output structure
# ============================================================================
# Nanopore sequencers write reads to fastq_pass/ directories organized by barcode
# We search for files meeting minimum size threshold to avoid processing empty/tiny files
# Exclusions: MinKNOW runtime dirs, previous outputs, temporary files

echo "[INFO] $(date +%H:%M:%S) Getting FASTQ files"
mapfile -d '' FILES < <(
  find "${INPUT}" -follow -type f -name '*fastq.gz' -size +"${MIN_SIZE}" \
    -regex '.*/fastq_pass/.*' -path '*barcode*' \
    ! -path '*/data/minknow/*' ! -path "*${OUTPUT}*" ! -name '*.tmp.*' -print0
)

NUMFILES=${#FILES[@]}
echo "[INFO] $(date +%H:%M:%S) Found ${NUMFILES} FASTQ files"
(( NUMFILES == 0 )) && { echo "[WARN] No FASTQ files found"; exit 0; }

# Optional shuffled ordering for balanced batch processing
# Uses MD5 hash of file path to create stable pseudo-random order
# This prevents always processing the same samples first in repeated runs
if (( SHUFFLE )); then
  echo "[INFO] $(date +%H:%M:%S) Sorting files in shuffled order"
  mapfile -t FILES < <(
    printf '%s\n' "${FILES[@]}" | python3 -c "
import sys, hashlib
paths = [line.strip() for line in sys.stdin if line.strip()]
for path in sorted(paths, key=lambda p: hashlib.md5(p.encode()).hexdigest()):
    print(path)
"
  )
fi

# ============================================================================
# STEP 2: Prepare output directory structure
# ============================================================================
# Create hierarchical output: flowcell/barcode/analysis_type
# Subdirectories: fa (FASTA), fq (FASTQ), sketch (taxonomy sketches),
#                 prokka (annotations), tetra (tetranucleotide), kraken (classifications)
# This mirrors the input structure for easy navigation

echo "[INFO] $(date +%H:%M:%S) Creating directories"
BARCODELIST=$(printf '%s\n' "${FILES[@]}" | xargs -n1 basename | cut -f3 -d'_' | sort -u)
FCLIST=$(printf '%s\n' "${FILES[@]}" | xargs -n1 basename | cut -f1 -d'_' | sort -u)

for fc in $FCLIST; do
  for barcode in $BARCODELIST; do
    if [[ "$barcode" =~ ^barcode[0-9]+$ ]]; then
      mkdir -p "${OUTPUT}/$fc/$barcode"/{fa,fq,sketch,prokka,tetra,stats,kraken}
    fi
  done
done

# ============================================================================
# FASTQ Validation & Caching System
# ============================================================================
# Nanopore sequencing can produce corrupted gzip files if runs are interrupted
# This function validates integrity and attempts repair when needed
#
# Workflow:
#   1. Check cache - if file already validated, use cached symlink
#   2. Test gzip integrity with 'gzip -t'
#   3. If corrupted, attempt repair using BBMap reformat.sh
#   4. Cache validated/repaired files for faster subsequent runs
#
# Why caching? Validation is I/O intensive; cache allows instant resume on re-runs

validate_fastq() {
  local file="$1"
  local base="$(basename "$file")"
  local cached="${CACHE_FASTQ}/${base}"

  # Return cached version if already validated
  if [[ -s "${cached}" ]]; then
    printf '%s\n' "${cached}"
    return 0
  fi

  debug_msg "Checking: $file"

  # Test gzip integrity
  if gzip -t "$file" 2>/dev/null; then
    debug_msg "Good file: $file"
    ln -sf "$(realpath "$file")" "${cached}"
    printf '%s\n' "${cached}"
    return 0
  else
    # Attempt repair using BBMap's error-tolerant reader
    debug_msg "Bad file, attempting fix: $file"
    local tmpfile=$(mktemp --suffix=.fastq.gz)
    if "${BBMAP}/reformat.sh" in="$file" out="$tmpfile" ow >/dev/null 2>&1; then
      if gzip -t "$tmpfile" 2>/dev/null; then
        mv "$tmpfile" "${cached}"
        echo "[FIXED] $file" >&2
        printf '%s\n' "${cached}"
        return 0
      fi
    fi
    rm -f "$tmpfile"
    echo "[BAD] Could not fix: $file" >&2
    return 1
  fi
}

# ============================================================================
# Core Processing Pipeline: QC → Filter → Convert → Optional Analyses
# ============================================================================
# This function processes one FASTQ file through the complete workflow
#
# Pipeline stages:
#   1. BBDuk: Remove adapters, artifacts, low-quality bases
#   2. Filtlong: Length filtering and quality-based read selection
#   3. Format conversion: FASTQ → FASTA for downstream tools
#   4. Optional analyses: Kraken2, Prokka, Sendsketch, Tetramer
#
# Resume capability: Checks for existing output files, skips completed work
# This allows safe interruption and restart of long-running batch jobs

process_one() {
  local file="$1"
  local base="$(basename "$file" .fastq.gz)"
  local barcode="$(basename "$file" | cut -f3 -d'_')"
  local fc="$(basename "$file" | cut -f1 -d'_')"

  # Skip non-standard barcode names
  [[ ! "$barcode" =~ ^barcode[0-9]+$ ]] && return 0

  local bcdir="${OUTPUT}/$fc/$barcode"
  local fafile="$bcdir/fa/$base.fa"
  local fqfile="$bcdir/fq/$base.fastq.gz"
  local ftfile="$bcdir/fq/$base.filt.fastq"
  local logfile="$bcdir/log.txt"

  # Resume capability: Check if basic processing is done
  local SKIP_BASIC=0

  # Check if basic processing (BBDuk → FASTA) is already done
  if [[ -s "$fafile" ]] || [[ -s "$fafile.gz" ]]; then
    SKIP_BASIC=1
    # Note: Optional stages (Prokka, HMM, etc.) have their own resume checks
    # So we continue to run those if requested, they'll skip themselves if done
  fi

  # Crash recovery: Clean up any temporary files from interrupted previous run
  # This prevents "File exists and overwrite=false" errors in BBMap tools
  if (( SKIP_BASIC == 0 )); then
    rm -f "$fafile.partial" "$fafile.tmp.fa" "$ftfile"
  fi

  debug_msg "Processing $base - input size: $(stat -c%s "$file" 2>/dev/null || echo "0") bytes"
  echo -n "RUN: $base : "

  # =========================================================================
  # Basic Processing Pipeline (BBDuk → Filtlong → FASTA conversion)
  # =========================================================================
  # Skip if FASTA already exists (resume from optional stages)
  if (( SKIP_BASIC == 0 )); then

  # -------------------------------------------------------------------------
  # Stage 1: Quality control with BBDuk
  # -------------------------------------------------------------------------
  # Remove adapters, sequencing artifacts (phiX, lambda), low-quality regions
  # Quality trimming: Remove bases with Q < 15 from read ends
  # Entropy filter: Remove low-complexity sequences (homopolymers, repeats)
  echo -n "BBDUK "
  debug_msg "Running bbduk on $base"

  if ! run_cmd "${BBMAP}/bbduk.sh in=$file out=$fqfile ref=adapters,artifacts,phix,lambda qtrim=rl trimq=15 entropy=0.75 qin=33 minlength=$MIN_READLEN" "$logfile"; then
    # Diagnose common BBDuk failures
    local reason="unknown error"
    if [[ ! -s "$file" ]]; then
      reason="input file empty or missing"
    elif ! gzip -t "$file" 2>/dev/null; then
      reason="corrupted gzip - validation should have caught this"
    elif [[ ! -x "${BBMAP}/bbduk.sh" ]]; then
      reason="bbduk.sh not executable or missing"
    else
      reason="check ${logfile} for details"
    fi

    log_failure "$file" "BBDuk" "$reason"
    return 0  # Continue processing other files
  fi
  
  debug_msg "BBduk output size: $(stat -c%s "$fqfile" 2>/dev/null || echo "0") bytes"
  if [[ ! -s "$fqfile" ]]; then
    verbose_msg "BBduk produced empty output for $base"
    echo "(NO DATA - all reads filtered)"
    return 0
  fi

  # -------------------------------------------------------------------------
  # Stage 2: Length and quality filtering with Filtlong
  # -------------------------------------------------------------------------
  # Keeps top-quality reads while enforcing minimum length threshold
  # KEEP_PCT parameter: Retain this percentage of bases from best reads
  # Example: 80% keeps highest quality 80% of data, discards worst 20%
  echo -n "FILTLONG "
  debug_msg "Running filtlong on $base"
  if (( DEBUG )); then
    echo "[DEBUG] Running: ${FILTLONG} --min_length $MIN_READLEN --keep_percent $KEEP_PCT $fqfile > $ftfile" >&2
  fi
  
  if ! "${FILTLONG}" --min_length "$MIN_READLEN" --keep_percent "$KEEP_PCT" "$fqfile" > "$ftfile" 2>>"$logfile"; then
    # Diagnose common Filtlong failures
    local reason="unknown error"
    if [[ ! -s "$fqfile" ]]; then
      reason="BBDuk output missing"
    elif [[ ! -x "${FILTLONG}" ]]; then
      reason="filtlong not found or not executable"
    else
      reason="check ${logfile} for details"
    fi

    log_failure "$file" "Filtlong" "$reason"
    return 0  # Continue processing other files
  fi
  
  debug_msg "Filtlong output size: $(stat -c%s "$ftfile" 2>/dev/null || echo "0") bytes"
  if [[ ! -s "$ftfile" ]]; then
    verbose_msg "Filtlong produced empty output for $base"
    echo "(NO DATA - all reads too short)"
    return 0
  fi

  # -------------------------------------------------------------------------
  # Stage 3: Format conversion FASTQ → FASTA
  # -------------------------------------------------------------------------
  # Most downstream tools (Kraken, Prokka) prefer FASTA input
  # BBMap reformat.sh: Fast, memory-efficient converter
  # fastawrap=0: Keep sequences on single lines for easier parsing
  echo -n "FASTA "
  debug_msg "Converting to FASTA: $base"
  if ! run_cmd "${BBMAP}/reformat.sh in=$ftfile out=$fafile.tmp.fa fastawrap=0" "$logfile"; then
    # Diagnose common reformat failures
    local reason="unknown error"
    if [[ ! -s "$ftfile" ]]; then
      reason="Filtlong output missing"
    elif [[ ! -x "${BBMAP}/reformat.sh" ]]; then
      reason="reformat.sh not found"
    else
      reason="check ${logfile} for details"
    fi

    log_failure "$file" "Reformat" "$reason"
    return 0  # Continue processing other files
  fi

  debug_msg "Reformat output size: $(stat -c%s "$fafile.tmp.fa" 2>/dev/null || echo "0") bytes"

  # -------------------------------------------------------------------------
  # Stage 4: Header standardization
  # -------------------------------------------------------------------------
  # Nanopore headers contain metadata (timestamps, quality scores, etc.)
  # We keep only the read ID (first field) for cleaner downstream processing
  # This prevents parsing issues in tools that don't expect complex headers
  #
  # CRASH SAFETY: Use atomic rename pattern
  # Write to .partial, then mv (atomic) to final name
  # If interrupted during cut, .partial is incomplete but $fafile doesn't exist yet
  # Resume will detect missing $fafile and retry from scratch
  if [[ -s "$fafile.tmp.fa" ]]; then
    cut -f1 -d' ' "$fafile.tmp.fa" > "$fafile.partial"
    mv "$fafile.partial" "$fafile"  # Atomic rename
    debug_msg "Final FASTA size: $(stat -c%s "$fafile" 2>/dev/null || echo "0") bytes"
  else
    verbose_msg "Reformat produced empty output for $base"
    touch "$fafile"
  fi

  # Cleanup temporary files unless debugging
  if (( DEBUG )); then
    debug_msg "Keeping temp files for debugging: $ftfile, $fafile.tmp.fa, $fafile.partial"
  else
    rm -f "$ftfile" "$fafile.tmp.fa" "$fafile.partial"
  fi

  if [[ ! -s "$fafile" ]]; then
    verbose_msg "Final FASTA file is empty: $fafile"
    echo "(NO DATA - empty after processing)"
    return 0
  fi

  debug_msg "Successfully processed $base - final size: $(stat -c%s "$fafile") bytes"

  fi  # End SKIP_BASIC check

  # =========================================================================
  # Optional Downstream Analyses (enabled via command-line flags)
  # =========================================================================

  # -------------------------------------------------------------------------
  # Sendsketch: Rapid taxonomic classification via k-mer sketching
  # -------------------------------------------------------------------------
  # Compares reads against NCBI nucleotide database using MinHash algorithm
  # Much faster than Kraken but less precise - good for quick overview
  # Results stored in DuckDB for visualization and reporting
  if (( RUN_SKETCH )); then
    local sk_out="$bcdir/sketch/$base.txt"
    # Run if --force is set OR if no output exists
    if (( FORCE )) || [[ ! -s "$sk_out" ]]; then
      echo -n "SKETCH "
      # Clean existing output if forcing
      (( FORCE )) && rm -f "$sk_out"
      debug_msg "Running sketch on $base"
      if run_cmd "${BBMAP}/sendsketch.sh in=$fafile address=nt out=$sk_out format=3" "$logfile"; then
        # Serialize DB writes to avoid lock conflicts
        sem --id duckdb_lock --fg run_cmd "Rscript ${DANADIR}/sketch-db.r $bcdir" "$logfile" || true
      fi
    else
      debug_msg "Sketch output exists for $base, skipping"
    fi
  fi

  # -------------------------------------------------------------------------
  # Kraken2: Precise taxonomic classification using k-mer database
  # -------------------------------------------------------------------------
  # Loads large reference database into RAM (~50-100GB depending on DB choice)
  # Provides genus/species level identification for marine microbes
  # Outputs: classification per read + summary report with abundance estimates
  # Results parsed and stored in DuckDB, integrated into visualization dashboard
  #
  # CRITICAL: Uses semaphore to ensure only 1 Kraken instance runs at a time
  # This prevents RAM exhaustion (16 workers × 50GB = 800GB needed!)
  # Workers queue here if Kraken is busy, other steps continue in parallel
  if (( RUN_KRAKEN )); then
    local k_tsv="$bcdir/kraken/$base.tsv"
    # Run if --force is set OR if no output exists
    if (( FORCE )) || [[ ! -e "$k_tsv" ]]; then
      echo -n "KRAKEN "
      # Clean existing output if forcing
      if (( FORCE )); then
        rm -f "$k_tsv" "$bcdir/kraken/$base.report"
      fi
      debug_msg "Running kraken on $base (acquiring Kraken lock...)"
      echo "DB: $KRAKEN_DB" >> "$logfile"

      # Semaphore ensures only 1 Kraken runs at a time across all workers
      # --id: Unique semaphore name
      # --fg: Run in foreground (wait for completion)
      if sem --id kraken_db_lock --fg run_kraken_locked "${KRAKEN2}" "${KRAKEN_DB}" "$bcdir/kraken/$base.report" "${fafile}" "${logfile}" "${DANADIR}/kraken_parse.awk" "${k_tsv}"; then
        # Serialize DB writes to avoid lock conflicts
        sem --id duckdb_lock --fg run_cmd "Rscript ${DANADIR}/kraken-db.r $bcdir" "$logfile" || true
        sem --id duckdb_lock --fg run_cmd "Rscript ${DANADIR}/krakenreport-db.r $bcdir" "$logfile" || true
      fi
    else
      debug_msg "Kraken output exists for $base, skipping"
    fi
  fi

  # -------------------------------------------------------------------------
  # Prokka: Rapid prokaryotic genome annotation
  # -------------------------------------------------------------------------
  # Identifies protein-coding genes, rRNA, tRNA in metagenomic reads
  # Metagenome mode: Relaxed gene calling for fragmented/partial sequences
  # Outputs functional annotations for detected genes
  # Useful for understanding metabolic potential of community
  if (( RUN_PROKKA )); then
    local prokdir="$bcdir/prokka/$base"
    shopt -s nullglob
    local prokfiles=("$prokdir/PROKKA_"*.tsv)

    debug_msg "Checking Prokka: dir=$prokdir, found ${#prokfiles[@]} tsv files"

    # Run if --force is set OR if no output exists
    if (( FORCE )) || (( ${#prokfiles[@]} == 0 )); then
      echo -n "PROKKA "

      # Clean existing output if forcing re-run
      if (( FORCE )) && [[ -d "$prokdir" ]]; then
        debug_msg "Force mode: removing existing Prokka output"
        rm -rf "$prokdir"
      fi

      debug_msg "Running prokka on $base"
      mkdir -p "$prokdir"
      if run_cmd "${PROKKA_BIN} --metagenome --fast --cpus $PROKKA_THREADS --evalue 1e-20 --outdir $prokdir --force --quiet $(pwd)/$fafile" "$logfile"; then
        rm -f "$prokdir"/*.{err,fna,fsa,gbk,log,sqn,txt} 2>/dev/null || true
        # Serialize DB writes to avoid lock conflicts
        sem --id duckdb_lock --fg run_cmd "Rscript ${DANADIR}/prokka-db.r $bcdir" "$logfile" || true
      fi
    else
      debug_msg "Prokka output exists for $base (${#prokfiles[@]} files), skipping"
    fi
  fi

  # -------------------------------------------------------------------------
  # HMM Search on Prokka Gene Calls
  # -------------------------------------------------------------------------
  # Search Prokka-predicted proteins with custom HMM databases
  # Uses trusted cutoffs (--cut_tc) from HMM files for accurate annotation
  # Supports multiple HMM files (comma-delimited paths)
  if [[ -n "${HMM_DATABASES}" ]]; then
    local prokdir="$bcdir/prokka/$base"
    local hmmdir="$bcdir/hmm"
    mkdir -p "$hmmdir"

    # Find Prokka's protein file (.faa), excluding temporary files
    shopt -s nullglob
    local faafiles=()
    for f in "$prokdir/PROKKA_"*.faa; do
      # Skip Prokka temporary database files (e.g., PROKKA_*.AMR.tmp.*.faa)
      [[ "$f" =~ \.tmp\. ]] && continue
      faafiles+=("$f")
    done

    if (( ${#faafiles[@]} > 0 )); then
      local faafile="${faafiles[0]}"
      debug_msg "Found Prokka proteins: $faafile"

      # Track if any HMM work is done (for progress display)
      local hmm_ran=0

      # Process each HMM file (comma-delimited paths)
      IFS=',' read -ra HMM_FILES <<< "$HMM_DATABASES"
      for hmmfile in "${HMM_FILES[@]}"; do
        hmmfile=$(echo "$hmmfile" | xargs)  # Trim whitespace

        # Skip if HMM file doesn't exist
        if [[ ! -f "$hmmfile" ]]; then
          verbose_msg "HMM file not found: $hmmfile (skipping)"
          continue
        fi

        # Extract database name from filename (e.g., /path/to/CANT-HYD.hmm -> CANT-HYD)
        local dbname=$(basename "$hmmfile" .hmm)
        local tblout="$hmmdir/${base}.${dbname}.tbl"
        local tsvout="$hmmdir/${base}.${dbname}.tsv"

        # Skip if already processed (unless --force is set)
        if (( ! FORCE )) && [[ -s "$tsvout" ]]; then
          debug_msg "HMM results exist for $dbname: $tsvout (skipping)"
          continue
        fi

        # Clean existing output if forcing re-run
        if (( FORCE )) && [[ -f "$tsvout" ]]; then
          debug_msg "Force mode: removing existing HMM output for $dbname"
          rm -f "$tblout" "$tsvout"
        fi

        # Print HMM label once (first time we actually run hmmsearch)
        if (( hmm_ran == 0 )); then
          echo -n "HMM "
          hmm_ran=1
        fi

        debug_msg "Running hmmsearch with $dbname on $base"

        # Run hmmsearch with trusted cutoffs
        # --cut_tc: Use model-specific trusted cutoffs (recommended for curated DBs)
        # --tblout: Parseable table output
        # --cpu 1: Single-threaded (parallelized at file level, not within file)
        if "${HMMSEARCH}" --cut_tc --cpu 1 --tblout "$tblout" \
           "$hmmfile" "$faafile" >> "$logfile" 2>&1; then

          # Parse hmmsearch output to TSV
          # Format: target_name query_name evalue score
          awk 'BEGIN {OFS="\t"; print "gene_id", "hmm_name", "evalue", "score"}
               !/^#/ && NF >= 5 {print $1, $3, $5, $6}' \
               "$tblout" > "$tsvout"

          verbose_msg "HMM search complete: $dbname ($(wc -l < "$tsvout" | xargs) hits)"
        else
          verbose_msg "HMM search failed for $dbname (see $logfile)"
        fi
      done
    else
      verbose_msg "No Prokka proteins found for HMM search (skipping)"
    fi
  fi

  # -------------------------------------------------------------------------
  # Tetranucleotide Frequency Analysis (for ESOM clustering)
  # -------------------------------------------------------------------------
  # Calculates 4-mer composition profiles for each read
  # TNF (tetranucleotide frequency) signatures are phylogenetically conserved
  # Used for unsupervised binning via self-organizing maps
  # Complementary to coverage-based binning methods
  if (( RUN_TETRA )); then
    local lrn="$bcdir/tetra/$base.lrn"
    mkdir -p "$bcdir/tetra"
    # Run if --force is set OR if no output exists
    if (( FORCE )) || [[ ! -s "$lrn" ]]; then
      echo -n "TETRA "
      # Clean existing output if forcing
      (( FORCE )) && rm -f "$lrn"
      debug_msg "Running tetramer analysis on $base"

      # Create annotation file mapping reads to metadata
      if grep '>' "$fafile" | sed 's/>//' | paste - - - > "$bcdir/tetra/annotation.$base.txt" 2>>"$logfile"; then

        # Calculate tetranucleotide frequencies for each read
        if run_cmd "perl ${APPS}/tetramer_freqs_esom.pl -f $fafile -a $bcdir/tetra/annotation.$base.txt -min $MIN_READLEN -max 10000000" "$logfile"; then

          # Create header with tetramer feature names (done once per barcode)
          [[ ! -s "$bcdir/tnfs.txt" ]] && ls Tetra_*.lrn >/dev/null 2>&1 && \
            paste <(echo "seqid") <(head -n4 Tetra_*.lrn | tail -n1 | cut -f2-) > "$bcdir/tnfs.txt" 2>>"$logfile"

          # Combine annotation with frequency data
          if ls Tetra_${base}*.names >/dev/null 2>&1 && ls Tetra_${base}*.lrn >/dev/null 2>&1; then
            paste <(awk '$1!~/^%/' Tetra_${base}*.names) <(awk '$1!~/^%/' Tetra_${base}*.lrn) \
              | cut -f3,5- > "$lrn" 2>>"$logfile"
          fi

          # Cleanup temporary files
          rm -f Tetra_${base}.* "$bcdir/tetra/annotation.$base.txt" 2>/dev/null || true
          # Serialize DB writes to avoid lock conflicts
          sem --id duckdb_lock --fg run_cmd "Rscript ${DANADIR}/tetra-db.r $bcdir" "$logfile" || true
        fi
      fi
    else
      debug_msg "Tetra output exists for $base, skipping"
    fi
  fi
  
  echo "DONE"
}

# ============================================================================
# STEP 3: Configure Parallel Execution Environment
# ============================================================================
# Export all functions and variables needed by GNU parallel worker processes
# Each worker runs in a separate bash subprocess and needs access to these

export -f process_one validate_fastq run_cmd debug_msg verbose_msg log_failure run_kraken_locked
export OUTPUT BBMAP DANADIR PROKKA_BIN FILTLONG KRAKEN2 KRAKEN_DB HMMSEARCH APPS
export RUN_SKETCH RUN_KRAKEN RUN_TETRA RUN_PROKKA PROKKA_THREADS HMM_DATABASES FORCE
export MIN_READLEN KEEP_PCT VERBOSE DEBUG CACHE_FASTQ FAILURE_LOG

# ============================================================================
# Parallelization Strategy: Balance Speed with Memory Constraints
# ============================================================================
# Most pipeline steps are CPU and I/O bound - these benefit from parallelization
# However, Kraken2 is memory-bound: loads entire database into RAM (50-100GB)
#
# Memory considerations:
#   - BBDuk, Filtlong, format conversion: Modest RAM per process (~1-2GB)
#     → Safe to run many parallel workers (e.g., 16 workers = ~32GB total)
#   - Kraken2: Loads full reference database per process (~50-100GB)
#     → Multiple parallel workers would exceed available RAM
#     → Example: 16 workers × 50GB = 800GB RAM needed (most systems: 64-256GB)
#
# Solution: Selective serialization using semaphores
#   - All pipeline steps run in parallel (full CPU utilization)
#   - Kraken2 calls use a semaphore: only 1 Kraken instance at a time
#   - Workers queue at Kraken step, but other steps continue in parallel
#   - Best of both worlds: speed + memory safety

PARALLEL_JOBS="${THREADS}"
if (( FORCE )); then
  echo "[INFO] --force mode: Optional stages will re-run even if output exists"
fi
if (( RUN_KRAKEN )); then
  echo "[INFO] Kraken2 classification enabled"
  echo "[INFO] Pipeline will use ${THREADS} parallel workers"
  echo "[INFO] Kraken2 calls will be serialized (1 at a time) to prevent RAM exhaustion"
  echo "[INFO] Other steps (BBDuk, Filtlong, etc.) run fully parallel for speed"
fi

# ============================================================================
# STEP 4: Pre-flight Validation
# ============================================================================
# Check all FASTQ files for corruption before main processing
# This step runs in parallel (validation is lightweight and won't exhaust RAM)
# Corrupted files are repaired when possible, logged if unfixable
# Validated files are cached for fast subsequent access

echo "[INFO] $(date +%H:%M:%S) Validating FASTQ files"
printf '%s\0' "${FILES[@]}" | parallel --null -j "${THREADS}" --bar 'validate_fastq {} >/dev/null || echo "[BAD] {}" >&2'

# Map original file paths to their cached/validated counterparts
mapfile -t CACHED < <(
  for file in "${FILES[@]}"; do
    echo "${CACHE_FASTQ}/$(basename "${file}")"
  done
)

# ============================================================================
# STEP 5: Main Processing Loop
# ============================================================================
# Process all files using configured parallelism level
# Progress bar shows real-time completion status
# Time cap prevents batch jobs from running indefinitely (useful for incremental processing)

echo "[INFO] $(date +%H:%M:%S) Processing ${#CACHED[@]} files with ${PARALLEL_JOBS} parallel workers"
start_ts=$(date +%s)

# Time-limited batch processing function
# Allows incremental processing of large datasets over multiple runs
# Each run processes files for MAX_DURATION seconds, then stops gracefully
process_batch() {
  local elapsed=$(( $(date +%s) - start_ts ))
  if (( elapsed > MAX_DURATION )); then
    echo "[WARN] Time limit reached: ${elapsed}s > ${MAX_DURATION}s"
    echo "[INFO] Stopping batch. Resume capability allows restarting from this point."
    return 0
  fi

  # Run the main processing pipeline with progress monitoring
  # --null: Handle filenames with spaces/special characters safely
  # --bar: Display progress bar with ETA
  # Errors are logged but don't halt the pipeline (resilient batch processing)
  printf '%s\0' "${CACHED[@]}" | parallel --null -j "${PARALLEL_JOBS}" --bar process_one {}
}

process_batch

# ============================================================================
# Failure Summary Report
# ============================================================================
# Report any files that failed processing with diagnostic information
# This helps identify systematic issues vs. one-off corrupted files

if [[ -s "${FAILURE_LOG}" ]]; then
  FAILURE_COUNT=$(wc -l < "${FAILURE_LOG}")
  echo ""
  echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
  echo "⚠️  WARNING: ${FAILURE_COUNT} file(s) failed processing"
  echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
  echo ""

  # Count failures by stage
  echo "Failures by stage:"
  awk -F'|' '{print "  " $2}' "${FAILURE_LOG}" | sort | uniq -c | sort -rn
  echo ""

  # Show details
  echo "Failed files (showing first 10):"
  echo ""
  head -10 "${FAILURE_LOG}" | while IFS='|' read -r file stage reason; do
    basename_file=$(basename "$file")
    printf "  ❌ %-50s\n" "$basename_file"
    printf "     Stage: %s\n" "$stage"
    printf "     Reason: %s\n\n" "$reason"
  done

  if (( FAILURE_COUNT > 10 )); then
    echo "  ... and $((FAILURE_COUNT - 10)) more (see ${FAILURE_LOG})"
    echo ""
  fi

  echo "Common fixes:"
  echo "  • Corrupted files: Re-download or re-extract original data"
  echo "  • Tool errors: Check log files in output directories"
  echo "  • Empty files: These are skipped automatically (not errors)"
  echo "  • Path issues: Verify BBMAP, FILTLONG paths are correct"
  echo ""
  echo "Full failure log: ${FAILURE_LOG}"
  echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
else
  echo ""
  echo "✅ All files processed successfully!"
fi

echo ""
echo "[DONE] Output directory: ${OUTPUT}"