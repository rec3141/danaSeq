#!/usr/bin/env bash
#
# run-illumina-amplicon.sh — Run the dānaSeq Illumina amplicon pipeline
#
# Auto-detects container runtime (Apptainer > Singularity > Docker) or
# falls back to local conda execution. Handles bind mounts, CA certificates,
# and user permissions transparently.
#
# Usage:
#   ./run-illumina-amplicon.sh --input /path/to/reads [pipeline options]
#
# Container options:
#   --runtime docker|apptainer|singularity|conda   Force a specific runtime
#   --image IMAGE       Override container image [default: ghcr.io/rec3141/danaseq-illumina-amplicon:latest]
#   --pull              Pull/update the container image before running
#   --sif PATH          Use a pre-built .sif file (Apptainer/Singularity)
#                       [default: .danaseq-illumina-amplicon.sif beside this
#                        script; falls back to pulling into --outdir]
#
# All other arguments are passed through to Nextflow.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
CONTAINER_IMAGE="ghcr.io/rec3141/danaseq-illumina-amplicon:latest"
CONTAINER_RUNTIME="auto"
SIF_PATH=""
DO_PULL=false

# ============================================================================
# Parse wrapper-specific arguments (strip them before passing to Nextflow)
# ============================================================================
NF_ARGS=()
EXTRA_RO_BINDS=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        --runtime)   CONTAINER_RUNTIME="$2"; shift 2 ;;
        --runtime=*) CONTAINER_RUNTIME="${1#--runtime=}"; shift ;;
        --image)     CONTAINER_IMAGE="$2"; shift 2 ;;
        --image=*)   CONTAINER_IMAGE="${1#--image=}"; shift ;;
        --pull)      DO_PULL=true; shift ;;
        --sif)       SIF_PATH="$2"; shift 2 ;;
        --sif=*)     SIF_PATH="${1#--sif=}"; shift ;;
        --bind)      EXTRA_RO_BINDS+=("$(realpath -m "$2")"); shift 2 ;;
        --bind=*)    EXTRA_RO_BINDS+=("$(realpath -m "${1#--bind=}")"); shift ;;
        *)           NF_ARGS+=("$1"); shift ;;
    esac
done

# ============================================================================
# Parse host path arguments from Nextflow args for container bind mounts
# ============================================================================
OUTDIR="results"
RW_PATHS=()
RO_PATHS=("${EXTRA_RO_BINDS[@]}")
SAMPLESHEET_PATH=""

RO_PATH_FLAGS=" --input --samplesheet --metadata --forward_bcs --reverse_bcs --primers_fwd --primers_rev --earthdna_db --dada2_db --dada2_species --dada2_species_db --sintax_db --globe_db "
RW_PATH_FLAGS=" --store_dir -w -work-dir --work-dir "

extract_ref_db_paths() {
    local spec="$1" entry rest db_path
    IFS=';' read -ra entries <<< "$spec"
    for entry in "${entries[@]}"; do
        [[ "$entry" == *:* ]] || continue
        rest="${entry#*:}"
        db_path="${rest%%:*}"
        [[ -n "$db_path" ]] && RO_PATHS+=("$db_path")
    done
}

for (( i=0; i<${#NF_ARGS[@]}; i++ )); do
    arg="${NF_ARGS[$i]}"
    case "$arg" in
        --outdir)
            [[ $((i+1)) -lt ${#NF_ARGS[@]} ]] && OUTDIR="${NF_ARGS[$((i+1))]}" ;;
        --outdir=*)
            OUTDIR="${arg#--outdir=}" ;;
        --ref_databases)
            [[ $((i+1)) -lt ${#NF_ARGS[@]} ]] && extract_ref_db_paths "${NF_ARGS[$((i+1))]}" ;;
        --ref_databases=*)
            extract_ref_db_paths "${arg#--ref_databases=}" ;;
        --*=*|-*=*)
            flag="${arg%%=*}"
            val="${arg#*=}"
            if [[ "$RO_PATH_FLAGS" == *" ${flag} "* && -n "$val" ]]; then
                RO_PATHS+=("$val")
                [[ "$flag" == "--samplesheet" ]] && SAMPLESHEET_PATH="$val"
            elif [[ "$RW_PATH_FLAGS" == *" ${flag} "* && -n "$val" ]]; then
                RW_PATHS+=("$val")
            fi ;;
        *)
            if [[ "$RO_PATH_FLAGS" == *" ${arg} "* && $((i+1)) -lt ${#NF_ARGS[@]} ]]; then
                val="${NF_ARGS[$((i+1))]}"
                RO_PATHS+=("$val")
                [[ "$arg" == "--samplesheet" ]] && SAMPLESHEET_PATH="$val"
            elif [[ "$RW_PATH_FLAGS" == *" ${arg} "* && $((i+1)) -lt ${#NF_ARGS[@]} ]]; then
                RW_PATHS+=("${NF_ARGS[$((i+1))]}")
            fi ;;
    esac
done

OUTDIR="$(mkdir -p "$OUTDIR" && cd "$OUTDIR" && pwd)"
RW_PATHS=("$OUTDIR" "${RW_PATHS[@]}")

# Ensure writable host directories (e.g. --store_dir, -w) exist before container launch
RESOLVED_RW=()
for p in "${RW_PATHS[@]}"; do
    [[ -n "$p" ]] || continue
    mkdir -p "$p"
    RESOLVED_RW+=("$(realpath -m "$p")")
done
if [[ ${#RESOLVED_RW[@]} -gt 0 ]]; then
    mapfile -t RESOLVED_RW < <(printf '%s\n' "${RESOLVED_RW[@]}" | awk 'NF' | sort -u)
fi

# Include read parent directories referenced inside --samplesheet (RFC-4180 CSV by column name)
if [[ -n "$SAMPLESHEET_PATH" && -f "$SAMPLESHEET_PATH" ]]; then
    PY_CSV=""
    for pycand in "$(command -v python3 2>/dev/null || true)" /usr/bin/python3 /opt/miniforge3/bin/python3; do
        if [[ -n "$pycand" && -x "$pycand" ]]; then PY_CSV="$pycand"; break; fi
    done
    if [[ -n "$PY_CSV" ]]; then
        mapfile -t SHEET_READ_DIRS < <(
            "$PY_CSV" - "$SAMPLESHEET_PATH" <<'PYEOF'
import csv, os, sys

sheet = sys.argv[1]
sheet_dir = os.path.dirname(os.path.abspath(sheet))
with open(sheet, newline="") as fh:
    lines = [ln for ln in fh if ln.strip() and not ln.lstrip().startswith("#")]
for row in csv.DictReader(lines):
    if not row:
        continue
    for key in ("r1_path", "r2_path", "r1", "r2"):
        val = (row.get(key) or "").strip()
        if not val:
            continue
        for cand in (val, os.path.join(sheet_dir, val)):
            if os.path.lexists(cand):
                print(os.path.dirname(os.path.abspath(cand)))
                print(os.path.dirname(os.path.realpath(cand)))
                break
PYEOF
        )
        for d in "${SHEET_READ_DIRS[@]}"; do
            [[ -n "$d" ]] && RO_PATHS+=("$d")
        done
    fi
fi

# Expand symlinks and file parent directories, deduplicate, and drop paths covered by a bound parent
ALL_RO_CANDIDATES=()
for p in "${RO_PATHS[@]}"; do
    [[ -n "$p" ]] || continue
    if [[ -L "$p" || -f "$p" ]]; then
        ALL_RO_CANDIDATES+=("$(dirname "$(realpath -m "$p")")")
        [[ -L "$p" ]] && ALL_RO_CANDIDATES+=("$(dirname "$p")")
    elif [[ -d "$p" ]]; then
        ALL_RO_CANDIDATES+=("$(cd "$p" && pwd)")
        r="$(realpath -m "$p")"
        [[ "$r" != "$p" && -d "$r" ]] && ALL_RO_CANDIDATES+=("$r")
    fi
done
if [[ ${#ALL_RO_CANDIDATES[@]} -gt 0 ]]; then
    mapfile -t ALL_RO_CANDIDATES < <(printf '%s\n' "${ALL_RO_CANDIDATES[@]}" | awk 'NF' | sort -u)
fi

RESOLVED_RO=()
for p in "${ALL_RO_CANDIDATES[@]}"; do
    [[ -d "$p" || -f "$p" ]] || continue
    covered=false
    for rw in "${RESOLVED_RW[@]}"; do
        [[ "$p" == "$rw" || "$p" == "$rw"/* ]] && { covered=true; break; }
    done
    if [[ "$covered" == false ]]; then
        for q in "${ALL_RO_CANDIDATES[@]}"; do
            [[ "$q" != "$p" && "$p" == "$q"/* ]] && { covered=true; break; }
        done
    fi
    [[ "$covered" == false ]] && RESOLVED_RO+=("$p")
done

# ============================================================================
# Auto-detect container runtime
# ============================================================================
if [[ "$CONTAINER_RUNTIME" == "auto" ]]; then
    if command -v apptainer &>/dev/null; then
        CONTAINER_RUNTIME=apptainer
    elif command -v singularity &>/dev/null; then
        CONTAINER_RUNTIME=singularity
    elif command -v docker &>/dev/null; then
        CONTAINER_RUNTIME=docker
    else
        CONTAINER_RUNTIME=conda
    fi
fi
echo "[INFO] Container runtime: $CONTAINER_RUNTIME"

# ============================================================================
# Find CA certificate bundle for HTTPS in containers (HPC often needs this)
# ============================================================================
container_ca=""
for ca in /etc/ssl/certs/ca-certificates.crt /etc/pki/tls/certs/ca-bundle.crt /etc/ssl/cert.pem; do
    if [[ -f "$ca" ]]; then
        container_ca="$ca"
        break
    fi
done

# ============================================================================
# Run with the detected runtime
# ============================================================================

case "$CONTAINER_RUNTIME" in

    docker)
        if $DO_PULL; then
            echo "[INFO] Pulling $CONTAINER_IMAGE"
            docker pull "$CONTAINER_IMAGE"
        fi

        BIND_ARGS=()
        for rw in "${RESOLVED_RW[@]}"; do
            BIND_ARGS+=("-v" "${rw}:${rw}")
        done
        for ro in "${RESOLVED_RO[@]}"; do
            BIND_ARGS+=("-v" "${ro}:${ro}:ro")
        done

        docker run --rm \
            --user "$(id -u):$(id -g)" \
            "${BIND_ARGS[@]}" \
            -w "$(pwd)" \
            ${container_ca:+--env "REQUESTS_CA_BUNDLE=${container_ca}"} \
            ${container_ca:+--env "SSL_CERT_FILE=${container_ca}"} \
            "$CONTAINER_IMAGE" \
            run /pipeline/main.nf \
            "${NF_ARGS[@]}" -resume
        ;;

    apptainer|singularity)
        RUNTIME_CMD="$CONTAINER_RUNTIME"

        # Resolve the .sif the same way every other danaSeq wrapper does: the
        # dot-file beside this script, which is what the rebuild scripts write.
        # Only when that is absent (a fresh checkout, or a read-only/over-quota
        # /project) do we fall back to pulling into the outdir — otherwise each
        # new outdir re-downloaded its own ~900 MB copy.
        if [[ -z "$SIF_PATH" ]]; then
            if [[ -f "${SCRIPT_DIR}/.danaseq-illumina-amplicon.sif" ]]; then
                SIF_PATH="${SCRIPT_DIR}/.danaseq-illumina-amplicon.sif"
            else
                SIF_PATH="${OUTDIR}/danaseq-illumina-amplicon.sif"
            fi
        fi

        if $DO_PULL || [[ ! -f "$SIF_PATH" ]]; then
            echo "[INFO] Pulling $CONTAINER_IMAGE -> $SIF_PATH"
            "$RUNTIME_CMD" pull --force "$SIF_PATH" "docker://${CONTAINER_IMAGE}"
        fi

        APPTAINER_BINDS=()
        for rw in "${RESOLVED_RW[@]}"; do
            APPTAINER_BINDS+=("--bind" "${rw}:${rw}")
        done
        for ro in "${RESOLVED_RO[@]}"; do
            APPTAINER_BINDS+=("--bind" "${ro}:${ro}:ro")
        done

        "$RUNTIME_CMD" run \
            "${APPTAINER_BINDS[@]}" \
            --pwd "$(pwd)" \
            ${container_ca:+--env "REQUESTS_CA_BUNDLE=${container_ca}"} \
            ${container_ca:+--env "SSL_CERT_FILE=${container_ca}"} \
            ${container_ca:+--env "CURL_CA_BUNDLE=${container_ca}"} \
            "$SIF_PATH" \
            run /pipeline/main.nf \
            "${NF_ARGS[@]}" -resume
        ;;

    conda)
        echo "[INFO] Running with local conda environments"
        echo "[INFO] Environments will be created in ${SCRIPT_DIR}/conda-envs/ on first run"

        # Use Java 17 if available and default is too new
        for jdir in /usr/lib/jvm/java-17-openjdk-amd64 /usr/lib/jvm/java-17-openjdk; do
            if [[ -d "$jdir" ]]; then
                export JAVA_HOME="$jdir"
                export PATH="$JAVA_HOME/bin:$PATH"
                break
            fi
        done

        nextflow run "${SCRIPT_DIR}/main.nf" \
            "${NF_ARGS[@]}" -resume
        ;;

    *)
        echo "[ERROR] Unknown runtime: $CONTAINER_RUNTIME" >&2
        exit 1
        ;;
esac
