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
while [[ $# -gt 0 ]]; do
    case "$1" in
        --runtime)  CONTAINER_RUNTIME="$2"; shift 2 ;;
        --image)    CONTAINER_IMAGE="$2"; shift 2 ;;
        --pull)     DO_PULL=true; shift ;;
        --sif)      SIF_PATH="$2"; shift 2 ;;
        *)          NF_ARGS+=("$1"); shift ;;
    esac
done

# ============================================================================
# Parse --outdir and --input from Nextflow args for bind mounts
# ============================================================================
OUTDIR="results"
INPUT_DIR=""
for i in "${!NF_ARGS[@]}"; do
    case "${NF_ARGS[$i]}" in
        --outdir) OUTDIR="${NF_ARGS[$((i+1))]}" ;;
        --input)  INPUT_DIR="${NF_ARGS[$((i+1))]}" ;;
    esac
done
OUTDIR="$(mkdir -p "$OUTDIR" && cd "$OUTDIR" && pwd)"

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

        BIND_ARGS=("-v" "${OUTDIR}:${OUTDIR}")
        [[ -n "$INPUT_DIR" ]] && BIND_ARGS+=("-v" "$(cd "$INPUT_DIR" && pwd):$(cd "$INPUT_DIR" && pwd):ro")

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

        BIND_ARGS="--bind ${OUTDIR}:${OUTDIR}"
        [[ -n "$INPUT_DIR" ]] && BIND_ARGS="${BIND_ARGS},$(cd "$INPUT_DIR" && pwd):$(cd "$INPUT_DIR" && pwd):ro"

        "$RUNTIME_CMD" run \
            $BIND_ARGS \
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
