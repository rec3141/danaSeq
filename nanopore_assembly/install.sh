#!/usr/bin/env bash
set -euo pipefail

# ============================================================================
# Dana Nanopore Assembly Pipeline - Conda Environment Installer
# ============================================================================
#
# Creates the conda environment for the nanopore assembly pipeline.
#
# Environment:
#   dana-mag-assembly  - Flye, metaMDBG, myloasm, BBMap, minimap2, samtools,
#                        CoverM, filtlong, Nextflow, OpenJDK
#
# Usage:
#   ./install.sh              # Install environment
#   ./install.sh --check      # Verify existing installation
#   ./install.sh --clean      # Remove environment
#
# Requirements: conda or mamba must be on $PATH
# ============================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ENV_DIR="${SCRIPT_DIR}/conda-envs"
MERGED_YAML_DIR="${SCRIPT_DIR}/envs/merged"

# Detect conda/mamba
if command -v mamba &>/dev/null; then
    CONDA_CMD="mamba"
elif command -v conda &>/dev/null; then
    CONDA_CMD="conda"
else
    echo "[ERROR] Neither conda nor mamba found on PATH." >&2
    exit 1
fi

echo "[INFO] Using: $CONDA_CMD ($(${CONDA_CMD} --version 2>&1))"

# Parse arguments
ACTION="install"
while (( $# )); do
    case "$1" in
        --prefix)  ENV_DIR="$2"; shift 2 ;;
        --check)   ACTION="check"; shift ;;
        --clean)   ACTION="clean"; shift ;;
        -h|--help) sed -n '/^# Usage:/,/^# ====/p' "$0" | head -n -1 | sed 's/^# //'; exit 0 ;;
        *) echo "[ERROR] Unknown argument: $1" >&2; exit 1 ;;
    esac
done

ENV_NAME="dana-mag-assembly"
ENV_PATH="${ENV_DIR}/${ENV_NAME}"
YAML_PATH="${MERGED_YAML_DIR}/assembly.yml"
CHECK_BINS=("flye" "filtlong" "minimap2" "samtools" "coverm" "bbduk.sh" "nextflow" "java")

do_install() {
    mkdir -p "${ENV_DIR}"

    if [[ ! -f "${YAML_PATH}" ]]; then
        echo "[ERROR] YAML file not found: ${YAML_PATH}" >&2
        exit 1
    fi

    if [[ -d "${ENV_PATH}/conda-meta" ]]; then
        echo "[INFO] ${ENV_NAME} already exists, skipping (use --clean to rebuild)"
        return 0
    fi

    echo "[INFO] Creating ${ENV_NAME} from assembly.yml..."
    ${CONDA_CMD} env create -y -p "${ENV_PATH}" -f "${YAML_PATH}"

    # Symlink conda into env for Nextflow activation
    local conda_bin
    conda_bin=$(which conda 2>/dev/null || echo "")
    if [[ -n "${conda_bin}" && ! -e "${ENV_PATH}/bin/conda" ]]; then
        ln -sf "${conda_bin}" "${ENV_PATH}/bin/conda"
        echo "[INFO] Symlinked conda into assembly env for Nextflow activation"
    fi

    # Compile C binaries
    compile_binaries

    echo "[SUCCESS] Environment installed to: ${ENV_PATH}"
}

compile_binaries() {
    local bin_dir="${SCRIPT_DIR}/bin"
    local cc="${CC:-gcc}"

    echo "[INFO] Compiling C binaries..."

    # fastq_filter (C++)
    if [[ -f "${bin_dir}/fastq_filter.cpp" ]]; then
        echo "  Compiling fastq_filter..."
        g++ -O2 -o "${bin_dir}/fastq_filter" "${bin_dir}/fastq_filter.cpp" -lz -lpthread
        chmod +x "${bin_dir}/fastq_filter"
        echo "  [OK] fastq_filter"
    fi

    # tetramer_freqs (C)
    if [[ -f "${bin_dir}/tetramer_freqs.c" ]]; then
        echo "  Compiling tetramer_freqs..."
        ${cc} -O2 -o "${bin_dir}/tetramer_freqs" "${bin_dir}/tetramer_freqs.c" -lz -lm
        chmod +x "${bin_dir}/tetramer_freqs"
        echo "  [OK] tetramer_freqs"
    fi

    echo "[INFO] C binaries compiled"
}

do_check() {
    local failed=0
    echo "[INFO] Checking ${ENV_NAME}..."

    if [[ ! -d "${ENV_PATH}" ]]; then
        echo "  [MISSING] Environment not found at ${ENV_PATH}"
        return 1
    fi

    for bin in "${CHECK_BINS[@]}"; do
        if [[ -x "${ENV_PATH}/bin/${bin}" ]]; then
            echo "  [OK] ${bin}"
        else
            echo "  [MISSING] ${bin}"
            failed=$((failed + 1))
        fi
    done

    # Check compiled C binaries in bin/
    local bin_dir="${SCRIPT_DIR}/bin"
    for cbin in fastq_filter tetramer_freqs; do
        if [[ -x "${bin_dir}/${cbin}" ]]; then
            echo "  [OK] bin/${cbin} (compiled)"
        else
            echo "  [MISSING] bin/${cbin} — run ./install.sh to compile"
            failed=$((failed + 1))
        fi
    done

    if (( failed > 0 )); then
        echo "[WARNING] ${failed} binary(ies) missing"
        return 1
    else
        echo "[OK] All binaries present"
    fi
}

do_clean() {
    if [[ -d "${ENV_PATH}" ]]; then
        echo "[INFO] Removing ${ENV_NAME}..."
        rm -rf "${ENV_PATH}"
        echo "[INFO] Removed"
    else
        echo "[INFO] ${ENV_NAME} not found, nothing to clean"
    fi
}

case "$ACTION" in
    install) do_install ;;
    check)   do_check ;;
    clean)   do_clean ;;
esac
