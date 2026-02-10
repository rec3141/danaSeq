#!/usr/bin/env bash
set -euo pipefail

# ============================================================================
# Dana Pipeline - Conda Environment Installer
# ============================================================================
#
# Creates isolated conda environments for each pipeline tool without
# modifying the host's base environment. All envs are stored under
# ./conda-envs/ (or a custom path via --prefix).
#
# Usage:
#   ./install.sh              # Install all environments
#   ./install.sh --prefix /custom/path
#   ./install.sh --check      # Verify existing installations
#   ./install.sh --clean      # Remove all pipeline environments
#
# Requirements: conda or mamba must be on $PATH
# ============================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ENV_DIR="${SCRIPT_DIR}/conda-envs"
TETRAMER_URL="https://raw.githubusercontent.com/tetramerFreqs/Binning/master/tetramer_freqs_esom.pl"

# Detect conda/mamba
if command -v mamba &>/dev/null; then
    CONDA_CMD="mamba"
elif command -v conda &>/dev/null; then
    CONDA_CMD="conda"
else
    echo "[ERROR] Neither conda nor mamba found on PATH." >&2
    echo "Install Miniforge: https://github.com/conda-forge/miniforge" >&2
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
        -h|--help)
            sed -n '/^# Usage:/,/^# ====/p' "$0" | head -n -1 | sed 's/^# //'
            exit 0 ;;
        *) echo "[ERROR] Unknown argument: $1" >&2; exit 1 ;;
    esac
done

# ============================================================================
# Environment definitions: name -> conda packages
# ============================================================================
# Each maps to one or more Nextflow processes. Separate environments prevent
# dependency conflicts between tools with incompatible requirements.

declare -A ENVS=(
    [dana-bbmap]="bioconda::bbmap"
    [dana-filtlong]="bioconda::filtlong"
    [dana-kraken2]="bioconda::kraken2 conda-forge::gawk"
    [dana-prokka]="bioconda::prokka"
    [dana-hmmer]="bioconda::hmmer"
    [dana-tetramer]="conda-forge::perl conda-forge::curl"
    [dana-r-duckdb]="conda-forge::r-base conda-forge::r-dbi conda-forge::r-duckdb conda-forge::r-readr"
)

# Which tool binary to check for each environment
declare -A ENV_CHECK=(
    [dana-bbmap]="bbduk.sh"
    [dana-filtlong]="filtlong"
    [dana-kraken2]="kraken2"
    [dana-prokka]="prokka"
    [dana-hmmer]="hmmsearch"
    [dana-tetramer]="perl"
    [dana-r-duckdb]="Rscript"
)

# ============================================================================
# Actions
# ============================================================================

do_install() {
    mkdir -p "${ENV_DIR}"
    local total=${#ENVS[@]}
    local count=0
    local failed=0

    echo ""
    echo "Installing ${total} conda environments into: ${ENV_DIR}"
    echo ""

    for env_name in $(printf '%s\n' "${!ENVS[@]}" | sort); do
        ((count++))
        local packages="${ENVS[$env_name]}"
        local env_path="${ENV_DIR}/${env_name}"

        echo "[$count/$total] ${env_name}"

        if [[ -d "${env_path}" ]] && [[ -f "${env_path}/bin/activate" ]]; then
            echo "  Already exists, skipping (use --clean to rebuild)"
            continue
        fi

        echo "  Creating: ${CONDA_CMD} create -y -p ${env_path} ${packages}"
        if ! ${CONDA_CMD} create -y -p "${env_path}" -c conda-forge -c bioconda ${packages} 2>&1 \
            | while IFS= read -r line; do echo "  $line"; done; then
            echo "  [ERROR] Failed to create ${env_name}" >&2
            ((failed++))
            continue
        fi

        echo "  Done"
    done

    # Download tetramer_freqs_esom.pl into the tetramer environment
    echo ""
    echo "Downloading tetramer_freqs_esom.pl..."
    local tetra_bin="${ENV_DIR}/dana-tetramer/bin/tetramer_freqs_esom.pl"
    if [[ -d "${ENV_DIR}/dana-tetramer" ]]; then
        curl -fsSL "${TETRAMER_URL}" -o "${tetra_bin}"
        chmod +x "${tetra_bin}"
        echo "  Installed to: ${tetra_bin}"
    fi

    echo ""
    if (( failed > 0 )); then
        echo "[WARNING] ${failed} environment(s) failed to install"
        echo "Run with --check to see which ones need attention"
        return 1
    else
        echo "[SUCCESS] All ${total} environments installed to: ${ENV_DIR}"
        echo ""
        echo "Nextflow will auto-detect these environments via conda.cacheDir."
        echo "Run the pipeline with:"
        echo "  nextflow run main.nf --input /path/to/data -resume"
    fi
}

do_check() {
    echo ""
    echo "Checking conda environments in: ${ENV_DIR}"
    echo ""

    local ok=0
    local missing=0

    for env_name in $(printf '%s\n' "${!ENVS[@]}" | sort); do
        local env_path="${ENV_DIR}/${env_name}"
        local check_bin="${ENV_CHECK[$env_name]}"

        printf "  %-20s" "${env_name}"

        if [[ ! -d "${env_path}" ]]; then
            echo "MISSING (not installed)"
            ((missing++))
            continue
        fi

        if [[ -x "${env_path}/bin/${check_bin}" ]]; then
            local version
            version=$("${env_path}/bin/${check_bin}" --version 2>&1 | head -1) || version="(installed)"
            echo "OK  ${version}"
            ((ok++))
        else
            echo "BROKEN (${check_bin} not found in env)"
            ((missing++))
        fi
    done

    # Special check for tetramer script
    printf "  %-20s" "tetramer_freqs"
    if [[ -x "${ENV_DIR}/dana-tetramer/bin/tetramer_freqs_esom.pl" ]]; then
        echo "OK"
        ((ok++))
    else
        echo "MISSING (run install.sh to download)"
        ((missing++))
    fi

    echo ""
    echo "Result: ${ok} OK, ${missing} missing"

    # Also check Nextflow
    echo ""
    printf "  %-20s" "nextflow"
    if command -v nextflow &>/dev/null; then
        echo "OK  $(nextflow -version 2>&1 | grep 'version' | head -1)"
    else
        echo "MISSING (install: curl -s https://get.nextflow.io | bash)"
    fi

    (( missing == 0 )) && return 0 || return 1
}

do_clean() {
    if [[ ! -d "${ENV_DIR}" ]]; then
        echo "Nothing to clean: ${ENV_DIR} does not exist"
        return 0
    fi

    echo "This will remove all conda environments in: ${ENV_DIR}"
    echo ""

    local total=0
    for env_name in "${!ENVS[@]}"; do
        if [[ -d "${ENV_DIR}/${env_name}" ]]; then
            echo "  ${env_name}"
            ((total++))
        fi
    done

    if (( total == 0 )); then
        echo "No pipeline environments found."
        return 0
    fi

    echo ""
    read -rp "Remove ${total} environments? [y/N] " confirm
    if [[ "${confirm}" =~ ^[Yy] ]]; then
        for env_name in "${!ENVS[@]}"; do
            if [[ -d "${ENV_DIR}/${env_name}" ]]; then
                echo "  Removing ${env_name}..."
                rm -rf "${ENV_DIR}/${env_name}"
            fi
        done
        echo "[SUCCESS] All pipeline environments removed"
    else
        echo "Cancelled."
    fi
}

# Run the requested action
case "${ACTION}" in
    install) do_install ;;
    check)   do_check ;;
    clean)   do_clean ;;
esac
