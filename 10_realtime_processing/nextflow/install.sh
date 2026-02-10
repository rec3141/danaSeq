#!/usr/bin/env bash
set -euo pipefail

# ============================================================================
# Dana Pipeline - Conda Environment Installer
# ============================================================================
#
# Creates isolated conda environments for the pipeline without modifying
# the host's base environment. All envs are prefix-installed under
# ./conda-envs/ (not in ~/.conda/envs/) so they never collide with
# existing environments.
#
# Three environments are needed because of dependency conflicts:
#   dana-bbmap   - BBMap (samtools/libdeflate conflicts with R)
#   dana-prokka  - Prokka (BioPerl pins perl to 5.26, conflicts with others)
#   dana-tools   - Everything else (filtlong, kraken2, hmmer, R/DuckDB, perl)
#
# Each is built from a pinned YAML file in envs/ for reproducibility.
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
ENVS_YAML_DIR="${SCRIPT_DIR}/envs"
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
# Environment definitions
# ============================================================================

# Ordered list of YAML files to install
ENV_YAMLS=(
    bbmap.yml
    prokka.yml
    tools.yml
)

# Which tool binary to check for each environment
declare -A ENV_CHECK=(
    [dana-bbmap]="bbduk.sh"
    [dana-prokka]="prokka"
    [dana-tools]="filtlong"
)

# Additional binaries to verify in dana-tools
TOOLS_CHECK=(filtlong kraken2 hmmsearch gawk perl Rscript)

yaml_to_envname() {
    local yaml="$1"
    echo "dana-${yaml%.yml}"
}

# ============================================================================
# Actions
# ============================================================================

do_install() {
    mkdir -p "${ENV_DIR}"
    local total=${#ENV_YAMLS[@]}
    local count=0
    local failed=0

    echo ""
    echo "Installing ${total} conda environments into: ${ENV_DIR}"
    echo "Using pinned YAML files from: ${ENVS_YAML_DIR}"
    echo ""

    for yaml in "${ENV_YAMLS[@]}"; do
        count=$((count + 1))
        local env_name
        env_name=$(yaml_to_envname "$yaml")
        local env_path="${ENV_DIR}/${env_name}"
        local yaml_path="${ENVS_YAML_DIR}/${yaml}"

        echo "[$count/$total] ${env_name} (from ${yaml})"

        if [[ ! -f "${yaml_path}" ]]; then
            echo "  [ERROR] YAML file not found: ${yaml_path}" >&2
            failed=$((failed + 1))
            continue
        fi

        if [[ -d "${env_path}/conda-meta" ]]; then
            echo "  Already exists, skipping (use --clean to rebuild)"
            continue
        fi

        echo "  Creating environment from ${yaml}..."
        local log_file
        log_file=$(mktemp)
        if ! ${CONDA_CMD} env create -y -p "${env_path}" -f "${yaml_path}" \
            > "${log_file}" 2>&1; then
            echo "  [ERROR] Failed to create ${env_name}" >&2
            sed 's/^/  /' "${log_file}" >&2
            rm -f "${log_file}"
            failed=$((failed + 1))
            continue
        fi
        rm -f "${log_file}"

        echo "  Done"
    done

    # Download tetramer_freqs_esom.pl into the tools environment
    echo ""
    echo "Downloading tetramer_freqs_esom.pl..."
    local tetra_bin="${ENV_DIR}/dana-tools/bin/tetramer_freqs_esom.pl"
    if [[ -d "${ENV_DIR}/dana-tools" ]]; then
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

    # Check each environment's primary binary
    for yaml in "${ENV_YAMLS[@]}"; do
        local env_name
        env_name=$(yaml_to_envname "$yaml")
        local env_path="${ENV_DIR}/${env_name}"
        local check_bin="${ENV_CHECK[$env_name]}"

        printf "  %-20s" "${env_name}"

        if [[ ! -d "${env_path}" ]]; then
            echo "MISSING (not installed)"
            missing=$((missing + 1))
            continue
        fi

        if [[ -x "${env_path}/bin/${check_bin}" ]]; then
            local version
            version=$("${env_path}/bin/${check_bin}" --version 2>&1 | head -1) || version="(installed)"
            echo "OK  ${version}"
            ok=$((ok + 1))
        else
            echo "BROKEN (${check_bin} not found in env)"
            missing=$((missing + 1))
        fi
    done

    # Detailed check of dana-tools binaries
    echo ""
    echo "  dana-tools contents:"
    for bin in "${TOOLS_CHECK[@]}"; do
        printf "    %-20s" "${bin}"
        if [[ -x "${ENV_DIR}/dana-tools/bin/${bin}" ]]; then
            echo "OK"
        else
            echo "MISSING"
            missing=$((missing + 1))
        fi
    done

    # Check tetramer script
    printf "    %-20s" "tetramer_freqs"
    if [[ -x "${ENV_DIR}/dana-tools/bin/tetramer_freqs_esom.pl" ]]; then
        echo "OK"
    else
        echo "MISSING (run install.sh to download)"
        missing=$((missing + 1))
    fi

    echo ""
    echo "Result: ${ok} envs OK, ${missing} issue(s)"

    # Check Nextflow
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
    for yaml in "${ENV_YAMLS[@]}"; do
        local env_name
        env_name=$(yaml_to_envname "$yaml")
        if [[ -d "${ENV_DIR}/${env_name}" ]]; then
            echo "  ${env_name}"
            total=$((total + 1))
        fi
    done

    if (( total == 0 )); then
        echo "No pipeline environments found."
        return 0
    fi

    echo ""
    read -rp "Remove ${total} environments? [y/N] " confirm
    if [[ "${confirm}" =~ ^[Yy] ]]; then
        for yaml in "${ENV_YAMLS[@]}"; do
            local env_name
            env_name=$(yaml_to_envname "$yaml")
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
