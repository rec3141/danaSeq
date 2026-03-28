#!/usr/bin/env bash
set -euo pipefail

# ============================================================================
# Dana MAG Assembly Pipeline - Conda Environment Installer
# ============================================================================
#
# Creates conda environments for the MAG assembly pipeline.
# All envs are prefix-installed under ./conda-envs/.
#
# Merged environments (9 envs, matching the slim Docker image):
#   dana-mag-assembly  - Flye, metaMDBG, myloasm, BBMap, minimap2, samtools, filtlong
#   dana-mag-binning   - MetaBAT2, MaxBin2, SemiBin2, COMEBin, LorBin, DAS_Tool
#   dana-mag-quality   - geNomad, CheckV, CheckM2, Tiara, Whokaryote
#   dana-mag-annotate  - Bakta, eggNOG-mapper, MarFERReT/DIAMOND
#   dana-mag-classify  - Kaiju, Kraken2, barrnap/vsearch, KofamScan, MetaEuk
#   dana-mag-genomic   - IntegronFinder, IslandPath, MacSyFinder, DefenseFinder, dbCAN
#   dana-mag-pathviz   - MinPath, KEGG-Decoder, viz dashboard (Node.js + Python)
#   dana-mag-strain    - Strainy, Floria, InStrain
#   dana-mag-derep     - galah, skani, sourmash
#
# Standalone environments (too large or conflicting to merge):
#   dana-mag-prokka    - Prokka gene annotation
#   dana-mag-gtdbtk    - GTDB-Tk phylogenetic classification
#   dana-mag-vamb      - VAMB variational autoencoder binning
#   dana-mag-binette   - Binette consensus refinement
#   dana-mag-antismash - antiSMASH secondary metabolite prediction
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
MERGED_YAML_DIR="${SCRIPT_DIR}/envs/merged"
STANDALONE_YAML_DIR="${SCRIPT_DIR}/envs"

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

# Merged environments — installed from envs/merged/*.yml
MERGED_YAMLS=(
    assembly.yml
    binning.yml
    quality.yml
    annotate.yml
    classify.yml
    genomic.yml
    pathviz.yml
    strain.yml
    derep.yml
)

# Standalone environments — installed from envs/*.yml
STANDALONE_YAMLS=(
    prokka.yml
    gtdbtk.yml
    vamb.yml
    binette.yml
    antismash.yml
)

# Which tool binary to check for each environment
declare -A ENV_CHECK=(
    # Merged envs
    [dana-mag-assembly]="flye"
    [dana-mag-binning]="metabat2"
    [dana-mag-quality]="genomad"
    [dana-mag-annotate]="bakta"
    [dana-mag-classify]="kaiju"
    [dana-mag-genomic]="integron_finder"
    [dana-mag-pathviz]="KEGG-decoder"
    [dana-mag-strain]="inStrain"
    [dana-mag-derep]="galah"
    # Standalone envs
    [dana-mag-prokka]="prokka"
    [dana-mag-gtdbtk]="gtdbtk"
    [dana-mag-vamb]="vamb"
    [dana-mag-binette]="binette"
    [dana-mag-antismash]="antismash"
)

# Additional binaries to verify per env
declare -A ENV_EXTRAS=(
    [dana-mag-assembly]="filtlong minimap2 samtools bbduk.sh metaMDBG myloasm nextflow java"
    [dana-mag-binning]="SemiBin2 LorBin run_MaxBin.pl DAS_Tool jgi_summarize_bam_contig_depths"
    [dana-mag-quality]="checkm2 checkv tiara"
    [dana-mag-annotate]="emapper.py diamond"
    [dana-mag-classify]="kraken2 barrnap exec_annotation metaeuk"
    [dana-mag-genomic]="hmmscan macsyfinder defense-finder run_dbcan"
    [dana-mag-pathviz]="node npm python3"
    [dana-mag-strain]="strainy floria"
    [dana-mag-derep]="skani sourmash"
)

# ============================================================================
# Actions
# ============================================================================

do_install() {
    mkdir -p "${ENV_DIR}"
    local total=$(( ${#MERGED_YAMLS[@]} + ${#STANDALONE_YAMLS[@]} ))
    local count=0
    local failed=0

    echo ""
    echo "Installing ${total} conda environments into: ${ENV_DIR}"
    echo "  Merged envs from: ${MERGED_YAML_DIR}"
    echo "  Standalone envs from: ${STANDALONE_YAML_DIR}"
    echo ""

    # --- Install merged environments ---
    for yaml in "${MERGED_YAMLS[@]}"; do
        count=$((count + 1))
        local env_name="dana-mag-${yaml%.yml}"
        local env_path="${ENV_DIR}/${env_name}"
        local yaml_path="${MERGED_YAML_DIR}/${yaml}"

        echo "[$count/$total] ${env_name} (merged, from ${yaml})"

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

        # Post-install hooks for merged envs
        post_install_merged "${env_name}" "${env_path}"

        echo "  Done"
    done

    # --- Install standalone environments ---
    for yaml in "${STANDALONE_YAMLS[@]}"; do
        count=$((count + 1))
        local env_name="dana-mag-${yaml%.yml}"
        local env_path="${ENV_DIR}/${env_name}"
        local yaml_path="${STANDALONE_YAML_DIR}/${yaml}"

        echo "[$count/$total] ${env_name} (standalone, from ${yaml})"

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

        # Post-install hooks for standalone envs
        post_install_standalone "${env_name}" "${env_path}"

        echo "  Done"
    done

    echo ""
    if (( failed > 0 )); then
        echo "[WARNING] ${failed} environment(s) failed to install"
        echo "Run with --check to see which ones need attention"
        return 1
    else
        echo "[SUCCESS] All ${total} environments installed to: ${ENV_DIR}"
        echo ""
        echo "Run the pipeline:"
        echo "  nextflow run main.nf --input /path/to/reads -resume"
    fi
}

post_install_merged() {
    local env_name="$1"
    local env_path="$2"

    # Ensure conda is on PATH for Nextflow env activation.
    # Nextflow's .command.run calls `conda info --json` to activate prefix envs.
    if [[ "${env_name}" == "dana-mag-assembly" ]]; then
        local conda_bin
        conda_bin=$(which conda 2>/dev/null || echo "")
        if [[ -n "${conda_bin}" && ! -e "${env_path}/bin/conda" ]]; then
            ln -sf "${conda_bin}" "${env_path}/bin/conda"
            echo "  Symlinked conda into assembly env for Nextflow activation"
        fi
    fi

    # Install pip packages that can't be in the conda YAML (build isolation issues)
    if [[ "${env_name}" == "dana-mag-binning" ]]; then
        echo "  Installing LorBin via pip..."
        "${env_path}/bin/pip" install --no-deps \
            'lorbin @ git+https://github.com/rec3141/LorBin.git' hnswlib \
            > /dev/null 2>&1
        echo "  LorBin installed"

        # Clone COMEBin fork and wire up bin/ scripts
        echo "  Installing COMEBin from fork..."
        local comebin_repo="https://github.com/rec3141/COMEBin.git"
        local comebin_branch="codex/recent-dependencies-fix"
        local comebin_dir="${env_path}/share/COMEBin"
        git clone --depth 1 -b "${comebin_branch}" "${comebin_repo}" "${comebin_dir}" \
            > /dev/null 2>&1
        cat > "${env_path}/bin/run_comebin.sh" <<'WRAPPER'
#!/usr/bin/env bash
# Resolve path args (-a/-o/-p) in caller's CWD before cd into source tree
COMEBIN_ROOT="$(dirname "$(dirname "$(readlink -f "$0")")")/share/COMEBin/COMEBin"
ORIG_CWD="$(pwd)"
RESOLVED_ARGS=()
while (( $# )); do
    case "$1" in
        -a|-o|-p) RESOLVED_ARGS+=("$1"); shift
                  RESOLVED_ARGS+=("$(cd "${ORIG_CWD}" && realpath "$1")"); shift ;;
        *) RESOLVED_ARGS+=("$1"); shift ;;
    esac
done
cd "${COMEBIN_ROOT}" && exec bash run_comebin.sh "${RESOLVED_ARGS[@]}"
WRAPPER
        chmod +x "${env_path}/bin/run_comebin.sh"
        if [[ -f "${comebin_dir}/COMEBin/scripts/gen_cov_file.sh" ]]; then
            ln -sf "${comebin_dir}/COMEBin/scripts/gen_cov_file.sh" "${env_path}/bin/gen_cov_file.sh"
            chmod +x "${env_path}/bin/gen_cov_file.sh"
        fi
        echo "  COMEBin installed from ${comebin_repo}@${comebin_branch}"
    fi

    if [[ "${env_name}" == "dana-mag-quality" ]]; then
        echo "  Installing whokaryote + tiara via pip..."
        "${env_path}/bin/pip" install --no-deps \
            'whokaryote @ git+https://github.com/LottePronk/whokaryote.git' tiara \
            > /dev/null 2>&1
        echo "  whokaryote + tiara installed"
    fi

    if [[ "${env_name}" == "dana-mag-strain" ]]; then
        echo "  Installing inStrain via pip..."
        "${env_path}/bin/pip" install --no-deps \
            'instrain @ git+https://github.com/rec3141/inStrain.git@fix-biopython-compat' \
            > /dev/null 2>&1
        echo "  inStrain installed"
    fi

    if [[ "${env_name}" == "dana-mag-pathviz" ]]; then
        # KEGGDecoder via pip (--no-deps: matplotlib already from conda)
        echo "  Installing KEGGDecoder via pip..."
        "${env_path}/bin/pip" install --no-deps KEGGDecoder > /dev/null 2>&1
        echo "  KEGGDecoder installed"

        # Clone MinPath
        local minpath_dir="${env_path}/share/minpath"
        if [[ ! -d "${minpath_dir}" ]]; then
            echo "  Cloning MinPath..."
            git clone --depth 1 https://github.com/mgtools/MinPath.git "${minpath_dir}" \
                > /dev/null 2>&1
            chmod +x "${minpath_dir}/MinPath.py" 2>/dev/null || true
            echo "  MinPath installed to ${minpath_dir}"
        fi
    fi
}

post_install_standalone() {
    local env_name="$1"
    local env_path="$2"

    # (magscot standalone env removed — Python port in bin/magscot.py uses binning env)
}

do_check() {
    echo ""
    echo "Checking conda environments in: ${ENV_DIR}"
    echo ""

    local ok=0
    local missing=0

    # Check all envs (merged + standalone)
    local all_envs=()
    for yaml in "${MERGED_YAMLS[@]}"; do
        all_envs+=("dana-mag-${yaml%.yml}")
    done
    for yaml in "${STANDALONE_YAMLS[@]}"; do
        all_envs+=("dana-mag-${yaml%.yml}")
    done

    for env_name in "${all_envs[@]}"; do
        local env_path="${ENV_DIR}/${env_name}"
        local check_bin="${ENV_CHECK[$env_name]:-}"

        printf "  %-25s" "${env_name}"

        if [[ ! -d "${env_path}" ]]; then
            echo "MISSING (not installed)"
            missing=$((missing + 1))
            continue
        fi

        if [[ -n "${check_bin}" && -x "${env_path}/bin/${check_bin}" ]]; then
            local version
            version=$("${env_path}/bin/${check_bin}" --version 2>&1 | head -1) || version="(installed)"
            echo "OK  ${version}"
            ok=$((ok + 1))
        elif [[ -n "${check_bin}" ]]; then
            echo "BROKEN (${check_bin} not found in env)"
            missing=$((missing + 1))
        else
            echo "OK  (no check binary defined)"
            ok=$((ok + 1))
        fi
    done

    # Check extra binaries
    echo ""
    echo "  Extra tool checks:"
    for env_name in "${!ENV_EXTRAS[@]}"; do
        local env_path="${ENV_DIR}/${env_name}"
        for bin in ${ENV_EXTRAS[$env_name]}; do
            printf "    %-30s" "${env_name}/${bin}"
            if [[ -x "${env_path}/bin/${bin}" ]]; then
                echo "OK"
            else
                echo "MISSING"
                missing=$((missing + 1))
            fi
        done
    done

    echo ""
    echo "Result: ${ok} envs OK, ${missing} issue(s)"
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
    local all_envs=()
    for yaml in "${MERGED_YAMLS[@]}"; do
        all_envs+=("dana-mag-${yaml%.yml}")
    done
    for yaml in "${STANDALONE_YAMLS[@]}"; do
        all_envs+=("dana-mag-${yaml%.yml}")
    done

    for env_name in "${all_envs[@]}"; do
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
        for env_name in "${all_envs[@]}"; do
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

# Install visualization dashboard dependencies (Node.js)
install_viz() {
    local viz_dir="${SCRIPT_DIR}/viz"
    if [[ ! -d "${viz_dir}" ]]; then
        echo "[INFO] No viz/ directory found, skipping dashboard setup"
        return 0
    fi
    # Use npm from pathviz conda env (nodejs>=18 is a conda dependency)
    local npm_bin="${ENV_DIR}/dana-mag-pathviz/bin/npm"
    if [[ ! -x "${npm_bin}" ]]; then
        echo "[WARNING] npm not found in dana-mag-pathviz env — skipping viz dashboard install"
        echo "  Re-run install to create the pathviz env first"
        return 0
    fi
    echo ""
    echo "Installing visualization dashboard dependencies ..."
    local node_dir
    node_dir="$(dirname "${npm_bin}")"
    (cd "${viz_dir}" && PATH="${node_dir}:${PATH}" "${npm_bin}" install --no-audit --no-fund) || {
        echo "[WARNING] viz npm install failed (non-fatal)"
        return 0
    }
    echo "[SUCCESS] Viz dashboard ready"
}

# Initialize ECOSSDB submodule (ecosystem services database)
install_ecossdb() {
    local ecossdb_dir="${SCRIPT_DIR}/ecossdb"
    if [[ -f "${ecossdb_dir}/bin/map_to_es.py" ]]; then
        echo "[OK] ECOSSDB submodule already initialized"
        return 0
    fi
    echo ""
    echo "Initializing ECOSSDB submodule (ecosystem services database) ..."
    (cd "${SCRIPT_DIR}" && git submodule update --init --depth 1 ecossdb) || {
        echo "[WARNING] ECOSSDB submodule init failed — ecosystem services profiling will be unavailable"
        echo "  To fix: cd ${SCRIPT_DIR} && git submodule update --init ecossdb"
        return 0
    }
    echo "[SUCCESS] ECOSSDB ready ($(wc -l < "${ecossdb_dir}/db/mappings/es_gene_mapping.tsv") gene-ES mappings)"
}

# Run the requested action
case "${ACTION}" in
    install) do_install; install_viz; install_ecossdb ;;
    check)   do_check ;;
    clean)   do_clean ;;
esac
