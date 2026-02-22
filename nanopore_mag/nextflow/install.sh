#!/usr/bin/env bash
set -euo pipefail

# ============================================================================
# Dana MAG Assembly Pipeline - Conda Environment Installer
# ============================================================================
#
# Creates isolated conda environments for the MAG assembly pipeline.
# All envs are prefix-installed under ./conda-envs/.
#
# Twenty-one environments are needed because of dependency conflicts:
#   dana-mag-flye    - Flye + Filtlong (Python version conflicts)
#   dana-mag-mapping - minimap2, samtools (universal, no conflicts)
#   dana-mag-semibin - SemiBin2, LorBin (ML dependencies: PyTorch isolated)
#   dana-mag-comebin - COMEBin (cloned from fork; modern Python 3.11 + numpy/torch)
#   dana-mag-binning - MetaBAT2, MaxBin2, DAS_Tool (binning suite)
#   dana-mag-prokka  - Prokka gene annotation
#   dana-mag-kofamscan - KofamScan (KEGG Orthology via adaptive HMM thresholds)
#   dana-mag-emapper  - eggNOG-mapper (COG/GO/EC/KEGG/Pfam annotation)
#   dana-mag-dbcan    - run_dbcan / dbCAN3 (CAZyme annotation)
#   dana-mag-genomad  - geNomad (virus + plasmid + provirus detection)
#   dana-mag-checkv   - CheckV (viral quality assessment)
#   dana-mag-integron  - IntegronFinder (integron + gene cassette detection)
#   dana-mag-islandpath  - IslandPath-DIMOB (genomic island detection)
#   dana-mag-macsyfinder - MacSyFinder (secretion systems + conjugation)
#   dana-mag-defensefinder - DefenseFinder (anti-phage defense systems)
#   dana-mag-kaiju    - Kaiju (protein-level taxonomy via Prokka)
#   dana-mag-kraken2  - Kraken2 (k-mer contig-level taxonomy)
#   dana-mag-checkm2  - CheckM2 (quality assessment)
#   dana-mag-tiara    - Tiara (deep learning eukaryotic contig classification)
#   dana-mag-whokaryote - Whokaryote (gene structure-based eukaryotic classification)
#   dana-mag-metaeuk  - MetaEuk (eukaryotic gene prediction, multi-exon)
#   dana-mag-rrna     - barrnap + vsearch (rRNA gene detection + SILVA classification)
#   dana-mag-pathway  - MinPath + KEGG-Decoder (parsimony pathways + biogeochemical heatmaps)
#   dana-mag-marferret - DIAMOND + Python/pandas (MarFERReT marine eukaryotic classification)
#   dana-mag-viz      - Node.js + Python/pandas/scipy (viz dashboard preprocessing + build)
#
#   dana-mag-metamdbg  - metaMDBG / nanoMDBG (de Bruijn graph ONT metagenome assembler)
#   dana-mag-myloasm   - myloasm (high-resolution ONT strain assembler)
#   dana-mag-derep     - galah, skani, sourmash (fast MAG dereplication + ANI + sketching)
#   dana-mag-drep      - dRep (established MAG dereplication; heavy deps, isolated)
#   dana-mag-instrain  - InStrain (strain-level population genomics; ancient biopython pin)
#   dana-mag-strainy   - Strainy (strain-aware assembly from long reads)
#   dana-mag-floria    - Floria (strain-aware phasing from long reads)
#   dana-mag-skder     - skDER (fast skani-based dereplication; Python 3.10 pin)
#
# BBMap (for optional dedupe) is shared with the realtime pipeline via
# symlinked YAML; its env is named dana-bbmap.
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
    flye.yml
    mapping.yml
    semibin.yml
    binning.yml
    prokka.yml
    bakta.yml
    kofamscan.yml
    emapper.yml
    dbcan.yml
    genomad.yml
    checkv.yml
    integron.yml
    islandpath.yml
    macsyfinder.yml
    defensefinder.yml
    kaiju.yml
    kraken2.yml
    checkm2.yml
    tiara.yml
    whokaryote.yml
    metaeuk.yml
    rrna.yml
    pathway.yml
    marferret.yml
    bbmap.yml
    viz.yml
    metamdbg.yml
    myloasm.yml
    derep.yml
    drep.yml
    instrain.yml
    strainy.yml
    floria.yml
    skder.yml
    comebin.yml
)

# Which tool binary to check for each environment
declare -A ENV_CHECK=(
    [dana-mag-flye]="flye"
    [dana-mag-mapping]="minimap2"
    [dana-mag-semibin]="SemiBin2"
    [dana-mag-binning]="metabat2"
    [dana-mag-comebin]="run_comebin.sh"
    [dana-mag-prokka]="prokka"
    [dana-mag-bakta]="bakta"
    [dana-mag-kofamscan]="exec_annotation"
    [dana-mag-emapper]="emapper.py"
    [dana-mag-dbcan]="run_dbcan"
    [dana-mag-genomad]="genomad"
    [dana-mag-checkv]="checkv"
    [dana-mag-integron]="integron_finder"
    [dana-mag-islandpath]="hmmscan"
    [dana-mag-macsyfinder]="macsyfinder"
    [dana-mag-defensefinder]="defense-finder"
    [dana-mag-kaiju]="kaiju"
    [dana-mag-kraken2]="kraken2"
    [dana-mag-checkm2]="checkm2"
    [dana-mag-tiara]="tiara"
    [dana-mag-whokaryote]="whokaryote.py"
    [dana-mag-metaeuk]="metaeuk"
    [dana-mag-rrna]="barrnap"
    [dana-mag-pathway]="KEGG-decoder"
    [dana-mag-marferret]="diamond"
    [dana-bbmap]="bbduk.sh"
    [dana-mag-viz]="node"
    [dana-mag-metamdbg]="metaMDBG"
    [dana-mag-myloasm]="myloasm"
    [dana-mag-derep]="galah"
    [dana-mag-drep]="dRep"
    [dana-mag-instrain]="inStrain"
    [dana-mag-strainy]="strainy"
    [dana-mag-floria]="floria"
    [dana-mag-skder]="skder"
)

# Additional binaries to verify
declare -A ENV_EXTRAS=(
    [dana-mag-flye]="filtlong nextflow java"
    [dana-mag-mapping]="samtools"
    [dana-mag-semibin]="LorBin"
    [dana-mag-comebin]="gen_cov_file.sh"
    [dana-mag-binning]="run_MaxBin.pl DAS_Tool jgi_summarize_bam_contig_depths"
    [dana-mag-viz]="npm python3"
    [dana-mag-derep]="skani sourmash"
)

yaml_to_envname() {
    local yaml="$1"
    local base="${yaml%.yml}"
    # bbmap.yml -> dana-bbmap (matches realtime pipeline)
    if [[ "$base" == "bbmap" ]]; then
        echo "dana-bbmap"
    else
        echo "dana-mag-${base}"
    fi
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

        # Post-install: ensure conda is on PATH for Nextflow env activation.
        # Nextflow's .command.run calls `conda info --json` to activate prefix envs.
        # The flye env hosts Nextflow and its bin/ is on PATH for all spawned processes.
        if [[ "${env_name}" == "dana-mag-flye" ]]; then
            local conda_bin
            conda_bin=$(which conda 2>/dev/null || echo "")
            if [[ -n "${conda_bin}" && ! -e "${env_path}/bin/conda" ]]; then
                ln -sf "${conda_bin}" "${env_path}/bin/conda"
                echo "  Symlinked conda into flye env for Nextflow activation"
            fi
        fi

        # Post-install: clone COMEBin fork and wire up bin/ scripts
        if [[ "${env_name}" == "dana-mag-comebin" ]]; then
            echo "  Installing COMEBin from fork..."
            local comebin_repo="https://github.com/rec3141/COMEBin.git"
            local comebin_branch="codex/recent-dependencies-fix"
            local comebin_dir="${env_path}/share/COMEBin"
            git clone --depth 1 -b "${comebin_branch}" "${comebin_repo}" "${comebin_dir}" \
                > /dev/null 2>&1
            # Create wrapper that cd's into the COMEBin source dir and runs the
            # inner run_comebin.sh (which uses pwd-relative paths to ../auxiliary/).
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

        # Post-install: install KEGGDecoder + clone MinPath into the pathway env
        if [[ "${env_name}" == "dana-mag-pathway" ]]; then
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
            else
                echo "  MinPath already present at ${minpath_dir}"
            fi
        fi

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

do_check() {
    echo ""
    echo "Checking conda environments in: ${ENV_DIR}"
    echo ""

    local ok=0
    local missing=0

    for yaml in "${ENV_YAMLS[@]}"; do
        local env_name
        env_name=$(yaml_to_envname "$yaml")
        local env_path="${ENV_DIR}/${env_name}"
        local check_bin="${ENV_CHECK[$env_name]}"

        printf "  %-25s" "${env_name}"

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

# Install visualization dashboard dependencies (Node.js)
install_viz() {
    local viz_dir="${SCRIPT_DIR}/viz"
    if [[ ! -d "${viz_dir}" ]]; then
        echo "[INFO] No viz/ directory found, skipping dashboard setup"
        return 0
    fi
    if ! command -v npm &>/dev/null; then
        echo "[WARNING] npm not found — skipping viz dashboard install"
        echo "  Install Node.js >= 18 to enable the MAG dashboard"
        return 0
    fi
    echo ""
    echo "Installing visualization dashboard dependencies ..."
    (cd "${viz_dir}" && npm install --no-audit --no-fund) || {
        echo "[WARNING] viz npm install failed (non-fatal)"
        return 0
    }
    echo "[SUCCESS] Viz dashboard ready. Run: cd viz && npm run dev"
}

# Run the requested action
case "${ACTION}" in
    install) do_install; install_viz ;;
    check)   do_check ;;
    clean)   do_clean ;;
esac
