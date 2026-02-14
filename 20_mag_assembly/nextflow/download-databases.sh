#!/usr/bin/env bash
set -euo pipefail

# ============================================================================
# Dana MAG Assembly Pipeline - Database Downloader
# ============================================================================
#
# Downloads reference databases required by optional pipeline modules.
# Databases are large (1-10 GB each) and are NOT included in the conda
# environments. Run this script once before using the tools that need them.
#
# Usage:
#   ./download-databases.sh                    # Interactive: choose databases
#   ./download-databases.sh --all              # Download all databases
#   ./download-databases.sh --genomad          # Download geNomad database only
#   ./download-databases.sh --checkv           # Download CheckV database only
#   ./download-databases.sh --checkm2          # Download CheckM2 database only
#   ./download-databases.sh --kaiju            # Download Kaiju RefSeq database only
#   ./download-databases.sh --dir /custom/path # Custom database directory
#   ./download-databases.sh --list             # Show available databases
#
# Default directory: ./databases/
# ============================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DB_DIR="${SCRIPT_DIR}/databases"
ENV_DIR="${SCRIPT_DIR}/conda-envs"

# Parse arguments
DOWNLOAD_GENOMAD=false
DOWNLOAD_CHECKV=false
DOWNLOAD_CHECKM2=false
DOWNLOAD_KAIJU=false
DOWNLOAD_ALL=false
LIST_ONLY=false
INTERACTIVE=true

while (( $# )); do
    case "$1" in
        --dir)       DB_DIR="$2"; shift 2 ;;
        --all)       DOWNLOAD_ALL=true; INTERACTIVE=false; shift ;;
        --genomad)   DOWNLOAD_GENOMAD=true; INTERACTIVE=false; shift ;;
        --checkv)    DOWNLOAD_CHECKV=true; INTERACTIVE=false; shift ;;
        --checkm2)   DOWNLOAD_CHECKM2=true; INTERACTIVE=false; shift ;;
        --kaiju)     DOWNLOAD_KAIJU=true; INTERACTIVE=false; shift ;;
        --list)      LIST_ONLY=true; INTERACTIVE=false; shift ;;
        -h|--help)
            sed -n '/^# Usage:/,/^# ====/p' "$0" | head -n -1 | sed 's/^# //'
            exit 0 ;;
        *) echo "[ERROR] Unknown argument: $1" >&2; exit 1 ;;
    esac
done

if $DOWNLOAD_ALL; then
    DOWNLOAD_GENOMAD=true
    DOWNLOAD_CHECKV=true
    DOWNLOAD_CHECKM2=true
    DOWNLOAD_KAIJU=true
fi

# ============================================================================
# Database info
# ============================================================================

show_databases() {
    echo ""
    echo "Available databases:"
    echo ""
    printf "  %-12s %-8s  %s\n" "Database" "Size" "Description"
    printf "  %-12s %-8s  %s\n" "--------" "----" "-----------"
    printf "  %-12s %-8s  %s\n" "genomad"  "~3.5 GB" "geNomad marker profiles + MMseqs2 (virus + plasmid detection)"
    printf "  %-12s %-8s  %s\n" "checkv"   "~1.4 GB" "CheckV reference genomes (viral quality assessment)"
    printf "  %-12s %-8s  %s\n" "checkm2"  "~3.5 GB" "CheckM2 DIAMOND db (MAG quality assessment)"
    printf "  %-12s %-8s  %s\n" "kaiju"    "~47 GB"  "Kaiju RefSeq protein db (contig-level taxonomy)"
    echo ""
    echo "Default download directory: ${DB_DIR}"
    echo ""
    echo "After downloading, pass database paths to the pipeline:"
    echo "  --genomad_db ${DB_DIR}/genomad_db"
    echo "  --checkv_db  ${DB_DIR}/checkv_db"
    echo "  --checkm2_db ${DB_DIR}/checkm2_db"
    echo "  --kaiju_db   ${DB_DIR}/kaiju_db"
    echo ""
}

if $LIST_ONLY; then
    show_databases
    exit 0
fi

# ============================================================================
# Interactive selection
# ============================================================================

if $INTERACTIVE; then
    show_databases
    echo "Select databases to download (space-separated, or 'all'):"
    echo "  1) genomad   - geNomad (virus + plasmid detection)"
    echo "  2) checkv    - CheckV (viral quality assessment)"
    echo "  3) checkm2   - CheckM2 (MAG quality assessment)"
    echo "  4) kaiju     - Kaiju RefSeq proteins (contig taxonomy, ~47 GB)"
    echo "  5) all       - All databases"
    echo ""
    read -rp "Choice [1-5, or names]: " choice

    case "$choice" in
        1|genomad)  DOWNLOAD_GENOMAD=true ;;
        2|checkv)   DOWNLOAD_CHECKV=true ;;
        3|checkm2)  DOWNLOAD_CHECKM2=true ;;
        4|kaiju)    DOWNLOAD_KAIJU=true ;;
        5|all)      DOWNLOAD_GENOMAD=true; DOWNLOAD_CHECKV=true; DOWNLOAD_CHECKM2=true; DOWNLOAD_KAIJU=true ;;
        *)
            # Parse space-separated names
            for item in $choice; do
                case "$item" in
                    1|genomad)  DOWNLOAD_GENOMAD=true ;;
                    2|checkv)   DOWNLOAD_CHECKV=true ;;
                    3|checkm2)  DOWNLOAD_CHECKM2=true ;;
                    4|kaiju)    DOWNLOAD_KAIJU=true ;;
                    all)        DOWNLOAD_GENOMAD=true; DOWNLOAD_CHECKV=true; DOWNLOAD_CHECKM2=true; DOWNLOAD_KAIJU=true ;;
                    *) echo "[WARNING] Unknown selection: $item" >&2 ;;
                esac
            done
            ;;
    esac

    if ! $DOWNLOAD_GENOMAD && ! $DOWNLOAD_CHECKV && ! $DOWNLOAD_CHECKM2 && ! $DOWNLOAD_KAIJU; then
        echo "No databases selected. Exiting."
        exit 0
    fi
fi

# ============================================================================
# Download functions
# ============================================================================

mkdir -p "${DB_DIR}"

download_genomad() {
    local db_path="${DB_DIR}/genomad_db"
    if [ -d "${db_path}" ] && [ -f "${db_path}/genomad_db" ]; then
        echo "[INFO] geNomad database already exists at ${db_path}"
        echo "  Delete ${db_path} and re-run to force re-download."
        return 0
    fi

    echo ""
    echo "[INFO] Downloading geNomad database (~3.5 GB)..."
    echo "  Destination: ${db_path}"

    local genomad_bin="${ENV_DIR}/dana-mag-genomad/bin/genomad"
    if [ ! -x "${genomad_bin}" ]; then
        echo "[ERROR] geNomad not installed. Run ./install.sh first." >&2
        return 1
    fi

    "${genomad_bin}" download-database "${DB_DIR}"
    echo "[SUCCESS] geNomad database downloaded to ${db_path}"
    echo "  Use with: --genomad_db ${db_path}"
}

download_checkv() {
    local db_path="${DB_DIR}/checkv_db"
    if [ -d "${db_path}" ] && [ -f "${db_path}/README.txt" ]; then
        echo "[INFO] CheckV database already exists at ${db_path}"
        echo "  Delete ${db_path} and re-run to force re-download."
        return 0
    fi

    echo ""
    echo "[INFO] Downloading CheckV database (~1.4 GB)..."
    echo "  Destination: ${db_path}"

    local checkv_bin="${ENV_DIR}/dana-mag-checkv/bin/checkv"
    if [ ! -x "${checkv_bin}" ]; then
        echo "[ERROR] CheckV not installed. Run ./install.sh first." >&2
        return 1
    fi

    "${checkv_bin}" download_database "${db_path}"
    # CheckV may download to a versioned subdirectory (checkv-db-v*)
    # If so, move contents up to the expected path
    local downloaded_dir
    downloaded_dir=$(ls -d "${db_path}"/checkv-db-v* 2>/dev/null | head -1)
    if [ -n "${downloaded_dir}" ] && [ -d "${downloaded_dir}" ]; then
        mv "${downloaded_dir}"/* "${db_path}"/
        rmdir "${downloaded_dir}"
    fi
    echo "[SUCCESS] CheckV database downloaded to ${db_path}"
    echo "  Use with: --checkv_db ${db_path}"
}

download_checkm2() {
    local db_path="${DB_DIR}/checkm2_db"
    if [ -d "${db_path}" ] && ls "${db_path}"/*.dmnd 1>/dev/null 2>&1; then
        echo "[INFO] CheckM2 database already exists at ${db_path}"
        echo "  Delete ${db_path} and re-run to force re-download."
        return 0
    fi

    echo ""
    echo "[INFO] Downloading CheckM2 database (~3.5 GB)..."
    echo "  Destination: ${db_path}"

    local checkm2_bin="${ENV_DIR}/dana-mag-checkm2/bin/checkm2"
    if [ ! -x "${checkm2_bin}" ]; then
        echo "[ERROR] CheckM2 not installed. Run ./install.sh first." >&2
        return 1
    fi

    mkdir -p "${db_path}"
    "${checkm2_bin}" database --download --path "${db_path}"
    echo "[SUCCESS] CheckM2 database downloaded to ${db_path}"
    echo "  Use with: --checkm2_db ${db_path}"
}

download_kaiju() {
    local db_path="${DB_DIR}/kaiju_db"
    if [ -d "${db_path}" ] && ls "${db_path}"/*.fmi 1>/dev/null 2>&1; then
        echo "[INFO] Kaiju database already exists at ${db_path}"
        echo "  Delete ${db_path} and re-run to force re-download."
        return 0
    fi

    echo ""
    echo "[INFO] Downloading Kaiju RefSeq reference database (~47 GB)..."
    echo "  This is a large download and will take a while."
    echo "  Destination: ${db_path}"

    local kaiju_makedb="${ENV_DIR}/dana-mag-kaiju/bin/kaiju-makedb"
    if [ ! -x "${kaiju_makedb}" ]; then
        echo "[ERROR] Kaiju not installed. Run ./install.sh first." >&2
        return 1
    fi

    mkdir -p "${db_path}"
    # kaiju-makedb downloads sequences and builds the FM-index
    # -s refseq_ref: RefSeq reference genomes (bacteria, archaea, viruses)
    (cd "${db_path}" && "${kaiju_makedb}" -s refseq_ref)
    echo "[SUCCESS] Kaiju database downloaded to ${db_path}"
    echo "  Use with: --kaiju_db ${db_path}"
}

# ============================================================================
# Execute downloads
# ============================================================================

echo ""
echo "Database directory: ${DB_DIR}"

failed=0

if $DOWNLOAD_GENOMAD; then
    download_genomad || failed=$((failed + 1))
fi

if $DOWNLOAD_CHECKV; then
    download_checkv || failed=$((failed + 1))
fi

if $DOWNLOAD_CHECKM2; then
    download_checkm2 || failed=$((failed + 1))
fi

if $DOWNLOAD_KAIJU; then
    download_kaiju || failed=$((failed + 1))
fi

echo ""
if (( failed > 0 )); then
    echo "[WARNING] ${failed} database download(s) failed"
    exit 1
else
    echo "[SUCCESS] All selected databases downloaded to: ${DB_DIR}"
    echo ""
    echo "Run the pipeline with database paths:"
    $DOWNLOAD_GENOMAD && echo "  --genomad_db ${DB_DIR}/genomad_db"
    $DOWNLOAD_CHECKV  && echo "  --checkv_db  ${DB_DIR}/checkv_db"
    $DOWNLOAD_CHECKM2 && echo "  --checkm2_db ${DB_DIR}/checkm2_db"
    $DOWNLOAD_KAIJU   && echo "  --kaiju_db   ${DB_DIR}/kaiju_db"
fi
