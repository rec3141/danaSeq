#!/usr/bin/env bash
set -euo pipefail

# ============================================================================
# danaSeq - Database Downloader
# ============================================================================
#
# Downloads reference databases required by the nanopore_mag and illumina_mag
# pipelines. Databases are large (1-10 GB each) and are NOT included in the
# conda environments. Run this script once before using the tools that need them.
#
# Usage:
#   ./download-databases.sh                    # Interactive: choose databases
#   ./download-databases.sh --all              # Download all databases
#   ./download-databases.sh --human            # Download BBTools human reference (~4 GB, for illumina_mag)
#   ./download-databases.sh --genomad          # Download geNomad database only
#   ./download-databases.sh --checkv           # Download CheckV database only
#   ./download-databases.sh --checkm            # Download CheckM v1 data (needed by dRep/COMEBin)
#   ./download-databases.sh --checkm2          # Download CheckM2 database only
#   ./download-databases.sh --kaiju            # Download Kaiju RefSeq database only
#   ./download-databases.sh --macsyfinder      # Download MacSyFinder models (TXSScan + CONJScan)
#   ./download-databases.sh --defensefinder    # Download DefenseFinder models (~100 MB)
#   ./download-databases.sh --bakta            # Download Bakta annotation database (~37 GB)
#   ./download-databases.sh --bakta-light      # Download Bakta light database (~1.4 GB)
#   ./download-databases.sh --kofam            # Download KOfam profiles (~4 GB)
#   ./download-databases.sh --eggnog           # Download eggNOG diamond database (~12 GB)
#   ./download-databases.sh --dbcan            # Download dbCAN databases (~2 GB)
#   ./download-databases.sh --metaeuk          # Download MetaEuk OrthoDB Eukaryota (~23 GB)
#   ./download-databases.sh --kraken2         # Download Kraken2 PlusPFP-8 (~8 GB)
#   ./download-databases.sh --silva           # Download SILVA SSU + LSU NR99 (~900 MB)
#   ./download-databases.sh --marferret       # Download MarFERReT marine eukaryotic database (~9 GB)
#   ./download-databases.sh --gtdbtk          # Download GTDB-Tk r226 reference data (~132 GB)
#   ./download-databases.sh --antismash      # Download antiSMASH databases (~2 GB)
#   ./download-databases.sh --docker            # Use Docker to run tool CLIs
#   ./download-databases.sh --apptainer         # Use Apptainer/Singularity (auto-pulls SIF)
#   ./download-databases.sh --container         # Auto-detect: apptainer > singularity > docker
#   ./download-databases.sh --docker --image IMG # Use a custom Docker/Apptainer image
#   ./download-databases.sh --db_dir /custom/path # Custom database directory
#   ./download-databases.sh --list             # Show available databases
#
# Default directory: ./databases/
# ============================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DB_DIR="${SCRIPT_DIR}/databases"
# Conda envs may live under either pipeline's directory
NANOPORE_ENV_DIR="${SCRIPT_DIR}/nanopore_mag/conda-envs"
ILLUMINA_ENV_DIR="${SCRIPT_DIR}/illumina_mag/conda-envs"
# Legacy default — nanopore_mag envs (used by most download functions)
ENV_DIR="${NANOPORE_ENV_DIR}"
USE_CONTAINER=false
CONTAINER_RUNTIME=""   # docker, apptainer, or singularity (set by --docker/--apptainer/--container)
CONTAINER_IMAGE="ghcr.io/rec3141/danaseq-mag:latest"
SIF_PATH=""            # resolved path to .sif file (apptainer/singularity only)

# Container helper — runs a command inside Docker or Apptainer with DB_DIR mounted at /data/db
container_run() {
    case "$CONTAINER_RUNTIME" in
        docker)
            docker run --rm --user "$(id -u):$(id -g)" --entrypoint "" \
                -v "${DB_DIR}:/data/db" "${CONTAINER_IMAGE}" "$@"
            ;;
        apptainer|singularity)
            # Ensure SIF image exists
            if [[ -z "$SIF_PATH" ]]; then
                SIF_PATH="${SCRIPT_DIR}/.danaseq-mag.sif"
                if [[ ! -f "$SIF_PATH" ]]; then
                    echo "[INFO] Pulling container image (one-time download)..."
                    echo "  Source: docker://${CONTAINER_IMAGE}"
                    echo "  Destination: ${SIF_PATH}"
                    "$CONTAINER_RUNTIME" pull "$SIF_PATH" "docker://${CONTAINER_IMAGE}"
                fi
            fi
            # Override host SSL env vars (RHEL sets REQUESTS_CA_BUNDLE to paths
            # that don't exist inside the Debian-based container) — point to the
            # container's own CA bundle so Python requests/urllib and wget all work
            local container_ca="/etc/ssl/certs/ca-certificates.crt"
            "$CONTAINER_RUNTIME" exec \
                --writable-tmpfs \
                --env "REQUESTS_CA_BUNDLE=${container_ca}" \
                --env "SSL_CERT_FILE=${container_ca}" \
                --env "CURL_CA_BUNDLE=${container_ca}" \
                --bind "${DB_DIR}:/data/db" "$SIF_PATH" "$@"
            ;;
    esac
}

# Detect conda/mamba — prefer pre-built envs, fall back to PATH (warning deferred to after arg parsing)
if [ -x "${ENV_DIR}/dana-mag-assembly/bin/mamba" ]; then
    CONDA_CMD="${ENV_DIR}/dana-mag-assembly/bin/mamba"
elif [ -x "${ENV_DIR}/dana-mag-assembly/bin/conda" ]; then
    CONDA_CMD="${ENV_DIR}/dana-mag-assembly/bin/conda"
elif command -v mamba &>/dev/null; then
    CONDA_CMD="mamba"
elif command -v conda &>/dev/null; then
    CONDA_CMD="conda"
else
    CONDA_CMD=""
fi

# Parse arguments
DOWNLOAD_HUMAN=false
DOWNLOAD_GENOMAD=false
DOWNLOAD_CHECKV=false
DOWNLOAD_CHECKM=false
DOWNLOAD_CHECKM2=false
DOWNLOAD_KAIJU=false
DOWNLOAD_MACSYFINDER=false
DOWNLOAD_DEFENSEFINDER=false
DOWNLOAD_BAKTA=false
DOWNLOAD_BAKTA_LIGHT=false
DOWNLOAD_KOFAM=false
DOWNLOAD_EGGNOG=false
DOWNLOAD_DBCAN=false
DOWNLOAD_METAEUK=false
DOWNLOAD_KRAKEN2=false
DOWNLOAD_SILVA=false
DOWNLOAD_MARFERRET=false
DOWNLOAD_GTDBTK=false
DOWNLOAD_ANTISMASH=false
DOWNLOAD_MAGSCOT=false
DOWNLOAD_ALL=false
LIST_ONLY=false
INTERACTIVE=true

while (( $# )); do
    case "$1" in
        --docker)    USE_CONTAINER=true; CONTAINER_RUNTIME=docker; shift ;;
        --apptainer|--singularity)
                     USE_CONTAINER=true; CONTAINER_RUNTIME="${1#--}"; shift ;;
        --container) USE_CONTAINER=true; CONTAINER_RUNTIME=auto; shift ;;
        --image)     CONTAINER_IMAGE="$2"; shift 2 ;;
        --sif)       SIF_PATH="$2"; shift 2 ;;
        --db_dir|--dir) DB_DIR="$2"; shift 2 ;;
        --all)       DOWNLOAD_ALL=true; INTERACTIVE=false; shift ;;
        --human)     DOWNLOAD_HUMAN=true; INTERACTIVE=false; shift ;;
        --genomad)   DOWNLOAD_GENOMAD=true; INTERACTIVE=false; shift ;;
        --checkv)    DOWNLOAD_CHECKV=true; INTERACTIVE=false; shift ;;
        --checkm)    DOWNLOAD_CHECKM=true; INTERACTIVE=false; shift ;;
        --checkm2)   DOWNLOAD_CHECKM2=true; INTERACTIVE=false; shift ;;
        --kaiju)     DOWNLOAD_KAIJU=true; INTERACTIVE=false; shift ;;
        --macsyfinder) DOWNLOAD_MACSYFINDER=true; INTERACTIVE=false; shift ;;
        --defensefinder) DOWNLOAD_DEFENSEFINDER=true; INTERACTIVE=false; shift ;;
        --bakta)       DOWNLOAD_BAKTA=true; INTERACTIVE=false; shift ;;
        --bakta-light) DOWNLOAD_BAKTA_LIGHT=true; INTERACTIVE=false; shift ;;
        --kofam)     DOWNLOAD_KOFAM=true; INTERACTIVE=false; shift ;;
        --eggnog)    DOWNLOAD_EGGNOG=true; INTERACTIVE=false; shift ;;
        --dbcan)     DOWNLOAD_DBCAN=true; INTERACTIVE=false; shift ;;
        --metaeuk)   DOWNLOAD_METAEUK=true; INTERACTIVE=false; shift ;;
        --kraken2)   DOWNLOAD_KRAKEN2=true; INTERACTIVE=false; shift ;;
        --silva)     DOWNLOAD_SILVA=true; INTERACTIVE=false; shift ;;
        --marferret) DOWNLOAD_MARFERRET=true; INTERACTIVE=false; shift ;;
        --gtdbtk)    DOWNLOAD_GTDBTK=true; INTERACTIVE=false; shift ;;
        --antismash) DOWNLOAD_ANTISMASH=true; INTERACTIVE=false; shift ;;
        --magscot)   DOWNLOAD_MAGSCOT=true; INTERACTIVE=false; shift ;;
        --list)      LIST_ONLY=true; INTERACTIVE=false; shift ;;
        -h|--help)
            sed -n '/^# Usage:/,/^# ====/p' "$0" | head -n -1 | sed 's/^# //'
            exit 0 ;;
        *) echo "[ERROR] Unknown argument: $1" >&2; exit 1 ;;
    esac
done

if $DOWNLOAD_ALL; then
    DOWNLOAD_HUMAN=true
    DOWNLOAD_GENOMAD=true
    DOWNLOAD_CHECKV=true
    DOWNLOAD_CHECKM=true
    DOWNLOAD_CHECKM2=true
    DOWNLOAD_KAIJU=true
    DOWNLOAD_MACSYFINDER=true
    DOWNLOAD_DEFENSEFINDER=true
    DOWNLOAD_BAKTA=true
    DOWNLOAD_BAKTA_LIGHT=true
    DOWNLOAD_KOFAM=true
    DOWNLOAD_EGGNOG=true
    DOWNLOAD_DBCAN=true
    DOWNLOAD_METAEUK=true
    DOWNLOAD_KRAKEN2=true
    DOWNLOAD_SILVA=true
    DOWNLOAD_MARFERRET=true
    DOWNLOAD_GTDBTK=true
    DOWNLOAD_ANTISMASH=true
    DOWNLOAD_MAGSCOT=true
fi

# Auto-detect container runtime if --container was used
if [[ "$CONTAINER_RUNTIME" == "auto" ]]; then
    if command -v apptainer &>/dev/null; then
        CONTAINER_RUNTIME=apptainer
    elif command -v singularity &>/dev/null; then
        CONTAINER_RUNTIME=singularity
    elif command -v docker &>/dev/null; then
        CONTAINER_RUNTIME=docker
    else
        echo "[ERROR] --container requires docker, apptainer, or singularity on PATH" >&2
        exit 1
    fi
    echo "[INFO] Auto-detected container runtime: ${CONTAINER_RUNTIME}"
fi

# Validate chosen runtime is available
if $USE_CONTAINER && [[ "$CONTAINER_RUNTIME" != "auto" ]]; then
    if ! command -v "$CONTAINER_RUNTIME" &>/dev/null; then
        echo "[ERROR] ${CONTAINER_RUNTIME} not found on PATH" >&2
        exit 1
    fi
fi

# Warn if no conda and not using a container
if [[ -z "$CONDA_CMD" ]] && ! $USE_CONTAINER; then
    echo "[WARNING] conda/mamba not found — tool-based downloads will require --docker or --apptainer" >&2
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
    printf "  %-12s %-8s  %s\n" "human"   "~4 GB"  "BBTools masked human reference (illumina_mag decontamination)"
    printf "  %-12s %-8s  %s\n" "genomad"  "~3.5 GB" "geNomad marker profiles + MMseqs2 (virus + plasmid detection)"
    printf "  %-12s %-8s  %s\n" "checkv"   "~1.4 GB" "CheckV reference genomes (viral quality assessment)"
    printf "  %-12s %-8s  %s\n" "checkm"   "~280 MB" "CheckM v1 reference data (needed by dRep + COMEBin)"
    printf "  %-12s %-8s  %s\n" "checkm2"  "~3.5 GB" "CheckM2 DIAMOND db (MAG quality assessment)"
    printf "  %-12s %-8s  %s\n" "kaiju"    "~47 GB"  "Kaiju RefSeq protein db (contig-level taxonomy)"
    printf "  %-12s %-8s  %s\n" "macsyfinder" "~50 MB" "MacSyFinder models: TXSScan + CONJScan (secretion + conjugation)"
    printf "  %-12s %-8s  %s\n" "defensefinder" "~100 MB" "DefenseFinder models: ~280 defense system HMM profiles"
    printf "  %-12s %-8s  %s\n" "bakta"       "~37 GB"  "Bakta full annotation db (UniProt, AMRFinderPlus, Pfam, etc.)"
    printf "  %-12s %-8s  %s\n" "bakta-light" "~1.4 GB" "Bakta light annotation db (faster downloads, reduced annotation)"
    printf "  %-12s %-8s  %s\n" "kofam"    "~4 GB"   "KOfam profiles + ko_list (KEGG Orthology via HMM)"
    printf "  %-12s %-8s  %s\n" "eggnog"   "~12 GB"  "eggNOG-mapper DIAMOND db (COG/GO/EC/KEGG/Pfam)"
    printf "  %-12s %-8s  %s\n" "dbcan"    "~2 GB"   "dbCAN HMM + DIAMOND + substrate db (CAZyme annotation)"
    printf "  %-12s %-8s  %s\n" "metaeuk"  "~8.5 GB" "MetaEuk OrthoDB v11 Eukaryota (eukaryotic gene prediction)"
    printf "  %-12s %-8s  %s\n" "kraken2"  "~8 GB"   "Kraken2 PlusPFP-8 (k-mer contig-level taxonomy)"
    printf "  %-12s %-8s  %s\n" "silva"    "~900 MB" "SILVA 138.2 SSU + LSU NR99 (rRNA gene classification)"
    printf "  %-12s %-8s  %s\n" "marferret" "~9 GB"  "MarFERReT v1.1.1 marine eukaryotic proteins (DIAMOND + taxonomy + Pfam)"
    printf "  %-12s %-8s  %s\n" "gtdbtk"   "~132 GB" "GTDB-Tk r226 reference data (phylogenetic MAG classification)"
    echo ""
    echo "  Note: Tiara, Whokaryote, and IslandPath HMM profiles are bundled with"
    echo "  their conda packages (no separate database download needed)."
    echo "  ECOSSDB gene-ES mappings and ontologies are bundled in the ecossdb/"
    echo "  submodule (initialized by install.sh, no download needed)."
    echo ""
    echo "Default download directory: ${DB_DIR}"
    echo ""
    echo "After downloading, pass database paths to the pipeline:"
    echo "  --human_ref  ${DB_DIR}/human_ref   (illumina_mag pipeline)"
    echo "  --genomad_db ${DB_DIR}/genomad_db"
    echo "  --checkv_db  ${DB_DIR}/checkv_db"
    echo "  --checkm_data ${DB_DIR}/checkm_data"
    echo "  --checkm2_db ${DB_DIR}/checkm2_db"
    echo "  --kaiju_db   ${DB_DIR}/kaiju_db"
    echo "  --macsyfinder_models ${DB_DIR}/macsyfinder_models"
    echo "  --defensefinder_models ${DB_DIR}/defensefinder_models"
    echo "  --bakta_db         ${DB_DIR}/bakta/db"
    echo "  --bakta_light_db   ${DB_DIR}/bakta/db-light"
    echo "  --kofam_db   ${DB_DIR}/kofam_db"
    echo "  --eggnog_db  ${DB_DIR}/eggnog_db"
    echo "  --dbcan_db   ${DB_DIR}/dbcan_db"
    echo "  --metaeuk_db ${DB_DIR}/metaeuk_db/metaeuk_db"
    echo "  --kraken2_db ${DB_DIR}/kraken2_db"
    echo "  --silva_ssu_db ${DB_DIR}/silva_db/SILVA_138.2_SSURef_NR99.fasta"
    echo "  --silva_lsu_db ${DB_DIR}/silva_db/SILVA_138.2_LSURef_NR99.fasta"
    echo "  --marferret_db ${DB_DIR}/marferret_db"
    echo "  --gtdbtk_db  ${DB_DIR}/gtdbtk_db"
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
    echo "  1) human     - BBTools masked human reference (illumina_mag, ~4 GB)"
    echo "  2) genomad   - geNomad (virus + plasmid detection)"
    echo "  3) checkv    - CheckV (viral quality assessment)"
    echo "  4) checkm    - CheckM v1 reference data (needed by dRep + COMEBin, ~280 MB)"
    echo "  5) checkm2   - CheckM2 (MAG quality assessment)"
    echo "  6) kaiju     - Kaiju RefSeq proteins (contig taxonomy, ~47 GB)"
    echo "  7) macsyfinder - MacSyFinder models: TXSScan + CONJScan (~50 MB)"
    echo "  8) defensefinder - DefenseFinder models: anti-phage defense (~100 MB)"
    echo "  9) bakta       - Bakta full annotation database (~37 GB)"
    echo " 10) bakta-light - Bakta light annotation database (~1.4 GB)"
    echo " 11) kofam     - KOfam profiles (KEGG Orthology, ~4 GB)"
    echo " 12) eggnog    - eggNOG-mapper DIAMOND db (~12 GB)"
    echo " 13) dbcan     - dbCAN HMM + DIAMOND db (CAZyme, ~2 GB)"
    echo " 14) metaeuk   - MetaEuk OrthoDB v12 Eukaryota (~23 GB)"
    echo " 15) kraken2   - Kraken2 PlusPFP-8 (contig taxonomy, ~8 GB)"
    echo " 16) silva     - SILVA SSU + LSU NR99 (rRNA classification, ~900 MB)"
    echo " 17) marferret - MarFERReT marine eukaryotic proteins (~9 GB)"
    echo " 18) all       - All databases"
    echo ""
    read -rp "Choice [1-18, or names]: " choice

    case "$choice" in
        1|human)        DOWNLOAD_HUMAN=true ;;
        2|genomad)      DOWNLOAD_GENOMAD=true ;;
        3|checkv)       DOWNLOAD_CHECKV=true ;;
        4|checkm)       DOWNLOAD_CHECKM=true ;;
        5|checkm2)      DOWNLOAD_CHECKM2=true ;;
        6|kaiju)        DOWNLOAD_KAIJU=true ;;
        7|macsyfinder)  DOWNLOAD_MACSYFINDER=true ;;
        8|defensefinder) DOWNLOAD_DEFENSEFINDER=true ;;
        9|bakta)        DOWNLOAD_BAKTA=true ;;
        10|bakta-light)  DOWNLOAD_BAKTA_LIGHT=true ;;
        11|kofam)       DOWNLOAD_KOFAM=true ;;
        12|eggnog)      DOWNLOAD_EGGNOG=true ;;
        13|dbcan)       DOWNLOAD_DBCAN=true ;;
        14|metaeuk)     DOWNLOAD_METAEUK=true ;;
        15|kraken2)     DOWNLOAD_KRAKEN2=true ;;
        16|silva)       DOWNLOAD_SILVA=true ;;
        17|marferret)   DOWNLOAD_MARFERRET=true ;;
        18|all)         DOWNLOAD_HUMAN=true; DOWNLOAD_GENOMAD=true; DOWNLOAD_CHECKV=true; DOWNLOAD_CHECKM=true; DOWNLOAD_CHECKM2=true; DOWNLOAD_KAIJU=true; DOWNLOAD_MACSYFINDER=true; DOWNLOAD_DEFENSEFINDER=true; DOWNLOAD_BAKTA=true; DOWNLOAD_BAKTA_LIGHT=true; DOWNLOAD_KOFAM=true; DOWNLOAD_EGGNOG=true; DOWNLOAD_DBCAN=true; DOWNLOAD_METAEUK=true; DOWNLOAD_KRAKEN2=true; DOWNLOAD_SILVA=true; DOWNLOAD_MARFERRET=true ;;
        *)
            # Parse space-separated names
            for item in $choice; do
                case "$item" in
                    1|human)        DOWNLOAD_HUMAN=true ;;
                    2|genomad)      DOWNLOAD_GENOMAD=true ;;
                    3|checkv)       DOWNLOAD_CHECKV=true ;;
                    4|checkm)       DOWNLOAD_CHECKM=true ;;
                    5|checkm2)      DOWNLOAD_CHECKM2=true ;;
                    6|kaiju)        DOWNLOAD_KAIJU=true ;;
                    7|macsyfinder)  DOWNLOAD_MACSYFINDER=true ;;
                    8|defensefinder) DOWNLOAD_DEFENSEFINDER=true ;;
                    9|bakta)        DOWNLOAD_BAKTA=true ;;
                    10|bakta-light)  DOWNLOAD_BAKTA_LIGHT=true ;;
                    11|kofam)       DOWNLOAD_KOFAM=true ;;
                    12|eggnog)      DOWNLOAD_EGGNOG=true ;;
                    13|dbcan)       DOWNLOAD_DBCAN=true ;;
                    14|metaeuk)     DOWNLOAD_METAEUK=true ;;
                    15|kraken2)     DOWNLOAD_KRAKEN2=true ;;
                    16|silva)       DOWNLOAD_SILVA=true ;;
                    17|marferret)   DOWNLOAD_MARFERRET=true ;;
                    all)            DOWNLOAD_HUMAN=true; DOWNLOAD_GENOMAD=true; DOWNLOAD_CHECKV=true; DOWNLOAD_CHECKM=true; DOWNLOAD_CHECKM2=true; DOWNLOAD_KAIJU=true; DOWNLOAD_MACSYFINDER=true; DOWNLOAD_DEFENSEFINDER=true; DOWNLOAD_BAKTA=true; DOWNLOAD_BAKTA_LIGHT=true; DOWNLOAD_KOFAM=true; DOWNLOAD_EGGNOG=true; DOWNLOAD_DBCAN=true; DOWNLOAD_METAEUK=true; DOWNLOAD_KRAKEN2=true; DOWNLOAD_SILVA=true; DOWNLOAD_MARFERRET=true ;;
                    *) echo "[WARNING] Unknown selection: $item" >&2 ;;
                esac
            done
            ;;
    esac

    if ! $DOWNLOAD_HUMAN && ! $DOWNLOAD_GENOMAD && ! $DOWNLOAD_CHECKV && ! $DOWNLOAD_CHECKM && ! $DOWNLOAD_CHECKM2 && ! $DOWNLOAD_KAIJU && ! $DOWNLOAD_MACSYFINDER && ! $DOWNLOAD_DEFENSEFINDER && ! $DOWNLOAD_BAKTA && ! $DOWNLOAD_BAKTA_LIGHT && ! $DOWNLOAD_KOFAM && ! $DOWNLOAD_EGGNOG && ! $DOWNLOAD_DBCAN && ! $DOWNLOAD_METAEUK && ! $DOWNLOAD_KRAKEN2 && ! $DOWNLOAD_SILVA && ! $DOWNLOAD_MARFERRET; then
        echo "No databases selected. Exiting."
        exit 0
    fi
fi

# ============================================================================
# Download functions
# ============================================================================

mkdir -p "${DB_DIR}"

download_human_ref() {
    local db_path="${DB_DIR}/human_ref"
    if [ -d "${db_path}/ref" ] && [ -f "${db_path}/ref/genome/1/summary.txt" ]; then
        echo "[INFO] BBTools human reference already exists at ${db_path}"
        echo "  Delete ${db_path} and re-run to force re-download."
        return 0
    fi

    echo ""
    echo "[INFO] Downloading BBTools masked human reference (~4 GB)..."
    echo "  Destination: ${db_path}"

    mkdir -p "${db_path}"

    # Find removehuman.sh from the illumina_mag bbmap env, or fall back to PATH
    local removehuman_bin=""
    if [ -x "${ILLUMINA_ENV_DIR}/dana-illumina-mag-bbmap/bin/removehuman.sh" ]; then
        removehuman_bin="${ILLUMINA_ENV_DIR}/dana-illumina-mag-bbmap/bin/removehuman.sh"
    elif [ -x "${NANOPORE_ENV_DIR}/dana-mag-assembly/bin/removehuman.sh" ]; then
        removehuman_bin="${NANOPORE_ENV_DIR}/dana-mag-assembly/bin/removehuman.sh"
    elif command -v removehuman.sh &>/dev/null; then
        removehuman_bin="removehuman.sh"
    fi

    if [ -z "${removehuman_bin}" ]; then
        echo "[ERROR] removehuman.sh not found. Install illumina_mag conda envs first:" >&2
        echo "  cd illumina_mag/nextflow && ./install.sh" >&2
        return 1
    fi

    # removehuman.sh auto-downloads the masked hg19 reference on first run.
    # We run it with empty input to trigger the download without processing reads.
    echo "  Using: ${removehuman_bin}"
    echo "  This will download and index the masked human genome (hg19)..."
    "${removehuman_bin}" path="${db_path}" build=1 2>&1 | tail -5 || true

    # Verify the BBTools index was built
    if ! [ -d "${db_path}/ref" ] || ! [ -f "${db_path}/ref/genome/1/summary.txt" ]; then
        echo "[ERROR] Human reference download/build may have failed" >&2
        echo "  Check ${db_path} for partial files" >&2
        return 1
    fi
    echo "[SUCCESS] BBTools human reference downloaded to ${db_path}"

    # Also download GRCh38 FASTA for minimap2 (nanopore pipeline).
    # minimap2 will auto-index from the FASTA on first use.
    local human_fa="${db_path}/GRCh38_noalt_as.fa.gz"
    if [ -f "${human_fa}" ]; then
        echo "[INFO] Human FASTA already exists at ${human_fa}"
    else
        echo "  Downloading GRCh38 no-alt FASTA for minimap2 (~900 MB)..."
        local fasta_url="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz"
        if wget -q --show-progress -O "${human_fa}" "${fasta_url}" 2>&1; then
            echo "[SUCCESS] Human FASTA downloaded to ${human_fa}"
        else
            echo "[WARNING] Could not download GRCh38 FASTA — nanopore human removal unavailable" >&2
            rm -f "${human_fa}"
        fi
    fi

    echo "  illumina_mag: --human_ref ${db_path}"
    echo "  nanopore_mag: --human_ref ${human_fa}"
}

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

    if $USE_CONTAINER; then
        container_run genomad download-database /data/db || return 1
    else
        local genomad_bin="${ENV_DIR}/dana-mag-quality/bin/genomad"
        if [ ! -x "${genomad_bin}" ]; then
            echo "[ERROR] geNomad not installed. Run ./install.sh first or use --docker/--apptainer." >&2
            return 1
        fi
        "${genomad_bin}" download-database "${DB_DIR}" || return 1
    fi
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

    if $USE_CONTAINER; then
        container_run checkv download_database /data/db/checkv_db || return 1
    else
        local checkv_bin="${ENV_DIR}/dana-mag-quality/bin/checkv"
        if [ ! -x "${checkv_bin}" ]; then
            echo "[ERROR] CheckV not installed. Run ./install.sh first or use --docker/--apptainer." >&2
            return 1
        fi
        "${checkv_bin}" download_database "${db_path}" || return 1
    fi
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

download_checkm() {
    local db_path="${DB_DIR}/checkm_data"
    if [ -d "${db_path}" ] && [ -f "${db_path}/genome_tree/genome_tree.derep.txt" ] \
            && [ -s "${db_path}/genome_tree/genome_tree.derep.txt" ]; then
        echo "[INFO] CheckM v1 data already exists at ${db_path}"
        echo "  Delete ${db_path} and re-run to force re-download."
        return 0
    fi

    echo ""
    echo "[INFO] Downloading CheckM v1 reference data (~280 MB)..."
    echo "  Needed by dRep and COMEBin (transitive dependency on checkm-genome)"
    echo "  Destination: ${db_path}"

    mkdir -p "${db_path}"

    local url="https://zenodo.org/record/7401545/files/checkm_data_2015_01_16.tar.gz"
    local tarball="${db_path}/checkm_data_2015_01_16.tar.gz"
    local expected_sha="971ec469348bd6c3d9eb96142f567f12443310fa06c1892643940f35f86ac92c"

    wget -q --show-progress -O "${tarball}" "${url}"

    # Verify checksum
    local actual_sha
    actual_sha=$(sha256sum "${tarball}" | cut -d' ' -f1)
    if [ "${actual_sha}" != "${expected_sha}" ]; then
        echo "[ERROR] SHA256 mismatch for ${tarball}" >&2
        echo "  Expected: ${expected_sha}" >&2
        echo "  Got:      ${actual_sha}" >&2
        rm -f "${tarball}"
        return 1
    fi

    tar -xzf "${tarball}" -C "${db_path}"
    rm -f "${tarball}"

    # Point each env that has checkm-genome to the shared data directory
    for env_name in dana-mag-binning dana-mag-derep; do
        local env_path="${ENV_DIR}/${env_name}"
        local checkm_bin="${env_path}/bin/checkm"
        if [ -x "${checkm_bin}" ]; then
            echo "  Setting checkm data root for ${env_name}..."
            "${checkm_bin}" data setRoot "${db_path}" || true
        fi
    done

    echo "[SUCCESS] CheckM v1 data downloaded to ${db_path}"
    echo "  Use with: --checkm_data ${db_path}"
}

download_checkm2() {
    local db_path="${DB_DIR}/checkm2_db"
    if [ -d "${db_path}" ] && find -L "${db_path}" -maxdepth 2 -name '*.dmnd' 2>/dev/null | grep -q .; then
        echo "[INFO] CheckM2 database already exists at ${db_path}"
        echo "  Delete ${db_path} and re-run to force re-download."
        return 0
    fi

    echo ""
    echo "[INFO] Downloading CheckM2 database (~3.5 GB)..."
    echo "  Destination: ${db_path}"

    mkdir -p "${db_path}"
    if $USE_CONTAINER; then
        container_run checkm2 database --download --path /data/db/checkm2_db || return 1
    else
        local checkm2_bin="${ENV_DIR}/dana-mag-quality/bin/checkm2"
        if [ ! -x "${checkm2_bin}" ]; then
            echo "[ERROR] CheckM2 not installed. Run ./install.sh first or use --docker/--apptainer." >&2
            return 1
        fi
        "${checkm2_bin}" database --download --path "${db_path}" || return 1
    fi
    # Create a stable top-level symlink to the DIAMOND db file
    # (Binette and other tools need the .dmnd path, not the directory)
    local dmnd_file
    dmnd_file=$(find -L "${db_path}" -maxdepth 2 -name '*.dmnd' -type f | head -1)
    if [ -n "${dmnd_file}" ]; then
        ln -sf "$(realpath --relative-to="${db_path}" "${dmnd_file}")" "${db_path}/checkm2.dmnd"
    fi

    echo "[SUCCESS] CheckM2 database downloaded to ${db_path}"
    echo "  Use with: --checkm2_db ${db_path}"
}

download_kaiju() {
    local db_path="${DB_DIR}/kaiju_db"
    if [ -d "${db_path}" ] && find -L "${db_path}" -maxdepth 2 -name '*.fmi' 2>/dev/null | grep -q .; then
        echo "[INFO] Kaiju database already exists at ${db_path}"
        echo "  Delete ${db_path} and re-run to force re-download."
        return 0
    fi

    echo ""
    echo "[INFO] Downloading Kaiju pre-built RefSeq database..."
    echo "  Destination: ${db_path}"

    mkdir -p "${db_path}"

    # Download pre-built index from kaiju server (NOT kaiju-makedb which rebuilds from GenBank)
    local base_url="https://kaiju-idx.s3.eu-central-1.amazonaws.com/2024"
    local tarball="kaiju_db_refseq_ref_2024-08-14.tgz"
    if [ ! -f "${db_path}/${tarball}" ]; then
        echo "  Downloading pre-built index: ${tarball}"
        wget -c -q --show-progress -O "${db_path}/${tarball}" \
            "${base_url}/${tarball}" || return 1
    fi

    echo "[INFO] Extracting ${tarball}..."
    tar xzf "${db_path}/${tarball}" -C "${db_path}" || return 1
    rm -f "${db_path}/${tarball}"

    echo "[SUCCESS] Kaiju database downloaded to ${db_path}"
    echo "  Use with: --kaiju_db ${db_path}"
}

download_macsyfinder() {
    local db_path="${DB_DIR}/macsyfinder_models"
    if [ -d "${db_path}/TXSScan" ] && [ -d "${db_path}/CONJScan" ]; then
        echo "[INFO] MacSyFinder models already exist at ${db_path}"
        echo "  Delete ${db_path} and re-run to force re-download."
        return 0
    fi

    echo ""
    echo "[INFO] Downloading MacSyFinder models (~50 MB)..."
    echo "  Models: TXSScan (20 secretion systems), CONJScan (17 conjugation systems)"
    echo "  Destination: ${db_path}"

    mkdir -p "${db_path}"
    if $USE_CONTAINER; then
        container_run macsydata install --target /data/db/macsyfinder_models TXSScan || return 1
        container_run macsydata install --target /data/db/macsyfinder_models CONJScan || return 1
    else
        local macsydata_bin="${ENV_DIR}/dana-mag-genomic/bin/macsydata"
        if [ ! -x "${macsydata_bin}" ]; then
            echo "[ERROR] MacSyFinder not installed. Run ./install.sh first or use --docker/--apptainer." >&2
            return 1
        fi
        "${macsydata_bin}" install --target "${db_path}" TXSScan || return 1
        "${macsydata_bin}" install --target "${db_path}" CONJScan || return 1
    fi
    echo "[SUCCESS] MacSyFinder models downloaded to ${db_path}"
    echo "  Use with: --macsyfinder_models ${db_path}"
}

download_defensefinder() {
    local db_path="${DB_DIR}/defensefinder_models"
    if [ -d "${db_path}/defense-finder-models" ] && [ -d "${db_path}/CasFinder" ]; then
        echo "[INFO] DefenseFinder models already exist at ${db_path}"
        echo "  Delete ${db_path} and re-run to force re-download."
        return 0
    fi

    echo ""
    echo "[INFO] Downloading DefenseFinder models (~100 MB)..."
    echo "  Models: ~280 anti-phage defense system HMM profiles"
    echo "  Destination: ${db_path}"

    mkdir -p "${db_path}"
    if $USE_CONTAINER; then
        container_run defense-finder update --models-dir /data/db/defensefinder_models || return 1
        # Workaround: CasFinder 3.1.1 model version incompatibility
        local cf_ver
        cf_ver=$(cat "${db_path}/CasFinder/metadata.yml" 2>/dev/null | grep -oP 'vers: \K.*' | head -1)
        if [ "${cf_ver}" = "3.1.1" ]; then
            echo "  Downgrading CasFinder 3.1.1 → 3.1.0 (model version compatibility fix)"
            rm -rf "${db_path}/CasFinder"
            container_run msf_data install --target /data/db/defensefinder_models CasFinder==3.1.0 || return 1
        fi
    else
        local df_bin="${ENV_DIR}/dana-mag-genomic/bin/defense-finder"
        if [ ! -x "${df_bin}" ]; then
            echo "[ERROR] DefenseFinder not installed. Run ./install.sh first or use --docker/--apptainer." >&2
            return 1
        fi
        "${df_bin}" update --models-dir "${db_path}"
        # Workaround: CasFinder 3.1.1 has model definition version 2.1 which is
        # incompatible with macsyfinder 2.1.4 bundled with defense-finder 2.0.1
        # (https://github.com/mdmparis/defense-finder/issues/91)
        # Downgrade to 3.1.0 which uses compatible model definition version 2.0
        local msf_data_bin="${ENV_DIR}/dana-mag-genomic/bin/macsydata"
        if [ -x "${msf_data_bin}" ]; then
            local cf_ver
            cf_ver=$(cat "${db_path}/CasFinder/metadata.yml" 2>/dev/null | grep -oP 'vers: \K.*' | head -1)
            if [ "${cf_ver}" = "3.1.1" ]; then
                echo "  Downgrading CasFinder 3.1.1 → 3.1.0 (model version compatibility fix)"
                rm -rf "${db_path}/CasFinder"
                "${msf_data_bin}" install --target "${db_path}" CasFinder==3.1.0
            fi
        fi
    fi
    echo "[SUCCESS] DefenseFinder models downloaded to ${db_path}"
    echo "  Use with: --defensefinder_models ${db_path}"
}

download_bakta() {
    local db_path="${DB_DIR}/bakta"
    if [ -d "${db_path}/db" ] && [ -f "${db_path}/db/version.json" ]; then
        echo "[INFO] Bakta full database already exists at ${db_path}/db"
        echo "  Delete ${db_path}/db and re-run to force re-download."
        return 0
    fi

    echo ""
    echo "[INFO] Downloading Bakta full database (~37 GB)..."
    echo "  Destination: ${db_path}/db"

    mkdir -p "${db_path}"
    if $USE_CONTAINER; then
        container_run bakta_db download --output /data/db/bakta --type full || return 1
    else
        local bakta_env="${ENV_DIR}/dana-mag-annotate"
        if [ ! -x "${bakta_env}/bin/bakta_db" ]; then
            echo "[ERROR] Bakta not installed. Run ./install.sh first or use --docker/--apptainer." >&2
            return 1
        fi
        # bakta_db requires AMRFinderPlus on PATH, so run inside the full conda env
        ${CONDA_CMD} run -p "${bakta_env}" bakta_db download --output "${db_path}" --type full || return 1
    fi
    echo "[SUCCESS] Bakta full database downloaded to ${db_path}/db"
    echo "  Use with: --bakta_db ${db_path}/db"
}

download_bakta_light() {
    local db_path="${DB_DIR}/bakta"
    if [ -d "${db_path}/db-light" ] && [ -f "${db_path}/db-light/version.json" ]; then
        echo "[INFO] Bakta light database already exists at ${db_path}/db-light"
        echo "  Delete ${db_path}/db-light and re-run to force re-download."
        return 0
    fi

    echo ""
    echo "[INFO] Downloading Bakta light database (~1.4 GB)..."
    echo "  Destination: ${db_path}/db-light"

    mkdir -p "${db_path}"
    if $USE_CONTAINER; then
        container_run bakta_db download --output /data/db/bakta --type light || return 1
    else
        local bakta_env="${ENV_DIR}/dana-mag-annotate"
        if [ ! -x "${bakta_env}/bin/bakta_db" ]; then
            echo "[ERROR] Bakta not installed. Run ./install.sh first or use --docker/--apptainer." >&2
            return 1
        fi
        # bakta_db requires AMRFinderPlus on PATH, so run inside the full conda env
        ${CONDA_CMD} run -p "${bakta_env}" bakta_db download --output "${db_path}" --type light || return 1
    fi
    echo "[SUCCESS] Bakta light database downloaded to ${db_path}/db-light"
    echo "  Use with: --bakta_light_db ${db_path}/db-light"
}

download_kofam() {
    local db_path="${DB_DIR}/kofam_db"
    if [ -d "${db_path}/profiles" ] && [ -f "${db_path}/ko_list" ]; then
        echo "[INFO] KOfam database already exists at ${db_path}"
        echo "  Delete ${db_path} and re-run to force re-download."
        return 0
    fi

    echo ""
    echo "[INFO] Downloading KOfam profiles (~4 GB)..."
    echo "  Destination: ${db_path}"

    mkdir -p "${db_path}"

    # Download ko_list (adaptive score thresholds)
    wget -q -O "${db_path}/ko_list.gz" "ftp://ftp.genome.jp/pub/db/kofam/ko_list.gz"
    gunzip -f "${db_path}/ko_list.gz"

    # Download HMM profiles
    wget -q -O "${db_path}/profiles.tar.gz" "ftp://ftp.genome.jp/pub/db/kofam/profiles.tar.gz"
    tar -xzf "${db_path}/profiles.tar.gz" -C "${db_path}"
    rm -f "${db_path}/profiles.tar.gz"

    echo "[SUCCESS] KOfam database downloaded to ${db_path}"
    echo "  Use with: --kofam_db ${db_path}"
}

download_eggnog() {
    local db_path="${DB_DIR}/eggnog_db"
    if [ -d "${db_path}" ] && [ -f "${db_path}/eggnog.db" ]; then
        echo "[INFO] eggNOG database already exists at ${db_path}"
        echo "  Delete ${db_path} and re-run to force re-download."
        return 0
    fi

    echo ""
    echo "[INFO] Downloading eggNOG-mapper database (~12 GB)..."
    echo "  Destination: ${db_path}"

    mkdir -p "${db_path}"
    if $USE_CONTAINER; then
        # Workaround: eggnog-mapper <=2.1.13 has broken download URLs
        # (eggnogdb.embl.de returns 404; eggnog5.embl.de is the working host)
        # Must patch all files — URL is in common.py, db_sqlite.py, and the download script
        container_run sh -c "find /opt/conda/envs/dana-mag-annotate -name '*.py' -exec grep -l 'eggnogdb.embl.de' {} \; | xargs sed -i 's|eggnogdb.embl.de|eggnog5.embl.de|g'; download_eggnog_data.py --data_dir /data/db/eggnog_db -y" || return 1
    else
        local emapper_bin="${ENV_DIR}/dana-mag-annotate/bin/download_eggnog_data.py"
        if [ ! -x "${emapper_bin}" ]; then
            echo "[ERROR] eggNOG-mapper not installed. Run ./install.sh first or use --docker/--apptainer." >&2
            return 1
        fi
        # Workaround: eggnog-mapper <=2.1.13 has broken download URLs
        # (eggnogdb.embl.de returns 404; eggnog5.embl.de is the working host)
        # https://github.com/eggnogdb/eggnog-mapper/issues/571
        sed -i "s|eggnogdb.embl.de|eggnog5.embl.de|g" "${emapper_bin}"
        "${emapper_bin}" --data_dir "${db_path}" -y || return 1
    fi
    echo "[SUCCESS] eggNOG database downloaded to ${db_path}"
    echo "  Use with: --eggnog_db ${db_path}"
}

download_dbcan() {
    local db_path="${DB_DIR}/dbcan_db"
    if [ -d "${db_path}" ] && [ -f "${db_path}/dbCAN.hmm" ] && [ -f "${db_path}/dbCAN.hmm.h3i" ]; then
        echo "[INFO] dbCAN database already exists at ${db_path}"
        echo "  Delete ${db_path} and re-run to force re-download."
        return 0
    fi

    echo ""
    echo "[INFO] Downloading dbCAN databases (~5 GB)..."
    echo "  Destination: ${db_path}"

    mkdir -p "${db_path}"
    if $USE_CONTAINER; then
        # Use AWS S3 mirror — bcb.unl.edu has persistent SSL cert issues
        container_run run_dbcan database --db_dir /data/db/dbcan_db --aws_s3 || return 1
        # Press HMM databases for HMMER (required before first run_dbcan use)
        echo "[INFO] Pressing HMM databases with hmmpress..."
        for hmm_name in dbCAN.hmm dbCAN-sub.hmm STP.hmm TF.hmm; do
            if [ -f "${db_path}/${hmm_name}" ] && [ ! -f "${db_path}/${hmm_name}.h3i" ]; then
                echo "  Pressing ${hmm_name}..."
                container_run /opt/conda/envs/dana-mag-genomic/bin/hmmpress "/data/db/dbcan_db/${hmm_name}" || return 1
            fi
        done
    else
        local run_dbcan="${ENV_DIR}/dana-mag-genomic/bin/run_dbcan"
        local hmmpress="${ENV_DIR}/dana-mag-genomic/bin/hmmpress"
        if [ ! -x "${run_dbcan}" ]; then
            echo "[ERROR] dbCAN not installed. Run ./install.sh first or use --docker/--apptainer." >&2
            return 1
        fi
        # Use AWS S3 mirror — bcb.unl.edu has persistent SSL cert issues
        # (https://github.com/bcb-unl/run_dbcan/issues/53)
        "${run_dbcan}" database --db_dir "${db_path}" --aws_s3 || return 1
        # Press HMM databases for HMMER (required before first run_dbcan use)
        echo "[INFO] Pressing HMM databases with hmmpress..."
        for hmm in "${db_path}"/dbCAN.hmm "${db_path}"/dbCAN-sub.hmm "${db_path}"/STP.hmm "${db_path}"/TF.hmm; do
            if [ -f "${hmm}" ] && [ ! -f "${hmm}.h3i" ]; then
                echo "  Pressing $(basename "${hmm}")..."
                "${hmmpress}" "${hmm}" || return 1
            fi
        done
    fi

    echo "[SUCCESS] dbCAN database downloaded to ${db_path}"
    echo "  Use with: --dbcan_db ${db_path}"
}

download_metaeuk() {
    local db_path="${DB_DIR}/metaeuk_db"
    if [ -f "${db_path}/metaeuk_db" ] && [ -f "${db_path}/metaeuk_db.dbtype" ]; then
        echo "[INFO] MetaEuk database already exists at ${db_path}"
        echo "  Delete ${db_path} and re-run to force re-download."
        return 0
    fi

    echo ""
    echo "[INFO] Downloading MetaEuk OrthoDB v11 Eukaryota database (~8.5 GB compressed)..."
    echo "  Source: bioinf.uni-greifswald.de (OrthoDB v11, ~2,000 eukaryotic genomes)"
    echo "  Destination: ${db_path}"
    echo ""
    echo "  Note: OrthoDB v11 is used instead of v12 (~23 GB) because v12 decompresses"
    echo "  to ~70-80 GB and exceeds 64 GB RAM even with --split-memory-limit."
    echo "  For systems with >128 GB RAM, download v12 manually from:"
    echo "  https://bioinf.uni-greifswald.de/bioinf/partitioned_odb12/Eukaryota.fa.gz"

    mkdir -p "${db_path}"

    # Download OrthoDB v11 Eukaryota protein FASTA (~2,000 eukaryotic genomes)
    local fasta_gz="${db_path}/Eukaryota.fa.gz"
    local fasta="${db_path}/Eukaryota.fa"
    if [ ! -f "${fasta}" ]; then
        echo "[INFO] Downloading Eukaryota.fa.gz..."
        wget -q --show-progress -O "${fasta_gz}" \
            "https://bioinf.uni-greifswald.de/bioinf/partitioned_odb11/Eukaryota.fa.gz"
        echo "[INFO] Decompressing..."
        gunzip -f "${fasta_gz}"
    else
        echo "[INFO] FASTA already downloaded: ${fasta}"
    fi

    # Convert to MMseqs2 database format
    echo "[INFO] Creating MMseqs2 database (metaeuk createdb)..."
    if $USE_CONTAINER; then
        container_run metaeuk createdb /data/db/metaeuk_db/Eukaryota.fa /data/db/metaeuk_db/metaeuk_db || return 1
    else
        local metaeuk_bin="${ENV_DIR}/dana-mag-classify/bin/metaeuk"
        if [ ! -x "${metaeuk_bin}" ]; then
            echo "[ERROR] MetaEuk not installed. Run ./install.sh first or use --docker/--apptainer." >&2
            return 1
        fi
        "${metaeuk_bin}" createdb "${fasta}" "${db_path}/metaeuk_db" || return 1
    fi

    echo "[SUCCESS] MetaEuk database created at ${db_path}/metaeuk_db"
    echo "  Use with: --metaeuk_db ${db_path}/metaeuk_db"
}

download_kraken2() {
    local db_path="${DB_DIR}/kraken2_db"
    if [ -d "${db_path}" ] && [ -f "${db_path}/hash.k2d" ]; then
        echo "[INFO] Kraken2 database already exists at ${db_path}"
        echo "  Delete ${db_path} and re-run to force re-download."
        return 0
    fi

    echo ""
    echo "[INFO] Downloading Kraken2 PlusPFP-8 database (~8 GB)..."
    echo "  Source: genome-idx S3 index (bacteria, archaea, viruses, fungi, protozoa, UniVec_Core)"
    echo "  Destination: ${db_path}"

    mkdir -p "${db_path}"

    local url="https://genome-idx.s3.amazonaws.com/kraken/k2_pluspfp_08gb_20240605.tar.gz"
    local tarball="${db_path}/k2_pluspfp_08gb.tar.gz"

    wget -q --show-progress -O "${tarball}" "${url}"
    tar -xzf "${tarball}" -C "${db_path}"
    rm -f "${tarball}"

    echo "[SUCCESS] Kraken2 database downloaded to ${db_path}"
    echo "  Use with: --kraken2_db ${db_path}"
}

download_silva() {
    local db_path="${DB_DIR}/silva_db"
    local ssu_fasta="${db_path}/SILVA_138.2_SSURef_NR99.fasta"
    local lsu_fasta="${db_path}/SILVA_138.2_LSURef_NR99.fasta"

    if [ -f "${ssu_fasta}" ] && [ -f "${lsu_fasta}" ]; then
        echo "[INFO] SILVA databases already exist at ${db_path}"
        echo "  Delete ${db_path} and re-run to force re-download."
        return 0
    fi

    echo ""
    echo "[INFO] Downloading SILVA 138.2 SSU + LSU NR99 databases (~900 MB)..."
    echo "  Source: arb-silva.de (SILVA 138.2, non-redundant 99% identity)"
    echo "  Destination: ${db_path}"

    mkdir -p "${db_path}"

    # Download SSU NR99 (~600 MB compressed)
    if [ ! -f "${ssu_fasta}" ]; then
        local ssu_gz="${db_path}/SILVA_138.2_SSURef_NR99_tax_silva.fasta.gz"
        echo "[INFO] Downloading SILVA SSU NR99..."
        wget -q --show-progress -O "${ssu_gz}" \
            "https://www.arb-silva.de/fileadmin/silva_databases/release_138.2/Exports/SILVA_138.2_SSURef_NR99_tax_silva.fasta.gz"
        echo "[INFO] Decompressing and converting U→T for vsearch compatibility..."
        gunzip -c "${ssu_gz}" | sed '/^[^>]/s/U/T/g' > "${ssu_fasta}"
        rm -f "${ssu_gz}"
        echo "[INFO] SSU NR99 ready: ${ssu_fasta}"
    fi

    # Download LSU NR99 (~300 MB compressed)
    if [ ! -f "${lsu_fasta}" ]; then
        local lsu_gz="${db_path}/SILVA_138.2_LSURef_NR99_tax_silva.fasta.gz"
        echo "[INFO] Downloading SILVA LSU NR99..."
        wget -q --show-progress -O "${lsu_gz}" \
            "https://www.arb-silva.de/fileadmin/silva_databases/release_138.2/Exports/SILVA_138.2_LSURef_NR99_tax_silva.fasta.gz"
        echo "[INFO] Decompressing and converting U→T for vsearch compatibility..."
        gunzip -c "${lsu_gz}" | sed '/^[^>]/s/U/T/g' > "${lsu_fasta}"
        rm -f "${lsu_gz}"
        echo "[INFO] LSU NR99 ready: ${lsu_fasta}"
    fi

    echo "[SUCCESS] SILVA databases downloaded to ${db_path}"
    echo "  Use with: --silva_ssu_db ${ssu_fasta}"
    echo "            --silva_lsu_db ${lsu_fasta}"
}

download_marferret() {
    local db_path="${DB_DIR}/marferret_db"
    if [ -d "${db_path}" ] && ls "${db_path}"/*.dmnd 1>/dev/null 2>&1; then
        echo "[INFO] MarFERReT database already exists at ${db_path}"
        echo "  Delete ${db_path} and re-run to force re-download."
        return 0
    fi

    echo ""
    echo "[INFO] Downloading MarFERReT v1.1.1 marine eukaryotic database (~9 GB total)..."
    echo "  Source: Zenodo (Carradec et al. 2023, Scientific Data 10:901)"
    echo "  Contents: ~28M protein sequences from 800 marine eukaryotic taxa"
    echo "  Destination: ${db_path}"

    mkdir -p "${db_path}"

    local zenodo_base="https://zenodo.org/records/10553848/files"

    # Pre-built DIAMOND database (~8.9 GB)
    if [ ! -f "${db_path}/MarFERReT.v1.1.1.dmnd" ]; then
        echo "[INFO] Downloading MarFERReT.v1.1.1.dmnd (~8.9 GB)..."
        wget -q --show-progress -O "${db_path}/MarFERReT.v1.1.1.dmnd" \
            "${zenodo_base}/MarFERReT.v1.1.1.dmnd"
    else
        echo "[INFO] DIAMOND database already present"
    fi

    # Taxonomy mapping (~90.6 MB)
    if [ ! -f "${db_path}/MarFERReT.v1.1.1.taxonomies.tab.gz" ]; then
        echo "[INFO] Downloading MarFERReT.v1.1.1.taxonomies.tab.gz (~90.6 MB)..."
        wget -q --show-progress -O "${db_path}/MarFERReT.v1.1.1.taxonomies.tab.gz" \
            "${zenodo_base}/MarFERReT.v1.1.1.taxonomies.tab.gz"
    else
        echo "[INFO] Taxonomy file already present"
    fi

    # Pfam annotations (~121.9 MB)
    if [ ! -f "${db_path}/MarFERReT.v1.1.1.best_pfam_annotations.csv.gz" ]; then
        echo "[INFO] Downloading MarFERReT.v1.1.1.best_pfam_annotations.csv.gz (~121.9 MB)..."
        wget -q --show-progress -O "${db_path}/MarFERReT.v1.1.1.best_pfam_annotations.csv.gz" \
            "${zenodo_base}/MarFERReT.v1.1.1.best_pfam_annotations.csv.gz"
    else
        echo "[INFO] Pfam annotation file already present"
    fi

    echo "[SUCCESS] MarFERReT database downloaded to ${db_path}"
    echo "  Use with: --marferret_db ${db_path}"
}

download_gtdbtk() {
    local db_path="${DB_DIR}/gtdbtk_db"
    if [ -d "${db_path}/taxonomy" ] && [ -d "${db_path}/markers" ]; then
        echo "[INFO] GTDB-Tk database already exists at ${db_path}"
        echo "  Delete ${db_path} and re-run to force re-download."
        return 0
    fi

    echo ""
    echo "[INFO] Downloading GTDB-Tk r226 reference data..."
    echo "  Destination: ${db_path}"
    echo ""
    echo "  The full package is ~132 GB (~85 GB extracted)."
    echo "  The skani/ directory (~100 GB compressed) is required for species-level"
    echo "  ANI classification even when using --skip_ani_screen."
    echo ""
    echo "  Mirrors (Denmark is fastest worldwide):"
    echo "    EU:  https://data.gtdb.aau.ecogenomic.org/releases/"
    echo "    AU1: https://data.ace.uq.edu.au/public/gtdb/data/releases/"
    echo "    AU2: https://data.gtdb.ecogenomic.org/releases/"

    mkdir -p "${db_path}"

    local tarball="gtdbtk_r226_data.tar.gz"
    local rel_path="release226/226.0/auxillary_files/gtdbtk_package/full_package/${tarball}"
    # Australian mirrors are fastest from most locations; Denmark as fallback
    local urls=(
        "https://data.gtdb.ecogenomic.org/releases/${rel_path}"
        "https://data.ace.uq.edu.au/public/gtdb/data/releases/${rel_path}"
        "https://data.gtdb.aau.ecogenomic.org/releases/${rel_path}"
    )

    if [ ! -f "${db_path}/${tarball}" ]; then
        echo "[INFO] Downloading ${tarball} (~132 GB)..."
        local downloaded=false
        if command -v aria2c &>/dev/null; then
            # aria2 supports multi-source parallel download + auto-resume
            echo "  Using aria2c (multi-source parallel download)"
            aria2c -x 4 -s 4 -c -d "${db_path}" -o "${tarball}" \
                "${urls[@]}" && downloaded=true
        fi
        if ! $downloaded; then
            # wget fallback — try mirrors in order, -c for resume
            for url in "${urls[@]}"; do
                echo "  Trying: ${url}"
                if wget -c -q --show-progress -O "${db_path}/${tarball}" "${url}"; then
                    downloaded=true
                    break
                fi
                echo "  Mirror failed, trying next..."
            done
        fi
        $downloaded || { echo "[ERROR] All download mirrors failed" >&2; return 1; }
    else
        echo "[INFO] Tarball already present, skipping download"
    fi

    echo "[INFO] Extracting ${tarball}..."
    tar xzf "${db_path}/${tarball}" -C "${db_path}" --strip-components=1 || return 1

    echo "[INFO] Removing tarball to free disk space..."
    rm -f "${db_path}/${tarball}"

    echo "[SUCCESS] GTDB-Tk database downloaded to ${db_path}"
    echo "  Use with: --gtdbtk_db ${db_path}"
}

download_antismash() {
    local db_path="${DB_DIR}/antismash_db"
    if [ -d "${db_path}" ] && [ -f "${db_path}/clusterblast/proteins.dmnd" ]; then
        echo "[INFO] antiSMASH database already exists at ${db_path}"
        echo "  Delete ${db_path} and re-run to force re-download."
        return 0
    fi

    echo ""
    echo "[INFO] Downloading antiSMASH databases..."
    echo "  Destination: ${db_path}"

    local antismash_bin="${ENV_DIR}/dana-mag-antismash/bin/download-antismash-databases"
    if [ ! -x "$antismash_bin" ]; then
        echo "[ERROR] antiSMASH conda env not found. Install it first:" >&2
        echo "  mamba create -p ${ENV_DIR}/dana-mag-antismash -y -c bioconda -c conda-forge --override-channels --channel-priority flexible 'antismash>=8.0'" >&2
        return 1
    fi

    "$antismash_bin" --database-dir "${db_path}" || return 1

    echo "[SUCCESS] antiSMASH databases downloaded to ${db_path}"
    echo "  Use with: --antismash_db ${db_path}"
}

download_magscot() {
    local db_path="${DB_DIR}/magscot_hmm"
    if [ -d "${db_path}" ] && [ -f "${db_path}/gtdbtk_rel207_tigrfam.hmm" ] && [ -f "${db_path}/gtdbtk_rel207_Pfam-A.hmm" ]; then
        echo "[INFO] MAGScoT HMMs already exist at ${db_path}"
        echo "  Delete ${db_path} and re-run to force re-download."
        return 0
    fi

    echo ""
    echo "[INFO] Downloading MAGScoT GTDB rel207 marker HMMs..."
    echo "  Destination: ${db_path}"

    mkdir -p "${db_path}"

    local base_url="https://raw.githubusercontent.com/ikmb/MAGScoT/main/hmm"
    curl -fSL "${base_url}/gtdbtk_rel207_tigrfam.hmm" -o "${db_path}/gtdbtk_rel207_tigrfam.hmm" || return 1
    curl -fSL "${base_url}/gtdbtk_rel207_Pfam-A.hmm"  -o "${db_path}/gtdbtk_rel207_Pfam-A.hmm"  || return 1

    echo "[SUCCESS] MAGScoT HMMs downloaded to ${db_path}"
    echo "  Use with: --magscot_hmm_dir ${db_path}"
}

# ============================================================================
# Execute downloads
# ============================================================================

echo ""
echo "Database directory: ${DB_DIR}"

failed=0

if $DOWNLOAD_HUMAN; then
    ( download_human_ref ) || failed=$((failed + 1))
fi

if $DOWNLOAD_GENOMAD; then
    ( download_genomad ) || failed=$((failed + 1))
fi

if $DOWNLOAD_CHECKV; then
    ( download_checkv ) || failed=$((failed + 1))
fi

if $DOWNLOAD_CHECKM; then
    ( download_checkm ) || failed=$((failed + 1))
fi

if $DOWNLOAD_CHECKM2; then
    ( download_checkm2 ) || failed=$((failed + 1))
fi

if $DOWNLOAD_KAIJU; then
    ( download_kaiju ) || failed=$((failed + 1))
fi

if $DOWNLOAD_MACSYFINDER; then
    ( download_macsyfinder ) || failed=$((failed + 1))
fi

if $DOWNLOAD_DEFENSEFINDER; then
    ( download_defensefinder ) || failed=$((failed + 1))
fi

if $DOWNLOAD_BAKTA; then
    ( download_bakta ) || failed=$((failed + 1))
fi

if $DOWNLOAD_BAKTA_LIGHT; then
    ( download_bakta_light ) || failed=$((failed + 1))
fi

if $DOWNLOAD_KOFAM; then
    ( download_kofam ) || failed=$((failed + 1))
fi

if $DOWNLOAD_EGGNOG; then
    ( download_eggnog ) || failed=$((failed + 1))
fi

if $DOWNLOAD_DBCAN; then
    ( download_dbcan ) || failed=$((failed + 1))
fi

if $DOWNLOAD_METAEUK; then
    ( download_metaeuk ) || failed=$((failed + 1))
fi

if $DOWNLOAD_KRAKEN2; then
    ( download_kraken2 ) || failed=$((failed + 1))
fi

if $DOWNLOAD_SILVA; then
    ( download_silva ) || failed=$((failed + 1))
fi

if $DOWNLOAD_MARFERRET; then
    ( download_marferret ) || failed=$((failed + 1))
fi

if $DOWNLOAD_GTDBTK; then
    ( download_gtdbtk ) || failed=$((failed + 1))
fi

if $DOWNLOAD_ANTISMASH; then
    ( download_antismash ) || failed=$((failed + 1))
fi

if $DOWNLOAD_MAGSCOT; then
    ( download_magscot ) || failed=$((failed + 1))
fi

echo ""
if (( failed > 0 )); then
    echo "[WARNING] ${failed} database download(s) failed"
    exit 1
else
    echo "[SUCCESS] All selected databases downloaded to: ${DB_DIR}"
    echo ""
    echo "Run the pipeline with database paths:"
    ! $DOWNLOAD_HUMAN  || echo "  --human_ref  ${DB_DIR}/human_ref   (illumina_mag)"
    ! $DOWNLOAD_GENOMAD || echo "  --genomad_db ${DB_DIR}/genomad_db"
    ! $DOWNLOAD_CHECKV  || echo "  --checkv_db  ${DB_DIR}/checkv_db"
    ! $DOWNLOAD_CHECKM  || echo "  --checkm_data ${DB_DIR}/checkm_data"
    ! $DOWNLOAD_CHECKM2 || echo "  --checkm2_db ${DB_DIR}/checkm2_db"
    ! $DOWNLOAD_KAIJU   || echo "  --kaiju_db   ${DB_DIR}/kaiju_db"
    ! $DOWNLOAD_MACSYFINDER || echo "  --macsyfinder_models ${DB_DIR}/macsyfinder_models"
    ! $DOWNLOAD_DEFENSEFINDER || echo "  --defensefinder_models ${DB_DIR}/defensefinder_models"
    ! $DOWNLOAD_BAKTA         || echo "  --bakta_db         ${DB_DIR}/bakta/db"
    ! $DOWNLOAD_BAKTA_LIGHT   || echo "  --bakta_light_db   ${DB_DIR}/bakta/db-light"
    ! $DOWNLOAD_KOFAM   || echo "  --kofam_db   ${DB_DIR}/kofam_db"
    ! $DOWNLOAD_EGGNOG  || echo "  --eggnog_db  ${DB_DIR}/eggnog_db"
    ! $DOWNLOAD_DBCAN   || echo "  --dbcan_db   ${DB_DIR}/dbcan_db"
    ! $DOWNLOAD_METAEUK || echo "  --metaeuk_db ${DB_DIR}/metaeuk_db/metaeuk_db"
    ! $DOWNLOAD_KRAKEN2 || echo "  --kraken2_db ${DB_DIR}/kraken2_db"
    ! $DOWNLOAD_SILVA   || echo "  --silva_ssu_db ${DB_DIR}/silva_db/SILVA_138.2_SSURef_NR99.fasta"
    ! $DOWNLOAD_SILVA   || echo "  --silva_lsu_db ${DB_DIR}/silva_db/SILVA_138.2_LSURef_NR99.fasta"
    ! $DOWNLOAD_MARFERRET || echo "  --marferret_db ${DB_DIR}/marferret_db"
    ! $DOWNLOAD_GTDBTK   || echo "  --gtdbtk_db  ${DB_DIR}/gtdbtk_db"
    ! $DOWNLOAD_ANTISMASH || echo "  --antismash_db ${DB_DIR}/antismash_db"
    ! $DOWNLOAD_MAGSCOT   || echo "  --magscot_hmm_dir ${DB_DIR}/magscot_hmm"
fi
