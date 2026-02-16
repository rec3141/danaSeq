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
#   ./download-databases.sh --macsyfinder      # Download MacSyFinder models (TXSScan + CONJScan)
#   ./download-databases.sh --defensefinder    # Download DefenseFinder models (~100 MB)
#   ./download-databases.sh --bakta            # Download Bakta annotation database (~37 GB)
#   ./download-databases.sh --kofam            # Download KOfam profiles (~4 GB)
#   ./download-databases.sh --eggnog           # Download eggNOG diamond database (~12 GB)
#   ./download-databases.sh --dbcan            # Download dbCAN databases (~2 GB)
#   ./download-databases.sh --metaeuk          # Download MetaEuk OrthoDB Eukaryota (~23 GB)
#   ./download-databases.sh --kraken2         # Download Kraken2 PlusPFP-8 (~8 GB)
#   ./download-databases.sh --silva           # Download SILVA SSU + LSU NR99 (~900 MB)
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
DOWNLOAD_MACSYFINDER=false
DOWNLOAD_DEFENSEFINDER=false
DOWNLOAD_BAKTA=false
DOWNLOAD_KOFAM=false
DOWNLOAD_EGGNOG=false
DOWNLOAD_DBCAN=false
DOWNLOAD_METAEUK=false
DOWNLOAD_KRAKEN2=false
DOWNLOAD_SILVA=false
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
        --macsyfinder) DOWNLOAD_MACSYFINDER=true; INTERACTIVE=false; shift ;;
        --defensefinder) DOWNLOAD_DEFENSEFINDER=true; INTERACTIVE=false; shift ;;
        --bakta)     DOWNLOAD_BAKTA=true; INTERACTIVE=false; shift ;;
        --kofam)     DOWNLOAD_KOFAM=true; INTERACTIVE=false; shift ;;
        --eggnog)    DOWNLOAD_EGGNOG=true; INTERACTIVE=false; shift ;;
        --dbcan)     DOWNLOAD_DBCAN=true; INTERACTIVE=false; shift ;;
        --metaeuk)   DOWNLOAD_METAEUK=true; INTERACTIVE=false; shift ;;
        --kraken2)   DOWNLOAD_KRAKEN2=true; INTERACTIVE=false; shift ;;
        --silva)     DOWNLOAD_SILVA=true; INTERACTIVE=false; shift ;;
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
    DOWNLOAD_MACSYFINDER=true
    DOWNLOAD_DEFENSEFINDER=true
    DOWNLOAD_BAKTA=true
    DOWNLOAD_KOFAM=true
    DOWNLOAD_EGGNOG=true
    DOWNLOAD_DBCAN=true
    DOWNLOAD_METAEUK=true
    DOWNLOAD_KRAKEN2=true
    DOWNLOAD_SILVA=true
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
    printf "  %-12s %-8s  %s\n" "macsyfinder" "~50 MB" "MacSyFinder models: TXSScan + CONJScan (secretion + conjugation)"
    printf "  %-12s %-8s  %s\n" "defensefinder" "~100 MB" "DefenseFinder models: ~280 defense system HMM profiles"
    printf "  %-12s %-8s  %s\n" "bakta"    "~37 GB"  "Bakta annotation db (UniProt, AMRFinderPlus, Pfam, etc.)"
    printf "  %-12s %-8s  %s\n" "kofam"    "~4 GB"   "KOfam profiles + ko_list (KEGG Orthology via HMM)"
    printf "  %-12s %-8s  %s\n" "eggnog"   "~12 GB"  "eggNOG-mapper DIAMOND db (COG/GO/EC/KEGG/Pfam)"
    printf "  %-12s %-8s  %s\n" "dbcan"    "~2 GB"   "dbCAN HMM + DIAMOND + substrate db (CAZyme annotation)"
    printf "  %-12s %-8s  %s\n" "metaeuk"  "~8.5 GB" "MetaEuk OrthoDB v11 Eukaryota (eukaryotic gene prediction)"
    printf "  %-12s %-8s  %s\n" "kraken2"  "~8 GB"   "Kraken2 PlusPFP-8 (k-mer contig-level taxonomy)"
    printf "  %-12s %-8s  %s\n" "silva"    "~900 MB" "SILVA 138.2 SSU + LSU NR99 (rRNA gene classification)"
    echo ""
    echo "  Note: Tiara and Whokaryote models are bundled with their conda packages"
    echo "  (no separate database download needed)."
    echo ""
    echo "Default download directory: ${DB_DIR}"
    echo ""
    echo "After downloading, pass database paths to the pipeline:"
    echo "  --genomad_db ${DB_DIR}/genomad_db"
    echo "  --checkv_db  ${DB_DIR}/checkv_db"
    echo "  --checkm2_db ${DB_DIR}/checkm2_db"
    echo "  --kaiju_db   ${DB_DIR}/kaiju_db"
    echo "  --macsyfinder_models ${DB_DIR}/macsyfinder_models"
    echo "  --defensefinder_models ${DB_DIR}/defensefinder_models"
    echo "  --bakta_db   ${DB_DIR}/bakta_db"
    echo "  --kofam_db   ${DB_DIR}/kofam_db"
    echo "  --eggnog_db  ${DB_DIR}/eggnog_db"
    echo "  --dbcan_db   ${DB_DIR}/dbcan_db"
    echo "  --metaeuk_db ${DB_DIR}/metaeuk_db/metaeuk_db"
    echo "  --kraken2_db ${DB_DIR}/kraken2_db"
    echo "  --silva_ssu_db ${DB_DIR}/silva_db/SILVA_138.2_SSURef_NR99.fasta"
    echo "  --silva_lsu_db ${DB_DIR}/silva_db/SILVA_138.2_LSURef_NR99.fasta"
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
    echo "  5) macsyfinder - MacSyFinder models: TXSScan + CONJScan (~50 MB)"
    echo "  6) defensefinder - DefenseFinder models: anti-phage defense (~100 MB)"
    echo "  7) bakta     - Bakta annotation database (~37 GB)"
    echo "  8) kofam     - KOfam profiles (KEGG Orthology, ~4 GB)"
    echo "  9) eggnog    - eggNOG-mapper DIAMOND db (~12 GB)"
    echo " 10) dbcan     - dbCAN HMM + DIAMOND db (CAZyme, ~2 GB)"
    echo " 11) metaeuk   - MetaEuk OrthoDB v12 Eukaryota (~23 GB)"
    echo " 12) kraken2   - Kraken2 PlusPFP-8 (contig taxonomy, ~8 GB)"
    echo " 13) silva     - SILVA SSU + LSU NR99 (rRNA classification, ~900 MB)"
    echo " 14) all       - All databases"
    echo ""
    read -rp "Choice [1-14, or names]: " choice

    case "$choice" in
        1|genomad)  DOWNLOAD_GENOMAD=true ;;
        2|checkv)   DOWNLOAD_CHECKV=true ;;
        3|checkm2)  DOWNLOAD_CHECKM2=true ;;
        4|kaiju)    DOWNLOAD_KAIJU=true ;;
        5|macsyfinder) DOWNLOAD_MACSYFINDER=true ;;
        6|defensefinder) DOWNLOAD_DEFENSEFINDER=true ;;
        7|bakta)    DOWNLOAD_BAKTA=true ;;
        8|kofam)    DOWNLOAD_KOFAM=true ;;
        9|eggnog)   DOWNLOAD_EGGNOG=true ;;
        10|dbcan)   DOWNLOAD_DBCAN=true ;;
        11|metaeuk) DOWNLOAD_METAEUK=true ;;
        12|kraken2) DOWNLOAD_KRAKEN2=true ;;
        13|silva)   DOWNLOAD_SILVA=true ;;
        14|all)     DOWNLOAD_GENOMAD=true; DOWNLOAD_CHECKV=true; DOWNLOAD_CHECKM2=true; DOWNLOAD_KAIJU=true; DOWNLOAD_MACSYFINDER=true; DOWNLOAD_DEFENSEFINDER=true; DOWNLOAD_BAKTA=true; DOWNLOAD_KOFAM=true; DOWNLOAD_EGGNOG=true; DOWNLOAD_DBCAN=true; DOWNLOAD_METAEUK=true; DOWNLOAD_KRAKEN2=true; DOWNLOAD_SILVA=true ;;
        *)
            # Parse space-separated names
            for item in $choice; do
                case "$item" in
                    1|genomad)  DOWNLOAD_GENOMAD=true ;;
                    2|checkv)   DOWNLOAD_CHECKV=true ;;
                    3|checkm2)  DOWNLOAD_CHECKM2=true ;;
                    4|kaiju)    DOWNLOAD_KAIJU=true ;;
                    5|macsyfinder) DOWNLOAD_MACSYFINDER=true ;;
                    6|defensefinder) DOWNLOAD_DEFENSEFINDER=true ;;
                    7|bakta)    DOWNLOAD_BAKTA=true ;;
                    8|kofam)    DOWNLOAD_KOFAM=true ;;
                    9|eggnog)   DOWNLOAD_EGGNOG=true ;;
                    10|dbcan)   DOWNLOAD_DBCAN=true ;;
                    11|metaeuk) DOWNLOAD_METAEUK=true ;;
                    12|kraken2) DOWNLOAD_KRAKEN2=true ;;
                    13|silva)   DOWNLOAD_SILVA=true ;;
                    all)        DOWNLOAD_GENOMAD=true; DOWNLOAD_CHECKV=true; DOWNLOAD_CHECKM2=true; DOWNLOAD_KAIJU=true; DOWNLOAD_MACSYFINDER=true; DOWNLOAD_DEFENSEFINDER=true; DOWNLOAD_BAKTA=true; DOWNLOAD_KOFAM=true; DOWNLOAD_EGGNOG=true; DOWNLOAD_DBCAN=true; DOWNLOAD_METAEUK=true; DOWNLOAD_KRAKEN2=true; DOWNLOAD_SILVA=true ;;
                    *) echo "[WARNING] Unknown selection: $item" >&2 ;;
                esac
            done
            ;;
    esac

    if ! $DOWNLOAD_GENOMAD && ! $DOWNLOAD_CHECKV && ! $DOWNLOAD_CHECKM2 && ! $DOWNLOAD_KAIJU && ! $DOWNLOAD_MACSYFINDER && ! $DOWNLOAD_DEFENSEFINDER && ! $DOWNLOAD_BAKTA && ! $DOWNLOAD_KOFAM && ! $DOWNLOAD_EGGNOG && ! $DOWNLOAD_DBCAN && ! $DOWNLOAD_METAEUK && ! $DOWNLOAD_KRAKEN2 && ! $DOWNLOAD_SILVA; then
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

    local msf_data_bin="${ENV_DIR}/dana-mag-macsyfinder/bin/msf_data"
    if [ ! -x "${msf_data_bin}" ]; then
        echo "[ERROR] MacSyFinder not installed. Run ./install.sh first." >&2
        return 1
    fi

    mkdir -p "${db_path}"
    "${msf_data_bin}" install --target "${db_path}" TXSScan
    "${msf_data_bin}" install --target "${db_path}" CONJScan
    echo "[SUCCESS] MacSyFinder models downloaded to ${db_path}"
    echo "  Use with: --macsyfinder_models ${db_path}"
}

download_defensefinder() {
    local db_path="${DB_DIR}/defensefinder_models"
    if [ -d "${db_path}" ] && ls "${db_path}"/*.hmm 1>/dev/null 2>&1; then
        echo "[INFO] DefenseFinder models already exist at ${db_path}"
        echo "  Delete ${db_path} and re-run to force re-download."
        return 0
    fi

    echo ""
    echo "[INFO] Downloading DefenseFinder models (~100 MB)..."
    echo "  Models: ~280 anti-phage defense system HMM profiles"
    echo "  Destination: ${db_path}"

    local df_bin="${ENV_DIR}/dana-mag-defensefinder/bin/defense-finder"
    if [ ! -x "${df_bin}" ]; then
        echo "[ERROR] DefenseFinder not installed. Run ./install.sh first." >&2
        return 1
    fi

    mkdir -p "${db_path}"
    "${df_bin}" update --models-dir "${db_path}"
    # Workaround: CasFinder 3.1.1 has model definition version 2.1 which is
    # incompatible with macsyfinder 2.1.4 bundled with defense-finder 2.0.1
    # (https://github.com/mdmparis/defense-finder/issues/91)
    # Downgrade to 3.1.0 which uses compatible model definition version 2.0
    local msf_data_bin="${ENV_DIR}/dana-mag-defensefinder/bin/macsydata"
    if [ -x "${msf_data_bin}" ]; then
        local cf_ver
        cf_ver=$(cat "${db_path}/CasFinder/metadata.yml" 2>/dev/null | grep -oP 'vers: \K.*' | head -1)
        if [ "${cf_ver}" = "3.1.1" ]; then
            echo "  Downgrading CasFinder 3.1.1 → 3.1.0 (model version compatibility fix)"
            rm -rf "${db_path}/CasFinder"
            "${msf_data_bin}" install --target "${db_path}" CasFinder==3.1.0
        fi
    fi
    echo "[SUCCESS] DefenseFinder models downloaded to ${db_path}"
    echo "  Use with: --defensefinder_models ${db_path}"
}

download_bakta() {
    local db_path="${DB_DIR}/bakta_db"
    if [ -d "${db_path}" ] && [ -f "${db_path}/version.json" ]; then
        echo "[INFO] Bakta database already exists at ${db_path}"
        echo "  Delete ${db_path} and re-run to force re-download."
        return 0
    fi

    echo ""
    echo "[INFO] Downloading Bakta database (~37 GB full, or ~1.4 GB light)..."
    echo "  Destination: ${db_path}"

    local bakta_env="${ENV_DIR}/dana-mag-bakta"
    if [ ! -x "${bakta_env}/bin/bakta_db" ]; then
        echo "[ERROR] Bakta not installed. Run ./install.sh first." >&2
        return 1
    fi

    mkdir -p "${db_path}"
    # bakta_db requires AMRFinderPlus on PATH, so run inside the full conda env
    ${CONDA_CMD} run -p "${bakta_env}" bakta_db download --output "${db_path}" --type full
    echo "[SUCCESS] Bakta database downloaded to ${db_path}"
    echo "  Use with: --bakta_db ${db_path}/db"
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

    local emapper_bin="${ENV_DIR}/dana-mag-emapper/bin/download_eggnog_data.py"
    if [ ! -x "${emapper_bin}" ]; then
        echo "[ERROR] eggNOG-mapper not installed. Run ./install.sh first." >&2
        return 1
    fi

    mkdir -p "${db_path}"
    # Workaround: eggnog-mapper <=2.1.13 has broken download URLs
    # (eggnogdb.embl.de returns 404; eggnog5.embl.de is the working host)
    # https://github.com/eggnogdb/eggnog-mapper/issues/571
    sed -i "s|eggnogdb.embl.de|eggnog5.embl.de|g" "${emapper_bin}"
    "${emapper_bin}" --data_dir "${db_path}" -y
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

    local run_dbcan="${ENV_DIR}/dana-mag-dbcan/bin/run_dbcan"
    local hmmpress="${ENV_DIR}/dana-mag-dbcan/bin/hmmpress"
    if [ ! -x "${run_dbcan}" ]; then
        echo "[ERROR] dbCAN not installed. Run ./install.sh first." >&2
        return 1
    fi

    mkdir -p "${db_path}"
    # Use AWS S3 mirror — bcb.unl.edu has persistent SSL cert issues
    # (https://github.com/bcb-unl/run_dbcan/issues/53)
    "${run_dbcan}" database --db_dir "${db_path}" --aws_s3

    # Press HMM databases for HMMER (required before first run_dbcan use)
    echo "[INFO] Pressing HMM databases with hmmpress..."
    for hmm in "${db_path}"/dbCAN.hmm "${db_path}"/dbCAN-sub.hmm "${db_path}"/STP.hmm "${db_path}"/TF.hmm; do
        if [ -f "${hmm}" ] && [ ! -f "${hmm}.h3i" ]; then
            echo "  Pressing $(basename "${hmm}")..."
            "${hmmpress}" "${hmm}"
        fi
    done

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

    local metaeuk_bin="${ENV_DIR}/dana-mag-metaeuk/bin/metaeuk"
    if [ ! -x "${metaeuk_bin}" ]; then
        echo "[ERROR] MetaEuk not installed. Run ./install.sh first." >&2
        return 1
    fi

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
    "${metaeuk_bin}" createdb "${fasta}" "${db_path}/metaeuk_db"

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

if $DOWNLOAD_MACSYFINDER; then
    download_macsyfinder || failed=$((failed + 1))
fi

if $DOWNLOAD_DEFENSEFINDER; then
    download_defensefinder || failed=$((failed + 1))
fi

if $DOWNLOAD_BAKTA; then
    download_bakta || failed=$((failed + 1))
fi

if $DOWNLOAD_KOFAM; then
    download_kofam || failed=$((failed + 1))
fi

if $DOWNLOAD_EGGNOG; then
    download_eggnog || failed=$((failed + 1))
fi

if $DOWNLOAD_DBCAN; then
    download_dbcan || failed=$((failed + 1))
fi

if $DOWNLOAD_METAEUK; then
    download_metaeuk || failed=$((failed + 1))
fi

if $DOWNLOAD_KRAKEN2; then
    download_kraken2 || failed=$((failed + 1))
fi

if $DOWNLOAD_SILVA; then
    download_silva || failed=$((failed + 1))
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
    $DOWNLOAD_MACSYFINDER && echo "  --macsyfinder_models ${DB_DIR}/macsyfinder_models"
    $DOWNLOAD_DEFENSEFINDER && echo "  --defensefinder_models ${DB_DIR}/defensefinder_models"
    $DOWNLOAD_BAKTA   && echo "  --bakta_db   ${DB_DIR}/bakta_db/db"
    $DOWNLOAD_KOFAM   && echo "  --kofam_db   ${DB_DIR}/kofam_db"
    $DOWNLOAD_EGGNOG  && echo "  --eggnog_db  ${DB_DIR}/eggnog_db"
    $DOWNLOAD_DBCAN   && echo "  --dbcan_db   ${DB_DIR}/dbcan_db"
    $DOWNLOAD_METAEUK && echo "  --metaeuk_db ${DB_DIR}/metaeuk_db/metaeuk_db"
    $DOWNLOAD_KRAKEN2 && echo "  --kraken2_db ${DB_DIR}/kraken2_db"
    $DOWNLOAD_SILVA   && echo "  --silva_ssu_db ${DB_DIR}/silva_db/SILVA_138.2_SSURef_NR99.fasta"
    $DOWNLOAD_SILVA   && echo "  --silva_lsu_db ${DB_DIR}/silva_db/SILVA_138.2_LSURef_NR99.fasta"
fi
