#!/usr/bin/env bash
set -euo pipefail
testing=${1:-""} #use --testing for testing

# --- paths (edit these two if needed) ---
INPUT_DIR="/mnt/Amundsen/Data/FULL_CSV/2025_LEG_04/"   # where ACSD_YYYYMMDD.csv live
WEBROOT="/home/cryomics/Desktop/Amundsen-Collins/underway_dashboard/"                                     # your bind-mounted webroot
# this gets mounted to /var/www/html/underway
OUTPUT_DIR=$(mktemp -d)
STATION_CSV="/mnt/Amundsen/Data/Rosette/2025_LEG_04/Logs/2025_04_CTD_logbook.csv"
RSCRIPT="/usr/bin/Rscript"
DASH_R="/home/cryomics/Desktop/underway_dev/underway_dashboard_v6.R"

# persistent cache lives beside output
CACHE_DIR="/home/cryomics/Desktop/underway_dev/cache"
mkdir -p "$CACHE_DIR" "$OUTPUT_DIR"


# avoid overlapping runs (unit is also set to be non-parallel, this is belt & suspenders)
exec 9>"$CACHE_DIR/.run.lock"
flock -n 9 || { echo "Another run is in progress; exiting."; exit 0; }

# Be nice to the box + ensure files end up world-readable
umask 002
export TZ=America/Toronto
export GCAL_SERVICE_JSON="/home/cryomics/underway_dev/.noble-vortex-471516-d6-7c784035a438.json"

# Generate dashboard (rolling windows + surprise panel)
CMD="nice -n 19 ionice -c3 $RSCRIPT $DASH_R --indir $INPUT_DIR --outdir $OUTPUT_DIR $testing"

eval $CMD
rsync -av ${OUTPUT_DIR}/ ${WEBROOT}

echo "[OK] $(date -u) report generated"
