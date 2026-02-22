#!/usr/bin/env bash
# repair-bins-tsv.sh — Regenerate *_bins.tsv from FASTA headers
#
# Fixes the bin numbering mismatch where the TSV was built from the binner's
# internal numbering but the FASTAs were renamed via a lexicographic glob,
# causing the two to diverge for bins >= 10.
#
# The FASTA files (bins/*.fa) are the ground truth — each file is named
# e.g. semibin_319.fa and contains the actual contigs. This script rebuilds
# the TSV by extracting contig names from the FASTA headers.
#
# Usage:
#   ./repair-bins-tsv.sh <results_dir> [--dry-run]
#
# Also repairs storeDir copies if --store-dir is given:
#   ./repair-bins-tsv.sh <results_dir> --store-dir /data/scratch/mag_store/results_Ebb_Flow_20260218

set -euo pipefail

RESULTS_DIR="${1:?Usage: $0 <results_dir> [--store-dir <path>] [--dry-run]}"
shift

DRY_RUN=false
STORE_DIR=""
while [[ $# -gt 0 ]]; do
    case "$1" in
        --dry-run) DRY_RUN=true; shift ;;
        --store-dir) STORE_DIR="$2"; shift 2 ;;
        *) echo "[ERROR] Unknown argument: $1" >&2; exit 1 ;;
    esac
done

BINNERS=(semibin metabat maxbin lorbin comebin)

repair_binner() {
    local dir="$1"
    local binner="$2"
    local binner_dir="$dir/binning/$binner"
    local tsv="$binner_dir/${binner}_bins.tsv"
    local bins_dir="$binner_dir/bins"

    if [ ! -d "$bins_dir" ]; then
        echo "  [SKIP] $binner — no bins/ directory in $dir"
        return
    fi

    local fa_files=()
    for f in "$bins_dir"/*.fa; do
        [ -e "$f" ] && fa_files+=("$f")
    done
    local n_fastas=${#fa_files[@]}
    if [ "$n_fastas" -eq 0 ]; then
        echo "  [SKIP] $binner — no .fa files in $bins_dir"
        return
    fi

    # Count contigs in existing TSV vs FASTAs
    local old_count=0
    if [ -f "$tsv" ]; then
        old_count=$(wc -l < "$tsv")
    fi

    local new_count=0
    for fa in "$bins_dir"/*.fa; do
        local c
        c=$(grep -c '>' "$fa" 2>/dev/null || echo 0)
        new_count=$((new_count + c))
    done

    if $DRY_RUN; then
        echo "  [DRY-RUN] $binner: $n_fastas FASTAs, TSV has $old_count lines, FASTAs have $new_count contigs"
        return
    fi

    # Backup existing TSV
    if [ -f "$tsv" ]; then
        cp "$tsv" "${tsv}.bak"
    fi

    # Rebuild TSV from FASTA headers (awk for speed)
    > "$tsv"
    for fa in "$bins_dir"/*.fa; do
        local bin_name
        bin_name=$(basename "$fa" .fa)
        awk -v bin="$bin_name" '/^>/{sub(/^>/,""); sub(/ .*/,""); print $0 "\t" bin}' "$fa" >> "$tsv"
    done

    local new_lines
    new_lines=$(wc -l < "$tsv")
    echo "  [FIXED] $binner: $old_count → $new_lines contig assignments ($n_fastas bins)"
}

echo "Repairing bins TSV files from FASTA headers..."
echo "Results: $RESULTS_DIR"
$DRY_RUN && echo "Mode: DRY RUN"

echo ""
echo "=== Results directory ==="
for binner in "${BINNERS[@]}"; do
    repair_binner "$RESULTS_DIR" "$binner"
done

if [ -n "$STORE_DIR" ]; then
    echo ""
    echo "=== Store directory ==="
    echo "Store: $STORE_DIR"
    for binner in "${BINNERS[@]}"; do
        repair_binner "$STORE_DIR" "$binner"
    done
fi

echo ""
echo "Done. Original TSVs backed up as *.bak"
