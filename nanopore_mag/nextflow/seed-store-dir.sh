#!/usr/bin/env bash
set -euo pipefail

# ============================================================================
# Seed a storeDir from existing pipeline results
# ============================================================================
#
# Creates a storeDir structure from an existing results directory.
# Once seeded, the pipeline can use --store_dir to skip all completed processes.
#
# Usage:
#   ./seed-store-dir.sh [--mode hardlink|symlink|copy] <results_dir> <store_dir>
#
# Modes:
#   hardlink  - Hard links (no extra disk space; same filesystem required)
#   symlink   - Symbolic links (cross-filesystem; breaks if results move)
#   copy      - Full copies (cross-filesystem; uses extra disk space)
#   auto      - Default: hardlink if same filesystem, else symlink
#
# Example:
#   ./seed-store-dir.sh full_test_20260216 /data/scratch/mag_store
#   ./seed-store-dir.sh --mode symlink full_test_20260216 /data/scratch/mag_store
#
# Notes:
#   - Handles legacy contig_bins.tsv → {name}_bins.tsv rename
#   - Handles legacy PROKKA_* → annotation.* rename
#   - Skips pipeline_info/ (not a process output)
#   - Safe to re-run (skips existing files)
#
# ============================================================================

die() { echo "[ERROR] $1" >&2; exit 1; }

# Parse --mode flag
MODE="auto"
if [[ "${1:-}" == "--mode" ]]; then
    MODE="${2:-}"
    shift 2
    [[ "$MODE" =~ ^(hardlink|symlink|copy|auto)$ ]] || die "Invalid mode: $MODE (use hardlink, symlink, copy, or auto)"
fi

[[ $# -eq 2 ]] || die "Usage: $0 [--mode hardlink|symlink|copy|auto] <results_dir> <store_dir>"

RESULTS="$(realpath "$1")"
STORE="$(realpath -m "$2")"

[[ -d "$RESULTS" ]] || die "Results directory not found: $1"

# Auto-detect mode: hardlink if same filesystem, else symlink
if [[ "$MODE" == "auto" ]]; then
    mkdir -p "$STORE"
    src_dev=$(stat -c %d "$RESULTS")
    dst_dev=$(stat -c %d "$STORE")
    if [[ "$src_dev" == "$dst_dev" ]]; then
        MODE="hardlink"
    else
        MODE="symlink"
    fi
    echo "[INFO] Auto-detected mode: $MODE"
fi

echo "[INFO] Results: $RESULTS"
echo "[INFO] Store:   $STORE"
echo "[INFO] Mode:    $MODE"

# Track stats
linked=0
skipped=0
renamed=0

# ----------------------------------------------------------------------------
# Helper: link/copy a single file (create parent dirs as needed)
# ----------------------------------------------------------------------------
link_file() {
    local src="$1" dst="$2"
    if [[ -e "$dst" ]]; then
        ((skipped++)) || true
        return
    fi
    mkdir -p "$(dirname "$dst")"
    case "$MODE" in
        hardlink) ln "$src" "$dst" ;;
        symlink)  ln -s "$src" "$dst" ;;
        copy)     cp "$src" "$dst" ;;
    esac
    ((linked++)) || true
}

# ----------------------------------------------------------------------------
# Helper: link/copy all files in a directory recursively
# ----------------------------------------------------------------------------
link_dir() {
    local src_dir="$1" dst_dir="$2"
    [[ -d "$src_dir" ]] || return 0
    while IFS= read -r -d '' f; do
        local rel="${f#"$src_dir"/}"
        link_file "$f" "$dst_dir/$rel"
    done < <(find "$src_dir" -type f -print0)
}

# ----------------------------------------------------------------------------
# 1. Direct mappings (filenames match exactly)
# ----------------------------------------------------------------------------

# Assembly
link_dir "$RESULTS/assembly" "$STORE/assembly"

# Mapping (per-sample BAMs + depths)
link_dir "$RESULTS/mapping" "$STORE/mapping"

# Concat (per-sample fastq.gz)
link_dir "$RESULTS/concat" "$STORE/concat"

# DAS Tool consensus
link_dir "$RESULTS/binning/dastool" "$STORE/binning/dastool"

# CheckM2
link_dir "$RESULTS/binning/checkm2" "$STORE/binning/checkm2"

# NCLB
link_dir "$RESULTS/binning/nclb" "$STORE/binning/nclb"

# Annotation (Bakta already uses annotation.* names)
link_dir "$RESULTS/annotation/bakta/basic" "$STORE/annotation/bakta/basic"
link_dir "$RESULTS/annotation/bakta/extra" "$STORE/annotation/bakta/extra"

# Annotation (Prokka — may have PROKKA_* or annotation.* names)
if [[ -d "$RESULTS/annotation/prokka" ]]; then
    for src in "$RESULTS/annotation/prokka"/*; do
        [[ -f "$src" ]] || continue
        base="$(basename "$src")"
        # Rename PROKKA_<timestamp>.ext → annotation.ext
        if [[ "$base" == PROKKA_* ]]; then
            ext="${base#*.}"  # everything after first dot
            dst="$STORE/annotation/prokka/annotation.$ext"
            link_file "$src" "$dst"
            ((renamed++)) || true
        else
            link_file "$src" "$STORE/annotation/prokka/$base"
        fi
    done
fi

# Taxonomy
link_dir "$RESULTS/taxonomy/kaiju" "$STORE/taxonomy/kaiju"
link_dir "$RESULTS/taxonomy/kraken2" "$STORE/taxonomy/kraken2"
link_dir "$RESULTS/taxonomy/sendsketch" "$STORE/taxonomy/sendsketch"
link_dir "$RESULTS/taxonomy/rrna" "$STORE/taxonomy/rrna"

# Eukaryotic
link_dir "$RESULTS/eukaryotic/tiara" "$STORE/eukaryotic/tiara"
link_dir "$RESULTS/eukaryotic/whokaryote" "$STORE/eukaryotic/whokaryote"
link_dir "$RESULTS/eukaryotic/metaeuk" "$STORE/eukaryotic/metaeuk"
link_dir "$RESULTS/eukaryotic/marferret" "$STORE/eukaryotic/marferret"

# MGE
link_dir "$RESULTS/mge/genomad" "$STORE/mge/genomad"
link_dir "$RESULTS/mge/checkv" "$STORE/mge/checkv"
link_dir "$RESULTS/mge/integrons" "$STORE/mge/integrons"
link_dir "$RESULTS/mge/islandpath" "$STORE/mge/islandpath"
link_dir "$RESULTS/mge/genomic_islands" "$STORE/mge/islandpath"  # legacy name
link_dir "$RESULTS/mge/macsyfinder" "$STORE/mge/macsyfinder"
link_dir "$RESULTS/mge/defensefinder" "$STORE/mge/defensefinder"

# Metabolism (most are direct)
link_dir "$RESULTS/metabolism/kofamscan" "$STORE/metabolism/kofamscan"
link_dir "$RESULTS/metabolism/emapper" "$STORE/metabolism/emapper"
link_dir "$RESULTS/metabolism/dbcan" "$STORE/metabolism/dbcan"
link_dir "$RESULTS/metabolism/merged" "$STORE/metabolism/merged"
link_dir "$RESULTS/metabolism/modules" "$STORE/metabolism/modules"
link_dir "$RESULTS/metabolism/minpath" "$STORE/metabolism/minpath"
link_dir "$RESULTS/metabolism/kegg_decoder" "$STORE/metabolism/kegg_decoder"

# Metabolism: MAP_TO_BINS — storeDir puts per_mag/ + community_annotations.tsv together
link_dir "$RESULTS/metabolism/per_mag" "$STORE/metabolism/per_mag/per_mag"
if [[ -f "$RESULTS/metabolism/community/community_annotations.tsv" ]]; then
    link_file "$RESULTS/metabolism/community/community_annotations.tsv" \
                  "$STORE/metabolism/per_mag/community_annotations.tsv"
fi

# ----------------------------------------------------------------------------
# 2. Binner TSV renames (contig_bins.tsv → {name}_bins.tsv)
# ----------------------------------------------------------------------------

for binner in semibin metabat maxbin lorbin comebin; do
    binner_dir="$RESULTS/binning/$binner"
    [[ -d "$binner_dir" ]] || continue

    # Link bins/ directory as-is
    link_dir "$binner_dir/bins" "$STORE/binning/$binner/bins"

    # Rename contig_bins.tsv → {binner}_bins.tsv
    if [[ -f "$binner_dir/contig_bins.tsv" ]]; then
        link_file "$binner_dir/contig_bins.tsv" "$STORE/binning/$binner/${binner}_bins.tsv"
        ((renamed++)) || true
    elif [[ -f "$binner_dir/${binner}_bins.tsv" ]]; then
        # Already in new naming convention
        link_file "$binner_dir/${binner}_bins.tsv" "$STORE/binning/$binner/${binner}_bins.tsv"
    fi
done

# ----------------------------------------------------------------------------
# Summary
# ----------------------------------------------------------------------------

echo ""
echo "[SUCCESS] Store directory seeded"
echo "  Linked:   $linked files ($MODE)"
echo "  Skipped:  $skipped files (already existed)"
echo "  Renamed:  $renamed files (legacy → current naming)"
echo ""
echo "Usage:"
echo "  ./run-mag.sh --input <reads> --outdir <output> --store_dir $STORE"
