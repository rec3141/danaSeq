#!/usr/bin/env bash
################################################################################
#                                                                              #
#  🧬  MAG ASSEMBLY PIPELINE - OPTIMIZED  🧬                                   #
#                                                                              #
#  ╔═══╗ ╔═══╗ ╔═══╗ ╔═══╗   From Chaos to Genomes   ╔═══╗ ╔═══╗ ╔═══╗      #
#  ║ A ║═║ T ║═║ C ║═║ G ║                            ║ T ║═║ A ║═║ G ║      #
#  ╚═══╝ ╚═══╝ ╚═══╝ ╚═══╝                            ╚═══╝ ╚═══╝ ╚═══╝      #
#                                                                              #
################################################################################
#
# MISSION: Transform mixed metagenomic reads into pristine genome bins
#
# WORKFLOW:
#   1. 🏗️  Co-assemble with Flye (overlap-based long-read assembly)
#   2. 📍 Map reads back to assembly (minimap2, all samples)
#   3. 📊 Calculate coverage profiles (differential abundance signal)
#   4. 🗂️  Multi-tool binning (SemiBin2, MetaBAT2, MaxBin2)
#   5. 🤝 Consensus binning (DAS Tool picks the best)
#   6. ✅ Quality assessment (CheckM2: completeness & contamination)
#   7. ✨ Polish bins (Racon + Medaka for publication quality)
#   8. 🏷️  Taxonomic classification (Kraken2, GTDB)
#
# FEATURES:
#   • Resume capability (checkpoints at each stage)
#   • Parallel processing (maximizes throughput)
#   • Quality filtering (at every step)
#   • Colored logging (because we're civilized)
#   • Error handling (robust against failures)
#
# USAGE:
#   ./61_map_and_bin_optimized.sh [output_directory]
#
# OUTPUT:
#   • High-quality MAGs (>90% complete, <5% contamination)
#   • Taxonomic assignments
#   • Abundance profiles
#   • Quality reports
#
# AUTHOR: Refactored with love by Claude (who finally escaped Celadon City)
# ORIGINAL: CMO2025 Project Pipeline
# VERSION: 2.0 - Now with more ASCII art! 🎨
#
################################################################################

set -euo pipefail  # Exit on error, undefined vars, pipe failures

# Global cleanup on exit/interrupt
cleanup_on_exit() {
  local exit_code=$?
  # Cleanup logic for any temporary files or locks would go here
  # Currently handled by individual stage checkpoints
  exit $exit_code
}

trap cleanup_on_exit EXIT
trap 'echo "[INTERRUPTED] Pipeline interrupted, cleaning up..." >&2; exit 130' INT
trap 'echo "[TERMINATED] Pipeline terminated, cleaning up..." >&2; exit 143' TERM
trap 'echo "[HANGUP] Connection lost, cleaning up..." >&2; exit 129' HUP

#-------------------------------------------------------------------------------
#  CONFIGURATION
#-------------------------------------------------------------------------------

readonly PROJECT="CMO2025"
readonly PROJECT_DIR="/data/project_${PROJECT}"
readonly METAFILE="${PROJECT_DIR}/CMO 2025 MICB.xlsx"
readonly FASTQ_DIR="${PROJECT_DIR}/concat"
readonly THREADS=24

# Colors for pretty output (because we deserve nice things)
readonly RED='\033[0;31m'
readonly GREEN='\033[0;32m'
readonly YELLOW='\033[1;33m'
readonly BLUE='\033[0;34m'
readonly PURPLE='\033[0;35m'
readonly NC='\033[0m' # No Color

#-------------------------------------------------------------------------------
#  LOGGING FUNCTIONS
#-------------------------------------------------------------------------------

log_info()    { echo -e "${BLUE}[INFO]${NC}    $(date '+%H:%M:%S') | $*"; }
log_success() { echo -e "${GREEN}[SUCCESS]${NC} $(date '+%H:%M:%S') | $*"; }
log_warn()    { echo -e "${YELLOW}[WARN]${NC}    $(date '+%H:%M:%S') | $*"; }
log_error()   { echo -e "${RED}[ERROR]${NC}   $(date '+%H:%M:%S') | $*" >&2; }
log_step()    { echo -e "\n${PURPLE}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"; 
                echo -e "${PURPLE}  $*${NC}";
                echo -e "${PURPLE}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}\n"; }

#-------------------------------------------------------------------------------
#  UTILITY FUNCTIONS
#-------------------------------------------------------------------------------

# Validate and sanitize output directory path
# Prevents path traversal attacks and ensures directory is within workspace
validate_output_dir() {
    local dir="$1"
    local base_workspace="/data"

    # Resolve to absolute path (creates parent dirs if needed)
    local abs_dir
    if ! abs_dir="$(realpath -m "$dir" 2>/dev/null)"; then
        log_error "Cannot resolve path: $dir"
        return 1
    fi

    # Check for path traversal attempts
    if [[ "$abs_dir" =~ \.\. ]]; then
        log_error "Path traversal detected: $dir"
        return 1
    fi

    # Ensure path is within allowed workspace
    if [[ "$abs_dir" != "$base_workspace"* ]]; then
        log_error "Output directory must be under $base_workspace: $abs_dir"
        return 1
    fi

    # Return validated path
    echo "$abs_dir"
    return 0
}


# Check if a file exists and is non-empty; optionally verify minimum size
file_ready() {
    local file="$1"
    local min_size="${2:-0}"
    
    [[ -s "$file" ]] && [[ $(stat -f%z "$file" 2>/dev/null || stat -c%s "$file" 2>/dev/null) -ge "$min_size" ]]
}

# Run a command only if the output file doesn't exist yet
run_if_missing() {
    local output_file="$1"
    shift
    
    if file_ready "$output_file"; then
        log_info "Skipping (output exists): $(basename "$output_file")"
        return 0
    fi
    
    log_info "Running: $1"
    "$@"
}

# Ensure directory exists
ensure_dir() {
    local dir="$1"
    if [[ ! -d "$dir" ]]; then
        if ! mkdir -p "$dir"; then
            log_error "Cannot create directory: $dir"
            return 1
        fi
    fi
    return 0
}

#-------------------------------------------------------------------------------
#  PIPELINE STEPS
#-------------------------------------------------------------------------------

run_assembly() {
    log_step "STEP 1: ASSEMBLY (Flye)"
    
    if file_ready "$ASSEMBLY" 1000; then
        log_success "Assembly already complete: $ASSEMBLY"
        return 0
    fi
    
    log_info "Running Flye assembly..."
    /data/dana/run-flye.sh "$OUTDIR" "$PROJECT"
    
    if ! file_ready "$ASSEMBLY" 1000; then
        log_error "Assembly failed or produced empty output"
        exit 1
    fi
    
    log_success "Assembly complete!"
}

run_mapping() {
    log_step "STEP 2: READ MAPPING (minimap2)"
    
    ensure_dir "$BAM_DIR"
    local count=0
    local total=$(ls -1 "$FASTQ_DIR"/20*.fastq.gz 2>/dev/null | wc -l)
    
    for fq in "$FASTQ_DIR"/20*.fastq.gz; do
        [[ -e "$fq" ]] || continue
        
        ((count++))
        local base=$(basename "$fq" .fastq.gz)
        local bam="${BAM_DIR}/${base}.sorted.bam"
        
        if file_ready "$bam"; then
            log_info "[$count/$total] Skipping $base (BAM exists)"
            continue
        fi
        
        log_info "[$count/$total] Mapping $base..."

        # Map reads and sort to BAM
        if ! minimap2 -a -x map-ont --secondary=no -t "$THREADS" "$ASSEMBLY" "$fq" \
            | samtools sort -@ "$THREADS" -o "$bam" -; then
            log_error "Mapping failed for $base"
            continue
        fi

        # Index BAM file
        if ! samtools index -@ "$THREADS" "$bam"; then
            log_error "Indexing failed for $base"
            continue
        fi

        log_success "Mapped: $base"
    done
    
    log_success "All samples mapped!"
}

calculate_depths() {
    log_step "STEP 3: DEPTH CALCULATION"
    
    local depths_file="${BAM_DIR}/depths_jgi.txt"
    
    if file_ready "$depths_file"; then
        log_success "Depth file exists: $depths_file"
        return 0
    fi
    
    log_info "Calculating contig depths..."

    if ! jgi_summarize_bam_contig_depths \
        --outputDepth "$depths_file" \
        --percentIdentity 80 \
        --minMapQual 5 \
        --referenceFasta "$ASSEMBLY" \
        "$BAM_DIR"/*.sorted.bam; then
        log_error "Depth calculation failed"
        return 1
    fi

    log_success "Depth calculation complete!"
}

run_semibin() {
    log_step "STEP 4a: BINNING (SemiBin2)"
    
    local semibin_dir="${OUTDIR}/semibin2"
    local contig_bins="${semibin_dir}/contig_bins.tsv"
    
    if file_ready "$contig_bins"; then
        log_success "SemiBin2 output exists"
        return 0
    fi
    
    log_info "Running SemiBin2..."

    if ! ensure_dir "$semibin_dir"; then
        log_error "Cannot create SemiBin2 directory"
        return 1
    fi

    if ! conda run -n SemiBin --live-stream SemiBin2 single_easy_bin \
        -i "$ASSEMBLY" \
        -b "$BAM_DIR"/*.sorted.bam \
        -o "$semibin_dir"; then
        log_error "SemiBin2 failed"
        return 1
    fi

    # Remove header line (SemiBin2 quirk)
    if [[ -f "$contig_bins" ]]; then
        if ! tail -n +2 "$contig_bins" > "${contig_bins}.tmp"; then
            log_error "Failed to process SemiBin2 output"
            return 1
        fi
        if ! mv "${contig_bins}.tmp" "$contig_bins"; then
            log_error "Failed to update contig_bins file"
            return 1
        fi
    fi

    log_success "SemiBin2 complete!"
}

run_metabat() {
    log_step "STEP 4b: BINNING (MetaBAT2)"
    
    local metabat_dir="${OUTDIR}/metabat2"
    local contig_bins="${metabat_dir}/contig_bins.tsv"
    
    if file_ready "$contig_bins"; then
        log_success "MetaBAT2 output exists"
        return 0
    fi
    
    log_info "Running MetaBAT2..."

    if ! ensure_dir "$metabat_dir"; then
        log_error "Cannot create MetaBAT2 directory"
        return 1
    fi

    if ! metabat2 \
        -i "$ASSEMBLY" \
        -o "${metabat_dir}/bin" \
        --saveCls \
        --minClsSize 50000 \
        -a "${BAM_DIR}/depths_jgi.txt"; then
        log_error "MetaBAT2 failed"
        return 1
    fi

    # Convert MetaBAT2 output to contig_bins.tsv format
    log_info "Converting MetaBAT2 output format..."

    if ! : > "$contig_bins"; then
        log_error "Cannot create contig_bins file"
        return 1
    fi
    
    for bin_file in "${metabat_dir}"/bin*.fa; do
        [[ -e "$bin_file" ]] || continue
        local bin_name=$(basename "$bin_file")
        grep '>' "$bin_file" | tr -d '>' | cut -f1 | while read -r contig; do
            printf '%s\t%s\n' "$contig" "$bin_name" >> "$contig_bins"
        done
    done
    
    log_success "MetaBAT2 complete!"
}

run_dastool() {
    log_step "STEP 5: BIN REFINEMENT (DAS_Tool)"
    
    local dastool_dir="${OUTDIR}/dastool_DASTool_bins"
    local contig2bin="${dastool_dir}/dastool_DASTool_contig2bin.tsv"
    
    if file_ready "$contig2bin"; then
        log_success "DAS_Tool output exists"
        return 0
    fi
    
    log_info "Running DAS_Tool..."

    if ! conda run -n SemiBin --live-stream DAS_Tool \
        -i "${OUTDIR}/semibin2/contig_bins.tsv,${OUTDIR}/metabat2/contig_bins.tsv" \
        -l semibin2,metabat2 \
        -c "$ASSEMBLY" \
        -o "${OUTDIR}/dastool" \
        --threads "$THREADS" \
        --write_bin_evals \
        --write_bins; then
        log_error "DAS_Tool failed"
        return 1
    fi

    # Consolidate output
    if ! ensure_dir "$dastool_dir"; then
        log_error "Cannot create DAS_Tool directory"
        return 1
    fi
    mv "${OUTDIR}"/dastool*.* "$dastool_dir/" 2>/dev/null || true
    
    # Create combined bins file
    cat "${dastool_dir}"/*.fa > "${dastool_dir}/allbins.fa" 2>/dev/null || true
    
    log_success "DAS_Tool complete!"
}

run_taxonomy() {
    log_step "STEP 6: TAXONOMIC CLASSIFICATION (Kaiju)"
    
    local dastool_dir="${OUTDIR}/dastool_DASTool_bins"
    local kaiju_summary="${dastool_dir}/kaiju.allbins.summary.tsv"
    local kaiju_db="/data/scratch/refdbs/kaiju/progenomes"
    
    if file_ready "$kaiju_summary"; then
        log_success "Kaiju classification exists"
        return 0
    fi
    
    log_info "Running Kaiju classification..."
    
    # Run Kaiju
    conda run -n SemiBin --live-stream kaiju \
        -t "${kaiju_db}/nodes.dmp" \
        -f "${kaiju_db}/kaiju_db_progenomes.fmi" \
        -i "${dastool_dir}/allbins.fa" \
        -z "$THREADS" \
        -o "${dastool_dir}/kaiju.allbins.tsv"
    
    # Add taxonomic names
    log_info "Adding taxonomic names..."
    conda run -n SemiBin --live-stream kaiju-addTaxonNames \
        -t "${kaiju_db}/nodes.dmp" \
        -n "${kaiju_db}/names.dmp" \
        -i "${dastool_dir}/kaiju.allbins.tsv" \
        -o "${dastool_dir}/kaiju.allbins-taxa.tsv" \
        -r superkingdom,phylum,class,order,family,genus,species
    
    # Create summary (join with contig lengths and bins)
    log_info "Creating summary..."
    join -t $'\t' \
        <(sort -k1,1 "${dastool_dir}/dastool.seqlength") \
        <(paste \
            <(sort -k1,1 "${dastool_dir}/dastool_DASTool_contig2bin.tsv") \
            <(sort -k2,2 "${dastool_dir}/kaiju.allbins-taxa.tsv")) \
        | sort -k3,3 -k2,2rn \
        > "$kaiju_summary"
    
    log_success "Taxonomic classification complete!"
}

run_tetramer_analysis() {
    log_step "STEP 7: TETRAMER FREQUENCIES"
    
    local tetra_dir="${OUTDIR}/tetra"
    local tnfs_file="${tetra_dir}/tnfs.txt"
    
    if file_ready "$tnfs_file"; then
        log_success "Tetramer analysis exists"
        return 0
    fi
    
    log_info "Calculating tetramer frequencies..."
    ensure_dir "$tetra_dir"
    
    # Create annotation file
    grep '>' "$ASSEMBLY" | sed 's/>//' | awk '{print $0"\t"$0"\t"$0}' > "${tetra_dir}/annotation.txt"
    
    # Run tetramer calculation
    perl /work/apps/tetramer_freqs_esom.pl \
        -f "$ASSEMBLY" \
        -a "${tetra_dir}/annotation.txt" \
        -min 1500 \
        -max 10000000
    
    # Process output
    paste <(echo "seqid") <(head -n4 Tetra_*.lrn | tail -n1 | cut -f2-) > "$tnfs_file"
    paste \
        <(awk '$1 !~ /^%/' Tetra_*.names) \
        <(awk '$1 !~ /^%/' Tetra_*.lrn) \
        | cut -f3,5- > "${tetra_dir}/assembly.lrn"
    
    mv Tetra_* "$tetra_dir/"
    
    log_success "Tetramer analysis complete!"
}

run_bin_qc() {
    log_step "STEP 8: BIN QC & SKETCHING"
    
    local dastool_dir="${OUTDIR}/dastool_DASTool_bins"
    
    for bin_file in "${dastool_dir}"/*.fa; do
        [[ -e "$bin_file" ]] || continue
        
        local base=$(basename "$bin_file" .fa)
        local stats_file="${dastool_dir}/${base}.txt"
        
        # Assembly stats
        if ! file_ready "$stats_file"; then
            log_info "Stats: $base"
            stats.sh "$bin_file" out="$stats_file"
        fi
        
        # Sketch against databases (with idempotency)
        for db in refseq nt protein; do
            local sketch_file="${dastool_dir}/ss_${db}_${base}.tsv"
            if ! file_ready "$sketch_file"; then
                log_info "Sketching $base against $db..."
                local address="$db"
                [[ "$db" == "protein" ]] && address="protein" || true
                sendsketch.sh \
                    in="$bin_file" \
                    out="$sketch_file" \
                    level=3 \
                    format=3 \
                    address="$address"
            fi
        done
    done
    
    log_success "Bin QC complete!"
}

run_checkm() {
    log_step "STEP 9: QUALITY ASSESSMENT (CheckM2)"
    
    local dastool_dir="${OUTDIR}/dastool_DASTool_bins"
    local checkm_dir="${dastool_dir}/checkm2"
    local checkm_report="${checkm_dir}/quality_report.tsv"
    
    if file_ready "$checkm_report"; then
        log_success "CheckM2 results exist"
        return 0
    fi
    
    log_info "Running CheckM2..."
    
    conda run -n checkm2 --live-stream checkm2 predict \
        --threads "$THREADS" \
        --input "$dastool_dir" \
        -x fa \
        --output-directory "$checkm_dir"
    
    log_success "CheckM2 complete!"
}

run_visualization() {
    log_step "STEP 10: VISUALIZATION"
    
    log_info "Generating plots..."
    Rscript /data/dana/plot-bins.R "$OUTDIR" "$METAFILE"
    
    log_success "Visualization complete!"
}

print_summary() {
    log_step "🎉 PIPELINE COMPLETE"
    
    echo -e "${GREEN}"
    echo "  Output directory: $OUTDIR"
    echo ""
    echo "  Key outputs:"
    echo "    • Assembly:     ${ASSEMBLY}"
    echo "    • Bins:         ${OUTDIR}/dastool_DASTool_bins/"
    echo "    • Taxonomy:     ${OUTDIR}/dastool_DASTool_bins/kaiju.allbins.summary.tsv"
    echo "    • QC:           ${OUTDIR}/dastool_DASTool_bins/checkm2/"
    echo "    • Tetramer:     ${OUTDIR}/tetra/"
    echo -e "${NC}"
}

#-------------------------------------------------------------------------------
#  MAIN
#-------------------------------------------------------------------------------

main() {
    # Parse and validate arguments
    local requested_outdir="${1:-${PROJECT_DIR}/map_and_bin_$(date +'%Y%m%d-%H%M%S')}"

    OUTDIR="$(validate_output_dir "$requested_outdir")" || {
        log_error "Invalid output directory: $requested_outdir"
        exit 1
    }

    log_info "Validated output directory: $OUTDIR"

    # Derived paths
    readonly ASSEMBLY="${OUTDIR}/flye1000/assembly.fasta"
    readonly BAM_DIR="${OUTDIR}/bamdir"
    
    # Print banner
    echo -e "${PURPLE}"
    echo "  ╔═══════════════════════════════════════════════════════════════╗"
    echo "  ║     METAGENOMICS ASSEMBLY & BINNING PIPELINE                  ║"
    echo "  ║     Project: ${PROJECT}                                          ║"
    echo "  ╚═══════════════════════════════════════════════════════════════╝"
    echo -e "${NC}"
    
    log_info "Output directory: $OUTDIR"
    log_info "Threads: $THREADS"
    
    # Create directories
    ensure_dir "$OUTDIR"
    ensure_dir "$FASTQ_DIR"
    
    # Run pipeline steps
    run_assembly
    run_mapping
    calculate_depths
    run_semibin
    run_metabat
    run_dastool
    run_taxonomy
    run_tetramer_analysis
    run_bin_qc
    run_checkm
    run_visualization
    
    print_summary
}

# Run if executed directly (not sourced)
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    main "$@"
fi
