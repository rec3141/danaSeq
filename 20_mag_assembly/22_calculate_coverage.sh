#!/usr/bin/env bash
set -euo pipefail

# ===========================================================
# ONT MAG post-polish helper
# 1) Build contig-depth tables for binning (CoverM; optional jgi)
# 2) (Optional) Run CheckM2 and GUNC on bins if provided
# ===========================================================

source activate checkm2

# Default values
THREADS=32
ASSEMBLY=""
READS=""
OUTDIR="coverage-out"
COVERM_IDENTITY=0.85
COVERM_ALIGNED_FRAC=0.50
COVERM_METHODS="mean covered_bases variance tpm rpkm"
BAMS=""
RUN_JGI=0
JGI_MIN_MAPQ=5
JGI_MIN_PID=80
BINS_DIR=""
BINS_EXT="fa"
CHECKM2_DB=""
GUNC_DB=""

# Check if terminal supports colors
if [ -t 1 ] && [ -n "${TERM:-}" ] && [ "${TERM:-}" != "dumb" ] && command -v tput >/dev/null 2>&1; then
    ncolors=$(tput colors 2>/dev/null || echo 0)
    if [ "$ncolors" -ge 8 ]; then
        # Use tput for better compatibility
        RED=$(tput setaf 1)
        GREEN=$(tput setaf 2)
        YELLOW=$(tput setaf 3)
        BLUE=$(tput setaf 4)
        MAGENTA=$(tput setaf 5)
        CYAN=$(tput setaf 6)
        BOLD=$(tput bold)
        NC=$(tput sgr0)  # Reset
    else
        # No color support
        RED=""
        GREEN=""
        YELLOW=""
        BLUE=""
        MAGENTA=""
        CYAN=""
        BOLD=""
        NC=""
    fi
else
    # No color support
    RED=""
    GREEN=""
    YELLOW=""
    BLUE=""
    MAGENTA=""
    CYAN=""
    BOLD=""
    NC=""
fi

# Function to print colored output
cprint() {
    printf "%b%s%b\n" "$1" "$2" "${NC}"
}

# Function to display help
show_help() {
    printf "%b\n" "${CYAN}"
    cat << 'EOF'
    ╔═══════════════════════════════════════════════════════════════╗
    ║                                                               ║
    ║     ██████╗ ███╗   ██╗████████╗    ███╗   ███╗ █████╗  ██████╗   ║
    ║    ██╔═══██╗████╗  ██║╚══██╔══╝    ████╗ ████║██╔══██╗██╔════╝   ║
    ║    ██║   ██║██╔██╗ ██║   ██║       ██╔████╔██║███████║██║  ███╗  ║
    ║    ██║   ██║██║╚██╗██║   ██║       ██║╚██╔╝██║██╔══██║██║   ██║  ║
    ║    ╚██████╔╝██║ ╚████║   ██║       ██║ ╚═╝ ██║██║  ██║╚██████╔╝  ║
    ║     ╚═════╝ ╚═╝  ╚═══╝   ╚═╝       ╚═╝     ╚═╝╚═╝  ╚═╝ ╚═════╝   ║
    ║                                                               ║
    ║           Post-Polish Coverage & Quality Control              ║
    ╚═══════════════════════════════════════════════════════════════╝
EOF
    printf "%b\n\n" "${NC}"

    printf "%b%bNAME%b\n" "${BOLD}" "" "${NC}"
    printf "    coverage.sh - ONT MAG post-polish helper for computing coverage and QC\n\n"
    
    printf "%b%bSYNOPSIS%b\n" "${BOLD}" "" "${NC}"
    printf "    %bcoverage.sh%b [%bOPTIONS%b] %b-a%b <assembly> %b-r%b <reads>\n\n" "${GREEN}" "${NC}" "${YELLOW}" "${NC}" "${CYAN}" "${NC}" "${CYAN}" "${NC}"
    
    printf "%b%bDESCRIPTION%b\n" "${BOLD}" "" "${NC}"
    printf "    This script performs post-polishing analysis on ONT assemblies:\n"
    printf "    1. Computes per-contig coverage using CoverM (with ONT-optimized settings)\n"
    printf "    2. Optionally generates JGI-style depth tables for binning\n"
    printf "    3. Optionally runs CheckM2 and GUNC quality assessment on bins\n\n"

    printf "%b%bREQUIRED OPTIONS%b\n" "${BOLD}" "" "${NC}"
    printf "    %b-a, --assembly%b FILE\n" "${CYAN}" "${NC}"
    printf "        Polished contigs file (e.g., medaka_out/consensus.fasta)\n\n"
    
    printf "    %b-r, --reads%b FILE(S)\n" "${CYAN}" "${NC}"
    printf "        ONT read files (space-separated if multiple)\n"
    printf "        Example: -r \"reads1.fastq.gz reads2.fastq.gz\"\n\n"

    printf "%b%bOPTIONAL OPTIONS%b\n" "${BOLD}" "" "${NC}"
    printf "    %b-o, --output%b DIR\n" "${CYAN}" "${NC}"
    printf "        Output directory (default: coverage-out)\n\n"
    
    printf "    %b-t, --threads%b NUM\n" "${CYAN}" "${NC}"
    printf "        Number of threads to use (default: 32)\n\n"
    
    printf "    %b-b, --bams%b FILE(S)\n" "${CYAN}" "${NC}"
    printf "        Use precomputed BAM files instead of raw reads\n"
    printf "        BAMs must be primary-alignments only, sorted & indexed\n\n"
    
    printf "    %b--identity%b FLOAT\n" "${CYAN}" "${NC}"
    printf "        Minimum read percent identity for CoverM (default: 0.85)\n"
    printf "        Recommended: 0.85 for ONT SUP non-duplex\n\n"
    
    printf "    %b--aligned-frac%b FLOAT\n" "${CYAN}" "${NC}"
    printf "        Minimum fraction of read aligned to count (default: 0.50)\n\n"
    
    printf "    %b--methods%b STRING\n" "${CYAN}" "${NC}"
    printf "        CoverM methods to compute (default: \"mean covered_bases variance tpm rpkm\")\n\n"
    
    printf "    %b--jgi%b\n" "${CYAN}" "${NC}"
    printf "        Enable JGI-style depth table generation\n\n"
    
    printf "    %b--jgi-mapq%b NUM\n" "${CYAN}" "${NC}"
    printf "        Minimum MAPQ for JGI summarize (default: 5)\n\n"
    
    printf "    %b--jgi-pid%b NUM\n" "${CYAN}" "${NC}"
    printf "        Minimum percent identity for JGI (default: 80)\n\n"

    printf "%b%bBINNING QC OPTIONS%b\n" "${BOLD}" "" "${NC}"
    printf "    %b--bins-dir%b DIR\n" "${CYAN}" "${NC}"
    printf "        Directory containing MAG FASTA files for QC\n\n"
    
    printf "    %b--bins-ext%b EXT\n" "${CYAN}" "${NC}"
    printf "        Extension of bin files (default: fa)\n"
    printf "        Options: fa, fna, fasta\n\n"
    
    printf "    %b--checkm2-db%b PATH\n" "${CYAN}" "${NC}"
    printf "        Path to CheckM2 database (enables CheckM2 if set)\n"
    printf "        Example: /data/checkm2_db_2024_05\n\n"
    
    printf "    %b--gunc-db%b PATH\n" "${CYAN}" "${NC}"
    printf "        Path to GUNC database (enables GUNC if set)\n"
    printf "        Example: /data/gunc_db_progenomes2.1\n\n"
    
    printf "    %b--no-color%b\n" "${CYAN}" "${NC}"
    printf "        Disable colored output (useful for logs/pipes)\n\n"
    
    printf "    %b-h, --help%b\n" "${CYAN}" "${NC}"
    printf "        Show this help message and exit\n\n"

    printf "%b%bEXAMPLES%b\n" "${BOLD}" "" "${NC}"
    printf "    %b# Basic usage with single read file%b\n" "${GREEN}" "${NC}"
    printf "    ./coverage.sh -a consensus.fasta -r ont_reads.fastq.gz\n\n"

    printf "    %b# Multiple read files with custom output%b\n" "${GREEN}" "${NC}"
    printf "    ./coverage.sh -a consensus.fasta -r \"reads1.fq.gz reads2.fq.gz\" -o my_coverage\n\n"

    printf "    %b# Using precomputed BAMs with JGI depth%b\n" "${GREEN}" "${NC}"
    printf "    ./coverage.sh -a consensus.fasta -b \"sample1.bam sample2.bam\" --jgi\n\n"

    printf "    %b# Full pipeline with bin QC%b\n" "${GREEN}" "${NC}"
    printf "    ./coverage.sh -a consensus.fasta -r reads.fastq.gz \\\\\n"
    printf "        --jgi \\\\\n"
    printf "        --bins-dir bins_dastool/nonredundant_bins \\\\\n"
    printf "        --checkm2-db /data/checkm2_db \\\\\n"
    printf "        --gunc-db /data/gunc_db\n\n"

    printf "%b%bOUTPUT FILES%b\n" "${BOLD}" "" "${NC}"
    printf "    %bcoverage-out/%b\n" "${MAGENTA}" "${NC}"
    printf "    ├── %bcoverm_contig.tsv%b      # Full CoverM coverage table\n" "${BLUE}" "${NC}"
    printf "    ├── %bcontig_mean_coverage.tsv%b  # Simplified mean coverage\n" "${BLUE}" "${NC}"
    printf "    ├── %bjgi_depth.txt%b          # JGI-style depth (if --jgi)\n" "${BLUE}" "${NC}"
    printf "    ├── %bcheckm2/%b               # CheckM2 results (if enabled)\n" "${BLUE}" "${NC}"
    printf "    │   └── quality_report.tsv\n"
    printf "    └── %bgunc/%b                  # GUNC results (if enabled)\n\n" "${BLUE}" "${NC}"

    printf "%b%bDEPENDENCIES%b\n" "${BOLD}" "" "${NC}"
    printf "    Required: coverm, samtools, minimap2\n"
    printf "    Optional: checkm2, gunc, jgi_summarize_bam_contig_depths\n\n"

    printf "%b%bNOTES%b\n" "${BOLD}" "" "${NC}"
    printf "    • Optimized for Oxford Nanopore long reads\n"
    printf "    • Default settings tuned for SUP accuracy reads\n"
    printf "    • For duplex reads, consider increasing --identity to 0.95\n\n"

    printf "%b%bAUTHOR%b\n" "${BOLD}" "" "${NC}"
    printf "    ONT MAG Pipeline Contributors\n\n"

    printf "%b%bVERSION%b\n" "${BOLD}" "" "${NC}"
    printf "    2.0.0 (Flag-based interface)\n\n"
}

# Function to display short usage
show_usage() {
    cprint "${YELLOW}" "Usage: $0 -a <assembly> -r <reads> [options]"
    printf "Try '$0 --help' for more information.\n"
}

# Function to parse command line arguments
parse_args() {
    local args=("$@")
    local i=0
    
    while [[ $i -lt ${#args[@]} ]]; do
        case "${args[$i]}" in
            -h|--help)
                show_help
                exit 0
                ;;
            -a|--assembly)
                ASSEMBLY="${args[$((i+1))]}"
                i=$((i+2))
                ;;
            -r|--reads)
                READS="${args[$((i+1))]}"
                i=$((i+2))
                ;;
            -o|--output)
                OUTDIR="${args[$((i+1))]}"
                i=$((i+2))
                ;;
            -t|--threads)
                THREADS="${args[$((i+1))]}"
                i=$((i+2))
                ;;
            -b|--bams)
                BAMS="${args[$((i+1))]}"
                i=$((i+2))
                ;;
            --identity)
                COVERM_IDENTITY="${args[$((i+1))]}"
                i=$((i+2))
                ;;
            --aligned-frac)
                COVERM_ALIGNED_FRAC="${args[$((i+1))]}"
                i=$((i+2))
                ;;
            --methods)
                COVERM_METHODS="${args[$((i+1))]}"
                i=$((i+2))
                ;;
            --jgi)
                RUN_JGI=1
                i=$((i+1))
                ;;
            --jgi-mapq)
                JGI_MIN_MAPQ="${args[$((i+1))]}"
                i=$((i+2))
                ;;
            --jgi-pid)
                JGI_MIN_PID="${args[$((i+1))]}"
                i=$((i+2))
                ;;
            --bins-dir)
                BINS_DIR="${args[$((i+1))]}"
                i=$((i+2))
                ;;
            --bins-ext)
                BINS_EXT="${args[$((i+1))]}"
                i=$((i+2))
                ;;
            --checkm2-db)
                CHECKM2_DB="${args[$((i+1))]}"
                i=$((i+2))
                ;;
            --gunc-db)
                GUNC_DB="${args[$((i+1))]}"
                i=$((i+2))
                ;;
            --no-color)
                RED=""
                GREEN=""
                YELLOW=""
                BLUE=""
                MAGENTA=""
                CYAN=""
                BOLD=""
                NC=""
                i=$((i+1))
                ;;
            *)
                printf "%bError: Unknown option %s%b\n" "${RED}" "${args[$i]}" "${NC}" >&2
                show_usage
                exit 1
                ;;
        esac
    done
}

# Main script starts here
main() {
    # Parse command line arguments
    if [[ $# -eq 0 ]]; then
        show_usage
        exit 1
    fi
    
    parse_args "$@"
    
    # Validate required arguments
    if [[ -z "$ASSEMBLY" ]]; then
        printf "%bError: Assembly file (-a/--assembly) is required%b\n" "${RED}" "${NC}" >&2
        show_usage
        exit 1
    fi
    
    if [[ -z "$READS" ]] && [[ -z "$BAMS" ]]; then
        printf "%bError: Either reads (-r/--reads) or BAMs (-b/--bams) must be provided%b\n" "${RED}" "${NC}" >&2
        show_usage
        exit 1
    fi
    
    # Tool checks
    need() { 
        command -v "$1" >/dev/null 2>&1 || { 
            printf "%bERROR: '%s' not found in PATH%b\n" "${RED}" "$1" "${NC}" >&2
            exit 127
        }
    }
    
    for t in coverm samtools minimap2; do need "$t"; done
    
    [[ -f "$ASSEMBLY" ]] || { 
        printf "%bERROR: Assembly file not found: %s%b\n" "${RED}" "$ASSEMBLY" "${NC}" >&2
        exit 1
    }
    
    mkdir -p "$OUTDIR"
    LOG="$OUTDIR/run.log"
    
    # Pretty header
    printf "%b╔═══════════════════════════════════════════════════════════════╗%b\n" "${CYAN}" "${NC}" | tee "$LOG"
    printf "%b║         %bPost-Polish Coverage & QC Analysis%b%b          ║%b\n" "${CYAN}" "${BOLD}" "${NC}" "${CYAN}" "${NC}" | tee -a "$LOG"
    printf "%b╚═══════════════════════════════════════════════════════════════╝%b\n" "${CYAN}" "${NC}" | tee -a "$LOG"
    printf "\n" | tee -a "$LOG"
    
    printf "%bConfiguration:%b\n" "${GREEN}" "${NC}" | tee -a "$LOG"
    printf "  %bAssembly:%b %s\n" "${BOLD}" "${NC}" "$ASSEMBLY" | tee -a "$LOG"
    printf "  %bThreads:%b  %s\n" "${BOLD}" "${NC}" "$THREADS" | tee -a "$LOG"
    printf "  %bReads:%b    %s %s %s\n" "${BOLD}" "${NC}" "${BAMS:+<using BAMs>}" "${BAMS:+"$BAMS"}" "${READS:+"$READS"}" | tee -a "$LOG"
    printf "  %bOutput:%b   %s\n" "${BOLD}" "${NC}" "$OUTDIR" | tee -a "$LOG"
    printf "\n" | tee -a "$LOG"
    
    # ---------- 1) Contig coverage with CoverM ----------
    COVERM_OUT="$OUTDIR/coverm_contig.tsv"
    printf "%b[CoverM]%b Computing per-contig depths → %s\n" "${BLUE}" "${NC}" "$COVERM_OUT" | tee -a "$LOG"
    
    if [[ -n "$BAMS" ]]; then
        # Use precomputed BAMs
        coverm contig \
            --bam-files $BAMS \
            -r "$ASSEMBLY" \
            --methods $COVERM_METHODS \
            --threads "$THREADS" \
            --min-read-percent-identity "$COVERM_IDENTITY" \
            --min-read-aligned-percent "$COVERM_ALIGNED_FRAC" \
            --output-file "$COVERM_OUT"
    else
        # Map internally with minimap2 using ONT preset
        coverm contig \
            --reads $READS \
            -r "$ASSEMBLY" \
            --mapper minimap2 \
            --minimap2-args "-x map-ont --secondary=no" \
            --methods $COVERM_METHODS \
            --threads "$THREADS" \
            --min-read-percent-identity "$COVERM_IDENTITY" \
            --min-read-aligned-percent "$COVERM_ALIGNED_FRAC" \
            --output-file "$COVERM_OUT"
    fi
    
    # Create simplified coverage matrix
    awk 'BEGIN{FS=OFS="\t"} NR==1{print $1,$(NF-3)} NR>1{print $1,$(NF-3)}' "$COVERM_OUT" > "$OUTDIR/contig_mean_coverage.tsv"
    printf "  %b✔%b Coverage tables generated\n" "${GREEN}" "${NC}" | tee -a "$LOG"
    
    # ---------- (optional) 1b) JGI-style depth ----------
    if [[ "$RUN_JGI" == "1" ]]; then
        need jgi_summarize_bam_contig_depths
        printf "%b[JGI]%b Generating ONT-tuned depth table (MAPQ≥%s, PID≥%s%%)\n" "${BLUE}" "${NC}" "$JGI_MIN_MAPQ" "$JGI_MIN_PID" | tee -a "$LOG"
        
        # Build/collect BAMs if not provided
        if [[ -z "$BAMS" ]]; then
            BAMDIR="$OUTDIR/jgi_bams"
            mkdir -p "$BAMDIR"
            i=1
            for R in $READS; do
                BAM="$BAMDIR/sample${i}.bam"
                printf "  Mapping %s → %s\n" "$R" "$BAM" | tee -a "$LOG"
                minimap2 -t "$THREADS" -x map-ont --secondary=no "$ASSEMBLY" "$R" \
                    | samtools sort -@ "$THREADS" -o "$BAM" -
                samtools index -@ "$THREADS" "$BAM"
                i=$((i+1))
            done
            BAMS="$(printf '%s ' $BAMDIR/*.bam)"
        fi
        
        # JGI summarize
        JGI_OUT="$OUTDIR/jgi_depth.txt"
        jgi_summarize_bam_contig_depths \
            --outputDepth "$JGI_OUT" \
            --minMapQual "$JGI_MIN_MAPQ" \
            --percentIdentity "$JGI_MIN_PID" \
            $BAMS
        
        printf "  %b✔%b JGI depth table → %s\n" "${GREEN}" "${NC}" "$JGI_OUT" | tee -a "$LOG"
    fi
    
    # ---------- 2) Optional QC on bins ----------
    if [[ -n "$BINS_DIR" ]]; then
        printf "\n" | tee -a "$LOG"
        printf "%b[QC]%b Running quality control on bins\n" "${MAGENTA}" "${NC}" | tee -a "$LOG"
        [[ -d "$BINS_DIR" ]] || { 
            printf "%bERROR: BINS_DIR is not a directory: %s%b\n" "${RED}" "$BINS_DIR" "${NC}" >&2
            exit 1
        }
        
        # CheckM2
        if [[ -n "$CHECKM2_DB" ]]; then
            need checkm2
            printf "%b[CheckM2]%b Assessing genome quality (*.%s)\n" "${BLUE}" "${NC}" "$BINS_EXT" | tee -a "$LOG"
            CHECKM2_OUT="$OUTDIR/checkm2"
            mkdir -p "$CHECKM2_OUT"
            checkm2 predict \
                --threads "$THREADS" \
                --database_path "$CHECKM2_DB" \
                --extension "$BINS_EXT" \
                --input "$BINS_DIR" \
                --output_directory "$CHECKM2_OUT" \
                --remove_intermediates
            printf "  %b✔%b CheckM2 summary → %s/quality_report.tsv\n" "${GREEN}" "${NC}" "$CHECKM2_OUT" | tee -a "$LOG"
        else
            printf "  %b⚠%b CheckM2 skipped (set --checkm2-db to enable)\n" "${YELLOW}" "${NC}" | tee -a "$LOG"
        fi
        
        # GUNC
        if [[ -n "$GUNC_DB" ]]; then
            need gunc
            printf "%b[GUNC]%b Assessing chimerism/contamination\n" "${BLUE}" "${NC}" | tee -a "$LOG"
            GUNC_OUT="$OUTDIR/gunc"
            mkdir -p "$GUNC_OUT"
            gunc run \
                --input_dir "$BINS_DIR" \
                --input_pattern "*.$BINS_EXT" \
                --threads "$THREADS" \
                --db_file "$GUNC_DB" \
                --output_dir "$GUNC_OUT" \
                --detailed_output
            printf "  %b✔%b GUNC outputs → %s\n" "${GREEN}" "${NC}" "$GUNC_OUT" | tee -a "$LOG"
        else
            printf "  %b⚠%b GUNC skipped (set --gunc-db to enable)\n" "${YELLOW}" "${NC}" | tee -a "$LOG"
        fi
    else
        printf "%b[QC]%b No bins directory provided; skipping CheckM2/GUNC\n" "${YELLOW}" "${NC}" | tee -a "$LOG"
    fi
    
    # Summary
    printf "\n" | tee -a "$LOG"
    printf "%b╔═══════════════════════════════════════════════════════════════╗%b\n" "${CYAN}" "${NC}" | tee -a "$LOG"
    printf "%b║                  %bPipeline Complete!%b%b                  ║%b\n" "${CYAN}" "${BOLD}" "${NC}" "${CYAN}" "${NC}" | tee -a "$LOG"
    printf "%b╚═══════════════════════════════════════════════════════════════╝%b\n" "${CYAN}" "${NC}" | tee -a "$LOG"
    printf "\n" | tee -a "$LOG"
    printf "%bOutput Summary:%b\n" "${GREEN}" "${NC}" | tee -a "$LOG"
    printf "  %bCoverM table:%b        %s\n" "${BOLD}" "${NC}" "$COVERM_OUT" | tee -a "$LOG"
    printf "  %bMean coverage:%b      %s/contig_mean_coverage.tsv\n" "${BOLD}" "${NC}" "$OUTDIR" | tee -a "$LOG"
    [[ "$RUN_JGI" == "1" ]] && printf "  %bJGI depth:%b          %s\n" "${BOLD}" "${NC}" "$JGI_OUT" | tee -a "$LOG"
    [[ -n "$BINS_DIR" ]] && printf "  %bQC results:%b         %s/{checkm2,gunc}\n" "${BOLD}" "${NC}" "$OUTDIR" | tee -a "$LOG"
    printf "\n" | tee -a "$LOG"
}

# Run main function
main "$@"
