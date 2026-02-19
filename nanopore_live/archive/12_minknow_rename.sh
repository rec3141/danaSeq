#!/bin/bash

# Configuration
SOURCE_BASE="/var/lib/minknow/data/basecalling/pass"
DEST_BASE="/data/minknow"  # Will be updated with project directory

# Arrays to store mappings
declare -A RUNID_TO_FLOWCELL
declare -A RUNID_TO_DESTDIR

# Function to build the mapping of run IDs to flowcell IDs
build_mapping() {
    echo "Scanning destination directories to build run ID/flowcell mapping..."
    
    # Find all run directories
    for run_dir in "$DEST_BASE"/*; do
        if [ -d "$run_dir" ]; then
            local run_name=$(basename "$run_dir")
            
            # Look for subdirectories that match the expected pattern
            for subdir in "$run_dir"/*; do
                if [ -d "$subdir" ]; then
                    local subdir_name=$(basename "$subdir")
                    
                    # Extract flowcell ID and run ID from directory name
                    # Expected pattern: YYYYMMDD_HHMM_MD-XXXXXX_FLOWCELLID_RUNID
                    if [[ "$subdir_name" =~ ([^_]+)_([^_]+)_([^_]+)_([^_]+)_([^_]+)$ ]]; then
                        local flowcell_id="${BASH_REMATCH[4]}"
                        local full_run_id="${BASH_REMATCH[5]}"
                        
                        # Extract short run ID (first 8 characters)
                        local short_run_id="${full_run_id:0:8}"
                        
                        RUNID_TO_FLOWCELL["$short_run_id"]="$flowcell_id"
                        RUNID_TO_DESTDIR["$short_run_id"]="$subdir/fastq_pass"
                        
                        echo "  Found: $short_run_id -> $flowcell_id (in $run_name)"
                    fi
                fi
            done
        fi
    done
    
    if [ ${#RUNID_TO_FLOWCELL[@]} -eq 0 ]; then
        echo "Error: No valid destination directories found!"
        exit 1
    fi
    
    echo "Built mapping for ${#RUNID_TO_FLOWCELL[@]} run ID/flowcell pairs."
    echo ""
}

# Function to extract run ID from filename
extract_run_id_from_filename() {
    local filename="$1"
    # Extract run ID from pattern: fastq_runid_XXXXXXXX-XXXX-XXXX-XXXX-XXXXXXXXXXXX_NUMBER_0.fastq.gz
    if [[ "$filename" =~ fastq_runid_([a-f0-9]{8})-[a-f0-9]{4}-[a-f0-9]{4}-[a-f0-9]{4}-[a-f0-9]{12}_([0-9]+)_0\.fastq\.gz ]]; then
        echo "${BASH_REMATCH[1]}"  # Return the 8-character run ID
    else
        echo ""
    fi
}

# Function to get all available barcodes in source
get_source_barcodes() {
    for barcode_dir in "$SOURCE_BASE"/barcode*; do
        if [ -d "$barcode_dir" ]; then
            basename "$barcode_dir"
        fi
    done
}

# Function to process files for a specific barcode
process_barcode() {
    local barcode="$1"
    local dry_run="$2"
    
    local source_dir="$SOURCE_BASE/$barcode"
    
    if [ ! -d "$source_dir" ]; then
        echo "  Warning: Source directory not found: $source_dir"
        return
    fi
    
    local file_count=0
    local processed_count=0
    
    echo "  Processing $barcode..."
    
    # Process each fastq.gz file in the barcode directory
    for file in "$source_dir"/fastq_runid_*_*_0.fastq.gz; do
        if [ ! -f "$file" ]; then
            continue
        fi
        
        ((file_count++))
        
        local filename=$(basename "$file")
        local run_id=$(extract_run_id_from_filename "$filename")
        
        if [ -z "$run_id" ]; then
            echo "    Warning: Could not extract run ID from $filename"
            continue
        fi
        
        local flowcell_id="${RUNID_TO_FLOWCELL[$run_id]}"
        
        if [ -z "$flowcell_id" ]; then
            echo "    Warning: No matching flowcell found for run ID $run_id in $filename"
            continue
        fi
        
        local dest_base_dir="${RUNID_TO_DESTDIR[$run_id]}"
        local dest_dir="$dest_base_dir/$barcode"
        
        # Extract number from filename
        if [[ "$filename" =~ _([0-9]+)_0\.fastq\.gz$ ]]; then
            local number="${BASH_REMATCH[1]}"
        else
            echo "    Warning: Could not extract number from $filename"
            continue
        fi
        
        # Construct new filename
        local new_filename="${flowcell_id}_pass_${barcode}_${run_id}_skipped_${number}.fastq.gz"
        local dest_file="$dest_dir/$new_filename"
        
        if [ "$dry_run" = true ]; then
            echo "    [DRY RUN] $filename -> $new_filename"
            echo "              Destination: $dest_dir"
        else
            # Create destination directory if it doesn't exist
            mkdir -p "$dest_dir"
            
            # Copy the file
            if cp "$file" "$dest_file"; then
                ((processed_count++))
                echo "    ✓ $filename -> $new_filename"
            else
                echo "    ✗ Failed to copy $filename"
            fi
        fi
    done
    
    if [ "$dry_run" = true ]; then
        echo "    Found $file_count files in $barcode"
    else
        echo "    Processed $processed_count/$file_count files in $barcode"
    fi
}

# Main script
main() {
    local dry_run=false
    local specific_barcode=""
    local project_dir=""
    
    # Parse command line arguments
    while [[ $# -gt 0 ]]; do
        case $1 in
            --dry-run|-n)
                dry_run=true
                shift
                ;;
            --barcode|-b)
                specific_barcode="$2"
                shift 2
                ;;
            --project|-p)
                project_dir="$2"
                shift 2
                ;;
            --help|-h)
                echo "Usage: $0 [OPTIONS]"
                echo ""
                echo "Options:"
                echo "  --project PROJECT, -p   Specify the project directory (e.g., QEI2025)"
                echo "  --dry-run, -n           Preview what will be done without copying files"
                echo "  --barcode BARCODE, -b   Process only the specified barcode (e.g., barcode17)"
                echo "  --help, -h              Show this help message"
                echo ""
                echo "Examples:"
                echo "  $0 -p QEI2025 --dry-run        # Preview all barcodes for QEI2025 project"
                echo "  $0 -p QEI2025 -b barcode17     # Process only barcode17 for QEI2025"
                echo "  $0 -p QEI2026 -n -b barcode17  # Preview barcode17 for QEI2026 project"
                exit 0
                ;;
            *)
                echo "Unknown option: $1"
                echo "Use --help for usage information"
                exit 1
                ;;
        esac
    done
    
    # Check if project directory is specified
    if [ -z "$project_dir" ]; then
        echo "Error: Project directory is required!"
        echo "Use --project PROJECT or -p PROJECT to specify the project directory"
        echo "Example: $0 -p QEI2025"
        echo "Use --help for more information"
        exit 1
    fi
    
    # Update the destination base path with the project directory
    DEST_BASE="$DEST_BASE/$project_dir"
    
    # Check if the project directory exists
    if [ ! -d "$DEST_BASE" ]; then
        echo "Error: Project directory does not exist: $DEST_BASE"
        echo "Available projects in /data/minknow/:"
        ls -1 /data/minknow/ 2>/dev/null | grep -E '^[A-Z0-9]+[0-9]{4}$' || echo "  (none found)"
        exit 1
    fi
    
    if [ "$dry_run" = true ]; then
        echo "=== DRY RUN MODE - No files will be copied ==="
        echo ""
    fi
    
    echo "Project: $project_dir"
    echo "Destination base: $DEST_BASE"
    echo ""
    
    # Build the flowcell/run ID mapping
    build_mapping
    
    # Show the mapping
    echo "Run ID -> Flowcell mapping:"
    for run_id in "${!RUNID_TO_FLOWCELL[@]}"; do
        echo "  $run_id -> ${RUNID_TO_FLOWCELL[$run_id]}"
    done
    echo ""
    
    # Get barcodes to process
    local barcodes_to_process
    if [ -n "$specific_barcode" ]; then
        barcodes_to_process=("$specific_barcode")
        echo "Processing specific barcode: $specific_barcode"
    else
        readarray -t barcodes_to_process < <(get_source_barcodes)
        echo "Processing all available barcodes: ${barcodes_to_process[*]}"
    fi
    echo ""
    
    # Process each barcode
    local total_barcodes=${#barcodes_to_process[@]}
    local barcode_num=0
    
    for barcode in "${barcodes_to_process[@]}"; do
        ((barcode_num++))
        echo "[$barcode_num/$total_barcodes] Processing $barcode"
        process_barcode "$barcode" "$dry_run"
        echo ""
    done
    
    if [ "$dry_run" = true ]; then
        echo "=== DRY RUN COMPLETE ==="
        echo "Project: $project_dir"
        echo "Run without --dry-run to actually copy the files"
    else
        echo "=== BATCH COPY COMPLETE ==="
        echo "Project: $project_dir"
    fi
}

# Run the main function with all arguments
main "$@"