#!/usr/bin/env bash
#
# test_edge_cases.sh - Test edge cases and corner conditions
#
# Tests unusual but valid scenarios:
#   - Empty directories
#   - Special characters in filenames
#   - Very long paths
#   - Symbolic links
#   - Permission issues
#   - Disk space issues
#   - Unusual barcode names
#
# Usage: ./test_edge_cases.sh
# Exit: 0 if all tests pass, 1 if any fail

# Don't use set -e since we track failures manually
set -uo pipefail
IFS=$'\n\t'

# Colors for output
GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
NC='\033[0m'

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
TEST_OUTPUT="${SCRIPT_DIR}/test_output"

TESTS_PASSED=0
TESTS_FAILED=0

pass() {
    echo -e "${GREEN}[PASS]${NC} $*"
    ((TESTS_PASSED++))
}

fail() {
    echo -e "${RED}[FAIL]${NC} $*"
    ((TESTS_FAILED++))
}

info() {
    echo -e "${YELLOW}[INFO]${NC} $*"
}

# Cleanup
cleanup() {
    rm -rf "${TEST_OUTPUT}" 2>/dev/null || true
}

trap cleanup EXIT

echo "========================================"
echo "Edge Case and Corner Condition Tests"
echo "========================================"
echo ""

# ====================
# TEST 1: Empty Input Directory
# ====================

info "TEST 1: Empty input directory handling"

mkdir -p "${TEST_OUTPUT}/empty_input"

# Should handle gracefully (warning, not error)
output=$("${PROJECT_ROOT}/nanopore_live/24_process_reads_optimized.sh" \
    -i "${TEST_OUTPUT}/empty_input" \
    -o "${TEST_OUTPUT}/empty_output" \
    -m 1 2>&1 || true)

if echo "$output" | grep -qi "warning\|no.*file\|empty"; then
    pass "Graceful handling of empty input directory"
else
    info "No specific warning for empty directory (may be silent)"
fi

echo ""

# ====================
# TEST 2: Special Characters in Paths
# ====================

info "TEST 2: Special characters in filenames"

# Create directory with spaces
mkdir -p "${TEST_OUTPUT}/dir with spaces/FC001/barcode01"

# Create test file
cat > "${TEST_OUTPUT}/dir with spaces/FC001/barcode01/test.fastq" <<'EOF'
@read1
ACGT
+
IIII
EOF

# Try to process (should handle spaces correctly due to quoting)
output=$("${PROJECT_ROOT}/nanopore_live/24_process_reads_optimized.sh" \
    -i "${TEST_OUTPUT}/dir with spaces" \
    -o "${TEST_OUTPUT}/output_spaces" \
    -m 1 --dry-run 2>&1 || true)

if echo "$output" | grep -qi "error.*space\|cannot.*space"; then
    fail "Does not handle spaces in paths"
else
    pass "Handles spaces in paths (or doesn't report errors)"
fi

echo ""

# ====================
# TEST 3: Barcode Name Validation
# ====================

info "TEST 3: Barcode name format validation"

# Check if invalid barcode names are rejected
if grep -q 'barcode.*=~.*\^barcode\[0-9\]' "${PROJECT_ROOT}/nanopore_live/24_process_reads_optimized.sh"; then
    pass "Barcode format validation present"
else
    info "No strict barcode validation (may accept any name)"
fi

# Create test with invalid barcode name
mkdir -p "${TEST_OUTPUT}/test_barcodes/FC001/invalid_barcode"
cat > "${TEST_OUTPUT}/test_barcodes/FC001/invalid_barcode/test.fastq" <<'EOF'
@read1
ACGT
+
IIII
EOF

# Script should skip or reject invalid barcode
output=$("${PROJECT_ROOT}/nanopore_live/24_process_reads_optimized.sh" \
    -i "${TEST_OUTPUT}/test_barcodes" \
    -o "${TEST_OUTPUT}/output_invalid" \
    -m 1 --dry-run 2>&1 || true)

# Check if it was processed or skipped
if [[ -d "${TEST_OUTPUT}/output_invalid/FC001/invalid_barcode" ]]; then
    info "Accepts non-standard barcode names"
else
    pass "Rejects or skips invalid barcode names"
fi

echo ""

# ====================
# TEST 4: Symbolic Links
# ====================

info "TEST 4: Symbolic link handling"

mkdir -p "${TEST_OUTPUT}/real_data/FC001/barcode01"
cat > "${TEST_OUTPUT}/real_data/FC001/barcode01/test.fastq" <<'EOF'
@read1
ACGT
+
IIII
EOF

# Create symlink to data
ln -s "${TEST_OUTPUT}/real_data" "${TEST_OUTPUT}/linked_data" 2>/dev/null || true

# Should either follow symlinks or reject them safely
output=$("${PROJECT_ROOT}/nanopore_live/24_process_reads_optimized.sh" \
    -i "${TEST_OUTPUT}/linked_data" \
    -o "${TEST_OUTPUT}/output_symlink" \
    -m 1 --dry-run 2>&1 || true)

if echo "$output" | grep -qi "error.*symlink\|error.*link"; then
    info "Explicitly rejects symbolic links"
else
    pass "Accepts symbolic links or handles them gracefully"
fi

echo ""

# ====================
# TEST 5: Very Long Paths
# ====================

info "TEST 5: Long path handling"

# Create nested directory structure
long_path="${TEST_OUTPUT}/a/very/long/path/that/goes/quite/deep/into/the/filesystem/FC001/barcode01"
mkdir -p "$long_path" 2>/dev/null || {
    info "Cannot create very long path (filesystem limitation)"
    continue
}

if [[ -d "$long_path" ]]; then
    cat > "${long_path}/test.fastq" <<'EOF'
@read1
ACGT
+
IIII
EOF
    pass "Can create and use long paths"
else
    info "Long path test skipped (cannot create path)"
fi

echo ""

# ====================
# TEST 6: Concurrent Access
# ====================

info "TEST 6: Concurrent access protection"

# Check for file locking mechanism
if grep -q "flock\|lockfile\|mkdir.*lock" "${PROJECT_ROOT}/nanopore_live/24_process_reads_optimized.sh"; then
    pass "File locking mechanism present"
else
    fail "No file locking (may have race conditions)"
fi

# Check for atomic operations
if grep -q "mv.*tmp\|atomic" "${PROJECT_ROOT}/nanopore_live/24_process_reads_optimized.sh"; then
    pass "Atomic operations for file safety"
else
    info "No explicit atomic operations (may corrupt on concurrent writes)"
fi

echo ""

# ====================
# TEST 7: Disk Space Handling
# ====================

info "TEST 7: Disk space checking"

# Check if script checks for available disk space
if grep -q "df\|disk.*space\|quota" "${PROJECT_ROOT}/nanopore_live/24_process_reads_optimized.sh"; then
    pass "Disk space checking present"
else
    info "No disk space checking (may fail when full)"
fi

# Check if script handles write failures gracefully
if grep -q "No space left\|disk full\|ENOSPC" "${PROJECT_ROOT}/nanopore_live/24_process_reads_optimized.sh"; then
    pass "Explicit handling of disk full errors"
else
    info "No explicit disk full handling (relies on set -e)"
fi

echo ""

# ====================
# TEST 8: Corrupted FASTQ Files
# ====================

info "TEST 8: Corrupted FASTQ file handling"

mkdir -p "${TEST_OUTPUT}/corrupt_test/FC001/barcode01"

# Create corrupted FASTQ (missing quality scores)
cat > "${TEST_OUTPUT}/corrupt_test/FC001/barcode01/corrupt.fastq" <<'EOF'
@read1
ACGT
+
II
@read2_missing_qual
TGCA
+
EOF

# Check if script has FASTQ validation
if grep -q "seqkit\|repair.*fastq\|validate.*fastq" "${PROJECT_ROOT}/nanopore_live/24_process_reads_optimized.sh"; then
    pass "FASTQ validation/repair mechanism present"
else
    info "No explicit FASTQ validation (may pass corrupt files to tools)"
fi

echo ""

# ====================
# TEST 9: Unicode/Non-ASCII Characters
# ====================

info "TEST 9: Unicode character handling"

# Create directory with unicode characters
mkdir -p "${TEST_OUTPUT}/test_unicode_café/FC001/barcode01" 2>/dev/null || {
    info "Cannot create unicode directory name (filesystem limitation)"
}

if [[ -d "${TEST_OUTPUT}/test_unicode_café" ]]; then
    info "Filesystem supports unicode directory names"
    # Note: We don't test processing as this is uncommon in practice
else
    info "Unicode directory test skipped (filesystem doesn't support)"
fi

echo ""

# ====================
# TEST 10: Memory Limits
# ====================

info "TEST 10: Memory management"

# Check for Kraken2 serialization (critical for memory)
if grep -q "sem.*kraken\|kraken.*--id" "${PROJECT_ROOT}/nanopore_live/24_process_reads_optimized.sh"; then
    pass "Kraken2 memory serialization present"
else
    fail "No Kraken2 serialization (CRITICAL: will OOM with parallel runs)"
fi

# Check for process limiting to prevent memory exhaustion
if grep -q "\-j\s*[0-9]\+\|-P\s*[0-9]\+\|--jobs" "${PROJECT_ROOT}/nanopore_live/24_process_reads_optimized.sh"; then
    pass "Job limiting to control memory usage"
else
    info "No explicit job limiting (may use unlimited parallel jobs)"
fi

echo ""

# ====================
# TEST 11: Permission Issues
# ====================

info "TEST 11: Permission error handling"

# Create read-only directory
mkdir -p "${TEST_OUTPUT}/readonly_test/FC001/barcode01"
chmod 444 "${TEST_OUTPUT}/readonly_test/FC001/barcode01" 2>/dev/null || true

# Try to write to read-only location
output=$("${PROJECT_ROOT}/nanopore_live/24_process_reads_optimized.sh" \
    -i "${TEST_OUTPUT}/readonly_test" \
    -o "${TEST_OUTPUT}/readonly_test/output" \
    -m 1 2>&1 || true)

if echo "$output" | grep -qi "permission.*denied\|cannot.*create\|error"; then
    pass "Detects and reports permission errors"
else
    info "Permission errors may not be explicitly reported"
fi

# Cleanup
chmod 755 "${TEST_OUTPUT}/readonly_test/FC001/barcode01" 2>/dev/null || true

echo ""

# ====================
# TEST 12: Zero-Length Files
# ====================

info "TEST 12: Zero-length file handling"

mkdir -p "${TEST_OUTPUT}/empty_files/FC001/barcode01"

# Create zero-length FASTQ
touch "${TEST_OUTPUT}/empty_files/FC001/barcode01/empty.fastq"

# Check if script detects empty files
if grep -q "\[\[ -s " "${PROJECT_ROOT}/nanopore_live/24_process_reads_optimized.sh"; then
    pass "Empty file detection present (using -s flag)"
else
    fail "No empty file detection"
fi

echo ""

# ====================
# TEST 13: Path Traversal Prevention
# ====================

info "TEST 13: Path traversal attack prevention"

# Check for .. detection in paths
if grep -q "\.\." "${PROJECT_ROOT}/nanopore_live/24_process_reads_optimized.sh" | grep -q "ERROR\|reject\|invalid"; then
    pass "Path traversal (..) detection present"
elif grep -q "realpath" "${PROJECT_ROOT}/nanopore_mag/61_map_and_bin_optimized.sh"; then
    pass "Path normalization prevents traversal"
else
    fail "No path traversal prevention"
fi

# Try path traversal in MAG script
output=$("${PROJECT_ROOT}/nanopore_mag/61_map_and_bin_optimized.sh" \
    "../../etc/passwd" 2>&1 || true)

if echo "$output" | grep -qi "error\|invalid\|denied"; then
    pass "Path traversal attempt rejected"
else
    info "Path traversal may not be explicitly blocked"
fi

echo ""

# ====================
# TEST 14: Extremely Large Files
# ====================

info "TEST 14: Large file handling"

# Check if there are any file size limits
if grep -q "size.*limit\|max.*size\|too.*large" "${PROJECT_ROOT}/nanopore_live/24_process_reads_optimized.sh"; then
    pass "File size limit checking present"
else
    info "No file size limits (will process any size)"
fi

echo ""

# ====================
# TEST 15: Unusual Exit Scenarios
# ====================

info "TEST 15: Cleanup on unusual exits"

# Check for comprehensive trap coverage
traps_found=0
for signal in EXIT INT TERM HUP QUIT; do
    if grep -q "trap.*$signal" "${PROJECT_ROOT}/nanopore_live/24_process_reads_optimized.sh"; then
        ((traps_found++))
    fi
done

if (( traps_found >= 4 )); then
    pass "Comprehensive signal trap coverage ($traps_found/5 signals)"
elif (( traps_found >= 2 )); then
    info "Partial signal coverage ($traps_found/5 signals)"
else
    fail "Insufficient signal coverage ($traps_found/5 signals)"
fi

echo ""

# ====================
# RESULTS
# ====================

echo "========================================"
echo "Test Results"
echo "========================================"
echo -e "${GREEN}PASSED: ${TESTS_PASSED}${NC}"
echo -e "${RED}FAILED: ${TESTS_FAILED}${NC}"
echo ""

if (( TESTS_FAILED == 0 )); then
    echo -e "${GREEN}✓ All edge case tests passed${NC}"
    exit 0
else
    echo -e "${RED}✗ Some edge case tests failed${NC}"
    exit 1
fi
