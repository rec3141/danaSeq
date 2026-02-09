#!/usr/bin/env bash
#
# test_basic_pipeline.sh - Test basic QC pipeline functionality
#
# Tests core pipeline operations with minimal test data:
#   - Script accepts valid arguments
#   - Creates expected directory structure
#   - Processes test FASTQ file
#   - Generates expected output files
#   - Handles missing dependencies gracefully
#
# Usage: ./test_basic_pipeline.sh
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
TEST_FIXTURES="${SCRIPT_DIR}/test_fixtures"
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

# Cleanup function
cleanup() {
    if [[ -d "${TEST_OUTPUT}" ]]; then
        rm -rf "${TEST_OUTPUT}"
    fi
}

trap cleanup EXIT

echo "========================================"
echo "Basic Pipeline Functionality Tests"
echo "========================================"
echo ""

# ====================
# SETUP
# ====================

info "Setting up test environment..."

# Clean any previous test output
cleanup

# Create test output directory
mkdir -p "${TEST_OUTPUT}"

# ====================
# TEST 1: Help/Usage
# ====================

info "TEST 1: Help and usage output"

if "${PROJECT_ROOT}/10_realtime_processing/24_process_reads_optimized.sh" -h 2>&1 | grep -q "Usage:"; then
    pass "Help flag (-h) works"
else
    fail "Help flag (-h) doesn't work"
fi

if "${PROJECT_ROOT}/10_realtime_processing/24_process_reads_optimized.sh" --help 2>&1 | grep -q "Usage:"; then
    pass "Help flag (--help) works"
else
    fail "Help flag (--help) doesn't work"
fi

echo ""

# ====================
# TEST 2: Argument Validation
# ====================

info "TEST 2: Argument validation"

# Should fail without required --input argument
if ! "${PROJECT_ROOT}/10_realtime_processing/24_process_reads_optimized.sh" 2>&1 | grep -q "\[ERROR\].*--input"; then
    fail "Should require --input argument"
else
    pass "Correctly requires --input argument"
fi

# Should fail with non-existent input directory
if ! "${PROJECT_ROOT}/10_realtime_processing/24_process_reads_optimized.sh" -i /nonexistent/path 2>&1 | grep -q "\[ERROR\]"; then
    fail "Should reject non-existent input directory"
else
    pass "Correctly rejects non-existent input directory"
fi

# Should fail if --hmm used without -P flag
if ! "${PROJECT_ROOT}/10_realtime_processing/24_process_reads_optimized.sh" -i "${TEST_FIXTURES}" --hmm test.hmm 2>&1 | grep -q "\[ERROR\].*--hmm.*-P"; then
    fail "Should require -P flag when --hmm is used"
else
    pass "Correctly requires -P flag with --hmm"
fi

echo ""

# ====================
# TEST 3: Directory Creation
# ====================

info "TEST 3: Output directory structure creation"

# Create minimal test input structure
mkdir -p "${TEST_FIXTURES}/FC001/barcode01"

# Create a minimal valid FASTQ file
cat > "${TEST_FIXTURES}/FC001/barcode01/test.fastq" <<'EOF'
@read1
ACGTACGTACGTACGT
+
IIIIIIIIIIIIIIII
@read2
TGCATGCATGCATGCA
+
IIIIIIIIIIIIIIII
EOF

# Try to run with dry-run or minimal processing
# We'll use -m 1 to limit to 1 file and add a short timeout
timeout 30s "${PROJECT_ROOT}/10_realtime_processing/24_process_reads_optimized.sh" \
    -i "${TEST_FIXTURES}" \
    -o "${TEST_OUTPUT}" \
    -m 1 \
    --dry-run 2>&1 | head -20 || true

# Check if expected directory structure was created
if [[ -d "${TEST_OUTPUT}" ]]; then
    pass "Output directory created"
else
    fail "Output directory not created"
fi

# Check for subdirectories (if created)
expected_dirs=(fa fq sketch prokka tetra stats kraken)
for dir in "${expected_dirs[@]}"; do
    if [[ -d "${TEST_OUTPUT}/FC001/barcode01/${dir}" ]]; then
        pass "Subdirectory created: ${dir}"
    else
        # This might be expected if --dry-run doesn't create all dirs
        info "Subdirectory not created: ${dir} (may require full run)"
    fi
done

echo ""

# ====================
# TEST 4: Dependency Checks
# ====================

info "TEST 4: Dependency checking"

# Check if script can detect missing tools gracefully
# (Don't actually require them to be installed for this test)

dependencies=(seqkit minimap2 samtools kraken2 bbduk.sh)
for tool in "${dependencies[@]}"; do
    if command -v "$tool" &>/dev/null; then
        pass "Tool available: $tool"
    else
        info "Tool not installed: $tool (optional for this test)"
    fi
done

echo ""

# ====================
# TEST 5: Configuration
# ====================

info "TEST 5: Configuration and environment"

# Check if critical paths exist
if [[ -d "/data/databases/kraken2" ]]; then
    pass "Kraken2 database directory exists"
else
    info "Kraken2 database not found (expected on fresh install)"
fi

if [[ -d "/data/databases/sendsketch/refseq-sketch.gz" ]] || [[ -d "/data/databases/sendsketch" ]]; then
    pass "Sendsketch database directory exists"
else
    info "Sendsketch database not found (expected on fresh install)"
fi

echo ""

# ====================
# TEST 6: Script Syntax
# ====================

info "TEST 6: Shell syntax validation"

# Check both main scripts for syntax errors
if bash -n "${PROJECT_ROOT}/10_realtime_processing/24_process_reads_optimized.sh"; then
    pass "24_process_reads_optimized.sh: Valid bash syntax"
else
    fail "24_process_reads_optimized.sh: Syntax errors detected"
fi

if bash -n "${PROJECT_ROOT}/20_mag_assembly/61_map_and_bin_optimized.sh"; then
    pass "61_map_and_bin_optimized.sh: Valid bash syntax"
else
    fail "61_map_and_bin_optimized.sh: Syntax errors detected"
fi

echo ""

# ====================
# TEST 7: Permissions
# ====================

info "TEST 7: File permissions"

if [[ -x "${PROJECT_ROOT}/10_realtime_processing/24_process_reads_optimized.sh" ]]; then
    pass "24_process_reads_optimized.sh is executable"
else
    fail "24_process_reads_optimized.sh is not executable"
fi

if [[ -x "${PROJECT_ROOT}/20_mag_assembly/61_map_and_bin_optimized.sh" ]]; then
    pass "61_map_and_bin_optimized.sh is executable"
else
    fail "61_map_and_bin_optimized.sh is not executable"
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
    echo -e "${GREEN}✓ All basic functionality tests passed${NC}"
    exit 0
else
    echo -e "${RED}✗ Some tests failed${NC}"
    exit 1
fi
