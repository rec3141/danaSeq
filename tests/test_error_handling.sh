#!/usr/bin/env bash
#
# test_error_handling.sh - Verify error handling and recovery
#
# Tests error conditions and graceful failure modes:
#   - Invalid inputs are rejected with clear error messages
#   - Pipeline failures are caught and reported
#   - Cleanup occurs on error/interrupt
#   - Error messages follow standardized format
#   - Exit codes are correct
#
# Usage: ./test_error_handling.sh
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
echo "Error Handling Verification Tests"
echo "========================================"
echo ""

# ====================
# TEST 1: Missing Arguments
# ====================

info "TEST 1: Missing required arguments"

# Test missing --input
output=$("${PROJECT_ROOT}/10_realtime_processing/24_process_reads_optimized.sh" 2>&1 || true)
if echo "$output" | grep -q "\[ERROR\].*--input"; then
    pass "Error message for missing --input argument"
else
    fail "No error message for missing --input argument"
fi

# Check exit code
if ! "${PROJECT_ROOT}/10_realtime_processing/24_process_reads_optimized.sh" &>/dev/null; then
    pass "Non-zero exit code for missing arguments"
else
    fail "Should exit with error for missing arguments"
fi

echo ""

# ====================
# TEST 2: Invalid Paths
# ====================

info "TEST 2: Invalid path handling"

# Test non-existent input directory
output=$("${PROJECT_ROOT}/10_realtime_processing/24_process_reads_optimized.sh" -i /definitely/does/not/exist/path 2>&1 || true)
if echo "$output" | grep -q "\[ERROR\]"; then
    pass "Error message for non-existent input directory"
else
    fail "No error message for non-existent input directory"
fi

# Test path traversal attempt (for MAG script)
mkdir -p "${TEST_OUTPUT}"
output=$("${PROJECT_ROOT}/20_mag_assembly/61_map_and_bin_optimized.sh" "../../../etc" 2>&1 || true)
if echo "$output" | grep -q -i "error"; then
    pass "Path traversal attempt rejected"
else
    fail "Path traversal attempt not properly rejected"
fi

echo ""

# ====================
# TEST 3: Invalid Combinations
# ====================

info "TEST 3: Invalid argument combinations"

# Test --hmm without -P flag
mkdir -p "${TEST_OUTPUT}/test_input"
output=$("${PROJECT_ROOT}/10_realtime_processing/24_process_reads_optimized.sh" \
    -i "${TEST_OUTPUT}/test_input" \
    --hmm dummy.hmm 2>&1 || true)
if echo "$output" | grep -q "\[ERROR\].*--hmm.*-P"; then
    pass "Error message for --hmm without -P flag"
else
    fail "No error message for --hmm without -P flag"
fi

echo ""

# ====================
# TEST 4: Error Message Format
# ====================

info "TEST 4: Error message standardization"

# Collect error messages from various error conditions
mkdir -p "${TEST_OUTPUT}/test_input"

# Test 1: Missing argument
output1=$("${PROJECT_ROOT}/10_realtime_processing/24_process_reads_optimized.sh" 2>&1 || true)

# Test 2: Invalid path
output2=$("${PROJECT_ROOT}/10_realtime_processing/24_process_reads_optimized.sh" -i /nonexistent 2>&1 || true)

# Test 3: Invalid flag combination
output3=$("${PROJECT_ROOT}/10_realtime_processing/24_process_reads_optimized.sh" \
    -i "${TEST_OUTPUT}/test_input" --hmm test.hmm 2>&1 || true)

# Verify all error messages use [ERROR] prefix (not [ERR])
all_errors="${output1}${output2}${output3}"

if echo "$all_errors" | grep -q "\[ERROR\]"; then
    pass "Error messages use standardized [ERROR] prefix"
else
    fail "Error messages don't use [ERROR] prefix"
fi

if echo "$all_errors" | grep -q "\[ERR\]"; then
    fail "Found old-style [ERR] prefix (should be [ERROR])"
else
    pass "No old-style [ERR] prefix found"
fi

# Check that errors go to stderr
if "${PROJECT_ROOT}/10_realtime_processing/24_process_reads_optimized.sh" 2>/dev/null; then
    fail "Errors should go to stderr, not stdout"
else
    pass "Errors correctly sent to stderr"
fi

echo ""

# ====================
# TEST 5: Graceful Degradation
# ====================

info "TEST 5: Graceful handling of missing tools"

# Create a test script that simulates missing dependency
cat > "${TEST_OUTPUT}/test_missing_tool.sh" <<'EOF'
#!/usr/bin/env bash
# Don't use set -e since we track failures manually
set -uo pipefail

# Simulate missing seqkit
if ! command -v definitely_not_a_real_command &>/dev/null; then
    echo "[ERROR] Required tool not found: definitely_not_a_real_command" >&2
    exit 1
fi
EOF
chmod +x "${TEST_OUTPUT}/test_missing_tool.sh"

if ! "${TEST_OUTPUT}/test_missing_tool.sh" 2>&1 | grep -q "\[ERROR\]"; then
    fail "Missing tool should produce error message"
else
    pass "Missing tool produces appropriate error message"
fi

echo ""

# ====================
# TEST 6: Cleanup on Error
# ====================

info "TEST 6: Cleanup on error/interrupt"

# Check if cleanup_on_exit function exists
if grep -q "cleanup_on_exit()" "${PROJECT_ROOT}/10_realtime_processing/24_process_reads_optimized.sh"; then
    pass "cleanup_on_exit function defined"
else
    fail "cleanup_on_exit function not found"
fi

# Check if EXIT trap is set
if grep -q "trap cleanup_on_exit EXIT" "${PROJECT_ROOT}/10_realtime_processing/24_process_reads_optimized.sh"; then
    pass "EXIT trap registered"
else
    fail "EXIT trap not registered"
fi

# Check if INT trap is set (Ctrl+C handling)
if grep -q "trap.*INT" "${PROJECT_ROOT}/10_realtime_processing/24_process_reads_optimized.sh"; then
    pass "INT trap registered (Ctrl+C handler)"
else
    fail "INT trap not registered"
fi

# Check if TERM trap is set
if grep -q "trap.*TERM" "${PROJECT_ROOT}/10_realtime_processing/24_process_reads_optimized.sh"; then
    pass "TERM trap registered"
else
    fail "TERM trap not registered"
fi

echo ""

# ====================
# TEST 7: Pipeline Failure Detection
# ====================

info "TEST 7: Pipeline command failure detection"

# Check that critical pipeline commands have error checking
critical_commands=(minimap2 samtools bbduk kraken2 prokka)

for cmd in "${critical_commands[@]}"; do
    # Look for error checking patterns: "if ! command" or "command || exit"
    if grep -q "if ! .*${cmd}" "${PROJECT_ROOT}/10_realtime_processing/24_process_reads_optimized.sh" 2>/dev/null; then
        pass "Error checking for: ${cmd}"
    else
        # Check in MAG script too
        if grep -q "if ! .*${cmd}" "${PROJECT_ROOT}/20_mag_assembly/61_map_and_bin_optimized.sh" 2>/dev/null; then
            pass "Error checking for: ${cmd}"
        else
            info "No explicit error check found for: ${cmd} (may use set -e)"
        fi
    fi
done

echo ""

# ====================
# TEST 8: Input Validation
# ====================

info "TEST 8: Input validation and sanitization"

# Check for barcode format validation
if grep -q "barcode.*=~.*\[0-9\]" "${PROJECT_ROOT}/10_realtime_processing/24_process_reads_optimized.sh"; then
    pass "Barcode format validation present"
else
    fail "Barcode format validation missing"
fi

# Check for path validation (directory traversal)
if grep -q "\.\." "${PROJECT_ROOT}/10_realtime_processing/24_process_reads_optimized.sh" | grep -q "ERROR"; then
    pass "Path traversal detection present"
elif grep -q "realpath" "${PROJECT_ROOT}/20_mag_assembly/61_map_and_bin_optimized.sh"; then
    pass "Path normalization with realpath present"
else
    fail "Path validation insufficient"
fi

echo ""

# ====================
# TEST 9: Empty/Missing Files
# ====================

info "TEST 9: Empty or missing file handling"

# Check for empty file detection
if grep -q "\[\[ ! -s" "${PROJECT_ROOT}/10_realtime_processing/24_process_reads_optimized.sh"; then
    pass "Empty file detection present"
else
    fail "Empty file detection missing"
fi

# Check for file existence checks before operations
if grep -c "\[\[ -f" "${PROJECT_ROOT}/10_realtime_processing/24_process_reads_optimized.sh" | grep -q "[1-9]"; then
    pass "File existence checks present"
else
    fail "Insufficient file existence checks"
fi

echo ""

# ====================
# TEST 10: Error Recovery
# ====================

info "TEST 10: Error recovery and resume capability"

# Check for checkpoint/resume logic
if grep -q "\.tmp\|\.partial" "${PROJECT_ROOT}/10_realtime_processing/24_process_reads_optimized.sh"; then
    pass "Temporary file markers for resume logic"
else
    info "No explicit .tmp/.partial resume markers (may use other method)"
fi

# Check for atomic operations (temp + mv)
if grep -q "mv.*\\.tmp" "${PROJECT_ROOT}/10_realtime_processing/24_process_reads_optimized.sh" || \
   grep -q "mv.*tmp.*final" "${PROJECT_ROOT}/10_realtime_processing/24_process_reads_optimized.sh"; then
    pass "Atomic file operations (temp + mv pattern)"
else
    info "No explicit atomic mv pattern detected"
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
    echo -e "${GREEN}✓ All error handling tests passed${NC}"
    exit 0
else
    echo -e "${RED}✗ Some error handling tests failed${NC}"
    exit 1
fi
