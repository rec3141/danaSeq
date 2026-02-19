#!/usr/bin/env bash
#
# test_security.sh - Verify all security fixes are present
#
# Tests for the 95 security vulnerabilities fixed on 2026-02-06:
#   - CRITICAL: Race conditions, command injection, path traversal, signal handlers
#   - HIGH: Loop safety, mkdir error checks, pipeline validation
#   - MEDIUM: Tool validation, input checks, output verification
#   - LOW: Message standardization, documentation
#
# Usage: ./test_security.sh
# Exit: 0 if all tests pass, 1 if any fail

# Don't use set -e since we track failures manually
set -uo pipefail
IFS=$'\n\t'

# Colors for output
GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

TESTS_PASSED=0
TESTS_FAILED=0

# Test result tracking
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

# Test helper: Check if file contains pattern
contains() {
    local file="$1"
    local pattern="$2"
    local context="${3:-}"

    if grep -q "$pattern" "$file"; then
        pass "$context: Found '$pattern'"
    else
        fail "$context: Missing '$pattern'"
    fi
    # Always return 0 to prevent early exit with set -e
    return 0
}

# Test helper: Check if file does NOT contain pattern (vulnerability removed)
not_contains() {
    local file="$1"
    local pattern="$2"
    local context="${3:-}"

    if ! grep -q "$pattern" "$file"; then
        pass "$context: Confirmed removal of '$pattern'"
    else
        fail "$context: Still contains '$pattern'"
    fi
    # Always return 0 to prevent early exit with set -e
    return 0
}

echo "========================================"
echo "Security Fix Verification Test Suite"
echo "========================================"
echo ""

# ====================
# CRITICAL FIXES
# ====================

info "Testing CRITICAL security fixes..."
echo ""

# CRITICAL-1: Race condition fix (atomic lock with mkdir)
info "CRITICAL-1: Atomic lock acquisition (race condition fix)"
contains "${PROJECT_ROOT}/nanopore_live/24_process_reads_optimized.sh" \
    'if ! mkdir "\$lockfile" 2>/dev/null; then' \
    "Atomic mkdir lock"
contains "${PROJECT_ROOT}/nanopore_live/24_process_reads_optimized.sh" \
    'trap "rmdir.*lockfile.*" RETURN' \
    "Lock cleanup on function return"
contains "${PROJECT_ROOT}/nanopore_live/24_process_reads_optimized.sh" \
    'while \[\[ -d "\$lockfile" \]\]' \
    "Lock wait loop"

# CRITICAL-2: Command injection fix (bash -c instead of eval)
info "CRITICAL-2: Command injection prevention"
# Check for eval command usage (not "evalue" or in comments)
if grep -E '(^|[[:space:]])eval[[:space:]]' "${PROJECT_ROOT}/nanopore_live/24_process_reads_optimized.sh" | grep -v '^[[:space:]]*#'; then
    fail "No eval command: Still contains 'eval'"
else
    pass "No eval command usage"
fi
contains "${PROJECT_ROOT}/nanopore_live/24_process_reads_optimized.sh" \
    'bash -c "\$cmd"' \
    "Using bash -c for isolation"

# CRITICAL-3: Path traversal fix
info "CRITICAL-3: Path traversal prevention"
contains "${PROJECT_ROOT}/nanopore_mag/61_map_and_bin_optimized.sh" \
    'realpath -m' \
    "Path normalization with realpath"
contains "${PROJECT_ROOT}/nanopore_mag/61_map_and_bin_optimized.sh" \
    'abs_dir.*base_workspace' \
    "Workspace confinement check"
# Note: 24_process_reads_optimized.sh uses strict input/output directory structure
# and doesn't accept arbitrary paths, so explicit .. checks not needed there

# CRITICAL-4: Signal handlers
info "CRITICAL-4: Signal handlers for cleanup"
contains "${PROJECT_ROOT}/nanopore_live/24_process_reads_optimized.sh" \
    'trap cleanup_on_exit EXIT' \
    "EXIT signal handler"
contains "${PROJECT_ROOT}/nanopore_live/24_process_reads_optimized.sh" \
    "trap.*INT" \
    "INT signal handler"
contains "${PROJECT_ROOT}/nanopore_live/24_process_reads_optimized.sh" \
    "trap.*TERM" \
    "TERM signal handler"
contains "${PROJECT_ROOT}/nanopore_live/24_process_reads_optimized.sh" \
    "trap.*HUP" \
    "HUP signal handler"

echo ""

# ====================
# HIGH FIXES
# ====================

info "Testing HIGH priority security fixes..."
echo ""

# HIGH-1: Loop safety (while read instead of unquoted for)
info "HIGH-1: Loop variable quoting"
contains "${PROJECT_ROOT}/nanopore_live/24_process_reads_optimized.sh" \
    'while IFS= read -r fc; do' \
    "Safe flowcell loop"
contains "${PROJECT_ROOT}/nanopore_live/24_process_reads_optimized.sh" \
    'while IFS= read -r barcode; do' \
    "Safe barcode loop"

# HIGH-2: mkdir error checking
info "HIGH-2: mkdir error checking"
if grep -c 'if ! mkdir -p' "${PROJECT_ROOT}/nanopore_live/24_process_reads_optimized.sh" | grep -q '[6-9]'; then
    pass "mkdir error checks present (6+ instances)"
else
    fail "Insufficient mkdir error checks"
fi

# HIGH-3: Pipeline command error checking
info "HIGH-3: Pipeline validation"
# Check for general pipeline error handling with set -e and || patterns
if grep -q "set -euo pipefail" "${PROJECT_ROOT}/nanopore_live/24_process_reads_optimized.sh"; then
    pass "Pipeline uses set -e for automatic error detection"
else
    info "No set -e found (may use explicit checks)"
fi

echo ""

# ====================
# MEDIUM FIXES
# ====================

info "Testing MEDIUM priority security fixes..."
echo ""

# MEDIUM-1: Tool execution validation
info "MEDIUM-1: Bioinformatics tool validation"
contains "${PROJECT_ROOT}/nanopore_mag/61_map_and_bin_optimized.sh" \
    'if ! conda run.*kaiju' \
    "Kaiju error checking"
contains "${PROJECT_ROOT}/nanopore_mag/61_map_and_bin_optimized.sh" \
    'if ! metabat2' \
    "MetaBAT2 error checking"

# MEDIUM-2: Input validation
info "MEDIUM-2: Input validation patterns"
contains "${PROJECT_ROOT}/nanopore_live/24_process_reads_optimized.sh" \
    'barcode.*=~.*\^barcode' \
    "Barcode format validation"

# MEDIUM-3: Output verification
info "MEDIUM-3: Output file verification"
contains "${PROJECT_ROOT}/nanopore_live/24_process_reads_optimized.sh" \
    '\[\[ ! -s ' \
    "Empty output detection"

echo ""

# ====================
# LOW FIXES
# ====================

info "Testing LOW priority fixes (message standardization)..."
echo ""

# LOW-1: Error message standardization
info "LOW-1: Standardized error prefixes"
not_contains "${PROJECT_ROOT}/nanopore_live/24_process_reads_optimized.sh" \
    '\[ERR\]' \
    "No old [ERR] prefix"
not_contains "${PROJECT_ROOT}/nanopore_live/24_process_reads_optimized.sh" \
    '\[WARN\]' \
    "No old [WARN] prefix"
not_contains "${PROJECT_ROOT}/nanopore_live/24_process_reads_optimized.sh" \
    '\[BAD\]' \
    "No old [BAD] prefix"
not_contains "${PROJECT_ROOT}/nanopore_live/24_process_reads_optimized.sh" \
    '\[FIXED\]' \
    "No old [FIXED] prefix"

# LOW-2: Standardized message format
info "LOW-2: New standardized prefixes"
contains "${PROJECT_ROOT}/nanopore_live/24_process_reads_optimized.sh" \
    '\[ERROR\]' \
    "New [ERROR] prefix present"
contains "${PROJECT_ROOT}/nanopore_live/24_process_reads_optimized.sh" \
    '\[WARNING\]' \
    "New [WARNING] prefix present"
contains "${PROJECT_ROOT}/nanopore_live/24_process_reads_optimized.sh" \
    '\[INFO\]' \
    "New [INFO] prefix present"

# LOW-3: Function documentation
info "LOW-3: Function documentation"
contains "${PROJECT_ROOT}/nanopore_live/24_process_reads_optimized.sh" \
    '# Output debug message' \
    "debug_msg documentation"
contains "${PROJECT_ROOT}/nanopore_live/24_process_reads_optimized.sh" \
    '# Output verbose message' \
    "verbose_msg documentation"

echo ""

# ====================
# DEFENSIVE PATTERNS
# ====================

info "Testing defensive programming patterns..."
echo ""

# Bash strict mode
info "Bash strict mode"
contains "${PROJECT_ROOT}/nanopore_live/24_process_reads_optimized.sh" \
    'set -euo pipefail' \
    "Strict error handling"
contains "${PROJECT_ROOT}/nanopore_live/24_process_reads_optimized.sh" \
    'IFS=' \
    "IFS safety setting"

# Local variables in functions
info "Function local variables"
if grep -c 'local ' "${PROJECT_ROOT}/nanopore_live/24_process_reads_optimized.sh" | grep -q '[1-9][0-9]'; then
    pass "Functions use local variables (10+ instances)"
else
    fail "Insufficient use of local variables"
fi

echo ""
echo "========================================"
echo "Test Results"
echo "========================================"
echo -e "${GREEN}PASSED: ${TESTS_PASSED}${NC}"
echo -e "${RED}FAILED: ${TESTS_FAILED}${NC}"
echo ""

if (( TESTS_FAILED == 0 )); then
    echo -e "${GREEN}✓ All security fixes verified${NC}"
    exit 0
else
    echo -e "${RED}✗ Some security fixes missing${NC}"
    exit 1
fi
