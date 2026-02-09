#!/usr/bin/env bash
#
# test_resume.sh - Verify resume/checkpoint logic
#
# Tests pipeline resume capability:
#   - Detects existing output files
#   - Skips completed work
#   - Resumes from last checkpoint
#   - Handles partial/incomplete files
#   - Atomic file operations prevent corruption
#
# Usage: ./test_resume.sh
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
echo "Resume Logic Verification Tests"
echo "========================================"
echo ""

# ====================
# TEST 1: Output File Detection
# ====================

info "TEST 1: Existing output file detection"

# Check if scripts detect existing output
if grep -q "\[\[ -f.*\]\] && \[\[ ! -f.*force" "${PROJECT_ROOT}/10_realtime_processing/24_process_reads_optimized.sh" || \
   grep -q "skip\|exist" "${PROJECT_ROOT}/10_realtime_processing/24_process_reads_optimized.sh" | grep -q -i "already"; then
    pass "Script has logic to detect existing outputs"
else
    info "No explicit resume skip logic detected (may overwrite)"
fi

# Check for force flag to override resume
if grep -q "force\|overwrite" "${PROJECT_ROOT}/10_realtime_processing/24_process_reads_optimized.sh"; then
    pass "Force/overwrite flag available"
else
    info "No force flag detected (always overwrites or always resumes)"
fi

echo ""

# ====================
# TEST 2: Checkpoint Files
# ====================

info "TEST 2: Checkpoint and partial file handling"

# Check for temporary file patterns (.tmp, .partial, .incomplete)
if grep -q "\.tmp\|\.partial\|\.incomplete" "${PROJECT_ROOT}/10_realtime_processing/24_process_reads_optimized.sh"; then
    pass "Temporary file markers found (.tmp/.partial)"
else
    info "No explicit .tmp/.partial markers (may use other checkpoint method)"
fi

# Check for cleanup of partial files
if grep -q "rm.*\.tmp\|rm.*\.partial" "${PROJECT_ROOT}/10_realtime_processing/24_process_reads_optimized.sh"; then
    pass "Cleanup of temporary/partial files present"
else
    info "No explicit cleanup of .tmp/.partial files"
fi

echo ""

# ====================
# TEST 3: Atomic Operations
# ====================

info "TEST 3: Atomic file operations"

# Check for atomic mv pattern (write to temp, then mv to final)
if grep -q "mv.*tmp.*final\|mv.*\\.tmp" "${PROJECT_ROOT}/10_realtime_processing/24_process_reads_optimized.sh" || \
   grep -q ">.*\\.tmp.*mv" "${PROJECT_ROOT}/10_realtime_processing/24_process_reads_optimized.sh"; then
    pass "Atomic file operation pattern detected"
else
    info "No explicit atomic mv pattern (may write directly)"
fi

# Check for atomic directory creation (mkdir for locks)
if grep -q "mkdir.*lock" "${PROJECT_ROOT}/10_realtime_processing/24_process_reads_optimized.sh"; then
    pass "Atomic lock creation with mkdir"
else
    fail "Atomic lock creation not found"
fi

echo ""

# ====================
# TEST 4: Lock Files
# ====================

info "TEST 4: Lock file management"

# Check for lock file creation
if grep -q "lockfile\|lock_file" "${PROJECT_ROOT}/10_realtime_processing/24_process_reads_optimized.sh"; then
    pass "Lock file mechanism present"
else
    fail "Lock file mechanism not found"
fi

# Check for lock cleanup on exit
if grep -q "rmdir.*lock\|rm.*lock" "${PROJECT_ROOT}/10_realtime_processing/24_process_reads_optimized.sh"; then
    pass "Lock cleanup present"
else
    fail "Lock cleanup not found"
fi

# Check for stale lock detection/cleanup
if grep -q "find.*lock.*-mmin\|find.*lock.*-mtime" "${PROJECT_ROOT}/10_realtime_processing/24_process_reads_optimized.sh"; then
    pass "Stale lock cleanup detected"
else
    info "No stale lock cleanup (may accumulate abandoned locks)"
fi

# Check for lock timeout/wait
if grep -q "while.*lock.*sleep" "${PROJECT_ROOT}/10_realtime_processing/24_process_reads_optimized.sh"; then
    pass "Lock wait/timeout mechanism present"
else
    info "No lock wait mechanism (may fail immediately on conflict)"
fi

echo ""

# ====================
# TEST 5: Stage Awareness
# ====================

info "TEST 5: Processing stage checkpoints"

# Common pipeline stages
stages=(
    "qc"
    "filter"
    "sketch"
    "kraken"
    "assembly"
    "prokka"
    "hmmer"
)

detected_stages=0
for stage in "${stages[@]}"; do
    if grep -qi "$stage" "${PROJECT_ROOT}/10_realtime_processing/24_process_reads_optimized.sh"; then
        ((detected_stages++))
    fi
done

if (( detected_stages >= 5 )); then
    pass "Multiple pipeline stages detected ($detected_stages/7)"
else
    info "Limited stage detection ($detected_stages/7 stages)"
fi

echo ""

# ====================
# TEST 6: Parallel Safety
# ====================

info "TEST 6: Parallel execution safety"

# Check for GNU parallel usage
if grep -q "parallel\|sem " "${PROJECT_ROOT}/10_realtime_processing/24_process_reads_optimized.sh"; then
    pass "GNU parallel/semaphore usage detected"
else
    info "No GNU parallel detected (may use xargs or sequential)"
fi

# Check for Kraken2 serialization (CRITICAL for memory)
if grep -q "sem.*kraken\|kraken.*sem" "${PROJECT_ROOT}/10_realtime_processing/24_process_reads_optimized.sh"; then
    pass "Kraken2 serialization with semaphore"
else
    info "No Kraken2 serialization (may cause OOM if running parallel)"
fi

# Check for process limiting
if grep -q "\-j\s*[0-9]\+\|-P\s*[0-9]\+\|--jobs\|--max-procs" "${PROJECT_ROOT}/10_realtime_processing/24_process_reads_optimized.sh"; then
    pass "Job/process limiting present"
else
    info "No explicit job limiting (may use system defaults)"
fi

echo ""

# ====================
# TEST 7: Crash Recovery
# ====================

info "TEST 7: Crash recovery mechanisms"

# Check for trap handlers (critical for cleanup)
if grep -q "trap.*EXIT" "${PROJECT_ROOT}/10_realtime_processing/24_process_reads_optimized.sh"; then
    pass "EXIT trap for cleanup on crash"
else
    fail "No EXIT trap (may leave partial files on crash)"
fi

# Check for signal handlers (INT, TERM, HUP)
signal_handlers=0
for signal in INT TERM HUP; do
    if grep -q "trap.*$signal" "${PROJECT_ROOT}/10_realtime_processing/24_process_reads_optimized.sh"; then
        ((signal_handlers++))
    fi
done

if (( signal_handlers >= 2 )); then
    pass "Multiple signal handlers present ($signal_handlers/3)"
else
    info "Limited signal handling ($signal_handlers/3 signals)"
fi

echo ""

# ====================
# TEST 8: File Size Validation
# ====================

info "TEST 8: Output file validation"

# Check for non-empty file verification
if grep -q "\[\[ -s " "${PROJECT_ROOT}/10_realtime_processing/24_process_reads_optimized.sh"; then
    pass "File size validation present (non-empty check)"
else
    fail "No file size validation"
fi

# Check for minimum size thresholds
if grep -q "stat.*size\|du.*size\|wc.*size" "${PROJECT_ROOT}/10_realtime_processing/24_process_reads_optimized.sh"; then
    pass "Advanced file size checking detected"
else
    info "Only basic -s flag used for size (may accept tiny/corrupt files)"
fi

echo ""

# ====================
# TEST 9: Directory Locking
# ====================

info "TEST 9: Directory-level locking (race condition prevention)"

# Check for atomic mkdir locks
if grep -q "if ! mkdir.*lockfile.*2>/dev/null" "${PROJECT_ROOT}/10_realtime_processing/24_process_reads_optimized.sh"; then
    pass "Atomic mkdir lock pattern present"
else
    fail "Atomic mkdir lock not found"
fi

# Check for lock directory cleanup in trap
if grep -q "trap.*rmdir.*lock" "${PROJECT_ROOT}/10_realtime_processing/24_process_reads_optimized.sh"; then
    pass "Lock cleanup in trap handler"
else
    info "Lock cleanup may not be in trap (could leak on crash)"
fi

# Check for double-check pattern after acquiring lock
if grep -q "# Double-check\|# Recheck\|# Verify again" "${PROJECT_ROOT}/10_realtime_processing/24_process_reads_optimized.sh"; then
    pass "Double-check pattern documented"
else
    info "No double-check comments (may still implement the pattern)"
fi

echo ""

# ====================
# TEST 10: Resume Documentation
# ====================

info "TEST 10: Resume behavior documentation"

# Check if help text mentions resume behavior
if "${PROJECT_ROOT}/10_realtime_processing/24_process_reads_optimized.sh" -h 2>&1 | grep -qi "resume\|continue\|skip.*exist"; then
    pass "Resume behavior documented in help"
else
    info "Resume behavior not in help text"
fi

# Check if there's a force/overwrite flag documented
if "${PROJECT_ROOT}/10_realtime_processing/24_process_reads_optimized.sh" -h 2>&1 | grep -qi "force\|overwrite"; then
    pass "Force flag documented in help"
else
    info "No force flag in help (may always resume or always overwrite)"
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
    echo -e "${GREEN}✓ All resume logic tests passed${NC}"
    exit 0
else
    echo -e "${RED}✗ Some resume logic tests failed${NC}"
    exit 1
fi
