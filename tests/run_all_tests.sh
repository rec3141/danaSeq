#!/usr/bin/env bash
#
# run_all_tests.sh - Master test runner for dānaSeq pipeline
#
# Runs all test suites and generates summary report:
#   1. test_security.sh - Security fix verification (CRITICAL)
#   2. test_basic_pipeline.sh - Basic functionality
#   3. test_error_handling.sh - Error detection and recovery
#   4. test_resume.sh - Resume/checkpoint logic
#   5. test_edge_cases.sh - Edge cases and corner conditions
#
# Usage: ./run_all_tests.sh [--verbose] [--stop-on-fail]
#
# Options:
#   --verbose       Show detailed test output
#   --stop-on-fail  Stop immediately if any test suite fails
#   --help          Show this help message
#
# Exit: 0 if all tests pass, 1 if any fail

set -euo pipefail
IFS=$'\n\t'

# Colors for output
GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
BOLD='\033[1m'
NC='\033[0m'

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Options
VERBOSE=0
STOP_ON_FAIL=0

# Parse arguments
while [[ $# -gt 0 ]]; do
    case "$1" in
        -v|--verbose)
            VERBOSE=1
            shift
            ;;
        -s|--stop-on-fail)
            STOP_ON_FAIL=1
            shift
            ;;
        -h|--help)
            grep '^#' "$0" | sed 's/^# \?//'
            exit 0
            ;;
        *)
            echo "Unknown option: $1" >&2
            echo "Use --help for usage information" >&2
            exit 1
            ;;
    esac
done

# Test suite definitions
declare -a TEST_SUITES=(
    "test_security.sh|Security Fixes|CRITICAL"
    "test_basic_pipeline.sh|Basic Functionality|HIGH"
    "test_error_handling.sh|Error Handling|HIGH"
    "test_resume.sh|Resume Logic|MEDIUM"
    "test_edge_cases.sh|Edge Cases|MEDIUM"
)

# Results tracking
SUITES_TOTAL=${#TEST_SUITES[@]}
SUITES_PASSED=0
SUITES_FAILED=0
SUITES_SKIPPED=0

declare -a FAILED_SUITES=()
declare -a SKIPPED_SUITES=()

# Helper functions
print_header() {
    echo ""
    echo -e "${BLUE}${BOLD}================================================================${NC}"
    echo -e "${BLUE}${BOLD}  $1${NC}"
    echo -e "${BLUE}${BOLD}================================================================${NC}"
    echo ""
}

print_suite_header() {
    local suite_name="$1"
    local priority="$2"

    echo ""
    echo -e "${YELLOW}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
    echo -e "${YELLOW}Running: ${BOLD}${suite_name}${NC} ${YELLOW}[Priority: ${priority}]${NC}"
    echo -e "${YELLOW}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
}

run_test_suite() {
    local script="$1"
    local name="$2"
    local priority="$3"

    local script_path="${SCRIPT_DIR}/${script}"

    # Check if test script exists
    if [[ ! -f "$script_path" ]]; then
        echo -e "${RED}[SKIP]${NC} Test script not found: ${script}"
        SKIPPED_SUITES+=("$name")
        ((SUITES_SKIPPED++))
        return 1
    fi

    # Make executable if needed
    chmod +x "$script_path" 2>/dev/null || true

    print_suite_header "$name" "$priority"

    # Run test
    local output
    local exit_code=0

    if (( VERBOSE )); then
        # Show full output
        if ! "$script_path"; then
            exit_code=$?
        fi
    else
        # Capture output, show summary only
        if ! output=$("$script_path" 2>&1); then
            exit_code=$?
        fi

        # Show last 10 lines (summary)
        if [[ -n "${output:-}" ]]; then
            echo "$output" | tail -20
        fi
    fi

    echo ""

    # Handle result
    if (( exit_code == 0 )); then
        echo -e "${GREEN}✓ ${BOLD}${name}${NC} ${GREEN}PASSED${NC}"
        ((SUITES_PASSED++))
        return 0
    else
        echo -e "${RED}✗ ${BOLD}${name}${NC} ${RED}FAILED${NC}"
        FAILED_SUITES+=("$name")
        ((SUITES_FAILED++))

        if (( STOP_ON_FAIL )); then
            echo ""
            echo -e "${RED}Stopping due to --stop-on-fail${NC}"
            return 1
        fi
        return 1
    fi
}

# Cleanup on exit
cleanup() {
    # Clean up any test artifacts
    rm -rf "${SCRIPT_DIR}/test_output" 2>/dev/null || true
    rm -rf "${SCRIPT_DIR}/test_fixtures/FC001" 2>/dev/null || true
}

trap cleanup EXIT

# Main execution
print_header "dānaSeq Pipeline - Test Suite Runner"

echo "Test Configuration:"
echo "  • Verbose output: $( (( VERBOSE )) && echo "ENABLED" || echo "DISABLED" )"
echo "  • Stop on failure: $( (( STOP_ON_FAIL )) && echo "ENABLED" || echo "DISABLED" )"
echo "  • Test suites: ${SUITES_TOTAL}"
echo ""

# Check if we're in the right directory
if [[ ! -f "${SCRIPT_DIR}/../nanopore_live/24_process_reads_optimized.sh" ]]; then
    echo -e "${RED}[ERROR]${NC} Cannot find pipeline scripts. Run from tests/ directory." >&2
    exit 1
fi

# Run each test suite
for suite in "${TEST_SUITES[@]}"; do
    IFS='|' read -r script name priority <<< "$suite"

    if ! run_test_suite "$script" "$name" "$priority"; then
        if (( STOP_ON_FAIL )); then
            break
        fi
    fi
done

# Generate final report
print_header "Test Results Summary"

echo -e "${BOLD}Results:${NC}"
echo -e "  ${GREEN}✓ Passed:  ${SUITES_PASSED}/${SUITES_TOTAL}${NC}"
echo -e "  ${RED}✗ Failed:  ${SUITES_FAILED}/${SUITES_TOTAL}${NC}"
echo -e "  ${YELLOW}⊘ Skipped: ${SUITES_SKIPPED}/${SUITES_TOTAL}${NC}"
echo ""

# Show failed suites if any
if (( SUITES_FAILED > 0 )); then
    echo -e "${RED}${BOLD}Failed Test Suites:${NC}"
    for suite in "${FAILED_SUITES[@]}"; do
        echo -e "  ${RED}✗${NC} $suite"
    done
    echo ""
fi

# Show skipped suites if any
if (( SUITES_SKIPPED > 0 )); then
    echo -e "${YELLOW}${BOLD}Skipped Test Suites:${NC}"
    for suite in "${SKIPPED_SUITES[@]}"; do
        echo -e "  ${YELLOW}⊘${NC} $suite"
    done
    echo ""
fi

# Calculate pass rate
if (( SUITES_TOTAL > 0 )); then
    PASS_RATE=$(( (SUITES_PASSED * 100) / SUITES_TOTAL ))
    echo -e "${BOLD}Pass Rate: ${PASS_RATE}%${NC}"
    echo ""
fi

# Final verdict
if (( SUITES_FAILED == 0 && SUITES_SKIPPED == 0 )); then
    echo -e "${GREEN}${BOLD}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
    echo -e "${GREEN}${BOLD}  ✓ ALL TESTS PASSED - Pipeline Ready for Production${NC}"
    echo -e "${GREEN}${BOLD}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
    echo ""
    exit 0
elif (( SUITES_FAILED == 0 )); then
    echo -e "${YELLOW}${BOLD}⚠ WARNING: Some tests were skipped${NC}"
    echo ""
    exit 0
else
    echo -e "${RED}${BOLD}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
    echo -e "${RED}${BOLD}  ✗ TESTS FAILED - Review failures before deployment${NC}"
    echo -e "${RED}${BOLD}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
    echo ""
    echo "Recommendations:"
    echo "  1. Review failed test output above"
    echo "  2. Run individual tests with: ./test_<name>.sh"
    echo "  3. Use --verbose flag for detailed output"
    echo "  4. Check SECURITY_AUDIT.md for remediation guidance"
    echo ""
    exit 1
fi
