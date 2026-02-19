# dānaSeq Pipeline Test Suite

Automated test suite for validating the dānaSeq metagenomics pipeline security fixes and functionality.

## Quick Start

```bash
# Run all tests
cd /data/danav2/tests
./run_all_tests.sh

# Run with verbose output
./run_all_tests.sh --verbose

# Stop on first failure
./run_all_tests.sh --stop-on-fail
```

## Test Suites

### 1. Security Fixes (`test_security.sh`) - CRITICAL

Verifies all 95 security vulnerabilities fixed on 2026-02-06:
- ✅ **CRITICAL (4):** Race conditions, command injection, path traversal, signal handlers
- ✅ **HIGH (17):** Loop safety, mkdir error checks, pipeline validation
- ✅ **MEDIUM (30):** Tool validation, input checks, output verification
- ✅ **LOW (44):** Message standardization, documentation

**Run individually:**
```bash
./test_security.sh
```

**What it checks:**
- Atomic lock acquisition (race condition fix)
- No `eval` command (command injection fix)
- Path normalization with `realpath` (path traversal fix)
- Signal handlers (EXIT, INT, TERM, HUP)
- Loop variable quoting
- Error checking on critical operations
- Standardized error message format

### 2. Basic Functionality (`test_basic_pipeline.sh`) - HIGH

Tests core pipeline operations:
- Script accepts valid arguments
- Creates expected directory structure
- Help/usage output works
- Detects missing dependencies gracefully
- Validates syntax

**Run individually:**
```bash
./test_basic_pipeline.sh
```

### 3. Error Handling (`test_error_handling.sh`) - HIGH

Verifies error detection and recovery:
- Invalid inputs rejected with clear messages
- Pipeline failures caught and reported
- Cleanup occurs on error/interrupt
- Exit codes are correct
- Error messages follow standardized format

**Run individually:**
```bash
./test_error_handling.sh
```

### 4. Resume Logic (`test_resume.sh`) - MEDIUM

Tests checkpoint and resume capability:
- Detects existing output files
- Skips completed work
- Handles partial/incomplete files
- Atomic file operations prevent corruption
- Lock files prevent race conditions

**Run individually:**
```bash
./test_resume.sh
```

### 5. Edge Cases (`test_edge_cases.sh`) - MEDIUM

Tests unusual but valid scenarios:
- Empty directories
- Special characters in filenames
- Very long paths
- Symbolic links
- Permission issues
- Corrupted FASTQ files
- Concurrent access protection

**Run individually:**
```bash
./test_edge_cases.sh
```

## Test Results Interpretation

### Exit Codes

- **0** = All tests passed
- **1** = One or more tests failed

### Output Format

Each test suite reports:
- `[PASS]` - Test passed ✓
- `[FAIL]` - Test failed ✗
- `[INFO]` - Informational message

### Example Output

```
======================================
Security Fix Verification Test Suite
======================================

[INFO] Testing CRITICAL security fixes...

[PASS] Atomic mkdir lock: Found 'if ! mkdir "$lockfile"'
[PASS] Lock cleanup on function return: Found 'trap "rmdir.*lockfile.*" RETURN'
[PASS] Lock wait loop: Found 'while [[ -d "$lockfile" ]]'

...

======================================
Test Results
======================================
PASSED: 42
FAILED: 0

✓ All security fixes verified
```

## Continuous Integration

### Running Tests in CI/CD

```bash
#!/bin/bash
# ci-test.sh

set -euo pipefail

cd /data/danav2/tests

# Run all tests with verbose output
./run_all_tests.sh --verbose --stop-on-fail

# Check exit code
if [[ $? -eq 0 ]]; then
    echo "✓ All tests passed - deployment approved"
    exit 0
else
    echo "✗ Tests failed - blocking deployment"
    exit 1
fi
```

### Pre-Commit Hook

Add to `.git/hooks/pre-commit`:

```bash
#!/bin/bash
# Run security tests before committing

echo "Running security tests..."
/data/danav2/tests/test_security.sh

if [[ $? -ne 0 ]]; then
    echo "❌ Security tests failed - commit blocked"
    exit 1
fi

echo "✅ Security tests passed"
```

## Test Fixtures

The `test_fixtures/` directory contains minimal test data:
- Small FASTQ files (1000 reads)
- Mock directory structures
- Example configuration files

**Note:** Test fixtures are generated on-demand by test scripts to minimize repository size.

## Troubleshooting

### All Tests Fail Immediately

**Problem:** Scripts not executable

**Solution:**
```bash
chmod +x *.sh
```

### test_security.sh Fails

**Problem:** Security fixes not applied

**Solution:** Review the relevant patch documentation:
- `SECURITY_PATCH_2026-02-06.md` (CRITICAL fixes)
- `HIGH_PRIORITY_FIXES_2026-02-06.md`
- `MEDIUM_PRIORITY_FIXES_2026-02-06.md`
- `LOW_PRIORITY_FIXES_2026-02-06.md`

### test_basic_pipeline.sh Fails

**Problem:** Missing dependencies or incorrect paths

**Solution:**
```bash
# Check dependencies
./status.sh

# Verify paths in scripts
grep -n "^[A-Z_]*=" ../nanopore_live/24_process_reads_optimized.sh
```

### test_edge_cases.sh Has Many Failures

**Problem:** This is expected if testing on unusual filesystem

**Solution:** Edge case tests are informational. Most failures are acceptable if:
1. Security tests pass
2. Basic functionality tests pass
3. Error handling tests pass

## Development

### Adding New Tests

1. Create new test script: `test_<name>.sh`
2. Follow the template structure:
   ```bash
   #!/usr/bin/env bash
   set -euo pipefail
   IFS=$'\n\t'

   # Test tracking
   TESTS_PASSED=0
   TESTS_FAILED=0

   pass() { echo "[PASS] $*"; ((TESTS_PASSED++)); }
   fail() { echo "[FAIL] $*"; ((TESTS_FAILED++)); }

   # Your tests here

   # Exit with appropriate code
   (( TESTS_FAILED == 0 )) && exit 0 || exit 1
   ```

3. Add to `run_all_tests.sh`:
   ```bash
   declare -a TEST_SUITES=(
       ...
       "test_<name>.sh|Description|PRIORITY"
   )
   ```

4. Make executable:
   ```bash
   chmod +x test_<name>.sh
   ```

### Test Coverage Goals

- **Security:** 100% of security fixes verified
- **Functionality:** 80% of core features tested
- **Error Handling:** All error paths covered
- **Edge Cases:** Common edge cases handled

## Integration with Deployment

Before deploying to production:

```bash
# 1. Run full test suite
./run_all_tests.sh --verbose > test_results.log 2>&1

# 2. Review results
less test_results.log

# 3. If all pass, proceed with deployment
if [[ $? -eq 0 ]]; then
    echo "✅ Tests passed - ready for deployment"
    # Continue with deployment steps
else
    echo "❌ Tests failed - review test_results.log"
    exit 1
fi
```

## Maintenance

### Regular Testing Schedule

- **After code changes:** Run `test_security.sh` and `test_error_handling.sh`
- **Before commits:** Run `test_security.sh`
- **Before deployment:** Run `run_all_tests.sh --verbose`
- **Monthly:** Run full test suite with fresh data

### Updating Tests

When adding new security fixes:
1. Add test case to `test_security.sh`
2. Document in test comments
3. Update this README
4. Run full test suite to verify no regressions

## References

- **Security Audit:** `../SECURITY_AUDIT.md`
- **Security Patches:** `../SECURITY_PATCH_2026-02-06.md` and related files
- **Deployment Guide:** `../DEPLOYMENT_GUIDE.md`
- **Methods Documentation:** `../METHODS.md`

## Support

For issues with tests:
1. Check test output for specific failure messages
2. Review relevant documentation (listed above)
3. Run individual test with verbose output: `bash -x test_<name>.sh`
4. Check system logs if testing infrastructure issues

---

**Test Suite Version:** 1.0
**Created:** 2026-02-06
**Last Updated:** 2026-02-06
**Maintained By:** dānaSeq Development Team
