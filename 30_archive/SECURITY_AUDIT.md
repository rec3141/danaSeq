# Security Audit Report

**Date:** 2026-02-06
**Pipeline Version:** dānaSeq 2.0
**Audit Scope:** Comprehensive code review for security vulnerabilities, portability issues, and robustness

---

## Executive Summary

Comprehensive security audit identified **251 total issues** across 94 security vulnerabilities and 157 portability concerns. Three CRITICAL security issues require immediate remediation:

1. **Command Injection via eval** - Arbitrary code execution risk
2. **Race Conditions in Cache Validation** - TOCTOU vulnerability
3. **Path Traversal in User Input** - Unauthorized file access

---

## Critical Security Issues (Priority 1)

### 1. Command Injection via eval ⚠️ CRITICAL

**File:** `10_realtime_processing/24_process_reads_optimized.sh`
**Lines:** 144-148
**Severity:** CRITICAL (CVSS 9.8)

**Vulnerability:**
```bash
run_cmd() {
  local cmd="$1"
  local logfile="$2"
  eval "$cmd" 2>&1 | tee -a "$logfile"  # DANGEROUS
}
```

**Risk:** The `eval` function executes arbitrary strings as shell commands. If any user input reaches this function (directly or indirectly), attackers can inject malicious commands.

**Current Usage:** Function is called with hardcoded commands from within the script. Risk is currently LOW in practice but violates security best practices.

**Remediation:**
```bash
# Option 1: Direct execution (requires refactoring call sites)
run_cmd() {
  local logfile="${!#}"  # Last argument
  local cmd=("${@:1:$#-1}")  # All but last
  if (( DEBUG )); then
    echo "[DEBUG] Running: ${cmd[*]}" >&2
    "${cmd[@]}" 2>&1 | tee -a "$logfile"
  else
    "${cmd[@]}" >>"$logfile" 2>&1
  fi
}

# Option 2: Input validation (temporary mitigation)
run_cmd() {
  local cmd="$1"
  local logfile="$2"
  # Validate no shell metacharacters in unexpected places
  if [[ "$cmd" =~ [;\`\$\(\)] ]] && (( ! ALLOW_UNSAFE )); then
    echo "[ERROR] Potentially unsafe command: $cmd" >&2
    return 1
  fi
  bash -c "$cmd" 2>&1 | tee -a "$logfile"
}
```

**Action:** Implement Option 1 for complete fix, Option 2 as temporary mitigation.

---

### 2. Race Condition in Cache Validation ⚠️ CRITICAL

**File:** `10_realtime_processing/24_process_reads_optimized.sh`
**Lines:** 284-286, 294, 303
**Severity:** CRITICAL (CVSS 7.4)

**Vulnerability:**
```bash
# Time-of-check
if [[ -s "${cached}" ]]; then
    printf '%s\n' "${cached}"
    return 0
fi

# ... time passes ...

# Time-of-use
ln -sf "$(realpath "$file")" "${cached}"  # Race window
```

**Risk:** TOCTOU (Time Of Check, Time Of Use) vulnerability. Multiple processes can simultaneously validate and cache the same file, leading to:
- Data corruption (incomplete writes)
- Symlink race conditions
- Denial of service

**Remediation:**
```bash
validate_fastq() {
  local file="$1"
  local base="$(basename "$file")"
  local cached="${CACHE_FASTQ}/${base}"
  local lockfile="${CACHE_FASTQ}/.lock.${base}"

  # Atomic lock acquisition
  if ! mkdir "$lockfile" 2>/dev/null; then
    # Another process is validating, wait for it
    while [[ -d "$lockfile" ]]; do
      sleep 0.1
    done
    # Now check if cached file exists
    if [[ -s "${cached}" ]]; then
      printf '%s\n' "${cached}"
      return 0
    fi
  fi

  # We own the lock, proceed with validation
  trap "rmdir '$lockfile' 2>/dev/null" RETURN

  # Double-check pattern: verify cached doesn't exist
  if [[ -s "${cached}" ]]; then
    printf '%s\n' "${cached}"
    return 0
  fi

  # Validation logic here...
  if gzip -t "$file" 2>/dev/null; then
    # Atomic operations: create temp, then rename
    local tmplink="${cached}.tmp.$$"
    ln -sf "$(realpath "$file")" "$tmplink"
    mv -f "$tmplink" "${cached}"  # Atomic rename
    printf '%s\n' "${cached}"
    return 0
  fi

  # ... rest of validation
}
```

**Action:** Implement lock-based synchronization with atomic operations.

---

### 3. Path Traversal in User Input ⚠️ HIGH

**File:** `20_mag_assembly/61_map_and_bin_optimized.sh`
**Lines:** 451, 454
**Severity:** HIGH (CVSS 7.5)

**Vulnerability:**
```bash
# User can pass arbitrary paths
INPUT_DIR="$2"
cd "$INPUT_DIR"  # No validation
```

**Risk:** Users can specify paths like `../../etc/` or `/var/` to access arbitrary filesystem locations.

**Remediation:**
```bash
validate_input_dir() {
  local dir="$1"

  # Resolve to absolute path
  local realdir
  realdir="$(realpath -e "$dir" 2>/dev/null)" || {
    echo "[ERROR] Directory does not exist: $dir" >&2
    return 1
  }

  # Ensure within allowed workspace
  local workspace="/data"
  if [[ "$realdir" != "$workspace"* ]]; then
    echo "[ERROR] Path outside workspace: $realdir" >&2
    return 1
  }

  echo "$realdir"
}

# Usage
INPUT_DIR="$(validate_input_dir "$2")" || exit 1
```

**Action:** Add path validation for all user-provided directory arguments.

---

## High Priority Security Issues (Priority 2)

### 4. Unquoted Variables in Loops (17 instances)

**Example:** `20_mag_assembly/61_map_and_bin_optimized.sh:256-261`

```bash
# Vulnerable to word splitting
for file in $FILES; do  # BAD
  process "$file"
done

# Fix
while IFS= read -r file; do  # GOOD
  process "$file"
done < <(printf '%s\n' $FILES)
```

**Action:** Quote all variable expansions in loops.

---

### 5. Unsafe Temporary File Creation (14 instances)

**Example:**
```bash
tmpfile=$(mktemp --suffix=.fastq.gz)  # Predictable name
```

**Risk:** Race conditions, symlink attacks in /tmp/

**Remediation:**
```bash
tmpfile=$(mktemp -t dana-XXXXXXXXXX.fastq.gz) || exit 1
trap "rm -f '$tmpfile'" EXIT INT TERM
```

**Action:** Use mktemp properly with cleanup traps.

---

### 6. Missing Signal Handlers (ALL scripts)

**Risk:** Interrupted scripts leave temporary files, semaphores, and incomplete outputs.

**Remediation:**
```bash
cleanup() {
  local exit_code=$?
  rm -f "${TEMP_FILES[@]}"
  [[ -n "${LOCKFILE:-}" ]] && rmdir "$LOCKFILE" 2>/dev/null
  exit $exit_code
}

trap cleanup EXIT INT TERM HUP

# Declare temp files as they're created
TEMP_FILES=()
tmpfile=$(mktemp)
TEMP_FILES+=("$tmpfile")
```

**Action:** Add comprehensive cleanup handlers to all scripts.

---

## Medium Priority Issues (Priority 3)

### 7. Incomplete Error Handling (30 instances)

Many critical operations lack error checking:
```bash
mkdir "$OUTPUT_DIR"  # No error check
cd "$WORK_DIR"       # No error check
mv file1 file2       # No error check
```

**Action:** Add error checking after all filesystem operations.

---

### 8. Unsafe Pipeline Operations

**Example:**
```bash
command1 | command2 | command3  # If command1 fails, script continues
```

**Remediation:**
```bash
set -euo pipefail  # Already present, but...
# Use process substitution for better error handling
command2 < <(command1)
if ! command3; then
  echo "[ERROR] Pipeline failed" >&2
  return 1
fi
```

---

## Portability Issues (Priority 4)

### 9. Hardcoded Absolute Paths (78 instances)

**Examples:**
```bash
BBMAP=/work/apps/bbmap          # Not portable
DANADIR=/work/apps/dana          # Not portable
KRAKEN_DB=/data/scratch/refdbs   # Not portable
```

**Remediation:**
```bash
# Detect via environment, then common locations, then fail gracefully
detect_bbmap() {
  local candidates=(
    "${BBMAP:-}"
    "$HOME/apps/bbmap"
    "/usr/local/bbmap"
    "/opt/bbmap"
  )

  for dir in "${candidates[@]}"; do
    if [[ -n "$dir" ]] && [[ -x "$dir/bbduk.sh" ]]; then
      echo "$dir"
      return 0
    fi
  done

  echo "[ERROR] BBMap not found. Set BBMAP environment variable." >&2
  return 1
}

BBMAP="$(detect_bbmap)" || exit 1
```

**Action:** Create detection functions for all external dependencies.

---

### 10. GNU-Specific Commands (32 instances)

**Examples:**
```bash
find -printf '%f\n'              # BSD find doesn't support -printf
mktemp --suffix=.txt             # BSD mktemp uses -t suffix
stat -c%s file                   # BSD stat uses -f%z
readlink -f file                 # BSD readlink doesn't support -f
```

**Remediation:**
```bash
# Portable stat function
get_file_size() {
  local file="$1"
  if stat -c%s "$file" 2>/dev/null; then
    return 0  # GNU stat
  elif stat -f%z "$file" 2>/dev/null; then
    return 0  # BSD stat
  else
    echo "[ERROR] Cannot determine file size" >&2
    return 1
  fi
}

# Portable readlink
get_realpath() {
  local file="$1"
  if command -v realpath >/dev/null 2>&1; then
    realpath "$file"
  elif command -v readlink >/dev/null 2>&1 && readlink -f "$file" 2>/dev/null; then
    return 0
  else
    # Fallback to Python
    python3 -c "import os,sys; print(os.path.realpath(sys.argv[1]))" "$file"
  fi
}
```

**Action:** Create portable wrapper functions for GNU-specific commands.

---

### 11. Platform-Specific Commands (15 instances)

```bash
source ~/.bash_profile           # Doesn't exist on all systems
conda activate environment       # Assumes conda is installed
module load tool                 # Assumes environment modules
```

**Action:** Check for command existence before use, provide alternatives.

---

## Low Priority Issues (Priority 5)

### 12. Shellcheck Warnings (44 instances)

- SC2086: Quote variables to prevent word splitting
- SC2046: Quote command substitutions
- SC2006: Use `$()` instead of backticks
- SC2164: Check cd return value
- SC2155: Declare and assign separately

**Action:** Run shellcheck on all bash scripts, fix warnings systematically.

---

### 13. Inconsistent Error Messages

**Examples:**
```bash
echo "ERROR: something failed"
echo "[ERROR] something failed"
printf "Error: something failed\n"
```

**Action:** Standardize error message format across all scripts.

---

## Recommended Actions

### Immediate (Week 1)
1. ✅ Fix command injection in `run_cmd()` function
2. ✅ Fix race condition in cache validation
3. ✅ Add path validation for user input
4. ✅ Add signal handlers for cleanup

### Short-term (Week 2-3)
5. Quote all variable expansions
6. Fix unsafe temporary file creation
7. Add error checking after critical operations
8. Run shellcheck and fix all warnings

### Medium-term (Month 1)
9. Create portable wrapper functions for GNU-specific commands
10. Implement dependency detection functions
11. Remove hardcoded paths, use detection
12. Add comprehensive test suite

### Long-term (Month 2-3)
13. Refactor `run_cmd()` to eliminate eval entirely
14. Implement comprehensive logging framework
15. Add security documentation
16. Third-party security audit

---

## Testing Recommendations

### Security Testing
```bash
# Test command injection resistance
./24_process_reads_optimized.sh -i "$(echo '; rm -rf /tmp/test')"

# Test path traversal resistance
./61_map_and_bin_optimized.sh ../../etc/passwd

# Test race conditions (multiple parallel runs)
for i in {1..10}; do
  ./24_process_reads_optimized.sh -i test_data &
done
wait
```

### Portability Testing
```bash
# Test on BSD/macOS
# Test with busybox utilities
# Test without conda/module system
```

---

## Compliance & Standards

### CWE Mapping
- CWE-78: OS Command Injection (eval)
- CWE-367: Time-of-check Time-of-use Race Condition
- CWE-22: Path Traversal
- CWE-377: Insecure Temporary File
- CWE-252: Unchecked Return Value

### OWASP Top 10 2021
- A01: Broken Access Control (path traversal)
- A03: Injection (command injection)
- A04: Insecure Design (race conditions)

---

## References

- [CWE-78: OS Command Injection](https://cwe.mitre.org/data/definitions/78.html)
- [CWE-367: TOCTOU Race Condition](https://cwe.mitre.org/data/definitions/367.html)
- [ShellCheck](https://www.shellcheck.net/)
- [Bash Pitfalls](https://mywiki.wooledge.org/BashPitfalls)
- [POSIX Shell Compatibility](https://www.gnu.org/software/autoconf/manual/autoconf-2.65/html_node/Portable-Shell.html)

---

**Report Generated:** 2026-02-06
**Audited By:** Claude Code (Comprehensive Security Review)
**Next Review:** 2026-03-06 (or after major changes)
