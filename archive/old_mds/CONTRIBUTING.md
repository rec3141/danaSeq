# Contributing to dānaSeq

Thank you for your interest in contributing to dānaSeq. This document provides guidelines for contributing to the project.

---

## Getting Started

### Prerequisites

1. Familiarity with bash scripting and bioinformatics workflows
2. Understanding of metagenomic analysis principles
3. Access to test datasets (Oxford Nanopore sequencing data)
4. Development environment with required dependencies (see `./status.sh`)

### Setting Up Development Environment

```bash
# Clone repository
git clone https://github.com/rec3141/danaSeq.git
cd danaSeq

# Check dependencies
./status.sh

# Test with small dataset
cd nanopore_live
./24_process_reads_optimized.sh -i /path/to/test/data
```

---

## Development Guidelines

### Code Style

**Bash Scripts:**
- Use descriptive variable names in `UPPER_CASE`
- Include comprehensive header comments
- Add error checking after critical operations
- Use `set -euo pipefail` for robust error handling
- Prefer `[[` over `[` for conditionals

**R Scripts:**
- Follow tidyverse style guide
- Use descriptive variable names in `snake_case`
- Comment complex data transformations
- Handle missing data explicitly

**Python:**
- Follow PEP 8 style guidelines
- Use type hints for function signatures
- Document functions with docstrings
- Prefer pathlib over os.path

### File Naming Convention

Scripts must follow the numbered prefix system:
- Use increments of 10 (10, 20, 30...)
- Include descriptive names: `<number>_<action>_<object>.sh`
- Examples: `24_process_reads_optimized.sh`, `61_map_and_bin_optimized.sh`

### Safety Requirements

**Critical:** Use `trash` or `trash-put` instead of `rm` for all file operations:

```bash
# Good
trash old_file.txt
trash-put deprecated_directory/

# Bad
rm -rf directory/  # Never use this
```

Before any deletion operation:
1. Create a manifest: `ls -laR path/ > manifest_$(date +%Y%m%d_%H%M%S).txt`
2. Ask for user confirmation if interactive
3. Document in commit message what was removed and why

### Resume Logic

All long-running operations must support resume capability:

```bash
# Check for existing output
if [[ -f "$OUTPUT_FILE" ]] && (( ! FORCE )); then
    echo "Output exists, skipping. Use --force to override."
    return 0
fi

# Use atomic operations
command > "${OUTPUT_FILE}.tmp"
mv "${OUTPUT_FILE}.tmp" "$OUTPUT_FILE"
```

### Error Handling

Implement robust error handling:

```bash
set -euo pipefail

# Function-level error handling
process_file() {
    local input=$1
    local output=$2

    if [[ ! -f "$input" ]]; then
        echo "ERROR: Input file not found: $input" >&2
        return 1
    fi

    if ! tool --input "$input" --output "$output" 2>&1 | tee -a "$LOG"; then
        echo "ERROR: Processing failed for $input" >&2
        log_failure "$input" "tool" "see $LOG for details"
        return 1
    fi

    return 0
}
```

### Logging

Consistent logging format:

```bash
# Color-coded output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m'

log_info() {
    echo -e "${GREEN}[$(date +'%Y-%m-%d %H:%M:%S')]${NC} INFO: $*" | tee -a "$LOG"
}

log_error() {
    echo -e "${RED}[$(date +'%Y-%m-%d %H:%M:%S')]${NC} ERROR: $*" | tee -a "$LOG" >&2
}

log_warn() {
    echo -e "${YELLOW}[$(date +'%Y-%m-%d %H:%M:%S')]${NC} WARN: $*" | tee -a "$LOG"
}
```

---

## Testing

### Before Submitting

1. Test with small dataset (1000-10000 reads)
2. Verify resume logic works correctly
3. Check error messages are informative
4. Confirm resource usage is reasonable
5. Run on clean system to catch missing dependencies

### Test Data

Maintain test datasets in `/data/test/`:
```
/data/test/
├── small/          # 1K reads, fast testing
├── medium/         # 10K reads, integration testing
└── large/          # 100K reads, performance testing
```

---

## Documentation

### Code Comments

```bash
################################################################################
# Script: 35_new_analysis.sh
# Purpose: Brief description of what this script does
#
# Author: Your Name
# Date: YYYY-MM-DD
# Version: 1.0
#
# Usage: ./35_new_analysis.sh -i <input> -o <output> [options]
#
# Options:
#   -i    Input directory (required)
#   -o    Output directory (required)
#   -f    Force overwrite existing files
#   -h    Display this help message
#
# Dependencies:
#   - tool1 (v1.0+)
#   - tool2 (v2.0+)
#
# Notes:
#   - Implements resume logic for long-running operations
#   - Requires 32GB RAM minimum
################################################################################
```

### README Updates

When adding new features:
1. Update relevant README.md in subdirectory
2. Update main README.md if user-facing
3. Update METHODS.md for scientific methodology changes
4. Update CLAUDE.md for architectural changes

---

## Pull Request Process

### Before Submitting PR

1. **Test thoroughly** on representative dataset
2. **Update documentation** (README, METHODS, CLAUDE.md)
3. **Check for hardcoded paths** - document in DEPLOYMENT_ISSUES.md
4. **Verify no secrets** in code (API keys, passwords)
5. **Clean commit history** - squash WIP commits

### PR Description

Include:
```markdown
## Summary
Brief description of changes

## Motivation
Why is this change needed?

## Testing
- [ ] Tested with small dataset
- [ ] Tested resume logic
- [ ] Verified error handling
- [ ] Checked resource usage

## Documentation
- [ ] Updated README
- [ ] Updated METHODS (if applicable)
- [ ] Updated CLAUDE.md (if applicable)
- [ ] Added inline comments

## Breaking Changes
List any breaking changes and migration path
```

### Review Process

1. Automated checks (if implemented)
2. Manual code review
3. Testing on independent dataset
4. Documentation review
5. Merge after approval

---

## Reporting Issues

### Bug Reports

Include:
- **dānaSeq version** (git commit hash)
- **OS and version**
- **Tool versions** (output of `./status.sh`)
- **Input data description** (read count, quality, source)
- **Command executed** (exact command with arguments)
- **Error message** (full error output)
- **Log files** (relevant portions)
- **Expected vs actual behavior**

### Feature Requests

Include:
- **Use case description**
- **Proposed solution** (if any)
- **Alternative solutions considered**
- **Impact on existing workflow**
- **References** (papers, tools, methods)

---

## Communication

### Channels

- **Issues**: Bug reports and feature requests
- **Pull Requests**: Code contributions
- **Discussions**: General questions and ideas
- **Email**: rec3141@gmail.com for private matters

### Response Time

This is research software under active development. Response times may vary based on expedition schedules and academic commitments. Please be patient.

---

## Code of Conduct

### Our Standards

- **Be respectful** of differing viewpoints and experiences
- **Accept constructive criticism** gracefully
- **Focus on what is best** for the community and science
- **Show empathy** towards other community members

### Unacceptable Behavior

- Harassment, discrimination, or exclusionary behavior
- Publishing others' private information without permission
- Trolling, insulting comments, or personal attacks
- Other conduct inappropriate in a professional setting

---

## License

By contributing to dānaSeq, you agree that your contributions will be licensed under the same license as the project (to be specified).

---

## Attribution

Contributors will be acknowledged in:
- `CONTRIBUTORS.md` file
- Release notes for significant contributions
- Publications resulting from contributed features (as appropriate)

---

## Questions?

If you have questions about contributing, please:
1. Check existing documentation
2. Search closed issues for similar questions
3. Open a new issue with the "question" label
4. Contact maintainers directly for sensitive matters

---

Thank you for helping improve dānaSeq! Your contributions advance open science and enable better understanding of aquatic ecosystems.
