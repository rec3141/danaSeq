#!/bin/bash
# Check for corrupted FASTA files with NULL bytes or malformed headers

OUTPUT_DIR="${1:-.}"

echo "Checking for FASTA corruption in: $OUTPUT_DIR"
echo ""

CORRUPTED=0
TRUNCATED=0
MISSING_HEADERS=0
NULL_BYTES=0

find "$OUTPUT_DIR" -name "*.fa" -type f | while read fafile; do
  # Check for NULL bytes
  if grep -q $'\x00' "$fafile" 2>/dev/null; then
    echo "⚠️  NULL BYTES: $fafile"
    ((NULL_BYTES++))

    # Show first occurrence
    echo "   First NULL byte at line:"
    head -20 "$fafile" | cat -v | head -5
    echo ""
    continue
  fi

  # Check if file is empty
  if [[ ! -s "$fafile" ]]; then
    echo "⚠️  EMPTY: $fafile"
    continue
  fi

  # Check if ends with header (incomplete)
  if tail -1 "$fafile" | grep -q '^>'; then
    echo "⚠️  TRUNCATED: $fafile (ends with header, missing sequence)"
    ((TRUNCATED++))
    continue
  fi

  # Check if has any headers
  if ! grep -q '^>' "$fafile"; then
    echo "⚠️  NO HEADERS: $fafile (not a valid FASTA)"
    ((MISSING_HEADERS++))
    continue
  fi

  # Check for headers without sequence
  if awk '/^>/ {if (seq) print seq; seq=0; next} {seq=1} END {if (!seq) exit 1}' "$fafile"; then
    : # OK
  else
    echo "⚠️  HEADER WITHOUT SEQUENCE: $fafile"
    continue
  fi

  # Check for sequence before first header
  first_line=$(head -1 "$fafile")
  if [[ ! "$first_line" =~ ^\> ]]; then
    echo "⚠️  SEQUENCE BEFORE HEADER: $fafile"
    echo "   First line: ${first_line:0:80}"
    ((CORRUPTED++))
    continue
  fi
done

echo ""
echo "Summary:"
echo "  NULL bytes: $NULL_BYTES"
echo "  Truncated: $TRUNCATED"
echo "  Missing headers: $MISSING_HEADERS"
echo "  Other corruption: $CORRUPTED"

if (( NULL_BYTES + TRUNCATED + MISSING_HEADERS + CORRUPTED == 0 )); then
  echo "✅ All FASTA files appear intact!"
else
  echo ""
  echo "To fix corrupted files:"
  echo "  1. Identify the source file"
  echo "  2. Remove the corrupted FASTA: rm path/to/bad.fa"
  echo "  3. Re-run pipeline: ./24_process_reads_optimized.sh -i data -K"
fi
