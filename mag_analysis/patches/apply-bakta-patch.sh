#!/bin/bash
# Patch Bakta's pseudogene detection so it finishes on whole metagenome assemblies.
#
#   apply-bakta-patch.sh <conda env prefix>
#
# Bakta 1.12.2 checks each pseudogene for a paralog by rebuilding a set of every
# CDS's UniRef90 cluster, and parses the DIAMOND blastx XML into memory whole. On
# a 1.5M-contig co-assembly that ran single-threaded for over four days in 320 GB
# and did not finish. The patch counts the clusters once and streams the XML; the
# results are unchanged. Already-patched installs are left alone, and a Bakta the
# patch does not apply to fails here rather than silently running unpatched.
set -euo pipefail
ENV=${1:?usage: apply-bakta-patch.sh <conda env prefix>}
PATCH="$(cd "$(dirname "$0")" && pwd)/bakta-1.12.2-pseudogenes.patch"
SITE=$("$ENV/bin/python" -c "import importlib.util, os; print(os.path.dirname(importlib.util.find_spec('bakta').submodule_search_locations[0]))")
CDS="$SITE/bakta/features/cds.py"
if grep -q "uniref90_counts" "$CDS"; then
    echo "  Bakta pseudogene patch already applied"
    exit 0
fi
patch --dry-run -p1 -d "$SITE" < "$PATCH" > /dev/null || {
    echo "[ERROR] the Bakta pseudogene patch does not apply to $CDS; update $PATCH" >&2
    exit 1
}
patch -p1 -d "$SITE" < "$PATCH" > /dev/null
echo "  Bakta pseudogene patch applied"
