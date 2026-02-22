#!/usr/bin/env bash
# patches.sh — Apply upstream bug-fix patches to installed conda environments
#
# These are fixes for open PRs on tool repositories that haven't been merged
# upstream yet. Run once after install.sh. Safe to re-run (idempotent).
#
# Tracked PRs:
#   LorMeBioAI/LorBin#6       - Fix crashes: sklearn n_neighbors=0, NameError in
#                               mcluster, missing argparse type= annotations
#   LorMeBioAI/LorBin#7       - Fix: pass --cuda flag to train_vae in bin subcommand
#   LottePronk/whokaryote#13  - Fix GFF parsing to support standard GFF3 (Bakta/PGAP)
#   oschwengers/bakta#429     - Make AMRFinderPlus failure non-fatal (warn + return)
#
# Usage:
#   ./patches.sh
#   ./patches.sh --prefix /custom/conda-envs

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ENV_DIR="${SCRIPT_DIR}/conda-envs"

while (( $# )); do
    case "$1" in
        --prefix) ENV_DIR="$2"; shift 2 ;;
        *) echo "[ERROR] Unknown argument: $1" >&2; exit 1 ;;
    esac
done

die()  { echo "[ERROR] $1" >&2; exit 1; }
info() { echo "[PATCH] $1"; }
skip() { echo "[SKIP]  $1"; }
ok()   { echo "  Done"; }

# ============================================================================
# LorBin PRs #6 + #7
# PR #6: sklearn.InvalidParameterError crash (n_neighbors=0), NameError in
#         mcluster (samplename → samplenames), argparse type= missing → TypeError
# PR #7: --cuda flag not forwarded to train_vae() → GPU training silently disabled
#
# Both PRs are on separate branches of rec3141/LorBin. We merge them into a
# single combined branch and pip-install from it into the semibin conda env.
# ============================================================================
patch_lorbin() {
    local env="${ENV_DIR}/dana-mag-semibin"
    [[ -x "${env}/bin/pip" ]] || { skip "dana-mag-semibin not installed"; return; }

    # Idempotency check: PR #7's change is the last applied (cuda flag)
    local installed_lorbin
    installed_lorbin=$(find "${env}" -name "lorbin.py" -path "*/lorbin/*" 2>/dev/null | head -1 || true)
    if [[ -n "$installed_lorbin" ]] && grep -q 'is_cuda=args.cuda' "$installed_lorbin" 2>/dev/null; then
        skip "lorbin already patched (PRs #6 + #7)"
        return
    fi

    info "lorbin PRs #6 + #7: sklearn crash + cuda flag + argparse types..."

    # Create combined fix branch in rec3141/LorBin if it doesn't exist
    if ! git ls-remote "https://github.com/rec3141/LorBin.git" \
            "refs/heads/fix/combined-fixes" 2>/dev/null | grep -q combined-fixes; then
        echo "  Creating rec3141/LorBin@fix/combined-fixes..."
        local tmp
        tmp=$(mktemp -d)
        # shellcheck disable=SC2064
        trap "rm -rf '$tmp'" RETURN
        git clone "https://github.com/rec3141/LorBin.git" "$tmp/LorBin" \
            > /dev/null 2>&1
        git -C "$tmp/LorBin" fetch origin \
            fix/sklearn-n-neighbors-zero fix/bin-cmd-cuda-flag > /dev/null 2>&1
        git -C "$tmp/LorBin" checkout -b fix/combined-fixes \
            origin/fix/sklearn-n-neighbors-zero > /dev/null 2>&1
        git -C "$tmp/LorBin" merge origin/fix/bin-cmd-cuda-flag \
            --no-edit > /dev/null 2>&1
        git -C "$tmp/LorBin" push origin fix/combined-fixes > /dev/null 2>&1
        echo "  Created rec3141/LorBin@fix/combined-fixes"
    fi

    "${env}/bin/pip" install --force-reinstall --no-deps \
        "lorbin @ git+https://github.com/rec3141/LorBin.git@fix/combined-fixes" \
        > /dev/null 2>&1
    ok
}

# ============================================================================
# Whokaryote PR #13
# The pipeline uses Bakta for annotation, whose GFF3 output uses standard
# ##sequence-region pragmas and Name= attributes — not Prodigal's custom
# comment headers. Without this fix, Whokaryote silently produces no output
# when called with --gff from Bakta.
# ============================================================================
patch_whokaryote() {
    local env="${ENV_DIR}/dana-mag-whokaryote"
    [[ -x "${env}/bin/pip" ]] || { skip "dana-mag-whokaryote not installed"; return; }

    # Idempotency check: patched version parses ##sequence-region
    local installed_feat
    installed_feat=$(find "${env}" -name "calculate_features.py" \
        -path "*/whokaryote_scripts/*" 2>/dev/null | head -1 || true)
    if [[ -n "$installed_feat" ]] && grep -q 'sequence-region' "$installed_feat" 2>/dev/null; then
        skip "whokaryote already patched (PR #13)"
        return
    fi

    info "whokaryote PR #13: Fix GFF3 parsing for Bakta/PGAP output..."
    "${env}/bin/pip" install --force-reinstall --no-deps \
        "git+https://github.com/rec3141/whokaryote.git@fix/gff-attribute-parsing" \
        > /dev/null 2>&1
    ok
}

# ============================================================================
# Bakta PR #429
# AMRFinderPlus often has database/version mismatch issues. Without this patch,
# any AMRFinderPlus failure raises an exception and crashes the entire Bakta
# annotation run. The fix: log a warning and return an empty set instead.
#
# Applied as a direct in-place patch to the installed amrfinder.py, since
# Bakta is conda-installed and pip-reinstalling the full package risks
# dependency version skew.
# ============================================================================
patch_bakta() {
    local env="${ENV_DIR}/dana-mag-bakta"
    [[ -x "${env}/bin/python" ]] || { skip "dana-mag-bakta not installed"; return; }

    local target
    target=$(find "${env}" -name "amrfinder.py" -path "*/bakta/expert/*" 2>/dev/null | head -1 || true)
    [[ -n "$target" ]] || { skip "bakta amrfinder.py not found in ${env}"; return; }

    # Idempotency check
    if grep -q 'AMR annotation will be skipped' "$target" 2>/dev/null; then
        skip "bakta already patched (PR #429)"
        return
    fi

    info "bakta PR #429: Make AMRFinderPlus failure non-fatal..."

    python3 - "$target" <<'PYEOF'
import sys
path = sys.argv[1]
with open(path) as f:
    lines = f.readlines()

# Find and replace:
#   raise Exception(f"amrfinder error! ...")
# with:
#   log.warning('AMRFinderPlus failed - AMR annotation will be skipped.')
#   return set()
old = '        raise Exception(f"amrfinder error!'
new_lines = [
    "        log.warning('AMRFinderPlus failed - AMR annotation will be skipped. "
    "Try: amrfinder_update --force_update --database <db_path>')\n",
    "        return set()\n",
]

patched = False
out = []
for line in lines:
    if line.startswith(old):
        out.extend(new_lines)
        patched = True
    else:
        out.append(line)

if not patched:
    print(f"  WARNING: target line not found in {path} — check bakta version", file=sys.stderr)
    sys.exit(1)

with open(path, 'w') as f:
    f.writelines(out)
print(f"  Patched {path}")
PYEOF
    ok
}

# ============================================================================
# Main
# ============================================================================

echo ""
echo "Applying upstream patches to conda environments in: ${ENV_DIR}"
echo ""

patch_lorbin
patch_whokaryote
patch_bakta

echo ""
echo "Done. Run './install.sh --check' to verify environments."
