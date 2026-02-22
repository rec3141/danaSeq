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
# Resolve conda-provided tools (no host dependencies required)
# ============================================================================

# git: present in comebin, macsyfinder, defensefinder — use whichever is installed
GIT=""
for _env in dana-mag-comebin dana-mag-macsyfinder dana-mag-defensefinder; do
    if [[ -x "${ENV_DIR}/${_env}/bin/git" ]]; then
        GIT="${ENV_DIR}/${_env}/bin/git"
        break
    fi
done
[[ -n "$GIT" ]] || die "No conda git found in ${ENV_DIR} — run install.sh first"

# python3: present in bakta and flye
PYTHON3=""
for _env in dana-mag-bakta dana-mag-flye; do
    if [[ -x "${ENV_DIR}/${_env}/bin/python3" ]]; then
        PYTHON3="${ENV_DIR}/${_env}/bin/python3"
        break
    fi
done
[[ -n "$PYTHON3" ]] || die "No conda python3 found in ${ENV_DIR} — run install.sh first"

# ============================================================================
# LorBin PRs #6 + #7
# PR #6: sklearn.InvalidParameterError crash (n_neighbors=0), NameError in
#         mcluster (samplename → samplenames), argparse type= missing → TypeError
# PR #7: --cuda flag not forwarded to train_vae() → GPU training silently disabled
#
# Both PRs are on separate branches of rec3141/LorBin. We clone, merge them
# locally, and pip-install from the local clone into the semibin conda env.
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

    # Merge both fix branches locally and pip-install from the local clone.
    # (We don't push because we don't have write access to rec3141/LorBin.)
    local tmp
    tmp=$(mktemp -d)
    # shellcheck disable=SC2064
    trap "rm -rf '$tmp'" RETURN
    echo "  Cloning rec3141/LorBin..."
    "$GIT" clone "https://github.com/rec3141/LorBin.git" "$tmp/LorBin" \
        > /dev/null 2>&1
    "$GIT" -C "$tmp/LorBin" fetch origin \
        fix/sklearn-n-neighbors-zero fix/bin-cmd-cuda-flag > /dev/null 2>&1
    "$GIT" -C "$tmp/LorBin" checkout -b fix/combined-fixes \
        origin/fix/sklearn-n-neighbors-zero > /dev/null 2>&1
    GIT_AUTHOR_NAME="patches" GIT_AUTHOR_EMAIL="patches@localhost" \
    GIT_COMMITTER_NAME="patches" GIT_COMMITTER_EMAIL="patches@localhost" \
    "$GIT" -C "$tmp/LorBin" merge origin/fix/bin-cmd-cuda-flag \
        --no-edit > /dev/null 2>&1

    "${env}/bin/pip" install --force-reinstall --no-deps "$tmp/LorBin" \
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
    # Prepend conda git dir to PATH so pip can resolve the git+https:// URL
    PATH="$(dirname "$GIT"):${PATH}" \
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

    "$PYTHON3" - "$target" <<'PYEOF'
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
echo "  git:     ${GIT}"
echo "  python3: ${PYTHON3}"
echo ""

patch_lorbin
patch_whokaryote
patch_bakta

echo ""
echo "Done. Run './install.sh --check' to verify environments."
