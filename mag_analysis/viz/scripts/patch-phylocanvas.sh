#!/usr/bin/env bash
# patch-phylocanvas.sh — Install and patch Phylocanvas.gl for per-node label colors
#
# Phylocanvas.gl renders all leaf labels with a single global fontColour.
# This patch makes the leaf-label TextLayer use each node's fillColour
# (set via the `styles` prop) when available, falling back to fontColour.
# This enables gold labels for named reference species vs dim labels for
# unnamed placeholders, and colored labels for user MAG bins.
#
# The patch is two sed replacements on the minified UMD bundle:
#   1. Change getColor from static fontColour to per-node accessor
#   2. Add getColor to updateTriggers so colors update when styles change
#
# Usage: bash scripts/patch-phylocanvas.sh
#   Run from the viz/ directory after npm install.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
VIZ_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"

PHYLOCANVAS_VERSION="1.61.0"
BUNDLE_SRC="$VIZ_DIR/node_modules/@phylocanvas/phylocanvas.gl/dist/bundle.min.js"
BUNDLE_DST="$VIZ_DIR/public/phylocanvas-bundle.js"

# Install the package if not present
if [[ ! -f "$BUNDLE_SRC" ]]; then
    echo "Installing @phylocanvas/phylocanvas.gl@${PHYLOCANVAS_VERSION}..."
    (cd "$VIZ_DIR" && npm install --no-save "@phylocanvas/phylocanvas.gl@${PHYLOCANVAS_VERSION}" 2>/dev/null)
fi

if [[ ! -f "$BUNDLE_SRC" ]]; then
    echo "[ERROR] Phylocanvas bundle not found at $BUNDLE_SRC" >&2
    exit 1
fi

echo "Copying bundle to public/..."
cp "$BUNDLE_SRC" "$BUNDLE_DST"

echo "Applying per-node label color patch..."

# Patch 1: Replace static getColor with per-node accessor on leaf label TextLayers.
# Original:  getText:nr,getAngle:Xe.default,getColor:this.props.fontColour
# Patched:   getText:nr,getAngle:Xe.default,getColor:(function(fc){return function(d){return d.fillColour||fc}})(this.props.fontColour)
#
# This affects exactly 2 TextLayer instances: leaf-labels-text and leaf-labels-selected.
# Other TextLayers (metadata headers, scalebar) use different getText functions and are unaffected.
sed -i 's|getText:nr,getAngle:Xe.default,getColor:this\.props\.fontColour|getText:nr,getAngle:Xe.default,getColor:(function(fc){return function(d){return d.fillColour\|\|fc}})(this.props.fontColour)|g' "$BUNDLE_DST"

# Patch 2: Add getColor to updateTriggers so deck.gl re-evaluates colors when styles change.
# Original:  updateTriggers:{getPosition:this.props.updateTriggers.getTextPosition}
# Patched:   updateTriggers:{getPosition:this.props.updateTriggers.getTextPosition,getColor:this.props.styles}
sed -i 's|updateTriggers:{getPosition:this\.props\.updateTriggers\.getTextPosition}|updateTriggers:{getPosition:this.props.updateTriggers.getTextPosition,getColor:this.props.styles}|g' "$BUNDLE_DST"

# Strip sourcemap reference (the .map file is not included)
sed -i 's|//# sourceMappingURL=bundle.min.js.map||g' "$BUNDLE_DST"

# Verify patches applied
PATCH1_COUNT=$(grep -c 'function(fc){return function(d)' "$BUNDLE_DST" || true)
PATCH2_COUNT=$(grep -c 'getColor:this.props.styles' "$BUNDLE_DST" || true)

if [[ "$PATCH1_COUNT" -lt 1 ]]; then
    echo "[WARNING] Label color patch may not have applied (expected >=1 getColor accessor, found $PATCH1_COUNT)" >&2
    echo "  The bundle format may have changed. Check @phylocanvas/phylocanvas.gl version." >&2
else
    echo "  Patched $PATCH1_COUNT leaf-label getColor accessor(s)"
fi

if [[ "$PATCH2_COUNT" -lt 1 ]]; then
    echo "[WARNING] updateTriggers patch may not have applied (found $PATCH2_COUNT)" >&2
else
    echo "  Patched $PATCH2_COUNT updateTrigger(s)"
fi

echo "Done. Patched bundle: $BUNDLE_DST"
