/**
 * Two-finger gestures that keep the axes apart.
 *
 * A pinch is normally reported as one number, because on a photograph the two
 * directions mean the same thing. On these plots they do not. A phylogeny is
 * long in branch length and crowded in tips, and the two wants are opposite:
 * spread the tips apart without stretching the branches, or the reverse. A t-SNE
 * has no meaningful aspect ratio at all. So this reports the horizontal and
 * vertical scale separately and lets the caller decide whether to use both.
 *
 * Everything is measured against the state at gesture start rather than the
 * previous move, so a slow pinch and a fast one of the same span agree, and
 * rounding cannot accumulate across a long gesture.
 *
 * attachPinch(el, {
 *   onStart(),
 *   onPinch({ sx, sy, cx, cy }),   // scale factors ≥0, centre in element coords
 *   onEnd(),
 * }) -> detach()
 */

// Below this a "pinch" is a fumbled two-finger tap or a drag with a resting
// thumb, and reacting to it makes the plot jump.
const MIN_SPAN_PX = 24;

export function attachPinch(el, { onStart, onPinch, onEnd, capture = false } = {}) {
  if (!el) return () => {};

  let active = false;
  let startX = 0;   // horizontal span between the touches at gesture start
  let startY = 0;

  function spans(touches) {
    const [a, b] = touches;
    return {
      x: Math.abs(a.clientX - b.clientX),
      y: Math.abs(a.clientY - b.clientY),
      cx: (a.clientX + b.clientX) / 2,
      cy: (a.clientY + b.clientY) / 2,
    };
  }

  function onTouchStart(e) {
    if (e.touches.length !== 2) return;
    const s = spans(e.touches);
    // A pinch that is purely horizontal has no vertical span to measure from,
    // and dividing by it later would give infinity. Keep the axis pinned
    // instead by remembering that it started too small to judge.
    startX = s.x;
    startY = s.y;
    active = true;
    onStart?.();
  }

  function onTouchMove(e) {
    if (!active || e.touches.length !== 2) return;
    // The browser would otherwise scroll the page or run its own page zoom.
    e.preventDefault();
    // With `capture`, take the gesture before anything below sees it. The tree
    // is drawn by deck.gl, whose own controller reads a pinch as one uniform
    // zoom; letting both act would apply the gesture twice and lose the axes.
    if (capture) e.stopPropagation();
    const s = spans(e.touches);
    const rect = el.getBoundingClientRect();
    onPinch?.({
      sx: startX >= MIN_SPAN_PX ? s.x / startX : 1,
      sy: startY >= MIN_SPAN_PX ? s.y / startY : 1,
      cx: s.cx - rect.left,
      cy: s.cy - rect.top,
    });
  }

  function onTouchEnd(e) {
    if (!active || e.touches.length >= 2) return;
    active = false;
    onEnd?.();
  }

  // Not passive: onTouchMove has to be able to preventDefault, or the browser
  // takes the gesture for itself.
  const opts = { passive: false, capture };
  el.addEventListener('touchstart', onTouchStart, { passive: true, capture });
  el.addEventListener('touchmove', onTouchMove, opts);
  el.addEventListener('touchend', onTouchEnd, { passive: true, capture });
  el.addEventListener('touchcancel', onTouchEnd, { passive: true, capture });

  return () => {
    el.removeEventListener('touchstart', onTouchStart, { capture });
    el.removeEventListener('touchmove', onTouchMove, opts);
    el.removeEventListener('touchend', onTouchEnd, { capture });
    el.removeEventListener('touchcancel', onTouchEnd, { capture });
  };
}

/**
 * Rescale one axis range about a fixed point, given as a 0..1 position within
 * the plotting area. Zooming in means the visible range shrinks, so the scale
 * factor divides.
 */
export function zoomRange(range, scale, anchor01) {
  if (!Array.isArray(range) || range.length !== 2 || !(scale > 0)) return range;
  const [lo, hi] = range;
  const at = lo + (hi - lo) * anchor01;
  const span = (hi - lo) / scale;
  return [at - span * anchor01, at + span * (1 - anchor01)];
}
