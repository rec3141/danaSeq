// Exact t-SNE (van der Maaten & Hinton 2008) on a precomputed distance matrix.
// O(n^2) per iteration — meant for sample-level embeddings (tens to a few
// hundred points), not contigs. Seeded so the layout is stable across renders.

function mulberry32(seed) {
  return () => {
    seed |= 0; seed = (seed + 0x6d2b79f5) | 0;
    let t = Math.imul(seed ^ (seed >>> 15), 1 | seed);
    t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

// Joint probabilities P from distances, binary-searching each row's
// precision so its entropy matches log(perplexity).
function jointP(dist, n, perplexity) {
  const P = new Float64Array(n * n);
  const target = Math.log(perplexity);
  for (let i = 0; i < n; i++) {
    let beta = 1, lo = -Infinity, hi = Infinity;
    for (let iter = 0; iter < 64; iter++) {
      let sum = 0, h = 0;
      for (let j = 0; j < n; j++) {
        if (j === i) continue;
        const d2 = dist[i * n + j] ** 2;
        const p = Math.exp(-d2 * beta);
        P[i * n + j] = p; sum += p; h += d2 * p;
      }
      if (sum === 0) sum = 1e-12;
      const H = Math.log(sum) + beta * h / sum;
      for (let j = 0; j < n; j++) P[i * n + j] /= sum;
      const diff = H - target;
      if (Math.abs(diff) < 1e-5) break;
      if (diff > 0) { lo = beta; beta = hi === Infinity ? beta * 2 : (beta + hi) / 2; }
      else { hi = beta; beta = lo === -Infinity ? beta / 2 : (beta + lo) / 2; }
    }
  }
  // Symmetrize
  for (let i = 0; i < n; i++) {
    for (let j = i + 1; j < n; j++) {
      const v = Math.max((P[i * n + j] + P[j * n + i]) / (2 * n), 1e-12);
      P[i * n + j] = v; P[j * n + i] = v;
    }
  }
  return P;
}

/**
 * @param {number[][]} dist  symmetric n x n distance matrix
 * @param {{perplexity?: number, iterations?: number, seed?: number}} opts
 * @returns {number[][]} n x 2 coordinates
 */
export function tsne(dist, { perplexity = 30, iterations = 1000, seed = 42 } = {}) {
  const n = dist.length;
  if (n < 2) return dist.map(() => [0, 0]);
  const perp = Math.max(1, Math.min(perplexity, (n - 1) / 3));
  const flat = new Float64Array(n * n);
  for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) flat[i * n + j] = dist[i][j];
  const P = jointP(flat, n, perp);

  const rand = mulberry32(seed);
  const Y = new Float64Array(n * 2);
  for (let k = 0; k < n * 2; k++) Y[k] = (rand() - 0.5) * 1e-4;
  const vel = new Float64Array(n * 2);
  const gains = new Float64Array(n * 2).fill(1);
  const num = new Float64Array(n * n);
  const grad = new Float64Array(n * 2);
  const eta = Math.max(n / 12, 10);

  for (let it = 0; it < iterations; it++) {
    const exag = it < 250 ? 12 : 1;
    const mom = it < 250 ? 0.5 : 0.8;
    let sumQ = 0;
    for (let i = 0; i < n; i++) {
      for (let j = i + 1; j < n; j++) {
        const dx = Y[2 * i] - Y[2 * j], dy = Y[2 * i + 1] - Y[2 * j + 1];
        const q = 1 / (1 + dx * dx + dy * dy);
        num[i * n + j] = q; num[j * n + i] = q; sumQ += 2 * q;
      }
    }
    grad.fill(0);
    for (let i = 0; i < n; i++) {
      for (let j = 0; j < n; j++) {
        if (i === j) continue;
        const q = num[i * n + j];
        const m = 4 * (exag * P[i * n + j] - q / sumQ) * q;
        grad[2 * i] += m * (Y[2 * i] - Y[2 * j]);
        grad[2 * i + 1] += m * (Y[2 * i + 1] - Y[2 * j + 1]);
      }
    }
    for (let k = 0; k < n * 2; k++) {
      gains[k] = Math.sign(grad[k]) !== Math.sign(vel[k]) ? gains[k] + 0.2 : Math.max(gains[k] * 0.8, 0.01);
      vel[k] = mom * vel[k] - eta * gains[k] * grad[k];
      Y[k] += vel[k];
    }
    // Recenter
    let mx = 0, my = 0;
    for (let i = 0; i < n; i++) { mx += Y[2 * i]; my += Y[2 * i + 1]; }
    mx /= n; my /= n;
    for (let i = 0; i < n; i++) { Y[2 * i] -= mx; Y[2 * i + 1] -= my; }
  }
  return Array.from({ length: n }, (_, i) => [Y[2 * i], Y[2 * i + 1]]);
}

/** Bray-Curtis dissimilarity matrix between rows of a non-negative matrix. */
export function brayCurtis(rows) {
  const n = rows.length;
  const D = Array.from({ length: n }, () => new Array(n).fill(0));
  for (let i = 0; i < n; i++) {
    for (let j = i + 1; j < n; j++) {
      let num = 0, den = 0;
      const a = rows[i], b = rows[j];
      for (let k = 0; k < a.length; k++) { num += Math.abs(a[k] - b[k]); den += a[k] + b[k]; }
      const d = den > 0 ? num / den : 0;
      D[i][j] = d; D[j][i] = d;
    }
  }
  return D;
}
