// Ward hierarchical clustering → leaf order, for heatmap rows and columns.
//
// Nearest-neighbour chain with Lance–Williams updates on squared Euclidean
// distances: O(n²) time and memory, so a thousand rows cluster in well under
// a second. Rows are clustered as given; callers pass log-transformed values.

function sqDist(a, b) {
  let s = 0;
  for (let k = 0; k < a.length; k++) { const d = a[k] - b[k]; s += d * d; }
  return s;
}

/**
 * @param {number[][]} rows  n vectors of equal length
 * @returns {number[]}       permutation of 0..n-1 (dendrogram leaf order)
 */
export function wardOrder(rows) {
  const n = rows.length;
  if (n < 3) return rows.map((_, i) => i);

  // Condensed distance matrix (full square: simple indexing, n ≤ ~2000)
  const D = new Float64Array(n * n);
  for (let i = 0; i < n; i++) {
    for (let j = i + 1; j < n; j++) {
      const d = sqDist(rows[i], rows[j]);
      D[i * n + j] = d; D[j * n + i] = d;
    }
  }
  const size = new Array(n).fill(1);
  const active = new Uint8Array(n).fill(1);
  // Each active cluster's leaf list, concatenated as clusters merge
  const leaves = rows.map((_, i) => [i]);

  const chain = [];
  let remaining = n;
  while (remaining > 1) {
    if (!chain.length) { for (let i = 0; i < n; i++) if (active[i]) { chain.push(i); break; } }
    const a = chain[chain.length - 1];
    // nearest active neighbour of a (prefer the previous chain element on ties)
    let b = -1, best = Infinity;
    const prev = chain.length > 1 ? chain[chain.length - 2] : -1;
    for (let j = 0; j < n; j++) {
      if (j === a || !active[j]) continue;
      const d = D[a * n + j];
      if (d < best || (d === best && j === prev)) { best = d; b = j; }
    }
    if (b === prev) {
      // a and b are reciprocal nearest neighbours: merge b into a
      chain.pop(); chain.pop();
      const na = size[a], nb = size[b];
      for (let k = 0; k < n; k++) {
        if (!active[k] || k === a || k === b) continue;
        const nk = size[k];
        const d = ((na + nk) * D[a * n + k] + (nb + nk) * D[b * n + k] - nk * D[a * n + b]) / (na + nb + nk);
        D[a * n + k] = d; D[k * n + a] = d;
      }
      size[a] = na + nb;
      leaves[a] = leaves[a].concat(leaves[b]);
      leaves[b] = null;
      active[b] = 0;
      remaining--;
    } else {
      chain.push(b);
    }
  }
  for (let i = 0; i < n; i++) if (active[i]) return leaves[i];
  return rows.map((_, i) => i);
}

/** Row and column leaf orders for a matrix of non-negative values (log1p). */
export function clusterMatrix(matrix) {
  const logged = matrix.map(r => r.map(v => Math.log1p(v)));
  const rowOrder = wardOrder(logged);
  const cols = logged.length ? logged[0].map((_, j) => logged.map(r => r[j])) : [];
  const colOrder = wardOrder(cols);
  return { rowOrder, colOrder };
}
