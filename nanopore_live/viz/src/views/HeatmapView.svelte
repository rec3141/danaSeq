<script>
  import PlotlyChart from '../components/charts/PlotlyChart.svelte';
  import { samples, sampleTaxonomy, metadata } from '../stores/data.js';
  import { cartItems, cartActive } from '../stores/cart.js';
  import { sampleClusters, sampleClusterK } from '../stores/clusters.js';
  import { paletteColor } from '../stores/taxonomy.js';

  const RANK_LEVELS = [
    { code: 'P', label: 'Phylum' },
    { code: 'C', label: 'Class' },
    { code: 'O', label: 'Order' },
    { code: 'F', label: 'Family' },
    { code: 'G', label: 'Genus' },
    { code: 'S', label: 'Species' },
  ];

  // Start at Genus (index 4) — most interpretable granularity for site comparison.
  let rankIdx = $state(4);
  let topN = $state(30);
  let valueMode = $state('pct');  // 'pct' | 'log'
  // Number of clusters cut from the Ward dendrogram. Slider range adapts to
  // the active sample count; the actual k used is min(kClusters, N).
  let kClusters = $state(4);

  let currentRank = $derived(RANK_LEVELS[rankIdx]);

  function cycleRank() { rankIdx = (rankIdx + 1) % RANK_LEVELS.length; }
  function cycleValue() { valueMode = valueMode === 'pct' ? 'log' : 'pct'; }

  let activeSamples = $derived.by(() => {
    if (!$samples) return [];
    if ($cartActive && $cartItems.size > 0) return $samples.filter(s => $cartItems.has(s.id));
    return $samples;
  });

  // Children map (inverted parents) for fast descendant walk at any rank.
  let childrenMap = $derived.by(() => {
    const parents = $sampleTaxonomy?.parents;
    if (!parents) return {};
    const out = {};
    for (const [child, parent] of Object.entries(parents)) {
      (out[parent] ||= []).push(child);
    }
    return out;
  });

  // Taxa at the current rank keyed by name.
  let taxaAtRank = $derived.by(() => {
    const ranks = $sampleTaxonomy?.ranks;
    if (!ranks) return [];
    const out = [];
    for (const [name, r] of Object.entries(ranks)) if (r === currentRank.code) out.push(name);
    return out;
  });

  // Descendant set per taxon at the current rank, used to roll up per-sample
  // direct counts into this rank's buckets.
  let descendantSets = $derived.by(() => {
    const kids = childrenMap;
    const out = {};
    function walk(taxon, acc) {
      acc.push(taxon);
      const cs = kids[taxon];
      if (cs) for (const c of cs) walk(c, acc);
    }
    for (const t of taxaAtRank) {
      const acc = [];
      walk(t, acc);
      out[t] = acc;
    }
    return out;
  });

  // Per-sample rollup: { sampleId: { taxon_at_rank: aggregated_count } }.
  let perSampleRankCounts = $derived.by(() => {
    const out = {};
    const samplesData = $sampleTaxonomy?.samples || {};
    for (const sid in samplesData) {
      const direct = samplesData[sid]?.direct || {};
      if (currentRank.code === 'P') { out[sid] = { ...(samplesData[sid].phylum || {}) }; continue; }
      if (currentRank.code === 'C') { out[sid] = { ...(samplesData[sid].class || {}) }; continue; }
      const row = {};
      for (const [taxon, descendants] of Object.entries(descendantSets)) {
        let n = 0;
        for (const d of descendants) if (direct[d]) n += direct[d];
        if (n > 0) row[taxon] = n;
      }
      out[sid] = row;
    }
    return out;
  });

  // Top-N taxa across active samples, ranked by summed counts.
  let topTaxa = $derived.by(() => {
    const totals = {};
    for (const s of activeSamples) {
      const row = perSampleRankCounts[s.id] || {};
      for (const [taxon, count] of Object.entries(row)) totals[taxon] = (totals[taxon] || 0) + count;
    }
    return Object.entries(totals)
      .sort((a, b) => b[1] - a[1])
      .slice(0, topN)
      .map(([name, count]) => ({ name, count }));
  });

  // Per-sample denominator: sum of classified reads at current rank (so the
  // heatmap cells are "% of reads classified into this rank" — a within-sample
  // fraction that's directly comparable across sites without depth bias).
  function sampleDenom(sid) {
    const row = perSampleRankCounts[sid] || {};
    let s = 0;
    for (const v of Object.values(row)) s += v;
    return s || 1;
  }

  // ---- Column ordering ---------------------------------------------------
  //
  // The order button cycles through: default (sample id) → each metadata
  // column → Ward clustering. Ward runs on squared Euclidean distances over
  // 4th-root-transformed relative abundances — a standard downweighting of
  // dominant taxa used in community ecology (Anderson & Willis 2003; Chao
  // et al., equivalent to the γ=¼ Tukey ladder). Leaves come out in the
  // dendrogram's traversal order so columns that cluster together sit
  // side-by-side in the heatmap.

  // Detect metadata columns present across active samples.
  let metaColumns = $derived.by(() => {
    if (!$metadata) return [];
    const keys = new Set();
    for (const s of activeSamples) {
      const m = $metadata[s.id];
      if (!m) continue;
      for (const k of Object.keys(m)) {
        if (k === 'lat' || k === 'lon') continue;
        keys.add(k);
      }
    }
    return [...keys].sort();
  });

  // Order modes: 'default' (sample id) | 'meta:<col>' | 'cluster'
  let orderModes = $derived([
    { id: 'default', label: 'Default' },
    ...metaColumns.map(c => ({ id: `meta:${c}`, label: c })),
    { id: 'cluster', label: 'Ward' },
  ]);
  let orderIdx = $state(0);
  let currentOrder = $derived(orderModes[orderIdx % Math.max(1, orderModes.length)]);

  function cycleOrder() { orderIdx = (orderIdx + 1) % orderModes.length; }

  // Label modes are independent from Order so sorting by a metadata column
  // doesn't force that column onto the x-axis tick (and vice versa). Default
  // label is FLOWCELL:barcodeNN.
  let labelModes = $derived([
    { id: 'default', label: 'FC:bc' },
    ...metaColumns.map(c => ({ id: `meta:${c}`, label: c })),
  ]);
  let labelIdx = $state(0);
  let currentLabel = $derived(labelModes[labelIdx % Math.max(1, labelModes.length)]);

  function cycleLabel() { labelIdx = (labelIdx + 1) % labelModes.length; }

  // Ward linkage on 4th-root-transformed relative abundances. The 4th-root
  // transform (x^¼) compresses the contribution of dominant taxa so rare
  // taxa inform the geometry more strongly — standard practice for
  // clustering community composition. Ward's Lance-Williams update
  // operates on squared Euclidean distances between cluster centroids
  // (in sample space); we carry those through the merges directly.
  //
  // Returns:
  //   {
  //     order:   [sampleId, ...] in dendrogram-leaf order,
  //     linkage: [{ left, right, distance, size, leaves: [sampleId,...] }, ...]
  //              one entry per merge (n-1 total). `left`/`right` are indices:
  //              0..n-1 reference original leaves (in the input sampleIds
  //              order); n..n+merges-1 reference earlier merges in this list.
  //   }
  function wardClustering(sampleIds) {
    const n = sampleIds.length;
    if (n === 0) return { order: [], linkage: [] };
    if (n === 1) return { order: sampleIds.slice(), linkage: [] };

    // Build per-sample vectors in a shared taxon order. Count-of-zero taxa
    // still contribute (x^¼ = 0), just without cost — no need to sparsify.
    const taxaSet = new Set();
    for (const sid of sampleIds) {
      for (const t of Object.keys(perSampleRankCounts[sid] || {})) taxaSet.add(t);
    }
    const taxa = [...taxaSet];
    const vecs = sampleIds.map(sid => {
      const counts = perSampleRankCounts[sid] || {};
      const denom = sampleDenom(sid);
      return taxa.map(t => Math.pow((counts[t] || 0) / denom, 0.25));
    });

    const sqEuc = (a, b) => {
      let s = 0;
      for (let k = 0; k < a.length; k++) { const d = a[k] - b[k]; s += d * d; }
      return s;
    };

    // `clusters[i]` is the live cluster currently sitting at slot i. As we
    // merge, the surviving slot keeps the smaller index; the other becomes
    // dead. `nodeId[i]` is the linkage-array node id (0..n-1 = leaves;
    // n+m = the m'th merge we record below) that this slot represents.
    const clusters = sampleIds.map(sid => ({ size: 1, leaves: [sid] }));
    const nodeId = sampleIds.map((_, i) => i);
    const linkage = [];

    const D = new Map();
    const key = (i, j) => i < j ? `${i},${j}` : `${j},${i}`;
    for (let i = 0; i < n; i++) {
      for (let j = i + 1; j < n; j++) D.set(key(i, j), sqEuc(vecs[i], vecs[j]));
    }
    const live = new Set([...Array(n).keys()]);
    while (live.size > 1) {
      let minD = Infinity, pi = -1, pj = -1;
      const arr = [...live];
      for (let a = 0; a < arr.length; a++) {
        for (let b = a + 1; b < arr.length; b++) {
          const d = D.get(key(arr[a], arr[b])) ?? Infinity;
          if (d < minD) { minD = d; pi = arr[a]; pj = arr[b]; }
        }
      }
      const ni = clusters[pi].size, nj = clusters[pj].size;
      const dij = D.get(key(pi, pj)) ?? 0;
      // Lance-Williams coefficients for Ward on squared Euclidean:
      //   d(uv, k) = ((n_u + n_k) d(u,k) + (n_v + n_k) d(v,k) - n_k d(u,v))
      //             / (n_u + n_v + n_k)
      for (const k of live) {
        if (k === pi || k === pj) continue;
        const nk = clusters[k].size;
        const dik = D.get(key(pi, k)) ?? 0;
        const djk = D.get(key(pj, k)) ?? 0;
        const dnew = ((ni + nk) * dik + (nj + nk) * djk - nk * dij) / (ni + nj + nk);
        D.set(key(pi, k), dnew);
      }
      const mergedLeaves = [...clusters[pi].leaves, ...clusters[pj].leaves];
      const newNodeId = n + linkage.length;
      // Record the merge (sqrt for a nicer height scale — squared-Euclidean
      // grows fast and dwarfs early merges otherwise).
      linkage.push({
        left: nodeId[pi],
        right: nodeId[pj],
        distance: Math.sqrt(Math.max(0, dij)),
        size: ni + nj,
        leaves: mergedLeaves,
      });
      clusters[pi] = { size: ni + nj, leaves: mergedLeaves };
      nodeId[pi] = newNodeId;
      live.delete(pj);
    }
    return { order: clusters[[...live][0]].leaves, linkage };
  }

  // Cut a linkage tree into exactly k clusters. Strategy: undo the last
  // (k-1) merges — i.e. consider the top (k-1) merges "broken", which
  // splits the single root into k disjoint sub-trees. Returns
  // { [sampleId]: 'C1'|'C2'|... } numbered in dendrogram leaf order so
  // adjacent clusters on the heatmap have adjacent labels.
  function cutLinkage(linkage, leafOrder, k) {
    const n = leafOrder.length;
    if (n === 0) return {};
    if (k <= 1 || linkage.length === 0) {
      const out = {};
      for (const sid of leafOrder) out[sid] = 'C1';
      return out;
    }
    const kk = Math.min(k, n);

    // Walk down: the top merge (linkage[L-1]) is the root; "breaking" it
    // gives 2 clusters (its left/right sub-trees). Breaking the next
    // highest gives 3, etc. After breaking the top (kk-1) merges, the
    // surviving merges are linkage[0 .. L-kk]. Roots of the kk sub-trees
    // are the nodes that never appear as a child within those surviving
    // merges (and aren't the broken-merges' parents either).
    const L = linkage.length;
    const survivingLast = L - (kk - 1);  // exclusive upper bound
    const childIds = new Set();
    for (let i = 0; i < survivingLast; i++) {
      childIds.add(linkage[i].left);
      childIds.add(linkage[i].right);
    }
    // Candidate roots: every node id (leaf 0..n-1, plus merge n..n+L-1)
    // referenced by surviving merges OR (for k=n) every leaf, that does
    // NOT appear as a child within surviving merges. Equivalent: the
    // children of broken merges that aren't themselves children of any
    // surviving merge.
    const subRoots = [];
    if (kk >= n) {
      // Each leaf is its own cluster.
      for (let i = 0; i < n; i++) subRoots.push(i);
    } else {
      // Children of broken merges, filtered to those not absorbed into a
      // surviving merge. Walk broken merges top-down (highest first) and
      // collect children that aren't subsumed.
      const seen = new Set();
      const candidates = [];
      for (let i = L - 1; i >= survivingLast; i--) {
        candidates.push(linkage[i].left);
        candidates.push(linkage[i].right);
      }
      for (const c of candidates) {
        if (childIds.has(c)) continue;  // absorbed by a surviving merge
        if (seen.has(c)) continue;
        seen.add(c);
        subRoots.push(c);
      }
      // Edge case: if subRoots is short of kk (shouldn't happen with
      // well-formed binary linkage but be defensive), pad with the first
      // unassigned leaves. Keeps cluster count == kk.
      if (subRoots.length < kk) {
        const usedLeaves = new Set();
        for (const r of subRoots) {
          if (r < n) usedLeaves.add(r);
          else for (const sid of linkage[r - n].leaves) usedLeaves.add(leafOrder.indexOf(sid));
        }
        for (let i = 0; i < n && subRoots.length < kk; i++) {
          if (!usedLeaves.has(i)) subRoots.push(i);
        }
      }
    }

    // Resolve each sub-root to its leaf set (sample ids).
    function leavesOf(nodeIdx) {
      if (nodeIdx < n) return [leafOrder[nodeIdx]];
      const m = linkage[nodeIdx - n];
      // m.leaves is already in dendrogram order — just return it.
      return m.leaves;
    }

    // Sort sub-roots by the dendrogram position of their first leaf so
    // cluster numbering goes left-to-right on the heatmap.
    const orderIndex = new Map();
    leafOrder.forEach((sid, i) => orderIndex.set(sid, i));
    const subWithLeaves = subRoots.map(r => {
      const leaves = leavesOf(r);
      const minIdx = Math.min(...leaves.map(sid => orderIndex.get(sid) ?? 1e9));
      return { leaves, minIdx };
    });
    subWithLeaves.sort((a, b) => a.minIdx - b.minIdx);

    const out = {};
    subWithLeaves.forEach((entry, i) => {
      const label = `C${i + 1}`;
      for (const sid of entry.leaves) out[sid] = label;
    });
    return out;
  }

  // Cached Ward linkage on the full active-sample set, keyed implicitly by
  // the (sample set, current rank) pair via Svelte's reactivity. Reused by
  // the dendrogram, the cluster color bar, and the (mode='cluster')
  // ordering branch — so the heavy O(n^3) merge runs once per rank/sample
  // change rather than once per consumer.
  let wardResult = $derived.by(() => {
    if (!activeSamples.length) return { order: [], linkage: [] };
    const ids = activeSamples.map(s => s.id);
    return wardClustering(ids);
  });

  let orderedSamples = $derived.by(() => {
    if (!activeSamples.length) return [];
    const meta = $metadata || {};
    const mode = currentOrder.id;
    if (mode === 'default') {
      return [...activeSamples].sort((a, b) => (a.id || '').localeCompare(b.id || ''));
    }
    if (mode === 'cluster') {
      const order = wardResult.order;
      const byId = new Map(activeSamples.map(s => [s.id, s]));
      return order.map(id => byId.get(id)).filter(Boolean);
    }
    // meta:<col> — alphanumeric sort with stable (id) tiebreak; numeric
    // values compared as numbers so temperature/depth sort naturally.
    const col = mode.slice(5);
    const getVal = (s) => {
      const v = meta[s.id]?.[col];
      if (v == null) return '';
      return typeof v === 'number' ? v : String(v).replace(/^"|"$/g, '');
    };
    return [...activeSamples].sort((a, b) => {
      const va = getVal(a), vb = getVal(b);
      if (typeof va === 'number' && typeof vb === 'number') {
        if (va !== vb) return va - vb;
      } else {
        const sa = String(va), sb = String(vb);
        const c = sa.localeCompare(sb);
        if (c !== 0) return c;
      }
      return (a.id || '').localeCompare(b.id || '');
    });
  });

  // ---- Cluster assignments + dendrogram coords ---------------------------

  // Slider's upper bound: 2..min(N, 20). Clamps the user-set k to whatever
  // the current sample count permits.
  let kMax = $derived(Math.min(Math.max(2, activeSamples.length), 20));
  let effectiveK = $derived(Math.max(1, Math.min(kClusters, activeSamples.length)));

  // Cluster id (e.g. 'C1') per sample, keyed by sample id. Shared with the
  // global sampleClusters store so other views can color-by-cluster.
  let clusterAssign = $derived.by(() => {
    if (!wardResult.order.length) return {};
    return cutLinkage(wardResult.linkage, wardResult.order, effectiveK);
  });

  // Sync the cluster assignment into the global store whenever it changes.
  // Other views (SampleExplorerView's metadata cycle, /map, /reads) can
  // import sampleClusters and color by it.
  $effect(() => {
    sampleClusters.set(clusterAssign);
    const ks = new Set(Object.values(clusterAssign));
    sampleClusterK.set(ks.size);
  });

  // Build dendrogram line segments in (x = leaf-position, y = distance)
  // coordinates. Returns an array of {x: [...], y: [...]} polylines, each
  // an inverted-U connecting two child cluster centers up at the merge's
  // distance. We build this in *original* sample order, so the consumer
  // re-projects x onto the column index of the displayed (orderedSamples)
  // sequence.
  let dendroSegments = $derived.by(() => {
    const { order, linkage } = wardResult;
    if (order.length < 2) return [];
    // Position of each leaf along the dendrogram x-axis = its index in
    // wardResult.order (the natural, Ward-imposed leaf ordering).
    const leafPos = new Map();
    order.forEach((sid, i) => leafPos.set(sid, i));
    // Each merge's center x = midpoint of its leaves' positions (mean of
    // first + last works for the dendrogram-leaf order, because each
    // merge's leaves are contiguous in that ordering).
    const nodeX = new Map();
    order.forEach((sid, i) => nodeX.set(i, i));  // leaf nodes 0..n-1
    const n = order.length;
    const segs = [];
    for (let m = 0; m < linkage.length; m++) {
      const merge = linkage[m];
      const lx = nodeX.get(merge.left);
      const rx = nodeX.get(merge.right);
      // Use leaves of each child cluster's parent merge (or the leaf
      // itself) to compute the merge center as the mean of its first +
      // last leaf positions — robust to non-binary mergers (we don't
      // produce those, but it's a clean general formula).
      const leftLeaves = merge.left < n ? [order[merge.left]] : linkage[merge.left - n].leaves;
      const rightLeaves = merge.right < n ? [order[merge.right]] : linkage[merge.right - n].leaves;
      const leftCenter = (leafPos.get(leftLeaves[0]) + leafPos.get(leftLeaves[leftLeaves.length - 1])) / 2;
      const rightCenter = (leafPos.get(rightLeaves[0]) + leafPos.get(rightLeaves[rightLeaves.length - 1])) / 2;
      // Vertical drops to each child's height + horizontal bar at this
      // merge's height. Children's heights are linkage[child-n].distance
      // for merges, or 0 for leaves.
      const lh = merge.left < n ? 0 : linkage[merge.left - n].distance;
      const rh = merge.right < n ? 0 : linkage[merge.right - n].distance;
      const h = merge.distance;
      // U-shape: (lx, lh) → (lx, h) → (rx, h) → (rx, rh)
      segs.push({ x: [leftCenter, leftCenter, rightCenter, rightCenter], y: [lh, h, h, rh] });
      // This merge's center for its parent reference.
      nodeX.set(n + m, (leftCenter + rightCenter) / 2);
    }
    return segs;
  });

  let dendroMaxHeight = $derived.by(() => {
    let h = 0;
    for (const m of wardResult.linkage) if (m.distance > h) h = m.distance;
    return h || 1;
  });

  // ---- Heatmap traces + layout -------------------------------------------

  // Single-line x-tick label controlled by `currentLabel` — never stacks two
  // strings so Plotly doesn't render overlapping primary/subline pairs. Falls
  // back to FLOWCELL:barcodeNN when a chosen metadata column is missing for
  // a given sample so no tick is ever blank.
  function xDisplayLabel(s) {
    const fc = s.flowcell && s.barcode ? `${s.flowcell}:${s.barcode}` : s.id;
    const mode = currentLabel.id;
    if (mode.startsWith('meta:')) {
      const col = mode.slice(5);
      const v = $metadata?.[s.id]?.[col];
      if (v != null && v !== '') return String(v).replace(/^"|"$/g, '');
    }
    return fc;
  }

  let heatmapTrace = $derived.by(() => {
    if (!topTaxa.length || !orderedSamples.length) return null;
    // Numeric indices for x — avoids Plotly's category-axis quirk where
    // both the auto category label AND our tickmode-array ticktext render
    // as overlapping (white + grey) labels. Ticks get their display name
    // via tickvals/ticktext only.
    const xVals = orderedSamples.map((_, i) => i);
    const xDisplay = orderedSamples.map(xDisplayLabel);
    const xHoverText = orderedSamples.map(s => xDisplayLabel(s).replace(/<[^>]*>/g, ''));
    const yLabels = topTaxa.map(t => t.name);
    const z = topTaxa.map(taxon =>
      orderedSamples.map(s => {
        const count = (perSampleRankCounts[s.id] || {})[taxon.name] ?? 0;
        if (valueMode === 'log') return count > 0 ? Math.log10(count) : null;
        return (count / sampleDenom(s.id)) * 100;
      })
    );
    // customdata: one value per column, broadcast per row by Plotly.
    const customdata = topTaxa.map(() => xHoverText);

    const traces = [{
      type: 'heatmap',
      x: xVals,
      y: yLabels,
      z,
      xaxis: 'x',
      yaxis: 'y',
      xgap: 0,
      ygap: 0,
      colorscale: 'Viridis',
      hoverongaps: false,
      colorbar: {
        title: { text: valueMode === 'log' ? 'log10(reads)' : '% of classified', font: { color: '#94a3b8', size: 11 } },
        tickfont: { color: '#94a3b8', size: 10 },
        thickness: 12,
        // Pin the colorbar to the heatmap subplot only (not the top
        // dendrogram subplot) by anchoring its y to the heatmap domain.
        len: 0.65, y: 0, yanchor: 'bottom',
      },
      customdata,
      hovertemplate: (valueMode === 'log'
        ? '<b>%{y}</b><br>%{customdata}<br>log10(reads): %{z:.2f}<extra></extra>'
        : '<b>%{y}</b><br>%{customdata}<br>%{z:.2f}%%<extra></extra>'),
    }];

    // Cluster color bar: one row of cells (yaxis2) showing each column's
    // cluster id. Use a categorical Plotly colorscale built from the
    // shared paletteColor() so the bar's hues match every other view.
    const orderedClusterLabels = orderedSamples.map(s => clusterAssign[s.id] || 'C1');
    const uniqueClusters = [];
    const seenC = new Set();
    for (const c of orderedClusterLabels) if (!seenC.has(c)) { seenC.add(c); uniqueClusters.push(c); }
    const clusterIdx = new Map();
    uniqueClusters.forEach((c, i) => clusterIdx.set(c, i));
    const clusterZ = [orderedClusterLabels.map(c => clusterIdx.get(c))];
    // Build a discrete colorscale for the cluster bar. Each cluster gets a
    // [t0, t1] band of one color so adjacent cells aren't blended.
    const nC = Math.max(uniqueClusters.length, 1);
    const clusterColorscale = [];
    for (let i = 0; i < nC; i++) {
      const t0 = i / nC;
      const t1 = (i + 1) / nC;
      const color = paletteColor(i);
      clusterColorscale.push([t0, color]);
      clusterColorscale.push([t1, color]);
    }
    traces.push({
      type: 'heatmap',
      x: xVals,
      y: ['cluster'],
      z: clusterZ,
      xaxis: 'x',
      yaxis: 'y2',
      xgap: 0,
      ygap: 0,
      zmin: -0.5,
      zmax: nC - 0.5,
      colorscale: clusterColorscale,
      showscale: false,
      customdata: [orderedClusterLabels],
      hovertemplate: '%{customdata}<extra></extra>',
    });

    // Dendrogram traces. We project Ward-leaf positions onto the
    // displayed (orderedSamples) column index so the U-bars line up with
    // the heatmap cells regardless of the user's chosen Order mode. If
    // the displayed order doesn't match Ward's natural leaf order
    // (because the user picked a metadata sort), the dendrogram becomes
    // crossed — which is the honest visual signal that "your sort
    // disagrees with what Ward thinks is similar".
    if (wardResult.order.length >= 2) {
      const colByDisplay = new Map();
      orderedSamples.forEach((s, i) => colByDisplay.set(s.id, i));
      const wardLeafToCol = wardResult.order.map(sid => colByDisplay.get(sid) ?? -1);
      const xs = [];
      const ys = [];
      // For each merge, project both child positions through the ordered-
      // column mapping. dendroSegments uses Ward's leaf positions
      // directly; we need to remap each x value to its displayed column.
      const { linkage, order } = wardResult;
      const n = order.length;
      const nodeCol = new Map();  // merge node id → mean column of its leaves in displayed order
      order.forEach((sid, i) => nodeCol.set(i, colByDisplay.get(sid) ?? i));
      for (let m = 0; m < linkage.length; m++) {
        const merge = linkage[m];
        const leftLeaves = merge.left < n ? [order[merge.left]] : linkage[merge.left - n].leaves;
        const rightLeaves = merge.right < n ? [order[merge.right]] : linkage[merge.right - n].leaves;
        const lCols = leftLeaves.map(sid => colByDisplay.get(sid) ?? 0);
        const rCols = rightLeaves.map(sid => colByDisplay.get(sid) ?? 0);
        const lc = lCols.reduce((a, b) => a + b, 0) / lCols.length;
        const rc = rCols.reduce((a, b) => a + b, 0) / rCols.length;
        const lh = merge.left < n ? 0 : linkage[merge.left - n].distance;
        const rh = merge.right < n ? 0 : linkage[merge.right - n].distance;
        const h = merge.distance;
        // U-shape with `null` separators between segments so Plotly
        // doesn't connect adjacent merges.
        xs.push(lc, lc, rc, rc, null);
        ys.push(lh, h, h, rh, null);
        nodeCol.set(n + m, (lc + rc) / 2);
      }
      traces.push({
        type: 'scatter',
        mode: 'lines',
        x: xs,
        y: ys,
        xaxis: 'x',
        yaxis: 'y3',
        line: { color: '#64748b', width: 1 },
        hoverinfo: 'skip',
        showlegend: false,
      });
    }

    return { traces, xVals, xDisplay, yLabels };
  });

  let heatmapLayout = $derived.by(() => {
    if (!heatmapTrace) return {};
    const { xVals, xDisplay, yLabels } = heatmapTrace;
    // Layout: top 30% dendrogram → 4% cluster color bar → 65% heatmap,
    // separated by a 1% gap each. Sample labels live above the dendrogram
    // (xaxis with side: 'top' and matches: 'x'-driven).
    const heatmapDomain = [0, 0.65];
    const clusterDomain = [0.66, 0.70];
    const dendroDomain = [0.72, 1.0];
    return {
      height: Math.max(520, 160 + yLabels.length * 18),
      margin: { t: 110, r: 20, b: 30, l: 220 },
      xaxis: {
        side: 'top',
        tickmode: 'array',
        tickvals: xVals,
        ticktext: xDisplay,
        tickangle: -30,
        tickfont: { size: 10, color: '#cbd5e1' },
        automargin: true,
        showgrid: false,
        zeroline: false,
        range: [-0.5, xVals.length - 0.5],
        anchor: 'y3',  // labels sit at the top of the dendrogram subplot
      },
      yaxis: {
        type: 'category',
        tickmode: 'array',
        tickvals: yLabels,
        ticktext: yLabels,
        tickfont: { size: 10, color: '#cbd5e1' },
        automargin: true,
        autorange: 'reversed',
        domain: heatmapDomain,
        anchor: 'x',
      },
      yaxis2: {
        type: 'category',
        tickvals: ['cluster'],
        ticktext: [`k=${effectiveK}`],
        tickfont: { size: 9, color: '#94a3b8' },
        domain: clusterDomain,
        anchor: 'x',
        showgrid: false,
        zeroline: false,
      },
      yaxis3: {
        domain: dendroDomain,
        anchor: 'x',
        showgrid: false,
        zeroline: false,
        showticklabels: false,
        range: [0, dendroMaxHeight * 1.05],
        fixedrange: true,
      },
    };
  });
</script>

<div class="space-y-4">
  <!-- Controls -->
  <div class="flex items-center gap-3 flex-wrap text-xs">
    <button
      class="px-3 py-1 rounded-md border border-cyan-400 bg-cyan-400/10 text-cyan-400 hover:bg-cyan-400/20 transition-colors"
      style="min-width: 5.5rem"
      onclick={cycleRank}
      title={`Click to cycle rank: ${RANK_LEVELS.map(r => r.label).join(' → ')}`}
    >
      {currentRank.label} &#x25BE;
    </button>

    <button
      class="px-3 py-1 rounded-md border border-cyan-400 bg-cyan-400/10 text-cyan-400 hover:bg-cyan-400/20 transition-colors"
      style="min-width: 6rem"
      onclick={cycleValue}
      title="Click to cycle: % of classified ↔ log10(reads)"
    >
      {valueMode === 'pct' ? '% classified' : 'log10(reads)'} &#x25BE;
    </button>

    <label class="flex items-center gap-2 text-slate-400">
      <span>Order:</span>
      <select
        bind:value={orderIdx}
        class="px-2 py-1 rounded-md border border-cyan-400 bg-slate-900 text-cyan-400 hover:bg-cyan-400/10 transition-colors focus:outline-none focus:ring-1 focus:ring-cyan-400"
      >
        {#each orderModes as m, i}
          <option value={i}>{m.label}</option>
        {/each}
      </select>
    </label>

    <label class="flex items-center gap-2 text-slate-400">
      <span>Label:</span>
      <select
        bind:value={labelIdx}
        class="px-2 py-1 rounded-md border border-cyan-400 bg-slate-900 text-cyan-400 hover:bg-cyan-400/10 transition-colors focus:outline-none focus:ring-1 focus:ring-cyan-400"
      >
        {#each labelModes as m, i}
          <option value={i}>{m.label}</option>
        {/each}
      </select>
    </label>

    <div class="flex items-center gap-2 text-slate-400">
      <span>Top N</span>
      <input type="range" min="5" max="100" step="1" bind:value={topN}
        class="w-28 accent-cyan-400" />
      <span class="text-slate-500 w-10 font-mono tabular-nums">{topN}</span>
    </div>

    <!-- k= cluster cut. Range adapts to N (max 20). Color bar above the
         heatmap repaints whenever this slider moves; the cluster
         assignment is also pushed to the global sampleClusters store so
         /samples can color-by-cluster. -->
    <div class="flex items-center gap-2 text-slate-400">
      <span>k=</span>
      <input type="range" min="2" max={kMax} step="1" bind:value={kClusters}
        class="w-24 accent-cyan-400"
        title="Cut Ward dendrogram into k clusters" />
      <span class="text-slate-500 w-6 font-mono tabular-nums">{effectiveK}</span>
    </div>

    <span class="text-slate-500">
      {activeSamples.length} samples
      {#if $cartActive && $cartItems.size > 0}<span class="text-cyan-400">(cart filtered)</span>{/if}
    </span>
  </div>

  <!-- Heatmap -->
  <div class="rounded-lg border border-slate-700 bg-slate-900/40 p-2">
    {#if heatmapTrace}
      <PlotlyChart
        traces={heatmapTrace.traces}
        layout={heatmapLayout}
        exportName={`danaseq_heatmap_${currentRank.label.toLowerCase()}_${currentOrder.label.toLowerCase().replace(/[^a-z0-9]+/g, '-')}_by-${currentLabel.label.toLowerCase().replace(/[^a-z0-9]+/g, '-')}`}
      />
    {:else}
      <div class="p-8 text-center text-slate-500 text-sm">
        {#if !$sampleTaxonomy}
          Loading taxonomy data…
        {:else if !activeSamples.length}
          No samples in the current view.
        {:else}
          No classifications at {currentRank.label} rank for the active samples.
        {/if}
      </div>
    {/if}
  </div>

  <p class="text-[11px] text-slate-500 italic">
    Cells show each sample's {currentRank.label}-level composition
    {#if valueMode === 'log'}(log10 raw read counts){:else}(% of classified reads){/if}.
    Columns ordered by <code class="bg-slate-800 px-1 rounded">{currentOrder.label}</code>; taxa ranked by total count across the current sample set.
    Ward dendrogram cut at <code class="bg-slate-800 px-1 rounded">k={effectiveK}</code> clusters
    (color bar above heatmap; assignments shared with other views).
  </p>
</div>
