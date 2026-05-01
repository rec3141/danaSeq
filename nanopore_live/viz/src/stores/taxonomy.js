/**
 * Taxonomy drill-down helpers for MapView.
 *
 * Works against the sample_taxonomy.json produced by viz/preprocess/preprocess.py:
 *   {
 *     parents: { taxon: parent_taxon, ... },          // full lineage map
 *     ranks:   { taxon: 'K'|'P'|'C'|'O'|'F'|'G'|'S'|'R2', ... },
 *     samples: { sample_id: { phylum: {...}, class: {...}, direct: {...} }, ... }
 *   }
 *
 * The `direct` dict holds per-sample counts at whatever rank each leaf taxon sits
 * at (filtered to counts ≥ 7 by preprocess). We walk up the parent chain from
 * each leaf to aggregate counts to any requested rank.
 */
import { writable, derived, get } from 'svelte/store';
import { taxonomySource } from './taxonomySource.js';
import { sampleTaxonomy } from './data.js';

// Rank codes in top-down order.
//
// Both classifiers now start at Domain:
//   - GTDB uses 'D' (its native prefix).
//   - Kraken uses 'R2' — Kraken's superkingdom rank that holds Bacteria /
//     Archaea / Eukaryota / Viruses. Kraken's per-read `*_domain` field is
//     populated from R2 (see compute_read_tsne.build_lineage_lookup), so
//     aligning the nav with R2 lets Reads filter/color correctly. Kraken's
//     'K' (Kingdom — e.g. Pseudomonadati) would not line up with any
//     per-read field and is intentionally skipped from the drill chain.
export const RANK_ORDER_KRAKEN = ['R2', 'P', 'C', 'O', 'F', 'G', 'S'];
export const RANK_ORDER_GTDB   = ['D',  'P', 'C', 'O', 'F', 'G', 'S'];

// Reactive rank order — flips with the global taxonomy source toggle.
// .svelte files should read $rankOrder; plain-JS helpers fall back to get().
export const rankOrder = derived(
  taxonomySource,
  ($src) => $src === 'gtdb' ? RANK_ORDER_GTDB : RANK_ORDER_KRAKEN,
);

// Legacy export kept so existing `import { RANK_ORDER }` call sites compile,
// but it's a snapshot of whatever source was active at module load. Prefer
// $rankOrder (Svelte) or get(rankOrder) (plain JS) for anything reactive.
export const RANK_ORDER = get(rankOrder);

export const RANK_LABELS = {
  D: 'Domain',
  K: 'Kingdom',
  P: 'Phylum',
  C: 'Class',
  O: 'Order',
  F: 'Family',
  G: 'Genus',
  S: 'Species',
  R2: 'Root',
};

export const RANK_LABELS_PLURAL = {
  D: 'Domains',
  K: 'Kingdoms',
  P: 'Phyla',
  C: 'Classes',
  O: 'Orders',
  F: 'Families',
  G: 'Genera',
  S: 'Species',
  R2: 'Roots',
};

// High-contrast palette reused from MapView. Cycles if there are more taxa
// than colors — rare at Kingdom/Phylum, possible at deeper ranks.
const PALETTE = [
  '#22d3ee', '#34d399', '#fbbf24', '#f87171', '#a78bfa', '#fb923c',
  '#2dd4bf', '#818cf8', '#f472b6', '#4ade80', '#e879f9', '#38bdf8',
  '#fde047', '#f43f5e', '#a3e635', '#6366f1', '#ec4899', '#14b8a6',
  '#eab308', '#06b6d4', '#d946ef', '#84cc16', '#facc15', '#fb7185',
];

export function paletteColor(i) {
  return PALETTE[((i % PALETTE.length) + PALETTE.length) % PALETTE.length];
}

// Get the next rank below `rank`, or null if already at the leaf (species).
// Callers in reactive contexts should pass the current $rankOrder explicitly;
// plain-JS callers fall through to the source-aware default.
export function nextRank(rank, order = get(rankOrder)) {
  const i = order.indexOf(rank);
  if (i < 0 || i === order.length - 1) return null;
  return order[i + 1];
}

// Get the rank above `rank`, or null if already at the top (Kingdom/Domain).
export function prevRank(rank, order = get(rankOrder)) {
  const i = order.indexOf(rank);
  if (i <= 0) return null;
  return order[i - 1];
}

// Walk up the parent chain and collect the nearest ancestor at each rank.
// Memoized per (parents, ranks) identity so repeated queries are cheap.
let _cache = new WeakMap();
export function ancestorInfo(taxon, parents, ranks) {
  let perParents = _cache.get(parents);
  if (!perParents) { perParents = new Map(); _cache.set(parents, perParents); }
  if (perParents.has(taxon)) return perParents.get(taxon);

  const info = {};
  let cur = taxon;
  const seen = new Set();
  while (cur && !seen.has(cur)) {
    seen.add(cur);
    const r = ranks[cur];
    // Prefer the nearest-to-leaf ancestor at each rank (first hit wins).
    if (r && !(r in info)) info[r] = cur;
    cur = parents[cur];
  }
  perParents.set(taxon, info);
  return info;
}

// True if `ancestor` appears in the lineage of `descendant` (inclusive).
export function isAncestorOf(ancestor, descendant, parents) {
  let cur = descendant;
  const seen = new Set();
  while (cur && !seen.has(cur)) {
    if (cur === ancestor) return true;
    seen.add(cur);
    cur = parents[cur];
  }
  return false;
}

/**
 * Aggregate a sample's `direct` counts to the given `targetRank`, keeping only
 * lineages that pass through `filterTaxon` (if provided).
 *
 * Returns: { taxonAtTargetRank: count, ... } — omits taxa with zero count.
 */
export function aggregateAtRank(sampleDirect, parents, ranks, targetRank, filterTaxon) {
  const out = {};
  if (!sampleDirect) return out;
  for (const leaf in sampleDirect) {
    const cnt = sampleDirect[leaf];
    if (!cnt) continue;
    // If a filter is active, leaf must descend from it.
    if (filterTaxon && !isAncestorOf(filterTaxon, leaf, parents)) continue;
    const info = ancestorInfo(leaf, parents, ranks);
    const atTarget = info[targetRank];
    if (!atTarget) continue;  // leaf isn't classified deep enough
    out[atTarget] = (out[atTarget] || 0) + cnt;
  }
  return out;
}

/**
 * Build { colorMap, ranked } for the union of sub-taxa across all samples at
 * (targetRank, filterTaxon). Deterministic ordering by global count descending
 * so the same taxon gets the same color across re-renders.
 */
export function buildSubTaxaIndex(samples, parents, ranks, targetRank, filterTaxon) {
  const globalCounts = {};
  const perSample = {};
  for (const sid in samples) {
    const direct = samples[sid]?.direct;
    const agg = aggregateAtRank(direct, parents, ranks, targetRank, filterTaxon);
    perSample[sid] = agg;
    for (const t in agg) globalCounts[t] = (globalCounts[t] || 0) + agg[t];
  }
  const ranked = Object.entries(globalCounts)
    .map(([name, count]) => ({ name, count }))
    .sort((a, b) => b.count - a.count || a.name.localeCompare(b.name));
  const colorMap = {};
  ranked.forEach((r, i) => {
    r.color = paletteColor(i);
    colorMap[r.name] = r.color;
  });
  return { colorMap, ranked, perSample };
}

/**
 * Nav state store for the map's taxonomy drill-down. Persisted across tab
 * switches within a session; reset on page load.
 *
 *   level:    currently-displayed rank (e.g. 'P')
 *   filter:   the taxon whose descendants are shown (or null = all)
 *   stack:    breadcrumb trail; each entry is { level, filter, isolated }
 *   isolated: if set, only this taxon's ring is drawn per sample (across-sample
 *             comparison mode). Set by clicking at the leaf rank where there's
 *             nothing below to drill into.
 */
function emptyNavFor(source) {
  const order = source === 'gtdb' ? RANK_ORDER_GTDB : RANK_ORDER_KRAKEN;
  return { level: order[0], filter: null, isolated: null, stack: [] };
}

function makeTaxNav() {
  const { subscribe, update, set } = writable(emptyNavFor(get(taxonomySource)));

  // Reset on classifier flip — rank sets differ (Kraken top is R2, GTDB top
  // is D), so stale level/filter/stack state would produce empty views with
  // no matching taxa. Read the new source directly out of the callback
  // rather than going through get(rankOrder), which can be stale on the
  // first firing of a chained derived.
  taxonomySource.subscribe((src) => set(emptyNavFor(src)));

  return {
    subscribe,
    // Click action: drill into the next rank if one exists, otherwise (at leaf
    // rank) isolate this taxon so only its ring is rendered on the map.
    drillInto(taxon) {
      update(s => {
        const next = nextRank(s.level);
        if (!next) {
          // At leaf rank — toggle isolation for cross-sample comparison.
          return { ...s, isolated: s.isolated === taxon ? null : taxon };
        }
        return {
          level: next,
          filter: taxon,
          isolated: null,
          stack: [...s.stack, { level: s.level, filter: s.filter, isolated: s.isolated }],
        };
      });
    },
    // Go up one step: first clear isolation if any, then pop the stack.
    up() {
      update(s => {
        if (s.isolated) return { ...s, isolated: null };
        if (!s.stack.length) return emptyNav();
        const prev = s.stack[s.stack.length - 1];
        return {
          level: prev.level,
          filter: prev.filter,
          isolated: prev.isolated || null,
          stack: s.stack.slice(0, -1),
        };
      });
    },
    reset() {
      set(emptyNavFor(get(taxonomySource)));
    },
    // Jump to a specific rank while preserving the current filter. Used by
    // the cycle-rank button on /reads so clicking through Phylum → Class →
    // Order … overrides the drill-down picker without wiping the context
    // the user already drilled into.
    setLevel(level) {
      update(s => (s.level === level ? s : { ...s, level, isolated: null, stack: [] }));
    },
    // Jump straight to a picked taxon (used by the search autocomplete).
    // Replaces the breadcrumb stack entirely with the picked taxon's ancestor
    // chain so that "Up" walks back through its lineage (Species → Genus →
    // Family → ... → Kingdom) rather than back to wherever the user happened
    // to be drilled when they ran the search.
    //
    // ancestorStack must contain one entry per drillable rank strictly above
    // the post-jump `level`, each `{ level, filter, isolated: null }`. The
    // caller is responsible for building it from the picked taxon's lineage.
    jumpTo(taxonRank, taxonName, parentForLeaf, ancestorStack) {
      update(_s => {
        const next = nextRank(taxonRank);
        if (!next) {
          return {
            level: taxonRank,
            filter: parentForLeaf,
            isolated: taxonName,
            stack: ancestorStack,
          };
        }
        return {
          level: next,
          filter: taxonName,
          isolated: null,
          stack: ancestorStack,
        };
      });
    },
  };
}

export const taxNav = makeTaxNav();

// Shared derived state consumed by TaxonomyDrillNav.svelte and every view that
// colors/sizes by sub-taxon (Map, Reads, Samples). Both views read the same
// index so their colors stay consistent.
export const activeSubTaxa = derived(
  [sampleTaxonomy, taxNav],
  ([$tax, $nav]) => {
    if (!$tax?.samples || !$tax?.parents || !$tax?.ranks) return null;
    return buildSubTaxaIndex($tax.samples, $tax.parents, $tax.ranks, $nav.level, $nav.filter);
  },
);

// Full-tree search candidates for the taxonomy nav's autocomplete: every
// drillable-rank taxon with a non-zero aggregated count across samples.
// Empty while the taxonomy isn't loaded.
export const taxSearchCandidates = derived(
  [sampleTaxonomy, rankOrder],
  ([$tax, $order]) => {
    if (!$tax?.parents || !$tax?.ranks || !$tax?.samples) return [];
    const drillable = new Set($order);
    const counts = {};
    for (const sid in $tax.samples) {
      const direct = $tax.samples[sid]?.direct || {};
      for (const leaf in direct) {
        const cnt = direct[leaf];
        if (!cnt) continue;
        let cur = leaf;
        const seen = new Set();
        while (cur && !seen.has(cur)) {
          seen.add(cur);
          const r = $tax.ranks[cur];
          if (r && drillable.has(r)) counts[cur] = (counts[cur] || 0) + cnt;
          cur = $tax.parents[cur];
        }
      }
    }
    const out = [];
    for (const name in counts) {
      if (counts[name] <= 0) continue;
      const chain = [];
      let cur = name;
      const seen = new Set();
      while (cur && !seen.has(cur)) {
        seen.add(cur);
        chain.unshift(cur);
        cur = $tax.parents[cur];
      }
      out.push({ name, rank: $tax.ranks[name], lineage: chain.join(';'), count: counts[name] });
    }
    out.sort((a, b) => b.count - a.count || a.name.localeCompare(b.name));
    return out;
  },
);

// Jump the drill-down to a search-picked taxon. Reconstructs the breadcrumb
// stack from the taxon's lineage so "Up" walks back through its ancestors.
export function jumpToTaxon(candidate) {
  const tax = get(sampleTaxonomy);
  const order = get(rankOrder);
  if (!candidate || !tax?.parents || !tax?.ranks) return;
  const info = ancestorInfo(candidate.name, tax.parents, tax.ranks);
  const pRank = prevRank(candidate.rank, order);
  const parentForLeaf = pRank ? (info[pRank] || null) : null;
  const postJumpLevel = nextRank(candidate.rank, order) || candidate.rank;
  const stack = [];
  for (const r of order) {
    if (r === postJumpLevel) break;
    const prR = prevRank(r, order);
    stack.push({ level: r, filter: prR ? (info[prR] || null) : null, isolated: null });
  }
  taxNav.jumpTo(candidate.rank, candidate.name, parentForLeaf, stack);
}
