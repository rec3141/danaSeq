<script>
  /**
   * Shared taxonomy drill-down sidebar used across Map / Reads / Samples.
   *
   * Reads its state from global stores (taxNav + activeSubTaxa + taxSearchCandidates)
   * so every view showing this component operates on the same drill-down cursor —
   * flipping tabs preserves where you are in the tree.
   *
   * Visual language mirrors microscape's AutocompleteInput + a scrollable ranked
   * sub-taxa list; clicking an entry drills into it, clicking at a leaf rank
   * toggles isolation.
   */
  import {
    taxNav, activeSubTaxa, taxSearchCandidates, jumpToTaxon,
    RANK_LABELS, RANK_LABELS_PLURAL, nextRank, rankOrder,
  } from '../../stores/taxonomy.js';

  let { heightClass = 'h-[500px]' } = $props();

  let taxSearch = $state('');
  let taxSearchFocused = $state(false);
  let taxSearchSelectedIdx = $state(-1);

  // Reset search whenever the drill-down cursor moves — stale queries confuse.
  $effect(() => { [$taxNav.level, $taxNav.filter, $taxNav.isolated]; taxSearch = ''; taxSearchSelectedIdx = -1; });

  let taxSuggestions = $derived.by(() => {
    if (!taxSearch.trim() || !taxSearchFocused) return [];
    let re = null;
    try { re = new RegExp(taxSearch, 'i'); } catch { /* invalid regex → substring */ }
    const q = taxSearch.toLowerCase();
    const out = [];
    for (const c of $taxSearchCandidates) {
      const hit = re ? re.test(c.name) : c.name.toLowerCase().includes(q);
      if (hit) {
        out.push(c);
        if (out.length >= 20) break;
      }
    }
    return out;
  });

  function pick(c) {
    jumpToTaxon(c);
    taxSearch = '';
    taxSearchFocused = false;
    taxSearchSelectedIdx = -1;
  }

  function onKeydown(e) {
    if (e.key === 'Escape') { taxSearchFocused = false; return; }
    if (!taxSuggestions.length) return;
    if (e.key === 'ArrowDown') {
      e.preventDefault();
      taxSearchSelectedIdx = Math.min(taxSearchSelectedIdx + 1, taxSuggestions.length - 1);
    } else if (e.key === 'ArrowUp') {
      e.preventDefault();
      taxSearchSelectedIdx = Math.max(taxSearchSelectedIdx - 1, -1);
    } else if (e.key === 'Enter' && taxSearchSelectedIdx >= 0) {
      e.preventDefault();
      pick(taxSuggestions[taxSearchSelectedIdx]);
    } else {
      taxSearchSelectedIdx = -1;
    }
  }

  let isLeaf = $derived(!nextRank($taxNav.level, $rankOrder));
</script>

{#if $activeSubTaxa}
<aside class="w-56 shrink-0 {heightClass} flex flex-col bg-slate-900/60 border border-slate-700 rounded-lg">
  <div class="px-3 py-2 border-b border-slate-700 bg-slate-900/80 rounded-t-lg">
    <div class="flex items-center justify-between">
      <p class="text-xs font-semibold uppercase tracking-wider text-slate-300">
        {RANK_LABELS[$taxNav.level] ?? $taxNav.level}
        <span class="text-slate-500 font-normal">({$activeSubTaxa.ranked.length})</span>
      </p>
      {#if $taxNav.filter || $taxNav.stack.length > 0 || $taxNav.isolated}
        <button
          class="text-[10px] text-slate-500 hover:text-slate-300"
          onclick={() => taxNav.reset()}
          title="Reset to top rank"
        >reset</button>
      {/if}
    </div>
    <div class="relative mt-1.5">
      <input
        type="text"
        bind:value={taxSearch}
        onfocus={() => taxSearchFocused = true}
        onblur={() => setTimeout(() => taxSearchFocused = false, 150)}
        onkeydown={onKeydown}
        placeholder="search full taxonomy…"
        class="w-full bg-slate-800 border border-slate-700 rounded px-2 py-1 pr-6 text-[11px] text-slate-200 placeholder-slate-500 focus:outline-none focus:border-cyan-500/50"
      />
      {#if taxSearch}
        <button
          class="absolute right-1 top-1/2 -translate-y-1/2 text-slate-500 hover:text-slate-300 px-1"
          onmousedown={() => { taxSearch = ''; taxSearchSelectedIdx = -1; }}
          title="Clear search"
        >&times;</button>
      {/if}
      {#if taxSuggestions.length > 0 && taxSearchFocused}
        <ul class="absolute z-[9999] left-0 right-0 top-full mt-1 bg-slate-800 border border-slate-700 rounded shadow-lg max-h-60 overflow-y-auto">
          {#each taxSuggestions as s, i}
            <li>
              <button
                class="w-full text-left px-2 py-1 text-[11px] flex items-baseline gap-1.5 {i === taxSearchSelectedIdx ? 'bg-cyan-500/20 text-cyan-200' : 'text-slate-300 hover:bg-slate-700'}"
                onmousedown={() => pick(s)}
                title={s.lineage}
              >
                <span class="truncate">{s.name}</span>
                <span class="ml-auto shrink-0 text-[9px] uppercase tracking-wider text-slate-500">{RANK_LABELS[s.rank] ?? s.rank}</span>
              </button>
            </li>
          {/each}
        </ul>
      {:else if taxSearch.trim() && taxSearchFocused}
        <div class="absolute z-[9999] left-0 right-0 top-full mt-1 bg-slate-800 border border-slate-700 rounded px-2 py-1 text-[11px] text-slate-500 italic">
          No taxa match "{taxSearch}"
        </div>
      {/if}
    </div>
  </div>
  <div class="flex-1 overflow-y-auto p-1 space-y-0.5">
    {#if $taxNav.isolated}
      <button
        class="flex items-center gap-1.5 w-full text-left text-xs hover:bg-slate-800 rounded px-2 py-1 text-cyan-400 border-b border-slate-700 mb-1"
        onclick={() => taxNav.up()}
      >
        &#x25B4; Show all {$taxNav.filter ? `${$taxNav.filter} ${(RANK_LABELS_PLURAL[$taxNav.level] ?? $taxNav.level).toLowerCase()}` : (RANK_LABELS_PLURAL[$taxNav.level] ?? $taxNav.level)}
      </button>
    {:else if $taxNav.stack.length > 0}
      {@const prev = $taxNav.stack[$taxNav.stack.length - 1]}
      <button
        class="flex items-center gap-1.5 w-full text-left text-xs hover:bg-slate-800 rounded px-2 py-1 text-cyan-400 border-b border-slate-700 mb-1"
        onclick={() => taxNav.up()}
      >
        &#x25B4; Up to {prev.filter ?? (RANK_LABELS[prev.level] ?? prev.level)}
      </button>
    {/if}
    {#if $activeSubTaxa.ranked.length === 0}
      <p class="text-[11px] text-slate-500 px-2 py-1 italic">No classifications at this rank.</p>
    {/if}
    {#each $activeSubTaxa.ranked as item}
      {@const isIsolated = $taxNav.isolated === item.name}
      <button
        class="flex items-center gap-1.5 w-full text-left text-xs rounded px-2 py-1 group transition-colors
          {isIsolated ? 'bg-amber-400/10 text-amber-300 hover:bg-amber-400/20' : 'hover:bg-slate-800 text-slate-200'}"
        onclick={() => taxNav.drillInto(item.name)}
        title={isLeaf
          ? (isIsolated ? `Clear isolation of ${item.name}` : `Isolate ${item.name}`)
          : `Drill into ${item.name}`}
      >
        <span class="inline-block h-2.5 w-2.5 rounded-full shrink-0 border border-slate-900" style="background:{item.color}"></span>
        <span class="truncate group-hover:text-slate-100">{item.name}</span>
        <span class="ml-auto shrink-0 font-mono text-[10px] text-slate-500">{item.count.toLocaleString()}</span>
      </button>
    {/each}
  </div>
</aside>
{/if}
