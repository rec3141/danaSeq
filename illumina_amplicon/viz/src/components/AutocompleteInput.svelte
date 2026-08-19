<script>
  let {
    value = $bindable(''),
    placeholder = '',
    label = '',
    candidates = [],
    maxSuggestions = 8,
    onPick = null,
    onType = null,   // fired on manual edits only — pick() sets the value directly
    // Extra strings to search for besides what was typed, so a name that has
    // been retired can still lead to the one that replaced it. Returns the
    // query first; anything after it is a synonym and ranks below a direct hit.
    expand = null,
  } = $props();

  let focused = $state(false);
  let selectedIdx = $state(-1);
  let inputEl;

  // Ranked, because the cut to maxSuggestions decides what the reader believes
  // exists. Alphabetical order put Alphaproteobacteria and Gammaproteobacteria
  // above the phylum a search for "Proteo" is actually after, and the phylum
  // fell off the end — the list then reads as proof it is not in the data.
  // Rank: what the query starts, then what merely contains it, then names
  // reached through a synonym; alphabetical within each.
  let suggestions = $derived.by(() => {
    if (!value || !focused) return [];
    const needles = expand ? expand(value) : [value.toLowerCase()];
    if (!needles.length) return [];
    const [query, ...synonyms] = needles;

    const rank = (c) => {
      const lc = c.toLowerCase();
      if (lc.startsWith(query)) return 0;
      if (lc.includes(query)) return 1;
      if (synonyms.some(n => lc.startsWith(n))) return 2;
      if (synonyms.some(n => lc.includes(n))) return 3;
      return -1;
    };

    return candidates
      .map(c => ({ c, r: rank(c) }))
      .filter(x => x.r >= 0)
      .sort((a, b) => a.r - b.r || a.c.localeCompare(b.c))
      .slice(0, maxSuggestions)
      .map(x => x.c);
  });

  function pick(s) {
    value = s;
    focused = false;
    selectedIdx = -1;
    if (onPick) onPick(s);
  }

  function handleKeydown(e) {
    if (!suggestions.length) return;
    if (e.key === 'ArrowDown') {
      e.preventDefault();
      selectedIdx = Math.min(selectedIdx + 1, suggestions.length - 1);
    } else if (e.key === 'ArrowUp') {
      e.preventDefault();
      selectedIdx = Math.max(selectedIdx - 1, -1);
    } else if (e.key === 'Enter' && selectedIdx >= 0) {
      e.preventDefault();
      pick(suggestions[selectedIdx]);
    } else if (e.key === 'Escape') {
      focused = false;
    } else {
      selectedIdx = -1;
    }
  }
</script>

<label class="block relative">
  {#if label}
    <span class="text-xs text-muted">{label}</span>
  {/if}
  <input
    type="text"
    bind:value
    bind:this={inputEl}
    {placeholder}
    oninput={() => onType?.()}
    onfocus={() => focused = true}
    onblur={() => setTimeout(() => focused = false, 150)}
    onkeydown={handleKeydown}
    class="mt-1 w-full rounded border border-line bg-raised px-2 py-1 text-sm text-fg placeholder-faint focus:border-blue-500 focus:outline-none"
  />
  {#if suggestions.length > 0 && focused}
    <ul class="absolute z-50 mt-1 w-full rounded border border-line bg-raised shadow-lg max-h-48 overflow-y-auto">
      {#each suggestions as s, i}
        <li>
          <button
            class="w-full text-left px-2 py-1 text-xs truncate {i === selectedIdx ? 'bg-blue-600 text-white' : 'text-fg2 hover:bg-raised2'}"
            onmousedown={() => pick(s)}
          >
            {s}
          </button>
        </li>
      {/each}
    </ul>
  {/if}
</label>
