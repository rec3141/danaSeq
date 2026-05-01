<script>
  import {
    cartItems, cartActive, removeFromCart, clearCart,
    readSelections, activeReadSelection, removeReadSelection, clearReadSelections,
  } from '../../stores/cart.js';
  import { samples } from '../../stores/data.js';

  let { open = false } = $props();

  function getSampleLabel(id) {
    const s = $samples?.find(s => s.id === id);
    if (s?.station) return `${id} (${s.station})`;
    return id;
  }

  function toggleActive(name) {
    activeReadSelection.update(cur => cur === name ? null : name);
  }
</script>

{#if open}
  <div class="fixed right-0 top-14 bottom-0 w-80 bg-slate-900 border-l border-slate-700 z-40 flex flex-col shadow-2xl">
    <!-- Samples section -->
    <div class="p-4 border-b border-slate-700 flex items-center justify-between">
      <h3 class="text-sm font-semibold text-slate-200">
        Sample Cart ({$cartItems.size})
      </h3>
      {#if $cartItems.size > 0}
        <button
          class="text-xs px-2 py-1 rounded text-rose-400 hover:text-rose-300 border border-slate-600 hover:border-rose-400 transition-colors"
          onclick={clearCart}
        >
          Clear
        </button>
      {/if}
    </div>
    <div class="flex-1 overflow-y-auto flex flex-col">
      <div class="p-2">
        {#if $cartItems.size === 0}
          <p class="text-slate-500 text-xs p-2">No samples in cart. Click samples in the explorer or map to add them.</p>
        {:else}
          <ul class="space-y-1">
            {#each [...$cartItems] as id}
              <li class="flex items-center justify-between px-3 py-2 rounded bg-slate-800 border border-slate-700 group">
                <span class="text-xs text-slate-300 font-mono truncate">{getSampleLabel(id)}</span>
                <button
                  class="text-slate-500 hover:text-rose-400 opacity-0 group-hover:opacity-100 transition-all ml-2"
                  onclick={() => removeFromCart(id)}
                >
                  <svg class="w-3.5 h-3.5" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                    <path stroke-linecap="round" stroke-linejoin="round" stroke-width="2" d="M6 18L18 6M6 6l12 12"/>
                  </svg>
                </button>
              </li>
            {/each}
          </ul>
        {/if}
      </div>

      <!-- Read selections section -->
      <div class="border-t border-slate-700 mt-2">
        <div class="p-4 flex items-center justify-between">
          <h3 class="text-sm font-semibold text-slate-200">
            Read Selections ({$readSelections.size})
          </h3>
          {#if $readSelections.size > 0}
            <button
              class="text-xs px-2 py-1 rounded text-rose-400 hover:text-rose-300 border border-slate-600 hover:border-rose-400 transition-colors"
              onclick={clearReadSelections}
            >
              Clear
            </button>
          {/if}
        </div>
        <div class="px-2 pb-3">
          {#if $readSelections.size === 0}
            <p class="text-slate-500 text-xs px-2">No read selections. Lasso reads on the Reads view and click "+ Cart" to save one.</p>
          {:else}
            <ul class="space-y-1">
              {#each [...$readSelections.entries()] as [name, sel]}
                {@const isActive = $activeReadSelection === name}
                <li class="px-3 py-2 rounded bg-slate-800 border group transition-colors
                  {isActive ? 'border-amber-400/60' : 'border-slate-700'}">
                  <div class="flex items-center justify-between gap-2">
                    <button
                      class="flex items-center gap-2 flex-1 min-w-0 text-left"
                      onclick={() => toggleActive(name)}
                      title={isActive ? 'Click to deactivate (no map overlay)' : 'Click to visualize on /map'}
                    >
                      <span class="inline-block w-3 h-3 rounded-full border flex-shrink-0
                        {isActive ? 'bg-amber-400 border-amber-300' : 'border-slate-500'}"></span>
                      <span class="text-xs font-mono truncate
                        {isActive ? 'text-amber-300' : 'text-slate-300'}">{name}</span>
                    </button>
                    <button
                      class="text-slate-500 hover:text-rose-400 opacity-0 group-hover:opacity-100 transition-all flex-shrink-0"
                      onclick={() => removeReadSelection(name)}
                      title="Remove this selection"
                    >
                      <svg class="w-3.5 h-3.5" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                        <path stroke-linecap="round" stroke-linejoin="round" stroke-width="2" d="M6 18L18 6M6 6l12 12"/>
                      </svg>
                    </button>
                  </div>
                  <div class="text-[10px] text-slate-500 mt-1 ml-5 font-mono">
                    {sel.readIds.size.toLocaleString()} reads · {sel.samples.size} samples
                  </div>
                </li>
              {/each}
            </ul>
          {/if}
        </div>
      </div>
    </div>
  </div>
{/if}
