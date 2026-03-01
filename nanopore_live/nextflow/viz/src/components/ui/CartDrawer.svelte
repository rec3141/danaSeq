<script>
  import { cartItems, cartActive, removeFromCart, clearCart } from '../../stores/cart.js';
  import { samples } from '../../stores/data.js';

  let { open = false } = $props();

  function getSampleLabel(id) {
    const s = $samples?.find(s => s.id === id);
    if (s?.station) return `${id} (${s.station})`;
    return id;
  }
</script>

{#if open}
  <div class="fixed right-0 top-14 bottom-0 w-80 bg-slate-900 border-l border-slate-700 z-40 flex flex-col shadow-2xl">
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
    <div class="flex-1 overflow-y-auto p-2">
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
  </div>
{/if}
