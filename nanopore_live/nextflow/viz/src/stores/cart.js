import { writable, derived } from 'svelte/store';

// Set of sample keys ("flowcell/barcode")
export const cartItems = writable(new Set());

// When true, all views filter to only cart items
export const cartActive = writable(false);

export function addToCart(sampleId) {
  cartItems.update(s => { s.add(sampleId); return new Set(s); });
}

export function removeFromCart(sampleId) {
  cartItems.update(s => { s.delete(sampleId); return new Set(s); });
}

export function toggleCart(sampleId) {
  cartItems.update(s => {
    if (s.has(sampleId)) s.delete(sampleId);
    else s.add(sampleId);
    return new Set(s);
  });
}

export function clearCart() {
  cartItems.set(new Set());
}

export const cartCount = derived(cartItems, $c => $c.size);
