import { writable } from 'svelte/store';

// Currently selected MAG name (drives detail panels)
export const selectedMag = writable(null);
