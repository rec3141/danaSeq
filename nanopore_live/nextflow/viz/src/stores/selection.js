import { writable } from 'svelte/store';

// Currently focused sample key (drives detail panels)
export const selectedSample = writable(null);

// Currently focused read ID
export const selectedRead = writable(null);
