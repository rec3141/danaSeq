/**
 * Sample cluster assignments — published by HeatmapView when a Ward
 * clustering runs and a k= cut is applied. Other views (e.g. /samples)
 * can read this store to color samples by their heatmap cluster.
 *
 * Shape: { sampleId: 'C1' | 'C2' | ... }
 *   - empty object when no clustering has been computed (e.g. heatmap
 *     hasn't been visited yet, or fewer than 2 active samples).
 *   - labels are 1-indexed in left-to-right dendrogram order so adjacent
 *     clusters in the heatmap also sit next to each other in the legend.
 */
import { writable } from 'svelte/store';

export const sampleClusters = writable({});

// Total number of clusters in the most recent assignment (helpful for UIs
// that want to size a categorical legend or pick a palette range).
export const sampleClusterK = writable(0);

// Shared color map for cluster labels, e.g. { 'C1': '#22d3ee', 'C2': '#34d399' }.
// Published by HeatmapView so /samples and /map paint the same color per cluster.
// Empty until a clustering has been computed.
export const sampleClusterColors = writable({});
