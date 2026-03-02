<script>
  import LeafletMap from '../components/charts/LeafletMap.svelte';
  import DataTable from '../components/ui/DataTable.svelte';
  import { samples, metadata } from '../stores/data.js';
  import { cartItems, cartActive, toggleCart } from '../stores/cart.js';
  import { selectedSample } from '../stores/selection.js';

  let colorBy = $state('dominant_phylum');

  const PALETTE = ['#22d3ee','#34d399','#fbbf24','#f87171','#a78bfa','#fb923c',
                   '#2dd4bf','#818cf8','#f472b6','#4ade80','#e879f9','#38bdf8'];

  // Build markers from samples + metadata
  let markers = $derived.by(() => {
    if (!$samples || !$metadata) return [];
    return $samples.map(s => {
      const m = $metadata[s.id];
      if (!m || m.lat == null || m.lon == null) return null;
      return {
        id: s.id,
        lat: m.lat,
        lon: m.lon,
        station: m.station,
        date: m.date,
        depth_m: m.depth_m,
        read_count: s.read_count,
        dominant_phylum: s.dominant_phylum,
        diversity: s.diversity,
      };
    }).filter(Boolean);
  });

  let markerColorMap = $derived.by(() => {
    const map = {};
    const values = [...new Set(markers.map(m => m[colorBy] || 'Unknown'))].filter(v => v !== 'Unknown').sort();
    values.forEach((v, i) => { map[v] = PALETTE[i % PALETTE.length]; });
    return map;
  });

  // Highlight cart items on map
  let highlightIds = $derived.by(() => {
    if ($cartItems.size === 0) return null;
    return $cartItems;
  });

  let noGeoData = $derived(!$metadata || markers.length === 0);

  const tableColumns = [
    { key: 'id', label: 'Sample' },
    { key: 'station', label: 'Station' },
    { key: 'lat', label: 'Lat', render: v => v?.toFixed(4) ?? '-' },
    { key: 'lon', label: 'Lon', render: v => v?.toFixed(4) ?? '-' },
    { key: 'depth_m', label: 'Depth (m)' },
    { key: 'date', label: 'Date' },
    { key: 'read_count', label: 'Reads', render: v => v?.toLocaleString() ?? '-' },
    { key: 'dominant_phylum', label: 'Dominant' },
  ];

  function handleMarkerClick(id) {
    selectedSample.set(id);
    toggleCart(id);
  }
</script>

<div class="space-y-6">
  {#if noGeoData}
    <div class="bg-amber-900/30 border border-amber-700 rounded-lg p-4 text-amber-300">
      <h3 class="font-semibold mb-1">No geographic metadata</h3>
      <p class="text-sm">Provide a metadata TSV file with <code class="bg-slate-800 px-1 rounded">lat</code> and <code class="bg-slate-800 px-1 rounded">lon</code> columns to enable the map view.</p>
      <p class="text-xs mt-2 text-amber-400">Re-run preprocess with: <code class="bg-slate-800 px-1 rounded">--metadata your_metadata.tsv</code></p>
    </div>
  {:else}
    <!-- Controls -->
    <div class="flex items-center gap-4">
      <label class="text-xs text-slate-400">
        Marker color:
        <select bind:value={colorBy} class="ml-1 bg-slate-800 border border-slate-600 rounded px-2 py-1 text-xs text-slate-200">
          <option value="dominant_phylum">Dominant Phylum</option>
          <option value="station">Station</option>
        </select>
      </label>
      <span class="text-xs text-slate-500">{markers.length} samples with GPS coordinates</span>
    </div>

    <!-- Map -->
    <div class="h-[500px] rounded-lg overflow-hidden border border-slate-700">
      <LeafletMap
        {markers}
        colorMap={markerColorMap}
        onMarkerClick={handleMarkerClick}
        {highlightIds}
      />
    </div>

    <!-- Location table -->
    <DataTable
      columns={tableColumns}
      rows={markers}
      onRowClick={(row) => { selectedSample.set(row.id); toggleCart(row.id); }}
      selectedId={$selectedSample}
      idKey="id"
      maxHeight="300px"
      exportFilename="sample_locations"
    />
  {/if}
</div>
