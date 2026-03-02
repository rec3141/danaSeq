<script>
  import { onMount } from 'svelte';

  let { markers = [], colorMap = {}, colorBy = 'dominant_phylum', sizeBy = 'fixed', sizeScale = 1.0, nudgeMeters = 100, onMarkerClick = null, highlightIds = null } = $props();
  let mapContainer;
  let map = null;
  let markerLayer = null;
  let hasFitted = false;

  const PALETTE = ['#22d3ee','#34d399','#fbbf24','#f87171','#a78bfa','#fb923c',
                   '#2dd4bf','#818cf8','#f472b6','#4ade80','#e879f9','#38bdf8'];

  function getColor(marker) {
    // Check color map by the active colorBy field value
    const key = marker[colorBy] || marker.dominant_phylum || marker.station || marker.id;
    if (colorMap[key]) return colorMap[key];
    if (colorMap[marker.id]) return colorMap[marker.id];
    if (marker.color) return marker.color;
    return PALETTE[Math.abs(hashCode(String(key))) % PALETTE.length];
  }

  function hashCode(s) {
    let hash = 0;
    for (let i = 0; i < s.length; i++) hash = ((hash << 5) - hash) + s.charCodeAt(i);
    return hash;
  }

  function markerRadius(m) {
    if (sizeBy === 'fixed') return 8 * sizeScale;
    const val = typeof m[sizeBy] === 'number' ? m[sizeBy] : 0;
    if (val === 0) return 3 * sizeScale;
    // Log-based sizing: good dynamic range from small (0-10) to large (1e6+)
    const r = Math.max(3, Math.min(30, 3 + Math.log1p(val) * 1.8));
    return r * sizeScale;
  }

  function tooltipContent(m) {
    const parts = [`<b>${m.id}</b>`];
    if (m.station) parts.push(`Station: ${m.station}`);
    if (m.date) parts.push(`Date: ${m.date}`);
    if (m.depth_m != null) parts.push(`Depth: ${m.depth_m}m`);
    if (m.read_count != null) parts.push(`Reads: ${m.read_count.toLocaleString()}`);
    if (m.dominant_phylum) parts.push(`Dominant: ${m.dominant_phylum}`);
    if (m.diversity != null) parts.push(`Shannon H: ${m.diversity.toFixed(2)}`);
    return `<div style="color:#e2e8f0;font-family:Inter,sans-serif;font-size:11px;line-height:1.4">${parts.join('<br>')}</div>`;
  }

  async function initMap() {
    const L = await import('leaflet');
    await import('leaflet/dist/leaflet.css');

    map = L.map(mapContainer, {
      center: [0, 0], zoom: 2,
      zoomControl: true,
      attributionControl: true,
    });

    L.tileLayer('https://{s}.tile.openstreetmap.org/{z}/{x}/{y}.png', {
      attribution: '&copy; OpenStreetMap',
      maxZoom: 18,
    }).addTo(map);

    markerLayer = L.layerGroup().addTo(map);
    updateMarkers(L);
  }

  function updateMarkers(L) {
    if (!markerLayer || !L) return;
    markerLayer.clearLayers();

    const validMarkers = markers.filter(m => m.lat != null && m.lon != null && !isNaN(m.lat) && !isNaN(m.lon));
    if (!validMarkers.length) return;

    // Offset overlapping markers (same lat/lon) so depth profiles are visible
    const posGroups = {};
    for (const m of validMarkers) {
      const key = `${m.lat.toFixed(4)}_${m.lon.toFixed(4)}`;
      if (!posGroups[key]) posGroups[key] = [];
      posGroups[key].push(m);
    }

    // Convert nudge meters to approximate degrees (1° lat ≈ 111,320m)
    const spreadDeg = nudgeMeters / 111320;

    for (const group of Object.values(posGroups)) {
      if (group.length <= 1 || spreadDeg === 0) continue;
      // Sort by date then depth (earliest+shallowest first = northernmost offset)
      group.sort((a, b) => {
        const da = a.date || '', db = b.date || '';
        if (da !== db) return da.localeCompare(db);
        return (a.depth_m ?? 0) - (b.depth_m ?? 0);
      });
      // Arrange in a small circle around the point
      const step = (2 * Math.PI) / group.length;
      group.forEach((m, i) => {
        m._offsetLat = m.lat + Math.cos(step * i) * spreadDeg;
        m._offsetLon = m.lon + Math.sin(step * i) * spreadDeg;
      });
    }

    for (const m of validMarkers) {
      const color = getColor(m);
      const isHighlighted = highlightIds ? highlightIds.has(m.id) : false;
      const radius = markerRadius(m);
      const weight = isHighlighted ? 1.5 : 1;
      const borderColor = isHighlighted ? '#22d3ee' : '#0f172a';
      const lat = m._offsetLat ?? m.lat;
      const lon = m._offsetLon ?? m.lon;

      const circle = L.circleMarker([lat, lon], {
        radius, fillColor: color, color: borderColor,
        weight, opacity: 1, fillOpacity: 0.8,
      });

      // Tooltip on hover
      circle.bindTooltip(tooltipContent(m), {
        className: 'dark-tooltip',
        direction: 'top', offset: [0, -radius],
      });

      if (onMarkerClick) {
        circle.on('click', (e) => {
          L.DomEvent.stopPropagation(e);
          onMarkerClick(m.id);
        });
      }

      circle.addTo(markerLayer);
    }

    // Only fit bounds on first load, not on every update
    if (!hasFitted) {
      const bounds = L.latLngBounds(validMarkers.map(m => [m.lat, m.lon]));
      map.fitBounds(bounds, { padding: [30, 30], maxZoom: 10 });
      hasFitted = true;
    }
  }

  onMount(() => {
    initMap();
    return () => { if (map) { map.remove(); map = null; } };
  });

  $effect(() => {
    const _deps = [markers, colorMap, colorBy, sizeBy, sizeScale, nudgeMeters, highlightIds];
    if (map) {
      import('leaflet').then(L => updateMarkers(L));
    }
  });
</script>

<div bind:this={mapContainer} class="w-full h-full rounded" style="min-height: 400px; z-index: 0;"></div>

<style>
  :global(.dark-tooltip) {
    background: #1e293b !important;
    border: 1px solid #475569 !important;
    color: #e2e8f0 !important;
    border-radius: 6px !important;
    padding: 6px 10px !important;
    box-shadow: 0 4px 12px rgba(0,0,0,0.5) !important;
  }
  :global(.dark-tooltip::before) {
    border-top-color: #475569 !important;
  }
</style>
