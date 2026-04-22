<script>
  import { onMount } from 'svelte';

  let { markers = [], colorMap = {}, colorBy = 'flowcell', sizeBy = 'fixed', sizeScale = 1.0, nudgeMeters = 100, onMarkerClick = null, onSelect = null, highlightIds = null, exportName = 'danaseq_map' } = $props();
  let mapContainer;
  let map = null;
  let markerLayer = null;
  let hasFitted = false;

  // Box-select state
  let selectBox = null;   // Leaflet rectangle overlay
  let selectStart = null; // {x, y} container pixel
  let renderedMarkers = []; // [{id, lat, lon}] for hit-testing

  const PALETTE = ['#22d3ee','#34d399','#fbbf24','#f87171','#a78bfa','#fb923c',
                   '#2dd4bf','#818cf8','#f472b6','#4ade80','#e879f9','#38bdf8'];

  function getColor(marker) {
    // Check color map by the active colorBy field value
    const key = marker[colorBy] || marker.station || marker.id;
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

  // Precomputed sizing normalization (set in updateMarkers)
  let sizeNorm = null;

  function markerRadius(m) {
    if (sizeBy === 'fixed') return 8 * sizeScale;
    const val = typeof m[sizeBy] === 'number' ? m[sizeBy] : 0;
    if (val <= 0 || !sizeNorm) return 3 * sizeScale;
    const t = (Math.log1p(val) - sizeNorm.logMin) / sizeNorm.range;
    return (3 + t * 27) * sizeScale; // 3px to 30px full range
  }

  function fmtVal(v) {
    if (v == null) return '-';
    if (typeof v === 'number') return Number.isInteger(v) ? v.toLocaleString() : v.toFixed(2);
    return String(v);
  }

  function fmtLabel(key) {
    return key.replace(/_/g, ' ').replace(/\b\w/g, l => l.toUpperCase());
  }

  function tooltipContent(m) {
    const label = m.station || m.id;
    const parts = [`<b>${label}</b>`];
    if (m.station) parts.push(`<span style="color:#64748b">${m.id}</span>`);
    if (m.date) parts.push(`Date: ${m.date}`);
    if (m.depth_m != null) parts.push(`Depth: ${m.depth_m}m`);
    // Show current color and size selections
    if (colorBy && m[colorBy] != null) parts.push(`${fmtLabel(colorBy)}: ${fmtVal(m[colorBy])}`);
    if (sizeBy !== 'fixed' && sizeBy !== colorBy && m[sizeBy] != null) parts.push(`${fmtLabel(sizeBy)}: ${fmtVal(m[sizeBy])}`);
    return `<div style="color:#e2e8f0;font-family:Inter,sans-serif;font-size:11px;line-height:1.4">${parts.join('<br>')}</div>`;
  }

  async function exportPng() {
    if (!map || !mapContainer) return;
    const size = map.getSize();
    const canvas = document.createElement('canvas');
    canvas.width = size.x * 2; canvas.height = size.y * 2;   // 2x for retina
    const ctx = canvas.getContext('2d');
    ctx.scale(2, 2);

    // 1. Draw tile images
    const tiles = mapContainer.querySelectorAll('.leaflet-tile');
    const tilePromises = [...tiles].map(tile => {
      return new Promise(resolve => {
        if (tile.complete) return resolve(tile);
        tile.onload = () => resolve(tile);
        tile.onerror = () => resolve(null);
      });
    });
    await Promise.all(tilePromises);

    const mapPane = mapContainer.querySelector('.leaflet-map-pane');
    const mapTransform = mapPane ? getComputedStyle(mapPane).transform : 'none';
    let mapOffsetX = 0, mapOffsetY = 0;
    if (mapTransform && mapTransform !== 'none') {
      const m = mapTransform.match(/matrix.*\((.+)\)/);
      if (m) {
        const vals = m[1].split(',').map(Number);
        mapOffsetX = vals[4] || 0; mapOffsetY = vals[5] || 0;
      }
    }

    for (const tile of tiles) {
      if (!tile || !tile.src) continue;
      try {
        const tileTransform = getComputedStyle(tile).transform;
        let tx = 0, ty = 0;
        if (tileTransform && tileTransform !== 'none') {
          const tm = tileTransform.match(/matrix.*\((.+)\)/);
          if (tm) { const v = tm[1].split(',').map(Number); tx = v[4] || 0; ty = v[5] || 0; }
        }
        const x = tile.offsetLeft + tx + mapOffsetX;
        const y = tile.offsetTop + ty + mapOffsetY;
        ctx.drawImage(tile, x, y, tile.width, tile.height);
      } catch (_) { /* CORS tiles may fail */ }
    }

    // 2. Draw circle markers
    const L = await import('leaflet');
    markerLayer.eachLayer(layer => {
      if (!layer.getLatLng) return;
      const pt = map.latLngToContainerPoint(layer.getLatLng());
      const opts = layer.options;
      ctx.beginPath();
      ctx.arc(pt.x, pt.y, opts.radius || 6, 0, Math.PI * 2);
      ctx.fillStyle = opts.fillColor || '#22d3ee';
      ctx.globalAlpha = opts.fillOpacity ?? 0.8;
      ctx.fill();
      ctx.globalAlpha = 1;
      ctx.strokeStyle = opts.color || '#0f172a';
      ctx.lineWidth = opts.weight || 1;
      ctx.stroke();
    });

    canvas.toBlob(blob => {
      const url = URL.createObjectURL(blob);
      const a = document.createElement('a');
      a.href = url; a.download = `${exportName}.png`; a.click();
      URL.revokeObjectURL(url);
    });
  }

  async function initMap() {
    const L = await import('leaflet');
    await import('leaflet/dist/leaflet.css');

    map = L.map(mapContainer, {
      center: [0, 0], zoom: 2,
      zoomControl: true,
      attributionControl: true,
      boxZoom: !onSelect, // disable default shift-drag zoom when select is enabled
    });

    L.tileLayer('https://{s}.tile.openstreetmap.org/{z}/{x}/{y}.png', {
      attribution: '&copy; OpenStreetMap',
      maxZoom: 18,
    }).addTo(map);

    markerLayer = L.layerGroup().addTo(map);

    // Shift-drag box select
    if (onSelect) {
      const container = map.getContainer();

      container.addEventListener('mousedown', (e) => {
        if (!e.shiftKey) return;
        e.stopPropagation();
        map.dragging.disable();
        const rect = container.getBoundingClientRect();
        selectStart = { x: e.clientX - rect.left, y: e.clientY - rect.top };

        selectBox = L.rectangle([[0,0],[0,0]], {
          color: '#22d3ee', weight: 2, fillColor: '#22d3ee',
          fillOpacity: 0.15, dashArray: '6 3', interactive: false,
        }).addTo(map);
      });

      container.addEventListener('mousemove', (e) => {
        if (!selectStart || !selectBox) return;
        const rect = container.getBoundingClientRect();
        const cur = { x: e.clientX - rect.left, y: e.clientY - rect.top };
        const sw = map.containerPointToLatLng([Math.min(selectStart.x, cur.x), Math.max(selectStart.y, cur.y)]);
        const ne = map.containerPointToLatLng([Math.max(selectStart.x, cur.x), Math.min(selectStart.y, cur.y)]);
        selectBox.setBounds([sw, ne]);
      });

      container.addEventListener('mouseup', (e) => {
        if (!selectStart || !selectBox) return;
        const rect = container.getBoundingClientRect();
        const cur = { x: e.clientX - rect.left, y: e.clientY - rect.top };
        const sw = map.containerPointToLatLng([Math.min(selectStart.x, cur.x), Math.max(selectStart.y, cur.y)]);
        const ne = map.containerPointToLatLng([Math.max(selectStart.x, cur.x), Math.min(selectStart.y, cur.y)]);
        const bounds = L.latLngBounds(sw, ne);

        // Find markers within the selection rectangle
        const selected = renderedMarkers.filter(m => bounds.contains([m.lat, m.lon])).map(m => m.id);
        if (selected.length > 0) {
          onSelect(selected);
        } else {
          onSelect(null);
        }

        map.removeLayer(selectBox);
        selectBox = null;
        selectStart = null;
        map.dragging.enable();
      });
    }

    updateMarkers(L);
  }

  function updateMarkers(L) {
    if (!markerLayer || !L) return;
    markerLayer.clearLayers();

    const validMarkers = markers.filter(m => m.lat != null && m.lon != null && !isNaN(m.lat) && !isNaN(m.lon));
    if (!validMarkers.length) return;

    // Precompute size normalization from data range
    sizeNorm = null;
    if (sizeBy !== 'fixed') {
      const vals = validMarkers.map(m => m[sizeBy]).filter(v => typeof v === 'number' && v > 0);
      if (vals.length > 1) {
        const logMin = Math.log1p(Math.min(...vals));
        const logMax = Math.log1p(Math.max(...vals));
        sizeNorm = { logMin, logMax, range: logMax - logMin || 1 };
      }
    }

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

    renderedMarkers = [];

    for (const m of validMarkers) {
      const color = getColor(m);
      const isHighlighted = highlightIds ? highlightIds.has(m.id) : false;
      const radius = markerRadius(m);
      const weight = isHighlighted ? 1.5 : 1;
      const borderColor = isHighlighted ? '#22d3ee' : '#0f172a';
      const lat = m._offsetLat ?? m.lat;
      const lon = m._offsetLon ?? m.lon;

      renderedMarkers.push({ id: m.id, lat, lon });

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

<div class="relative w-full h-full" style="min-height: 400px;">
  <div bind:this={mapContainer} class="w-full h-full rounded" style="z-index: 0;"></div>
  <button
    class="absolute top-2 right-2 z-[1000] p-1.5 rounded bg-white/90 text-slate-600 hover:text-cyan-600 hover:bg-white shadow transition-colors"
    onclick={exportPng}
    title="Download PNG"
  >
    <svg class="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
      <path stroke-linecap="round" stroke-linejoin="round" stroke-width="2" d="M4 16v1a3 3 0 003 3h10a3 3 0 003-3v-1m-4 0l-4 4m0 0l-4-4m4 4V4"/>
    </svg>
  </button>
</div>

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
