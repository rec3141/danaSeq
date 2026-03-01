<script>
  import { onMount } from 'svelte';

  let { markers = [], colorMap = {}, onMarkerClick = null, highlightIds = null } = $props();
  let mapContainer;
  let map = null;
  let markerLayer = null;

  const PALETTE = ['#22d3ee','#34d399','#fbbf24','#f87171','#a78bfa','#fb923c',
                   '#2dd4bf','#818cf8','#f472b6','#4ade80','#e879f9','#38bdf8'];

  function getColor(marker) {
    if (colorMap[marker.id]) return colorMap[marker.id];
    if (marker.color) return marker.color;
    return PALETTE[Math.abs(hashCode(marker.id)) % PALETTE.length];
  }

  function hashCode(s) {
    let hash = 0;
    for (let i = 0; i < s.length; i++) hash = ((hash << 5) - hash) + s.charCodeAt(i);
    return hash;
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

    for (const m of validMarkers) {
      const color = getColor(m);
      const isHighlighted = highlightIds ? highlightIds.has(m.id) : false;
      const radius = m.radius || 8;
      const weight = isHighlighted ? 3 : 1;
      const borderColor = isHighlighted ? '#22d3ee' : '#0f172a';

      const circle = L.circleMarker([m.lat, m.lon], {
        radius, fillColor: color, color: borderColor,
        weight, opacity: 1, fillOpacity: 0.8,
      });

      const popupContent = [
        `<div style="color: #e2e8f0; font-family: Inter, sans-serif; font-size: 12px;">`,
        `<strong>${m.id}</strong>`,
        m.station ? `<br>Station: ${m.station}` : '',
        m.date ? `<br>Date: ${m.date}` : '',
        m.depth_m != null ? `<br>Depth: ${m.depth_m}m` : '',
        m.read_count != null ? `<br>Reads: ${m.read_count.toLocaleString()}` : '',
        `</div>`,
      ].join('');

      circle.bindPopup(popupContent, {
        className: 'dark-popup',
        closeButton: true,
      });

      if (onMarkerClick) {
        circle.on('click', () => onMarkerClick(m.id));
      }

      circle.addTo(markerLayer);
    }

    // Fit bounds to markers
    const bounds = L.latLngBounds(validMarkers.map(m => [m.lat, m.lon]));
    map.fitBounds(bounds, { padding: [30, 30], maxZoom: 10 });
  }

  onMount(() => {
    initMap();
    return () => { if (map) { map.remove(); map = null; } };
  });

  $effect(() => {
    const _deps = [markers, colorMap, highlightIds];
    if (map) {
      import('leaflet').then(L => updateMarkers(L));
    }
  });
</script>

<div bind:this={mapContainer} class="w-full h-full rounded" style="min-height: 400px; z-index: 0;"></div>
