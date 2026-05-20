import { writable, derived } from 'svelte/store';

export const overview     = writable(null);
export const samples      = writable([]);
export const references   = writable([]);
export const readFlow     = writable(null);
export const expression   = writable(null);  // { <ref>: {genes, samples, matrix} }

export const loading = writable(true);
export const error   = writable(null);

async function fetchJSON(url) {
  // Resolve absolute /data/* against Vite's BASE_URL so the SPA works under a subdir
  if (url.startsWith('/data/')) {
    const base = (import.meta.env.BASE_URL || '/').replace(/\/$/, '');
    url = base + url;
  }
  try {
    const gzRes = await fetch(url + '.gz');
    if (gzRes.ok) {
      const ct = gzRes.headers.get('Content-Type') || '';
      if (ct.includes('json') || ct.includes('gzip') || ct.includes('octet-stream')) {
        const ce = gzRes.headers.get('Content-Encoding');
        if (ce && ce.includes('gzip')) return gzRes.json();
        const ds = new DecompressionStream('gzip');
        const text = await new Response(gzRes.body.pipeThrough(ds)).text();
        return JSON.parse(text);
      }
    }
  } catch (_) { /* fall through */ }
  const res = await fetch(url);
  if (!res.ok) throw new Error(`Failed to load ${url}: ${res.status}`);
  const ct = res.headers.get('Content-Type') || '';
  if (ct.includes('html')) return null;
  return res.json();
}

export async function loadAllData() {
  loading.set(true);
  error.set(null);
  try {
    const [ov, sm, rf, rd, ex] = await Promise.all([
      fetchJSON('./data/overview.json').catch(() => null),
      fetchJSON('./data/samples.json').catch(() => []),
      fetchJSON('./data/references.json').catch(() => []),
      fetchJSON('./data/read_flow.json').catch(() => null),
      fetchJSON('./data/expression.json').catch(() => null),
    ]);
    overview.set(ov);
    samples.set(sm || []);
    references.set(rf || []);
    readFlow.set(rd);
    expression.set(ex);
  } catch (e) {
    error.set(e.message);
    console.error('Data load error:', e);
  } finally {
    loading.set(false);
  }
}
