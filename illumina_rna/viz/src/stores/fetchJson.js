// Fetch a data file, preferring the pre-compressed <url>.gz.
//
// dataBase is prepended to every URL. The SPA leaves it empty (URLs resolve
// against the page); scripts running outside a browser page set an absolute
// origin with setDataBase().
let dataBase = '';
export function setDataBase(base) { dataBase = base; }

export async function fetchJSON(url) {
  url = dataBase + url;
  // Try pre-compressed .json.gz first — works with any static server, no server config needed.
  // If the server sets Content-Encoding: gzip (e.g. vite preview), the browser auto-decompresses
  // and we can just parse directly. Otherwise we decompress manually via DecompressionStream.
  const gzRes = await fetch(url + '.gz?t=' + Date.now());
  const gzContentType = (gzRes.headers.get('Content-Type') || '').toLowerCase();
  if (gzRes.ok && !gzContentType.includes('text/html')) {
    const ce = (gzRes.headers.get('Content-Encoding') || '').toLowerCase();
    if (ce.includes('gzip') || ce.includes('br') || ce.includes('deflate')) {
      // Browser already decompressed transparently
      return gzRes.json();
    }
    // Sniff first two bytes — if they match the gzip magic number (1f 8b),
    // decompress manually; otherwise the server already decoded for us.
    const buf = await gzRes.arrayBuffer();
    const header = new Uint8Array(buf, 0, 2);
    if (header[0] === 0x1f && header[1] === 0x8b) {
      const ds = new DecompressionStream('gzip');
      const text = await new Response(
        new Blob([buf]).stream().pipeThrough(ds)
      ).text();
      return JSON.parse(text);
    }
    // Not actually gzipped (server decoded it) — parse as plain JSON
    const text = new TextDecoder().decode(buf);
    return JSON.parse(text);
  }
  // Fall back to plain JSON (older preprocess runs, dev environments)
  const res = await fetch(url + '?t=' + Date.now());
  if (!res.ok) throw new Error(`Failed to load ${url}: ${res.status}`);
  return res.json();
}
