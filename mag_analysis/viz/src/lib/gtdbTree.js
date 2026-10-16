// GTDB-Tk trees carry the reference taxonomy on their internal nodes, as
// labels like '100.0:o__Brevinematales' or '100.0:p__Foo; c__Bar' (support,
// then one or more ranks). Reading a leaf's ancestors therefore gives its
// lineage down to whatever rank the tree labels, without the GTDB database.

const RANK_ORDER = ['d', 'p', 'c', 'o', 'f', 'g', 's'];

/** Taxa in an internal label: "100.0:p__Foo; c__Bar" -> ['p__Foo', 'c__Bar']. */
function labelTaxa(label) {
  const taxa = [];
  const body = label.includes(':') ? label.slice(label.indexOf(':') + 1) : label;
  for (const seg of body.split(';')) {
    const t = seg.trim();
    if (/^[dpcofgs]__./.test(t)) taxa.push(t);
  }
  return taxa;
}

/**
 * Parse a GTDB-Tk Newick string once.
 *
 * Returns
 *   newick       the same tree with each internal label reduced to its taxa
 *                (p__Foo / c__Bar), unquoted because Phylocanvas's parser
 *                doesn't read quotes, and support-only labels dropped;
 *   leafLineage  Map leaf name -> ['d__Bacteria', 'p__…', …], root first.
 */
export function parseGtdbTree(nwk) {
  const parent = [];   // node index -> parent index (-1 for the root)
  const taxa = [];     // node index -> taxa on that node's label
  const leaves = [];   // [name, node index]
  const out = [];
  const stack = [];
  let cur = -1;
  let i = 0;

  const newNode = (p) => { parent.push(p); taxa.push(null); return parent.length - 1; };

  // Read a (possibly quoted) label starting at i; returns [text, raw].
  function readLabel() {
    if (nwk[i] === "'") {
      let j = i + 1, text = '';
      while (j < nwk.length) {
        if (nwk[j] === "'") {
          if (nwk[j + 1] === "'") { text += "'"; j += 2; continue; }
          break;
        }
        text += nwk[j++];
      }
      const raw = nwk.slice(i, j + 1);
      i = j + 1;
      return [text, raw];
    }
    let j = i;
    while (j < nwk.length && !':,();'.includes(nwk[j])) j++;
    const raw = nwk.slice(i, j);
    i = j;
    return [raw.trim(), raw];
  }

  while (i < nwk.length) {
    const ch = nwk[i];
    if (ch === '(') {
      stack.push(cur);
      cur = newNode(cur);
      out.push(ch);
      i++;
    } else if (ch === ')') {
      out.push(ch);
      i++;
      const [text, raw] = readLabel();
      const t = labelTaxa(text);
      if (t.length) {
        taxa[cur] = t;
        out.push(t.map(x => x.replace(/[\s'(),:;[\]]/g, '_')).join(' / '));
      } else if (text && !/^[\d.eE+-]+$/.test(text)) {
        out.push(raw); // a non-taxon, non-support label: keep it
      }
      cur = stack.pop();
    } else if (ch === ',' || ch === ';' || /\s/.test(ch)) {
      out.push(ch);
      i++;
    } else if (ch === ':') {
      let j = i + 1;
      while (j < nwk.length && !',();'.includes(nwk[j])) j++;
      out.push(nwk.slice(i, j));
      i = j;
    } else {
      const [name, raw] = readLabel();
      leaves.push([name, newNode(cur)]);
      out.push(raw);
    }
  }

  const leafLineage = new Map();
  for (const [name, idx] of leaves) {
    const lineage = [];
    for (let n = parent[idx]; n >= 0; n = parent[n]) {
      if (taxa[n]) lineage.unshift(...taxa[n]);
    }
    leafLineage.set(name, lineage);
  }
  return { newick: out.join(''), leafLineage };
}

/** Lineage shared by every given lineage (the LCA's taxonomy). */
export function commonLineage(lineages) {
  let common = null;
  for (const l of lineages) {
    if (!l) continue;
    if (common === null) { common = l.slice(); continue; }
    let k = 0;
    while (k < common.length && k < l.length && common[k] === l[k]) k++;
    common.length = k;
    if (!k) break;
  }
  return common ?? [];
}

/** ['d__Bacteria', …] -> 'd__Bacteria;…', the form parseLineage() reads. */
export function lineageString(lineage) {
  return [...lineage]
    .sort((a, b) => RANK_ORDER.indexOf(a[0]) - RANK_ORDER.indexOf(b[0]))
    .join(';');
}
