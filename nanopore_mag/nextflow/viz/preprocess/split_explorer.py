#!/usr/bin/env python3
"""Split old-format contig_explorer.json into two files:
  - contig_explorer.json: contigs with PCA + metadata (no embeddings)
  - contig_embeddings.json: {contig_id: {tsne_x, tsne_y, umap_x, umap_y}}

Usage: python split_explorer.py <contig_explorer.json> [output_dir]
"""
import json
import os
import sys


def main():
    src = sys.argv[1]
    out_dir = sys.argv[2] if len(sys.argv) > 2 else os.path.dirname(src)

    with open(src) as f:
        data = json.load(f)

    embeddings = {}
    embedding_keys = {'tsne_x', 'tsne_y', 'umap_x', 'umap_y'}

    for c in data.get('contigs', []):
        emb = {}
        for k in list(c.keys()):
            if k in embedding_keys:
                emb[k] = c.pop(k)
        if emb:
            embeddings[c['id']] = emb

    # Update flags
    n_tsne = sum(1 for v in embeddings.values() if 'tsne_x' in v)
    n_umap = sum(1 for v in embeddings.values() if 'umap_x' in v)
    data['has_tsne'] = n_tsne > 0
    data['has_umap'] = n_umap > 0

    # Write explorer (without embeddings)
    explorer_path = os.path.join(out_dir, 'contig_explorer.json')
    with open(explorer_path, 'w') as f:
        json.dump(data, f, separators=(',', ':'))
    print(f"  Wrote {explorer_path}: {len(data['contigs'])} contigs, {os.path.getsize(explorer_path)/1e6:.1f} MB")

    # Write embeddings
    if embeddings:
        emb_path = os.path.join(out_dir, 'contig_embeddings.json')
        with open(emb_path, 'w') as f:
            json.dump(embeddings, f, separators=(',', ':'))
        print(f"  Wrote {emb_path}: t-SNE={n_tsne}, UMAP={n_umap}, {os.path.getsize(emb_path)/1e6:.1f} MB")
    else:
        print("  No embeddings found in source file")


if __name__ == '__main__':
    main()
