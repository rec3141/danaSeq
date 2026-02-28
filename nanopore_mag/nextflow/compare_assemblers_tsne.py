#!/usr/bin/env python3
"""
Merge TNF profiles from multiple assembler runs, compute t-SNE,
and produce an interactive Plotly HTML for exploring assembly space.

Usage:
    python3 compare_assemblers_tsne.py \
        --assemblers flye:test_results_flye \
                     metamdbg:test_results_metamdbg \
                     myloasm:test_results_myloasm \
        --out assembler_tsne.html
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.manifold import TSNE
import plotly.graph_objects as go
from sklearn.cluster import KMeans


def load_assembler(label, results_dir):
    """Load TNF + assembly_info for one assembler run."""
    rd = Path(results_dir)
    tnf_path = rd / "assembly" / "tnf.tsv"
    info_path = rd / "assembly" / "assembly_info.txt"

    if not tnf_path.exists():
        sys.exit(f"[ERROR] TNF file not found: {tnf_path}")
    if not info_path.exists():
        sys.exit(f"[ERROR] Assembly info not found: {info_path}")

    # TNF: no header, first col is contig name, rest are 136 freq values
    tnf = pd.read_csv(tnf_path, sep="\t", header=None)
    tnf.columns = ["contig"] + [f"tnf_{i}" for i in range(1, tnf.shape[1])]

    # Assembly info: has header
    info = pd.read_csv(info_path, sep="\t")
    info.rename(columns={"#seq_name": "contig"}, inplace=True)

    # Merge
    df = tnf.merge(info[["contig", "length", "cov.", "circ."]], on="contig", how="left")

    # Prefix contig names with assembler label
    df = df.copy()
    df["contig"] = label + "_" + df["contig"]
    df["assembler"] = label

    return df


def main():
    parser = argparse.ArgumentParser(description="t-SNE comparison of assembler TNF profiles")
    parser.add_argument(
        "--assemblers", nargs="+", required=True,
        help="label:results_dir pairs, e.g. flye:test_results_flye"
    )
    parser.add_argument("--out", default="assembler_tsne.html", help="Output HTML file")
    parser.add_argument("--perplexity", type=float, default=30, help="t-SNE perplexity")
    parser.add_argument("--min-length", type=int, default=5000,
                        help="Minimum contig length to include (default: 5000)")
    parser.add_argument("--kmeans", type=int, default=128, help="Number of k-means clusters")
    args = parser.parse_args()

    # Load all assemblers
    frames = []
    for spec in args.assemblers:
        if ":" not in spec:
            sys.exit(f"[ERROR] Expected label:path format, got: {spec}")
        label, path = spec.split(":", 1)
        print(f"Loading {label} from {path}...")
        frames.append(load_assembler(label, path))

    df = pd.concat(frames, ignore_index=True)
    print(f"Total contigs loaded: {len(df)}")

    # Filter by minimum length
    df = df[df["length"] >= args.min_length].reset_index(drop=True)
    print(f"Contigs >= {args.min_length} bp: {len(df)}")

    if len(df) < 10:
        sys.exit("[ERROR] Too few contigs after filtering")

    # Extract TNF matrix
    tnf_cols = [c for c in df.columns if c.startswith("tnf_")]
    X = df[tnf_cols].values

    # t-SNE
    perp = min(args.perplexity, (len(df) - 1) / 3)  # can't exceed (n-1)/3
    print(f"Running t-SNE (perplexity={perp:.0f}, n={len(df)})...")
    tsne = TSNE(n_components=2, perplexity=perp, random_state=42, max_iter=1000)
    embedding = tsne.fit_transform(X)
    df["tsne_x"] = embedding[:, 0]
    df["tsne_y"] = embedding[:, 1]

    # Convert coverage to numeric
    df["cov."] = pd.to_numeric(df["cov."], errors="coerce").fillna(0)

    # Log-transform length for sizing
    df["log_length"] = np.log10(df["length"])
    sizes = (df["log_length"] - df["log_length"].min()) / (df["log_length"].max() - df["log_length"].min()) * 5 + 1.5

    # k-means clustering on t-SNE embedding
    k = args.kmeans
    print(f"Running k-means (k={k})...")
    km = KMeans(n_clusters=k, random_state=42, n_init=10)
    df["cluster"] = km.fit_predict(embedding)

    # Hover text (shared by both layers)
    hover = [
        f"<b>{r.contig}</b><br>assembler: {r.assembler}<br>cluster: {r.cluster}"
        f"<br>length: {r.length:,}<br>cov: {r['cov.']:.1f}<br>circ: {r['circ.']}"
        for _, r in df.iterrows()
    ]

    # --- Layer 1: color by assembler ---
    asm_colors = {"flye": "#636EFA", "metamdbg": "#EF553B", "myloasm": "#00CC96"}
    asm_traces = []
    for asm in ["flye", "metamdbg", "myloasm"]:
        mask = df["assembler"] == asm
        asm_traces.append(go.Scattergl(
            x=df.loc[mask, "tsne_x"], y=df.loc[mask, "tsne_y"],
            mode="markers",
            name=asm,
            marker=dict(size=sizes[mask], color=asm_colors[asm], opacity=0.6,
                        line=dict(width=0.3, color="DarkSlateGrey")),
            text=[hover[i] for i in mask[mask].index],
            hoverinfo="text",
            visible=True,
        ))

    # --- Layer 2: color by cluster ---
    # Generate distinct colors via a cyclical palette
    palette = [
        '#e6194b','#3cb44b','#ffe119','#4363d8','#f58231','#911eb4','#42d4f4',
        '#f032e6','#bfef45','#fabed4','#469990','#dcbeff','#9A6324','#fffac8',
        '#800000','#aaffc3','#808000','#ffd8b1','#000075','#a9a9a9','#000000',
        '#e6beff','#1abc9c','#e74c3c','#2ecc71','#3498db','#9b59b6','#f39c12',
    ]
    cluster_traces = []
    for ci in sorted(df["cluster"].unique()):
        mask = df["cluster"] == ci
        color = palette[ci % len(palette)]
        cluster_traces.append(go.Scattergl(
            x=df.loc[mask, "tsne_x"], y=df.loc[mask, "tsne_y"],
            mode="markers",
            name=f"cluster {ci}",
            marker=dict(size=sizes[mask], color=color, opacity=0.6,
                        line=dict(width=0.3, color="DarkSlateGrey")),
            text=[hover[i] for i in mask[mask].index],
            hoverinfo="text",
            visible=False,
            legendgroup=f"c{ci}",
        ))

    fig = go.Figure(data=asm_traces + cluster_traces)

    n_asm = len(asm_traces)
    n_clust = len(cluster_traces)

    fig.update_layout(
        width=1200, height=800,
        template="plotly_white",
        dragmode="pan",
        title=f"t-SNE of Tetranucleotide Frequencies — {len(df)} contigs >= {args.min_length} bp",
        xaxis_title="t-SNE 1", yaxis_title="t-SNE 2",
        legend_title_text="Assembler",
        updatemenus=[dict(
            type="buttons",
            direction="left",
            x=0.0, y=1.12,
            xanchor="left",
            buttons=[
                dict(
                    label="Color: Assembler",
                    method="update",
                    args=[
                        {"visible": [True]*n_asm + [False]*n_clust},
                        {"legend_title_text": "Assembler", "showlegend": True},
                    ],
                ),
                dict(
                    label="Color: Cluster",
                    method="update",
                    args=[
                        {"visible": [False]*n_asm + [True]*n_clust},
                        {"legend_title_text": f"Cluster (k={k})", "showlegend": True},
                    ],
                ),
            ],
        )],
    )

    fig.write_html(args.out, include_plotlyjs=True,
                   config={"scrollZoom": True})
    print(f"Saved interactive plot: {args.out}")

    # Summary stats
    for label in df["assembler"].unique():
        sub = df[df["assembler"] == label]
        print(f"  {label}: {len(sub)} contigs, median length {sub['length'].median():.0f}")


if __name__ == "__main__":
    main()
