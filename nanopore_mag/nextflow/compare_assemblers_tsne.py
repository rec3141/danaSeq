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
import plotly.express as px


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

    # Build interactive plot
    fig = px.scatter(
        df,
        x="tsne_x", y="tsne_y",
        color="assembler",
        size="log_length",
        size_max=12,
        hover_data={
            "contig": True,
            "length": ":,",
            "cov.": ":.1f",
            "circ.": True,
            "assembler": True,
            "tsne_x": False,
            "tsne_y": False,
            "log_length": False,
        },
        title=f"t-SNE of Tetranucleotide Frequencies — {len(df)} contigs >= {args.min_length} bp",
        labels={"tsne_x": "t-SNE 1", "tsne_y": "t-SNE 2"},
        opacity=0.6,
        color_discrete_map={
            "flye": "#636EFA",
            "metamdbg": "#EF553B",
            "myloasm": "#00CC96",
        },
    )
    fig.update_layout(
        width=1200, height=800,
        legend_title_text="Assembler",
        template="plotly_white",
    )
    fig.update_traces(marker=dict(line=dict(width=0.3, color="DarkSlateGrey")))

    fig.write_html(args.out, include_plotlyjs=True)
    print(f"Saved interactive plot: {args.out}")

    # Summary stats
    for label in df["assembler"].unique():
        sub = df[df["assembler"] == label]
        print(f"  {label}: {len(sub)} contigs, median length {sub['length'].median():.0f}")


if __name__ == "__main__":
    main()
