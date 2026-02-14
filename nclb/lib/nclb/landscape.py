"""Landscape — UMAP embedding, HDBSCAN clustering, and community naming.

Every contig lives somewhere in composition-abundance space. The landscape
module projects that high-dimensional space into a 2D map, finds natural
clusters via HDBSCAN, and gives each contig its coordinates plus the
centroid of its cluster. Contigs can judge distance directly from coordinates.

Communities get playful names (adjective + noun) instead of ugly bin IDs.
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import numpy as np
from sklearn.preprocessing import StandardScaler


# ---------------------------------------------------------------------------
# Community naming — playful, memorable, deterministic
# ---------------------------------------------------------------------------

ADJECTIVES = [
    "amber", "arctic", "astral", "azure", "blazing", "bold", "bright",
    "bronze", "celestial", "cobalt", "coral", "crimson", "crystal", "dappled",
    "dawn", "deep", "distant", "drifting", "dusky", "electric", "emerald",
    "fading", "fierce", "flickering", "frosty", "gilded", "glinting",
    "golden", "granite", "hollow", "hushed", "iron", "ivory", "jade",
    "keen", "kindled", "lavender", "leaden", "luminous", "lunar", "misty",
    "molten", "moss", "neon", "nimble", "noble", "obsidian", "opal",
    "pale", "phantom", "plum", "polar", "proud", "quiet", "radiant",
    "roaming", "ruby", "russet", "sable", "scarlet", "shadow", "silent",
    "silver", "slate", "solar", "spectral", "stark", "steady", "storm",
    "swift", "tawny", "tidal", "twilight", "umbral", "vast", "velvet",
    "verdant", "violet", "wandering", "wild", "woven", "zinc",
]

NOUNS = [
    "anchor", "arrow", "atlas", "beacon", "blade", "bolt", "bow",
    "cairn", "cinder", "comet", "compass", "corona", "crest", "crow",
    "dagger", "dart", "drake", "drift", "dusk", "echo", "ember",
    "falcon", "fern", "flare", "forge", "fox", "frost", "gale",
    "gate", "ghost", "grove", "halo", "hare", "hawk", "hearth",
    "helm", "heron", "hollow", "horn", "ibis", "isle", "jade",
    "keel", "kestrel", "lantern", "lark", "ledge", "lynx", "mantle",
    "marsh", "mesa", "moth", "nebula", "nexus", "oar", "onyx",
    "orbit", "osprey", "otter", "peak", "petal", "pike", "plume",
    "quill", "raven", "reef", "ridge", "rill", "rover", "rune",
    "sage", "shard", "shoal", "sigil", "spark", "spire", "stag",
    "stone", "surge", "talon", "tarn", "thistle", "torch", "vale",
    "vane", "vigil", "vole", "warden", "wave", "wren", "zenith",
]


def generate_community_name(seed: str, index: int = 0) -> str:
    """Generate a deterministic playful name from a seed string.

    Uses a hash to pick adjective + noun, so the same seed always
    produces the same name. The index parameter breaks ties if
    multiple communities hash to the same pair.
    """
    h = hashlib.sha256(f"{seed}:{index}".encode()).hexdigest()
    adj_idx = int(h[:8], 16) % len(ADJECTIVES)
    noun_idx = int(h[8:16], 16) % len(NOUNS)
    return f"{ADJECTIVES[adj_idx].title()} {NOUNS[noun_idx].title()}"


def name_communities(community_names: list[str]) -> dict[str, str]:
    """Assign playful names to a list of communities.

    Returns {internal_name: display_name}, guaranteed unique.
    """
    used_names = set()
    mapping = {}
    for internal in community_names:
        index = 0
        while True:
            display = generate_community_name(internal, index)
            if display not in used_names:
                break
            index += 1
        used_names.add(display)
        mapping[internal] = display
    return mapping


# ---------------------------------------------------------------------------
# Core computation
# ---------------------------------------------------------------------------

@dataclass
class LandscapeResult:
    """Result of landscape computation for a single contig."""
    x: float
    y: float
    cluster: int                     # -1 = unclustered
    cluster_cx: float = 0.0          # cluster centroid x
    cluster_cy: float = 0.0          # cluster centroid y


def compute_landscape(
    identities: dict,
    min_contig_size: int = 2000,
    umap_neighbors: int = 30,
    umap_min_dist: float = 0.1,
    hdbscan_min_cluster: int = 30,
    hdbscan_min_samples: int = 10,
    random_state: int = 42,
) -> dict[str, LandscapeResult]:
    """Compute UMAP embedding and HDBSCAN clustering for all contigs.

    Each contig gets its (x, y) position and, if clustered, the centroid
    of its cluster. Contigs can judge distance from raw coordinates.

    Returns:
        {contig_name: LandscapeResult}
    """
    import umap
    from sklearn.cluster import HDBSCAN

    # Build feature matrix for contigs above size threshold
    names = []
    tnf_list = []
    cov_list = []

    for name, c in identities.items():
        if c.size < min_contig_size:
            continue
        names.append(name)
        tnf_list.append(c.tnf)
        cov_list.append(c.coverage)

    if len(names) < 20:
        return {}

    tnf_matrix = np.array(tnf_list)
    cov_matrix = np.array(cov_list)

    # Scale: TNF weighted 2x (composition matters more than abundance)
    scaler_tnf = StandardScaler()
    scaler_cov = StandardScaler()
    tnf_scaled = scaler_tnf.fit_transform(tnf_matrix)
    cov_scaled = scaler_cov.fit_transform(cov_matrix)
    features = np.hstack([tnf_scaled * 2.0, cov_scaled])

    # UMAP embedding
    reducer = umap.UMAP(
        n_components=2,
        n_neighbors=umap_neighbors,
        min_dist=umap_min_dist,
        random_state=random_state,
        metric="euclidean",
    )
    embedding = reducer.fit_transform(features)

    # HDBSCAN on UMAP embedding (2D preserves local structure)
    clusterer = HDBSCAN(
        min_cluster_size=hdbscan_min_cluster,
        min_samples=hdbscan_min_samples,
    )
    labels = clusterer.fit_predict(embedding)

    # Compute cluster centroids
    centroids: dict[int, tuple[float, float]] = {}
    for label in set(labels):
        if label == -1:
            continue
        mask = labels == label
        centroids[label] = (
            float(embedding[mask, 0].mean()),
            float(embedding[mask, 1].mean()),
        )

    # Build results
    results: dict[str, LandscapeResult] = {}
    for i, name in enumerate(names):
        label = int(labels[i])
        cx, cy = centroids.get(label, (0.0, 0.0))
        results[name] = LandscapeResult(
            x=float(embedding[i, 0]),
            y=float(embedding[i, 1]),
            cluster=label,
            cluster_cx=cx,
            cluster_cy=cy,
        )

    return results


# ---------------------------------------------------------------------------
# Serialization
# ---------------------------------------------------------------------------

def save_landscape(results: dict[str, LandscapeResult], path: Path) -> None:
    """Save landscape data to JSON for reuse by converse."""
    data = {
        name: {
            "x": round(r.x, 4),
            "y": round(r.y, 4),
            "cluster": r.cluster,
            "cluster_cx": round(r.cluster_cx, 4),
            "cluster_cy": round(r.cluster_cy, 4),
        }
        for name, r in results.items()
    }
    with open(path, "w") as f:
        json.dump(data, f)


def load_landscape(path: Path) -> dict[str, LandscapeResult]:
    """Load landscape data from JSON."""
    with open(path) as f:
        data = json.load(f)
    return {
        name: LandscapeResult(**vals)
        for name, vals in data.items()
    }


# ---------------------------------------------------------------------------
# Visualization
# ---------------------------------------------------------------------------

def plot_landscape(
    results: dict[str, LandscapeResult],
    identities: dict,
    communities: dict,
    output_path: Path,
    community_names: dict[str, str] | None = None,
    title: str = "NCLB Contig Landscape",
) -> None:
    """Generate 4-panel landscape visualization."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from collections import Counter

    # Build arrays aligned with results
    names = list(results.keys())
    xs = np.array([results[n].x for n in names])
    ys = np.array([results[n].y for n in names])
    clusters = np.array([results[n].cluster for n in names])

    categories = []
    voices = []
    comm_labels = []
    sizes = []
    for n in names:
        c = identities.get(n)
        if not c:
            categories.append("unknown")
            voices.append(0)
            comm_labels.append("unhoused")
            sizes.append(0)
            continue
        sizes.append(c.size)
        voices.append(c.voice_strength)
        if c.community:
            categories.append("housed")
            comm_labels.append(c.community)
        elif c.voice_strength == 0:
            categories.append("voiceless")
            comm_labels.append("unhoused")
        else:
            categories.append("voiced-unhoused")
            comm_labels.append("unhoused")

    categories = np.array(categories)
    voices = np.array(voices)
    sizes = np.array(sizes)

    fig, axes = plt.subplots(2, 2, figsize=(20, 18))

    # --- Panel 1: Category ---
    ax = axes[0, 0]
    cat_config = {
        "housed": ("#2196F3", 3, 2, 0.4),
        "voiceless": ("#9E9E9E", 8, 1, 0.5),
        "voiced-unhoused": ("#FF9800", 5, 3, 0.4),
    }
    for cat in ["housed", "voiceless", "voiced-unhoused"]:
        color, s, z, alpha = cat_config[cat]
        mask = categories == cat
        ax.scatter(xs[mask], ys[mask], c=color, s=s, alpha=alpha,
                   label=f"{cat} ({mask.sum():,})", zorder=z, rasterized=True)
    ax.legend(loc="upper right", fontsize=10)
    ax.set_title("Category", fontsize=14)
    ax.set_xlabel("UMAP 1")
    ax.set_ylabel("UMAP 2")

    # --- Panel 2: All named communities ---
    ax = axes[0, 1]
    comm_counts = Counter(c for c in comm_labels if c != "unhoused")
    all_comms = [c for c, _ in comm_counts.most_common()]
    n_comms = len(all_comms)

    # Build enough distinct colors by cycling through multiple colormaps
    _cmaps = [plt.cm.tab20, plt.cm.tab20b, plt.cm.tab20c]
    comm_colors = []
    for i in range(n_comms):
        cmap_idx = min(i // 20, len(_cmaps) - 1)
        comm_colors.append(_cmaps[cmap_idx](i % 20))

    unhoused_mask = np.array([c == "unhoused" for c in comm_labels])
    ax.scatter(xs[unhoused_mask], ys[unhoused_mask],
               c="#E0E0E0", s=2, alpha=0.2, zorder=1, rasterized=True)
    for i, comm in enumerate(all_comms):
        mask = np.array([c == comm for c in comm_labels])
        display = community_names.get(comm, comm) if community_names else comm
        short = display if len(display) <= 20 else display[:18] + ".."
        ax.scatter(xs[mask], ys[mask], c=[comm_colors[i]], s=8, alpha=0.6,
                   label=f"{short} ({mask.sum()})", zorder=2, rasterized=True)
        # Label at centroid for larger communities
        if mask.sum() >= 15:
            cx, cy = xs[mask].mean(), ys[mask].mean()
            ax.text(cx, cy, display.split()[-1], fontsize=5, ha="center",
                    va="center", fontweight="bold", alpha=0.8,
                    bbox=dict(boxstyle="round,pad=0.1", facecolor="white",
                              alpha=0.6, edgecolor="none"))
    ncol = max(1, (n_comms + 14) // 15)
    ax.legend(loc="upper right", fontsize=4, ncol=ncol)
    ax.set_title(f"Communities ({n_comms} total)", fontsize=14)
    ax.set_xlabel("UMAP 1")
    ax.set_ylabel("UMAP 2")

    # --- Panel 3: Voice strength ---
    ax = axes[1, 0]
    sc = ax.scatter(xs, ys, c=voices, cmap="YlOrRd", s=3, alpha=0.4,
                    vmin=0, vmax=5, rasterized=True)
    plt.colorbar(sc, ax=ax, label="Voice strength (0-5 binners)")
    ax.set_title("Voice Strength (Binner Agreement)", fontsize=14)
    ax.set_xlabel("UMAP 1")
    ax.set_ylabel("UMAP 2")

    # --- Panel 4: HDBSCAN clusters with centroids ---
    ax = axes[1, 1]
    n_clusters = len(set(clusters)) - (1 if -1 in clusters else 0)
    n_noise = (clusters == -1).sum()
    n_clustered = len(clusters) - n_noise

    # Noise points
    noise_mask = clusters == -1
    ax.scatter(xs[noise_mask], ys[noise_mask], c="#FFE0B2", s=3, alpha=0.2,
               zorder=1, rasterized=True, label=f"unclustered ({n_noise:,})")

    # Cluster points (color by cluster, max 20 in legend)
    unique_labels = sorted(l for l in set(clusters) if l != -1)
    # Sort by size for legend (biggest first)
    label_sizes = [(l, (clusters == l).sum()) for l in unique_labels]
    label_sizes.sort(key=lambda x: -x[1])

    cluster_cmap = plt.cm.tab20
    for rank, (label, n_members) in enumerate(label_sizes):
        mask = clusters == label
        total_bp = sizes[mask].sum()
        color = cluster_cmap(rank % 20)
        show_label = rank < 15  # legend for top 15
        ax.scatter(xs[mask], ys[mask], c=[color], s=6, alpha=0.5,
                   zorder=2, rasterized=True,
                   label=f"cl-{label} ({n_members}, {total_bp/1e6:.1f}Mb)" if show_label else None)

        # Draw centroid marker for larger clusters
        if n_members >= 50:
            cx = xs[mask].mean()
            cy = ys[mask].mean()
            ax.plot(cx, cy, "x", color="black", markersize=6, markeredgewidth=1.5, zorder=3)

    ax.legend(loc="upper right", fontsize=6, ncol=2)
    ax.set_title(f"HDBSCAN ({n_clusters} clusters, {n_clustered:,} contigs)", fontsize=14)
    ax.set_xlabel("UMAP 1")
    ax.set_ylabel("UMAP 2")

    n_total = len(names)
    plt.suptitle(f"{title} — {n_total:,} contigs", fontsize=16, y=0.98)
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.savefig(output_path, dpi=450, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Summary statistics
# ---------------------------------------------------------------------------

def landscape_summary(
    results: dict[str, LandscapeResult],
    identities: dict,
) -> dict:
    """Compute summary statistics for the landscape."""
    from collections import Counter

    clusters = [r.cluster for r in results.values()]

    n_clusters = len(set(clusters)) - (1 if -1 in clusters else 0)
    n_noise = sum(1 for c in clusters if c == -1)
    n_clustered = len(clusters) - n_noise

    # Per-cluster stats
    cluster_stats = []
    for label in sorted(set(clusters)):
        if label == -1:
            continue
        members = [n for n, r in results.items() if r.cluster == label]
        member_sizes = [identities[n].size for n in members if n in identities]

        # Centroid
        member_xs = [results[n].x for n in members]
        member_ys = [results[n].y for n in members]
        cx = float(np.mean(member_xs))
        cy = float(np.mean(member_ys))

        # Taxonomy summary
        ancestries = [identities[n].ancestry for n in members
                      if n in identities and identities[n].ancestry]
        top_ancestry = ""
        if ancestries:
            phyla = []
            for a in ancestries:
                parts = a.split(";")
                if len(parts) >= 2:
                    phyla.append(parts[1].strip())
            if phyla:
                top_ancestry = Counter(phyla).most_common(1)[0][0]

        cluster_stats.append({
            "cluster": label,
            "n_members": len(members),
            "total_size": sum(member_sizes),
            "median_size": int(np.median(member_sizes)) if member_sizes else 0,
            "centroid": [round(cx, 2), round(cy, 2)],
            "top_phylum": top_ancestry,
        })

    return {
        "n_contigs_mapped": len(results),
        "n_clusters": n_clusters,
        "n_clustered": n_clustered,
        "n_unclustered": n_noise,
        "clusters": cluster_stats,
    }
