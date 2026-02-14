"""Contig identity cards and community profiles.

Every contig carries an innate identity derived from its physical nature.
Communities are emergent groups whose members resonate together.
"""

from __future__ import annotations

import csv
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import numpy as np


# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------

@dataclass
class ContigIdentity:
    """Everything a contig knows about itself."""

    name: str
    size: int                                         # bp
    tnf: np.ndarray                                   # 136-dim composition
    coverage: np.ndarray                              # per-sample energy
    gc: float                                         # GC spirit (derived from TNF)
    assembly_coverage: float = 0.0                    # Flye coverage estimate
    is_circular: bool = False
    is_repeat: bool = False
    multiplicity: int = 1

    # Ancestry (filled later if taxonomy available)
    ancestry: Optional[str] = None

    # Gifts — gene names from Prokka (filled later)
    gifts: list[str] = field(default_factory=list)
    coding_density: float = 0.0
    marker_genes: list[str] = field(default_factory=list)

    # Social connections — assembly graph neighbors (contig names)
    connections: list[str] = field(default_factory=list)

    # Oracle testimony — what each binner said
    testimony: dict[str, Optional[str]] = field(default_factory=dict)
    voice_strength: int = 0                           # how many oracles spoke

    # Mobile genetic element annotations (from geNomad/CheckV)
    is_viral: bool = False                            # geNomad virus detection
    virus_score: float = 0.0                          # geNomad virus confidence
    virus_hallmarks: int = 0                          # viral hallmark gene count
    virus_taxonomy: Optional[str] = None              # viral taxonomy from geNomad
    is_plasmid: bool = False                          # geNomad plasmid detection
    plasmid_score: float = 0.0                        # geNomad plasmid confidence
    plasmid_hallmarks: int = 0                        # plasmid hallmark gene count
    conjugation_genes: list[str] = field(default_factory=list)  # tra/vir genes
    amr_genes: list[str] = field(default_factory=list)          # antimicrobial resistance
    is_provirus: bool = False                         # CheckV provirus detection
    proviral_length: int = 0                          # length of proviral region
    checkv_quality: Optional[str] = None              # CheckV quality tier
    viral_genes: int = 0                              # CheckV viral gene count
    host_genes: int = 0                               # CheckV host gene count
    mge_type: Optional[str] = None                    # virus|plasmid|provirus|None

    # Integron annotations (from IntegronFinder)
    integrons: list[dict] = field(default_factory=list)  # [{id, type, n_attC, n_proteins}]
    has_integron: bool = False

    # Genomic island annotations (from IslandPath-DIMOB)
    genomic_islands: list[dict] = field(default_factory=list)  # [{start, end, length}]
    has_genomic_island: bool = False

    # Secretion / conjugation systems (from MacSyFinder)
    secretion_systems: list[dict] = field(default_factory=list)  # [{sys_id, model, wholeness, genes}]
    has_secretion_system: bool = False

    # Defense systems (from DefenseFinder)
    defense_systems: list[dict] = field(default_factory=list)  # [{sys_id, type, subtype, genes}]
    has_defense_system: bool = False

    # Landscape position (from UMAP embedding)
    landscape_x: float = 0.0                          # UMAP x coordinate
    landscape_y: float = 0.0                          # UMAP y coordinate
    landscape_cluster: int = -1                       # HDBSCAN cluster (-1 = unclustered)
    landscape_cluster_cx: float = 0.0                 # cluster centroid x
    landscape_cluster_cy: float = 0.0                 # cluster centroid y

    # Current state
    community: Optional[str] = None                   # current community name
    membership_type: str = "unhoused"                 # core|accessory|traveler|unhoused
    valence: float = 0.0                              # current valence score


@dataclass
class CommunityProfile:
    """Collective identity of a community (bin)."""

    name: str
    source_binner: str                                # which oracle created it
    members: list[str] = field(default_factory=list)  # contig names

    # Collective composition
    tnf_centroid: Optional[np.ndarray] = None         # mean TNF
    mean_coverage: Optional[np.ndarray] = None        # mean coverage profile
    mean_gc: float = 0.0
    gc_stdev: float = 0.0

    # Genome metrics
    total_size: int = 0
    n50: int = 0
    completeness: float = 0.0                         # SCG completeness %
    redundancy: float = 0.0                           # SCG redundancy %
    contamination: float = 0.0                        # CheckM2 contamination %
    checkm2_completeness: float = 0.0                 # CheckM2 completeness %

    # Harmony metrics (computed by valence module)
    tnf_coherence: float = 0.0                        # mean cosine to centroid
    coverage_correlation: float = 0.0                 # mean pairwise Pearson
    graph_connectivity: float = 0.0                   # fraction of pairs with edges

    # Gene inventory
    marker_gene_inventory: list[str] = field(default_factory=list)
    missing_markers: list[str] = field(default_factory=list)

    # Elder status
    elder_rank: str = "none"                          # none|apprentice|full|sage|contigsattva


# ---------------------------------------------------------------------------
# Canonical tetramers (same order as tetramer_freqs.py)
# ---------------------------------------------------------------------------

def _canonical_tetramers() -> list[str]:
    """Generate the 136 canonical (RC-collapsed) tetramers in alphabetical order."""
    from itertools import product
    bases = "ACGT"
    rc_map = str.maketrans("ACGT", "TGCA")
    seen = set()
    canonical = []
    for combo in product(bases, repeat=4):
        kmer = "".join(combo)
        rc = kmer.translate(rc_map)[::-1]
        pair = min(kmer, rc)
        if pair not in seen:
            seen.add(pair)
            canonical.append(pair)
    return sorted(canonical)


CANONICAL_TETRAMERS = _canonical_tetramers()

# GC weight vector: for each canonical tetramer, fraction of G+C bases
GC_WEIGHTS = np.array([
    sum(1 for b in kmer if b in "GC") / 4.0
    for kmer in CANONICAL_TETRAMERS
])


# ---------------------------------------------------------------------------
# Parsers
# ---------------------------------------------------------------------------

def load_tnf(path: Path) -> dict[str, np.ndarray]:
    """Load tetranucleotide frequencies from tnf.tsv.

    Format: no header, tab-separated: contig_name  freq1  freq2  ...  freq136
    """
    tnf = {}
    with open(path) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            name = parts[0]
            freqs = np.array([float(x) for x in parts[1:]], dtype=np.float64)
            tnf[name] = freqs
    return tnf


def load_depths(path: Path) -> tuple[list[str], dict[str, int], dict[str, np.ndarray]]:
    """Load CoverM depth table.

    Format: header row, then contigName  contigLen  totalAvgDepth  sample1.bam  sample1.bam-var  ...
    Returns: (sample_names, contig_lengths, contig_coverages)
    """
    lengths = {}
    coverages = {}
    sample_names = []

    with open(path) as f:
        header = f.readline().rstrip("\n").split("\t")
        # Coverage columns are every other column starting at index 3
        # (skip contigName, contigLen, totalAvgDepth, then bam, bam-var alternating)
        for i in range(3, len(header), 2):
            # Strip .sorted.bam suffix to get sample name
            name = header[i].replace(".sorted.bam", "")
            sample_names.append(name)

        for line in f:
            parts = line.rstrip("\n").split("\t")
            contig = parts[0]
            length = int(parts[1])
            lengths[contig] = length
            # Extract mean coverage per sample (skip variance columns)
            cov = []
            for i in range(3, len(parts), 2):
                cov.append(float(parts[i]))
            coverages[contig] = np.array(cov, dtype=np.float64)

    return sample_names, lengths, coverages


def load_assembly_info(path: Path) -> dict[str, dict]:
    """Load Flye assembly_info.txt.

    Format: #seq_name  length  cov.  circ.  repeat  mult.  alt_group  graph_path
    """
    info = {}
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            name = parts[0]
            info[name] = {
                "length": int(parts[1]),
                "coverage": float(parts[2]),
                "circular": parts[3] == "Y",
                "repeat": parts[4] == "Y",
                "multiplicity": int(parts[5]),
            }
    return info


def load_gfa_graph(path: Path) -> tuple[dict[str, list[str]], dict[str, list[str]]]:
    """Parse Flye GFA 1.0 to build contig adjacency.

    Returns: (contig_to_edges, edge_adjacency)
    where contig_to_edges maps contig→[edge_names] from P-lines
    and edge_adjacency maps edge→[neighbor_edges] from L-lines.

    Then contig adjacency is derived: two contigs are neighbors if they
    share an edge, or their edges are linked.
    """
    # P-lines: contig → list of edges
    contig_edges: dict[str, set[str]] = {}
    # L-lines: edge → set of neighbor edges
    edge_adj: dict[str, set[str]] = {}

    with open(path) as f:
        for line in f:
            if line.startswith("P\t"):
                parts = line.rstrip("\n").split("\t")
                contig = parts[1]
                edge_list = parts[2].split(",")
                edges = set()
                for e in edge_list:
                    edge_name = e.rstrip("+-")
                    edges.add(edge_name)
                contig_edges[contig] = edges

            elif line.startswith("L\t"):
                parts = line.rstrip("\n").split("\t")
                e1 = parts[1]
                e2 = parts[3]
                edge_adj.setdefault(e1, set()).add(e2)
                edge_adj.setdefault(e2, set()).add(e1)

    # Build contig adjacency: two contigs are neighbors if they share
    # an edge or if their edges are linked
    # First, reverse map: edge → contigs that use it
    edge_to_contigs: dict[str, set[str]] = {}
    for contig, edges in contig_edges.items():
        for edge in edges:
            edge_to_contigs.setdefault(edge, set()).add(contig)

    contig_adj: dict[str, set[str]] = {}
    for contig, edges in contig_edges.items():
        neighbors = set()
        for edge in edges:
            # Contigs sharing this edge
            for other in edge_to_contigs.get(edge, set()):
                if other != contig:
                    neighbors.add(other)
            # Contigs whose edges are linked to this edge
            for linked_edge in edge_adj.get(edge, set()):
                for other in edge_to_contigs.get(linked_edge, set()):
                    if other != contig:
                        neighbors.add(other)
        contig_adj[contig] = neighbors

    return {k: sorted(v) for k, v in contig_adj.items()}


def load_binner_assignments(paths: dict[str, Path]) -> dict[str, dict[str, Optional[str]]]:
    """Load contig→bin assignments from each binner.

    Args:
        paths: {binner_name: path_to_contig_bins.tsv}

    Returns:
        {contig_name: {binner_name: bin_label_or_None}}
    """
    all_contigs: dict[str, dict[str, Optional[str]]] = {}

    for binner, path in paths.items():
        with open(path) as f:
            for line in f:
                parts = line.rstrip("\n").split("\t")
                if len(parts) >= 2:
                    contig, bin_label = parts[0], parts[1]
                    if contig not in all_contigs:
                        all_contigs[contig] = {}
                    all_contigs[contig][binner] = bin_label

    # Fill in None for binners that didn't assign a contig
    binner_names = list(paths.keys())
    for contig in all_contigs:
        for binner in binner_names:
            if binner not in all_contigs[contig]:
                all_contigs[contig][binner] = None

    return all_contigs


def load_consensus(path: Path) -> dict[str, str]:
    """Load DAS Tool contig2bin.tsv.

    Format: no header, tab-separated: contig_name  bin_name
    """
    assignments = {}
    with open(path) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2:
                assignments[parts[0]] = parts[1]
    return assignments


def load_dastool_summary(path: Path) -> dict[str, dict]:
    """Load DAS Tool summary.tsv with quality metrics per bin."""
    bins = {}
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            bins[row["bin"]] = {
                "source": row["bin_set"],
                "size": int(row["size"]),
                "contigs": int(row["contigs"]),
                "n50": int(row["N50"]),
                "completeness": float(row["SCG_completeness"]),
                "redundancy": float(row["SCG_redundancy"]),
                "score": float(row["bin_score"]),
            }
    return bins


def load_genomad_viruses(path: Path) -> dict[str, dict]:
    """Load geNomad virus_summary.tsv."""
    viruses = {}
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            viruses[row["seq_name"]] = {
                "virus_score": float(row.get("virus_score", 0)),
                "n_hallmarks": int(row.get("n_hallmarks", 0)),
                "taxonomy": row.get("taxonomy", ""),
                "n_genes": int(row.get("n_genes", 0)),
                "topology": row.get("topology", ""),
            }
    return viruses


def load_genomad_plasmids(path: Path) -> dict[str, dict]:
    """Load geNomad plasmid_summary.tsv."""
    plasmids = {}
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            conj = row.get("conjugation_genes", "NA")
            conj_list = [g for g in conj.split(";") if g and g != "NA"]
            amr = row.get("amr_genes", "NA")
            amr_list = [g for g in amr.split(";") if g and g != "NA"]
            plasmids[row["seq_name"]] = {
                "plasmid_score": float(row.get("plasmid_score", 0)),
                "n_hallmarks": int(row.get("n_hallmarks", 0)),
                "conjugation_genes": conj_list,
                "amr_genes": amr_list,
                "n_genes": int(row.get("n_genes", 0)),
            }
    return plasmids


def load_checkv_quality(path: Path) -> dict[str, dict]:
    """Load CheckV quality_summary.tsv."""
    quality = {}
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            proviral_len = row.get("proviral_length", "NA")
            def _safe_float(val, default=0.0):
                try:
                    return float(val)
                except (ValueError, TypeError):
                    return default

            def _safe_int(val, default=0):
                try:
                    return int(val)
                except (ValueError, TypeError):
                    return default

            quality[row["contig_id"]] = {
                "provirus": row.get("provirus", "No") == "Yes",
                "proviral_length": int(proviral_len) if proviral_len != "NA" else 0,
                "viral_genes": _safe_int(row.get("viral_genes", 0)),
                "host_genes": _safe_int(row.get("host_genes", 0)),
                "checkv_quality": row.get("checkv_quality", ""),
                "completeness": _safe_float(row.get("completeness", 0)),
                "contamination": _safe_float(row.get("contamination", 0)),
            }
    return quality


def load_kaiju_taxonomy(path: Path) -> dict[str, dict]:
    """Load Kaiju per-contig taxonomy from kaiju_contigs.tsv.

    Format: contig_id  classified_genes  total_genes  fraction_classified  taxon_id  lineage
    """
    taxonomy = {}
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            lineage = row.get("lineage", "").strip().rstrip(";")
            if lineage == "Unclassified" or not lineage:
                continue
            taxonomy[row["contig_id"]] = {
                "lineage": lineage,
                "taxon_id": row.get("taxon_id", ""),
                "classified_genes": int(row.get("classified_genes", 0)),
                "total_genes": int(row.get("total_genes", 0)),
                "fraction_classified": float(row.get("fraction_classified", 0)),
            }
    return taxonomy


def load_checkm2(path: Path) -> dict[str, dict]:
    """Load CheckM2 quality_report.tsv."""
    quality = {}
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            quality[row["Name"]] = {
                "completeness": float(row["Completeness"]),
                "contamination": float(row["Contamination"]),
                "coding_density": float(row.get("Coding_Density", 0)),
                "gc": float(row.get("GC_Content", 0)),
                "genome_size": int(row.get("Genome_Size", 0)),
                "n50": int(row.get("Contig_N50", 0)),
            }
    return quality


def load_integrons(path: Path) -> dict[str, list[dict]]:
    """Load IntegronFinder results, grouped by contig.

    Returns {contig_name: [{id, type, n_attC, n_proteins}]}.
    """
    integrons: dict[str, list[dict]] = {}
    current: dict[str, dict] = {}  # integron_id -> accumulator

    with open(path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split("\t")
            if len(parts) < 8:
                continue
            integron_id = parts[0]
            contig = parts[1]
            element_type = parts[7]  # type_elt: protein, attC, attI, Pc, Pi

            key = f"{integron_id}:{contig}"
            if key not in current:
                current[key] = {
                    "id": integron_id, "contig": contig,
                    "type": parts[13] if len(parts) > 13 else "unknown",  # CALIN, complete, In0
                    "n_attC": 0, "n_proteins": 0,
                }
            if element_type == "attC":
                current[key]["n_attC"] += 1
            elif element_type == "protein":
                current[key]["n_proteins"] += 1

    for key, info in current.items():
        contig = info["contig"]
        entry = {"id": info["id"], "type": info["type"],
                 "n_attC": info["n_attC"], "n_proteins": info["n_proteins"]}
        integrons.setdefault(contig, []).append(entry)

    return integrons


def load_genomic_islands(path: Path) -> dict[str, list[dict]]:
    """Load IslandPath-DIMOB genomic island predictions.

    Input format: island_id\\tstart\\tend (tab-separated, header line).
    Island IDs are contig names (islands are regions within contigs).
    Returns {contig_name: [{start, end, length}]}.
    """
    islands: dict[str, list[dict]] = {}

    with open(path) as f:
        header = f.readline().strip().split("\t")
        # Detect column layout: may have island_id prefix column
        try:
            contig_col = header.index("contig")
            start_col = header.index("start")
            end_col = header.index("end")
        except ValueError:
            # Fallback: assume contig, start, end
            contig_col, start_col, end_col = 0, 1, 2
        n_cols = max(contig_col, start_col, end_col) + 1
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < n_cols:
                continue
            contig = parts[contig_col]
            start = int(parts[start_col])
            end = int(parts[end_col])
            islands.setdefault(contig, []).append({
                "start": start, "end": end, "length": end - start,
            })

    return islands


def load_macsyfinder(path: Path) -> dict[str, list[dict]]:
    """Load MacSyFinder secretion/conjugation system predictions.

    The hit_id column contains Prokka locus tags (e.g. HFMBNMNI_49108).
    We need a GFF to map locus tags → contigs, but since all hits are from
    the same replicon, we group by sys_id and extract the contig from the
    locus tag position in the assembly.

    Returns {prokka_locus_tag: {sys_id, gene_name, model, wholeness, hit_status}}.
    Also returns a sys_id-keyed summary for later contig mapping.
    """
    # Parse system-level data: group genes by sys_id
    systems: dict[str, dict] = {}  # sys_id -> {model, wholeness, genes: [gene_name]}
    gene_hits: dict[str, dict] = {}  # locus_tag -> {sys_id, gene_name, model, wholeness, status}

    with open(path) as f:
        for line in f:
            if line.startswith("#") or line.startswith("replicon\t") or not line.strip():
                continue
            parts = line.strip().split("\t")
            if len(parts) < 10:
                continue
            locus_tag = parts[1]   # hit_id (Prokka locus tag)
            gene_name = parts[2]   # gene_name (e.g. T4SS_MOBQ)
            model = parts[4]       # model_fqn (e.g. CONJScan/Chromosome/MOB)
            sys_id = parts[5]      # sys_id
            wholeness = parts[6]   # sys_wholeness

            if sys_id not in systems:
                systems[sys_id] = {
                    "model": model, "wholeness": float(wholeness),
                    "genes": [],
                }
            systems[sys_id]["genes"].append(gene_name)

            gene_hits[locus_tag] = {
                "sys_id": sys_id,
                "gene_name": gene_name,
                "model": model,
                "wholeness": float(wholeness),
                "status": parts[8] if len(parts) > 8 else "unknown",
            }

    return gene_hits, systems


def map_macsyfinder_to_contigs(
    gene_hits: dict[str, dict],
    systems: dict[str, dict],
    gff_path: Path,
) -> dict[str, list[dict]]:
    """Map MacSyFinder locus tags to contigs using Prokka GFF.

    Returns {contig_name: [{sys_id, model, wholeness, genes: [gene_name]}]}.
    """
    # Build locus_tag → contig mapping from GFF
    tag_to_contig: dict[str, str] = {}
    with open(gff_path) as f:
        for line in f:
            if line.startswith("##FASTA"):
                break
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9 or parts[2] != "CDS":
                continue
            contig = parts[0]
            attrs = parts[8]
            # Extract ID from attributes: ID=HFMBNMNI_00001;...
            for attr in attrs.split(";"):
                if attr.startswith("ID="):
                    tag_to_contig[attr[3:]] = contig
                    break

    # Map systems to contigs
    contig_systems: dict[str, dict[str, dict]] = {}  # contig -> sys_id -> system_info
    for locus_tag, hit in gene_hits.items():
        contig = tag_to_contig.get(locus_tag)
        if not contig:
            continue
        sys_id = hit["sys_id"]
        if contig not in contig_systems:
            contig_systems[contig] = {}
        if sys_id not in contig_systems[contig]:
            sys_info = systems.get(sys_id, {})
            contig_systems[contig][sys_id] = {
                "sys_id": sys_id,
                "model": sys_info.get("model", hit["model"]),
                "wholeness": sys_info.get("wholeness", hit["wholeness"]),
                "genes": [],
            }
        contig_systems[contig][sys_id]["genes"].append(hit["gene_name"])

    # Flatten: contig -> list of systems
    result: dict[str, list[dict]] = {}
    for contig, sys_dict in contig_systems.items():
        result[contig] = list(sys_dict.values())

    return result


def load_defensefinder(genes_path: Path) -> tuple[dict[str, dict], dict[str, dict]]:
    """Load DefenseFinder gene-level results.

    The hit_id column contains Prokka locus tags (e.g. HFMBNMNI_00542).
    We group genes by sys_id and track system type/subtype.

    Returns (gene_hits, systems):
        gene_hits: {locus_tag: {sys_id, gene_name, type, subtype}}
        systems: {sys_id: {type, subtype, genes: [gene_name]}}
    """
    systems: dict[str, dict] = {}
    gene_hits: dict[str, dict] = {}

    with open(genes_path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            locus_tag = row["hit_id"]
            gene_name = row["gene_name"]
            sys_id = row["sys_id"]
            sys_type = row.get("type", "")
            sys_subtype = row.get("subtype", "")

            if sys_id not in systems:
                systems[sys_id] = {
                    "type": sys_type,
                    "subtype": sys_subtype,
                    "genes": [],
                }
            systems[sys_id]["genes"].append(gene_name)

            gene_hits[locus_tag] = {
                "sys_id": sys_id,
                "gene_name": gene_name,
                "type": sys_type,
                "subtype": sys_subtype,
            }

    return gene_hits, systems


def map_defensefinder_to_contigs(
    gene_hits: dict[str, dict],
    systems: dict[str, dict],
    gff_path: Path,
) -> dict[str, list[dict]]:
    """Map DefenseFinder locus tags to contigs using Prokka GFF.

    Returns {contig_name: [{sys_id, type, subtype, genes: [gene_name]}]}.
    """
    # Build locus_tag -> contig mapping from GFF
    tag_to_contig: dict[str, str] = {}
    with open(gff_path) as f:
        for line in f:
            if line.startswith("##FASTA"):
                break
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9 or parts[2] != "CDS":
                continue
            contig = parts[0]
            attrs = parts[8]
            for attr in attrs.split(";"):
                if attr.startswith("ID="):
                    tag_to_contig[attr[3:]] = contig
                    break

    # Map systems to contigs
    contig_systems: dict[str, dict[str, dict]] = {}
    for locus_tag, hit in gene_hits.items():
        contig = tag_to_contig.get(locus_tag)
        if not contig:
            continue
        sys_id = hit["sys_id"]
        if contig not in contig_systems:
            contig_systems[contig] = {}
        if sys_id not in contig_systems[contig]:
            sys_info = systems.get(sys_id, {})
            contig_systems[contig][sys_id] = {
                "sys_id": sys_id,
                "type": sys_info.get("type", hit["type"]),
                "subtype": sys_info.get("subtype", hit["subtype"]),
                "genes": [],
            }
        contig_systems[contig][sys_id]["genes"].append(hit["gene_name"])

    result: dict[str, list[dict]] = {}
    for contig, sys_dict in contig_systems.items():
        result[contig] = list(sys_dict.values())

    return result


# ---------------------------------------------------------------------------
# Assembly: build all identity cards
# ---------------------------------------------------------------------------

def build_identities(
    tnf_path: Path,
    depths_path: Path,
    assembly_info_path: Path,
    gfa_path: Path,
    binner_paths: dict[str, Path],
    consensus_path: Path,
    summary_path: Path,
    checkm2_path: Optional[Path] = None,
    virus_summary_path: Optional[Path] = None,
    plasmid_summary_path: Optional[Path] = None,
    checkv_quality_path: Optional[Path] = None,
    kaiju_taxonomy_path: Optional[Path] = None,
    integron_path: Optional[Path] = None,
    genomic_island_path: Optional[Path] = None,
    macsyfinder_path: Optional[Path] = None,
    defensefinder_path: Optional[Path] = None,
    prokka_gff_path: Optional[Path] = None,
) -> tuple[dict[str, ContigIdentity], dict[str, CommunityProfile]]:
    """Build identity cards for all contigs and community profiles.

    This is Gathering 1: Self-Knowledge.
    """
    # Load everything
    tnf_data = load_tnf(tnf_path)
    sample_names, contig_lengths, coverages = load_depths(depths_path)
    assembly_info = load_assembly_info(assembly_info_path)
    graph = load_gfa_graph(gfa_path)
    binner_data = load_binner_assignments(binner_paths)
    consensus = load_consensus(consensus_path)
    dastool_summary = load_dastool_summary(summary_path)
    checkm2_data = load_checkm2(checkm2_path) if checkm2_path else {}

    # Load MGE data
    virus_data = load_genomad_viruses(virus_summary_path) if virus_summary_path else {}
    plasmid_data = load_genomad_plasmids(plasmid_summary_path) if plasmid_summary_path else {}
    checkv_data = load_checkv_quality(checkv_quality_path) if checkv_quality_path else {}

    # Load taxonomy
    kaiju_data = load_kaiju_taxonomy(kaiju_taxonomy_path) if kaiju_taxonomy_path else {}

    # Load integrons, genomic islands, secretion systems
    integron_data = load_integrons(integron_path) if integron_path else {}
    island_data = load_genomic_islands(genomic_island_path) if genomic_island_path else {}
    msf_data: dict[str, list[dict]] = {}
    if macsyfinder_path and prokka_gff_path:
        gene_hits, systems = load_macsyfinder(macsyfinder_path)
        msf_data = map_macsyfinder_to_contigs(gene_hits, systems, prokka_gff_path)

    df_data: dict[str, list[dict]] = {}
    if defensefinder_path and prokka_gff_path:
        df_gene_hits, df_systems = load_defensefinder(defensefinder_path)
        df_data = map_defensefinder_to_contigs(df_gene_hits, df_systems, prokka_gff_path)

    # Build identity cards for every contig that has TNF data
    identities: dict[str, ContigIdentity] = {}

    for name, tnf_vec in tnf_data.items():
        length = contig_lengths.get(name, 0)
        cov = coverages.get(name, np.zeros(len(sample_names)))
        gc = float(tnf_vec @ GC_WEIGHTS)

        info = assembly_info.get(name, {})

        # Oracle testimony
        testimony = binner_data.get(name, {})
        voice_strength = sum(1 for v in testimony.values() if v is not None)

        # Community assignment
        community = consensus.get(name)
        membership_type = "core" if community else "unhoused"

        identity = ContigIdentity(
            name=name,
            size=length,
            tnf=tnf_vec,
            coverage=cov,
            gc=gc,
            assembly_coverage=info.get("coverage", 0.0),
            is_circular=info.get("circular", False),
            is_repeat=info.get("repeat", False),
            multiplicity=info.get("multiplicity", 1),
            connections=graph.get(name, []),
            testimony=testimony,
            voice_strength=voice_strength,
            community=community,
            membership_type=membership_type,
        )

        # Annotate with taxonomy
        kaiju = kaiju_data.get(name)
        if kaiju:
            identity.ancestry = kaiju["lineage"]

        # Annotate with MGE data
        vir = virus_data.get(name)
        if vir:
            identity.is_viral = True
            identity.virus_score = vir["virus_score"]
            identity.virus_hallmarks = vir["n_hallmarks"]
            identity.virus_taxonomy = vir["taxonomy"] or None
            identity.mge_type = "virus"

        plas = plasmid_data.get(name)
        if plas:
            identity.is_plasmid = True
            identity.plasmid_score = plas["plasmid_score"]
            identity.plasmid_hallmarks = plas["n_hallmarks"]
            identity.conjugation_genes = plas["conjugation_genes"]
            identity.amr_genes = plas["amr_genes"]
            if not identity.is_viral:  # virus takes precedence
                identity.mge_type = "plasmid"

        cv = checkv_data.get(name)
        if cv:
            identity.is_provirus = cv["provirus"]
            identity.proviral_length = cv["proviral_length"]
            identity.checkv_quality = cv["checkv_quality"]
            identity.viral_genes = cv["viral_genes"]
            identity.host_genes = cv["host_genes"]
            if cv["provirus"]:
                identity.mge_type = "provirus"
                # Proviruses are Travelers — integrated in host but viral in nature
                if community:
                    identity.membership_type = "traveler"

        # Annotate with integron data
        intg = integron_data.get(name)
        if intg:
            identity.integrons = intg
            identity.has_integron = True

        # Annotate with genomic island data
        isld = island_data.get(name)
        if isld:
            identity.genomic_islands = isld
            identity.has_genomic_island = True

        # Annotate with secretion/conjugation system data
        msf = msf_data.get(name)
        if msf:
            identity.secretion_systems = msf
            identity.has_secretion_system = True

        # Annotate with defense system data
        dfs = df_data.get(name)
        if dfs:
            identity.defense_systems = dfs
            identity.has_defense_system = True

        identities[name] = identity

    # Build community profiles
    communities: dict[str, CommunityProfile] = {}

    for bin_name, summary in dastool_summary.items():
        members = [c for c, b in consensus.items() if b == bin_name]
        member_tnf = np.array([identities[c].tnf for c in members if c in identities])
        member_cov = np.array([identities[c].coverage for c in members if c in identities])
        member_gc = [identities[c].gc for c in members if c in identities]

        tnf_centroid = member_tnf.mean(axis=0) if len(member_tnf) > 0 else None
        mean_coverage = member_cov.mean(axis=0) if len(member_cov) > 0 else None

        checkm2 = checkm2_data.get(bin_name, {})

        completeness = summary["completeness"]
        # Determine elder rank
        if completeness >= 95:
            elder_rank = "contigsattva"
        elif completeness >= 90:
            elder_rank = "sage"
        elif completeness >= 70:
            elder_rank = "full"
        elif completeness >= 50:
            elder_rank = "apprentice"
        else:
            elder_rank = "none"

        profile = CommunityProfile(
            name=bin_name,
            source_binner=summary["source"],
            members=members,
            tnf_centroid=tnf_centroid,
            mean_coverage=mean_coverage,
            mean_gc=float(np.mean(member_gc)) if member_gc else 0.0,
            gc_stdev=float(np.std(member_gc)) if member_gc else 0.0,
            total_size=summary["size"],
            n50=summary["n50"],
            completeness=completeness,
            redundancy=summary["redundancy"],
            contamination=checkm2.get("contamination", 0.0),
            checkm2_completeness=checkm2.get("completeness", 0.0),
            elder_rank=elder_rank,
        )
        communities[bin_name] = profile

    return identities, communities


# ---------------------------------------------------------------------------
# Binner agreement seeding
# ---------------------------------------------------------------------------

def seed_communities_from_binner_agreement(
    identities: dict[str, ContigIdentity],
    max_bin_size: int = 500,
    min_members: int = 5,
    min_total_size: int = 100_000,
) -> tuple[dict[str, CommunityProfile], dict[str, str]]:
    """Seed new communities from binner co-assignment (3+ out of N majority vote).

    Contigs placed in the same bin by multiple binners are grouped into new
    seed communities, even if DAS Tool didn't select that bin as consensus.

    Algorithm:
        1. Start at maximum agreement (all binners agree), work down to 3
        2. Build composite keys from binner labels for each combination
        3. Filter: skip mega-bins (>max_bin_size), require ≥min_members and ≥min_total_size
        4. Higher agreement level takes priority (no double-assignment)

    Returns:
        (new_communities, assignments) where assignments maps contig→community_name
    """
    from itertools import combinations

    # Reconstruct binner assignments from testimony
    binner_names = sorted(set().union(
        *(c.testimony.keys() for c in identities.values() if c.testimony)
    ))
    n_binners = len(binner_names)
    if n_binners < 3:
        return {}, {}

    # Compute bin sizes per binner to filter mega-bins
    bin_sizes: dict[tuple[str, str], int] = {}
    for c in identities.values():
        for binner, label in c.testimony.items():
            if label is not None:
                key = (binner, label)
                bin_sizes[key] = bin_sizes.get(key, 0) + 1

    mega_bins = {k for k, v in bin_sizes.items() if v > max_bin_size}

    # Only consider unhoused contigs
    unhoused = {name for name, c in identities.items() if c.community is None}

    assigned: dict[str, str] = {}
    community_members: dict[str, list[str]] = {}
    community_agreement: dict[str, int] = {}
    comm_counter = 0

    # Work from highest agreement to lowest
    for k in range(n_binners, 2, -1):  # e.g., 5, 4, 3
        remaining = unhoused - set(assigned.keys())
        if not remaining:
            break

        for combo in combinations(binner_names, k):
            # Build composite key for each remaining contig
            groups: dict[str, list[str]] = {}
            for contig in remaining:
                if contig in assigned:
                    continue
                testimony = identities[contig].testimony
                labels = []
                skip = False
                for binner in combo:
                    label = testimony.get(binner)
                    if label is None:
                        skip = True
                        break
                    if (binner, label) in mega_bins:
                        skip = True
                        break
                    labels.append(f"{binner}:{label}")
                if skip:
                    continue
                key = "|".join(labels)
                groups.setdefault(key, []).append(contig)

            # Create communities from qualifying groups
            for key, members in groups.items():
                new_members = [c for c in members if c not in assigned]
                if len(new_members) < min_members:
                    continue
                total_size = sum(
                    identities[c].size for c in new_members if c in identities
                )
                if total_size < min_total_size:
                    continue

                comm_name = f"binner-agree-{k}of{n_binners}_{comm_counter}"
                community_members[comm_name] = new_members
                community_agreement[comm_name] = k
                comm_counter += 1
                for c in new_members:
                    assigned[c] = comm_name

    # Build CommunityProfile objects for new communities
    new_communities: dict[str, CommunityProfile] = {}
    for comm_name, members in community_members.items():
        member_ids = [identities[c] for c in members if c in identities]
        if not member_ids:
            continue

        member_tnf = np.array([c.tnf for c in member_ids])
        member_cov = np.array([c.coverage for c in member_ids])
        member_gc = [c.gc for c in member_ids]

        tnf_centroid = member_tnf.mean(axis=0)
        mean_coverage = member_cov.mean(axis=0)
        total_size = sum(c.size for c in member_ids)

        # Compute N50
        sizes_sorted = sorted([c.size for c in member_ids], reverse=True)
        cumsum = 0
        n50 = 0
        for s in sizes_sorted:
            cumsum += s
            if cumsum >= total_size / 2:
                n50 = s
                break

        agreement = community_agreement[comm_name]

        profile = CommunityProfile(
            name=comm_name,
            source_binner=f"binner-agreement-{agreement}of{n_binners}",
            members=members,
            tnf_centroid=tnf_centroid,
            mean_coverage=mean_coverage,
            mean_gc=float(np.mean(member_gc)),
            gc_stdev=float(np.std(member_gc)),
            total_size=total_size,
            n50=n50,
            elder_rank="none",
        )
        new_communities[comm_name] = profile

    return new_communities, assigned
