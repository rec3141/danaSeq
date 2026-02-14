#!/usr/bin/env python3
"""Gathering 2: The Conversations (tool-use architecture).

Each contig is an agent that discovers itself through tool calls.
The LLM investigates by calling tools like who_am_i(), how_do_i_resonate_with(),
what_is_this_community() — each returns actual data, not pre-digested metrics.

Supports OpenAI-compatible servers (LM Studio, ollama, vLLM) and Anthropic API.

Usage:
    nclb_converse.py --results /path/to/results [--output proposals.json]
    nclb_converse.py --results /path/to/results --backend anthropic
    nclb_converse.py --results /path/to/results --base-url http://host:1234/v1
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "lib"))

from nclb.identity import (
    build_identities, load_gfa_graph, ContigIdentity, CommunityProfile,
)
from nclb.voices import (
    ContigToolkit, CONTIG_TOOLS_OPENAI, parse_json_response,
)
from nclb.valence import contig_valence, tnf_coherence, coverage_coherence


# ---------------------------------------------------------------------------
# System prompts (brief — the LLM discovers context through tools)
# ---------------------------------------------------------------------------

ROUND1_SYSTEM = """You are the voice of metagenome-assembled genome communities.

You have tools to investigate any contig or community. Use them to understand
why members are uneasy before deciding whether to release them.

IMPORTANT biological rules:
- PROVIRUS contigs are integrated phages — composition divergence is EXPECTED.
  Do NOT release a provirus just because its GC/TNF differs from the host.
- PLASMID contigs are mobile elements the organism carries — not contamination.
- Only release contigs with MULTIPLE lines of evidence: zero coverage AND
  graph isolation AND conflicting oracle testimony AND no MGE explanation.

After investigating, respond with JSON (no commentary):
{"community": "name", "assessment": "narrative", "release": [{"contig": "name", "reason": "evidence"}], "concerns": []}"""

ROUND2_SYSTEM = """You speak for unhoused contigs seeking community.

You have tools to investigate each contig's identity, graph connections,
and resonance with candidate communities. Call them to build evidence.

Strategy:
1. Call who_am_i() to learn the contig's full identity
2. Call find_graph_connections() to see which communities it's physically linked to
3. Call how_do_i_resonate_with() for promising communities
4. Compare coverage profiles — do they track together across samples?

For VIRAL/PLASMID/PROVIRUS contigs: graph connections to a host community
are strong evidence of integration. Place them with the host.

After investigating, respond with JSON (no commentary):
{"decisions": [{"contig": "name", "action": "join|wait|wander", "community": "name_or_null", "evidence": "specific signals", "valence": 0.0}]}"""

ROUND3_SYSTEM = """You evaluate candidate new communities formed from voiceless contigs.

Signs of a real community: TNF coherence >0.9, synchronized coverage across
samples, reasonable genome size, consistent ancestry.

Respond with JSON (no commentary):
{"evaluations": [{"cluster_id": 0, "verdict": "accept|reject|uncertain", "reason": "assessment", "suggested_name": "name"}]}"""


# ---------------------------------------------------------------------------
# Tool-use conversation loop
# ---------------------------------------------------------------------------

def run_tool_conversation(
    client,
    model: str,
    system: str,
    user_prompt: str,
    toolkit: ContigToolkit,
    max_rounds: int = 15,
    log_fn=None,
) -> dict:
    """Run a tool-use conversation loop with an OpenAI-compatible API.

    The LLM calls tools to investigate, we execute and return results,
    until it produces a final text response.

    When approaching the round limit, injects a nudge asking the model to
    produce its final JSON decision, then forces a tool-free final round.
    """
    messages = [
        {"role": "system", "content": system},
        {"role": "user", "content": user_prompt},
    ]

    tool_calls_total = 0
    nudged = False
    # Reserve last 3 rounds for winding down: nudge at -3, remind at -2,
    # force text at -1
    wind_down_start = max(max_rounds - 3, max_rounds // 2)

    for round_num in range(max_rounds):
        # Wind-down phase: nudge the model to produce decisions
        if round_num == wind_down_start and not nudged:
            messages.append({
                "role": "user",
                "content": (
                    "You have used most of your investigation budget. "
                    "Please finish your investigation and respond with "
                    "your final JSON decision now. Remember the required "
                    "format from your instructions."
                ),
            })
            nudged = True
        elif round_num == max_rounds - 2 and nudged:
            messages.append({
                "role": "user",
                "content": (
                    "This is your LAST chance to investigate. On the next "
                    "round you must produce your JSON response. Respond "
                    "with the JSON now if you are ready."
                ),
            })

        # On the final round, don't offer tools — force a text response
        is_final_round = (round_num == max_rounds - 1)
        call_kwargs = dict(
            model=model,
            messages=messages,
            max_tokens=4096,
            temperature=0.3,
        )
        if not is_final_round:
            call_kwargs["tools"] = CONTIG_TOOLS_OPENAI

        try:
            response = client.chat.completions.create(**call_kwargs)
        except Exception as e:
            if log_fn:
                log_fn(f"    [ERROR] API call failed: {e}")
            return {"error": str(e), "tool_calls": tool_calls_total}

        choice = response.choices[0]
        msg = choice.message

        # Check for tool calls
        if msg.tool_calls:
            # Append the assistant message with tool calls
            messages.append({
                "role": "assistant",
                "content": msg.content or "",
                "tool_calls": [
                    {
                        "id": tc.id,
                        "type": "function",
                        "function": {
                            "name": tc.function.name,
                            "arguments": tc.function.arguments,
                        },
                    }
                    for tc in msg.tool_calls
                ],
            })

            # Execute each tool call
            for tc in msg.tool_calls:
                tool_calls_total += 1
                try:
                    args = json.loads(tc.function.arguments)
                    result = toolkit.dispatch(tc.function.name, args)
                except Exception as e:
                    result = {"error": f"Tool execution failed: {e}"}

                messages.append({
                    "role": "tool",
                    "tool_call_id": tc.id,
                    "content": json.dumps(result, default=str),
                })

            if log_fn and tool_calls_total % 5 == 0:
                log_fn(f"    ({tool_calls_total} tool calls so far)")
            continue

        # No tool calls — this is the final response
        text = msg.content or ""
        if text:
            try:
                parsed = parse_json_response(text)
                parsed["_tool_calls"] = tool_calls_total
                return parsed
            except (json.JSONDecodeError, ValueError):
                # Model produced text but not valid JSON (e.g. XML tool
                # calls when tools were removed). If this is the forced
                # final round, retry once asking for pure JSON.
                if is_final_round:
                    messages.append({"role": "assistant", "content": text})
                    messages.append({
                        "role": "user",
                        "content": (
                            "That was not valid JSON. Please respond with "
                            "ONLY a JSON object matching the required format. "
                            "No tool calls, no XML, no commentary — just the "
                            "JSON object."
                        ),
                    })
                    try:
                        retry = client.chat.completions.create(
                            model=model, messages=messages,
                            max_tokens=4096, temperature=0.1,
                        )
                        retry_text = retry.choices[0].message.content or ""
                        parsed = parse_json_response(retry_text)
                        parsed["_tool_calls"] = tool_calls_total
                        return parsed
                    except Exception:
                        pass
                return {"raw_response": text, "_tool_calls": tool_calls_total}

        # Empty response — shouldn't happen
        if log_fn:
            log_fn(f"    [WARNING] Empty response at round {round_num}")
        break

    # Exhausted rounds — try one final forced-text round without tools
    if log_fn:
        log_fn(f"    [INFO] Forcing final response after {tool_calls_total} tool calls")
    messages.append({
        "role": "user",
        "content": (
            "Investigation complete. You MUST respond now with your final "
            "JSON decision. No more tool calls. Produce ONLY the JSON object "
            "matching the required format from your instructions."
        ),
    })
    try:
        final = client.chat.completions.create(
            model=model, messages=messages,
            max_tokens=4096, temperature=0.1,
        )
        text = final.choices[0].message.content or ""
        if text:
            try:
                parsed = parse_json_response(text)
                parsed["_tool_calls"] = tool_calls_total
                return parsed
            except (json.JSONDecodeError, ValueError):
                return {"raw_response": text, "_tool_calls": tool_calls_total}
    except Exception as e:
        if log_fn:
            log_fn(f"    [ERROR] Final forced response failed: {e}")

    return {"error": "max_rounds_exhausted", "_tool_calls": tool_calls_total}


# ---------------------------------------------------------------------------
# Round builders (brief prompts — tools provide the data)
# ---------------------------------------------------------------------------

def round1_prompt(comm_name: str, comm_data: dict, uneasy_names: list[str]) -> str:
    """Brief prompt for Round 1 — the LLM investigates via tools."""
    uneasy_section = "None — all members have positive valence."
    if uneasy_names:
        uneasy_section = ", ".join(uneasy_names)

    return f"""Examine community {comm_name} [{comm_data.get('elder_rank', 'none')}].

Quick overview:
  Members: {len(comm_data.get('members', []))}
  Size: {comm_data['total_size']:,} bp
  Wholeness: {comm_data['completeness']:.1f}% | Redundancy: {comm_data['redundancy']:.1f}%
  Mean valence: {comm_data['mean_valence']:+.3f} | Min valence: {comm_data['min_valence']:+.3f}

Uneasy members (negative valence): {uneasy_section}

Investigate each uneasy member using who_am_i() and how_do_i_resonate_with().
Check their coverage, graph connections, and oracle testimony.
Then decide which should be released."""


def round2_prompt(contig_names: list[str]) -> str:
    """Brief prompt for Round 2 — the LLM investigates via tools."""
    names = ", ".join(contig_names)
    return f"""You speak for {len(contig_names)} unhoused contigs seeking community: {names}

For each contig:
1. Call who_am_i() to learn its full identity
2. Call find_graph_connections() to see graph links to communities
3. Call how_do_i_resonate_with() for the most promising community
4. Decide: join, wait, or wander

Investigate each contig and recommend placement."""


def round3_prompt(clusters: list[dict]) -> str:
    """Prompt for Round 3 — data-in-prompt is fine for cluster summaries."""
    lines = []
    for cl in clusters:
        lines.append(
            f"  Cluster {cl['id']}: {cl['n_members']} contigs, "
            f"{cl['total_size']:,}bp, TNF coherence={cl['tnf_coherence']:.3f}, "
            f"Coverage correlation={cl.get('coverage_correlation', 0):.3f}, "
            f"Mean GC={cl.get('mean_gc', 0):.3f}"
        )
    return f"""{len(clusters)} clusters emerged from voiceless contigs.

CANDIDATES:
{chr(10).join(lines)}

Evaluate each: real organism or noise?"""


# ---------------------------------------------------------------------------
# HDBSCAN clustering for Round 3
# ---------------------------------------------------------------------------

def cluster_voiceless(identities: dict) -> list[dict]:
    """Cluster truly voiceless contigs using HDBSCAN."""
    try:
        from sklearn.cluster import HDBSCAN
        from sklearn.preprocessing import StandardScaler
    except ImportError:
        return []

    import numpy as np

    voiceless = []
    for c in identities.values():
        if c.community is None and c.voice_strength == 0 and not c.connections:
            if c.size >= 5000:
                voiceless.append(c)

    if len(voiceless) < 10:
        return []

    tnf_matrix = np.array([c.tnf for c in voiceless])
    cov_matrix = np.array([c.coverage for c in voiceless])

    scaler_tnf = StandardScaler()
    scaler_cov = StandardScaler()
    tnf_scaled = scaler_tnf.fit_transform(tnf_matrix)
    cov_scaled = scaler_cov.fit_transform(cov_matrix) if cov_matrix.shape[1] > 0 else np.zeros((len(voiceless), 0))

    features = np.hstack([tnf_scaled * 2.0, cov_scaled])
    clusterer = HDBSCAN(min_cluster_size=5, min_samples=3)
    labels = clusterer.fit_predict(features)

    clusters = []
    for label in set(labels):
        if label == -1:
            continue
        members = [voiceless[i] for i, l in enumerate(labels) if l == label]
        total_size = sum(c.size for c in members)
        if total_size < 100000:
            continue
        clusters.append({
            "id": int(label),
            "n_members": len(members),
            "total_size": total_size,
            "member_names": [c.name for c in members],
            "tnf_coherence": float(tnf_coherence(members)),
            "coverage_correlation": float(coverage_coherence(members)),
            "mean_gc": float(np.mean([c.gc for c in members])),
        })
    return clusters


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="NCLB Gathering 2: Tool-use conversations"
    )
    parser.add_argument(
        "--results", "-r", type=Path, required=True,
        help="Path to Nextflow results directory"
    )
    parser.add_argument(
        "--gathering", "-g", type=Path, default=None,
        help="Path to gathering.json (for pre-computed resonance data)"
    )
    parser.add_argument(
        "--output", "-o", type=Path, default=None,
        help="Output proposals.json path"
    )
    parser.add_argument(
        "--base-url", type=str, default="http://10.151.30.147:1234/v1",
        help="Base URL for OpenAI-compatible server"
    )
    parser.add_argument(
        "--model", "-m", type=str, default=None,
        help="Model name (default: auto-detect from server)"
    )
    parser.add_argument(
        "--max-tool-rounds", type=int, default=15,
        help="Max tool-use rounds per conversation (default: 15)"
    )
    parser.add_argument(
        "--batch-size", type=int, default=5,
        help="Contigs per Round 2 conversation (default: 5)"
    )
    parser.add_argument(
        "--max-round2", type=int, default=100,
        help="Max Round 2 contigs to process (default: 100)"
    )
    parser.add_argument(
        "--quiet", "-q", action="store_true",
        help="Suppress progress output"
    )
    args = parser.parse_args()

    results = args.results
    gathering_path = args.gathering or (results / "binning" / "nclb" / "gathering.json")
    output_path = args.output or (results / "binning" / "nclb" / "proposals.json")

    def log(msg: str):
        if not args.quiet:
            print(msg, file=sys.stderr, flush=True)

    t0 = time.time()

    # --- Load full identity data ---
    log("[INFO] Loading identity data from results directory...")
    assembly_dir = results / "assembly"
    mapping_dir = results / "mapping"
    binning_dir = results / "binning"

    binner_paths = {}
    for name in ["semibin", "metabat", "maxbin", "lorbin", "comebin"]:
        tsv = binning_dir / name / "contig_bins.tsv"
        if tsv.exists():
            binner_paths[name] = tsv

    checkm2_path = binning_dir / "checkm2" / "quality_report.tsv"
    if not checkm2_path.exists():
        checkm2_path = None

    # MGE data
    mge_dir = results / "mge"
    virus_path = mge_dir / "genomad" / "virus_summary.tsv"
    plasmid_path = mge_dir / "genomad" / "plasmid_summary.tsv"
    checkv_path = mge_dir / "checkv" / "quality_summary.tsv"

    # Taxonomy
    kaiju_path = results / "taxonomy" / "kaiju" / "kaiju_contigs.tsv"

    # Integrons, genomic islands, secretion systems
    integron_path = mge_dir / "integrons" / "integrons.tsv"
    island_path = mge_dir / "genomic_islands" / "genomic_islands.tsv"
    msf_path = mge_dir / "macsyfinder" / "all_systems.tsv"
    df_path = mge_dir / "defensefinder" / "genes.tsv"
    prokka_gff_path = None
    annotation_dir = results / "annotation" / "prokka"
    if annotation_dir.exists():
        gffs = sorted(annotation_dir.glob("*/*.gff")) + sorted(annotation_dir.glob("*.gff"))
        if gffs:
            prokka_gff_path = gffs[-1]

    identities, communities = build_identities(
        tnf_path=assembly_dir / "tnf.tsv",
        depths_path=mapping_dir / "depths.txt",
        assembly_info_path=assembly_dir / "assembly_info.txt",
        gfa_path=assembly_dir / "assembly_graph.gfa",
        binner_paths=binner_paths,
        consensus_path=binning_dir / "dastool" / "contig2bin.tsv",
        summary_path=binning_dir / "dastool" / "summary.tsv",
        checkm2_path=checkm2_path,
        virus_summary_path=virus_path if virus_path.exists() else None,
        plasmid_summary_path=plasmid_path if plasmid_path.exists() else None,
        checkv_quality_path=checkv_path if checkv_path.exists() else None,
        kaiju_taxonomy_path=kaiju_path if kaiju_path.exists() else None,
        integron_path=integron_path if integron_path.exists() else None,
        genomic_island_path=island_path if island_path.exists() else None,
        macsyfinder_path=msf_path if msf_path.exists() else None,
        defensefinder_path=df_path if df_path.exists() else None,
        prokka_gff_path=prokka_gff_path,
    )
    adjacency = load_gfa_graph(assembly_dir / "assembly_graph.gfa")
    log(f"[INFO] Loaded {len(identities):,} contigs, {len(communities)} DAS Tool communities")

    # Seed additional communities from binner agreement
    from nclb.identity import seed_communities_from_binner_agreement
    new_comms, binner_assigned = seed_communities_from_binner_agreement(identities)
    if new_comms:
        for contig, comm_name in binner_assigned.items():
            identities[contig].community = comm_name
            identities[contig].membership_type = "core"
        communities.update(new_comms)
        log(f"[INFO] Seeded {len(new_comms)} binner-agreement communities ({len(binner_assigned):,} contigs)")
    log(f"[INFO] Total: {len(identities):,} contigs, {len(communities)} communities")

    # Compute valence for housed contigs
    from nclb.valence import contig_valence as cv
    for name, c in identities.items():
        if c.community and c.community in communities:
            c.valence = cv(c, communities[c.community], adjacency)

    # Compute community harmony
    from nclb.valence import community_harmony
    from nclb.graph import graph_connectivity
    for comm in communities.values():
        members = [identities[n] for n in comm.members if n in identities]
        comm.tnf_coherence = tnf_coherence(members)
        comm.coverage_correlation = coverage_coherence(members)
        comm.graph_connectivity = graph_connectivity(comm.members, adjacency)

    # Load landscape data (UMAP coordinates) if available
    landscape_path = binning_dir / "nclb" / "landscape.json"
    if landscape_path.exists():
        from nclb.landscape import load_landscape
        landscape_data = load_landscape(landscape_path)
        for name, lr in landscape_data.items():
            if name in identities:
                identities[name].landscape_x = lr.x
                identities[name].landscape_y = lr.y
                identities[name].landscape_cluster = lr.cluster
                identities[name].landscape_cluster_cx = lr.cluster_cx
                identities[name].landscape_cluster_cy = lr.cluster_cy
        log(f"[INFO] Loaded landscape data for {len(landscape_data):,} contigs")

    # Load gathering.json for pre-computed resonance and uneasy member data
    gathering = None
    if gathering_path.exists():
        with open(gathering_path) as f:
            gathering = json.load(f)
        log(f"[INFO] Loaded gathering data from {gathering_path}")

    # --- Create toolkit ---
    toolkit = ContigToolkit(identities, communities, adjacency)

    # --- Configure OpenAI client ---
    from openai import OpenAI
    base_url = args.base_url
    model = args.model

    if not model:
        try:
            import urllib.request
            resp = urllib.request.urlopen(f"{base_url}/models", timeout=5)
            models_data = json.loads(resp.read())
            if models_data.get("data"):
                model = models_data["data"][0]["id"]
                log(f"[INFO] Auto-detected model: {model}")
        except Exception as e:
            log(f"[WARNING] Could not auto-detect model: {e}")
            model = "local-model"

    client = OpenAI(base_url=base_url, api_key="lm-studio")
    log(f"[INFO] Server: {base_url}, Model: {model}")
    log("")

    all_proposals = []

    # =====================================================================
    # Round 1: Community Health Check (tool-use)
    # =====================================================================
    log("=" * 70)
    log("ROUND 1: COMMUNITY HEALTH CHECK")
    log("=" * 70)

    # Build community data dicts for prompts
    comm_data_map = {}
    uneasy_map = {}
    for comm_name, comm in communities.items():
        harmony = community_harmony(comm, identities, adjacency)
        comm_data_map[comm_name] = {
            "name": comm_name,
            "elder_rank": comm.elder_rank,
            "members": comm.members,
            "total_size": comm.total_size,
            "completeness": comm.completeness,
            "redundancy": comm.redundancy,
            "tnf_coherence": comm.tnf_coherence,
            "coverage_correlation": comm.coverage_correlation,
            "graph_connectivity": comm.graph_connectivity,
            "mean_valence": harmony["mean_valence"],
            "min_valence": harmony["min_valence"],
        }
        # Find uneasy members
        uneasy_names = []
        for m in comm.members:
            c = identities.get(m)
            if c and c.valence < 0:
                uneasy_names.append(m)
        uneasy_map[comm_name] = uneasy_names

    for comm_name in sorted(communities.keys()):
        uneasy = uneasy_map[comm_name]
        comm_data = comm_data_map[comm_name]

        if not uneasy:
            log(f"  {comm_name}: 0 uneasy, skipping")
            all_proposals.append({
                "round": 1, "community": comm_name,
                "result": {"community": comm_name, "assessment": "All members content", "release": []},
            })
            continue

        prompt = round1_prompt(comm_name, comm_data, uneasy)
        try:
            result = run_tool_conversation(
                client, model, ROUND1_SYSTEM, prompt, toolkit,
                max_rounds=args.max_tool_rounds, log_fn=log,
            )
            n_releases = len(result.get("release", []))
            n_tc = result.pop("_tool_calls", 0)
            log(f"  {comm_name}: {n_releases} releases ({n_tc} tool calls)")
            all_proposals.append({
                "round": 1, "community": comm_name, "result": result,
            })
        except Exception as e:
            log(f"  [ERROR] {comm_name}: {e}")

    r1_releases = sum(
        len(p.get("result", {}).get("release", []))
        for p in all_proposals if p["round"] == 1
    )
    log(f"\nRound 1 complete: {r1_releases} releases from {len(communities)} communities")

    # =====================================================================
    # Round 2: Unhoused Contigs Speak (tool-use)
    # =====================================================================
    log("")
    log("=" * 70)
    log("ROUND 2: UNHOUSED CONTIGS SPEAK")
    log("=" * 70)

    # Select unhoused contigs with voice (sorted by voice strength, then size)
    unhoused = [
        c for c in identities.values()
        if c.community is None and c.voice_strength >= 2
    ]
    unhoused.sort(key=lambda c: (-c.voice_strength, -c.size))
    unhoused = unhoused[:args.max_round2]
    log(f"  Processing {len(unhoused)} unhoused contigs (voice >= 2)")

    # Batch contigs for conversation
    batches = []
    for i in range(0, len(unhoused), args.batch_size):
        batches.append([c.name for c in unhoused[i:i + args.batch_size]])

    for batch_idx, batch_names in enumerate(batches):
        prompt = round2_prompt(batch_names)
        try:
            result = run_tool_conversation(
                client, model, ROUND2_SYSTEM, prompt, toolkit,
                max_rounds=args.max_tool_rounds, log_fn=log,
            )
            n_decisions = len(result.get("decisions", []))
            n_tc = result.pop("_tool_calls", 0)
            n_joins = sum(1 for d in result.get("decisions", []) if d.get("action") == "join")
            log(f"  Batch {batch_idx+1}/{len(batches)}: {n_decisions} decisions ({n_joins} joins, {n_tc} tool calls)")
            all_proposals.append({
                "round": 2, "batch": batch_idx, "result": result,
            })
        except Exception as e:
            log(f"  [ERROR] Batch {batch_idx+1}: {e}")

    r2_decisions = sum(
        len(p.get("result", {}).get("decisions", []))
        for p in all_proposals if p["round"] == 2
    )
    r2_joins = sum(
        1 for p in all_proposals if p["round"] == 2
        for d in p.get("result", {}).get("decisions", [])
        if d.get("action") == "join"
    )
    log(f"\nRound 2 complete: {r2_decisions} decisions ({r2_joins} joins)")

    # =====================================================================
    # Round 3: Voiceless Clusters
    # =====================================================================
    log("")
    log("=" * 70)
    log("ROUND 3: THE VOICELESS SPEAK")
    log("=" * 70)

    r3_proposals = []
    clusters = cluster_voiceless(identities)
    if clusters:
        log(f"  Found {len(clusters)} voiceless clusters")
        prompt = round3_prompt(clusters)
        try:
            result = run_tool_conversation(
                client, model, ROUND3_SYSTEM, prompt, toolkit,
                max_rounds=3, log_fn=log,
            )
            result.pop("_tool_calls", None)
            r3_proposals.append({"round": 3, "clusters": clusters, "result": result})
        except Exception as e:
            log(f"  [ERROR] Round 3: {e}")
    else:
        log("  No viable voiceless clusters found")
    all_proposals.extend(r3_proposals)

    # --- Summary ---
    summary = {
        "round1": {"n_communities": len(communities), "n_releases": r1_releases},
        "round2": {"n_batches": len(batches), "n_decisions": r2_decisions, "n_joins": r2_joins},
        "round3": {
            "n_clusters": sum(len(p.get("clusters", [])) for p in r3_proposals),
            "n_accepted": sum(
                1 for p in r3_proposals
                for e in p.get("result", {}).get("evaluations", [])
                if e.get("verdict") == "accept"
            ),
        },
    }

    # --- Write output ---
    output = {
        "version": "0.3.0",
        "timestamp": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "backend": "openai",
        "model": model,
        "architecture": "tool-use",
        "summary": summary,
        "proposals": all_proposals,
    }

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w") as f:
        json.dump(output, f, indent=2, default=str)

    elapsed = time.time() - t0

    log("")
    log("=" * 70)
    log("GATHERING 2: COMPLETE")
    log("=" * 70)
    log(f"  Round 1: {r1_releases} releases")
    log(f"  Round 2: {r2_joins} joins out of {r2_decisions} decisions")
    log(f"  Round 3: {summary['round3']['n_accepted']} new communities")
    log(f"  Output: {output_path}")
    log(f"  Time:   {elapsed:.1f}s")
    log("=" * 70)


if __name__ == "__main__":
    main()
