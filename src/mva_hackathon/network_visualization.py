"""Phase 2: Protein-Protein Interaction (PPI) and Pathway Network Visualization.

Generates reproducible, deterministic interactive network visualizations (Plotly HTML)
and structural topological summaries (CSV):
- 09_phase2_network_nodes.csv: Complete node-level metrics (degree, component rank, clustering, MVA link)
- 10_phase2_component_summary.csv: Connected component ranking, sizes, densities, and MVA gene membership
- 11_phase2_largest_component.html: Interactive network graph for the largest connected component (Component 1)
- 12_phase2_mva_focus_network.html: Interactive network graph for the MVA-focused candidate component (Component 2, containing BUB1B and MAD1L1)
- 13_phase2_pathway_enrichment.html: Interactive pathway enrichment bubble chart (Reactome and g:Profiler)

Implementation Details:
- Uses a deterministic pure-Python Fruchterman-Reingold spring layout algorithm with fixed random seeds.
  (Avoids adding heavyweight C-graph dependencies like networkx or graph-tool while ensuring 100% reproducible node coordinates).
- Visualizes focused subgraphs to avoid uninterpretable hairballs of ~2,000 nodes.
- Uses include_plotlyjs='cdn' to generate lightweight HTML artifacts (< 200 KB each), complying with the < 8 MB git file size policy.
- Strictly non-diagnostic: interactive tooltips, annotations, and legends convey research exploration caveats.
"""

from __future__ import annotations

from collections import deque
from dataclasses import dataclass, field
import json
import logging
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np
import pandas as pd
import plotly.graph_objects as go

logger = logging.getLogger("mva_hackathon.network_visualization")

from mva_hackathon.mva_gene_sets import MVA_CANONICAL_GENES, MVA_EXTENDED_GENES


def fruchterman_reingold_layout(
    nodes: Sequence[str],
    edges: Sequence[tuple[str, str]],
    seed: int = 42,
    iterations: int = 150,
    k: float | None = None,
) -> dict[str, tuple[float, float]]:
    """Compute a deterministic 2D spring layout using Fruchterman-Reingold.

    Parameters:
        nodes: Sequence of node identifier strings.
        edges: Sequence of (node_u, node_v) pairs.
        seed: Random seed for initial placement (fixed for exact reproducibility).
        iterations: Number of cooling iterations.
        k: Optimal spring distance. If None, k = 1.0 / sqrt(len(nodes)).

    Returns:
        Mapping from node identifier to (x, y) float coordinates in [-1, 1].
    """
    n = len(nodes)
    if n == 0:
        return {}
    if n == 1:
        return {nodes[0]: (0.0, 0.0)}

    rng = np.random.RandomState(seed)
    pos = {node: rng.uniform(-1.0, 1.0, 2) for node in nodes}

    if k is None:
        k = 1.0 / np.sqrt(n)

    t = 1.0
    dt = t / (iterations + 1)

    for it in range(iterations):
        disp = {node: np.zeros(2) for node in nodes}

        # 1. Repulsive forces between all node pairs
        for i in range(n):
            u = nodes[i]
            for j in range(i + 1, n):
                v = nodes[j]
                delta = pos[u] - pos[v]
                dist = np.linalg.norm(delta)
                if dist < 1e-4:
                    delta = rng.uniform(-0.01, 0.01, 2)
                    dist = np.linalg.norm(delta)
                rep = (k * k / dist) * (delta / dist)
                disp[u] += rep
                disp[v] -= rep

        # 2. Attractive forces along existing edges
        for u, v in edges:
            if u in pos and v in pos and u != v:
                delta = pos[u] - pos[v]
                dist = np.linalg.norm(delta)
                if dist < 1e-4:
                    continue
                attr = (dist * dist / k) * (delta / dist)
                disp[u] -= attr
                disp[v] += attr

        # 3. Displace positions with temperature cooling
        for node in nodes:
            d = disp[node]
            d_norm = np.linalg.norm(d)
            if d_norm > 0:
                pos[node] += (d / d_norm) * min(d_norm, t)

        t -= dt

    return {node: (float(pos[node][0]), float(pos[node][1])) for node in nodes}


class NetworkVisualizer:
    """End-to-end visualizer for Phase 2 networks and pathway enrichments."""

    def __init__(self, seed: int = 42) -> None:
        self.seed = seed

    def build_network_node_table(
        self,
        nodes_df: pd.DataFrame,
        edges_df: pd.DataFrame,
    ) -> tuple[pd.DataFrame, pd.DataFrame, list[set[str]]]:
        """Compute comprehensive node table and component ranking summary.

        Returns:
            - network_nodes_df: All evaluated nodes (including singletons) with component ranks and degrees
            - component_summary_df: Connected components ranked by size with density and MVA annotations
            - components: List of sets of node names corresponding to ranked components
        """
        # Build adjacency graph from edges
        adj: dict[str, set[str]] = {}
        for _, row in edges_df.iterrows():
            u = str(row["PREFERRED_NAME_A"]).strip().upper()
            v = str(row["PREFERRED_NAME_B"]).strip().upper()
            if u and v and u != v:
                adj.setdefault(u, set()).add(v)
                adj.setdefault(v, set()).add(u)

        # Connected component traversal on edge subgraph
        visited = set()
        components: list[set[str]] = []
        for node in sorted(adj.keys()):
            if node not in visited:
                comp: set[str] = set()
                q = deque([node])
                visited.add(node)
                while q:
                    curr = q.popleft()
                    comp.add(curr)
                    for nbr in adj[curr]:
                        if nbr not in visited:
                            visited.add(nbr)
                            q.append(nbr)
                components.append(comp)

        # Sort components descending by size, then by smallest node name
        components.sort(key=lambda c: (len(c), sorted(list(c))[0]), reverse=True)

        # Build component rank mapping
        node_to_comp_rank: dict[str, int] = {}
        for rank, comp in enumerate(components, start=1):
            for n in comp:
                node_to_comp_rank[n] = rank

        # Component summary records
        comp_records = []
        for rank, comp in enumerate(components, start=1):
            comp_nodes = sorted(list(comp))
            comp_edges = edges_df[
                edges_df["PREFERRED_NAME_A"].isin(comp) & edges_df["PREFERRED_NAME_B"].isin(comp)
            ]
            n_nodes = len(comp)
            n_edges = len(comp_edges)
            density = (2.0 * n_edges) / (n_nodes * (n_nodes - 1)) if n_nodes > 1 else 0.0
            mva_hits = [g for g in comp_nodes if g in MVA_CANONICAL_GENES or g in MVA_EXTENDED_GENES]
            comp_records.append({
                "COMPONENT_RANK": rank,
                "COMPONENT_SIZE": n_nodes,
                "EDGE_COUNT": n_edges,
                "DENSITY": round(density, 4),
                "HAS_MVA_KNOWN_GENE": len(mva_hits) > 0,
                "MVA_GENES": ";".join(mva_hits) if mva_hits else "None",
                "GENE_LIST": ";".join(comp_nodes),
            })

        comp_summary_df = pd.DataFrame(comp_records)

        # Synthesize network nodes table
        # Merge input metadata from nodes_df (07_phase2_network_gene_context.csv)
        node_meta_by_sym: dict[str, dict[str, Any]] = {}
        if not nodes_df.empty:
            for _, r in nodes_df.iterrows():
                sym = str(r.get("HGNC_SYMBOL", "")).strip().upper()
                if sym:
                    node_meta_by_sym[sym] = r.to_dict()

        all_node_names = sorted(list(set(node_meta_by_sym.keys()).union(set(adj.keys()))))
        node_records = []
        for sym in all_node_names:
            meta = node_meta_by_sym.get(sym, {})
            deg = len(adj.get(sym, set()))
            comp_rank = node_to_comp_rank.get(sym, 0)  # 0 indicates isolated / singleton
            is_mva_canonical = sym in MVA_CANONICAL_GENES
            is_mva_extended = sym in MVA_EXTENDED_GENES

            node_records.append({
                "HGNC_SYMBOL": sym,
                "INPUT_SYMBOL": meta.get("INPUT_SYMBOL", sym),
                "NETWORK_DEGREE": deg,
                "COMPONENT_RANK": comp_rank,
                "IS_ISOLATED": deg == 0,
                "MVA_CANONICAL": is_mva_canonical,
                "MVA_EXTENDED": is_mva_extended,
                "MVA_MECHANISM_LINK": meta.get("MVA_MECHANISM_LINK", "Uncharacterized / Secondary finding"),
                "PATHWAY": meta.get("PATHWAY", "General Cellular / Genomic Maintenance"),
                "TOTAL_PASSING_CANDIDATE_VARIANTS": meta.get("TOTAL_PASSING_CANDIDATE_VARIANTS", 0),
                "COMPOUND_HET_ARCHETYPE_FIT": meta.get("COMPOUND_HET_ARCHETYPE_FIT", False),
                "CLUSTERING_COEFFICIENT": meta.get("CLUSTERING_COEFFICIENT", 0.0),
                "TOP_ENRICHED_PATHWAYS": meta.get("TOP_ENRICHED_PATHWAYS", "None"),
            })

        network_nodes_df = pd.DataFrame(node_records)
        network_nodes_df = network_nodes_df.sort_values(
            by=["COMPONENT_RANK", "NETWORK_DEGREE"], ascending=[True, False]
        ).reset_index(drop=True)

        return network_nodes_df, comp_summary_df, components

    def create_network_figure(
        self,
        nodes: Sequence[str],
        edges_df: pd.DataFrame,
        nodes_df: pd.DataFrame,
        title: str,
        subtitle: str,
        seed: int = 42,
    ) -> go.Figure:
        """Create an interactive Plotly network figure with spring layout."""
        node_set = set(nodes)
        sub_edges = edges_df[
            edges_df["PREFERRED_NAME_A"].isin(node_set) & edges_df["PREFERRED_NAME_B"].isin(node_set)
        ]

        # Edge pairs
        edge_pairs = [
            (str(r["PREFERRED_NAME_A"]).strip().upper(), str(r["PREFERRED_NAME_B"]).strip().upper())
            for _, r in sub_edges.iterrows()
        ]

        # Compute deterministic coordinates
        sorted_nodes = sorted(list(node_set))
        pos = fruchterman_reingold_layout(sorted_nodes, edge_pairs, seed=seed, iterations=180)

        # Meta lookup
        meta_by_sym = {}
        for _, r in nodes_df.iterrows():
            sym = str(r["HGNC_SYMBOL"]).strip().upper()
            meta_by_sym[sym] = r

        # Calculate degrees within this subgraph
        deg_map: dict[str, int] = {n: 0 for n in sorted_nodes}
        for u, v in edge_pairs:
            deg_map[u] = deg_map.get(u, 0) + 1
            deg_map[v] = deg_map.get(v, 0) + 1

        # 1. Edge traces
        # Collect edge segments
        edge_x = []
        edge_y = []
        hover_edge_x = []
        hover_edge_y = []
        hover_edge_text = []

        for _, r in sub_edges.iterrows():
            u = str(r["PREFERRED_NAME_A"]).strip().upper()
            v = str(r["PREFERRED_NAME_B"]).strip().upper()
            if u in pos and v in pos:
                x0, y0 = pos[u]
                x1, y1 = pos[v]
                edge_x.extend([x0, x1, None])
                edge_y.extend([y0, y1, None])
                # Midpoint for hover
                xm, ym = (x0 + x1) / 2.0, (y0 + y1) / 2.0
                hover_edge_x.append(xm)
                hover_edge_y.append(ym)
                score = r.get("SCORE", 0.0)
                curated = r.get("CURATED_SUPPORT", "NO")
                escore = r.get("EXPERIMENTAL_SCORE", 0.0)
                dscore = r.get("DATABASE_SCORE", 0.0)
                txt = (
                    f"Interaction: {u} - {v}<br>"
                    f"STRING Score: {score:.3f}<br>"
                    f"Curated Support: {curated}<br>"
                    f"Exp Score: {escore:.3f} | DB Score: {dscore:.3f}"
                )
                hover_edge_text.append(txt)

        edge_trace = go.Scatter(
            x=edge_x,
            y=edge_y,
            line=dict(width=1.2, color="#bdc3c7"),
            hoverinfo="none",
            mode="lines",
            showlegend=False,
        )

        edge_hover_trace = go.Scatter(
            x=hover_edge_x,
            y=hover_edge_y,
            mode="markers",
            marker=dict(size=4, color="#bdc3c7", opacity=0.01),
            hoverinfo="text",
            hovertext=hover_edge_text,
            showlegend=False,
        )

        # 2. Node traces grouped by biological category
        # Groups: Canonical MVA, Extended MVA/SAC, Compound-Het Candidate, Other Candidate
        groups = {
            "Canonical MVA Gene (BUB1B, CEP57, TRIP13)": {
                "color": "#e74c3c",  # Crimson Red
                "size_base": 24,
                "symbol": "diamond",
                "nodes": [],
            },
            "Reported Aneuploidy / SAC Gene (MAD1L1, etc.)": {
                "color": "#e67e22",  # Amber Orange
                "size_base": 20,
                "symbol": "hexagon",
                "nodes": [],
            },
            "Compound-Het / Multi-Hit Candidate": {
                "color": "#2980b9",  # Steel Blue
                "size_base": 16,
                "symbol": "circle",
                "nodes": [],
            },
            "Candidate Interactor (Single-Hit / Maintenance)": {
                "color": "#16a085",  # Teal Green
                "size_base": 14,
                "symbol": "circle",
                "nodes": [],
            },
        }

        for n in sorted_nodes:
            meta = meta_by_sym.get(n, {})
            is_canonical = n in MVA_CANONICAL_GENES
            is_extended = n in MVA_EXTENDED_GENES
            is_archetype = bool(meta.get("COMPOUND_HET_ARCHETYPE_FIT", False))

            if is_canonical:
                groups["Canonical MVA Gene (BUB1B, CEP57, TRIP13)"]["nodes"].append(n)
            elif is_extended:
                groups["Reported Aneuploidy / SAC Gene (MAD1L1, etc.)"]["nodes"].append(n)
            elif is_archetype:
                groups["Compound-Het / Multi-Hit Candidate"]["nodes"].append(n)
            else:
                groups["Candidate Interactor (Single-Hit / Maintenance)"]["nodes"].append(n)

        node_traces = []
        for grp_name, grp_info in groups.items():
            grp_nodes = grp_info["nodes"]
            if not grp_nodes:
                continue

            nx = [pos[n][0] for n in grp_nodes]
            ny = [pos[n][1] for n in grp_nodes]
            n_text = []
            n_size = []

            for n in grp_nodes:
                meta = meta_by_sym.get(n, {})
                deg = deg_map.get(n, 0)
                n_size.append(grp_info["size_base"] + min(deg * 2.0, 16.0))
                pathway = meta.get("PATHWAY", "General Maintenance")
                mva_link = meta.get("MVA_MECHANISM_LINK", "None")
                n_vars = meta.get("TOTAL_PASSING_CANDIDATE_VARIANTS", 0)
                arch = meta.get("COMPOUND_HET_ARCHETYPE_FIT", False)
                top_p = meta.get("TOP_ENRICHED_PATHWAYS", "None")

                hover_str = (
                    f"<b>Gene: {n}</b><br>"
                    f"Subgraph Degree: {deg}<br>"
                    f"Passing Variants: {n_vars}<br>"
                    f"Compound-Het Fit: {'YES' if arch else 'NO'}<br>"
                    f"Pathway: {pathway}<br>"
                    f"MVA Link: {mva_link}<br>"
                    f"Top Enriched: {top_p}"
                )
                n_text.append(hover_str)

            trace = go.Scatter(
                x=nx,
                y=ny,
                mode="markers+text",
                name=grp_name,
                marker=dict(
                    symbol=grp_info["symbol"],
                    size=n_size,
                    color=grp_info["color"],
                    line=dict(width=1.5, color="#2c3e50"),
                ),
                text=[f"<b>{n}</b>" for n in grp_nodes],
                textposition="top center",
                textfont=dict(size=10, family="sans-serif", color="#2c3e50"),
                hoverinfo="text",
                hovertext=n_text,
            )
            node_traces.append(trace)

        fig = go.Figure(data=[edge_trace, edge_hover_trace, *node_traces])
        fig.update_layout(
            title=dict(
                text=f"<b>{title}</b><br><span style='font-size:12px; color:#7f8c8d;'>{subtitle}</span>",
                x=0.03,
                y=0.96,
            ),
            showlegend=True,
            legend=dict(
                x=0.02,
                y=0.05,
                bgcolor="rgba(255, 255, 255, 0.85)",
                bordercolor="#dcdde1",
                borderwidth=1,
            ),
            hovermode="closest",
            margin=dict(b=40, l=40, r=40, t=80),
            xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
            yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
            plot_bgcolor="#f8f9fa",
            paper_bgcolor="#ffffff",
            width=1100,
            height=850,
            annotations=[
                dict(
                    text="<b>Disclaimer:</b> Research exploration only. High-confidence STRING PPI (score >= 0.700). Node topology does not infer disease causality or drug suitability.",
                    showarrow=False,
                    xref="paper",
                    yref="paper",
                    x=0.5,
                    y=-0.04,
                    font=dict(size=11, color="#7f8c8d"),
                )
            ],
        )

        return fig

    def create_pathway_enrichment_figure(
        self,
        enrichment_df: pd.DataFrame,
        top_n: int = 25,
    ) -> go.Figure:
        """Create an interactive Plotly bubble/bar chart of top enriched pathways."""
        if enrichment_df.empty:
            fig = go.Figure()
            fig.update_layout(title="No Pathway Enrichment Terms Available")
            return fig

        # Clean and sort terms by FDR
        df = enrichment_df.copy()
        df["FDR"] = pd.to_numeric(df["FDR"], errors="coerce").fillna(1.0)
        df["ENTITIES_FOUND"] = pd.to_numeric(df["ENTITIES_FOUND"], errors="coerce").fillna(1)
        df["NEG_LOG_FDR"] = -np.log10(np.clip(df["FDR"], 1e-15, 1.0))

        # Filter to top N terms
        df = df.sort_values(by=["FDR", "ENTITIES_FOUND"], ascending=[True, False]).head(top_n)
        # Reverse order for horizontal bar/bubble display
        df = df.iloc[::-1].reset_index(drop=True)

        # Truncate very long pathway names for labels
        df["DISPLAY_NAME"] = df["PATHWAY_NAME"].apply(
            lambda s: (str(s)[:55] + "...") if len(str(s)) > 55 else str(s)
        )

        # Database color palette
        color_map = {
            "Reactome": "#2980b9",
            "g:Profiler (GO:BP)": "#27ae60",
            "g:Profiler (KEGG)": "#e67e22",
            "g:Profiler (WP)": "#8e44ad",
            "g:Profiler (REAC)": "#16a085",
            "g:Profiler (HP)": "#d35400",
        }

        fig = go.Figure()
        for db in df["SOURCE_DATABASE"].unique():
            sub = df[df["SOURCE_DATABASE"] == db]
            hover_text = [
                f"<b>{r['PATHWAY_NAME']}</b><br>"
                f"Source: {r['SOURCE_DATABASE']}<br>"
                f"ID: {r['PATHWAY_ID']}<br>"
                f"Entities Found: {r['ENTITIES_FOUND']} / {r.get('ENTITIES_TOTAL', 'N/A')}<br>"
                f"FDR: {r['FDR']:.2e} (-log10: {r['NEG_LOG_FDR']:.2f})<br>"
                f"Hit Genes: {r.get('HIT_GENES', 'N/A')}"
                for _, r in sub.iterrows()
            ]

            fig.add_trace(
                go.Scatter(
                    x=sub["NEG_LOG_FDR"],
                    y=sub["DISPLAY_NAME"],
                    mode="markers",
                    name=str(db),
                    marker=dict(
                        size=np.clip(np.sqrt(sub["ENTITIES_FOUND"]) * 4.5 + 8, 8, 38),
                        color=color_map.get(str(db), "#34495e"),
                        line=dict(width=1.2, color="#2c3e50"),
                        opacity=0.85,
                    ),
                    hoverinfo="text",
                    hovertext=hover_text,
                )
            )

        fig.update_layout(
            title=dict(
                text="<b>Top Enriched Biological Pathways and Gene Ontologies</b><br><span style='font-size:12px; color:#7f8c8d;'>Ranked by False Discovery Rate (FDR); bubble size represents entity hit count</span>",
                x=0.03,
                y=0.97,
            ),
            xaxis=dict(
                title="<b>Significance [-log10(FDR)]</b>",
                gridcolor="#f1f2f6",
                zeroline=True,
                zerolinecolor="#dcdde1",
            ),
            yaxis=dict(
                title="",
                autorange=True,
                gridcolor="#f1f2f6",
            ),
            plot_bgcolor="#ffffff",
            paper_bgcolor="#ffffff",
            margin=dict(l=320, r=40, t=70, b=60),
            legend=dict(
                title=dict(text="<b>Database</b>"),
                x=0.78,
                y=0.05,
                bgcolor="rgba(255, 255, 255, 0.85)",
                bordercolor="#dcdde1",
                borderwidth=1,
            ),
            width=1100,
            height=850,
            annotations=[
                dict(
                    text="FDR threshold 0.05: -log10(0.05) = 1.30",
                    x=1.30,
                    xref="x",
                    y=len(df) - 1,
                    yref="y",
                    showarrow=True,
                    arrowhead=2,
                    ax=40,
                    ay=-20,
                    font=dict(size=10, color="#e74c3c"),
                )
            ],
            shapes=[
                dict(
                    type="line",
                    x0=1.301,
                    x1=1.301,
                    y0=-0.5,
                    y1=len(df) - 0.5,
                    line=dict(color="#e74c3c", width=1.5, dash="dash"),
                )
            ],
        )

        return fig

    def run_visualization_pipeline(
        self,
        nodes_context_path: Path,
        edges_path: Path,
        enrichment_path: Path,
        output_dir: Path,
        flat_prefix: bool = True,
    ) -> dict[str, Path]:
        """Execute end-to-end Phase 2 visualization generation and persist artifacts."""
        output_dir.mkdir(parents=True, exist_ok=True)
        paths: dict[str, Path] = {}

        # 1. Load inputs
        logger.info("Loading Phase 2 inputs for visualization: %s", nodes_context_path)
        nodes_df = pd.read_csv(nodes_context_path, low_memory=False)
        edges_df = pd.read_csv(edges_path, low_memory=False)
        enrichment_df = pd.read_csv(enrichment_path, low_memory=False)

        # 2. Build complete network node table and component ranking summary
        network_nodes_df, comp_summary_df, components = self.build_network_node_table(
            nodes_df, edges_df
        )

        # 3. Write tables: 09_phase2_network_nodes.csv & 10_phase2_component_summary.csv
        nodes_csv_name = "09_phase2_network_nodes.csv" if flat_prefix else "network_nodes.csv"
        comp_csv_name = "10_phase2_component_summary.csv" if flat_prefix else "component_summary.csv"

        p_nodes = output_dir / nodes_csv_name
        p_comp = output_dir / comp_csv_name
        network_nodes_df.to_csv(p_nodes, index=False)
        comp_summary_df.to_csv(p_comp, index=False)
        paths["network_nodes"] = p_nodes
        paths["component_summary"] = p_comp

        # 4. Generate Figure 11: Largest Component (Component 1)
        c1_nodes = sorted(list(components[0])) if components else []
        fig_largest = self.create_network_figure(
            nodes=c1_nodes,
            edges_df=edges_df,
            nodes_df=nodes_df,
            title=f"Phase 2: Largest Connected Interactome Component (Component 1, N={len(c1_nodes)})",
            subtitle="STRING v12.5 high-confidence PPI (score >= 0.700) among passing candidate genes; fixed seed Fruchterman-Reingold layout.",
            seed=self.seed,
        )
        largest_html_name = "11_phase2_largest_component.html" if flat_prefix else "largest_component.html"
        p_largest = output_dir / largest_html_name
        fig_largest.write_html(p_largest, include_plotlyjs="cdn")
        paths["largest_component_html"] = p_largest

        # 5. Generate Figure 12: MVA-Focused Network Component
        # Find component containing BUB1B or MAD1L1
        mva_comp_nodes = []
        for c in components:
            if "BUB1B" in c or "MAD1L1" in c or "TRIP13" in c:
                mva_comp_nodes = sorted(list(c))
                break
        if not mva_comp_nodes and components:
            # Fallback if no known MVA gene in components: use 2nd largest
            mva_comp_nodes = sorted(list(components[1])) if len(components) > 1 else c1_nodes

        fig_mva = self.create_network_figure(
            nodes=mva_comp_nodes,
            edges_df=edges_df,
            nodes_df=nodes_df,
            title=f"Phase 2: MVA-Focused Checkpoint Interactome Subgraph (N={len(mva_comp_nodes)})",
            subtitle="Connected component anchoring BUB1B and MAD1L1 interactors (STRING score >= 0.700); deterministic layout.",
            seed=self.seed,
        )
        mva_html_name = "12_phase2_mva_focus_network.html" if flat_prefix else "mva_focus_network.html"
        p_mva = output_dir / mva_html_name
        fig_mva.write_html(p_mva, include_plotlyjs="cdn")
        paths["mva_focus_html"] = p_mva

        # 6. Generate Figure 13: Pathway Enrichment Plot
        fig_enrichment = self.create_pathway_enrichment_figure(enrichment_df, top_n=25)
        enrichment_html_name = "13_phase2_pathway_enrichment.html" if flat_prefix else "pathway_enrichment.html"
        p_enrichment = output_dir / enrichment_html_name
        fig_enrichment.write_html(p_enrichment, include_plotlyjs="cdn")
        paths["pathway_enrichment_html"] = p_enrichment

        return paths
