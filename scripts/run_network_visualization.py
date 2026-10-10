#!/usr/bin/env python3
"""CLI script for Phase 2 PPI and Pathway Network Visualization.

Consumes Phase 2 processed artifacts:
- 07_phase2_network_gene_context.csv (or network_gene_context.csv)
- 04_phase2_string_edges.csv (or string_edges.csv)
- 06_phase2_pathway_enrichment.csv (or pathway_enrichment.csv)

Generates:
- 09_phase2_network_nodes.csv: Complete node-level metrics (degree, component rank, clustering, MVA link)
- 10_phase2_component_summary.csv: Connected component ranking, sizes, densities, and MVA gene membership
- 11_phase2_largest_component.html: Interactive network graph for the largest connected component (Component 1)
- 12_phase2_mva_focus_network.html: Interactive network graph for the MVA-focused candidate component (Component 2, containing BUB1B and MAD1L1)
- 13_phase2_pathway_enrichment.html: Interactive pathway enrichment bubble chart (Reactome and g:Profiler)

Supports deterministic fixed-seed layout, flat output naming, and pure-Python execution.
Does not print raw patient-derived candidate rows.
"""

from __future__ import annotations

import argparse
import logging
import os
from pathlib import Path
import sys

# Ensure repository src/ is accessible
sys.path.insert(0, str(Path(__file__).parents[1] / "src"))

import pandas as pd

from mva_hackathon.network_visualization import NetworkVisualizer

logging.basicConfig(
    level=logging.INFO,
    format="[%(asctime)s] %(levelname)s [%(name)s] %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger("run_network_visualization")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Phase 2: Generate Interactive Network Visualizations and Component Summaries."
    )
    parser.add_argument(
        "--gene-context",
        type=Path,
        required=True,
        help="Path to Phase 2 network gene context CSV (07_phase2_network_gene_context.csv)",
    )
    parser.add_argument(
        "--edges",
        type=Path,
        required=True,
        help="Path to Phase 2 STRING edges CSV (04_phase2_string_edges.csv)",
    )
    parser.add_argument(
        "--enrichment",
        type=Path,
        required=True,
        help="Path to Phase 2 pathway enrichment CSV (06_phase2_pathway_enrichment.csv)",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        required=True,
        help="Directory where visualization outputs will be written",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Random seed for deterministic spring layout (default: 42)",
    )
    parser.add_argument(
        "--flat-names",
        action="store_true",
        help="Write outputs using flat Track 2 ordered prefixes (09_phase2_network_nodes.csv, etc.)",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    print("=== Phase 2 PPI & Pathway Network Visualization Engine ===")
    print(f"Gene context input: {args.gene_context}")
    print(f"STRING edges input: {args.edges}")
    print(f"Pathway enrichment input: {args.enrichment}")
    print(f"Output directory:   {args.output_dir}")
    print(f"Layout seed:        {args.seed}")
    print(f"Flat naming mode:   {args.flat_names}")

    visualizer = NetworkVisualizer(seed=args.seed)
    written_paths = visualizer.run_visualization_pipeline(
        nodes_context_path=args.gene_context,
        edges_path=args.edges,
        enrichment_path=args.enrichment,
        output_dir=args.output_dir,
        flat_prefix=args.flat_names,
    )

    print("\nPhase 2 Visualization Artifacts Generated Successfully:")
    for key, path in written_paths.items():
        size_bytes = os.path.getsize(path)
        size_kb = size_bytes / 1024.0
        print(f"  - {key}: {path.name} ({size_kb:.1f} KB)")

    print("\n[RESEARCH DISCLAIMER] Network graphs reflect physical/functional association in STRING DB v12.5 and do not infer causality, clinical diagnosis, or drug suitability.")


if __name__ == "__main__":
    main()
