#!/usr/bin/env python3
"""CLI script for Phase 3 Query Gene Panel Generation.

Derives a graph-grounded, multi-tiered gene query panel from Phase 1 & 2 artifacts:
- Gate 0: Symbol hygiene & non-query quarantine
- Gate 1: Core MVA causal genes (7 genes)
- Gate 2: Graph-derived interactome neighbors (MVA-focus component)
- Gate 3: Context modifiers (ranked Phase 1/2 candidates + curated literature references)
- Gate 4: Positive and programmatically derived negative drug-target controls
- Gate 5 & 6: Deduplication, integrity assertions, and panel sizing caps

Outputs:
- 18_phase3_query_gene_panel.csv
- 18_phase3_query_panel_provenance.json
- 19_phase3_query_panel_citations.csv

Never logs or prints patient data rows to standard output.
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path
import sys

# Ensure repository src/ is accessible
sys.path.insert(0, str(Path(__file__).parents[1] / "src"))

from mva_hackathon.query_gene_panel import (
    QueryGenePanelBuilder,
    QueryGenePanelConfig,
    write_query_panel_outputs,
)

logging.basicConfig(
    level=logging.INFO,
    format="[%(asctime)s] %(levelname)s [%(name)s] %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger("run_query_gene_panel")


def parse_args(args: list[str] | None = None) -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Generate Phase 3 Query Gene Panel for drug repurposing inquiry."
    )
    parser.add_argument(
        "--track2-dir",
        type=Path,
        default=Path("data/processed/track2"),
        help="Directory containing Track 2 flat artifacts (02_ through 10_). Default: data/processed/track2",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Directory to write panel artifacts. Default: matches --track2-dir",
    )
    parser.add_argument(
        "--flat-names",
        dest="flat_names",
        action="store_true",
        default=True,
        help="Use flat numbered filenames (18_phase3_.. and 19_phase3_..). Default: True.",
    )
    parser.add_argument(
        "--no-flat-names",
        dest="flat_names",
        action="store_false",
        help="Use hierarchical phase subdirectories (phase3/..).",
    )
    parser.add_argument(
        "--verify-literature",
        action="store_true",
        default=False,
        help="Verify primary literature via live NCBI E-utilities queries.",
    )
    parser.add_argument(
        "--max-context-modifiers",
        type=int,
        default=25,
        help="Maximum candidate-derived context modifiers to include. Default: 25.",
    )
    parser.add_argument(
        "--max-panel-size",
        type=int,
        default=120,
        help="Maximum total genes in the query panel. Default: 120.",
    )
    parser.add_argument(
        "--no-mva-focus-only",
        action="store_true",
        default=False,
        help="Do not restrict graph neighbors to the MVA-focus component.",
    )
    return parser.parse_args(args)


def main(args: list[str] | None = None) -> int:
    """Execute query gene panel generation."""
    parsed = parse_args(args)

    track2_dir = parsed.track2_dir.resolve()
    output_dir = (parsed.output_dir or parsed.track2_dir).resolve()

    if not track2_dir.is_dir():
        logger.error("Track 2 directory not found at %s", track2_dir)
        return 1

    config = QueryGenePanelConfig(
        track2_dir=track2_dir,
        output_dir=output_dir,
        flat_names=parsed.flat_names,
        max_context_modifiers=parsed.max_context_modifiers,
        max_panel_size=parsed.max_panel_size,
        mva_focus_component_only=not parsed.no_mva_focus_only,
        verify_literature=parsed.verify_literature,
    )

    logger.info("Initializing QueryGenePanelBuilder...")
    builder = QueryGenePanelBuilder(config)

    logger.info("Building query gene panel across selection gates...")
    panel_df, citations_df, provenance = builder.build_panel()

    cat_counts = provenance.get("panel_category_counts", {})
    ctrl_counts = provenance.get("control_counts", {})
    comp_info = provenance.get("detected_mva_focus_component", {})

    logger.info(
        "Panel derived: %d genes total. Categories: %s. Controls: %s. Detected MVA component rank: %s (size %d, %d edges).",
        len(panel_df),
        cat_counts,
        ctrl_counts,
        comp_info.get("component_rank"),
        comp_info.get("size", 0),
        comp_info.get("edge_count", 0),
    )

    logger.info("Writing output artifacts to %s...", output_dir)
    written = write_query_panel_outputs(
        panel_df=panel_df,
        provenance=provenance,
        citations_df=citations_df,
        output_dir=output_dir,
        flat_names=parsed.flat_names,
    )

    for name, path in written.items():
        size_bytes = path.stat().st_size
        logger.info("Written %s -> %s (%d bytes)", name, path, size_bytes)

    print(f"Phase 3 Query Gene Panel complete: {len(panel_df)} genes, {len(citations_df)} citations across {len(cat_counts)} categories.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
