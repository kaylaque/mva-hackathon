#!/usr/bin/env python3
"""CLI script for Phase 4 Signed-Pathway/Reaction Adjudication & Literature Validation.

Evaluates directional regulatory reachability and reaction flux alterations:
- Ingests SIGNOR and OmniPath curated signed networks with cross-database conflict resolution.
- Adjudicates Reactome biochemical reactions from enriched Phase 2 pathways.
- Classifies candidate pairs into H1 (disease mechanism correction), H2 (cancer complication),
  H3 (supportive care), or H0 (no class).
- Verifies literature evidence via Crossref REST (DOIs) and NCBI Entrez E-utilities (PMIDs).
- Executes multi-gate safety review, applying pediatric hazard thresholds and directionality checks.
- Produces prioritized, rankable reports ensuring best-supported hypotheses sort first.

DATA-BLIND RULE: Outputs to standard output contain aggregate counts, statistics, and paths only;
never patient-derived records.
"""

from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path
import sys

# Ensure repository src/ is accessible
sys.path.insert(0, str(Path(__file__).parents[1] / "src"))

from mva_hackathon.adjudication_safety import (
    Phase4AdjudicationConfig,
    Phase4AdjudicationEngine,
)

logging.basicConfig(
    level=logging.INFO,
    format="[%(asctime)s] %(levelname)s [%(name)s] %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger("run_phase4_adjudication")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Phase 4: Signed-Pathway/Reaction Adjudication & Literature Validation."
    )
    parser.add_argument(
        "--gene-mechanism",
        type=Path,
        default=Path("data/processed/track2/02_phase1_gene_mechanism_matrix.csv"),
        help="Path to Phase 1 gene mechanism matrix CSV",
    )
    parser.add_argument(
        "--drug-evidence",
        type=Path,
        default=Path("data/processed/track2/14_phase3_drug_evidence_matrix.csv"),
        help="Path to Phase 3 drug evidence matrix CSV",
    )
    parser.add_argument(
        "--safety-ranking",
        type=Path,
        default=Path("data/processed/track2/15_phase3_safety_ranking.csv"),
        help="Path to Phase 3 safety ranking CSV",
    )
    parser.add_argument(
        "--query-panel",
        type=Path,
        default=Path("data/processed/track2/18_phase3_query_gene_panel.csv"),
        help="Path to Phase 3 query gene panel CSV",
    )
    parser.add_argument(
        "--pathway-enrichment",
        type=Path,
        default=Path("data/processed/track2/06_phase2_pathway_enrichment.csv"),
        help="Path to Phase 2 pathway enrichment CSV",
    )
    parser.add_argument(
        "--network-context",
        type=Path,
        default=Path("data/processed/track2/07_phase2_network_gene_context.csv"),
        help="Path to Phase 2 network gene context CSV",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("data/processed/track2"),
        help="Directory where Phase 4 flat tables will be written",
    )
    parser.add_argument(
        "--flat-names",
        action="store_true",
        default=True,
        help="Use flat numbered filenames (20_ through 24_). Default: True.",
    )
    parser.add_argument(
        "--no-flat-names",
        dest="flat_names",
        action="store_false",
        help="Use un-numbered filenames (reaction_level_matrix.csv, etc.)",
    )
    parser.add_argument(
        "--max-path-length",
        type=int,
        default=2,
        help="Maximum path length for signed regulatory graph traversal (default: 2)",
    )
    parser.add_argument(
        "--reaction-limit",
        type=int,
        default=500,
        help="Hard cap on reaction records extracted from Reactome ContentService (default: 500)",
    )
    parser.add_argument(
        "--enable-omnipath-inferred",
        action="store_true",
        default=False,
        help="Include inferred OmniPath datasets (pathwayextra, kinaseextra)",
    )
    parser.add_argument(
        "--enable-synlethdb",
        action="store_true",
        default=False,
        help="Query synthetic lethality databases (SynLethDB / DepMap; default: False)",
    )
    parser.add_argument(
        "--offline-mode",
        action="store_true",
        default=False,
        help="Run without outbound network requests using local/mock data",
    )
    parser.add_argument(
        "--verify-literature",
        action="store_true",
        default=True,
        help="Perform live Crossref/NCBI literature verification (default: True)",
    )
    parser.add_argument(
        "--no-verify-literature",
        dest="verify_literature",
        action="store_false",
        help="Disable outbound Crossref/NCBI literature verification",
    )
    parser.add_argument(
        "--targets",
        type=str,
        default="",
        help="Comma-separated list of target symbols to override entry selection",
    )
    parser.add_argument(
        "--drugs",
        type=str,
        default="",
        help="Comma-separated list of drug names to override entry selection",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    targets_override = [t.strip().upper() for t in args.targets.split(",") if t.strip()]
    drugs_override = [d.strip() for d in args.drugs.split(",") if d.strip()]

    config = Phase4AdjudicationConfig(
        gene_mechanism_path=args.gene_mechanism,
        drug_evidence_path=args.drug_evidence,
        safety_ranking_path=args.safety_ranking,
        query_panel_path=args.query_panel,
        pathway_enrichment_path=args.pathway_enrichment,
        network_context_path=args.network_context,
        output_dir=args.output_dir,
        flat_names=args.flat_names,
        max_path_length=args.max_path_length,
        reaction_limit=args.reaction_limit,
        enable_omnipath_inferred=args.enable_omnipath_inferred,
        enable_synlethdb=args.enable_synlethdb,
        offline_mode=args.offline_mode,
        verify_literature=args.verify_literature,
        targets_override=targets_override,
        drugs_override=drugs_override,
    )

    logger.info("Initializing Phase 4 Adjudication Engine...")
    engine = Phase4AdjudicationEngine(config)
    results = engine.run()

    # Data-blind summary reporting: print only aggregate counts and output paths
    print("\n" + "=" * 65)
    print("PHASE 4 ADJUDICATION SUMMARY (AGGREGATE COUNTS & METRICS ONLY)")
    print("=" * 65)
    print(f"Entry Candidate Pairs Evaluated:   {results['entry_pairs_count']}")
    print(f"Signed Directed Edges Reconciled:  {results['edge_count']}")
    print(f"Cross-DB Curated Conflicts:        {results['conflict_count']}")
    print(f"Reactome Reactions Adjudicated:    {results['reaction_count']}")
    print(f"Adjudicated Safety Rows Emitted:   {results['safety_rows_count']}")
    print("\nGate Validation Statistics:")
    for gate, count in results["gate_stats"].items():
        print(f"  {gate:25s}: {count}")

    print("\nGenerated Outputs (all <= 8MB verified):")
    for name, path in results["output_paths"].items():
        size_bytes = path.stat().st_size if path.is_file() else 0
        size_kb = size_bytes / 1024
        print(f"  {name:25s}: {path} ({size_kb:.1f} KB)")

    print("=" * 65)
    print("RESEARCH ONLY: Strictly not clinical advice. All recommendations PROHIBITED.")
    print("=" * 65 + "\n")


if __name__ == "__main__":
    main()
