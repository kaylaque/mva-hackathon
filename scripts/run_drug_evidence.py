#!/usr/bin/env python3
"""CLI script for Phase 3 Multi-Source Drug Repurposing Evidence and Safety Review.

Consumes Phase 1 gene mechanism matrix (and optional Phase 2 network context table).
Queries runtime authoritative biomedical pharmacological resources:
- Open Targets Platform v4 GraphQL (clinical candidates, tractability, safety liabilities)
- ChEMBL REST API v2.0 (target mechanisms of action, clinical phases)
- DGIdb GraphQL API (drug-gene interactions and evidence sources)

Applies strict variant mechanism -> pathway -> target action -> safety directionality
evaluations, explicitly flagging mechanistic hazards and prohibiting clinical recommendations.
"""

from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path
import sys

# Ensure repository src/ is accessible
sys.path.insert(0, str(Path(__file__).parents[1] / "src"))

import pandas as pd

from mva_hackathon.drug_evidence import (
    DrugEvidenceConfig,
    DrugEvidencePipeline,
    write_drug_evidence_outputs,
)

logging.basicConfig(
    level=logging.INFO,
    format="[%(asctime)s] %(levelname)s [%(name)s] %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger("run_drug_evidence")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Phase 3: Multi-Source Drug-Repurposing Evidence, Directionality & Safety Review."
    )
    parser.add_argument(
        "--gene-mechanism",
        type=Path,
        required=True,
        help="Path to Phase 1 gene mechanism matrix CSV (e.g. 02_phase1_gene_mechanism_matrix.csv)",
    )
    parser.add_argument(
        "--query-panel",
        type=Path,
        default=None,
        help="Optional path to Phase 3 query gene panel CSV (18_phase3_query_gene_panel.csv)",
    )
    parser.add_argument(
        "--network-context",
        type=Path,
        default=None,
        help="Optional path to Phase 2 network gene context CSV (e.g. 07_phase2_network_gene_context.csv)",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        required=True,
        help="Target directory where Phase 3 flat tables will be written",
    )
    parser.add_argument(
        "--seed-mode",
        choices=["query_panel", "mva_core", "mva_network", "passing_candidates", "custom"],
        default="mva_network",
        help="Candidate seed selection strategy (default: mva_network; auto-switches to query_panel if --query-panel provided)",
    )
    parser.add_argument(
        "--custom-genes",
        type=str,
        default="",
        help="Comma-separated list of gene symbols when --seed-mode=custom",
    )
    parser.add_argument(
        "--min-passing",
        type=int,
        default=1,
        help="Minimum passing candidate variants per gene for passing_candidates mode (default: 1)",
    )
    parser.add_argument(
        "--flat-names",
        action="store_true",
        default=True,
        help="Use flat numbered filenames (14_phase3_.. through 17_phase3_..). Default: True.",
    )
    parser.add_argument(
        "--no-flat-names",
        dest="flat_names",
        action="store_false",
        help="Use un-numbered filenames (drug_evidence_matrix.csv, etc.)",
    )
    parser.add_argument(
        "--disable-open-targets",
        action="store_true",
        default=False,
        help="Skip Open Targets queries",
    )
    parser.add_argument(
        "--disable-chembl",
        action="store_true",
        default=False,
        help="Skip ChEMBL queries",
    )
    parser.add_argument(
        "--disable-dgidb",
        action="store_true",
        default=False,
        help="Skip DGIdb queries",
    )
    parser.add_argument(
        "--offline",
        action="store_true",
        default=False,
        help="Run in offline mode without making network requests",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()

    if not args.gene_mechanism.exists():
        logger.error("Gene mechanism file does not exist: %s", args.gene_mechanism)
        return 1

    custom_gene_list = [
        g.strip().upper() for g in args.custom_genes.split(",") if g.strip()
    ]

    seed_mode = args.seed_mode
    if args.query_panel and args.query_panel.exists() and seed_mode == "mva_network":
        seed_mode = "query_panel"

    config = DrugEvidenceConfig(
        seed_mode=seed_mode,
        custom_genes=custom_gene_list,
        query_panel_path=args.query_panel,
        min_passing_variants=args.min_passing,
        enable_open_targets=not args.disable_open_targets and not args.offline,
        enable_chembl=not args.disable_chembl and not args.offline,
        enable_dgidb=not args.disable_dgidb and not args.offline,
        offline_mode=args.offline,
    )

    logger.info("Loading Phase 1 gene mechanism matrix from: %s", args.gene_mechanism)
    df_gene = pd.read_csv(args.gene_mechanism)

    df_net = None
    if args.network_context and args.network_context.exists():
        logger.info("Loading Phase 2 network gene context from: %s", args.network_context)
        df_net = pd.read_csv(args.network_context)

    pipeline = DrugEvidencePipeline(config)

    # 1. Select seed genes
    seed_genes = pipeline.select_seed_genes(df_gene, df_net)
    logger.info("Selected %d seed genes for Phase 3 evaluation: %s", len(seed_genes), seed_genes)

    # 2. Collect evidence and review safety/directionality
    logger.info("Executing Phase 3 pharmacological queries and directionality analysis...")
    matrix_df, ranking_df, action_map_df, provenance = pipeline.collect_evidence_for_seeds(
        seed_genes=seed_genes,
        gene_matrix_df=df_gene,
    )

    logger.info("Retrieved %d normalized drug-target evidence records.", len(matrix_df))
    logger.info("Synthesized %d ranked drug hypotheses.", len(ranking_df))
    logger.info("Contradiction/hazard count: %d", provenance.get("contradiction_count", 0))

    # 3. Write outputs
    written_paths = write_drug_evidence_outputs(
        evidence_matrix_df=matrix_df,
        ranking_df=ranking_df,
        action_map_df=action_map_df,
        provenance=provenance,
        output_dir=args.output_dir,
        flat_names=args.flat_names,
    )

    for k, p in written_paths.items():
        logger.info("Written [%s]: %s (%.1f KB)", k, p, p.stat().st_size / 1024)

    return 0


if __name__ == "__main__":
    sys.exit(main())
