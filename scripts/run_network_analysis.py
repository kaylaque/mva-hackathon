#!/usr/bin/env python3
"""CLI script for Phase 2 PPI and Pathway Network Context Analysis.

Consumes caller-supplied Phase 1 gene mechanism matrix (and optional variant matrix).
Queries external APIs at runtime:
- g:Profiler (g:Convert): Official HGNC mapping
- STRING Database v12.5: High-confidence PPI network (default score 700)
- Reactome Analysis Service: Pathway overrepresentation
- g:Profiler (g:GOSt): Functional enrichment

Writes downstream tables and JSON provenance to a caller-supplied directory.
Does NOT rerun upstream QC, VEP, snpEff, enrichment, or Exomiser.
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

from mva_hackathon.network_analysis import (
    NetworkAnalysisConfig,
    NetworkContextEngine,
    write_network_outputs,
)

logging.basicConfig(
    level=logging.INFO,
    format="[%(asctime)s] %(levelname)s [%(name)s] %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger("run_network_analysis")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Phase 2: Generate PPI and Pathway Network Context for Phase 1 Candidate Genes."
    )
    parser.add_argument(
        "--gene-mechanism",
        type=Path,
        required=True,
        help="Path to Phase 1 gene mechanism matrix (CSV/TSV)",
    )
    parser.add_argument(
        "--variant-matrix",
        type=Path,
        default=None,
        help="Optional path to Phase 1 variant mechanism matrix (CSV/TSV)",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        required=True,
        help="Caller-supplied directory where Phase 2 output tables will be written",
    )
    parser.add_argument(
        "--string-score",
        type=int,
        default=700,
        help="STRING interaction confidence threshold 0-1000 (default: 700 high confidence)",
    )
    parser.add_argument(
        "--string-additional-nodes",
        type=int,
        default=0,
        help="Number of network expansion interactors to add (default: 0, seeds only)",
    )
    parser.add_argument(
        "--min-passing",
        type=int,
        default=1,
        help="Keep only genes with >= N passing Phase 0 candidate variants (0 = use full gene matrix)",
    )
    parser.add_argument(
        "--string-batch-size",
        type=int,
        default=300,
        help="Max symbols per STRING GET request (batches large seed lists to avoid URL limits)",
    )
    parser.add_argument(
        "--min-escore",
        type=float,
        default=None,
        help="Optional minimum STRING experimental score channel threshold (0.0 - 1.0)",
    )
    parser.add_argument(
        "--min-dscore",
        type=float,
        default=None,
        help="Optional minimum STRING curated database score channel threshold (0.0 - 1.0)",
    )
    parser.add_argument(
        "--min-tscore",
        type=float,
        default=None,
        help="Optional minimum STRING textmining score channel threshold (0.0 - 1.0)",
    )
    parser.add_argument(
        "--min-ascore",
        type=float,
        default=None,
        help="Optional minimum STRING coexpression score channel threshold (0.0 - 1.0)",
    )
    parser.add_argument(
        "--timeout",
        type=float,
        default=20.0,
        help="HTTP request timeout in seconds (default: 20.0)",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    print("=== Phase 2 PPI & Pathway Network Context Engine ===")
    print(f"Reading gene mechanism matrix from: {args.gene_mechanism}")

    sep = "\t" if args.gene_mechanism.suffix.lower() in {".tsv", ".txt"} else ","
    gene_df = pd.read_csv(args.gene_mechanism, sep=sep, low_memory=False)

    var_df = None
    if args.variant_matrix and args.variant_matrix.exists():
        v_sep = "\t" if args.variant_matrix.suffix.lower() in {".tsv", ".txt"} else ","
        var_df = pd.read_csv(args.variant_matrix, sep=v_sep, low_memory=False)
        print(f"Loaded optional variant matrix: {len(var_df):,} rows")

    config = NetworkAnalysisConfig(
        string_required_score=args.string_score,
        string_additional_nodes=args.string_additional_nodes,
        min_passing_variants=args.min_passing,
        string_batch_size=args.string_batch_size,
        min_escore=args.min_escore,
        min_dscore=args.min_dscore,
        min_tscore=args.min_tscore,
        min_ascore=args.min_ascore,
        request_timeout=args.timeout,
    )

    engine = NetworkContextEngine(config=config)
    results = engine.run_analysis(gene_df, variant_matrix_df=var_df)

    # Write outputs
    written_paths = write_network_outputs(results, args.output_dir)

    print("\nPhase 2 Outputs Generated Successfully:")
    for key, path in written_paths.items():
        print(f"  - {key}: {path}")

    summary = results["network_summary_df"].iloc[0].to_dict()
    print("\nNetwork Summary Highlights:")
    print(f"  * Total seed genes evaluated: {summary.get('TOTAL_SEED_GENES')}")
    print(f"  * Total PPI edges (score >= {summary.get('STRING_SCORE_THRESHOLD')}): {summary.get('TOTAL_PPI_EDGES')}")
    print(f"  * Average degree: {summary.get('AVERAGE_DEGREE')}")
    print(f"  * Connected components: {summary.get('CONNECTED_COMPONENTS_COUNT')} (largest: {summary.get('LARGEST_COMPONENT_SIZE')})")
    print(f"  * Isolated genes: {summary.get('ISOLATED_GENES_COUNT')}")
    print(f"  * Enriched pathway terms: {len(results['pathway_enrichment_df']):,}")

    print("\n[RESEARCH DISCLAIMER] Network topology reflects physical/functional association and does not infer directionality or drug suitability.")


if __name__ == "__main__":
    main()
