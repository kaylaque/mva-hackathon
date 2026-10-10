#!/usr/bin/env python3
"""CLI script for Phase 0 Candidate Intake, Gating, and Consolidation.

Consumes caller-supplied existing evidence tables and optional Exomiser tables.
Does NOT rerun upstream QC, VEP, snpEff, or Exomiser.
Applies transparent research gates and outputs consolidated tables with full provenance.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

# Ensure repository src/ is accessible
sys.path.insert(0, str(Path(__file__).parents[1] / "src"))

from mva_hackathon.candidate_consolidation import (
    CandidateConsolidator,
    ConsolidationGates,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Phase 0: Consolidate existing variant evidence & Exomiser tables with transparent gates."
    )
    parser.add_argument(
        "--evidence",
        type=Path,
        required=True,
        help="Path to caller-supplied variant evidence table (e.g. variant_evidence_master.csv or VEP high/moderate table)",
    )
    parser.add_argument(
        "--exomiser-variants",
        type=Path,
        default=None,
        help="Optional path to caller-supplied Exomiser variants table (TSV/CSV)",
    )
    parser.add_argument(
        "--exomiser-genes",
        type=Path,
        default=None,
        help="Optional path to caller-supplied Exomiser genes table (TSV/CSV)",
    )
    parser.add_argument(
        "--output-candidates",
        type=Path,
        required=True,
        help="Destination path for consolidated candidate table (CSV/TSV)",
    )
    parser.add_argument(
        "--output-provenance",
        type=Path,
        default=None,
        help="Optional destination path for JSON provenance record",
    )
    parser.add_argument(
        "--max-af",
        type=float,
        default=0.001,
        help="Maximum population allele frequency threshold for rarity gate (default: 0.001)",
    )
    parser.add_argument(
        "--min-cadd",
        type=float,
        default=20.0,
        help="Minimum CADD PHRED score for deleteriousness gate (default: 20.0)",
    )
    parser.add_argument(
        "--min-spliceai",
        type=float,
        default=0.2,
        help="Minimum SpliceAI delta score for splicing gate (default: 0.2)",
    )
    parser.add_argument(
        "--min-revel",
        type=float,
        default=0.5,
        help="Minimum REVEL score for missense gate (default: 0.5)",
    )
    parser.add_argument(
        "--min-alphamissense",
        type=float,
        default=0.5,
        help="Minimum AlphaMissense score for missense gate (default: 0.5)",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    gates = ConsolidationGates(
        max_af=args.max_af,
        min_cadd_phred=args.min_cadd,
        min_spliceai_ds=args.min_spliceai,
        min_revel=args.min_revel,
        min_alphamissense=args.min_alphamissense,
    )

    consolidator = CandidateConsolidator(gates=gates)
    print(f"=== Phase 0 Candidate Consolidation ===")
    print(f"Reading evidence from: {args.evidence}")
    if args.exomiser_variants:
        print(f"Integrating Exomiser variants from: {args.exomiser_variants}")
    if args.exomiser_genes:
        print(f"Integrating Exomiser genes from: {args.exomiser_genes}")

    consolidated_df, provenance = consolidator.consolidate(
        evidence_path=args.evidence,
        exomiser_variants_path=args.exomiser_variants,
        exomiser_genes_path=args.exomiser_genes,
        output_path=args.output_candidates,
    )

    if args.output_provenance:
        args.output_provenance.parent.mkdir(parents=True, exist_ok=True)
        with args.output_provenance.open("w", encoding="utf-8") as f:
            json.dump(provenance, f, indent=2)
        print(f"Provenance recorded at: {args.output_provenance}")

    print(f"Consolidated table written to: {args.output_candidates}")
    print(f"Summary: {len(consolidated_df):,} total variants ingested, "
          f"{provenance['total_variants_passing_gates']:,} passing gates, "
          f"{provenance['total_multihit_genes']} multi-hit gene candidates.")
    print("Non-diagnostic research output. Complete.")


if __name__ == "__main__":
    main()
