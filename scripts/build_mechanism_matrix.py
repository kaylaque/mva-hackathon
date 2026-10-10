#!/usr/bin/env python3
"""CLI script for Phase 1 Mechanism & Directionality Matrix generation.

Consumes caller-supplied consolidated candidate tables from Phase 0 (or master evidence tables).
Does NOT rerun upstream QC, VEP, snpEff, or Exomiser.
Constructs variant- and gene-level mechanism/directionality matrices with clear gates and provenance.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

# Ensure repository src/ is accessible
sys.path.insert(0, str(Path(__file__).parents[1] / "src"))

import pandas as pd

from mva_hackathon.candidate_consolidation import (
    CandidateConsolidator,
    ConsolidationGates,
)
from mva_hackathon.mechanism_matrix import MechanismMatrixBuilder


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Phase 1: Generate Mechanism and Directionality Matrices for candidate genes and variants."
    )
    parser.add_argument(
        "--candidates",
        type=Path,
        required=True,
        help="Path to the Phase 0 consolidated candidate table (CSV/TSV)",
    )
    parser.add_argument(
        "--output-variant-matrix",
        type=Path,
        required=True,
        help="Destination path for variant-level mechanism & directionality matrix (CSV/TSV)",
    )
    parser.add_argument(
        "--output-gene-matrix",
        type=Path,
        required=True,
        help="Destination path for gene-level mechanism & directionality matrix (CSV/TSV)",
    )
    parser.add_argument(
        "--output-provenance",
        type=Path,
        default=None,
        help="Optional destination path for JSON provenance record",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    print("=== Phase 1 Mechanism & Directionality Matrix Engine ===")
    print(f"Reading candidates from: {args.candidates}")

    sep = "\t" if args.candidates.suffix.lower() in {".tsv", ".txt"} else ","
    df = pd.read_csv(args.candidates, sep=sep, low_memory=False)

    required_columns = {
        "VARIANT_KEY",
        "SYMBOL",
        "GATE_OVERALL_PASS",
        "GATE_RARITY",
        "GATE_IMPACT",
        "CONSOLIDATION_PRIORITY_SCORE",
    }
    missing = sorted(required_columns - set(df.columns))
    if missing:
        raise ValueError(
            "Phase 1 requires the output of consolidate_candidates.py. "
            f"Missing columns: {', '.join(missing)}"
        )

    matrix_builder = MechanismMatrixBuilder()
    variant_matrix, gene_matrix, provenance = matrix_builder.build_matrix(df)

    # Write variant matrix
    args.output_variant_matrix.parent.mkdir(parents=True, exist_ok=True)
    v_sep = "\t" if args.output_variant_matrix.suffix.lower() in {".tsv", ".txt"} else ","
    variant_matrix.to_csv(args.output_variant_matrix, sep=v_sep, index=False)
    print(f"Variant mechanism matrix written to: {args.output_variant_matrix} ({len(variant_matrix):,} rows)")

    # Write gene matrix
    args.output_gene_matrix.parent.mkdir(parents=True, exist_ok=True)
    g_sep = "\t" if args.output_gene_matrix.suffix.lower() in {".tsv", ".txt"} else ","
    gene_matrix.to_csv(args.output_gene_matrix, sep=g_sep, index=False)
    print(f"Gene mechanism matrix written to: {args.output_gene_matrix} ({len(gene_matrix):,} genes)")

    if args.output_provenance:
        args.output_provenance.parent.mkdir(parents=True, exist_ok=True)
        with args.output_provenance.open("w", encoding="utf-8") as f:
            json.dump(provenance, f, indent=2)
        print(f"Provenance recorded at: {args.output_provenance}")

    print("Non-diagnostic research output. Complete.")


if __name__ == "__main__":
    main()
