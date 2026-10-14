#!/usr/bin/env python3
"""Write aggregate-only scientific validation report for pathway artifacts."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).parents[1] / "src"))

from mva_hackathon.pathway_audit import audit_pathway_outputs


def main() -> int:
    parser = argparse.ArgumentParser(description="Aggregate-only pathway validation audit.")
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    summary, report_path = audit_pathway_outputs(args.output_dir)
    baseline = summary["baseline"]
    print(
        "audit completed "
        f"baseline_bh_significant={baseline['provider_bh_le_0_05']} "
        f"threshold_two_pathways={summary['threshold_two'].get('prioritized_pathways', 'UNAVAILABLE')} "
        f"report={report_path.name}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
