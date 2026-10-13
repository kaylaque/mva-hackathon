#!/usr/bin/env python3
"""Run Compass: saved enrichment evidence to auditable research decisions.

Exit codes: 0 completed processing (scientific nulls and unresolved biology included),
2 source-incomplete, 1 unsafe or invalid inputs.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).parents[1] / "src"))

from mva_hackathon.compass import run_compass


class _Parser(argparse.ArgumentParser):
    def error(self, message: str) -> None:
        raise ValueError(message)


def parse_args() -> argparse.Namespace:
    parser = _Parser(description="Consolidate saved evidence and public records into research decisions.")
    parser.add_argument("--postanalysis-dir", type=Path, required=True, help="verified post-analysis output directory")
    parser.add_argument("--network-dir", type=Path, default=None, help="saved network artifacts (projected columns only)")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--cache-dir", type=Path, required=True, help="external public-response cache")
    parser.add_argument("--reviewed-evidence", type=Path, default=None, help="schema-validated reviewed evidence CSV")
    parser.add_argument("--resume", action="store_true", help="continue an existing output directory")
    parser.add_argument("--offline", action="store_true", help="use only verified cached public responses")
    return parser.parse_args()


def main() -> int:
    try:
        args = parse_args()
        result = run_compass(
            args.postanalysis_dir,
            args.network_dir,
            args.output_dir,
            cache_dir=args.cache_dir,
            reviewed_evidence=args.reviewed_evidence,
            resume=args.resume,
            offline=args.offline,
        )
    except (OSError, ValueError) as error:
        print(f"compass failed: {error}", file=sys.stderr)
        return 1
    print(json.dumps({"status": result["status"], "output_dir": result["output_dir"]}, sort_keys=True))
    return 2 if result["status"] == "source_incomplete" else 0


if __name__ == "__main__":
    raise SystemExit(main())
