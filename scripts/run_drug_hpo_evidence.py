#!/usr/bin/env python3
"""Run cached post-Compass drug-HPO evidence enrichment."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).parents[1] / "src"))

from mva_hackathon.drug_hpo_evidence import run_drug_hpo_evidence


class _Parser(argparse.ArgumentParser):
    def error(self, message: str) -> None:
        raise ValueError(message)


def parse_args() -> argparse.Namespace:
    parser = _Parser(description="Enrich a saved Compass run with cached public drug-HPO evidence.")
    parser.add_argument("--compass-dir", type=Path, required=True)
    parser.add_argument("--hpo-codes", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--cache-dir", type=Path, required=True)
    parser.add_argument("--reviewed-claims", type=Path)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--offline", action="store_true")
    parser.add_argument("--max-workers", type=int, default=8)
    return parser.parse_args()


def main() -> int:
    try:
        args = parse_args()
        result = run_drug_hpo_evidence(
            args.compass_dir, args.hpo_codes, args.output_dir,
            cache_dir=args.cache_dir, reviewed_claims=args.reviewed_claims,
            resume=args.resume, offline=args.offline, max_workers=args.max_workers,
        )
    except (OSError, ValueError) as error:
        print(f"drug-HPO evidence failed: {error}", file=sys.stderr)
        return 1
    print(json.dumps(result, sort_keys=True))
    return 2 if result["status"] == "source_incomplete" else 0


if __name__ == "__main__":
    raise SystemExit(main())
