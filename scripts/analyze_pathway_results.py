#!/usr/bin/env python3
"""Verify saved pathway artifacts and write offline post-analysis tables."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).parents[1] / "src"))

from mva_hackathon.pathway_postanalysis import analyze_pathway_results

class _Parser(argparse.ArgumentParser):
    def error(self, message: str) -> None:
        raise ValueError(message)


def parse_args() -> argparse.Namespace:
    parser = _Parser(description="Offline, provenance-bound DAVID and g:Profiler post-analysis.")
    parser.add_argument("--reuse-pathway-dir", type=Path, required=True)
    parser.add_argument("--david-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def main() -> int:
    try:
        args = parse_args()
        result = analyze_pathway_results(args.reuse_pathway_dir, args.david_dir, args.output_dir)
    except (OSError, ValueError) as error:
        print(f"pathway post-analysis failed: {error}", file=sys.stderr)
        return 1
    print(json.dumps(result, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
