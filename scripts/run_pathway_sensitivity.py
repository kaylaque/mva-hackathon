#!/usr/bin/env python3
"""Run only the controlled pathway sensitivity analyses."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).parents[1] / "src"))

from mva_hackathon.pathway_sensitivity import run_sensitivities


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Prepare or run controlled pathway sensitivity analyses.")
    parser.add_argument("--baseline-dir", type=Path, required=True)
    parser.add_argument("--multihit-genes", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--execute", action="store_true", help="Submit the two prepared g:Profiler requests.")
    parser.add_argument("--offline-fixture", type=Path, help="Synthetic multihit.json and hp_excluded.json responses.")
    args = parser.parse_args()
    if args.execute and args.offline_fixture is not None:
        parser.error("--execute and --offline-fixture are mutually exclusive")
    return args


def main() -> int:
    args = parse_args()
    try:
        status = run_sensitivities(
            args.baseline_dir,
            args.multihit_genes,
            args.output_dir,
            execute=args.execute,
            offline_fixture=args.offline_fixture,
        )
    except (FileNotFoundError, ValueError, NotImplementedError) as error:
        print(f"INVALID_INPUT: {error}", file=sys.stderr)
        return 1
    report = args.output_dir / "sensitivity_report.md"
    if report.is_file():
        for line in report.read_text(encoding="utf-8").splitlines():
            if line.startswith("- "):
                print(line[2:])
    return status


if __name__ == "__main__":
    raise SystemExit(main())
