#!/usr/bin/env python3
"""Prepare or explicitly run Exomiser; no work starts without --execute."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).parents[1] / "src"))


from mva_hackathon.config import PipelinePaths
from mva_hackathon.enrichment import EvidencePipeline
from mva_hackathon.exomiser import ExomiserRun, ExomiserRunner


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).parents[1])
    parser.add_argument("--jar", type=Path, help="Path to exomiser-cli.jar")
    parser.add_argument("--hpo", action="append", help="Repeat for each HPO term, e.g. HP:0000001")
    parser.add_argument("--execute", action="store_true", help="Run Exomiser after creating inputs")
    parser.add_argument("--import-results", type=Path, help="Normalize a completed Exomiser TSV/CSV export")
    args = parser.parse_args()
    if not args.import_results and (args.jar is None or not args.hpo):
        parser.error("--jar and at least one --hpo are required unless --import-results is used")
    paths = PipelinePaths.from_root(args.root)
    runner = ExomiserRunner(paths, ExomiserRun(args.jar, tuple(args.hpo or ())))
    if args.import_results:
        runner.normalize_results(args.import_results)
        print(f"Imported Exomiser results into {paths.stages['exomiser']}")
        return
    variants = EvidencePipeline(paths).load_vep()
    _, config = runner.write_input(variants)
    print(f"Prepared Exomiser config: {config}")
    if args.execute:
        runner.execute(config)


if __name__ == "__main__":
    main()
