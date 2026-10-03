#!/usr/bin/env python3
"""Consolidate completed VEP enrichment stages without running annotations."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).parents[1] / "src"))


from mva_hackathon.config import PipelinePaths
from mva_hackathon.enrichment import EvidencePipeline


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).parents[1])
    parser.add_argument("--vep", type=Path, help="Override data/interim/vep/vep_high_moderate.csv")
    args = parser.parse_args()
    evidence = EvidencePipeline(PipelinePaths.from_root(args.root)).consolidate(args.vep)
    print(f"Consolidated {len(evidence):,} unique variants.")


if __name__ == "__main__":
    main()
