#!/usr/bin/env python3
"""Migration script for Track 2 result layout.

Copies hierarchical Phase 0, Phase 1, and Phase 2 result artifacts into flat,
ordered, descriptive files directly in data/processed/track2/ while keeping
the original phase directories intact (duplicates preserved).

Expected files and mapping:
- phase0/consolidated_candidates.csv      -> 00_phase0_consolidated_candidates.csv
- phase0/provenance.json                  -> 00_phase0_provenance.json
- phase1/variant_mechanism_matrix.csv     -> 01_phase1_variant_mechanism_matrix.csv
- phase1/gene_mechanism_matrix.csv        -> 02_phase1_gene_mechanism_matrix.csv
- phase1/provenance.json                  -> 02_phase1_provenance.json
- phase2/seed_genes.csv                   -> 03_phase2_seed_genes.csv
- phase2/string_edges.csv                 -> 04_phase2_string_edges.csv
- phase2/network_summary.csv              -> 05_phase2_network_summary.csv
- phase2/pathway_enrichment.csv           -> 06_phase2_pathway_enrichment.csv
- phase2/network_gene_context.csv         -> 07_phase2_network_gene_context.csv
- phase2/phase2_network_provenance.json   -> 08_phase2_network_provenance.json

Never reads or prints patient/genomic data rows into standard output or logs.
"""

from __future__ import annotations

import argparse
import hashlib
import logging
from pathlib import Path
import shutil
import sys
from typing import NamedTuple

logging.basicConfig(
    level=logging.INFO,
    format="[%(asctime)s] %(levelname)s: %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger("migrate_track2_layout")


class ArtifactMapping(NamedTuple):
    source_rel: str
    target_rel: str
    description: str
    optional: bool = False


TRACK2_ARTIFACT_MAPPINGS: list[ArtifactMapping] = [
    ArtifactMapping(
        source_rel="phase0/consolidated_candidates.csv",
        target_rel="00_phase0_consolidated_candidates.csv",
        description="Phase 0 consolidated candidate table with quality and research gates",
    ),
    ArtifactMapping(
        source_rel="phase0/provenance.json",
        target_rel="00_phase0_provenance.json",
        description="Phase 0 execution provenance and gating thresholds",
    ),
    ArtifactMapping(
        source_rel="phase1/variant_mechanism_matrix.csv",
        target_rel="01_phase1_variant_mechanism_matrix.csv",
        description="Phase 1 variant-level biological mechanism & directionality matrix",
    ),
    ArtifactMapping(
        source_rel="phase1/gene_mechanism_matrix.csv",
        target_rel="02_phase1_gene_mechanism_matrix.csv",
        description="Phase 1 gene-level mechanism matrix and biallelic archetype classification",
    ),
    ArtifactMapping(
        source_rel="phase1/provenance.json",
        target_rel="02_phase1_provenance.json",
        description="Phase 1 execution provenance and evaluated candidate counts",
    ),
    ArtifactMapping(
        source_rel="phase2/seed_genes.csv",
        target_rel="03_phase2_seed_genes.csv",
        description="Phase 2 seed genes extracted from passing Phase 1 candidates",
    ),
    ArtifactMapping(
        source_rel="phase2/string_edges.csv",
        target_rel="04_phase2_string_edges.csv",
        description="Phase 2 STRING DB v12.5 high-confidence PPI interactome edges",
    ),
    ArtifactMapping(
        source_rel="phase2/network_summary.csv",
        target_rel="05_phase2_network_summary.csv",
        description="Phase 2 global PPI network topology and graph statistics",
    ),
    ArtifactMapping(
        source_rel="phase2/pathway_enrichment.csv",
        target_rel="06_phase2_pathway_enrichment.csv",
        description="Phase 2 Reactome & g:Profiler functional pathway enrichment terms",
    ),
    ArtifactMapping(
        source_rel="phase2/network_gene_context.csv",
        target_rel="07_phase2_network_gene_context.csv",
        description="Phase 2 unified gene-level network context synthesizing graph metrics and pathways",
    ),
    ArtifactMapping(
        source_rel="phase2/phase2_network_provenance.json",
        target_rel="08_phase2_network_provenance.json",
        description="Phase 2 API queries, parameter configurations, and research caveats",
    ),
    ArtifactMapping(
        source_rel="phase2/network_nodes.csv",
        target_rel="09_phase2_network_nodes.csv",
        description="Phase 2 complete node-level interactome metrics, component ranks, and degrees",
        optional=True,
    ),
    ArtifactMapping(
        source_rel="phase2/component_summary.csv",
        target_rel="10_phase2_component_summary.csv",
        description="Phase 2 connected component ranking, sizes, densities, and MVA membership",
        optional=True,
    ),
    ArtifactMapping(
        source_rel="phase2/largest_component.html",
        target_rel="11_phase2_largest_component.html",
        description="Phase 2 interactive network graph for the largest connected component",
        optional=True,
    ),
    ArtifactMapping(
        source_rel="phase2/mva_focus_network.html",
        target_rel="12_phase2_mva_focus_network.html",
        description="Phase 2 interactive network graph for the MVA-focused candidate component",
        optional=True,
    ),
    ArtifactMapping(
        source_rel="phase2/pathway_enrichment.html",
        target_rel="13_phase2_pathway_enrichment.html",
        description="Phase 2 interactive pathway enrichment bubble chart",
        optional=True,
    ),
]


def check_sources_exist(base_dir: Path) -> tuple[bool, list[str]]:
    """Verify that all expected non-optional hierarchical source artifacts exist.

    Returns (all_present, missing_list).
    """
    missing: list[str] = []
    for mapping in TRACK2_ARTIFACT_MAPPINGS:
        if mapping.optional:
            continue
        src = base_dir / mapping.source_rel
        if not src.is_file():
            missing.append(f"{mapping.source_rel} (expected at {src})")
    return (len(missing) == 0, missing)


def sha256sum(path: Path) -> str:
    """Return a streaming SHA-256 digest without loading an artifact into memory."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def migrate_track2_artifacts(
    base_dir: Path,
    dry_run: bool = False,
    validate_only: bool = False,
) -> dict[str, dict[str, str | int | bool]]:
    """Perform safe migration or validation of Track 2 artifacts.

    Guarantees:
    - Never prints or reads patient/genomic data rows.
    - Preserves phase directories (phase0, phase1, phase2).
    - Checks that all 11 source artifacts exist before copying.
    - Validates file existence, byte sizes, and SHA-256 content hashes after copying.
    """
    if not base_dir.is_dir():
        raise FileNotFoundError(f"Track 2 base directory does not exist: {base_dir}")

    all_present, missing_sources = check_sources_exist(base_dir)
    if not all_present and not validate_only:
        missing_msg = "\n  - ".join(missing_sources)
        raise FileNotFoundError(
            f"Cannot migrate Track 2 artifacts; missing required source files in {base_dir}:\n  - {missing_msg}"
        )

    results: dict[str, dict[str, str | int | bool]] = {}

    for mapping in TRACK2_ARTIFACT_MAPPINGS:
        src_path = base_dir / mapping.source_rel
        dst_path = base_dir / mapping.target_rel

        if not src_path.is_file() and mapping.optional:
            continue

        if validate_only:
            src_exists = src_path.is_file()
            dst_exists = dst_path.is_file()
            src_bytes = src_path.stat().st_size if src_exists else -1
            dst_bytes = dst_path.stat().st_size if dst_exists else -1
            src_sha256 = sha256sum(src_path) if src_exists else ""
            dst_sha256 = sha256sum(dst_path) if dst_exists else ""
            size_match = (src_bytes == dst_bytes and src_bytes > 0) if (src_exists and dst_exists) else False
            content_match = bool(src_exists and dst_exists and src_sha256 == dst_sha256)

            results[mapping.target_rel] = {
                "source": str(mapping.source_rel),
                "target": str(mapping.target_rel),
                "source_exists": src_exists,
                "target_exists": dst_exists,
                "source_bytes": src_bytes,
                "target_bytes": dst_bytes,
                "source_sha256": src_sha256,
                "target_sha256": dst_sha256,
                "size_match": size_match,
                "content_match": content_match,
                "valid": bool(dst_exists and size_match and content_match),
            }
            continue

        src_bytes = src_path.stat().st_size

        if dry_run:
            results[mapping.target_rel] = {
                "action": "DRY_RUN_COPY",
                "source": str(mapping.source_rel),
                "target": str(mapping.target_rel),
                "source_bytes": src_bytes,
                "description": mapping.description,
            }
            continue

        # Execute safe copy preserving file metadata
        shutil.copy2(src_path, dst_path)

        # Immediate validation
        dst_exists = dst_path.is_file()
        dst_bytes = dst_path.stat().st_size if dst_exists else -1
        src_sha256 = sha256sum(src_path)
        dst_sha256 = sha256sum(dst_path) if dst_exists else ""
        size_match = src_bytes == dst_bytes
        content_match = dst_exists and src_sha256 == dst_sha256

        if not dst_exists or not size_match or not content_match:
            raise IOError(
                f"Validation failed after copying {src_path} to {dst_path}: "
                f"exists={dst_exists}, src_size={src_bytes}, dst_size={dst_bytes}, "
                f"content_match={content_match}"
            )

        results[mapping.target_rel] = {
            "action": "COPIED",
            "source": str(mapping.source_rel),
            "target": str(mapping.target_rel),
            "source_bytes": src_bytes,
            "target_bytes": dst_bytes,
            "source_sha256": src_sha256,
            "target_sha256": dst_sha256,
            "size_match": size_match,
            "content_match": content_match,
            "description": mapping.description,
        }

    return results


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Migrate Track 2 hierarchical phase artifacts into flat ordered descriptive filenames."
    )
    parser.add_argument(
        "--base-dir",
        type=Path,
        default=Path("data/processed/track2"),
        help="Base directory containing phase0, phase1, phase2 subdirectories (default: data/processed/track2)",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Inspect existing sources and planned target copies without copying files",
    )
    parser.add_argument(
        "--validate-only",
        action="store_true",
        help="Verify existing source and target artifacts without modifying files",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    base_dir = args.base_dir.resolve()

    mode_label = "VALIDATION" if args.validate_only else ("DRY-RUN" if args.dry_run else "MIGRATION")
    logger.info(f"Starting Track 2 result layout {mode_label} in: {base_dir}")

    results = migrate_track2_artifacts(
        base_dir=base_dir,
        dry_run=args.dry_run,
        validate_only=args.validate_only,
    )

    print(f"\n========================================================")
    print(f" Track 2 Artifact Migration Summary ({mode_label})")
    print(f" Target Directory: {base_dir}")
    print(f"========================================================")

    all_valid = True
    for target_name, info in results.items():
        if args.validate_only:
            status = "VALID" if info["valid"] else "FAILED"
            if not info["valid"]:
                all_valid = False
            print(f"[{status:6s}] {target_name:<42s} (src: {info['source_bytes']:,} B, tgt: {info['target_bytes']:,} B)")
        elif args.dry_run:
            print(f"[PLAN  ] {info['source']:<38s} -> {target_name:<42s} ({info['source_bytes']:,} B)")
        else:
            print(f"[COPIED] {info['source']:<38s} -> {target_name:<42s} ({info['target_bytes']:,} B)")

    print(f"========================================================")
    if args.validate_only:
        if all_valid:
            print("Validation successful: All 11 flat target artifacts match source files by SHA-256.")
            sys.exit(0)
        else:
            print("Validation failed: One or more artifacts are missing or differ from their source.", file=sys.stderr)
            sys.exit(1)
    elif args.dry_run:
        print(f"Dry-run complete: {len(results)} artifacts planned for copy. No files were modified.")
    else:
        print(f"Migration complete: {len(results)} artifacts copied and validated successfully.")
        print("Note: Phase directories (phase0, phase1, phase2) have been kept intact as requested.")


if __name__ == "__main__":
    main()
