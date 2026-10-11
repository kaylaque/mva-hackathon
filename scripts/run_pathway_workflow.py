#!/usr/bin/env python3
"""Run the independent research-only pathway workflow from gene-level input."""

from __future__ import annotations

import logging

import argparse
import json
from pathlib import Path
import sys
from typing import Any

sys.path.insert(0, str(Path(__file__).parents[1] / "src"))

import pandas as pd

from mva_hackathon.drug_evidence import (
    ChEMBLClient,
    DGIdbClient,
    DrugEvidenceConfig,
    OpenTargetsClient,
    PoliteHttpClient,
)
from mva_hackathon.network_analysis import HttpClient, NetworkAnalysisConfig
from mva_hackathon.pathway_workflow import PathwayWorkflow, PathwayWorkflowConfig

_STAGES = ("p0", "p1", "p2", "p3", "p4", "p5", "p6", "p7")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Research-only pathway workflow; accepts only gene-level Phase 1 input.")
    parser.add_argument("--gene-mechanism", type=Path, required=True)
    parser.add_argument("--network-context", type=Path)
    parser.add_argument("--network-panel", type=Path)
    parser.add_argument("--network-drugs", type=Path)
    parser.add_argument("--min-passing", type=int, default=1)
    parser.add_argument("--q-threshold", type=float, default=0.05)
    parser.add_argument("--top-pathways", type=int, default=15)
    parser.add_argument("--stage", choices=("all", *_STAGES), default="all")
    parser.add_argument("--offline-fixture", type=Path, help="Synthetic fixture directory; live network is prohibited.")
    parser.add_argument("--no-chembl", action="store_true")
    parser.add_argument("--no-open-targets", action="store_true")
    parser.add_argument("--no-dgidb", action="store_true")
    parser.add_argument("--david", action="store_true", help="Reserved; credentials must be provided through DAVID_EMAIL.")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--force", action="store_true")
    return parser.parse_args()


def _read_fixture(path: Path, name: str, default: Any) -> Any:
    fixture = path / name
    if not fixture.is_file():
        return default
    return json.loads(fixture.read_text(encoding="utf-8"))


def _drug_records(targets: pd.DataFrame, args: argparse.Namespace) -> tuple[list[dict[str, Any]], list[str]]:
    """Reuse drug clients; failures remain source gaps rather than no-hit claims."""
    config = DrugEvidenceConfig(
        enable_open_targets=not args.no_open_targets,
        enable_chembl=not args.no_chembl,
        enable_dgidb=not args.no_dgidb,
    )
    client = PoliteHttpClient(config)
    failures: list[str] = []
    records: list[dict[str, Any]] = []
    symbols = sorted(targets["SYMBOL"].dropna().astype(str).unique())
    if config.enable_chembl:
        chembl = ChEMBLClient(client, config)
        for symbol in symbols:
            try:
                for row in chembl.fetch_target_mechanisms(symbol):
                    records.append(
                        {
                            "SOURCE": "ChEMBL",
                            "TARGET_SYMBOL": symbol,
                            "ACTION_TYPE": row.get("action_type"),
                            "DRUG_NAME": row.get("molecule_chembl_id"),
                            "MOLECULE_CHEMBL_ID": row.get("molecule_chembl_id"),
                            "TARGET_CHEMBL_ID": row.get("target_chembl_id"),
                            "MECHANISM_ID": row.get("molecule_chembl_id"),
                        }
                    )
            except RuntimeError:
                failures.append("ChEMBL_UNAVAILABLE")
    if config.enable_dgidb:
        try:
            for node in DGIdbClient(client, config).fetch_gene_interactions(symbols):
                records.append(
                    {
                        "SOURCE": "DGIdb",
                        "TARGET_SYMBOL": (node.get("gene") or {}).get("name"),
                        "ACTION_TYPE": ";".join(item.get("type", "") for item in node.get("interactionTypes", [])),
                        "DRUG_NAME": (node.get("drug") or {}).get("name"),
                        "DGIDB_INTERACTION_ID": (node.get("drug") or {}).get("conceptId"),
                    }
                )
        except RuntimeError:
            failures.append("DGIdb_UNAVAILABLE")
    if config.enable_open_targets:
        ot = OpenTargetsClient(client, config)
        for ensembl_id, symbol in targets[["ENSEMBL_ID", "SYMBOL"]].dropna().drop_duplicates().itertuples(index=False):
            try:
                target = ot.fetch_target_drugs_and_safety(str(ensembl_id))
                for row in target.get("drugAndClinicalCandidates", {}).get("rows", []):
                    drug = row.get("drug") or {}
                    for mechanism in drug.get("mechanismsOfAction", {}).get("rows", []):
                        for mechanism_target in mechanism.get("targets", []):
                            records.append(
                                {
                                    "SOURCE": "OpenTargets",
                                    "TARGET_SYMBOL": mechanism_target.get("approvedSymbol", symbol),
                                    "ACTION_TYPE": mechanism.get("actionType"),
                                    "DRUG_NAME": drug.get("name"),
                                    "OPEN_TARGETS_DRUG_ID": drug.get("id"),
                                }
                            )
            except RuntimeError:
                failures.append("OPEN_TARGETS_UNAVAILABLE")
    return records, sorted(set(failures))


def main() -> int:
    args = parse_args()
    for logger_name in ("mva_hackathon.network_analysis", "mva_hackathon.drug_evidence"):
        logging.getLogger(logger_name).setLevel(logging.CRITICAL)
    output_dir = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    workflow = PathwayWorkflow(
        PathwayWorkflowConfig(
            min_passing=args.min_passing,
            q_threshold=args.q_threshold,
            top_pathways=args.top_pathways,
            offline=args.offline_fixture is not None,
        )
    )
    stages: dict[str, dict[str, Any]] = {}
    normalized, stages["p0"] = workflow.run_p0(args.gene_mechanism)
    if stages["p0"]["status"] == "failed":
        workflow.atomic_write_csv(normalized, output_dir / "p0_gene_qc.csv")
        workflow.atomic_write_json(workflow.provenance(stages), output_dir / "pathway_provenance.json")
        return 1

    if args.offline_fixture:
        response = _read_fixture(args.offline_fixture, "gprofiler.json", None)
        if response is None:
            raise ValueError("offline fixture is missing gprofiler.json")
        mapping = workflow.mapping_from_gprofiler(normalized, response)
        memberships = _read_fixture(args.offline_fixture, "memberships.json", {})
    else:
        http_client = HttpClient(NetworkAnalysisConfig())
        http_client.session.headers.update({"Accept": "application/json"})
        mapping = workflow.fetch_gprofiler_mapping(normalized, http_client)
        merged = normalized.merge(mapping, on="SYMBOL", how="left", validate="one_to_one")
        response = workflow.fetch_gprofiler(merged, http_client)
        memberships = {}
    normalized = normalized.merge(mapping, on="SYMBOL", how="left", validate="one_to_one")
    workflow.atomic_write_csv(
        normalized.rename(columns={"TOTAL_PASSING_CANDIDATE_VARIANTS": "PASSING_VARIANTS", "TOTAL_VARIANTS_EVALUATED": "TOTAL_VARIANTS"})[
            ["SYMBOL", "HGNC_ID", "ENSEMBL_ID", "STATUS", "IN_FOREGROUND", "PASSING_VARIANTS", "TOTAL_VARIANTS"]
        ],
        output_dir / "p0_gene_qc.csv",
    )

    p1 = workflow.run_p1(normalized, response)
    workflow.atomic_write_csv(p1, output_dir / "p1_enrichment_results.csv")
    stages["p1"] = workflow.stage_summary(p1, "NO_ENRICHMENT_TERMS")
    p2 = workflow.run_p2(p1)
    workflow.atomic_write_csv(p2, output_dir / "p2_prioritized_pathways.csv")
    stages["p2"] = workflow.stage_summary(p2, "NO_SIGNIFICANT_ENRICHMENT")

    if not args.offline_fixture:
        for term_id in p2.loc[p2["SOURCE"].eq("REAC"), "TERM_ID"].head(args.top_pathways):
            try:
                memberships[str(term_id)] = workflow.reactome_members(str(term_id), http_client)
            except RuntimeError:
                stages.setdefault("p4", {"status": "partial", "source_failures": []})["source_failures"].append("REACTOME_CONTENT_UNAVAILABLE")
    if args.offline_fixture:
        p3_evidence = _read_fixture(args.offline_fixture, "mechanism_evidence.json", [])
        p3_status: dict[str, Any] = workflow.stage_summary(pd.DataFrame(p3_evidence), "NO_ELIGIBLE_PATHWAYS")
    elif p2.empty:
        p3_evidence = []
        p3_status = workflow.stage_summary(pd.DataFrame(), "NO_ELIGIBLE_PATHWAYS")
    else:
        try:
            p3_evidence = workflow.mechanism_evidence(p2, http_client)
            p3_status = workflow.stage_summary(pd.DataFrame(p3_evidence), "NO_MECHANISM_EVIDENCE")
        except RuntimeError:
            p3_evidence = []
            p3_status = {"status": "partial", "count": 0, "reason": "UNIPROT_OR_QUICKGO_UNAVAILABLE"}
    p3 = workflow.run_p3(normalized, p3_evidence)
    workflow.atomic_write_csv(p3, output_dir / "p3_mechanism_evidence.csv")
    stages["p3"] = {**p3_status, "count": len(p3)}
    p4 = workflow.run_p4(p2, memberships, normalized)
    workflow.atomic_write_csv(p4, output_dir / "p4_expanded_targets.csv")
    p4_summary = workflow.stage_summary(p4, "NO_ELIGIBLE_PATHWAYS")
    stages["p4"] = {**p4_summary, **stages.get("p4", {}), "count": len(p4)}
    if p4.empty:
        source_records = []
        source_failures = []
        p5_outcome = "NO_ELIGIBLE_TARGETS"
    elif args.offline_fixture:
        source_records = _read_fixture(args.offline_fixture, "drug_records.json", [])
        source_failures = []
        p5_outcome = "NO_DRUG_EVIDENCE"
    else:
        source_records, source_failures = _drug_records(p4, args)
        p5_outcome = "NO_DRUG_EVIDENCE"
    p5 = workflow.run_p5(p4, source_records)
    workflow.atomic_write_csv(p5, output_dir / "p5_drug_evidence.csv")
    stages["p5"] = (
        {"status": "partial", "count": len(p5), "reason": "SOURCE_UNAVAILABLE", "source_failures": source_failures}
        if source_failures
        else workflow.stage_summary(p5, p5_outcome)
    )

    panel_path = args.network_panel or ((args.offline_fixture / "network_panel.csv") if args.offline_fixture else None)
    panel = pd.read_csv(panel_path) if panel_path and panel_path.is_file() else None
    p6 = workflow.run_p6(p4, panel)
    workflow.atomic_write_csv(pd.DataFrame([p6]), output_dir / "p6_evaluation_summary.csv")
    stages["p6"] = p6
    workflow.run_p7(output_dir / "p7_pathway_report.md", stages, p2, p4, p5)
    stages["p7"] = {"status": "completed", "count": 1}
    workflow.atomic_write_json(workflow.provenance(stages), output_dir / "pathway_provenance.json")
    return 2 if any(report.get("status") in {"partial", "failed"} for report in stages.values()) else 0


if __name__ == "__main__":
    raise SystemExit(main())
