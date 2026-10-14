"""Provenance-bound, sensitivity-only g:Profiler preparation."""

from __future__ import annotations

from datetime import UTC, datetime
import hashlib
import json
import math
from pathlib import Path
import subprocess
from typing import Any

import pandas as pd
from mva_hackathon.network_analysis import HttpClient, NetworkAnalysisConfig
from mva_hackathon.pathway_workflow import PathwayWorkflow, PathwayWorkflowConfig

_GPROFILER_URL = "https://biit.cs.ut.ee/gprofiler/api/gost/profile/"
_P0_COLUMNS = ("SYMBOL", "ENSEMBL_ID", "STATUS", "IN_FOREGROUND")
_P1_COLUMNS = ("TERM_ID", "TERM_NAME", "SOURCE", "RAW_P", "BH_Q", "FOLD_ENRICHMENT")
_ENRICHMENT_COLUMNS = (
    "SOURCE", "TERM_ID", "TERM_NAME", "PATHWAY_DB", "INTERSECTION_SIZE", "TERM_SIZE",
    "QUERY_SIZE", "BACKGROUND_SIZE", "RAW_P", "BH_Q", "PROVIDER_ADJUSTED_P",
    "PROVIDER_CORRECTION", "FOLD_ENRICHMENT", "CONTRIBUTING_GENES", "BACKGROUND_LABEL",
    "LOCAL_STATS_REASON",
)


def _canonical_hash(value: Any) -> str:
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"), default=str).encode()).hexdigest()


def _strict_bool(value: Any) -> bool:
    if isinstance(value, bool):
        return value
    if isinstance(value, str) and value.strip().lower() in {"true", "false"}:
        return value.strip().lower() == "true"
    raise ValueError("invalid IN_FOREGROUND boolean")


def _require_nonignored(path: Path) -> None:
    repo = Path(__file__).parents[2].resolve()
    resolved = path.resolve()
    if not resolved.is_relative_to(repo):
        return
    result = subprocess.run(
        ["git", "check-ignore", "--quiet", str(resolved.relative_to(repo))],
        cwd=repo,
        check=False,
        capture_output=True,
        text=True,
    )
    if result.returncode == 0:
        raise ValueError(f"ignored artifact path is forbidden: {path}")
    if result.returncode != 1:
        raise RuntimeError("cannot verify artifact ignore status")


def _read_baseline(baseline_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame, str]:
    p0_path = baseline_dir / "p0_gene_qc.csv"
    p1_path = baseline_dir / "p1_enrichment_results.csv"
    if not p0_path.is_file() or not p1_path.is_file():
        missing = [path.name for path in (p0_path, p1_path) if not path.is_file()]
        raise FileNotFoundError(f"missing pathway artifact: {', '.join(missing)}")
    resolved_baseline = baseline_dir.resolve()
    for path in (p0_path, p1_path):
        if not path.resolve().is_relative_to(resolved_baseline):
            raise ValueError("baseline artifact escapes baseline directory")
        _require_nonignored(path)
    headers = pd.read_csv(p0_path, nrows=0).columns
    missing = set(_P0_COLUMNS) - set(headers)
    if missing:
        raise ValueError(f"missing required P0 columns: {', '.join(sorted(missing))}")
    p0 = pd.read_csv(p0_path, usecols=list(_P0_COLUMNS), dtype="string")
    p0["SYMBOL"] = p0["SYMBOL"].fillna("").str.strip().str.upper()
    p0["STATUS"] = p0["STATUS"].fillna("").str.strip()
    p0["ENSEMBL_ID"] = p0["ENSEMBL_ID"].str.strip()
    if p0["SYMBOL"].eq("").any():
        raise ValueError("blank P0 SYMBOL")
    p0["IN_FOREGROUND"] = p0["IN_FOREGROUND"].map(_strict_bool)
    duplicate_columns = ["ENSEMBL_ID", "STATUS", "IN_FOREGROUND"]
    duplicates = p0.loc[p0.duplicated("SYMBOL", keep=False)]
    if not duplicates.empty and not duplicates.groupby("SYMBOL").apply(
        lambda rows: len(rows.loc[:, duplicate_columns].drop_duplicates()) == 1, include_groups=False
    ).all():
        raise ValueError("conflicting duplicate P0 rows")
    p0 = p0.drop_duplicates("SYMBOL", keep="first")
    mapped = p0.loc[p0["STATUS"].eq("mapped") & p0["ENSEMBL_ID"].notna() & p0["ENSEMBL_ID"].ne("")].copy()
    if p0.loc[p0["STATUS"].eq("mapped") & ~p0.index.isin(mapped.index)].shape[0]:
        raise ValueError("mapped P0 row has no ENSG identifier")
    p1_headers = pd.read_csv(p1_path, nrows=0).columns
    missing_p1 = set(_P1_COLUMNS) - set(p1_headers)
    if missing_p1:
        raise ValueError(f"missing required P1 columns: {', '.join(sorted(missing_p1))}")
    p1 = pd.read_csv(p1_path, usecols=list(_P1_COLUMNS))
    projection = {
        "p0": sorted(
            [{column: str(row[column]) for column in _P0_COLUMNS} for _, row in p0.iterrows()],
            key=lambda row: json.dumps(row, sort_keys=True),
        ),
        "p1": sorted(
            [{column: str(row[column]) for column in _P1_COLUMNS} for _, row in p1.iterrows()],
            key=lambda row: json.dumps(row, sort_keys=True),
        ),
    }
    return p0, mapped, _canonical_hash(projection)


def _read_multihit(path: Path, baseline: pd.DataFrame) -> set[str]:
    headers = pd.read_csv(path, nrows=0).columns.tolist()
    if headers != ["SYMBOL"]:
        raise ValueError("multi-hit CSV must contain exactly SYMBOL")
    symbols = pd.read_csv(path, dtype="string")["SYMBOL"].fillna("").str.strip().str.upper()
    if symbols.eq("").any():
        raise ValueError("blank multi-hit SYMBOL")
    result = set(symbols)
    known = set(baseline["SYMBOL"])
    if not result.issubset(known):
        raise ValueError("multi-hit symbol is unknown")
    foreground = set(baseline.loc[baseline["IN_FOREGROUND"], "SYMBOL"])
    if not result.issubset(foreground):
        raise ValueError("multi-hit symbols must be baseline foreground")
    return result


def _write_json(path: Path, payload: dict[str, Any]) -> None:
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def _prepare_analysis(
    output_dir: Path,
    name: str,
    mapped: pd.DataFrame,
    projection_hash: str,
    membership: set[str],
    sources: list[str],
) -> tuple[pd.DataFrame, PathwayWorkflow, dict[str, Any]]:
    run_dir = output_dir / name
    run_dir.mkdir()
    selected = mapped.copy()
    selected["IN_FOREGROUND"] = selected["SYMBOL"].isin(membership) if name == "multihit" else selected["IN_FOREGROUND"]
    workflow = PathwayWorkflow(PathwayWorkflowConfig(q_threshold=0.05))
    payload = workflow.gprofiler_payload(selected, domain_scope="custom", sources=sources)
    request = {"endpoint": _GPROFILER_URL, "payload": payload, "request_cache_key": workflow.request_cache_key(_GPROFILER_URL, payload)}
    membership_hash = _canonical_hash(sorted(membership))
    _write_json(run_dir / "request.json", request)
    _write_json(run_dir / "manifest.json", {
        "status": "PREPARED_NOT_EXECUTED",
        "analysis": name,
        "baseline_projection_sha256": projection_hash,
        "membership_sha256": membership_hash,
        "request_sha256": _canonical_hash(request),
        "sources": sources,
        "retrieved_at": None,
        "provider_version": "unknown",
    })
    pd.DataFrame(columns=_ENRICHMENT_COLUMNS).to_csv(run_dir / "enrichment.csv", index=False)
    pd.DataFrame(columns=_ENRICHMENT_COLUMNS).to_csv(run_dir / "prioritized.csv", index=False)
    return selected, workflow, request


def _valid_provider_response(response: dict[str, Any]) -> bool:
    for term in response["result"]:
        if not isinstance(term, dict):
            return False
        values = [term.get(name) for name in ("effective_domain_size", "term_size", "query_size", "intersection_size")]
        if any(isinstance(value, bool) or not isinstance(value, int) or value < 0 for value in values):
            return False
        domain, term_size, query_size, intersection = values
        if not domain or term_size > domain or query_size > domain or intersection > min(term_size, query_size):
            return False
        try:
            adjusted_p = float(term.get("p_value"))
        except (TypeError, ValueError):
            return False
        if not math.isfinite(adjusted_p) or not 0 <= adjusted_p <= 1:
            return False
    return True


def _complete_offline(run_dir: Path, selected: pd.DataFrame, workflow: PathwayWorkflow, response: dict[str, Any]) -> str:
    if not isinstance(response, dict) or not isinstance(response.get("result"), list):
        return "UNAVAILABLE_INVALID_PROVIDER_RESPONSE"
    metadata = response.get("meta", {})
    if not isinstance(metadata, dict) or not _valid_provider_response(response):
        return "UNAVAILABLE_INVALID_PROVIDER_RESPONSE"
    expected = {"organism": "hsapiens", "domain_scope": "custom", "significance_threshold_method": "fdr"}
    if any(key in metadata and metadata[key] != value for key, value in expected.items()):
        return "UNAVAILABLE_PROVIDER_METADATA_MISMATCH"
    enrichment = workflow.run_p1(selected, response).reindex(columns=_ENRICHMENT_COLUMNS)
    prioritized = workflow.run_p2(enrichment)
    enrichment.to_csv(run_dir / "enrichment.csv", index=False)
    prioritized.to_csv(run_dir / "prioritized.csv", index=False)
    manifest = json.loads((run_dir / "manifest.json").read_text(encoding="utf-8"))
    manifest.update({
        "status": "COMPLETED" if set(expected).issubset(metadata) else "UNVERIFIED_PROVIDER_METADATA",
        "retrieved_at": datetime.now(UTC).isoformat(),
        "provider_version": metadata.get("version", "unknown"),
        "enrichment_sha256": hashlib.sha256((run_dir / "enrichment.csv").read_bytes()).hexdigest(),
        "prioritized_sha256": hashlib.sha256((run_dir / "prioritized.csv").read_bytes()).hexdigest(),
    })
    _write_json(run_dir / "manifest.json", manifest)
    return str(manifest["status"])


def verify_sensitivity_artifacts(output_dir: Path) -> bool:
    """Validate hashes for completed sensitivity outputs before later consumption."""
    for name in ("multihit", "hp_excluded"):
        run_dir = output_dir / name
        try:
            manifest = json.loads((run_dir / "manifest.json").read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            return False
        if manifest.get("status") not in {"COMPLETED", "UNVERIFIED_PROVIDER_METADATA"}:
            continue
        for filename, field in (("enrichment.csv", "enrichment_sha256"), ("prioritized.csv", "prioritized_sha256")):
            path = run_dir / filename
            if not path.is_file() or manifest.get(field) != hashlib.sha256(path.read_bytes()).hexdigest():
                return False
    return True


def run_sensitivities(
    baseline_dir: Path,
    multihit_genes: Path,
    output_dir: Path,
    *,
    execute: bool = False,
    offline_fixture: Path | None = None,
) -> int:
    """Prepare or execute exactly two sensitivity analyses."""
    if execute and offline_fixture is not None:
        raise ValueError("--execute and --offline-fixture are mutually exclusive")
    if output_dir.exists():
        raise ValueError("output directory must not already exist")
    repo = Path(__file__).parents[2].resolve()
    if output_dir.resolve().is_relative_to(repo):
        raise ValueError("output directory must be outside the repository")
    if output_dir.resolve().is_relative_to(baseline_dir.resolve()) or baseline_dir.resolve().is_relative_to(output_dir.resolve()):
        raise ValueError("output and baseline directories must be disjoint")
    _require_nonignored(output_dir)
    baseline, mapped, projection_hash = _read_baseline(baseline_dir)
    multi_hit = _read_multihit(multihit_genes, baseline)
    mapped_multi_hit = multi_hit & set(mapped["SYMBOL"])
    if not mapped_multi_hit:
        raise ValueError("no mapped multi-hit genes remain")
    output_dir.mkdir()
    prepared = {
        "multihit": _prepare_analysis(output_dir, "multihit", mapped, projection_hash, mapped_multi_hit, ["GO:BP", "REAC", "WP", "HP"]),
        "hp_excluded": _prepare_analysis(output_dir, "hp_excluded", mapped, projection_hash, set(), ["GO:BP", "REAC", "WP"]),
    }
    multihit_manifest_path = output_dir / "multihit" / "manifest.json"
    multihit_manifest = json.loads(multihit_manifest_path.read_text(encoding="utf-8"))
    multihit_manifest["unmapped_multihit_symbols"] = sorted(multi_hit - mapped_multi_hit)
    _write_json(multihit_manifest_path, multihit_manifest)
    statuses = {name: "PREPARED_NOT_EXECUTED" for name in prepared}
    if offline_fixture is not None:
        for name, (selected, workflow, _) in prepared.items():
            fixture = offline_fixture / f"{name}.json"
            try:
                response = json.loads(fixture.read_text(encoding="utf-8"))
            except (OSError, json.JSONDecodeError):
                statuses[name] = "UNAVAILABLE_PROVIDER_RESPONSE"
            else:
                statuses[name] = _complete_offline(output_dir / name, selected, workflow, response)
    elif execute:
        client = HttpClient(NetworkAnalysisConfig())
        for name, (selected, workflow, request) in prepared.items():
            response = client.post(request["endpoint"], json_data=request["payload"])
            if response is None or response.status_code != 200:
                statuses[name] = "UNAVAILABLE_PROVIDER_RESPONSE"
                continue
            try:
                provider_response = response.json()
            except (TypeError, ValueError):
                statuses[name] = "UNAVAILABLE_PROVIDER_RESPONSE"
            else:
                statuses[name] = _complete_offline(output_dir / name, selected, workflow, provider_response)
    for name, status in statuses.items():
        if status.startswith("UNAVAILABLE"):
            manifest_path = output_dir / name / "manifest.json"
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
            manifest["status"] = status
            manifest["retrieved_at"] = datetime.now(UTC).isoformat()
            _write_json(manifest_path, manifest)
    report_lines = [
        "# Pathway sensitivity report",
        "",
        "- Baseline: immutable allowlisted projection; this is a separate artifact family.",
        *[f"- {name}: {status}" for name, status in statuses.items()],
        "",
        "Provider-adjusted q-values are not re-corrected.",
        "HP terms are phenotype annotations; they do not establish a biological mechanism.",
        "No automatic target or drug generation occurred.",
    ]
    (output_dir / "sensitivity_report.md").write_text("\n".join(report_lines) + "\n", encoding="utf-8")
    return 2 if any(status.startswith("UNAVAILABLE") for status in statuses.values()) else 0
