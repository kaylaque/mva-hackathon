"""Resumable import of human-operated native DAVID exports.

This module only prepares gene-only DAVID lists and validates/imports exports saved
by the native DAVID interface.  It never performs browser or service requests.
"""

from __future__ import annotations

from datetime import UTC, datetime
import hashlib
import json
import math
from pathlib import Path
import re
import shutil
from typing import Any

import pandas as pd

from mva_hackathon.pathway_audit import _require_nonignored_paths
from mva_hackathon.pathway_workflow import PathwayWorkflow

_NAMESPACE = "ENSEMBL_GENE_ID"
_TAXON = 9606
_FOREGROUND_LIMIT = 3000
_P0_COLUMNS = ("SYMBOL", "ENSEMBL_ID", "STATUS", "IN_FOREGROUND")
_P1_COLUMNS = (
    "SOURCE", "TERM_ID", "TERM_NAME", "BH_Q", "PROVIDER_ADJUSTED_P",
    "PROVIDER_CORRECTION", "RAW_P", "FOLD_ENRICHMENT", "CONTRIBUTING_GENES",
    "BACKGROUND_LABEL",
)
_REQUIRED_P1_COLUMNS = {"SOURCE", "TERM_ID", "BH_Q", "CONTRIBUTING_GENES"}
_CHART_COLUMNS = {"Category", "Term", "Count", "P-Value", "Benjamini"}
_CLUSTER_COLUMNS = {"Cluster", "Cluster Enrichment Score", "Category", "Term", "Count", "P-Value", "Benjamini"}
_ID_PATTERN = re.compile(r"^(GO:\d+|R-HSA-\d+|WP\d+|HP:\d+|hsa\d+)(?=~|\s|$)", re.IGNORECASE)


def _sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def _canonical(value: Any) -> Any:
    if value is None or value is pd.NA or (isinstance(value, float) and math.isnan(value)):
        return None
    if isinstance(value, bool):
        return value
    if isinstance(value, (int, str)):
        return value
    if isinstance(value, float):
        return value if math.isfinite(value) else None
    if isinstance(value, dict):
        return {str(key): _canonical(item) for key, item in sorted(value.items())}
    if isinstance(value, (list, tuple, set)):
        return [_canonical(item) for item in value]
    return str(value)


def _json_hash(value: Any) -> str:
    return _sha256_bytes(json.dumps(_canonical(value), sort_keys=True, separators=(",", ":")).encode())


def _read_json(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ValueError(f"invalid JSON artifact: {path.name}") from error
    if not isinstance(value, dict):
        raise ValueError(f"invalid JSON object: {path.name}")
    return value


def _bool(value: Any, name: str) -> bool:
    if isinstance(value, bool):
        return value
    if isinstance(value, str) and value.strip().lower() in {"true", "false"}:
        return value.strip().lower() == "true"
    raise ValueError(f"invalid boolean {name}")


def _nonempty_string(value: Any, name: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"missing nonempty {name}")
    return value.strip()


def _nonnegative_int(value: Any, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise ValueError(f"invalid nonnegative {name}")
    return value


def _finite(value: Any, name: str, *, unit_interval: bool = False) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError) as error:
        raise ValueError(f"invalid {name}") from error
    if not math.isfinite(number) or (unit_interval and not 0 <= number <= 1):
        raise ValueError(f"invalid {name}")
    return number


def _repo_processed_root() -> Path:
    return Path(__file__).parents[2].resolve() / "data" / "processed"


def _validate_baseline_path(baseline_dir: Path, paths: list[Path]) -> None:
    resolved_root = baseline_dir.resolve()
    processed_root = _repo_processed_root()
    for path in paths:
        resolved = path.resolve()
        if not resolved.is_relative_to(resolved_root):
            raise ValueError("baseline artifact escapes baseline directory")
        if resolved.is_relative_to(processed_root):
            continue
        _require_nonignored_paths([path])


def _validate_output_path(baseline_dir: Path, output_dir: Path) -> None:
    baseline = baseline_dir.resolve()
    output = output_dir.resolve()
    if output == baseline or output.is_relative_to(baseline) or baseline.is_relative_to(output):
        raise ValueError("DAVID output and baseline directories must be disjoint")
    _require_nonignored_paths([output_dir, output_dir / "david_manifest.json"])


def _project_csvs(baseline_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    p0_path = baseline_dir / "p0_gene_qc.csv"
    p1_path = baseline_dir / "p1_enrichment_results.csv"
    if not p0_path.is_file() or not p1_path.is_file():
        missing = [path.name for path in (p0_path, p1_path) if not path.is_file()]
        raise FileNotFoundError(f"missing pathway artifact: {', '.join(missing)}")
    _validate_baseline_path(baseline_dir, [p0_path, p1_path])
    p0 = pd.read_csv(p0_path, usecols=list(_P0_COLUMNS))
    p1_headers = pd.read_csv(p1_path, nrows=0).columns
    missing = _REQUIRED_P1_COLUMNS - set(p1_headers)
    if missing:
        raise ValueError(f"missing required P1 columns: {', '.join(sorted(missing))}")
    p1 = pd.read_csv(p1_path, usecols=lambda column: column in _P1_COLUMNS)
    return p0, p1


def _projection_hash(frame: pd.DataFrame, columns: tuple[str, ...]) -> str:
    rows = [
        {column: _canonical(row[column]) for column in columns if column in frame}
        for _, row in frame.loc[:, [column for column in columns if column in frame]].iterrows()
    ]
    return _json_hash(sorted(rows, key=lambda row: json.dumps(row, sort_keys=True, separators=(",", ":"))))


def _resolved_lists(p0: pd.DataFrame) -> tuple[list[str], list[str], dict[str, str]]:
    if p0[list(_P0_COLUMNS)].isna().any(axis=None):
        # ENSEMBL_ID may be absent only for unmapped rows.
        invalid = p0.loc[p0["STATUS"].astype(str).eq("mapped") & p0["ENSEMBL_ID"].isna()]
        if not invalid.empty:
            raise ValueError("mapped P0 row has no ENSG identifier")
    frame = p0.copy()
    frame["SYMBOL"] = frame["SYMBOL"].fillna("").astype(str).str.strip().str.upper()
    frame["STATUS"] = frame["STATUS"].fillna("").astype(str).str.strip()
    frame["ENSEMBL_ID"] = frame["ENSEMBL_ID"].astype("string").str.strip()
    frame["IN_FOREGROUND"] = [_bool(value, "IN_FOREGROUND") for value in frame["IN_FOREGROUND"]]
    if frame["SYMBOL"].eq("").any():
        raise ValueError("blank P0 SYMBOL")
    duplicate_columns = ["SYMBOL", "ENSEMBL_ID", "STATUS", "IN_FOREGROUND"]
    conflicts = frame.loc[frame.duplicated("SYMBOL", keep=False)].groupby("SYMBOL").filter(
        lambda rows: len(rows.loc[:, duplicate_columns].drop_duplicates()) > 1
    )
    if not conflicts.empty:
        raise ValueError("conflicting duplicate P0 rows")
    frame = frame.drop_duplicates("SYMBOL", keep="first")
    mapped = frame.loc[frame["STATUS"].eq("mapped") & frame["ENSEMBL_ID"].notna()].copy()
    mapped["ENSEMBL_ID"] = mapped["ENSEMBL_ID"].astype(str)
    genes = PathwayWorkflow().gprofiler_payload(mapped)
    foreground, background = list(genes["query"]), list(genes["background"])
    if not set(foreground).issubset(background):
        raise ValueError("foreground must be a subset of background")
    if len(foreground) > _FOREGROUND_LIMIT:
        raise ValueError(f"DAVID foreground query exceeds {_FOREGROUND_LIMIT} identifiers")
    symbol_to_ensembl = (
        mapped.groupby("SYMBOL")["ENSEMBL_ID"].agg(lambda values: sorted(set(values))).to_dict()
    )
    unambiguous = {symbol: values[0] for symbol, values in symbol_to_ensembl.items() if len(values) == 1}
    return foreground, background, unambiguous


def _prepared_manifest(p0: pd.DataFrame, p1: pd.DataFrame, foreground: list[str], background: list[str]) -> dict[str, Any]:
    foreground_bytes = ("\n".join(foreground) + "\n").encode()
    background_bytes = ("\n".join(background) + "\n").encode()
    projections = {
        "p0_gene_qc.csv": _projection_hash(p0, _P0_COLUMNS),
        "p1_enrichment_results.csv": _projection_hash(p1, _P1_COLUMNS),
    }
    settings = {
        "chart": {"Count": 2, "EASE": 0.1},
        "clustering": {"displayed": "Medium", "actual_categories_and_thresholds": None},
    }
    list_sha256 = _sha256_bytes(foreground_bytes + b"\0" + background_bytes)
    run_key = PathwayWorkflow.request_cache_key(
        "david-native-ui",
        {"projections": projections, "list_sha256": list_sha256, "namespace": _NAMESPACE, "taxon": _TAXON, "settings": settings},
    )
    return {
        "schema_version": 1,
        "run_key": run_key,
        "baseline_projection_sha256": _json_hash(projections),
        "baseline_projections": {"columns": {"p0_gene_qc.csv": list(_P0_COLUMNS), "p1_enrichment_results.csv": list(_P1_COLUMNS)}, "sha256": projections},
        "list_sha256": list_sha256,
        "namespace": _NAMESPACE,
        "taxon": _TAXON,
        "foreground_count": len(foreground),
        "background_count": len(background),
        "status": "awaiting_web_results",
        "settings": settings,
        "result_sha256": None,
    }


def _write_json(path: Path, payload: dict[str, Any]) -> None:
    PathwayWorkflow.atomic_write_json(payload, path)


def _matching_manifest(existing: dict[str, Any], expected: dict[str, Any]) -> bool:
    return all(existing.get(key) == expected.get(key) for key in ("schema_version", "run_key", "list_sha256", "namespace", "taxon", "baseline_projection_sha256"))


def _verify_completed_manifest(output_dir: Path, manifest: dict[str, Any]) -> bool:
    version = manifest.get("import_version")
    run_key = manifest.get("run_key")
    result_hashes = manifest.get("result_sha256")
    if not isinstance(version, str) or not isinstance(run_key, str) or not isinstance(result_hashes, dict):
        return False
    directory = output_dir / "imports" / run_key / version
    if not directory.is_dir():
        return False
    return all(
        isinstance(expected, str) and (directory / name).is_file() and PathwayWorkflow.file_hash(directory / name) == expected
        for name, expected in result_hashes.items()
    )


def _validate_mapping(receipt: dict[str, Any], expected_foreground: int) -> None:
    foreground = receipt.get("foreground_mapping")
    if not isinstance(foreground, dict):
        raise ValueError("missing foreground_mapping")
    submitted = _nonnegative_int(foreground.get("submitted"), "foreground submitted count")
    mapped = _nonnegative_int(foreground.get("mapped"), "foreground mapped count")
    unique = _nonnegative_int(foreground.get("unique"), "foreground unique count")
    if submitted != expected_foreground:
        raise ValueError("foreground submitted count does not match prepared query")
    if mapped > submitted or unique > submitted:
        raise ValueError("foreground mapped or unique count exceeds submitted count")
    if unique > mapped:
        raise ValueError("foreground unique count exceeds mapped count")
    background = receipt.get("background_mapping")
    if not isinstance(background, dict):
        raise ValueError("missing background_mapping")
    mapped_background = background.get("mapped")
    if mapped_background is not None:
        _nonnegative_int(mapped_background, "background mapped count")
    _nonempty_string(background.get("reason"), "background mapping reason")


def _validate_timestamps(value: Any) -> None:
    if not isinstance(value, dict) or not value:
        raise ValueError("missing nonempty receipt timestamps")
    for name, timestamp in value.items():
        text = _nonempty_string(timestamp, f"receipt timestamp {name}")
        try:
            datetime.fromisoformat(text.replace("Z", "+00:00"))
        except ValueError as error:
            raise ValueError(f"malformed receipt timestamp {name}") from error


def _validate_receipt(receipt: dict[str, Any], manifest: dict[str, Any]) -> dict[str, dict[str, Any]]:
    required = ("schema_version", "run_key", "list_sha256", "namespace", "taxon", "population_confirmed", "timestamps", "selected_population", "foreground_mapping", "background_mapping", "interface", "reports", "export_sha256")
    missing = [field for field in required if field not in receipt]
    if missing:
        raise ValueError(f"missing receipt fields: {', '.join(missing)}")
    if receipt["schema_version"] != 1:
        raise ValueError("unsupported receipt schema_version")
    for field in ("run_key", "list_sha256", "namespace", "taxon"):
        if receipt[field] != manifest[field]:
            raise ValueError(f"receipt {field} does not match prepared DAVID run")
    if not _bool(receipt["population_confirmed"], "population_confirmed"):
        raise ValueError("population_confirmed must be true")
    _nonempty_string(receipt["selected_population"], "selected_population")
    _validate_timestamps(receipt["timestamps"])
    _validate_mapping(receipt, int(manifest["foreground_count"]))
    interface = receipt["interface"]
    if not isinstance(interface, dict):
        raise ValueError("missing interface")
    categories = interface.get("selected_categories")
    if not isinstance(categories, list) or not categories or any(not isinstance(item, str) or not item.strip() for item in categories):
        raise ValueError("missing nonempty interface selected_categories")
    settings = interface.get("settings")
    expected_settings = manifest["settings"]
    if not isinstance(settings, dict):
        raise ValueError("missing nonempty interface settings")
    if settings.get("chart") != expected_settings["chart"]:
        raise ValueError("interface chart Count and EASE must match prepared settings")
    clustering = settings.get("clustering")
    expected_clustering = expected_settings["clustering"]
    if (
        not isinstance(clustering, dict)
        or clustering.get("displayed") != expected_clustering["displayed"]
        or not isinstance(clustering.get("actual_categories_and_thresholds"), dict)
        or not clustering["actual_categories_and_thresholds"]
    ):
        raise ValueError("interface clustering settings must record actual categories and thresholds")
    reports = receipt["reports"]
    if not isinstance(reports, dict):
        raise ValueError("missing reports")
    validated: dict[str, dict[str, Any]] = {}
    for name in ("chart", "clusters"):
        state = reports.get(name)
        if not isinstance(state, dict):
            raise ValueError(f"missing {name} report state")
        if not _bool(state.get("complete"), f"{name} complete"):
            raise ValueError(f"{name} report is incomplete")
        if not _bool(state.get("complete_unfiltered"), f"{name} complete_unfiltered"):
            raise ValueError(f"{name} report must be complete_unfiltered")
        state["empty"] = _bool(state.get("empty"), f"{name} empty")
        validated[name] = state
    exports = receipt["export_sha256"]
    if (
        not isinstance(exports, dict)
        or any(
            not isinstance(exports.get(f"{name}.csv"), str)
            or re.fullmatch(r"[0-9a-fA-F]{64}", exports[f"{name}.csv"]) is None
            for name in ("chart", "clusters")
        )
    ):
        raise ValueError("missing export SHA-256 values")
    return validated


def _read_export(path: Path, required: set[str], report: str, completed_empty: bool) -> pd.DataFrame:
    if not path.is_file() or path.stat().st_size == 0:
        raise ValueError(f"missing or zero-byte {report} export")
    try:
        frame = pd.read_csv(path, encoding="utf-8-sig")
    except (OSError, UnicodeDecodeError, pd.errors.EmptyDataError, pd.errors.ParserError) as error:
        raise ValueError(f"malformed {report} export") from error
    missing = required - set(frame.columns)
    if missing:
        raise ValueError(f"malformed {report} export missing: {', '.join(sorted(missing))}")
    if frame.empty and not completed_empty:
        raise ValueError(f"header-only {report} export requires completed-empty receipt")
    if not frame.empty and completed_empty:
        raise ValueError(f"nonempty {report} export conflicts with completed-empty receipt")
    for _, row in frame.iterrows():
        for column in ("Category", "Term"):
            _nonempty_string(row[column], f"{report} {column}")
        count = _finite(row["Count"], f"{report} Count")
        if count < 0 or not count.is_integer():
            raise ValueError(f"invalid nonnegative integer {report} Count")
        _finite(row["P-Value"], f"{report} P-Value", unit_interval=True)
        _finite(row["Benjamini"], f"{report} Benjamini", unit_interval=True)
        if report == "clusters":
            _finite(row["Cluster Enrichment Score"], "cluster enrichment score")
            if _finite(row["Cluster Enrichment Score"], "cluster enrichment score") < 0:
                raise ValueError("invalid cluster enrichment score")
    frame["PROVIDER"] = "DAVID"
    return frame


def _text(value: Any) -> str:
    if value is None or value is pd.NA or (isinstance(value, float) and math.isnan(value)):
        return ""
    return str(value).strip()


def _stable_id(value: Any) -> str | None:
    candidate = _text(value)
    if candidate.lower().startswith("path:"):
        candidate = candidate[5:].strip()
    match = _ID_PATTERN.match(candidate)
    if match is None:
        return None
    identifier = match.group(1)
    return identifier.lower() if identifier.lower().startswith("hsa") else identifier.upper()
def _validate_report_categories(frames: tuple[pd.DataFrame, pd.DataFrame], selected_categories: list[str]) -> None:
    selected = {category.strip() for category in selected_categories}
    for frame in frames:
        categories = {str(value).strip() for value in frame["Category"].dropna()}
        unexpected = categories - selected
        if unexpected:
            raise ValueError("export Category is absent from recorded selected_categories")




def _ontology(source: Any, identifier: str | None) -> str:
    text = _text(source).upper()
    if identifier and identifier.startswith("GO:"):
        return "GO"
    if identifier and identifier.startswith("R-HSA-"):
        return "REACTOME"
    if identifier and identifier.startswith("WP"):
        return "WIKIPATHWAYS"
    if identifier and identifier.startswith("HP:"):
        return "HP"
    if identifier and identifier.lower().startswith("HSA"):
        return "KEGG"
    if "GO" in text:
        return "GO"
    if "REACT" in text:
        return "REACTOME"
    if "WIKI" in text or text.startswith("WP"):
        return "WIKIPATHWAYS"
    if "KEGG" in text:
        return "KEGG"
    if "HP" in text or "HPO" in text:
        return "HP"
    return "UNMAPPED"


def _json_rows(frame: pd.DataFrame) -> list[dict[str, Any]]:
    return [_canonical(record) for record in frame.to_dict(orient="records")]


def _gprofiler_genes(value: Any, mapping: dict[str, str]) -> set[str] | None:
    if not isinstance(value, str) or not value.strip():
        return set()
    symbols = {item.strip().upper() for item in value.split(";") if item.strip()}
    if not symbols:
        return set()
    genes = {mapping[symbol] for symbol in symbols if symbol in mapping}
    return genes if len(genes) == len(symbols) else None


def _david_genes(value: Any, submitted: set[str]) -> set[str] | None:
    if not isinstance(value, str) or not value.strip():
        return set()
    values = {item.strip() for item in re.split(r"[;,]", value) if item.strip()}
    return values if values.issubset(submitted) else None


def _term_key(source: Any, term: Any) -> tuple[str, str | None]:
    identifier = _stable_id(term)
    stated_ontology = _ontology(source, None)
    identifier_ontology = _ontology("", identifier)
    if identifier is not None and stated_ontology != "UNMAPPED" and stated_ontology != identifier_ontology:
        return stated_ontology, None
    return (identifier_ontology if identifier is not None else stated_ontology), identifier


def _comparison(p1: pd.DataFrame, chart: pd.DataFrame, mapping: dict[str, str], submitted: list[str]) -> pd.DataFrame:
    gprof: dict[tuple[str, str], pd.DataFrame] = {}
    gprof_keys = p1.apply(lambda row: _term_key(row.get("SOURCE"), row.get("TERM_ID")), axis=1)
    for key, rows in p1.groupby(gprof_keys, dropna=False):
        if key[1] is not None:
            gprof[key] = rows.copy()
    david: dict[tuple[str, str], pd.DataFrame] = {}
    david_keys = chart.apply(lambda row: _term_key(row.get("Category"), row.get("Term")), axis=1)
    for key, rows in chart.groupby(david_keys, dropna=False):
        if key[1] is not None:
            david[key] = rows.copy()
    records: list[dict[str, Any]] = []
    for key in sorted(set(gprof) | set(david)):
        g_rows, d_rows = gprof.get(key), david.get(key)
        if g_rows is not None and d_rows is not None:
            status = "shared"
        elif g_rows is not None:
            status = "gprofiler_only"
        else:
            status = "david_only"
        g_genes = _gprofiler_genes(g_rows.iloc[0].get("CONTRIBUTING_GENES"), mapping) if g_rows is not None else set()
        d_column = "User Ids" if d_rows is not None and "User Ids" in d_rows else None
        d_genes = _david_genes(d_rows.iloc[0][d_column], set(submitted)) if d_column else set()
        note = "stable identifier match; provider statistics remain separate"
        jaccard: float | None = None
        if status == "shared":
            if g_genes is None or d_genes is None:
                note = "stable identifier match; contributor mapping unavailable"
            elif g_genes and d_genes:
                jaccard = len(g_genes & d_genes) / len(g_genes | d_genes)
            else:
                note = "stable identifier match; no nonempty comparable contributor sets"
        records.append({
            "ONTOLOGY": key[0], "TERM_ID": key[1], "MATCH_STATUS": status,
            "GPROFILER_ROWS_JSON": json.dumps(_json_rows(g_rows) if g_rows is not None else [], sort_keys=True),
            "DAVID_ROWS_JSON": json.dumps(_json_rows(d_rows) if d_rows is not None else [], sort_keys=True),
            "GENE_JACCARD": jaccard, "COMPARABILITY_NOTE": note,
        })
    for index, row in p1.loc[gprof_keys.map(lambda key: key[1] is None)].iterrows():
        records.append({"ONTOLOGY": _term_key(row.get("SOURCE"), row.get("TERM_ID"))[0], "TERM_ID": None, "MATCH_STATUS": "unmapped_term", "GPROFILER_ROWS_JSON": json.dumps([_canonical(row.to_dict())], sort_keys=True), "DAVID_ROWS_JSON": "[]", "GENE_JACCARD": None, "COMPARABILITY_NOTE": "g:Profiler term has no supported stable identifier"})
    for index, row in chart.loc[david_keys.map(lambda key: key[1] is None)].iterrows():
        records.append({"ONTOLOGY": _ontology(row.get("Category"), None), "TERM_ID": None, "MATCH_STATUS": "unmapped_term", "GPROFILER_ROWS_JSON": "[]", "DAVID_ROWS_JSON": json.dumps([_canonical(row.to_dict())], sort_keys=True), "GENE_JACCARD": None, "COMPARABILITY_NOTE": "DAVID title has no supported stable identifier"})
    return pd.DataFrame(records, columns=["ONTOLOGY", "TERM_ID", "MATCH_STATUS", "GPROFILER_ROWS_JSON", "DAVID_ROWS_JSON", "GENE_JACCARD", "COMPARABILITY_NOTE"])


def _write_report(path: Path, receipt: dict[str, Any], chart: pd.DataFrame, clusters: pd.DataFrame, comparison: pd.DataFrame) -> None:
    mapping = receipt["foreground_mapping"]
    background = receipt["background_mapping"]
    lines = [
        "# DAVID addendum", "", "DAVID is supplementary annotation evidence; it does not alter P0–P7 candidate or drug selection.", "",
        "## Mapping and background QC", "",
        f"- Foreground submitted/mapped/unique: {mapping['submitted']}/{mapping['mapped']}/{mapping['unique']}",
        f"- Background mapped: {background.get('mapped')}; reason: {background['reason']}",
        f"- Confirmed selected population: {receipt['selected_population']}", "",
        "## Provider-specific results", "",
        f"- DAVID chart terms: {len(chart)}; Benjamini ≤ 0.05: {int(pd.to_numeric(chart.get('Benjamini'), errors='coerce').le(0.05).sum())}",
        f"- DAVID clusters: {len(clusters)}", f"- Shared stable-ID terms: {int(comparison['MATCH_STATUS'].eq('shared').sum())}",
        f"- DAVID-only terms: {int(comparison['MATCH_STATUS'].eq('david_only').sum())}",
        f"- g:Profiler-only terms: {int(comparison['MATCH_STATUS'].eq('gprofiler_only').sum())}", "",
        "## Limitations", "", "- DAVID and g:Profiler statistics remain provider-specific; shared annotations are not independent replication.", "- Native FDR and cluster enrichment score are retained without treating either as Benjamini BH.", "- Existing P2–P7 artifacts remain unchanged context; this addendum does not invoke downstream target or drug selection.",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _import_exports(output_dir: Path, manifest: dict[str, Any], p0: pd.DataFrame, p1: pd.DataFrame, exports_dir: Path) -> tuple[str, dict[str, str]]:
    receipt_path = exports_dir / "receipt.json"
    receipt = _read_json(receipt_path)
    reports = _validate_receipt(receipt, manifest)
    receipt_hash = PathwayWorkflow.file_hash(receipt_path)
    version_dir = output_dir / "imports" / str(manifest["run_key"]) / receipt_hash
    if version_dir.exists():
        completed = _read_json(version_dir / "import_manifest.json") if (version_dir / "import_manifest.json").is_file() else {}
        hashes = completed.get("result_sha256") if isinstance(completed, dict) else None
        if isinstance(hashes, dict) and all((version_dir / name).is_file() and PathwayWorkflow.file_hash(version_dir / name) == expected for name, expected in hashes.items()):
            return receipt_hash, hashes
        raise ValueError("immutable DAVID import version is incomplete or corrupt")
    chart_path, clusters_path = exports_dir / "chart.csv", exports_dir / "clusters.csv"
    hashes = receipt["export_sha256"]
    if not chart_path.is_file() or not clusters_path.is_file() or PathwayWorkflow.file_hash(chart_path) != hashes["chart.csv"] or PathwayWorkflow.file_hash(clusters_path) != hashes["clusters.csv"]:
        raise ValueError("DAVID export SHA-256 does not match receipt")
    chart = _read_export(chart_path, _CHART_COLUMNS, "chart", reports["chart"]["empty"])
    clusters = _read_export(clusters_path, _CLUSTER_COLUMNS, "clusters", reports["clusters"]["empty"])
    _validate_report_categories((chart, clusters), receipt["interface"]["selected_categories"])
    foreground, _, symbol_mapping = _resolved_lists(p0)
    comparison = _comparison(p1, chart, symbol_mapping, foreground)
    version_dir.mkdir(parents=True)
    try:
        shutil.copyfile(receipt_path, version_dir / "receipt.json")
        PathwayWorkflow.atomic_write_csv(chart, version_dir / "david_chart.csv")
        PathwayWorkflow.atomic_write_csv(clusters, version_dir / "david_clusters.csv")
        PathwayWorkflow.atomic_write_csv(comparison, version_dir / "provider_term_comparison.csv")
        _write_report(version_dir / "david_comparison_report.md", receipt, chart, clusters, comparison)
        result_hashes = {name: PathwayWorkflow.file_hash(version_dir / name) for name in ("receipt.json", "david_chart.csv", "david_clusters.csv", "provider_term_comparison.csv", "david_comparison_report.md")}
        _write_json(version_dir / "import_manifest.json", {"schema_version": 1, "run_key": manifest["run_key"], "receipt_sha256": receipt_hash, "result_sha256": result_hashes})
    except Exception:
        # Leave no partially trusted immutable import directory.
        for path in version_dir.iterdir():
            path.unlink()
        version_dir.rmdir()
        raise
    return receipt_hash, result_hashes


def david_run_summary(output_dir: Path, exit_code: int, *, reason: str | None = None) -> dict[str, Any]:
    """Return the narrow DAVID status contract consumed by the workflow runner."""
    manifest_path = output_dir / "david_manifest.json"
    manifest = _read_json(manifest_path) if manifest_path.is_file() else {}
    status = {0: "completed", 1: "failed", 2: "awaiting_web_results"}.get(exit_code, "failed")
    result: dict[str, Any] = {
        "status": status,
        "run_key": manifest.get("run_key"),
        "import_version": manifest.get("import_version"),
        "external_output_dir": str(output_dir.resolve()),
    }
    if reason is not None:
        result["reason"] = reason
    return result


def run_david_addition(baseline_dir: Path, output_dir: Path, exports_dir: Path | None, *, force: bool = False) -> int:
    """Prepare or import one verified native DAVID addendum without network access."""
    baseline_dir, output_dir = Path(baseline_dir), Path(output_dir)
    _validate_output_path(baseline_dir, output_dir)
    p0, p1 = _project_csvs(baseline_dir)
    foreground, background, _ = _resolved_lists(p0)
    expected = _prepared_manifest(p0, p1, foreground, background)
    manifest_path = output_dir / "david_manifest.json"
    existing = _read_json(manifest_path) if manifest_path.is_file() else None
    if existing is not None and not _matching_manifest(existing, expected):
        if not force:
            raise ValueError("STALE_OR_INVALID_PROVENANCE")
        existing = None
    if existing is not None and existing.get("status") == "completed" and exports_dir is None:
        if _verify_completed_manifest(output_dir, existing):
            return 0
        if not force:
            raise ValueError("STALE_OR_INVALID_PROVENANCE")
    if existing is None:
        output_dir.mkdir(parents=True, exist_ok=True)
        (output_dir / "foreground.txt").write_text("\n".join(foreground) + "\n", encoding="utf-8")
        (output_dir / "background.txt").write_text("\n".join(background) + "\n", encoding="utf-8")
        _write_json(manifest_path, expected)
        existing = expected
    else:
        expected_lists = _sha256_bytes((output_dir / "foreground.txt").read_bytes() + b"\0" + (output_dir / "background.txt").read_bytes()) if (output_dir / "foreground.txt").is_file() and (output_dir / "background.txt").is_file() else None
        if expected_lists != existing.get("list_sha256"):
            raise ValueError("STALE_OR_INVALID_PROVENANCE")
    if exports_dir is None:
        return 2
    receipt_hash, result_hashes = _import_exports(output_dir, existing, p0, p1, Path(exports_dir))
    if existing.get("status") != "completed" or existing.get("import_version") != receipt_hash:
        completed = dict(existing)
        completed.update({"status": "completed", "import_version": receipt_hash, "result_sha256": result_hashes, "completed_at": datetime.now(UTC).isoformat()})
        _write_json(manifest_path, completed)
    return 0
