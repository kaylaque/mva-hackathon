"""Offline, provenance-bound comparison of saved DAVID and g:Profiler results."""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
import shutil
import subprocess
from typing import Any

import pandas as pd

from mva_hackathon.david_workflow import (
    _CHART_COLUMNS,
    _CLUSTER_COLUMNS,
    _P0_COLUMNS,
    _P1_COLUMNS,
    _comparison,
    _complete_contributor_union,
    _ontology,
    _prepared_manifest,
    _project_csvs,
    _resolved_lists,
    _term_key,
    _text,
    _validate_receipt,
)
from mva_hackathon.mva_gene_sets import (
    CHROMOSOME_COHESIN_GENES,
    CORE_MVA_CAUSAL_GENES,
    SAC_CORE_GENES,
)
from mva_hackathon.pathway_workflow import PathwayWorkflow


_REQUIRED_IMPORTS = (
    "receipt.json",
    "david_chart.csv",
    "david_clusters.csv",
    "provider_term_comparison.csv",
    "david_comparison_report.md",
)
_HEX64 = set("0123456789abcdef")
_MECHANISM_PANEL = CORE_MVA_CAUSAL_GENES | SAC_CORE_GENES | CHROMOSOME_COHESIN_GENES
_EVIDENCE_COLUMNS = (
    "EVIDENCE_ID", "ONTOLOGY", "TERM_ID", "MATCH_STATUS", "GPROFILER_ROWS_JSON", "DAVID_ROWS_JSON",
    "GPROFILER_SIGNIFICANCE", "DAVID_SIGNIFICANCE", "EVIDENCE_TIER", "GENE_JACCARD",
    "GPROFILER_CONTRIBUTOR_STATUS", "DAVID_CONTRIBUTOR_STATUS", "CONTRIBUTOR_STATUS", "MVA_CONTEXT_JSON",
)
_GENE_COLUMNS = ("EVIDENCE_ID", "ONTOLOGY", "TERM_ID", "PROVIDER", "GENE_SYMBOL", "MEMBERSHIP_STATUS", "MVA_GENE_SETS")
_CLUSTER_OUTPUT_COLUMNS = (
    "CLUSTER_ID", "CLUSTER_ENRICHMENT_SCORE", "ANNOTATION_ROW_COUNT", "MAPPED_TERM_COUNT",
    "UNMAPPED_TERM_COUNT", "SHARED_TERM_COUNT", "DAVID_SIGNIFICANT_TERM_COUNT", "TERM_IDS_JSON",
    "TERM_NAMES_JSON", "CONTEXT_GENES_JSON", "CONTEXT_MEMBERSHIP_SOURCE",
)
_DISCORDANCE_COLUMNS = (
    "ONTOLOGY", "TERM_ID", "MATCH_STATUS", "GPROFILER_SIGNIFICANCE", "DAVID_SIGNIFICANCE",
    "UNRESOLVED_MATCH_REASON", "REVIEW_NOTE",
)
_OUTPUT_FILES = (
    "integrated_pathway_evidence.csv",
    "gene_pathway_support.csv",
    "david_cluster_summary.csv",
    "provider_discordance_review.csv",
    "integrated_pathway_report.md",
    "postanalysis_provenance.json",
)


def _fail(message: str, error: Exception | None = None) -> None:
    raise ValueError(message) from error

def _analysis_threshold(value: float) -> float:
    try:
        threshold = float(value)
    except (TypeError, ValueError) as error:
        _fail("analysis threshold must be finite and in (0, 1]", error)
    if not math.isfinite(threshold) or not 0 < threshold <= 1:
        _fail("analysis threshold must be finite and in (0, 1]")
    return threshold


def _sha256(path: Path) -> str:
    try:
        return hashlib.sha256(path.read_bytes()).hexdigest()
    except OSError as error:
        _fail(f"unreadable artifact: {path.name}", error)


def _json(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as error:
        _fail(f"invalid JSON artifact: {path.name}", error)
    if not isinstance(value, dict):
        _fail(f"invalid JSON object: {path.name}")
    return value


def _is_hex64(value: Any) -> bool:
    return isinstance(value, str) and len(value) == 64 and set(value.lower()) <= _HEX64


def _git_root(path: Path) -> Path | None:
    completed = subprocess.run(
        ["git", "-C", str(path), "rev-parse", "--show-toplevel"],
        capture_output=True,
        text=True,
        check=False,
    )
    if completed.returncode != 0:
        return None
    return Path(completed.stdout.strip()).resolve()


def _assert_existing_file(path: Path, root: Path, label: str) -> Path:
    resolved = path.resolve()
    if not resolved.is_file() or not resolved.is_relative_to(root):
        _fail(f"invalid {label} path")
    return resolved


def _assert_nonignored(path: Path, repo_root: Path | None) -> None:
    if repo_root is None or not path.resolve().is_relative_to(repo_root):
        return
    completed = subprocess.run(
        ["git", "check-ignore", "--quiet", str(path.resolve().relative_to(repo_root))],
        cwd=repo_root,
        capture_output=True,
        text=True,
        check=False,
    )
    if completed.returncode == 0:
        _fail("ignored export path is forbidden")
    if completed.returncode != 1:
        _fail("cannot validate export path")


def _validate_output(output_dir: Path, baseline_dir: Path, david_dir: Path) -> Path:
    output = output_dir.resolve()
    if output.exists():
        _fail("output directory already exists")
    if any(output == item or output.is_relative_to(item) or item.is_relative_to(output) for item in (baseline_dir, david_dir)):
        _fail("output directory must be disjoint from inputs")
    parent = output.parent
    if not parent.is_dir():
        _fail("output parent does not exist")
    if _git_root(parent) is not None:
        _fail("output directory must be outside a Git repository")
    return output


def _read_baseline(baseline_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    root = baseline_dir.resolve()
    repo = _git_root(root)
    if repo is not None and root.is_relative_to(repo) and not root.is_relative_to(repo / "data" / "processed"):
        _fail("baseline inside a repository must be under data/processed")
    p0_path, p1_path = root / "p0_gene_qc.csv", root / "p1_enrichment_results.csv"
    for path in (p0_path, p1_path):
        _assert_existing_file(path, root, "baseline artifact")
        if repo is not None and path.resolve().is_relative_to(repo) and not path.resolve().is_relative_to(repo / "data" / "processed"):
            _fail("baseline artifact escapes processed root")
    try:
        p0, p1 = _project_csvs(root)
        foreground, background, mapping = _resolved_lists(p0)
    except (OSError, ValueError, pd.errors.ParserError) as error:
        _fail("invalid projected baseline artifacts", error)
    return p0, p1, _prepared_manifest(p0, p1, foreground, background) | {"_mapping": mapping, "_foreground": foreground}


def _read_import(david_dir: Path, baseline: dict[str, Any]) -> tuple[dict[str, Any], pd.DataFrame, pd.DataFrame, Path, dict[str, str]]:
    root = david_dir.resolve()
    repo = _git_root(root)
    manifest_path = _assert_existing_file(root / "david_manifest.json", root, "DAVID manifest")
    _assert_nonignored(manifest_path, repo)
    manifest = _json(manifest_path)
    if manifest.get("status") != "completed" or not _is_hex64(manifest.get("run_key")) or not _is_hex64(manifest.get("import_version")):
        _fail("DAVID manifest is not a completed immutable import")
    if manifest.get("baseline_projection_sha256") != baseline["baseline_projection_sha256"]:
        _fail("DAVID manifest baseline projection mismatch")
    result_hashes = manifest.get("result_sha256")
    if not isinstance(result_hashes, dict) or set(result_hashes) != set(_REQUIRED_IMPORTS):
        _fail("DAVID manifest has incomplete result hash map")
    version_dir = root / "imports" / manifest["run_key"] / manifest["import_version"]
    if not version_dir.is_dir():
        _fail("DAVID immutable import is missing")
    import_manifest_path = _assert_existing_file(version_dir / "import_manifest.json", version_dir, "import manifest")
    _assert_nonignored(import_manifest_path, repo)
    import_manifest = _json(import_manifest_path)
    hashes = import_manifest.get("result_sha256")
    if import_manifest.get("run_key") != manifest["run_key"] or import_manifest.get("receipt_sha256") != manifest["import_version"]:
        _fail("DAVID import identity mismatch")
    if not isinstance(hashes, dict) or hashes != result_hashes or set(hashes) != set(_REQUIRED_IMPORTS):
        _fail("DAVID import hash map mismatch")
    for name in _REQUIRED_IMPORTS:
        if not _is_hex64(hashes[name]):
            _fail("invalid DAVID artifact hash")
        artifact = _assert_existing_file(version_dir / name, version_dir, "DAVID artifact")
        _assert_nonignored(artifact, repo)
        if _sha256(artifact) != hashes[name]:
            _fail("DAVID artifact hash mismatch")
    receipt = _json(version_dir / "receipt.json")
    if receipt.get("run_key") != manifest["run_key"]:
        _fail("DAVID receipt run identity mismatch")
    try:
        _validate_receipt(receipt, manifest)
    except ValueError as error:
        _fail("DAVID receipt contract is not satisfied", error)
    try:
        chart = pd.read_csv(version_dir / "david_chart.csv", encoding="utf-8-sig")
        clusters = pd.read_csv(version_dir / "david_clusters.csv", encoding="utf-8-sig")
    except (OSError, UnicodeDecodeError, pd.errors.ParserError) as error:
        _fail("malformed DAVID report", error)
    if not _CHART_COLUMNS.issubset(chart.columns) or not _CLUSTER_COLUMNS.issubset(clusters.columns):
        _fail("DAVID report missing required columns")
    return receipt, chart, clusters, version_dir, {name: hashes[name] for name in _REQUIRED_IMPORTS}


def _state(rows: list[dict[str, Any]], column: str, threshold: float) -> str:
    """Summarize a provider statistic column without hiding malformed values.

    Every present value is validated first, so one missing value never masks a
    malformed value in another row.  Missing data stays unknown rather than
    being folded into not_significant.
    """
    if not rows:
        return "absent"
    values: list[float] = []
    missing = False
    for row in rows:
        value = row.get(column)
        if value is None or (isinstance(value, float) and math.isnan(value)) or str(value).strip() == "":
            missing = True
            continue
        try:
            number = float(value)
        except (TypeError, ValueError) as error:
            _fail(f"invalid {column}", error)
        if not math.isfinite(number) or not 0 <= number <= 1:
            _fail(f"invalid {column}")
        values.append(number)
    if missing or not values:
        return "unknown"
    significant = [value <= threshold for value in values]
    return "significant" if all(significant) else "not_significant" if not any(significant) else "mixed"


def _tier(match: str, g_state: str, d_state: str) -> str:
    if match == "unmapped_term":
        return "review_required"
    if "unknown" in (g_state, d_state) or "mixed" in (g_state, d_state):
        return "review_required"
    if g_state == d_state == "significant":
        return "both_significant"
    if match == "shared" and {g_state, d_state} == {"significant", "not_significant"}:
        return "shared_one_significant"
    if match == "gprofiler_only" and g_state == "significant":
        return "gprofiler_only_significant"
    if match == "david_only" and d_state == "significant":
        return "david_only_significant"
    return "context_only"


def _context(genes: set[str]) -> list[str]:
    sets = {
        "CORE_MVA_CAUSAL_GENES": CORE_MVA_CAUSAL_GENES,
        "SAC_CORE_GENES": SAC_CORE_GENES,
        "CHROMOSOME_COHESIN_GENES": CHROMOSOME_COHESIN_GENES,
    }
    return [name for name, members in sets.items() if genes & members]


def _rows(value: str) -> list[dict[str, Any]]:
    parsed = json.loads(value)
    if not isinstance(parsed, list) or not all(isinstance(row, dict) for row in parsed):
        _fail("invalid provider row payload")
    return parsed


def _membership(rows: list[dict[str, Any]], provider: str, mapping: dict[str, str], submitted: set[str]) -> tuple[set[str], str]:
    """Return observed contributor identifiers and the provider's coverage status.

    ``complete``    every contributing row resolved to a nonempty identifier set
    ``partial``     at least one row resolved while another did not
    ``unavailable`` rows exist but no contributor set could be resolved
    ``absent``      the provider produced no row for this term
    """
    if not rows:
        return set(), "absent"
    observed: set[str] = set()
    unresolved = False
    for row in rows:
        genes = _complete_contributor_union([row], provider, mapping, submitted)
        if genes:
            observed |= genes
        else:
            unresolved = True
    if not observed:
        return set(), "unavailable"
    return observed, "partial" if unresolved else "complete"


def _combined_status(g_status: str, d_status: str) -> str:
    """Combined coverage is complete only for two complete nonempty provider sets."""
    states = {g_status, d_status}
    if states == {"complete"}:
        return "complete"
    if states <= {"complete", "partial"}:
        return "partial"
    if states == {"absent"}:
        return "absent"
    return "unavailable"


def _reverse_symbols(mapping: dict[str, str]) -> dict[str, str]:
    """One-to-one reverse mapping; aliases collapsing onto one identifier are dropped."""
    grouped: dict[str, set[str]] = {}
    for symbol, ensembl in mapping.items():
        grouped.setdefault(ensembl, set()).add(symbol)
    return {ensembl: next(iter(symbols)) for ensembl, symbols in grouped.items() if len(symbols) == 1}


def _evidence_digest(record: dict[str, Any]) -> str:
    """Stable identity for one comparison record, distinct even when TERM_ID is null."""
    payload = json.dumps(
        {
            "ontology": record["ONTOLOGY"],
            "term_id": record["TERM_ID"],
            "match": record["MATCH_STATUS"],
            "g": _rows(record["GPROFILER_ROWS_JSON"]),
            "d": _rows(record["DAVID_ROWS_JSON"]),
        },
        sort_keys=True,
        default=str,
    )
    return "EV-" + hashlib.sha256(payload.encode()).hexdigest()[:24]


def _cluster_context(frame: pd.DataFrame, mapping: dict[str, str], submitted: set[str]) -> tuple[list[str], str]:
    """Mechanism-panel genes observed directly in a cluster's own annotation rows."""
    if "User Ids" not in frame.columns:
        return [], "unavailable"
    observed: set[str] = set()
    for value in frame["User Ids"]:
        genes = _complete_contributor_union([{"User Ids": value}], "DAVID", mapping, submitted)
        if genes:
            observed |= genes
    if not observed:
        return [], "unavailable"
    reverse = _reverse_symbols(mapping)
    symbols = sorted({reverse[gene] for gene in observed if gene in reverse and reverse[gene] in _MECHANISM_PANEL})
    return symbols, "cluster_observed"


def _build_evidence(
    p1: pd.DataFrame,
    chart: pd.DataFrame,
    mapping: dict[str, str],
    submitted: set[str],
    david_threshold: float,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    reverse_mapping = _reverse_symbols(mapping)
    comparison = _comparison(p1, chart, mapping, sorted(submitted))
    evidence: list[dict[str, Any]] = []
    genes: list[dict[str, Any]] = []
    seen: dict[str, int] = {}
    for _, source in comparison.iterrows():
        g_rows, d_rows = _rows(source["GPROFILER_ROWS_JSON"]), _rows(source["DAVID_ROWS_JSON"])
        g_set, g_status = _membership(g_rows, "g:Profiler", mapping, submitted)
        d_set, d_status = _membership(d_rows, "DAVID", mapping, submitted)
        g_state = _state(g_rows, "BH_Q", 0.05) if g_rows else "not_comparable"
        d_state = _state(d_rows, "Benjamini", david_threshold) if d_rows else "not_comparable"
        # Overlap is only defined for two complete nonempty provider sets.
        jaccard: float | None = None
        if source["MATCH_STATUS"] == "shared" and g_status == d_status == "complete":
            jaccard = len(g_set & d_set) / len(g_set | d_set)
        digest = _evidence_digest(source)
        seen[digest] = seen.get(digest, 0) + 1
        evidence_id = digest if seen[digest] == 1 else f"{digest}-{seen[digest]}"
        observed = {reverse_mapping[gene] for gene in (g_set | d_set) if gene in reverse_mapping}
        evidence.append({
            "EVIDENCE_ID": evidence_id,
            "ONTOLOGY": source["ONTOLOGY"], "TERM_ID": source["TERM_ID"], "MATCH_STATUS": source["MATCH_STATUS"],
            "GPROFILER_ROWS_JSON": source["GPROFILER_ROWS_JSON"], "DAVID_ROWS_JSON": source["DAVID_ROWS_JSON"],
            "GPROFILER_SIGNIFICANCE": g_state, "DAVID_SIGNIFICANCE": d_state,
            "EVIDENCE_TIER": _tier(source["MATCH_STATUS"], g_state, d_state), "GENE_JACCARD": jaccard,
            "GPROFILER_CONTRIBUTOR_STATUS": g_status, "DAVID_CONTRIBUTOR_STATUS": d_status,
            "CONTRIBUTOR_STATUS": _combined_status(g_status, d_status),
            "MVA_CONTEXT_JSON": json.dumps(_context(observed)),
        })
        for provider, gene_set, membership_status in (("g:Profiler", g_set, g_status), ("DAVID", d_set, d_status)):
            for gene in sorted(gene_set):
                symbol = reverse_mapping.get(gene)
                if symbol is None:
                    continue
                genes.append({
                    "EVIDENCE_ID": evidence_id, "ONTOLOGY": source["ONTOLOGY"], "TERM_ID": source["TERM_ID"],
                    "PROVIDER": provider, "GENE_SYMBOL": symbol,
                    "MEMBERSHIP_STATUS": "observed" if membership_status == "complete" else "partial",
                    "MVA_GENE_SETS": json.dumps(_context({symbol})),
                })
    return pd.DataFrame(evidence, columns=_EVIDENCE_COLUMNS), pd.DataFrame(genes, columns=_GENE_COLUMNS)


def _cluster_summary(clusters: pd.DataFrame, evidence: pd.DataFrame, mapping: dict[str, str], submitted: set[str]) -> pd.DataFrame:
    """One row per distinct nonempty DAVID cluster, using Full Term identifiers.

    Mapped term counts aggregate unique stable-ID/category pairs; significance is read
    from each cluster's own Benjamini-bearing rows, never from chart significance.
    """
    records: list[dict[str, Any]] = []
    mapped = evidence.loc[evidence["TERM_ID"].notna()].set_index(["ONTOLOGY", "TERM_ID"]) if not evidence.empty else None
    for cluster_id, frame in clusters.groupby("Cluster", dropna=True):
        label = str(cluster_id).strip()
        if not label or label.lower() == "nan":
            continue
        scores = pd.to_numeric(frame["Cluster Enrichment Score"], errors="coerce")
        if scores.isna().any() or not (scores >= 0).all() or not all(math.isfinite(value) for value in scores) or scores.nunique() != 1:
            _fail("inconsistent cluster enrichment score")
        term_column = "Full Term" if "Full Term" in frame.columns else "Term"
        keys = [_term_key(category, value) for category, value in zip(frame.get("Category"), frame.get(term_column))]
        pairs = sorted({(_text(category), key[1]) for category, key in zip(frame.get("Category"), keys) if key[1] is not None})
        terms = sorted({(key[0], key[1]) for key in keys if key[1] is not None})
        names = sorted({_text(value) for value in frame.get(term_column, pd.Series(dtype=object)) if _text(value)})
        shared = david_significant = 0
        if mapped is not None:
            for key in terms:
                if key not in mapped.index:
                    continue
                selected = mapped.loc[key]
                selected = selected if isinstance(selected, pd.DataFrame) else selected.to_frame().T
                shared += int(selected["MATCH_STATUS"].eq("shared").sum())
                david_significant += int(selected["DAVID_SIGNIFICANCE"].eq("significant").sum())
        context_genes, context_source = _cluster_context(frame, mapping, submitted)
        records.append({
            "CLUSTER_ID": label, "CLUSTER_ENRICHMENT_SCORE": float(scores.iloc[0]), "ANNOTATION_ROW_COUNT": len(frame),
            "MAPPED_TERM_COUNT": len(pairs), "UNMAPPED_TERM_COUNT": int(sum(key[1] is None for key in keys)),
            "SHARED_TERM_COUNT": shared, "DAVID_SIGNIFICANT_TERM_COUNT": david_significant,
            "TERM_IDS_JSON": json.dumps([term_id for _, term_id in pairs]), "TERM_NAMES_JSON": json.dumps(names),
            "CONTEXT_GENES_JSON": json.dumps(context_genes), "CONTEXT_MEMBERSHIP_SOURCE": context_source,
        })
    return pd.DataFrame(records, columns=_CLUSTER_OUTPUT_COLUMNS)


def _discordance(evidence: pd.DataFrame, categories: list[str]) -> pd.DataFrame:
    """Review rows for unmapped annotations and for providers that disagree."""
    rows: list[dict[str, Any]] = []
    for _, item in evidence.iterrows():
        states = {item["GPROFILER_SIGNIFICANCE"], item["DAVID_SIGNIFICANCE"]}
        unmapped = item["MATCH_STATUS"] == "unmapped_term"
        observed_states = states - {"not_comparable"}
        if not unmapped and item["EVIDENCE_TIER"] not in {"gprofiler_only_significant", "david_only_significant", "shared_one_significant"} and len(observed_states) <= 1:
            continue
        reason = (
            "unmapped: the provider term carries no supported stable identifier, so no counterpart state exists"
            if unmapped
            else "mapped on a stable identifier; the two providers differ in coverage, effective background, thresholds or mappings"
        )
        rows.append({
            "ONTOLOGY": item["ONTOLOGY"], "TERM_ID": item["TERM_ID"], "MATCH_STATUS": item["MATCH_STATUS"],
            "GPROFILER_SIGNIFICANCE": item["GPROFILER_SIGNIFICANCE"], "DAVID_SIGNIFICANCE": item["DAVID_SIGNIFICANCE"],
            "UNRESOLVED_MATCH_REASON": reason,
            "REVIEW_NOTE": "Unresolved: provider coverage, mappings, effective background, and export thresholds may differ; DAVID categories: " + ", ".join(categories),
        })
    return pd.DataFrame(rows, columns=_DISCORDANCE_COLUMNS)


def _write_report(
    path: Path,
    receipt: dict[str, Any],
    evidence: pd.DataFrame,
    clusters: pd.DataFrame,
    genes: pd.DataFrame,
    discordance: pd.DataFrame,
    provenance: dict[str, Any],
) -> None:
    interface = receipt.get("interface", {})
    categories = interface.get("selected_categories", []) if isinstance(interface, dict) else []
    lines = [
        "# Integrated pathway report", "",
        "Research-only post-analysis of saved g:Profiler and DAVID exports. Descriptive annotation is not",
        "efficacy evidence, not a diagnostic finding, and not independent replication.", "",
        "## Definitions used", "",
        "- **Pathway-supported**: a provider reported the term below its own adjusted-q threshold; each provider is judged on its own statistics.",
        f"- **DAVID sensitivity threshold**: Benjamini <= {provenance['analysis_threshold']}; g:Profiler BH_Q remains at 0.05.",
        "- **Nonspecific for an MVA mechanism**: the term is not itself a spindle-assembly-checkpoint, chromosome-segregation or cohesin process, so it is context rather than mechanism.",
        "- **Modification annotation**: membership in a post-translational-modification term names members, not measured activity.",
        "- **not_comparable**: the other provider produced no row for this term, so no cross-provider statement is made.",
        "- **Jaccard**: gene overlap between two complete nonempty contributor sets; unavailable sets report nothing rather than an overlap of zero.",
        "",
        "## Verified provenance identities", "",
        f"- run_key: `{provenance['run_key']}`",
        f"- import_version: `{provenance['import_version']}`",
        f"- baseline projection SHA-256: `{provenance['baseline_projections']}`",
        f"- analysis threshold: {provenance['analysis_threshold']}",
        f"- DAVID categories: {', '.join(str(item) for item in categories) or 'none recorded'}",
        f"- distinct nonempty clusters: {clusters['CLUSTER_ID'].nunique() if not clusters.empty else 0}",
        f"- cluster annotation rows: {int(clusters['ANNOTATION_ROW_COUNT'].sum()) if not clusters.empty else 0}",
        "",
        "## Native provider results", "",
        "### DAVID-significant rows",
        "",
    ]
    david_significant = evidence.loc[evidence["DAVID_SIGNIFICANCE"].eq("significant")]
    if david_significant.empty:
        lines.append("- none: no DAVID row met its own Benjamini threshold")
    for _, row in david_significant.iterrows():
        rows = _rows(row["DAVID_ROWS_JSON"])
        names = [str(entry.get("Term") or entry.get("Full Term") or "") for entry in rows]
        native = [str(entry.get("Benjamini")) for entry in rows]
        lines.append(f"- {row['TERM_ID']}: {'; '.join(names)} (native Benjamini {', '.join(native)}; g:Profiler {row['GPROFILER_SIGNIFICANCE']})")
    gprofiler_significant = evidence.loc[evidence["GPROFILER_SIGNIFICANCE"].eq("significant")]
    lines.extend(["", "### g:Profiler-significant unique terms", ""])
    if gprofiler_significant.empty:
        lines.append("- none: no g:Profiler row met its own adjusted-q threshold")
    for term_id in sorted({str(value) for value in gprofiler_significant["TERM_ID"].fillna("") if str(value)}):
        lines.append(f"- {term_id}")
    lines.extend(["", "### Unmapped significant annotations requiring review", ""])
    unmapped = evidence.loc[evidence["MATCH_STATUS"].eq("unmapped_term")]
    if unmapped.empty:
        lines.append("- none")
    for _, row in unmapped.iterrows():
        lines.append(f"- {row['ONTOLOGY']} {row['EVIDENCE_ID']}: {row['GPROFILER_SIGNIFICANCE']} / {row['DAVID_SIGNIFICANCE']}")
    lines.extend([
        "", "## Per-provider contributor coverage", "",
        "| Provider | complete | partial | unavailable | absent |",
        "| --- | --- | --- | --- | --- |",
    ])
    for provider, column in (("g:Profiler", "GPROFILER_CONTRIBUTOR_STATUS"), ("DAVID", "DAVID_CONTRIBUTOR_STATUS")):
        counts = evidence[column].value_counts() if not evidence.empty else pd.Series(dtype=int)
        lines.append(f"| {provider} | {int(counts.get('complete', 0))} | {int(counts.get('partial', 0))} | {int(counts.get('unavailable', 0))} | {int(counts.get('absent', 0))} |")
    lines.extend(["", f"- combined complete (shared term, two complete nonempty sets): {int(evidence['CONTRIBUTOR_STATUS'].eq('complete').sum()) if not evidence.empty else 0}", ""])
    lines.extend(["## Top clusters", "", "| Cluster | Score | Rows | Mapped terms | Shared | Context genes | Context source |", "| --- | --- | --- | --- | --- | --- | --- |"])
    top = clusters.head(20)
    if top.empty:
        lines.append("| none | | | | | | |")
    for _, row in top.iterrows():
        lines.append(
            f"| {row['CLUSTER_ID']} | {row['CLUSTER_ENRICHMENT_SCORE']} | {row['ANNOTATION_ROW_COUNT']} | "
            f"{row['MAPPED_TERM_COUNT']} | {row['SHARED_TERM_COUNT']} | {', '.join(json.loads(row['CONTEXT_GENES_JSON'])) or 'none'} | {row['CONTEXT_MEMBERSHIP_SOURCE']} |"
        )
    lines.extend(["", "## Observed gene membership", ""])
    for provider in sorted({str(value) for value in genes["PROVIDER"]}) if not genes.empty else []:
        subset = genes.loc[genes["PROVIDER"].eq(provider)]
        symbols = sorted({str(value) for value in subset["GENE_SYMBOL"]})
        lines.append(f"- {provider} ({len(symbols)} genes, {len(subset)} term-gene rows): {', '.join(symbols) if symbols else 'none'}")
    lines.extend(["", f"## Open review rows: {len(discordance)}", ""])
    for _, row in discordance.head(20).iterrows():
        lines.append(f"- {row['MATCH_STATUS']} {row['ONTOLOGY']} {row['TERM_ID']}: {row['UNRESOLVED_MATCH_REASON']}")
    lines.extend(["", "## Output links", ""])
    for name in _OUTPUT_FILES:
        lines.append(f"- [{name}](./{name})")
    lines.extend([
        "", "## Limitations", "",
        "- Provider statistics remain separate; no pooled hypothesis test was calculated.",
        "- Authoritative public definitions were not retrieved in this offline run, so term interpretation stays provisional.",
        "- File hashes for every generated artifact are recorded in `postanalysis_provenance.json`, which is excluded from its own hash map by convention.",
        "- Null Jaccard denotes unavailable complete contributor sets, not disjoint gene membership.",
        "- Gene-set membership is observed annotation context, not pathogenicity evidence.",
        "- A nonspecific term for an MVA mechanism is context; modification terms annotate membership, not measured activity.",
    ])
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _write_csv(frame: pd.DataFrame, path: Path) -> None:
    frame.to_csv(path, index=False)


def analyze_pathway_results(
    baseline_dir: Path,
    david_dir: Path,
    output_dir: Path,
    *,
    analysis_threshold: float = 0.05,
) -> dict[str, Any]:
    """Verify saved artifacts and write an offline, auditable post-analysis.

    Tables are generated in an invocation-owned staging directory beside the
    destination and published by a single rename only after every artifact and its
    hashed provenance are present.  A failure leaves the destination untouched.
    """
    david_threshold = _analysis_threshold(analysis_threshold)
    baseline_dir, david_dir = baseline_dir.resolve(), david_dir.resolve()
    output = _validate_output(output_dir, baseline_dir, david_dir)
    p0, p1, baseline = _read_baseline(baseline_dir)
    receipt, chart, clusters, import_dir, input_hashes = _read_import(david_dir, baseline)
    evidence, genes = _build_evidence(
        p1, chart, baseline["_mapping"], set(baseline["_foreground"]), david_threshold
    )
    cluster_summary = _cluster_summary(clusters, evidence, baseline["_mapping"], set(baseline["_foreground"]))
    categories = receipt.get("interface", {}).get("selected_categories", [])
    discordance = _discordance(evidence, categories if isinstance(categories, list) else [])
    staging = output.parent / f".{output.name}.staging"
    if staging.exists():
        _fail("staging directory already exists")
    staging.mkdir()
    try:
        frames = {
            "integrated_pathway_evidence.csv": evidence,
            "gene_pathway_support.csv": genes,
            "david_cluster_summary.csv": cluster_summary,
            "provider_discordance_review.csv": discordance,
        }
        for name, frame in frames.items():
            _write_csv(frame, staging / name)
        provenance = {
            "schema_version": 1, "run_key": baseline["run_key"], "import_version": import_dir.name,
            "baseline_projections": baseline["baseline_projections"]["sha256"], "external_input_sha256": input_hashes,
            "analysis_threshold": david_threshold, "status": "completed",
            "output_sha256": {},
        }
        _write_report(staging / "integrated_pathway_report.md", receipt, evidence, cluster_summary, genes, discordance, provenance)
        content = [name for name in _OUTPUT_FILES if name != "postanalysis_provenance.json"]
        missing = [name for name in content if not (staging / name).is_file() or (staging / name).stat().st_size == 0]
        if missing:
            _fail(f"incomplete staged artifacts: {', '.join(missing)}")
        provenance["output_sha256"] = {name: _sha256(staging / name) for name in content}
        (staging / "postanalysis_provenance.json").write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        unpublished = [name for name in _OUTPUT_FILES if not (staging / name).is_file()]
        if unpublished:
            _fail(f"incomplete staged artifacts: {', '.join(unpublished)}")
        staging.rename(output)
    except Exception:
        shutil.rmtree(staging, ignore_errors=True)
        raise
    return {
        "status": "completed",
        "output_dir": str(output),
        "analysis_threshold": david_threshold,
        "evidence_terms": len(evidence),
        "clusters": len(cluster_summary),
    }
