"""Compass: saved enrichment evidence to auditable research decisions.

Compass is an evidence-to-decision layer, not another enrichment or network engine.  It
reads saved post-analysis tables and safe network projections, retrieves public membership,
protein-role and pharmacology records through one cached requester, and produces tables
that always say what is supported, what is merely context, and what is unresolved.

Two research lanes stay explicitly separate:

``PATHWAY_SUPPORTED``
    A provider's own corrected statistic clears its threshold for a pathway or process the
    target is an observed contributor to.  This is enrichment annotation, not causality.
``MECHANISM_PRIOR``
    The finite MVA / spindle-assembly-checkpoint / cohesin panel, researched from public
    gene-function evidence independent of whether any enrichment was significant.

``NETWORK_CONTEXT`` is saved network membership only: it never becomes independent causal
confirmation, and one shared source is never counted as replication.  No patient direction,
variant support, genotype or phenotype value enters this module, and no clinical
recommendation is produced.
"""

from __future__ import annotations

import csv
from datetime import UTC, datetime
import hashlib
import json
import os
from pathlib import Path
import shutil
from typing import Any, Iterable, Sequence

import pandas as pd

from mva_hackathon import compass_sources as sources
from mva_hackathon.drug_evidence import (
    assess_target_action,
    desired_change_for,
    disease_direction_for,
)
from mva_hackathon.mva_gene_sets import (
    CHROMOSOME_COHESIN_GENES,
    CORE_MVA_CAUSAL_GENES,
    SAC_CORE_GENES,
)

SCHEMA_VERSION = 1
SIGNIFICANCE_THRESHOLD = 0.05
_REPORT_SAMPLE = 20

LABEL_PATHWAY_SUPPORTED = "PATHWAY_SUPPORTED"
LABEL_MECHANISM_PRIOR = "MECHANISM_PRIOR"
LABEL_NETWORK_CONTEXT = "NETWORK_CONTEXT"

MECHANISM_PRIOR_PANEL = frozenset(CORE_MVA_CAUSAL_GENES | SAC_CORE_GENES | CHROMOSOME_COHESIN_GENES)

ENTITY_PATHWAY = "PATHWAY"
ENTITY_PROCESS = "BIOLOGICAL_PROCESS"
ENTITY_CONTEXT = "ANNOTATION_CONTEXT"

STOP_UNMAPPED = "STOP_UNMAPPED"
STOP_NOT_SIGNIFICANT = "STOP_NOT_SIGNIFICANT"
HOLD_STATISTICS = "HOLD_STATISTICS"
HOLD_MEMBERSHIP = "HOLD_MEMBERSHIP"
HOLD_RELEVANCE = "HOLD_RELEVANCE"
ELIGIBLE = "ELIGIBLE"

DECISION_REJECT = "REJECT"
DECISION_HOLD = "HOLD"
DECISION_TEST_FIRST = "TEST_FIRST"
DECISION_ADVANCE = "ADVANCE_RESEARCH"

PROVIDER_BY_ONTOLOGY = {
    "REACTOME": "Reactome",
    "KEGG": "KEGG",
    "WIKIPATHWAYS": "WikiPathways",
    "GO": "GO",
    "HP": "HPO",
}

PATHWAY_ONTOLOGIES = {"REACTOME", "KEGG", "WIKIPATHWAYS"}
PROCESS_ONTOLOGIES = {"GO"}

ELIGIBILITY_COLUMNS = (
    "EVIDENCE_ID", "PROVIDER", "ONTOLOGY", "TERM_ID", "TERM_NAME", "ENTITY_TYPE", "SIGNIFICANCE",
    "CONTRIBUTOR_STATUS", "MEMBERSHIP_STATUS", "DECISION", "REASON_CODES", "EVIDENCE_IDS",
)
MEMBERSHIP_COLUMNS = (
    "EVIDENCE_ID", "PROVIDER", "ONTOLOGY", "TERM_ID", "TARGET_SYMBOL", "MEMBERSHIP_KIND",
    "MEMBERSHIP_STATUS", "SOURCE_RECORD_ID",
)
PROTEIN_ROLE_COLUMNS = ("TARGET_SYMBOL", "UNIPROT_ACCESSION", "PROTEIN_ROLE", "ROLE_CONTEXT", "SOURCE_STATUS", "EVIDENCE_IDS")
MECHANISM_EVIDENCE_COLUMNS = (
    "EVIDENCE_ID", "TARGET_SYMBOL", "CLAIM_TYPE", "CLAIM", "SOURCE", "SOURCE_RECORD_ID", "CITATION_URL",
    "PMID", "DOI", "MODEL_SYSTEM", "SPECIES", "CELL_TYPE", "DEVELOPMENTAL_CONTEXT", "PERTURBATION",
    "DOSE_TEXT", "DURATION_TEXT", "OBSERVED_EFFECT", "SUPPORT_DIRECTION", "EVIDENCE_CLASS",
    "REVIEW_STATUS", "SOURCE_STATUS", "SOURCE_SECTION", "CONTEXT_LIMITATION",
)
RATIONALE_COLUMNS = (
    "RATIONALE_ID", "TARGET_SYMBOL", "TERM_ID", "HYPOTHESIS_SOURCE", "PATHWAY_ENRICHMENT_SUPPORT",
    "OBSERVED_MEMBERSHIP", "EXPANDED_MEMBERSHIP", "NETWORK_MEMBERSHIP", "PROTEIN_ROLE",
    "DISEASE_STATE_DIRECTION", "DESIRED_TARGET_CHANGE", "CONTEXT", "EVIDENCE_IDS", "COUNTEREVIDENCE_IDS",
    "DECISION", "REASON_CODES",
)
DRUG_REVIEW_COLUMNS = (
    "RATIONALE_ID", "TARGET_SYMBOL", "DRUG_ID", "DRUG_NAME", "SOURCE_RECORD_IDS", "ACTION_TYPE",
    "ACTION_DIRECTION", "DIRECTIONALITY_FIT", "HAZARD_STATUS", "EVIDENCE_IDS", "DECISION", "REASON_CODES",
)
LEDGER_COLUMNS = ("ENTITY_TYPE", "ENTITY_ID", "STAGE", "DECISION", "REASON_CODES", "EVIDENCE_IDS")
TEST_PLAN_COLUMNS = (
    "RATIONALE_ID", "HYPOTHESIS", "MODEL_SYSTEM", "PERTURBATION_DIRECTION", "PRIMARY_READOUT",
    "RESCUE_CRITERION", "WORSENING_CRITERION", "SAFETY_READOUTS", "CONTROLS", "EVIDENCE_IDS", "PLAN_STATUS",
)
SOURCE_STATUS_COLUMNS = ("REQUEST_ID", "STAGE", "SOURCE", "SCOPE", "STATUS", "COMPLETE", "RETRIEVED_AT", "RESPONSE_SHA256", "REASON")

REVIEWED_CLAIM_TYPES = {"PROTEIN_ROLE", "DISEASE_MECHANISM", "PERTURBATION_EFFECT", "RESCUE", "HAZARD", "SAFETY"}
REVIEWED_SUPPORT_DIRECTIONS = {
    "REDUCED_FUNCTION", "INCREASED_FUNCTION", "ALTERED_FUNCTION", "UNKNOWN", "SUPPORTS", "REFUTES", "WORSENS",
}
REVIEWED_EVIDENCE_CLASSES = {
    "HUMAN_GENETICS", "HUMAN_CELL_EXPERIMENT", "ANIMAL_EXPERIMENT", "CELL_EXPERIMENT", "BIOCHEMICAL_ASSAY",
    "CURATED_ANNOTATION", "COMPUTATIONAL", "HYPOTHESIS",
}

_COMPASS_TABLES = (
    "pathway_eligibility.csv",
    "expanded_members.csv",
    "protein_roles.csv",
    "mechanism_evidence.csv",
    "target_mechanism_rationale.csv",
    "target_drug_review.csv",
    "target_rejection_ledger.csv",
    "target_test_plan.csv",
    "source_status.csv",
    "compass_report.md",
)


def _fail(message: str, error: Exception | None = None) -> None:
    raise ValueError(message) from error


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _digest(value: Any) -> str:
    return hashlib.sha256(json.dumps(value, sort_keys=True, default=str).encode()).hexdigest()


def _text(value: Any) -> str:
    if value is None or value is pd.NA:
        return ""
    if isinstance(value, float) and pd.isna(value):
        return ""
    return str(value).strip()


def _flag(value: Any) -> bool:
    return bool(value) if not isinstance(value, str) else value.strip().lower() in {"true", "1", "yes"}


def _column(frame: pd.DataFrame, name: str) -> pd.Series:
    if name in frame.columns:
        return frame[name]
    return pd.Series([""] * len(frame), index=frame.index, dtype=object)


def _table(frame: pd.DataFrame, columns: Sequence[str]) -> pd.DataFrame:
    """Return a table with exactly the contracted columns, filling absent ones."""
    out = pd.DataFrame(index=frame.index)
    for column in columns:
        out[column] = _column(frame, column)
    return out.loc[:, list(columns)]


def _write_csv(frame: pd.DataFrame, path: Path) -> None:
    frame.to_csv(path, index=False, quoting=csv.QUOTE_MINIMAL)


# --------------------------------------------------------------------------- input loading


def _project_csv(path: Path, columns: Sequence[str]) -> pd.DataFrame:
    """Read only the projected columns of one saved artifact; extra columns never propagate."""
    if not path.is_file():
        _fail(f"missing required post-analysis artifact: {path.name}")
    try:
        frame = pd.read_csv(path, dtype=str, keep_default_na=False)
    except (OSError, UnicodeDecodeError, pd.errors.ParserError) as error:
        _fail(f"malformed post-analysis artifact: {path.name}", error)
    for column in columns:
        if column not in frame.columns:
            _fail(f"post-analysis artifact {path.name} is missing column {column}")
    return frame.loc[:, list(columns)]


def load_postanalysis(postanalysis_dir: Path) -> dict[str, Any]:
    """Load verified post-analysis tables plus their provenance identities."""
    root = Path(postanalysis_dir).resolve()
    if not root.is_dir():
        _fail("post-analysis directory does not exist")
    evidence = _project_csv(root / "integrated_pathway_evidence.csv", (
        "EVIDENCE_ID", "ONTOLOGY", "TERM_ID", "MATCH_STATUS", "GPROFILER_ROWS_JSON", "DAVID_ROWS_JSON",
        "GPROFILER_SIGNIFICANCE", "DAVID_SIGNIFICANCE", "EVIDENCE_TIER",
        "GPROFILER_CONTRIBUTOR_STATUS", "DAVID_CONTRIBUTOR_STATUS", "CONTRIBUTOR_STATUS", "MVA_CONTEXT_JSON",
    ))
    support = _project_csv(root / "gene_pathway_support.csv", (
        "EVIDENCE_ID", "ONTOLOGY", "TERM_ID", "PROVIDER", "GENE_SYMBOL", "MEMBERSHIP_STATUS",
    ))
    provenance_path = root / "postanalysis_provenance.json"
    if not provenance_path.is_file():
        _fail("post-analysis provenance is missing")
    try:
        provenance = json.loads(provenance_path.read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as error:
        _fail("post-analysis provenance is unreadable", error)
    return {
        "dir": root,
        "evidence": evidence,
        "support": support,
        "provenance": provenance,
        "projection_hash": _digest({
            "evidence": _digest(evidence.to_dict(orient="records")),
            "support": _digest(support.to_dict(orient="records")),
            "baseline": provenance.get("baseline_projections"),
            "import_version": provenance.get("import_version"),
        }),
    }


def load_network_projections(network_dir: Path | None) -> dict[str, Any]:
    """Read only the permitted gene/membership columns of saved Phase 2 artifacts.

    Saved network output predates any strong provenance contract, so everything derived
    here is labelled ``LEGACY_UNVERIFIED`` and is never treated as independent replication.
    """
    if network_dir is None:
        return {"status": "NETWORK_UNAVAILABLE", "memberships": pd.DataFrame(), "projection_hash": "", "reason": "no network directory supplied"}
    root = Path(network_dir).resolve()
    if not root.is_dir():
        return {"status": "NETWORK_UNAVAILABLE", "memberships": pd.DataFrame(), "projection_hash": "", "reason": "network directory does not exist"}
    edges_path = root / "04_phase2_string_edges.csv"
    context_path = root / "07_phase2_network_gene_context.csv"
    if not edges_path.is_file() and not context_path.is_file():
        return {"status": "NETWORK_UNAVAILABLE", "memberships": pd.DataFrame(), "projection_hash": "", "reason": "saved network tables are absent"}
    edges = _project_csv(edges_path, ("PREFERRED_NAME_A", "PREFERRED_NAME_B")) if edges_path.is_file() else pd.DataFrame(columns=["PREFERRED_NAME_A", "PREFERRED_NAME_B"])
    context = _project_csv(context_path, ("HGNC_SYMBOL",)) if context_path.is_file() else pd.DataFrame(columns=["HGNC_SYMBOL"])
    component_path = root / "10_phase2_component_summary.csv"
    component_gene_lists: list[list[str]] = []
    if component_path.is_file():
        try:
            saved = pd.read_csv(component_path, dtype=str, keep_default_na=False)
        except (OSError, UnicodeDecodeError, pd.errors.ParserError):
            saved = pd.DataFrame()
        if "GENE_LIST" in saved.columns:
            component_gene_lists = [
                [gene.strip().upper() for gene in _text(value).split(";") if gene.strip()]
                for value in saved["GENE_LIST"]
            ]
    components = _reconstruct_components(edges, context, component_gene_lists)
    records: list[dict[str, Any]] = []
    for component_id, members in enumerate(components, start=1):
        for symbol in sorted(members):
            records.append({
                "TARGET_SYMBOL": symbol,
                "MEMBERSHIP_KIND": "NETWORK_COMPONENT",
                "COMPONENT_ID": str(component_id),
                "MEMBERSHIP_STATUS": "LEGACY_UNVERIFIED",
                "SOURCE_RECORD_ID": f"string-component:{component_id}",
            })
    projection_hash = _digest({
        "edges": _digest(edges.to_dict(orient="records")),
        "context": _digest(context.to_dict(orient="records")),
        "components": _digest(component_gene_lists),
    })
    sizes = {}
    for path in (edges_path, context_path, component_path):
        if path.is_file():
            stat = path.stat()
            sizes[path.name] = {"size": stat.st_size, "mtime": int(stat.st_mtime)}
    return {
        "status": "NETWORK_CONTEXT_AVAILABLE" if records else "NETWORK_UNAVAILABLE",
        "memberships": pd.DataFrame(records, columns=["TARGET_SYMBOL", "MEMBERSHIP_KIND", "COMPONENT_ID", "MEMBERSHIP_STATUS", "SOURCE_RECORD_ID"]),
        "projection_hash": projection_hash,
        "files": sizes,
        "reason": "" if records else "no projected network membership was reconstructed",
    }


def _reconstruct_components(
    edges: pd.DataFrame,
    context: pd.DataFrame,
    saved_components: list[list[str]],
) -> list[set[str]]:
    """Reconstruct connected components from projected edges and every projected gene node.

    Component identifiers here are reconstruction order, not the saved visualization rank.
    """
    adjacency: dict[str, set[str]] = {}
    for _, row in edges.iterrows():
        left, right = _text(row.get("PREFERRED_NAME_A")).upper(), _text(row.get("PREFERRED_NAME_B")).upper()
        if left and right and left != right:
            adjacency.setdefault(left, set()).add(right)
            adjacency.setdefault(right, set()).add(left)
    nodes = set(adjacency)
    nodes.update(_text(value).upper() for value in context.get("HGNC_SYMBOL", pd.Series(dtype=object)) if _text(value))
    for group in saved_components:
        nodes.update(group)
    components: list[set[str]] = []
    visited: set[str] = set()
    for node in sorted(nodes):
        if node in visited:
            continue
        stack = [node]
        visited.add(node)
        members: set[str] = set()
        while stack:
            current = stack.pop()
            members.add(current)
            for neighbour in sorted(adjacency.get(current, set())):
                if neighbour not in visited:
                    visited.add(neighbour)
                    stack.append(neighbour)
        components.append(members)
    components.sort(key=lambda group: (len(group), sorted(group)[0]), reverse=True)
    return components


# --------------------------------------------------------------------------- eligibility


def _entity_type(ontology: str) -> str:
    if ontology in PATHWAY_ONTOLOGIES:
        return ENTITY_PATHWAY
    if ontology in PROCESS_ONTOLOGIES:
        return ENTITY_PROCESS
    return ENTITY_CONTEXT


def _provider_rows(record: pd.Series, provider: str) -> int:
    column = "GPROFILER_ROWS_JSON" if provider == "g:Profiler" else "DAVID_ROWS_JSON"
    try:
        value = json.loads(_text(record.get(column)) or "[]")
    except json.JSONDecodeError:
        return 0
    return len(value) if isinstance(value, list) else 0


def build_eligibility(evidence: pd.DataFrame, support: pd.DataFrame) -> pd.DataFrame:
    """One row per saved annotation per observing provider, each judged on its own statistics."""
    observed_genes: dict[tuple[str, str], set[str]] = {}
    for _, row in support.iterrows():
        term_id = _text(row.get("TERM_ID"))
        if term_id:
            observed_genes.setdefault((_text(row.get("ONTOLOGY")), term_id), set()).add(_text(row.get("GENE_SYMBOL")).upper())
    records: list[dict[str, Any]] = []
    for _, row in evidence.iterrows():
        evidence_id = _text(row.get("EVIDENCE_ID"))
        ontology = _text(row.get("ONTOLOGY")).upper()
        term_id = _text(row.get("TERM_ID"))
        entity_type = _entity_type(ontology)
        gene_names = sorted(observed_genes.get((ontology, term_id), set()))
        for provider, significance_column, status_column in (
            ("g:Profiler", "GPROFILER_SIGNIFICANCE", "GPROFILER_CONTRIBUTOR_STATUS"),
            ("DAVID", "DAVID_SIGNIFICANCE", "DAVID_CONTRIBUTOR_STATUS"),
        ):
            if _provider_rows(row, provider) == 0:
                # A provider that produced no row for this term contributes nothing to the gate.
                continue
            significance = _text(row.get(significance_column)) or "unknown"
            contributor_status = _text(row.get(status_column)) or "unavailable"
            membership_status = "observed" if gene_names else "unavailable"
            reasons: list[str] = []
            if not term_id or _text(row.get("MATCH_STATUS")) == "unmapped_term":
                decision = STOP_UNMAPPED
                reasons.append("NO_SUPPORTED_STABLE_IDENTIFIER")
            elif significance == "unknown":
                decision = HOLD_STATISTICS
                reasons.append("CORRECTED_SIGNIFICANCE_UNKNOWN")
            elif significance == "mixed":
                decision = HOLD_STATISTICS
                reasons.append("REPEATED_CORRECTED_VALUES_DISAGREE")
            elif significance == "not_comparable":
                decision = STOP_NOT_SIGNIFICANT
                reasons.append("PROVIDER_PRODUCED_NO_ROW_FOR_THIS_TERM")
            elif significance != "significant" or not _clears_threshold(row, provider):
                decision = STOP_NOT_SIGNIFICANT
                reasons.append("ABOVE_PROVIDER_THRESHOLD")
            elif entity_type == ENTITY_CONTEXT:
                decision = STOP_NOT_SIGNIFICANT
                reasons.append("ANNOTATION_CONTEXT_NOT_A_PATHWAY")
            elif contributor_status in {"unavailable", "absent"}:
                decision = HOLD_MEMBERSHIP
                reasons.append("CONTRIBUTOR_SET_UNAVAILABLE")
            elif not gene_names:
                decision = HOLD_MEMBERSHIP
                reasons.append("NO_OBSERVED_MEMBER_GENE")
            elif not (set(gene_names) & MECHANISM_PRIOR_PANEL):
                decision = HOLD_RELEVANCE
                reasons.append("NO_CITED_MVA_ROLE_TITLE_ALONE_IS_NOT_A_GATE")
            else:
                decision = ELIGIBLE
                reasons.append("PROVIDER_SIGNIFICANCE_AND_OBSERVED_PANEL_MEMBER")
            records.append({
                "EVIDENCE_ID": f"{evidence_id}::{provider.replace(':', '').replace(' ', '_').lower()}",
                "PROVIDER": provider,
                "ONTOLOGY": ontology,
                "TERM_ID": term_id,
                "TERM_NAME": _term_name(row, provider),
                "ENTITY_TYPE": entity_type,
                "SIGNIFICANCE": significance,
                "CONTRIBUTOR_STATUS": contributor_status,
                "MEMBERSHIP_STATUS": membership_status,
                "DECISION": decision,
                "REASON_CODES": ";".join(reasons),
                "EVIDENCE_IDS": evidence_id,
            })
    return pd.DataFrame(records, columns=ELIGIBILITY_COLUMNS)


def _clears_threshold(row: pd.Series, provider: str) -> bool:
    """Check the provider's own native corrected statistic, never a pooled value."""
    column = "GPROFILER_ROWS_JSON" if provider == "g:Profiler" else "DAVID_ROWS_JSON"
    value_column = "BH_Q" if provider == "g:Profiler" else "Benjamini"
    try:
        rows = json.loads(_text(row.get(column)) or "[]")
    except json.JSONDecodeError:
        return False
    if not isinstance(rows, list) or not rows:
        return False
    for entry in rows:
        if not isinstance(entry, dict):
            return False
        raw = entry.get(value_column)
        try:
            number = float(raw)
        except (TypeError, ValueError):
            return False
        if number > SIGNIFICANCE_THRESHOLD:
            return False
    return True


def _term_name(row: pd.Series, provider: str) -> str:
    column = "GPROFILER_ROWS_JSON" if provider == "g:Profiler" else "DAVID_ROWS_JSON"
    try:
        rows = json.loads(_text(row.get(column)) or "[]")
    except json.JSONDecodeError:
        return ""
    for entry in rows:
        if isinstance(entry, dict):
            name = _text(entry.get("TERM_NAME")) or _text(entry.get("Term")) or _text(entry.get("Full Term"))
            if name:
                return name
    return ""


# --------------------------------------------------------------------------- membership


def fetch_expanded_members(
    eligibility: pd.DataFrame,
    support: pd.DataFrame,
    cache_dir: Path,
    *,
    offline: bool = False,
    transport: sources.Transport | None = None,
) -> tuple[pd.DataFrame, list[dict[str, Any]], list[str]]:
    """Observed enrichment members plus complete public database membership for eligible terms."""
    records: list[dict[str, Any]] = []
    status_rows: list[dict[str, Any]] = []
    failures: list[str] = []
    eligible: dict[tuple[str, str], str] = {}
    for _, row in eligibility.iterrows():
        if _text(row.get("DECISION")) == ELIGIBLE:
            eligible[(_text(row.get("ONTOLOGY")), _text(row.get("TERM_ID")))] = _text(row.get("EVIDENCE_ID"))
    for _, row in support.iterrows():
        term_id = _text(row.get("TERM_ID"))
        if not term_id or (_text(row.get("ONTOLOGY")), term_id) not in eligible:
            continue
        records.append({
            "EVIDENCE_ID": _text(row.get("EVIDENCE_ID")),
            "PROVIDER": _text(row.get("PROVIDER")),
            "ONTOLOGY": _text(row.get("ONTOLOGY")),
            "TERM_ID": term_id,
            "TARGET_SYMBOL": _text(row.get("GENE_SYMBOL")).upper(),
            "MEMBERSHIP_KIND": "OBSERVED_INTERSECTION",
            "MEMBERSHIP_STATUS": _text(row.get("MEMBERSHIP_STATUS")) or "observed",
            "SOURCE_RECORD_ID": f"enrichment:{_text(row.get('EVIDENCE_ID'))}",
        })
    for (ontology, term_id), eligibility_id in sorted(eligible.items()):
        envelope = sources.fetch_membership(ontology, term_id, cache_dir, offline=offline, transport=transport)
        provider = PROVIDER_BY_ONTOLOGY.get(ontology, ontology)
        status_rows.append(_source_status_row(envelope, "membership", provider, term_id))
        if envelope["status"] in {sources.STATUS_UNAVAILABLE, sources.STATUS_INVALID_RESPONSE, sources.STATUS_OFFLINE_MISS}:
            failures.append(f"{provider}:{term_id}")
        membership_status = "complete" if envelope["complete"] else ("partial" if envelope["records"] else envelope["status"].lower())
        for record in envelope["records"]:
            symbol = _text(record.get("SYMBOL")).upper()
            if not symbol:
                continue
            records.append({
                "EVIDENCE_ID": eligibility_id,
                "PROVIDER": provider,
                "ONTOLOGY": ontology,
                "TERM_ID": term_id,
                "TARGET_SYMBOL": symbol,
                "MEMBERSHIP_KIND": "EXPANDED_DATABASE_MEMBER",
                "MEMBERSHIP_STATUS": membership_status,
                "SOURCE_RECORD_ID": _text(record.get("SOURCE_RECORD_ID")),
            })
    frame = pd.DataFrame(records, columns=MEMBERSHIP_COLUMNS)
    if not frame.empty:
        frame = frame.drop_duplicates(["TERM_ID", "TARGET_SYMBOL", "MEMBERSHIP_KIND"]).reset_index(drop=True)
    return frame, status_rows, failures


def _source_status_row(envelope: dict[str, Any], stage: str, source: str, scope: str) -> dict[str, Any]:
    """Public request identity and outcome — never a private query or patient value."""
    return {
        "REQUEST_ID": _text(envelope.get("request_id")),
        "STAGE": stage,
        "SOURCE": source,
        "SCOPE": scope,
        "STATUS": _text(envelope.get("status")),
        "COMPLETE": bool(envelope.get("complete")),
        "RETRIEVED_AT": _text(envelope.get("retrieved_at")),
        "RESPONSE_SHA256": _text(envelope.get("response_sha256")),
        "REASON": _text(envelope.get("reason")),
    }


# --------------------------------------------------------------------------- protein and literature


def canonicalize_targets(
    symbols: Iterable[str],
    cache_dir: Path,
    *,
    offline: bool = False,
    transport: sources.Transport | None = None,
) -> tuple[dict[str, str], dict[str, str], list[dict[str, Any]], list[str]]:
    """Resolve reported symbols to authoritative HGNC symbols.

    Returns canonical-by-reported mapping, the identity status per reported symbol, the
    source-status ledger rows, and any source failures.  Ambiguous or unresolvable
    symbols keep their reported identity and are separated rather than silently merged.
    """
    canonical: dict[str, str] = {}
    statuses: dict[str, str] = {}
    status_rows: list[dict[str, Any]] = []
    failures: list[str] = []
    for symbol in sorted({_text(item).upper() for item in symbols if _text(item)}):
        envelope = sources.fetch_gene_identity(symbol, cache_dir, offline=offline, transport=transport)
        statuses[symbol] = envelope["status"]
        status_rows.append(_source_status_row(envelope, "gene_identity", "HGNC", symbol))
        if envelope["status"] in {sources.STATUS_UNAVAILABLE, sources.STATUS_INVALID_RESPONSE, sources.STATUS_OFFLINE_MISS}:
            failures.append(f"HGNC:{symbol}")
        if envelope["status"] == sources.STATUS_OK and envelope["records"]:
            record = envelope["records"][0]
            canonical[symbol] = _text(record.get("CANONICAL_SYMBOL")) or symbol
            if _text(record.get("RESOLUTION")) != "exact_symbol":
                status_rows.append({
                    "REQUEST_ID": _text(envelope.get("request_id")),
                    "STAGE": "gene_identity_alias",
                    "SOURCE": "HGNC",
                    "SCOPE": f"{symbol}->{canonical[symbol]}",
                    "STATUS": sources.STATUS_OK,
                    "COMPLETE": True,
                    "RETRIEVED_AT": _text(envelope.get("retrieved_at")),
                    "RESPONSE_SHA256": _text(envelope.get("response_sha256")),
                    "REASON": _text(envelope.get("reason")),
                })
        else:
            canonical[symbol] = symbol
    return canonical, statuses, status_rows, failures


def fetch_protein_roles(
    symbols: Iterable[str],
    cache_dir: Path,
    *,
    offline: bool = False,
    transport: sources.Transport | None = None,
) -> tuple[pd.DataFrame, list[dict[str, Any]], dict[str, str], list[str]]:
    """Curated UniProt roles per candidate, with the source status kept explicit."""
    rows: list[dict[str, Any]] = []
    status_rows: list[dict[str, Any]] = []
    statuses: dict[str, str] = {}
    failures: list[str] = []
    for symbol in sorted({_text(item).upper() for item in symbols if _text(item)}):
        envelope = sources.fetch_protein_evidence(symbol, cache_dir, offline=offline, transport=transport)
        statuses[symbol] = envelope["status"]
        status_rows.append(_source_status_row(envelope, "protein_role", "UniProt", symbol))
        if envelope["status"] in {sources.STATUS_UNAVAILABLE, sources.STATUS_INVALID_RESPONSE, sources.STATUS_OFFLINE_MISS}:
            failures.append(f"UniProt:{symbol}")
        for record in envelope["records"]:
            accession = _text(record.get("UNIPROT_ACCESSION"))
            if _text(record.get("PROTEIN_ROLE")):
                rows.append(_protein_row(symbol, accession, "curated_function", record["PROTEIN_ROLE"], envelope["status"]))
            if _text(record.get("ROLE_CONTEXT")):
                rows.append(_protein_row(symbol, accession, "subunit_and_localization", record["ROLE_CONTEXT"], envelope["status"]))
            if _text(record.get("MODIFICATION_ANNOTATION")):
                rows.append(_protein_row(
                    symbol, accession, "modification_annotation",
                    "Modification annotation names annotated modifications; it is not measured activity. " + _text(record["MODIFICATION_ANNOTATION"]),
                    envelope["status"],
                ))
            for annotation in record.get("GO_ANNOTATIONS") or []:
                if not isinstance(annotation, dict):
                    continue
                rows.append(_protein_row(
                    symbol, accession, "go_annotation_membership",
                    f"{_text(annotation.get('GO_ID'))} ({_text(annotation.get('EVIDENCE_CODE'))}, {_text(annotation.get('EVIDENCE_CLASS'))})",
                    envelope["status"],
                ))
    return pd.DataFrame(rows, columns=PROTEIN_ROLE_COLUMNS), status_rows, statuses, failures


def _protein_row(symbol: str, accession: str, role: str, context: str, status: str) -> dict[str, Any]:
    return {
        "TARGET_SYMBOL": symbol,
        "UNIPROT_ACCESSION": accession,
        "PROTEIN_ROLE": role,
        "ROLE_CONTEXT": context,
        "SOURCE_STATUS": status,
        "EVIDENCE_IDS": f"ME-{_digest(['uniprot', accession, role])[:16]}",
    }


def fetch_literature_evidence(
    symbols: Iterable[str],
    cache_dir: Path,
    *,
    offline: bool = False,
    transport: sources.Transport | None = None,
) -> tuple[pd.DataFrame, list[dict[str, Any]], list[str]]:
    """Bounded literature shortlist per candidate, explicitly unreviewed."""
    rows: list[dict[str, Any]] = []
    status_rows: list[dict[str, Any]] = []
    failures: list[str] = []
    for symbol in sorted({_text(item).upper() for item in symbols if _text(item)}):
        envelope = sources.fetch_literature(symbol, cache_dir, offline=offline, transport=transport)
        status_rows.append(_source_status_row(envelope, "literature", "PubMed", symbol))
        if envelope["status"] in {sources.STATUS_UNAVAILABLE, sources.STATUS_INVALID_RESPONSE, sources.STATUS_OFFLINE_MISS}:
            failures.append(f"PubMed:{symbol}")
        for record in envelope["records"]:
            pmid = _text(record.get("PMID"))
            rows.append({
                "EVIDENCE_ID": f"ME-{_digest(['pubmed', pmid, symbol])[:16]}",
                "TARGET_SYMBOL": symbol,
                "CLAIM_TYPE": "LITERATURE_SEARCH_HIT",
                "CLAIM": _text(record.get("TITLE")),
                "SOURCE": "PubMed",
                "SOURCE_RECORD_ID": _text(record.get("SOURCE_RECORD_ID")),
                "CITATION_URL": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/" if pmid else "",
                "PMID": pmid,
                "DOI": _text(record.get("DOI")),
                "MODEL_SYSTEM": "",
                "SPECIES": "",
                "CELL_TYPE": "",
                "DEVELOPMENTAL_CONTEXT": "",
                "PERTURBATION": "",
                "DOSE_TEXT": "",
                "DURATION_TEXT": "",
                "OBSERVED_EFFECT": "",
                "SUPPORT_DIRECTION": "UNKNOWN",
                "EVIDENCE_CLASS": "HYPOTHESIS",
                "REVIEW_STATUS": "UNREVIEWED_SEARCH_HIT",
                "SOURCE_STATUS": envelope["status"],
                "SOURCE_SECTION": "title",
                "CONTEXT_LIMITATION": "A search hit is not a claim; mechanism review must cite the source section.",
            })
    return pd.DataFrame(rows, columns=MECHANISM_EVIDENCE_COLUMNS), status_rows, failures


def load_reviewed_evidence(path: Path | None) -> tuple[pd.DataFrame, str]:
    """Load and schema-validate the reviewed evidence CSV, returning its content hash."""
    if path is None:
        return pd.DataFrame(columns=MECHANISM_EVIDENCE_COLUMNS), ""
    reviewed_path = Path(path).resolve()
    if not reviewed_path.is_file():
        _fail("reviewed evidence file does not exist")
    try:
        frame = pd.read_csv(reviewed_path, dtype=str, keep_default_na=False)
    except (OSError, UnicodeDecodeError, pd.errors.ParserError) as error:
        _fail("reviewed evidence file is unreadable", error)
    missing = [column for column in MECHANISM_EVIDENCE_COLUMNS if column not in frame.columns]
    if missing:
        _fail(f"reviewed evidence is missing columns: {', '.join(missing)}")
    frame = frame.loc[:, list(MECHANISM_EVIDENCE_COLUMNS)].fillna("")
    for _, row in frame.iterrows():
        claim_type = _text(row.get("CLAIM_TYPE")).upper()
        direction = _text(row.get("SUPPORT_DIRECTION")).upper()
        evidence_class = _text(row.get("EVIDENCE_CLASS")).upper()
        if claim_type not in REVIEWED_CLAIM_TYPES:
            _fail(f"reviewed evidence has an unsupported CLAIM_TYPE: {claim_type}")
        if direction not in REVIEWED_SUPPORT_DIRECTIONS:
            _fail(f"reviewed evidence has an unsupported SUPPORT_DIRECTION: {direction}")
        if evidence_class not in REVIEWED_EVIDENCE_CLASSES:
            _fail(f"reviewed evidence has an unsupported EVIDENCE_CLASS: {evidence_class}")
        if _text(row.get("REVIEW_STATUS")).upper() != "REVIEWED":
            _fail("reviewed evidence rows must carry REVIEW_STATUS=REVIEWED")
        if not _text(row.get("SOURCE_SECTION")):
            _fail("reviewed evidence rows must record the cited SOURCE_SECTION")
        if not _text(row.get("CLAIM")):
            _fail("reviewed evidence rows must record a faithful CLAIM")
        if claim_type != "SAFETY" and not (_text(row.get("CITATION_URL")) or _text(row.get("PMID")) or _text(row.get("DOI"))):
            _fail("reviewed evidence rows must carry a resolvable citation")
    return frame, hashlib.sha256(reviewed_path.read_bytes()).hexdigest()


def fetch_drug_review(
    rationales: pd.DataFrame,
    cache_dir: Path,
    *,
    offline: bool = False,
    transport: sources.Transport | None = None,
) -> tuple[pd.DataFrame, list[dict[str, Any]], list[str]]:
    """Query pharmacology only for rationales whose desired direction is known."""
    rows: list[dict[str, Any]] = []
    status_rows: list[dict[str, Any]] = []
    failures: list[str] = []
    eligible = rationales.loc[rationales["DECISION"].isin([DECISION_TEST_FIRST, DECISION_ADVANCE])]
    for _, rationale in eligible.iterrows():
        symbol = _text(rationale["TARGET_SYMBOL"])
        envelope = sources.fetch_drug_evidence(symbol, cache_dir, offline=offline, transport=transport)
        status_rows.append(_source_status_row(envelope, "pharmacology", "ChEMBL+Open Targets", symbol))
        if envelope["status"] in {sources.STATUS_UNAVAILABLE, sources.STATUS_INVALID_RESPONSE, sources.STATUS_OFFLINE_MISS}:
            failures.append(f"Pharmacology:{symbol}")
        hazard = _text(rationale["REASON_CODES"]).find("CITED_CONTEXTUAL_HAZARD") >= 0
        for record in envelope["records"]:
            assessment = assess_target_action(
                _text(rationale["DESIRED_TARGET_CHANGE"]),
                _text(record.get("ACTION_TYPE")) or None,
                disease_direction=_text(rationale["DISEASE_STATE_DIRECTION"]) or "UNKNOWN",
                hazard_present=True if hazard else None,
            )
            rows.append({
                "RATIONALE_ID": _text(rationale["RATIONALE_ID"]),
                "TARGET_SYMBOL": symbol,
                "DRUG_ID": _text(record.get("DRUG_ID")),
                "DRUG_NAME": _text(record.get("DRUG_NAME")) or _text(record.get("DRUG_ID")),
                "SOURCE_RECORD_IDS": _text(record.get("SOURCE_RECORD_ID")),
                "ACTION_TYPE": _text(record.get("ACTION_TYPE")) or "UNKNOWN",
                "ACTION_DIRECTION": assessment["action_direction"],
                "DIRECTIONALITY_FIT": assessment["directionality_fit"],
                "HAZARD_STATUS": assessment["hazard_status"],
                "EVIDENCE_IDS": f"ME-{_digest(['pharmacology', symbol, record.get('DRUG_ID'), record.get('SOURCE')])[:16]}",
                "DECISION": assessment["decision"],
                "REASON_CODES": ";".join(assessment["reason_codes"] + [f"ASSAY_STATUS={_text(record.get('ASSAY_STATUS')) or 'ASSAY_UNAVAILABLE'}"]),
            })
    return pd.DataFrame(rows, columns=DRUG_REVIEW_COLUMNS), status_rows, failures


# --------------------------------------------------------------------------- pure decisions


def _reviewed_for(evidence: pd.DataFrame, symbol: str) -> pd.DataFrame:
    if evidence.empty:
        return evidence
    mask = evidence["TARGET_SYMBOL"].str.upper().eq(symbol) & evidence["REVIEW_STATUS"].str.upper().eq("REVIEWED")
    return evidence.loc[mask]


def _disease_direction(claims: pd.DataFrame) -> tuple[str, bool]:
    values = {_text(value).upper() for value in claims.get("SUPPORT_DIRECTION", pd.Series(dtype=object))}
    directions = {value for value in values if value in {"REDUCED_FUNCTION", "INCREASED_FUNCTION", "ALTERED_FUNCTION"}}
    if len(directions) > 1:
        return "UNKNOWN", True
    if len(directions) == 1:
        return directions.pop(), False
    return "UNKNOWN", False


def build_rationales(
    eligibility: pd.DataFrame,
    memberships: pd.DataFrame,
    protein_roles: pd.DataFrame,
    evidence: pd.DataFrame,
    network_memberships: pd.DataFrame,
) -> pd.DataFrame:
    """Pure target rationales: what is supported, what is context, and what is unresolved."""
    eligible_terms = {
        (_text(row["ONTOLOGY"]), _text(row["TERM_ID"]))
        for _, row in eligibility.iterrows()
        if _text(row["DECISION"]) == ELIGIBLE
    }
    observed: dict[str, set[str]] = {}
    expanded: dict[str, set[str]] = {}
    term_for: dict[str, set[str]] = {}
    if not memberships.empty:
        for _, row in memberships.iterrows():
            symbol = _text(row["TARGET_SYMBOL"]).upper()
            key = (_text(row["ONTOLOGY"]), _text(row["TERM_ID"]))
            if key not in eligible_terms:
                continue
            term_for.setdefault(symbol, set()).add(_text(row["TERM_ID"]))
            if _text(row["MEMBERSHIP_KIND"]) == "OBSERVED_INTERSECTION":
                observed.setdefault(symbol, set()).add(_text(row["TERM_ID"]))
            else:
                expanded.setdefault(symbol, set()).add(_text(row["TERM_ID"]))
    network_symbols = {_text(value).upper() for value in network_memberships.get("TARGET_SYMBOL", pd.Series(dtype=object)) if _text(value)}
    role_status = {
        _text(row["TARGET_SYMBOL"]).upper(): _text(row["SOURCE_STATUS"])
        for _, row in protein_roles.iterrows()
        if _text(row["PROTEIN_ROLE"]) == "curated_function"
    }
    role_text: dict[str, str] = {}
    for _, row in protein_roles.iterrows():
        symbol = _text(row["TARGET_SYMBOL"]).upper()
        if _text(row["PROTEIN_ROLE"]) == "curated_function" and _text(row["ROLE_CONTEXT"]):
            role_text.setdefault(symbol, _text(row["ROLE_CONTEXT"])[:240])

    # Candidate targets: the finite mechanism panel, plus pathway-supported targets, plus any
    # network gene that a reviewed evidence row explicitly nominates.
    nominated = {
        _text(row["TARGET_SYMBOL"]).upper()
        for _, row in evidence.iterrows()
        if _text(row["REVIEW_STATUS"]).upper() == "REVIEWED" and _text(row["TARGET_SYMBOL"])
    }
    candidates = set(MECHANISM_PRIOR_PANEL) | set(observed) | set(expanded) | (network_symbols & nominated)

    records: list[dict[str, Any]] = []
    for symbol in sorted(candidates):
        claims = _reviewed_for(evidence, symbol)
        pathway_supported = bool(observed.get(symbol))
        mechanism_prior = symbol in MECHANISM_PRIOR_PANEL
        network_member = symbol in network_symbols
        labels = []
        if pathway_supported:
            labels.append(LABEL_PATHWAY_SUPPORTED)
        if mechanism_prior:
            labels.append(LABEL_MECHANISM_PRIOR)
        if network_member:
            labels.append(LABEL_NETWORK_CONTEXT)
        disease_direction, direction_conflict = _disease_direction(claims)
        hazards = claims.loc[claims["CLAIM_TYPE"].str.upper().isin(["HAZARD"]) | claims["SUPPORT_DIRECTION"].str.upper().eq("WORSENS")]
        supports = claims.loc[claims["SUPPORT_DIRECTION"].str.upper().eq("SUPPORTS")]
        rescues = supports.loc[supports["CLAIM_TYPE"].str.upper().eq("RESCUE")]
        perturbations = supports.loc[supports["CLAIM_TYPE"].str.upper().isin(["PERTURBATION_EFFECT", "RESCUE"])]
        safety = claims.loc[claims["CLAIM_TYPE"].str.upper().eq("SAFETY")]
        identity_status = role_status.get(symbol, "UNAVAILABLE")
        desired = "UNKNOWN"
        if not claims.empty and disease_direction in {"REDUCED_FUNCTION", "INCREASED_FUNCTION"}:
            desired = "RESTORE_NORMAL"
        reasons: list[str] = []
        if pathway_supported:
            reasons.append("PATHWAY_SUPPORTED_OBSERVED_CONTRIBUTOR")
        if mechanism_prior:
            reasons.append("MECHANISM_PRIOR_PANEL")
        if not pathway_supported:
            reasons.append("NO_PATHWAY_ENRICHMENT_SUPPORT")
        if network_member and not pathway_supported and not mechanism_prior:
            reasons.append("NETWORK_CONTEXT_ONLY")
        model = next((_text(value) for value in perturbations.get("MODEL_SYSTEM", pd.Series(dtype=object)) if _text(value)), "")
        if not hazards.empty:
            decision = DECISION_REJECT
            reasons.append("CITED_CONTEXTUAL_HAZARD")
        elif identity_status in {"AMBIGUOUS", "EMPTY"}:
            decision = DECISION_HOLD
            reasons.append("IDENTITY_UNRESOLVED")
        elif identity_status in {"UNAVAILABLE", "INVALID_RESPONSE", "OFFLINE_MISS"}:
            decision = DECISION_HOLD
            reasons.append("SOURCE_INCOMPLETE")
        elif direction_conflict:
            decision = DECISION_HOLD
            reasons.append("EVIDENCE_CONFLICT")
        elif claims.empty:
            decision = DECISION_HOLD
            reasons.append("AWAITING_ROLE_EVIDENCE")
        elif desired == "UNKNOWN":
            decision = DECISION_HOLD
            reasons.append("DIRECTION_UNKNOWN")
        elif not rescues.empty and not (model and not safety.empty):
            # A rescue claim without a specified model or recorded safety readouts cannot advance.
            decision = DECISION_HOLD
            reasons.append("MODEL_UNSPECIFIED" if not model else "SAFETY_UNRECORDED")
        elif not rescues.empty:
            decision = DECISION_ADVANCE
            reasons.append("REVIEWED_RESCUE_WITH_MODEL_AND_SAFETY_READOUT")
        elif not perturbations.empty and model:
            decision = DECISION_TEST_FIRST
            reasons.append("CONTEXTUAL_HYPOTHESIS_WITH_PERTURBATION")
        elif not model:
            decision = DECISION_HOLD
            reasons.append("MODEL_UNSPECIFIED")
        else:
            decision = DECISION_HOLD
            reasons.append("SAFETY_UNRECORDED")
        records.append({
            "RATIONALE_ID": f"RAT-{_digest(['rationale', symbol])[:16]}",
            "TARGET_SYMBOL": symbol,
            "TERM_ID": ";".join(sorted(term_for.get(symbol, set()))),
            "HYPOTHESIS_SOURCE": "|".join(labels),
            "PATHWAY_ENRICHMENT_SUPPORT": LABEL_PATHWAY_SUPPORTED if pathway_supported else "NOT_IN_OBSERVED_INTERSECTION",
            "OBSERVED_MEMBERSHIP": pathway_supported,
            "EXPANDED_MEMBERSHIP": bool(expanded.get(symbol)),
            "NETWORK_MEMBERSHIP": network_member,
            "PROTEIN_ROLE": role_text.get(symbol, ""),
            "DISEASE_STATE_DIRECTION": disease_direction,
            "DESIRED_TARGET_CHANGE": desired,
            "CONTEXT": next((_text(value) for value in claims.get("MODEL_SYSTEM", pd.Series(dtype=object)) if _text(value)), "") or "unreported",
            "EVIDENCE_IDS": ";".join(sorted(_text(value) for value in claims.get("EVIDENCE_ID", pd.Series(dtype=object)) if _text(value))),
            "COUNTEREVIDENCE_IDS": ";".join(sorted(_text(value) for value in hazards.get("EVIDENCE_ID", pd.Series(dtype=object)) if _text(value))),
            "DECISION": decision,
            "REASON_CODES": ";".join(reasons),
        })
    return pd.DataFrame(records, columns=RATIONALE_COLUMNS)


def build_test_plans(rationales: pd.DataFrame, evidence: pd.DataFrame) -> pd.DataFrame:
    """Test plans, or an explicit evidence-resolution plan when the rationale is unresolved."""
    records: list[dict[str, Any]] = []
    for _, rationale in rationales.iterrows():
        symbol = _text(rationale["TARGET_SYMBOL"])
        claims = _reviewed_for(evidence, symbol)
        model = next((_text(value) for value in claims.get("MODEL_SYSTEM", pd.Series(dtype=object)) if _text(value)), "")
        evidence_ids = _text(rationale["EVIDENCE_IDS"])
        decision = _text(rationale["DECISION"])
        desired = _text(rationale["DESIRED_TARGET_CHANGE"])
        if decision in {DECISION_HOLD, DECISION_REJECT}:
            records.append({
                "RATIONALE_ID": _text(rationale["RATIONALE_ID"]),
                "HYPOTHESIS": f"Resolve {_text(rationale['REASON_CODES'])} for {symbol} before proposing any experiment",
                "MODEL_SYSTEM": model,
                "PERTURBATION_DIRECTION": "UNKNOWN",
                "PRIMARY_READOUT": "not applicable until the cited evidence is resolved",
                "RESCUE_CRITERION": "not applicable until the cited evidence is resolved",
                "WORSENING_CRITERION": "record any cited worsening signal as counterevidence",
                "SAFETY_READOUTS": "unresolved; safety must be established before any intervention is considered",
                "CONTROLS": "not applicable until the cited evidence is resolved",
                "EVIDENCE_IDS": evidence_ids,
                "PLAN_STATUS": "EVIDENCE_RESOLUTION_REQUIRED",
            })
            continue
        records.append({
            "RATIONALE_ID": _text(rationale["RATIONALE_ID"]),
            "HYPOTHESIS": (
                f"{desired} {symbol} function restores checkpoint or segregation competence in the cited model; "
                "this is a research hypothesis, not a treatment proposal"
            ),
            "MODEL_SYSTEM": model,
            "PERTURBATION_DIRECTION": f"genetic perturbation first: {desired}; concentrations require experimental calibration",
            "PRIMARY_READOUT": "mitotic timing and checkpoint-arrest competence; thresholds require experimental calibration",
            "RESCUE_CRITERION": "improved chromosome segregation fidelity without attributing selective cell death to correction",
            "WORSENING_CRITERION": "increased micronuclei, aneuploidy or premature mitotic exit",
            "SAFETY_READOUTS": "viability, DNA-damage response and proliferation in non-transformed cells; limits require experimental calibration",
            "CONTROLS": "vehicle and non-targeting controls, orthogonal target engagement, isogenic baseline where available",
            "EVIDENCE_IDS": evidence_ids,
            "PLAN_STATUS": "PROPOSED" if model else "REQUIRES_MODEL_SELECTION",
        })
    return pd.DataFrame(records, columns=TEST_PLAN_COLUMNS)


def build_rejection_ledger(
    eligibility: pd.DataFrame,
    rationales: pd.DataFrame,
    drug_review: pd.DataFrame,
    network_memberships: pd.DataFrame,
) -> pd.DataFrame:
    """Every held, unsupported and rejected case, so nothing is dropped silently."""
    records: list[dict[str, Any]] = []
    for _, row in eligibility.iterrows():
        if _text(row["DECISION"]) == ELIGIBLE:
            continue
        records.append({
            "ENTITY_TYPE": "TERM",
            "ENTITY_ID": f"{_text(row['ONTOLOGY'])}:{_text(row['TERM_ID']) or _text(row['EVIDENCE_ID'])}",
            "STAGE": "ELIGIBILITY",
            "DECISION": _text(row["DECISION"]),
            "REASON_CODES": _text(row["REASON_CODES"]),
            "EVIDENCE_IDS": _text(row["EVIDENCE_IDS"]),
        })
    for _, row in rationales.iterrows():
        if _text(row["DECISION"]) in {DECISION_TEST_FIRST, DECISION_ADVANCE}:
            continue
        records.append({
            "ENTITY_TYPE": "TARGET",
            "ENTITY_ID": _text(row["TARGET_SYMBOL"]),
            "STAGE": "RATIONALE",
            "DECISION": _text(row["DECISION"]),
            "REASON_CODES": _text(row["REASON_CODES"]),
            "EVIDENCE_IDS": _text(row["EVIDENCE_IDS"]),
        })
    for _, row in drug_review.iterrows():
        if _text(row["DECISION"]) in {DECISION_TEST_FIRST, DECISION_ADVANCE}:
            continue
        records.append({
            "ENTITY_TYPE": "DRUG",
            "ENTITY_ID": f"{_text(row['TARGET_SYMBOL'])}:{_text(row['DRUG_ID'])}",
            "STAGE": "PHARMACOLOGY",
            "DECISION": _text(row["DECISION"]),
            "REASON_CODES": _text(row["REASON_CODES"]),
            "EVIDENCE_IDS": _text(row["EVIDENCE_IDS"]),
        })
    supported = set(rationales.loc[rationales["NETWORK_MEMBERSHIP"].eq(True), "TARGET_SYMBOL"])
    for symbol in sorted({_text(value) for value in network_memberships.get("TARGET_SYMBOL", pd.Series(dtype=object)) if _text(value)} - supported):
        records.append({
            "ENTITY_TYPE": "NETWORK_GENE",
            "ENTITY_ID": symbol,
            "STAGE": "NETWORK",
            "DECISION": DECISION_HOLD,
            "REASON_CODES": "NETWORK_CONTEXT_ONLY",
            "EVIDENCE_IDS": "",
        })
    return pd.DataFrame(records, columns=LEDGER_COLUMNS)


# --------------------------------------------------------------------------- report


def _write_report(
    path: Path,
    eligibility: pd.DataFrame,
    memberships: pd.DataFrame,
    protein_roles: pd.DataFrame,
    evidence: pd.DataFrame,
    rationales: pd.DataFrame,
    drug_review: pd.DataFrame,
    ledger: pd.DataFrame,
    test_plans: pd.DataFrame,
    source_status: pd.DataFrame,
    network: dict[str, Any],
    provenance: dict[str, Any],
) -> None:
    approved = rationales.loc[rationales["DECISION"].eq(DECISION_ADVANCE)]
    test_first = rationales.loc[rationales["DECISION"].eq(DECISION_TEST_FIRST)]
    held = rationales.loc[rationales["DECISION"].eq(DECISION_HOLD)]
    rejected = rationales.loc[rationales["DECISION"].eq(DECISION_REJECT)]
    reviewed = evidence.loc[evidence["REVIEW_STATUS"].str.upper().eq("REVIEWED")]
    lines = [
        "# Compass report", "",
        "Research-exploration evidence integration. Nothing here is a diagnosis, a clinical",
        "recommendation, or experimental proof. Annotation is not efficacy; membership is not",
        "activity; an approved drug is not an approved treatment for this disease.", "",
        "## Lanes", "",
        "| Lane | Meaning | Targets |",
        "| --- | --- | --- |",
        f"| PATHWAY_SUPPORTED | provider-adjusted significance with an observed panel contributor | {int(rationales['HYPOTHESIS_SOURCE'].str.contains(LABEL_PATHWAY_SUPPORTED).sum())} |",
        f"| MECHANISM_PRIOR | finite MVA / SAC / cohesin panel researched from public gene function | {int(rationales['HYPOTHESIS_SOURCE'].str.contains(LABEL_MECHANISM_PRIOR).sum())} |",
        f"| NETWORK_CONTEXT | saved network membership only, never independent confirmation | {int(rationales['HYPOTHESIS_SOURCE'].str.contains(LABEL_NETWORK_CONTEXT).sum())} |",
        "",
        "A target may carry several labels; one shared source is never independent replication.", "",
        "## Pathway evidence", "",
        f"- saved annotations evaluated: {len(eligibility)} eligibility rows",
        f"- eligible pathways or processes: {int(eligibility['DECISION'].eq(ELIGIBLE).sum())}",
        f"- observed enrichment members: {int(memberships['MEMBERSHIP_KIND'].eq('OBSERVED_INTERSECTION').sum()) if not memberships.empty else 0}",
        f"- expanded public database members: {int(memberships['MEMBERSHIP_KIND'].eq('EXPANDED_DATABASE_MEMBER').sum()) if not memberships.empty else 0}",
        "- expanded membership never implies the gene was part of the original enrichment intersection",
        "",
        "### Held and stopped annotations",
        "",
    ]
    stopped = eligibility.loc[~eligibility["DECISION"].eq(ELIGIBLE)]
    if stopped.empty:
        lines.append("- none")
    else:
        for decision, count in stopped["DECISION"].value_counts().items():
            lines.append(f"- {decision}: {count} annotation(s)")
        for _, row in stopped.head(_REPORT_SAMPLE).iterrows():
            lines.append(f"- example: {row['PROVIDER']} {row['ONTOLOGY']} {row['TERM_ID'] or row['EVIDENCE_ID']}: {row['DECISION']} ({row['REASON_CODES']})")
        if len(stopped) > _REPORT_SAMPLE:
            lines.append(f"- {len(stopped) - _REPORT_SAMPLE} further stopped annotation(s) are listed in `pathway_eligibility.csv`")
    lines.extend(["", "## Mechanism-prior research", ""])
    lines.append(f"- curated protein-role rows: {len(protein_roles)}")
    lines.append(f"- reviewed mechanism assertions: {len(reviewed)}")
    lines.append(f"- unreviewed literature search hits retained for audit: {int(evidence['REVIEW_STATUS'].str.upper().eq('UNREVIEWED_SEARCH_HIT').sum())}")
    lines.append("- symbols resolved through the HGNC alias index: " + (
        "; ".join(
            f"{_text(row['SCOPE'])}"
            for _, row in source_status.iterrows()
            if _text(row["STAGE"]) == "gene_identity_alias"
        ) or "none"
    ))
    ambiguous = sorted(
        _text(row["SCOPE"])
        for _, row in source_status.iterrows()
        if _text(row["STAGE"]) == "gene_identity" and _text(row["STATUS"]) == sources.STATUS_AMBIGUOUS
    )
    lines.append(f"- identities left unresolved rather than chosen: {', '.join(ambiguous) or 'none'}")
    lines.extend(["", "## Decisions", "", f"- ADVANCE_RESEARCH (experimental prioritisation only): {len(approved)}", f"- TEST_FIRST: {len(test_first)}", f"- HOLD: {len(held)}", f"- REJECT: {len(rejected)}", ""])
    for _, row in rationales.head(_REPORT_SAMPLE).iterrows():
        lines.append(
            f"- {row['TARGET_SYMBOL']}: {row['DECISION']} | source={row['HYPOTHESIS_SOURCE'] or 'none'} | "
            f"direction={row['DISEASE_STATE_DIRECTION']} | desired={row['DESIRED_TARGET_CHANGE']} | reasons={row['REASON_CODES']}"
        )
    if len(rationales) > _REPORT_SAMPLE:
        lines.append(f"- {len(rationales) - _REPORT_SAMPLE} further rationale(s) are listed in `target_mechanism_rationale.csv`")
    lines.extend(["", "## Counterevidence and safety unknowns", ""])
    hazards = drug_review.loc[drug_review["DECISION"].eq(DECISION_REJECT)]
    if hazards.empty and rejected.empty:
        lines.append("- no cited worsening hazard was identified in the reviewed evidence")
    for _, row in rejected.iterrows():
        lines.append(f"- {row['TARGET_SYMBOL']}: cited hazard or opposite action ({row['COUNTEREVIDENCE_IDS'] or 'see reason codes'})")
    for _, row in hazards.iterrows():
        lines.append(f"- {row['TARGET_SYMBOL']} / {row['DRUG_NAME']}: {row['DIRECTIONALITY_FIT']} ({row['REASON_CODES']})")
    lines.append("- Unknown safety is recorded as unknown; it is never reported as safe.")
    lines.extend(["", "## Test plans", ""])
    for _, row in test_plans.iterrows():
        lines.append(f"- {row['RATIONALE_ID']} [{row['PLAN_STATUS']}] {row['HYPOTHESIS']}")
    lines.extend(["", "## Network context", ""])
    lines.append(f"- status: {network.get('status')} ({network.get('reason') or 'projected membership reconstructed'})")
    lines.append("- saved network projections are LEGACY_UNVERIFIED and are never independent causal confirmation")
    lines.extend(["", "## Source completeness", "", "| Stage | Source | Scope | Status | Complete |", "| --- | --- | --- | --- | --- |"])
    if source_status.empty:
        lines.append("| none | | | | |")
    for _, row in source_status.iterrows():
        lines.append(f"| {row['STAGE']} | {row['SOURCE']} | {row['SCOPE']} | {row['STATUS']} | {row['COMPLETE']} |")
    lines.extend(["", "## Rejection ledger summary", ""])
    if ledger.empty:
        lines.append("- empty: nothing was withheld")
    for decision, count in ledger["DECISION"].value_counts().items():
        lines.append(f"- {decision}: {count}")
    lines.extend(["", "## Output links", ""])
    for name in provenance["output_sha256"]:
        lines.append(f"- [{name}](./{name})")
    lines.append("- [compass_provenance.json](./compass_provenance.json)")
    lines.extend([
        "", "## Limitations", "",
        "- Provider statistics remain provider-specific; a single provider is judged on its own corrected values.",
        "- Term relevance is never established from a title; a held term needs a cited pathway role.",
        "- Public protein-role and modification annotation describe annotation, not measured activity.",
        "- A literature search hit is unreviewed until a reviewed assertion cites its source section.",
        "- Insufficient or failing source access is reported as such and never converted into a negative result.",
        "- No patient direction, variant support or clinical recommendation is produced by this layer.",
    ])
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


# --------------------------------------------------------------------------- orchestration


def _fixed_identity(postanalysis: dict[str, Any], network: dict[str, Any]) -> dict[str, Any]:
    """Identity that must match for a resumed run to write into the same destination."""
    return {
        "schema_version": SCHEMA_VERSION,
        "postanalysis_projection_hash": postanalysis["projection_hash"],
        "network_projection_hash": network.get("projection_hash", ""),
        "panel_hash": _digest(sorted(MECHANISM_PRIOR_PANEL)),
        "significance_threshold": SIGNIFICANCE_THRESHOLD,
        "selection_rules": "per-provider-significance+observed-panel-member",
    }


def _validate_output(output_dir: Path) -> Path:
    output = Path(output_dir).resolve()
    parent = output.parent
    if not parent.is_dir():
        _fail("output parent directory does not exist")
    return output


def _publish(staging: Path, output: Path, *, resume: bool) -> None:
    if not output.exists():
        staging.rename(output)
        return
    for name in sorted(os.listdir(staging)):
        os.replace(staging / name, output / name)
    shutil.rmtree(staging, ignore_errors=True)


def run_compass(
    postanalysis_dir: Path,
    network_dir: Path | None,
    output_dir: Path,
    *,
    cache_dir: Path,
    reviewed_evidence: Path | None = None,
    resume: bool = False,
    offline: bool = False,
    transport: sources.Transport | None = None,
) -> dict[str, Any]:
    """Consolidate saved evidence and public records into auditable research decisions."""
    output = _validate_output(output_dir)
    provenance_path = output / "compass_provenance.json"
    if output.exists() and not resume:
        _fail("output directory already exists; use --resume to continue an earlier run")
    postanalysis = load_postanalysis(postanalysis_dir)
    network = load_network_projections(network_dir)
    identity = _fixed_identity(postanalysis, network)
    previous: dict[str, Any] = {}
    if resume and provenance_path.is_file():
        try:
            previous = json.loads(provenance_path.read_text(encoding="utf-8"))
        except (OSError, UnicodeDecodeError, json.JSONDecodeError) as error:
            _fail("existing Compass provenance is unreadable", error)
        if previous.get("identity") != identity:
            _fail("resumed inputs differ from the recorded run; use a new output destination")

    reviewed, reviewed_hash = load_reviewed_evidence(reviewed_evidence)
    eligibility = build_eligibility(postanalysis["evidence"], postanalysis["support"])
    memberships, membership_status, membership_failures = fetch_expanded_members(
        eligibility, postanalysis["support"], Path(cache_dir), offline=offline, transport=transport
    )
    reported = set(MECHANISM_PRIOR_PANEL) | (set(memberships["TARGET_SYMBOL"]) if not memberships.empty else set())
    canonical, identity_status, identity_status_rows, identity_failures = canonicalize_targets(
        reported, Path(cache_dir), offline=offline, transport=transport
    )
    targets = {canonical.get(symbol, symbol) for symbol in reported}
    protein_roles, protein_status, protein_state, protein_failures = fetch_protein_roles(
        targets, Path(cache_dir), offline=offline, transport=transport
    )
    literature, literature_status, literature_failures = fetch_literature_evidence(
        targets, Path(cache_dir), offline=offline, transport=transport
    )
    evidence = pd.concat([literature, reviewed], ignore_index=True) if not reviewed.empty else literature
    evidence = _table(evidence, MECHANISM_EVIDENCE_COLUMNS)
    rationales = build_rationales(eligibility, memberships, protein_roles, evidence, network["memberships"])
    drug_review, drug_status, drug_failures = fetch_drug_review(rationales, Path(cache_dir), offline=offline, transport=transport)
    ledger = build_rejection_ledger(eligibility, rationales, drug_review, network["memberships"])
    test_plans = build_test_plans(rationales, evidence)
    source_status = pd.DataFrame(
        identity_status_rows + membership_status + protein_status + literature_status + drug_status,
        columns=SOURCE_STATUS_COLUMNS,
    )
    failures = identity_failures + membership_failures + protein_failures + literature_failures + drug_failures

    staging = output.parent / f".{output.name}.staging"
    if staging.exists():
        _fail("staging directory already exists")
    staging.mkdir()
    try:
        frames = {
            "pathway_eligibility.csv": eligibility,
            "expanded_members.csv": memberships,
            "protein_roles.csv": protein_roles,
            "mechanism_evidence.csv": evidence,
            "target_mechanism_rationale.csv": rationales,
            "target_drug_review.csv": drug_review,
            "target_rejection_ledger.csv": ledger,
            "target_test_plan.csv": test_plans,
            "source_status.csv": source_status,
        }
        for name, frame in frames.items():
            _write_csv(_table(frame, {
                "pathway_eligibility.csv": ELIGIBILITY_COLUMNS,
                "expanded_members.csv": MEMBERSHIP_COLUMNS,
                "protein_roles.csv": PROTEIN_ROLE_COLUMNS,
                "mechanism_evidence.csv": MECHANISM_EVIDENCE_COLUMNS,
                "target_mechanism_rationale.csv": RATIONALE_COLUMNS,
                "target_drug_review.csv": DRUG_REVIEW_COLUMNS,
                "target_rejection_ledger.csv": LEDGER_COLUMNS,
                "target_test_plan.csv": TEST_PLAN_COLUMNS,
                "source_status.csv": SOURCE_STATUS_COLUMNS,
            }[name]), staging / name)
        generated_at = datetime.now(UTC).isoformat()
        provenance = {
            "schema_version": SCHEMA_VERSION,
            "generated_at": generated_at,
            "resumed_from": previous.get("generated_at", "") if previous else "",
            "identity": identity,
            "reviewed_evidence_sha256": reviewed_hash,
            "postanalysis_dir": str(postanalysis["dir"]),
            "network_status": network.get("status", ""),
            "network_projection_hash": network.get("projection_hash", ""),
            "network_files": network.get("files", {}),
            "offline": offline,
            "source_failures": sorted(failures),
            "public_requests": sorted({_text(value) for value in source_status["REQUEST_ID"] if _text(value)}),
            "output_sha256": {},
        }
        _write_report(
            staging / "compass_report.md", eligibility, memberships, protein_roles, evidence,
            rationales, drug_review, ledger, test_plans, source_status, network, provenance,
        )
        published = [name for name in _COMPASS_TABLES if name != "compass_report.md"]
        provenance["output_sha256"] = {name: _sha256(staging / name) for name in published + ["compass_report.md"]}
        (staging / "compass_provenance.json").write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        missing = [name for name in _COMPASS_TABLES if not (staging / name).is_file()]
        if missing:
            _fail(f"incomplete Compass artifacts: {', '.join(missing)}")
        _publish(staging, output, resume=resume)
    except Exception:
        shutil.rmtree(staging, ignore_errors=True)
        raise
    return {
        "status": "source_incomplete" if failures else "completed",
        "output_dir": str(output),
        "eligibility_rows": len(eligibility),
        "eligible_terms": int(eligibility["DECISION"].eq(ELIGIBLE).sum()),
        "rationale_rows": len(rationales),
        "test_first": int(rationales["DECISION"].eq(DECISION_TEST_FIRST).sum()),
        "advance_research": int(rationales["DECISION"].eq(DECISION_ADVANCE).sum()),
        "hold": int(rationales["DECISION"].eq(DECISION_HOLD).sum()),
        "reject": int(rationales["DECISION"].eq(DECISION_REJECT).sum()),
        "drug_rows": len(drug_review),
        "source_failures": sorted(failures),
    }
