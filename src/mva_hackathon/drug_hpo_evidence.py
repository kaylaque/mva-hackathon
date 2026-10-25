"""Input validation and stable identities for drug–HPO evidence enrichment."""

from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor, as_completed

import hashlib
import json
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Sequence
import unicodedata
import pandas as pd
import shutil
from mva_hackathon.compass import DRUG_REVIEW_COLUMNS, RATIONALE_COLUMNS
from mva_hackathon.drug_evidence import assess_target_action, classify_action_direction


_HPO_CODE = re.compile(r"HP:\d{7}\Z")
_CHEMBL_ID = re.compile(r"CHEMBL\d+\Z")
_ROW_ID_FIELDS = ("RATIONALE_ID", "TARGET_SYMBOL", "DRUG_ID", "ACTION_TYPE", "SOURCE_RECORD_IDS")
_REQUIRED_EVIDENCE_LEVEL = "RELAXED_DAVID_BENJAMINI_0_90"


def sha256_bytes(value: bytes) -> str:
    """Return lowercase SHA-256 for exact input bytes."""
    return hashlib.sha256(value).hexdigest()


def sha256_file(path: Path) -> str:
    """Hash a file without changing its byte representation."""
    return sha256_bytes(path.read_bytes())


def load_hpo_codes(path: Path) -> tuple[list[str], str]:
    """Load reviewed code-only HPO input, preserving original-byte provenance."""
    path = Path(path)
    if not path.is_file():
        raise ValueError("HPO code file does not exist")
    content = path.read_bytes()
    try:
        text = content.decode("utf-8")
    except UnicodeDecodeError:
        raise ValueError("HPO code file contains an invalid line: 1") from None
    codes: list[str] = []
    seen: set[str] = set()
    for line_number, line in enumerate(text.splitlines(), start=1):
        if not line.strip():
            continue
        if not _HPO_CODE.fullmatch(line):
            raise ValueError(f"HPO code file contains an invalid line: {line_number}")
        if line not in seen:
            seen.add(line)
            codes.append(line)
    if not codes:
        raise ValueError("HPO code file contains no valid codes")
    return codes, sha256_bytes(content)


def assign_input_row_ids(drug_rows: pd.DataFrame) -> pd.DataFrame:
    """Copy drug rows and assign deterministic IDs unique to each input position."""
    missing = [field for field in _ROW_ID_FIELDS if field not in drug_rows.columns]
    if missing:
        raise ValueError(f"drug input missing required field: {missing[0]}")
    result = drug_rows.copy()
    row_numbers: list[int] = []
    row_ids: list[str] = []
    for row_number, (_, row) in enumerate(result.iterrows()):
        values = [row_number, *[str(row[field]) for field in _ROW_ID_FIELDS]]
        digest = hashlib.sha256(json.dumps(values, ensure_ascii=False, separators=(",", ":")).encode("utf-8")).hexdigest()
        row_numbers.append(row_number)
        row_ids.append(f"DHE-{digest[:16]}")
    result.insert(0, "INPUT_ROW_NUMBER", row_numbers)
    result.insert(1, "ROW_ID", row_ids)
    return result


def load_compass_snapshot(
    compass_dir: Path,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any], dict[str, str]]:
    """Load and validate the required non-sensitive Compass snapshot artifacts."""
    compass_dir = Path(compass_dir)
    artifact_names = (
        "target_drug_review.csv",
        "target_mechanism_rationale.csv",
        "compass_provenance.json",
    )
    paths = {name: compass_dir / name for name in artifact_names}
    for name, path in paths.items():
        if not path.is_file():
            raise ValueError(f"required Compass artifact is missing: {name}")

    try:
        drugs = pd.read_csv(paths["target_drug_review.csv"], dtype=str, keep_default_na=False)
    except (OSError, UnicodeError, pd.errors.ParserError):
        raise ValueError("invalid Compass artifact: target_drug_review.csv") from None
    try:
        rationales = pd.read_csv(paths["target_mechanism_rationale.csv"], dtype=str, keep_default_na=False)
    except (OSError, UnicodeError, pd.errors.ParserError):
        raise ValueError("invalid Compass artifact: target_mechanism_rationale.csv") from None
    for frame, columns, name in (
        (drugs, DRUG_REVIEW_COLUMNS, "target_drug_review.csv"),
        (rationales, RATIONALE_COLUMNS, "target_mechanism_rationale.csv"),
    ):
        missing = [column for column in columns if column not in frame.columns]
        if missing:
            raise ValueError(f"Compass artifact missing required field: {name} {missing[0]}")
    try:
        provenance = json.loads(paths["compass_provenance.json"].read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError):
        raise ValueError("invalid Compass artifact: compass_provenance.json") from None
    if not isinstance(provenance, dict) or not isinstance(provenance.get("identity"), dict):
        raise ValueError("Compass artifact missing required field: compass_provenance.json identity")

    rationale_ids = set(rationales["RATIONALE_ID"])
    if not drugs["RATIONALE_ID"].isin(rationale_ids).all():
        raise ValueError("Compass artifact contains unresolved RATIONALE_ID")
    if not drugs["PATHWAY_EVIDENCE_LEVEL"].eq(_REQUIRED_EVIDENCE_LEVEL).all():
        raise ValueError("Compass artifact has invalid PATHWAY_EVIDENCE_LEVEL")
    identity = provenance["identity"]
    if identity.get("selection_profile") != "david_relaxed":
        raise ValueError("Compass artifact has invalid selection_profile")
    try:
        threshold = float(identity.get("david_benjamini_threshold"))
    except (TypeError, ValueError):
        threshold = float("nan")
    if threshold != 0.9:
        raise ValueError("Compass artifact has invalid david_benjamini_threshold")
    if identity.get("drug_scope") != "relaxed_pathway":
        raise ValueError("Compass artifact has invalid drug_scope")
    selected = set(drugs["RATIONALE_ID"])
    selected_rationales = rationales.loc[rationales["RATIONALE_ID"].isin(selected)]
    if not selected_rationales["EXPLORATORY_DAVID_OBSERVED_TERM_IDS"].map(
        lambda value: "hsa05170" in {term.strip() for term in value.split(";")}
    ).all():
        raise ValueError("Compass artifact missing hsa05170 in EXPLORATORY_DAVID_OBSERVED_TERM_IDS")

    hashes = {name: sha256_file(path) for name, path in paths.items()}
    return drugs, rationales, provenance, hashes
CLAIM_COLUMNS = [
    "EVIDENCE_ID", "CANONICAL_DRUG_ID", "REGIMEN_ID", "TARGET_SYMBOL", "HPO_CODE",
    "CLAIM_TYPE", "CLAIM_STATUS", "OBSERVED_OR_INFERRED", "OUTCOME_DIRECTION",
    "APPLICABILITY", "ATTRIBUTION", "STUDY_DESIGN", "MODEL_OR_POPULATION",
    "INTERVENTION", "COMPARATOR", "DOSE_TEXT", "DURATION_TEXT", "SOURCE_TYPE",
    "SOURCE_RECORD_ID", "PMID", "DOI", "NCT_ID", "FDA_SET_ID", "PASSAGE_LOCATOR",
    "PASSAGE_SHA256", "CHECKED_AT", "LIMITATIONS", "REVIEW_STATUS",
]
CLAIM_STATUSES = {"PROPOSED", "QUALIFYING", "CONTEXT_ONLY", "CONTRADICTED", "REJECTED_IDENTITY", "RETRACTED"}
SEARCH_STATUSES = {"not_searched", "retrieval_failed", "screening_pending", "reviewed_no_qualifying_evidence", "evidence_coded"}
_CLAIM_ID_FIELDS = ("CANONICAL_DRUG_ID", "REGIMEN_ID", "TARGET_SYMBOL", "HPO_CODE", "CLAIM_TYPE", "SOURCE_TYPE", "SOURCE_RECORD_ID", "PMID", "DOI", "NCT_ID", "FDA_SET_ID", "PASSAGE_LOCATOR")
_BENEFIT = re.compile(r"\b(improved|resolved|prevented|stabilized|benefited|responded)\b", re.I)
_WORSEN = re.compile(r"\b(worsened|caused|induced|exacerbated|contraindicated)\b", re.I)
_DIRECTIONAL = re.compile(r"\b(increased|decreased|reduced)\b", re.I)
_NEGATION = re.compile(r"\b(no|not|never|without|failed to|didn't|doesn't|cannot|unable to)\b", re.I)
_SENTENCE = re.compile(r"(?<=[.!?])\s+")


def _clean(value: Any) -> str:
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return ""
    return " ".join(str(value).split()).strip()


def _terms(value: Any) -> list[str]:
    if value is None:
        return []
    if isinstance(value, (list, tuple, set)):
        values = list(value)
    else:
        values = re.split(r"[;|]", str(value))
    return list(dict.fromkeys(_clean(item) for item in values if _clean(item)))


def _normalized_name(value: str) -> str:
    return " ".join(unicodedata.normalize("NFKC", value).casefold().split())


def _combination_id(drug_id: str, component_ids: list[str], name: str) -> str:
    identity = json.dumps([drug_id, component_ids, _normalized_name(name)], ensure_ascii=False, separators=(",", ":"))
    return "DHR-" + hashlib.sha256(identity.encode("utf-8")).hexdigest()[:16]


def normalize_drug_inventory(drug_rows: pd.DataFrame) -> pd.DataFrame:
    """Add stable canonical identities without merging source drug rows."""
    result = drug_rows.copy()
    if "ROW_ID" not in result:
        result = assign_input_row_ids(result)
    for column in ("DRUG_ID", "DRUG_NAME"):
        if column not in result:
            result[column] = ""
    names_by_id = {_clean(row["DRUG_ID"]): _clean(row["DRUG_NAME"]) for _, row in result.iterrows() if _clean(row["DRUG_ID"]) and _clean(row["DRUG_NAME"])}
    ids, aliases, regimens, components, component_names = [], [], [], [], []
    for _, row in result.iterrows():
        drug_id, name = _clean(row["DRUG_ID"]), _clean(row["DRUG_NAME"])
        source_aliases = _terms(row.get("ALIASES", []))
        source_aliases.extend(_terms(row.get("SOURCE_ALIASES", [])))
        # Aliases are admitted only when present as exact source-record names.
        names = [name, *source_aliases]
        aliases.append(list(dict.fromkeys(value for value in names if value)))
        canonical = drug_id or _normalized_name(name)
        ids.append(canonical)
        component_ids = _terms(row.get("COMPONENT_DRUG_IDS", []))
        components.append(component_ids)
        component_names.append([names_by_id[component] for component in component_ids if component in names_by_id])
        regimen = bool(component_ids) or len(_terms(row.get("COMPONENTS", []))) > 1
        regimens.append(_combination_id(drug_id, component_ids, name) if regimen else "")
    result["CANONICAL_DRUG_ID"] = ids
    result["ALIASES"] = aliases
    result["REGIMEN_ID"] = regimens
    result["COMPONENT_DRUG_IDS"] = components
    result["COMPONENT_NAMES"] = component_names
    return result


def _query(drug_names: Sequence[str], phenotype_names: Sequence[str] = ()) -> str:
    quote = lambda item: '"' + item.replace('"', '\\"') + '"'
    left = " OR ".join(quote(item) for item in drug_names if item)
    if not phenotype_names:
        return left
    right = " OR ".join(quote(item) for item in phenotype_names if item)
    return f"({left}) AND ({right})"


def build_search_inventory(drugs: pd.DataFrame, hpo_terms: pd.DataFrame) -> pd.DataFrame:
    """Build stable source scopes for every canonical drug and HPO combination."""
    normalized = normalize_drug_inventory(drugs) if not {"CANONICAL_DRUG_ID", "ALIASES"}.issubset(drugs.columns) else drugs
    output = []
    term_rows = list(hpo_terms.to_dict("records"))
    for drug_id, group in normalized.groupby("CANONICAL_DRUG_ID", sort=True, dropna=False):
        aliases = list(dict.fromkeys(alias for group_aliases in group["ALIASES"] for alias in _terms(group_aliases)))
        name = _clean(group.iloc[0].get("DRUG_NAME", ""))
        names = list(dict.fromkeys([name, *aliases]))
        output.append(["literature", drug_id, "", _query(names), "not_searched"])
        output.append(["trials", drug_id, "", _query(names), "not_searched"])
        for term in sorted(term_rows, key=lambda item: _clean(item.get("HPO_CODE"))):
            phenotype = list(dict.fromkeys([_clean(term.get("HPO_LABEL")), *_terms(term.get("HPO_SYNONYMS"))]))
            code = _clean(term.get("HPO_CODE"))
            output.append(["literature", drug_id, code, _query(names, phenotype), "not_searched"])
            output.append(["trials", drug_id, code, _query(names, phenotype), "not_searched"])
    columns = ["SOURCE", "CANONICAL_DRUG_ID", "HPO_CODE", "QUERY", "STATUS"]
    return pd.DataFrame(output, columns=columns).sort_values(columns[:4], kind="stable").reset_index(drop=True)


def _contains(text: str, candidates: Sequence[str]) -> bool:
    folded = _normalized_name(text)
    return any(re.search(r"(?<!\w)" + re.escape(_normalized_name(item)) + r"(?!\w)", folded) for item in candidates if item)


def _split_sentences(text: str, section: str) -> list[tuple[str, str]]:
    return [(sentence.strip(), f"{section}:p{paragraph + 1}:s{index + 1}") for paragraph, block in enumerate(re.split(r"\n+", text)) if block.strip() for index, sentence in enumerate(_SENTENCE.split(block.strip())) if sentence.strip()]


def _record_id(record: dict[str, Any], source: str, index: int) -> str:
    return _clean(record.get("SOURCE_RECORD_ID") or record.get("SPL_ID") or record.get("NCT_ID") or record.get("PMID") or record.get("DOI") or f"{source.lower()}-{index + 1}")


def _base_claim(drug: pd.Series, record: dict[str, Any], source: str, record_id: str, passage: str, locator: str, claim_type: str, status: str, *, hpo_code: str = "", direction: str = "", applicability: str = "unknown", attribution: str = "INDIVIDUAL") -> dict[str, Any]:
    passage_hash = sha256_bytes(passage.encode("utf-8"))
    claim = {column: "" for column in CLAIM_COLUMNS}
    claim.update({
        "CANONICAL_DRUG_ID": _clean(drug.get("CANONICAL_DRUG_ID")),
        "REGIMEN_ID": _clean(drug.get("REGIMEN_ID")),
        "TARGET_SYMBOL": _clean(drug.get("TARGET_SYMBOL")),
        "HPO_CODE": hpo_code, "CLAIM_TYPE": claim_type, "CLAIM_STATUS": status,
        "OBSERVED_OR_INFERRED": "OBSERVED" if claim_type in {"BENEFIT_OBSERVED", "WORSENING_OBSERVED", "MECHANISM_OBSERVED"} else "CONTEXT",
        "OUTCOME_DIRECTION": direction, "APPLICABILITY": applicability,
        "ATTRIBUTION": attribution, "STUDY_DESIGN": _clean(record.get("STUDY_DESIGN") or record.get("STUDY_TYPE") or record.get("PUBLICATION_TYPES")),
        "MODEL_OR_POPULATION": _clean(record.get("MODEL_OR_POPULATION")),
        "INTERVENTION": _clean(record.get("INTERVENTION") or record.get("INTERVENTIONS")),
        "COMPARATOR": _clean(record.get("COMPARATOR")), "DOSE_TEXT": _clean(record.get("DOSE_TEXT")),
        "DURATION_TEXT": _clean(record.get("DURATION_TEXT")), "SOURCE_TYPE": source,
        "SOURCE_RECORD_ID": record_id, "PMID": _clean(record.get("PMID")),
        "DOI": _clean(record.get("DOI")), "NCT_ID": _clean(record.get("NCT_ID")),
        "FDA_SET_ID": _clean(record.get("FDA_SET_ID")), "PASSAGE_LOCATOR": locator,
        "PASSAGE_SHA256": passage_hash, "CHECKED_AT": _clean(record.get("CHECKED_AT")),
        "LIMITATIONS": "", "REVIEW_STATUS": "AUTOMATED",
    })
    if not claim["CHECKED_AT"]:
        claim["CHECKED_AT"] = ""
    _set_evidence_id(claim)
    return claim


def _set_evidence_id(claim: dict[str, Any]) -> None:
    identity = [claim[field] for field in _CLAIM_ID_FIELDS] + [claim["PASSAGE_SHA256"]]
    claim["EVIDENCE_ID"] = "DHEV-" + hashlib.sha256(json.dumps(identity, ensure_ascii=False, separators=(",", ":")).encode("utf-8")).hexdigest()[:16]


def _passages(record: dict[str, Any], source: str) -> list[tuple[str, str]]:
    fields = (("ABSTRACT", "ABSTRACT"), ("FULL_TEXT_PARAGRAPHS", "FULL_TEXT")) if source == "LITERATURE" else ()
    passages = []
    for field, section in fields:
        value = record.get(field, [])
        blocks = value if isinstance(value, list) else [value]
        for block in blocks:
            passages.extend(_split_sentences(_clean(block), section))
    if source == "FDA":
        for field, section, typ in (
            ("BOXED_WARNINGS", "BOXED_WARNING", "HAZARD"), ("CONTRAINDICATIONS", "CONTRAINDICATION", "HAZARD"),
            ("WARNINGS", "WARNING", "HAZARD"), ("ADVERSE_REACTIONS", "ADVERSE_REACTION", "SAFETY_CONTEXT"),
            ("PEDIATRIC_USE", "PEDIATRIC_USE", "SAFETY_CONTEXT"),
        ):
            value = record.get(field, [])
            blocks = value if isinstance(value, list) else [value]
            for block in blocks:
                for text, locator in _split_sentences(_clean(block), section):
                    passages.append((text, f"{locator}|{typ}"))
    if source == "TRIAL":
        # Trial registration is context, not observed benefit.
        for field, section in (("BRIEF_SUMMARY", "TRIAL_SUMMARY"), ("OUTCOME_MEASURES", "TRIAL_OUTCOME")):
            value = record.get(field, [])
            blocks = value if isinstance(value, list) else [value]
            for block in blocks:
                passages.extend(_split_sentences(_clean(block), section))
        metadata = f"Interventions: {_clean(record.get('INTERVENTIONS'))}; conditions: {_clean(record.get('CONDITIONS'))}"
        if metadata != "Interventions: ; conditions: ":
            passages.extend(_split_sentences(metadata, "TRIAL_METADATA"))
    return passages


def extract_claims(drugs: pd.DataFrame, hpo_terms: pd.DataFrame, literature_records: Sequence[dict[str, Any]], fda_records: Sequence[dict[str, Any]], trial_records: Sequence[dict[str, Any]]) -> pd.DataFrame:
    """Extract deterministic sentence-linked claims; passages are never published."""
    normalized = normalize_drug_inventory(drugs) if not {"CANONICAL_DRUG_ID", "ALIASES"}.issubset(drugs.columns) else drugs
    names_by_drug = {
        drug_id: list(dict.fromkeys(name for _, drug_row in group.iterrows() for name in _terms(drug_row.get("ALIASES", []))))
        for drug_id, group in normalized.groupby("CANONICAL_DRUG_ID", sort=False)
    }
    term_rows = hpo_terms.to_dict("records")
    records = [("LITERATURE", literature_records), ("FDA", fda_records), ("TRIAL", trial_records)]
    claims: list[dict[str, Any]] = []
    for source, source_records in records:
        for record_index, record in enumerate(source_records):
            rid = _record_id(record, source, record_index)
            retracted = bool(record.get("RETRACTED")) or "retract" in _clean(record.get("PUBLICATION_STATUS")).casefold()
            for _, drug in normalized.iterrows():
                if _clean(record.get("_SCOPE_CANONICAL_DRUG_ID")) not in {"", _clean(drug.get("CANONICAL_DRUG_ID"))}:
                    continue
                names = _terms(drug.get("ALIASES", [])) or [_clean(drug.get("DRUG_NAME"))]
                names = list(dict.fromkeys([_clean(drug.get("DRUG_NAME")), *names]))
                if source == "FDA":
                    exact_names = [*_terms(record.get("GENERIC_NAMES")), *_terms(record.get("BRAND_NAMES")), *_terms(record.get("SUBSTANCE_NAMES"))]
                    if not any(_normalized_name(n) == _normalized_name(x) for n in names for x in exact_names):
                        continue
                passages = _passages(record, source)
                if source == "FDA":
                    for passage, locator in passages:
                        claim_type = locator.rsplit("|", 1)[-1]
                        matched_terms = [
                            (_clean(term.get("HPO_CODE")), _clean(term.get("HPO_LABEL")), _terms(term.get("HPO_SYNONYMS")))
                            for term in term_rows
                            if _contains(passage, [_clean(term.get("HPO_LABEL")), *_terms(term.get("HPO_SYNONYMS"))])
                        ]
                        matches = matched_terms or [("", "", [])]
                        for code, _, _ in matches:
                            status = "QUALIFYING" if claim_type == "HAZARD" and code else "CONTEXT_ONLY"
                            claim = _base_claim(drug, record, source, rid, passage, locator, claim_type, status, hpo_code=code, direction="HAZARD" if claim_type == "HAZARD" else "", applicability="direct" if code else "unknown")
                            claims.append(claim)
                    continue
                if source == "TRIAL" and not passages:
                    passages = [(f"Interventions: {_clean(record.get('INTERVENTIONS'))}; conditions: {_clean(record.get('CONDITIONS'))}", "TRIAL_METADATA")]
                for term in term_rows:
                    code = _clean(term.get("HPO_CODE"))
                    if _clean(record.get("_SCOPE_HPO_CODE")) not in {"", code}:
                        continue
                    phenotypes = list(dict.fromkeys([_clean(term.get("HPO_LABEL")), *_terms(term.get("HPO_SYNONYMS"))]))
                    pairs: list[tuple[str, str, str]] = []
                    for index, (sentence, locator) in enumerate(passages):
                        pair = [(sentence, locator)]
                        if index + 1 < len(passages):
                            pair.append(passages[index + 1])
                        text = " ".join(s for s, _ in pair)
                        if (source != "FDA" and not _contains(text, names)) or not _contains(text, phenotypes):
                            continue
                        direction = ""
                        claim_type, status = "PHENOTYPE_ASSOCIATION", "CONTEXT_ONLY"
                        benefit = bool(_BENEFIT.search(text)) and not bool(_NEGATION.search(text))
                        worsening = bool(_WORSEN.search(text)) and not bool(_NEGATION.search(text))
                        if benefit and worsening:
                            claim_type, status, direction = "PHENOTYPE_ASSOCIATION", "CONTRADICTED", "CONFLICT"
                        elif benefit:
                            claim_type, status, direction = "BENEFIT_OBSERVED", "QUALIFYING", "BENEFIT"
                        elif worsening:
                            claim_type, status, direction = "WORSENING_OBSERVED", "QUALIFYING", "WORSENING"
                        elif _DIRECTIONAL.search(text):
                            claim_type, status, direction = "PHENOTYPE_ASSOCIATION", "PROPOSED", "AMBIGUOUS"
                        elif _NEGATION.search(text):
                            claim_type, status = "PHENOTYPE_ASSOCIATION", "CONTEXT_ONLY"
                        else:
                            continue
                        matched_components = [component for component in _terms(drug.get("COMPONENT_NAMES", [])) if _contains(text, [component])]
                        explicit_drugs = [drug_id for drug_id, drug_names in names_by_drug.items() if _contains(text, drug_names)]
                        attribution = "REGIMEN_ONLY" if len(explicit_drugs) > 1 or (drug.get("REGIMEN_ID") and len(matched_components) != 1) else "INDIVIDUAL"
                        if retracted:
                            status = "RETRACTED"
                        if source == "TRIAL":
                            claim_type, status, direction = "TRIAL_CONTEXT", "CONTEXT_ONLY", ""
                        pairs.append((text, "|".join(locator for _, locator in pair), f"{claim_type}\x1f{status}\x1f{direction}\x1f{attribution}"))
                    if pairs:
                        states = {p[2] for p in pairs}
                        if not retracted and (any("\x1fCONFLICT" in state for state in states) or (any("\x1fBENEFIT" in state for state in states) and any("\x1fWORSENING" in state for state in states))):
                            pairs = [(text, loc, f"PHENOTYPE_ASSOCIATION\x1fCONTRADICTED\x1fCONFLICT\x1f{state.split(chr(31))[-1]}") for text, loc, state in pairs]
                        for passage, locator, state in pairs:
                            claim_type, status, direction, attribution = state.split("\x1f")
                            claim = _base_claim(drug, record, source, rid, passage, locator, claim_type, status, hpo_code=code, direction=direction, applicability="direct", attribution=attribution)
                            if source == "FDA":
                                section_type = locator.rsplit("|", 1)[-1]
                                claim["CLAIM_TYPE"] = section_type
                                claim["CLAIM_STATUS"] = "CONTEXT_ONLY" if section_type == "SAFETY_CONTEXT" else status
                                claim["OUTCOME_DIRECTION"] = "HAZARD" if section_type == "HAZARD" else ""
                                _set_evidence_id(claim)
                            claims.append(claim)
                if source == "LITERATURE":
                    for sentence, locator in passages:
                        target = _clean(drug.get("TARGET_SYMBOL"))
                        action = _clean(drug.get("ACTION_TYPE"))
                        action_direction = classify_action_direction(action)
                        if target and _contains(sentence, names) and _contains(sentence, [target]) and action_direction != "UNKNOWN":
                            action_terms = {"INCREASE": r"\b(activat\w*|agoniz\w*|stimulat\w*|increas\w*)\b", "DECREASE": r"\b(inhibit\w*|antagoniz\w*|block\w*|suppress\w*|decreas\w*|reduce\w*)\b"}
                            explicit_action = re.search(action_terms[action_direction], sentence, re.I)
                            observed_effect = re.search(r"\b(reduced|increased|decreased|altered|changed|suppressed|enhanced)\b", sentence, re.I)
                            if explicit_action and observed_effect:
                                claim = _base_claim(drug, record, source, rid, sentence, locator, "MECHANISM_OBSERVED", "RETRACTED" if retracted else "QUALIFYING", direction=action_direction)
                                claims.append(claim)
    output = pd.DataFrame(claims, columns=CLAIM_COLUMNS)
    if output.empty:
        return output
    return output.drop_duplicates("EVIDENCE_ID", keep="first").sort_values(["CANONICAL_DRUG_ID", "HPO_CODE", "SOURCE_TYPE", "SOURCE_RECORD_ID", "PASSAGE_LOCATOR", "EVIDENCE_ID"], kind="stable").reset_index(drop=True)


def load_reviewed_claims(path: Path | None, proposed: pd.DataFrame) -> tuple[pd.DataFrame, str]:
    """Apply only reviewer-approved edits to claims with an unchanged identity."""
    if path is None:
        return proposed.copy(), ""
    path = Path(path)
    if not path.is_file():
        raise ValueError("reviewed claims file does not exist")
    content = path.read_bytes()
    try:
        corrections = pd.read_csv(path, dtype=str, keep_default_na=False)
    except Exception:
        raise ValueError("reviewed claims file is invalid") from None
    if len(corrections.columns) != len(CLAIM_COLUMNS) or set(corrections.columns) != set(CLAIM_COLUMNS):
        raise ValueError("reviewed claims file has invalid claim columns")
    corrections = corrections[CLAIM_COLUMNS]
    if corrections["EVIDENCE_ID"].duplicated().any():
        raise ValueError("reviewed claims file contains duplicate EVIDENCE_ID values")
    if not set(corrections["CLAIM_STATUS"]).issubset(CLAIM_STATUSES):
        raise ValueError("reviewed claims file has invalid CLAIM_STATUS")
    if not set(corrections["REVIEW_STATUS"]).issubset({"AUTOMATED", "REVIEWED"}):
        raise ValueError("reviewed claims file has invalid REVIEW_STATUS")
    allowed = {"CLAIM_STATUS", "APPLICABILITY", "LIMITATIONS", "REVIEW_STATUS"}
    immutable = [column for column in CLAIM_COLUMNS if column not in allowed]
    indexed = proposed.set_index("EVIDENCE_ID", drop=False)
    result = proposed.copy()
    for _, correction in corrections.iterrows():
        evidence_id = correction["EVIDENCE_ID"]
        if evidence_id not in indexed.index:
            raise ValueError("reviewed claim does not match a proposed EVIDENCE_ID")
        original = indexed.loc[evidence_id]
        if any(correction[column] != _clean(original[column]) for column in immutable):
            raise ValueError("reviewed claim changed source identity or passage hash")
        row_index = result.index[result["EVIDENCE_ID"].eq(evidence_id)][0]
        for column in allowed:
            result.loc[row_index, column] = correction[column]
    return result, sha256_bytes(content)
def _reason_list(value: Any) -> list[str]:
    return list(dict.fromkeys(item for item in re.split(r"[;,|]", _clean(value)) if item))


def _claim_study_id(claim: pd.Series) -> str:
    return _clean(claim.get("NCT_ID")) or _clean(claim.get("PMID")) or _clean(claim.get("DOI")) or _clean(claim.get("SOURCE_RECORD_ID"))


def _priority_band(claims: pd.DataFrame) -> str:
    if claims.empty:
        return ""
    designs = " ".join(claims["STUDY_DESIGN"].map(_clean).str.casefold())
    human = bool(re.search(r"\b(human|clinical|patient|case report|randomi[sz]ed|controlled)\b", designs))
    if human:
        return "A" if len({_claim_study_id(row) for _, row in claims.iterrows() if _claim_study_id(row)}) >= 2 else "B"
    preclinical = bool(re.search(r"\b(animal|mouse|mice|rat|in vitro|cell|preclinical|model)\b", designs))
    return "C" if preclinical else ""


def assess_drug_hpo_rows(
    drug_rows: pd.DataFrame,
    rationales: pd.DataFrame,
    claims: pd.DataFrame,
    source_status: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Assess separate mechanism and phenotype routes, then rank eligible source rows."""
    drugs = drug_rows.copy()
    if "ROW_ID" not in drugs:
        drugs = assign_input_row_ids(drugs)
    normalized = normalize_drug_inventory(drugs)
    claims = claims.reindex(columns=CLAIM_COLUMNS)
    rationale_by_id = rationales.set_index("RATIONALE_ID", drop=False) if not rationales.empty and "RATIONALE_ID" in rationales else pd.DataFrame()
    records: list[dict[str, Any]] = []
    assessment_records: list[dict[str, Any]] = []
    rank_data: list[tuple[int, tuple[Any, ...]]] = []

    for _, row in normalized.iterrows():
        rationale_id = _clean(row.get("RATIONALE_ID"))
        rationale = rationale_by_id.loc[rationale_id] if not rationale_by_id.empty and rationale_id in rationale_by_id.index else pd.Series(dtype=object)
        if isinstance(rationale, pd.DataFrame):
            rationale = rationale.iloc[0]
        drug_id, target = _clean(row.get("CANONICAL_DRUG_ID")), _clean(row.get("TARGET_SYMBOL"))
        related = claims.loc[claims.get("CANONICAL_DRUG_ID", pd.Series(index=claims.index, dtype=str)).map(_clean).eq(drug_id)].copy() if not claims.empty else claims.copy()
        if "TARGET_SYMBOL" in related:
            related["_TARGET_MATCH"] = related["TARGET_SYMBOL"].map(_clean).isin({"", target})
        else:
            related["_TARGET_MATCH"] = False
        mechanisms = related.loc[
            related["_TARGET_MATCH"] & related["CLAIM_TYPE"].eq("MECHANISM_OBSERVED")
            & related["CLAIM_STATUS"].eq("QUALIFYING") & related["OBSERVED_OR_INFERRED"].eq("OBSERVED")
            & related["ATTRIBUTION"].ne("REGIMEN_ONLY")
            & related["OUTCOME_DIRECTION"].map(_clean).str.upper().eq(classify_action_direction(_clean(row.get("ACTION_TYPE"))))
        ]
        desired = _clean(rationale.get("DESIRED_TARGET_CHANGE"))
        disease = _clean(rationale.get("DISEASE_STATE_DIRECTION")) or "UNKNOWN"
        directional = assess_target_action(desired, _clean(row.get("ACTION_TYPE")) or None, disease_direction=disease)
        mechanism_fit = directional["directionality_fit"]
        mechanism_ok = bool(mechanisms.empty is False and mechanism_fit == "COMPATIBLE_DIRECTION")
        mechanism_reject = mechanism_fit == "OPPOSITE_ACTION" and not related.loc[
            related["_TARGET_MATCH"] & related["CLAIM_TYPE"].eq("MECHANISM_OBSERVED")
            & related["CLAIM_STATUS"].eq("QUALIFYING") & related["OBSERVED_OR_INFERRED"].eq("OBSERVED")
        ].empty

        benefits = related.loc[
            related["CLAIM_TYPE"].eq("BENEFIT_OBSERVED") & related["CLAIM_STATUS"].eq("QUALIFYING")
            & related["OBSERVED_OR_INFERRED"].eq("OBSERVED")
            & related["APPLICABILITY"].isin(["direct", "analogous"])
            & ~related["ATTRIBUTION"].eq("REGIMEN_ONLY")
        ].copy()
        worsening = related.loc[
            related["CLAIM_TYPE"].eq("WORSENING_OBSERVED") & related["CLAIM_STATUS"].eq("QUALIFYING")
            & related["APPLICABILITY"].isin(["direct", "analogous"])
        ]
        hazard = related.loc[
            related["CLAIM_TYPE"].eq("HAZARD") & related["CLAIM_STATUS"].eq("QUALIFYING")
        ]
        uncertain_hazard = related.loc[related["CLAIM_TYPE"].eq("HAZARD") & related["APPLICABILITY"].eq("unknown")]
        uncertain_worsening = related.loc[
            related["CLAIM_TYPE"].eq("WORSENING_OBSERVED") & related["APPLICABILITY"].eq("unknown")
        ]
        contradict = related.loc[related["CLAIM_STATUS"].eq("CONTRADICTED")]
        safety_complete = False
        if not source_status.empty and {"CANONICAL_DRUG_ID", "SOURCE", "COMPLETE"}.issubset(source_status.columns):
            safety_rows = source_status.loc[source_status["CANONICAL_DRUG_ID"].map(_clean).eq(drug_id) & source_status["SOURCE"].map(_clean).str.upper().eq("FDA")]
            safety_complete = not safety_rows.empty and safety_rows["COMPLETE"].map(_clean).str.lower().isin({"true", "1"}).all()
        conflict = not contradict.empty or (not benefits.empty and not worsening.empty) or not uncertain_hazard.empty or not uncertain_worsening.empty
        hard_hazard = not hazard.empty or not worsening.empty

        phenotype_ok = not benefits.empty and bool(benefits["HPO_CODE"].map(_clean).any()) and bool(benefits["INTERVENTION"].map(_clean).any()) and bool(benefits["MODEL_OR_POPULATION"].map(_clean).any())
        # No unknown direction or inferred mechanism can qualify the mechanism route.
        routes = []
        if mechanism_ok:
            routes.append(("MECHANISM_DIRECTED", mechanisms))
        if phenotype_ok:
            routes.append(("PHENOTYPE_DIRECTED", benefits))
        decision = _clean(row.get("DECISION")) or "HOLD"
        reasons = _reason_list(row.get("REASON_CODES"))
        if not hazard.empty:
            decision = "REJECT"
            reasons.append("CITED_CONTEXTUAL_HAZARD")
            routes = []
        elif conflict:
            decision = "HOLD"
            reasons.append("EVIDENCE_CONFLICT")
            routes = []
        elif hard_hazard:
            decision = "REJECT"
            reasons.append("CITED_CONTEXTUAL_HAZARD")
            routes = []
        elif mechanism_reject and not phenotype_ok:
            decision = "REJECT"
        elif routes and decision != "ADVANCE_RESEARCH":
            decision = "TEST_FIRST"
        route_labels = [label for label, _ in routes]
        selected = pd.concat([frame for _, frame in routes]).drop_duplicates("EVIDENCE_ID") if routes else pd.DataFrame(columns=claims.columns)
        band = _priority_band(selected)
        if routes and not safety_complete:
            reasons.append("SAFETY_UNASSESSED")
        if "PHENOTYPE_DIRECTED" in route_labels:
            reasons.extend(["PHENOTYPE_DIRECTED", "DIRECTION_UNKNOWN"] if directional["directionality_fit"] != "COMPATIBLE_DIRECTION" else ["PHENOTYPE_DIRECTED"])
        if "MECHANISM_DIRECTED" in route_labels and "PHENOTYPE_DIRECTED" not in route_labels:
            reasons.append("MECHANISM_ONLY")
        if not routes and decision not in {"ADVANCE_RESEARCH"}:
            band = ""
        eligible = decision in {"TEST_FIRST", "ADVANCE_RESEARCH"} and bool(routes or _clean(row.get("DECISION")) == "ADVANCE_RESEARCH")
        priority_route = "PHENOTYPE_DIRECTED" if "PHENOTYPE_DIRECTED" in route_labels else (route_labels[0] if route_labels else "")
        applicability_order = {"direct": 0, "analogous": 1, "unknown": 2}
        applicability = min((applicability_order.get(_clean(v), 2) for v in selected.get("APPLICABILITY", [])), default=2)
        study_ids = {_claim_study_id(claim) for _, claim in selected.iterrows() if _claim_study_id(claim)}
        fit_order = 0 if mechanism_fit == "COMPATIBLE_DIRECTION" else 1
        status_rows = source_status.loc[source_status.get("CANONICAL_DRUG_ID", pd.Series(index=source_status.index, dtype=str)).map(_clean).eq(drug_id)] if not source_status.empty else source_status
        completeness = "complete"
        if not status_rows.empty and "COMPLETE" in status_rows:
            completeness = "complete" if status_rows["COMPLETE"].map(_clean).str.lower().isin({"true", "1"}).all() else "incomplete"
        rank_data.append((len(records), (0 if eligible else 1, {"A": 0, "B": 1, "C": 2, "": 3}.get(band, 3), applicability, 0 if len(study_ids) > 1 else 1, fit_order, 0 if safety_complete else 1, 0 if completeness == "complete" else 1, _clean(row.get("ROW_ID")))))

        output = row.to_dict()
        for column in ("DECISION", "REASON_CODES", "DIRECTIONALITY_FIT", "HAZARD_STATUS", "EVIDENCE_GAPS"):
            output[f"ORIGINAL_{column}"] = _clean(row.get(column))
        output.update({
            "DECISION": decision,
            "REASON_CODES": ";".join(dict.fromkeys(reasons)),
            "DIRECTIONALITY_FIT": mechanism_fit if mechanism_fit != "UNKNOWN" else "UNKNOWN",
            "HAZARD_STATUS": "cited_hazard" if hard_hazard else ("no_cited_hazard" if safety_complete else "unknown"),
            "PROMOTION_ROUTE": ";".join(route_labels),
            "MECHANISM_EVIDENCE_IDS": ";".join(sorted(mechanisms["EVIDENCE_ID"].map(_clean))) if "MECHANISM_DIRECTED" in route_labels else "",
            "PHENOTYPE_EVIDENCE_IDS": ";".join(sorted(benefits["EVIDENCE_ID"].map(_clean))) if "PHENOTYPE_DIRECTED" in route_labels else "",
            "PRIORITY_BAND": band if eligible else "",
            "PRIORITY_APPLICABILITY": ("direct", "analogous", "unknown")[applicability],
            "PRIORITY_REPLICATED": len(study_ids) > 1,
            "PRIORITY_DIRECTIONAL_FIT": "matched" if mechanism_fit == "COMPATIBLE_DIRECTION" else "unknown",
            "PRIORITY_SAFETY": "assessed" if safety_complete else "unassessed",
            "PRIORITY_SOURCE_COMPLETENESS": completeness,
        })
        records.append(output)
        hpo_codes = sorted(set(related["HPO_CODE"].map(_clean)) - {""}) if "HPO_CODE" in related else [""]
        for hpo in hpo_codes or [""]:
            assessment_records.append({
                "ROW_ID": _clean(row.get("ROW_ID")), "RATIONALE_ID": rationale_id, "CANONICAL_DRUG_ID": drug_id,
                "HPO_CODE": hpo, "MECHANISM_DECISION": "TEST_FIRST" if mechanism_ok and not hard_hazard and not conflict else ("REJECT" if mechanism_reject else "HOLD"),
                "PHENOTYPE_DECISION": "TEST_FIRST" if phenotype_ok and not hard_hazard and not conflict else "HOLD",
                "MECHANISM_EVIDENCE_IDS": ";".join(sorted(mechanisms["EVIDENCE_ID"].map(_clean))),
                "PHENOTYPE_EVIDENCE_IDS": ";".join(sorted(benefits.loc[benefits["HPO_CODE"].map(_clean).eq(hpo), "EVIDENCE_ID"].map(_clean))) if "HPO_CODE" in benefits else "",
                "DECISION": decision, "REASON_CODES": ";".join(dict.fromkeys(reasons)),
            })
    enriched = pd.DataFrame(records)
    order = sorted(rank_data, key=lambda item: item[1])
    ranks = {index: rank + 1 for rank, (index, _) in enumerate(order)}
    enriched["PRIORITY_RANK"] = [ranks[index] for index in range(len(enriched))]
    assessments = pd.DataFrame(assessment_records, columns=["ROW_ID", "RATIONALE_ID", "CANONICAL_DRUG_ID", "HPO_CODE", "MECHANISM_DECISION", "PHENOTYPE_DECISION", "MECHANISM_EVIDENCE_IDS", "PHENOTYPE_EVIDENCE_IDS", "DECISION", "REASON_CODES"])
    return enriched, assessments
_OUTPUT_FILES = (
    "hpo_dictionary.csv", "drug_identity.csv", "search_runs.csv",
    "public_evidence_claims.csv", "drug_hpo_assessments.csv",
    "target_drug_review_enriched.csv", "target_drug_review_decision_ready.csv",
    "source_status.csv", "drug_hpo_report.md", "drug_hpo_provenance.json",
)
_EXTRACTION_RULE_VERSION = "drug-hpo-extraction-1"


def _frame_csv(frame: pd.DataFrame, path: Path) -> None:
    frame.to_csv(path, index=False, lineterminator="\n")


def _safe_envelope(envelope: dict[str, Any], *, identity: dict[str, str], source: str, hpo_code: str = "") -> dict[str, Any]:
    return {
        "SOURCE": source, "CANONICAL_DRUG_ID": identity.get("CANONICAL_DRUG_ID", ""),
        "HPO_CODE": hpo_code, "STATUS": str(envelope.get("status", "UNAVAILABLE")),
        "COMPLETE": bool(envelope.get("complete", False)),
        "REQUEST_ID": str(envelope.get("request_id", "")),
        "RESPONSE_SHA256": str(envelope.get("response_sha256", "")),
        "RETRIEVED_AT": str(envelope.get("retrieved_at", "")),
        "REASON": str(envelope.get("reason", "")),
    }



def _task_exception_envelope(source: str, error: Exception) -> dict[str, Any]:
    return {
        "status": "TASK_EXCEPTION",
        "records": [],
        "complete": False,
        "source": source,
        "request_id": "",
        "response_sha256": "",
        "retrieved_at": "",
        "reason": f"task exception: {type(error).__name__}",
    }


def _unresolved_identity_envelope(source: str) -> dict[str, Any]:
    return {
        "status": "IDENTITY_UNRESOLVED",
        "records": [],
        "complete": False,
        "source": source,
        "request_id": "",
        "response_sha256": "",
        "retrieved_at": "",
        "reason": "unresolved public drug identity",
    }


def _report(enriched: pd.DataFrame, claims: pd.DataFrame, statuses: pd.DataFrame,
            hpo_count: int, decision_ready_count: int) -> str:
    decisions = enriched.get("DECISION", pd.Series(dtype=str)).value_counts().to_dict()
    routes = enriched.get("PROMOTION_ROUTE", pd.Series(dtype=str)).fillna("")
    route_counts = {
        "mechanism_qualified": int(routes.str.contains("MECHANISM_DIRECTED", regex=False).sum()),
        "phenotype_qualified": int(routes.str.contains("PHENOTYPE_DIRECTED", regex=False).sum()),
        "both_routes": int((routes.str.contains("MECHANISM_DIRECTED", regex=False) & routes.str.contains("PHENOTYPE_DIRECTED", regex=False)).sum()),
    }
    bands = enriched.get("PRIORITY_BAND", pd.Series(dtype=str)).value_counts().to_dict()
    reasons = enriched.get("REASON_CODES", pd.Series(dtype=str)).fillna("")
    incomplete = int((~statuses.get("COMPLETE", pd.Series(dtype=bool)).map(lambda x: str(x).lower() in {"true", "1"})).sum()) if not statuses.empty else 0
    return "\n".join([
        "# Drug–HPO evidence enrichment report", "",
        "> Research-only prioritization; not diagnosis, treatment advice, or evidence of MVA efficacy.", "",
        "Automated decisions and claim extraction are deterministic screening aids, not a systematic review.",
        "Evidence review is targeted and may omit relevant sources or unpublished findings.", "",
        "## Counts", "",
        f"- Input drug rows: {len(enriched)}",
        f"- HPO terms: {hpo_count}",
        f"- Evidence claims: {len(claims)}",
        f"- Decision-ready rows: {decision_ready_count}",
        f"- Decision transitions: {sum(1 for _, row in enriched.iterrows() if _clean(row.get('DECISION')) != _clean(row.get('ORIGINAL_DECISION')))}",
        f"- Decisions: {json.dumps({str(k): int(v) for k, v in sorted(decisions.items())}, sort_keys=True)}",
        f"- Routes: {json.dumps(route_counts, sort_keys=True)}",
        f"- Priority bands: {json.dumps({str(k): int(v) for k, v in sorted(bands.items())}, sort_keys=True)}",
        f"- Hazard/conflict rows: {int(reasons.str.contains('CITED_CONTEXTUAL_HAZARD', regex=False).sum())}/{int(reasons.str.contains('EVIDENCE_CONFLICT', regex=False).sum())}",
        f"- Incomplete source scopes: {incomplete}/{len(statuses)}", "",
    ])


def run_drug_hpo_evidence(
    compass_dir: Path, hpo_codes: Path, output_dir: Path, *, cache_dir: Path,
    reviewed_claims: Path | None = None, resume: bool = False, offline: bool = False,
    transport: Any = None, max_workers: int = 8,
) -> dict[str, Any]:
    """Enrich an immutable Compass snapshot and atomically publish data-blind outputs."""
    from . import drug_hpo_sources as sources
    compass_dir, hpo_codes, output_dir, cache_dir = map(Path, (compass_dir, hpo_codes, output_dir, cache_dir))
    if max_workers < 1:
        raise ValueError("max_workers must be at least 1")
    drugs, rationales, compass_provenance, compass_hashes = load_compass_snapshot(compass_dir)
    codes, hpo_hash = load_hpo_codes(hpo_codes)
    normalized = normalize_drug_inventory(assign_input_row_ids(drugs))
    identity = {
        "compass_artifact_sha256": compass_hashes, "hpo_codes_sha256": hpo_hash,
        "reviewed_claims_sha256": sha256_file(Path(reviewed_claims)) if reviewed_claims else "",
        "offline": bool(offline), "extraction_rule_version": _EXTRACTION_RULE_VERSION,
    }
    fixed_identity = {k: v for k, v in identity.items() if k != "offline"}
    if output_dir.exists():
        if not resume:
            raise ValueError("output directory already exists; use --resume or choose a new destination")
        try:
            previous = json.loads((output_dir / "drug_hpo_provenance.json").read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            raise ValueError("existing output has invalid drug-HPO provenance") from None
        if previous.get("identity") != fixed_identity:
            raise ValueError("resumed inputs differ from the recorded drug-HPO run; use a new output destination")
        return previous["summary"]

    output_dir.parent.mkdir(parents=True, exist_ok=True)
    cache_dir.mkdir(parents=True, exist_ok=True)
    staging = output_dir.parent / f".{output_dir.name}.staging"
    if staging.exists():
        raise ValueError("staging directory already exists; remove it or choose a new output destination")
    staging.mkdir()
    try:
        hpo_records, hpo_status, ontology_hashes = [], [], []
        hpo_tasks = [
            (ordinal, "OLS4", "", code, lambda code=code: sources.resolve_hpo_term(code, cache_dir, offline=offline, transport=transport))
            for ordinal, code in enumerate(codes)
        ]

        def collect(tasks: list[tuple[int, str, str, str, Any]]) -> list[tuple[tuple[int, str, str, str], dict[str, Any]]]:
            completed: list[tuple[tuple[int, str, str, str], dict[str, Any]]] = []
            with ThreadPoolExecutor(max_workers=max_workers) as executor:
                futures = {executor.submit(task): metadata for *metadata, task in tasks}
                for future in as_completed(futures):
                    metadata = futures[future]
                    try:
                        envelope = future.result()
                    except Exception as error:
                        envelope = _task_exception_envelope(metadata[1], error)
                    completed.append((metadata, envelope))
            return sorted(completed, key=lambda item: item[0][0])

        for (_, source, _, code), env in collect(hpo_tasks):
            hpo_status.append(_safe_envelope(env, identity={}, source=source, hpo_code=code))
            for term in env.get("records", []):
                record = dict(term)
                record.setdefault("HPO_CODE", code)
                hpo_records.append(record)
                ontology_hashes.append({"HPO_CODE": code, "RESPONSE_SHA256": env.get("response_sha256", ""), "ONTOLOGY_VERSION": record.get("ONTOLOGY_VERSION", "")})
        hpo_terms = pd.DataFrame(hpo_records, columns=["HPO_CODE", "HPO_LABEL", "HPO_SYNONYMS", "IRI", "ONTOLOGY_VERSION", "ONTOLOGY_UPDATED"])
        records_by_source: dict[str, list[dict[str, Any]]] = {"literature": [], "fda": [], "trials": []}
        status_records = list(hpo_status)
        aliases_by_drug = {
            did: list(dict.fromkeys(alias for aliases in group["ALIASES"] for alias in _terms(aliases)))
            for did, group in normalized.groupby("CANONICAL_DRUG_ID", sort=False)
        }
        distinct = normalized.drop_duplicates("CANONICAL_DRUG_ID", keep="first").sort_values("CANONICAL_DRUG_ID", kind="stable")
        source_tasks: list[tuple[int, str, str, str, Any]] = []
        ordinal = len(hpo_tasks)
        for _, drug in distinct.iterrows():
            did = _clean(drug["CANONICAL_DRUG_ID"])
            names = [name for name in aliases_by_drug[did] if not _CHEMBL_ID.fullmatch(name)]
            if not names:
                source_tasks.append((ordinal, "FDA", did, "", lambda: _unresolved_identity_envelope("FDA")))
                ordinal += 1
                for _, term in hpo_terms.sort_values("HPO_CODE", kind="stable").iterrows():
                    code = _clean(term["HPO_CODE"])
                    source_tasks.append((ordinal, "LITERATURE", did, code, lambda: _unresolved_identity_envelope("LITERATURE")))
                    ordinal += 1
                    source_tasks.append((ordinal, "TRIALS", did, code, lambda: _unresolved_identity_envelope("TRIALS")))
                    ordinal += 1
                continue
            drug_name, drug_aliases = names[0], names[1:]
            source_tasks.append((ordinal, "FDA", did, "", lambda drug_name=drug_name, drug_aliases=drug_aliases: sources.fetch_fda_labels(drug_name, drug_aliases, cache_dir, offline=offline, transport=transport)))
            ordinal += 1
            for _, term in hpo_terms.sort_values("HPO_CODE", kind="stable").iterrows():
                code, label, synonyms = _clean(term["HPO_CODE"]), _clean(term["HPO_LABEL"]), _terms(term.get("HPO_SYNONYMS"))
                source_tasks.append((ordinal, "LITERATURE", did, code, lambda drug_name=drug_name, drug_aliases=drug_aliases, code=code, label=label, synonyms=synonyms: sources.search_drug_hpo_literature(drug_name, drug_aliases, code, label, synonyms, cache_dir, offline=offline, transport=transport)))
                ordinal += 1
                source_tasks.append((ordinal, "TRIALS", did, code, lambda drug_name=drug_name, drug_aliases=drug_aliases, label=label, synonyms=synonyms: sources.fetch_trials(drug_name, drug_aliases, label, synonyms, cache_dir, offline=offline, transport=transport)))
                ordinal += 1
        for (_, source, did, code), env in collect(source_tasks):
            status_records.append(_safe_envelope(env, identity={"CANONICAL_DRUG_ID": did}, source=source, hpo_code=code))
            record_key = source.lower() if source != "FDA" else "fda"
            records_by_source[record_key].extend(
                {**record, "_SCOPE_CANONICAL_DRUG_ID": did, "_SCOPE_HPO_CODE": code}
                for record in env.get("records", [])
            )
        source_status = pd.DataFrame(status_records, columns=["SOURCE", "CANONICAL_DRUG_ID", "HPO_CODE", "STATUS", "COMPLETE", "REQUEST_ID", "RESPONSE_SHA256", "RETRIEVED_AT", "REASON"])
        searches = build_search_inventory(normalized, hpo_terms)
        fda_searches = [
            ["fda", _clean(drug.get("CANONICAL_DRUG_ID")), "", _query(aliases_by_drug[_clean(drug.get("CANONICAL_DRUG_ID"))]), "not_searched"]
            for _, drug in distinct.iterrows()
        ]
        if fda_searches:
            searches = pd.concat([searches, pd.DataFrame(fda_searches, columns=searches.columns)], ignore_index=True)
            searches = searches.sort_values(["SOURCE", "CANONICAL_DRUG_ID", "HPO_CODE", "QUERY"], kind="stable").reset_index(drop=True)
        searches["STATUS"] = "not_searched"
        if not source_status.empty:
            failed_keys = {(r.CANONICAL_DRUG_ID, r.SOURCE, r.HPO_CODE) for r in source_status.itertuples() if not r.COMPLETE}
            for i, row in searches.iterrows():
                if (row.CANONICAL_DRUG_ID, row.SOURCE.upper() if row.SOURCE != "literature" else "LITERATURE", row.HPO_CODE) in failed_keys:
                    searches.at[i, "STATUS"] = "retrieval_failed"
                elif row.SOURCE in {"literature", "trials"}:
                    searches.at[i, "STATUS"] = "screening_pending"
        proposed = extract_claims(normalized, hpo_terms, records_by_source["literature"], records_by_source["fda"], records_by_source["trials"])
        claims, reviewed_hash = load_reviewed_claims(reviewed_claims, proposed)
        for index, scope in searches.iterrows():
            source_name = {"literature": "LITERATURE", "fda": "FDA", "trials": "TRIALS"}[scope["SOURCE"]]
            matching = source_status.loc[
                source_status["SOURCE"].eq(source_name)
                & source_status["CANONICAL_DRUG_ID"].eq(scope["CANONICAL_DRUG_ID"])
                & source_status["HPO_CODE"].eq(scope["HPO_CODE"])
            ]
            if matching.empty or not matching["COMPLETE"].astype(bool).all():
                searches.at[index, "STATUS"] = "retrieval_failed"
                continue
            source_type = {"LITERATURE": "LITERATURE", "FDA": "FDA", "TRIALS": "TRIAL"}[source_name]
            relevant = claims.loc[claims["CANONICAL_DRUG_ID"].eq(scope["CANONICAL_DRUG_ID"]) & claims["SOURCE_TYPE"].eq(source_type)]
            if scope["HPO_CODE"]:
                relevant = relevant.loc[relevant["HPO_CODE"].eq(scope["HPO_CODE"])]
            searches.at[index, "STATUS"] = "evidence_coded" if not relevant.empty else "reviewed_no_qualifying_evidence"
        enriched, assessments = assess_drug_hpo_rows(normalized, rationales, claims, source_status)
        decisions = enriched.get("DECISION", pd.Series(dtype=str))
        decision_ready = enriched.loc[decisions.isin(["TEST_FIRST", "ADVANCE_RESEARCH"])].copy()
        frames = {
            "hpo_dictionary.csv": hpo_terms, "drug_identity.csv": normalized,
            "search_runs.csv": searches, "public_evidence_claims.csv": claims,
            "drug_hpo_assessments.csv": assessments, "target_drug_review_enriched.csv": enriched,
            "target_drug_review_decision_ready.csv": decision_ready, "source_status.csv": source_status,
        }
        for filename, frame in frames.items():
            _frame_csv(frame, staging / filename)
        (staging / "drug_hpo_report.md").write_text(_report(enriched, claims, source_status, len(hpo_terms), len(decision_ready)), encoding="utf-8")
        output_hashes = {name: sha256_file(staging / name) for name in (*frames.keys(), "drug_hpo_report.md")}
        failures = source_status.loc[~source_status["COMPLETE"].astype(bool)] if not source_status.empty else source_status
        summary = {
            "status": "source_incomplete" if not failures.empty else "completed", "output_dir": str(output_dir),
            "input_drug_rows": len(drugs), "hpo_terms": len(hpo_terms), "evidence_claims": len(claims),
            "test_first": int(decisions.eq("TEST_FIRST").sum()), "advance_research": int(decisions.eq("ADVANCE_RESEARCH").sum()),
            "hold": int(decisions.eq("HOLD").sum()), "reject": int(decisions.eq("REJECT").sum()),
            "mechanism_qualified": int(enriched["PROMOTION_ROUTE"].str.contains("MECHANISM_DIRECTED", regex=False).sum()),
            "phenotype_qualified": int(enriched["PROMOTION_ROUTE"].str.contains("PHENOTYPE_DIRECTED", regex=False).sum()),
            "both_routes": int((enriched["PROMOTION_ROUTE"].str.contains("MECHANISM_DIRECTED", regex=False) & enriched["PROMOTION_ROUTE"].str.contains("PHENOTYPE_DIRECTED", regex=False)).sum()),
            "priority_bands": {band: int(enriched["PRIORITY_BAND"].eq(band).sum()) for band in ("A", "B", "C")},
            "source_failures": len(failures),
        }
        provenance = {
            "schema_version": "drug-hpo-evidence-1", "generated_at": datetime.now(timezone.utc).isoformat(),
            "identity": fixed_identity, "compass_artifact_sha256": compass_hashes,
            "hpo_codes_sha256": hpo_hash, "ontology_records": ontology_hashes,
            "reviewed_claims_sha256": reviewed_hash, "source_requests": source_status[["SOURCE", "REQUEST_ID", "RESPONSE_SHA256"]].to_dict("records"),
            "extraction_rule_version": _EXTRACTION_RULE_VERSION,
            "source_failures": failures[["SOURCE", "CANONICAL_DRUG_ID", "HPO_CODE", "STATUS", "REASON"]].to_dict("records") if not failures.empty else [],
            "output_sha256": output_hashes, "summary": {key: value for key, value in summary.items() if key != "output_dir"},
        }
        (staging / "drug_hpo_provenance.json").write_text(json.dumps(provenance, sort_keys=True, indent=2) + "\n", encoding="utf-8")
        staging.rename(output_dir)
        return summary
    except Exception:
        if staging.exists():
            shutil.rmtree(staging)
        raise
