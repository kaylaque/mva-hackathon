from __future__ import annotations

import hashlib
import sys
from pathlib import Path

import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).parents[1] / "src"))

from mva_hackathon.drug_hpo_evidence import (
    CLAIM_COLUMNS,
    build_search_inventory,
    extract_claims,
    load_reviewed_claims,
    normalize_drug_inventory,
)


def _drugs() -> pd.DataFrame:
    return pd.DataFrame([
        {"ROW_ID": "row-1", "DRUG_ID": "CHEMBL1", "DRUG_NAME": "Drug A", "TARGET_SYMBOL": "GENE1", "ACTION_TYPE": "INHIBITOR", "ACTION_DIRECTION": "DECREASE", "RATIONALE_ID": "R1"},
        {"ROW_ID": "row-2", "DRUG_ID": "CHEMBL2", "DRUG_NAME": "Drug B", "TARGET_SYMBOL": "GENE2", "ACTION_TYPE": "AGONIST", "ACTION_DIRECTION": "INCREASE", "RATIONALE_ID": "R2"},
        {"ROW_ID": "row-3", "DRUG_ID": "", "DRUG_NAME": "Drug C", "TARGET_SYMBOL": "GENE3", "ACTION_TYPE": "INHIBITOR", "ACTION_DIRECTION": "DECREASE", "RATIONALE_ID": "R3"},
    ])


def test_normalize_drug_inventory_preserves_rows_and_models_regimens() -> None:
    rows = pd.DataFrame([
        {"ROW_ID": "r1", "DRUG_ID": "CHEMBL1", "DRUG_NAME": "Drug A", "SOURCE_RECORD_IDS": "s1", "TARGET_SYMBOL": "G1"},
        {"ROW_ID": "r2", "DRUG_ID": "CHEMBL1", "DRUG_NAME": "Drug A alias", "SOURCE_RECORD_IDS": "s2", "TARGET_SYMBOL": "G2"},
        {"ROW_ID": "r3", "DRUG_ID": "", "DRUG_NAME": "  Drug C  ", "SOURCE_RECORD_IDS": "s3", "TARGET_SYMBOL": "G3", "COMPONENT_DRUG_IDS": "CHEMBL1;CHEMBL2"},
    ])

    result = normalize_drug_inventory(rows)

    assert result.ROW_ID.tolist() == ["r1", "r2", "r3"]
    assert result.CANONICAL_DRUG_ID.tolist() == ["CHEMBL1", "CHEMBL1", "drug c"]
    assert result.loc[0, "ALIASES"] == ["Drug A"]
    assert result.loc[2, "REGIMEN_ID"]
    assert result.loc[2, "COMPONENT_DRUG_IDS"] == ["CHEMBL1", "CHEMBL2"]


def test_search_inventory_has_stable_drug_and_hpo_scopes() -> None:
    drugs = normalize_drug_inventory(_drugs())
    terms = pd.DataFrame([{"HPO_CODE": "HP:0000002", "HPO_LABEL": "Ataxia", "HPO_SYNONYMS": ["Gait issue"]}, {"HPO_CODE": "HP:0000001", "HPO_LABEL": "Seizure", "HPO_SYNONYMS": []}])

    result = build_search_inventory(drugs, terms)

    assert result.columns.tolist() == ["SOURCE", "CANONICAL_DRUG_ID", "HPO_CODE", "QUERY", "STATUS"]
    assert result.STATUS.eq("not_searched").all()
    assert len(result[result.SOURCE.eq("literature")]) == 3 * 3
    assert len(result[result.SOURCE.eq("trials")]) == 3 * 3
    assert result.duplicated(["SOURCE", "CANONICAL_DRUG_ID", "HPO_CODE", "QUERY"]).sum() == 0


def test_extract_claims_links_exact_phenotype_passages_and_never_inferrs_from_title() -> None:
    drugs = normalize_drug_inventory(_drugs().iloc[:1])
    terms = pd.DataFrame([{"HPO_CODE": "HP:0000001", "HPO_LABEL": "Ataxia", "HPO_SYNONYMS": ["gait instability"]}])
    literature = [{"PMID": "123", "TITLE": "Drug A cures ataxia", "ABSTRACT": "Drug A was administered. Ataxia improved after treatment.", "PASSAGE_SOURCE": "ABSTRACT", "PUBLICATION_TYPES": ["Clinical Trial"], "RETRACTED": False}, {"PMID": "124", "TITLE": "", "ABSTRACT": "Drug A resolved ataxia.", "PUBLICATION_TYPES": [], "RETRACTED": True}]
    fda = [{"FDA_SET_ID": "set1", "GENERIC_NAMES": ["Drug A"], "CONTRAINDICATIONS": "Ataxia worsened with Drug A.", "WARNINGS": "Drug A may cause gait instability."}]
    trials = [{"NCT_ID": "NCT00000001", "INTERVENTIONS": ["Drug A"], "CONDITIONS": ["Ataxia"], "BRIEF_SUMMARY": "Drug A was studied for ataxia.", "RESULTS_POSTED": False}]

    result = extract_claims(drugs, terms, literature, fda, trials)

    assert result.columns.tolist() == CLAIM_COLUMNS
    benefit = result[(result.SOURCE_TYPE == "LITERATURE") & (result.CLAIM_TYPE == "BENEFIT_OBSERVED")].iloc[0]
    assert benefit.CLAIM_STATUS == "QUALIFYING"
    assert benefit.OUTCOME_DIRECTION == "BENEFIT"
    assert benefit.PASSAGE_LOCATOR.startswith("ABSTRACT:")
    assert len(benefit.PASSAGE_SHA256) == 64
    assert not ((result.SOURCE_TYPE == "TRIAL") & (result.CLAIM_TYPE == "BENEFIT_OBSERVED")).any()
    trial_context = result[result.SOURCE_TYPE.eq("TRIAL")]
    assert trial_context.CLAIM_TYPE.eq("TRIAL_CONTEXT").all()
    assert trial_context.CLAIM_STATUS.eq("CONTEXT_ONLY").all()
    assert result[result.SOURCE_TYPE == "FDA"].CLAIM_TYPE.isin(["HAZARD", "SAFETY_CONTEXT"]).all()
    assert not result[result.SOURCE_TYPE == "FDA"].CLAIM_TYPE.eq("BENEFIT_OBSERVED").any()
    assert "RETRACTED" in set(result.CLAIM_STATUS)


def test_unmatched_fda_hazards_remain_contextual_without_hpo_attribution() -> None:
    drugs = normalize_drug_inventory(_drugs().iloc[:1])
    terms = pd.DataFrame([{"HPO_CODE": "HP:0000001", "HPO_LABEL": "Ataxia", "HPO_SYNONYMS": []}])
    fda = [{"FDA_SET_ID": "set1", "GENERIC_NAMES": ["Drug A"], "CONTRAINDICATIONS": "Severe liver injury is a contraindication."}]

    result = extract_claims(drugs, terms, [], fda, [])

    assert len(result) == 1
    assert result.iloc[0].CLAIM_TYPE == "HAZARD"
    assert result.iloc[0].CLAIM_STATUS == "CONTEXT_ONLY"
    assert result.iloc[0].HPO_CODE == ""
    assert result.iloc[0].APPLICABILITY == "unknown"


def test_ambiguous_directions_and_negated_benefit_do_not_qualify() -> None:
    drugs = normalize_drug_inventory(_drugs().iloc[:1])
    terms = pd.DataFrame([{"HPO_CODE": "HP:0000001", "HPO_LABEL": "Ataxia", "HPO_SYNONYMS": []}])
    literature = [{"PMID": "1", "TITLE": "", "ABSTRACT": "Drug A increased ataxia. Drug A did not improve ataxia.", "PUBLICATION_TYPES": []}]

    result = extract_claims(drugs, terms, literature, [], [])

    assert set(result.CLAIM_STATUS) <= {"PROPOSED", "CONTEXT_ONLY", "CONTRADICTED"}
    assert not result.CLAIM_STATUS.eq("QUALIFYING").any()


def test_title_only_and_article_count_never_create_benefit_claims() -> None:
    drugs = normalize_drug_inventory(_drugs().iloc[:1])
    terms = pd.DataFrame([{"HPO_CODE": "HP:0000001", "HPO_LABEL": "Ataxia", "HPO_SYNONYMS": []}])
    literature = [{"PMID": "1", "TITLE": "Drug A improved ataxia", "ABSTRACT": "", "PUBLICATION_TYPES": ["Review"]}]
    result = extract_claims(drugs, terms, literature, [], [])
    assert result.empty


def test_mechanism_claim_requires_exact_target_and_action_in_observed_passage() -> None:
    drugs = normalize_drug_inventory(_drugs().iloc[:1])
    terms = pd.DataFrame(columns=["HPO_CODE", "HPO_LABEL", "HPO_SYNONYMS"])
    literature = [
        {"PMID": "42", "TITLE": "", "ABSTRACT": "Drug A inhibited GENE1 activity and reduced pathway activity in cells.", "PUBLICATION_TYPES": ["Journal Article"]},
        {"PMID": "43", "TITLE": "", "ABSTRACT": "Drug A inhibited GENE1 activity and reduced pathway activity in cells.", "PUBLICATION_TYPES": ["Retracted Publication"], "RETRACTED": True},
    ]

    result = extract_claims(drugs, terms, literature, [], [])

    mechanism = result[result.CLAIM_TYPE.eq("MECHANISM_OBSERVED")]
    assert len(mechanism) == 2
    assert mechanism[mechanism.PMID.eq("42")].iloc[0].CLAIM_STATUS == "QUALIFYING"
    assert mechanism[mechanism.PMID.eq("43")].iloc[0].CLAIM_STATUS == "RETRACTED"


def test_reviewed_corrections_only_change_allowed_claim_fields(tmp_path: Path) -> None:
    proposed = pd.DataFrame([{column: "" for column in CLAIM_COLUMNS}])
    proposed.loc[0, ["EVIDENCE_ID", "CLAIM_STATUS", "PASSAGE_SHA256", "REVIEW_STATUS"]] = ["DHEV-123", "PROPOSED", "abc", "AUTOMATED"]
    correction = proposed.copy()
    correction.loc[0, "CLAIM_STATUS"] = "QUALIFYING"
    correction.loc[0, "REVIEW_STATUS"] = "REVIEWED"
    path = tmp_path / "corrections.csv"
    correction.to_csv(path, index=False)

    result, digest = load_reviewed_claims(path, proposed)

    assert result.loc[0, "CLAIM_STATUS"] == "QUALIFYING"
    assert result.loc[0, "REVIEW_STATUS"] == "REVIEWED"
    assert len(digest) == 64
    correction[CLAIM_COLUMNS[::-1]].to_csv(path, index=False)
    reordered, _ = load_reviewed_claims(path, proposed)
    assert reordered.loc[0, "CLAIM_STATUS"] == "QUALIFYING"

    correction.loc[0, "PMID"] = "changed"
    correction.to_csv(path, index=False)
    with pytest.raises(ValueError, match="source identity"):
        load_reviewed_claims(path, proposed)

def test_full_text_claims_keep_locator_and_conflicts_link_both_directions() -> None:
    drugs = normalize_drug_inventory(_drugs().iloc[:1])
    terms = pd.DataFrame([{"HPO_CODE": "HP:0000001", "HPO_LABEL": "Ataxia", "HPO_SYNONYMS": []}])
    literature = [{"PMID": "9", "TITLE": "Drug A improved ataxia", "ABSTRACT": "", "FULL_TEXT_PARAGRAPHS": ["Drug A improved ataxia.", "Drug A worsened ataxia."], "PUBLICATION_TYPES": ["Case Reports"]}]
    result = extract_claims(drugs, terms, literature, [], [])
    phenotype = result[result.HPO_CODE.eq("HP:0000001")]
    assert len(phenotype) == 2
    assert phenotype.CLAIM_STATUS.eq("CONTRADICTED").all()
    assert phenotype.PASSAGE_LOCATOR.str.startswith("FULL_TEXT:p").all()
    assert phenotype.PASSAGE_SHA256.nunique() == 2


def test_regimen_outcome_is_not_attributed_to_component_without_independent_mention() -> None:
    rows = pd.DataFrame([
        {"ROW_ID": "r1", "DRUG_ID": "COMBO", "DRUG_NAME": "Drug A + Drug B", "COMPONENT_DRUG_IDS": ["CHEMBL1", "CHEMBL2"], "TARGET_SYMBOL": "GENE1", "ACTION_TYPE": "INHIBITOR"},
        {"ROW_ID": "r2", "DRUG_ID": "CHEMBL1", "DRUG_NAME": "Drug A", "TARGET_SYMBOL": "GENE2", "ACTION_TYPE": "INHIBITOR"},
        {"ROW_ID": "r3", "DRUG_ID": "CHEMBL2", "DRUG_NAME": "Drug B", "TARGET_SYMBOL": "GENE3", "ACTION_TYPE": "INHIBITOR"},
    ])
    drugs = normalize_drug_inventory(rows)
    terms = pd.DataFrame([{"HPO_CODE": "HP:0000001", "HPO_LABEL": "Ataxia", "HPO_SYNONYMS": []}])
    literature = [{"PMID": "9", "TITLE": "", "ABSTRACT": "Drug A + Drug B improved ataxia.", "PUBLICATION_TYPES": ["Case Reports"]}]
    result = extract_claims(drugs, terms, literature, [], [])
    combo = result[result.CANONICAL_DRUG_ID.eq("COMBO")]
    assert combo.REGIMEN_ID.ne("").all()
    assert combo.ATTRIBUTION.eq("REGIMEN_ONLY").all()
    assert result[result.CANONICAL_DRUG_ID.isin(["CHEMBL1", "CHEMBL2"])].ATTRIBUTION.eq("REGIMEN_ONLY").all()


def test_corrections_reject_unknown_claim_and_invalid_review_status(tmp_path: Path) -> None:
    proposed = pd.DataFrame([{column: "" for column in CLAIM_COLUMNS}])
    proposed.loc[0, ["EVIDENCE_ID", "CLAIM_STATUS", "PASSAGE_SHA256", "REVIEW_STATUS"]] = ["DHEV-123", "PROPOSED", "abc", "AUTOMATED"]
    correction = proposed.copy()
    correction.loc[0, "EVIDENCE_ID"] = "DHEV-unknown"
    path = tmp_path / "corrections.csv"
    correction.to_csv(path, index=False)
    with pytest.raises(ValueError, match="does not match"):
        load_reviewed_claims(path, proposed)
    correction = proposed.copy()
    correction.loc[0, "REVIEW_STATUS"] = "UNREVIEWED"
    correction.to_csv(path, index=False)
    with pytest.raises(ValueError, match="REVIEW_STATUS"):
        load_reviewed_claims(path, proposed)

