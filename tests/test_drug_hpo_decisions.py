from __future__ import annotations

import pandas as pd

from mva_hackathon.compass import DRUG_REVIEW_COLUMNS, RATIONALE_COLUMNS, filter_drug_candidates
from mva_hackathon.drug_hpo_evidence import CLAIM_COLUMNS, assess_drug_hpo_rows


def _inputs():
    rationale = {column: "" for column in RATIONALE_COLUMNS}
    rationale.update({"RATIONALE_ID": "R1", "TARGET_SYMBOL": "T1", "DESIRED_TARGET_CHANGE": "DECREASE", "DISEASE_STATE_DIRECTION": "REDUCED_FUNCTION"})
    row = {column: "" for column in DRUG_REVIEW_COLUMNS}
    row.update({"RATIONALE_ID": "R1", "TARGET_SYMBOL": "T1", "DRUG_ID": "CHEMBL1", "DRUG_NAME": "Drug", "ACTION_TYPE": "INHIBITOR", "ACTION_DIRECTION": "DECREASE", "DIRECTIONALITY_FIT": "UNKNOWN", "HAZARD_STATUS": "unknown", "DECISION": "HOLD", "REASON_CODES": "ORIGINAL_REASON", "EVIDENCE_GAPS": "original gap"})
    claim = {column: "" for column in CLAIM_COLUMNS}
    claim.update({"EVIDENCE_ID": "E1", "CANONICAL_DRUG_ID": "CHEMBL1", "TARGET_SYMBOL": "T1", "HPO_CODE": "HP:0000001", "CLAIM_TYPE": "MECHANISM_OBSERVED", "CLAIM_STATUS": "QUALIFYING", "OBSERVED_OR_INFERRED": "OBSERVED", "OUTCOME_DIRECTION": "DECREASE", "APPLICABILITY": "direct", "ATTRIBUTION": "INDIVIDUAL", "STUDY_DESIGN": "animal", "SOURCE_TYPE": "LITERATURE", "SOURCE_RECORD_ID": "S1", "PASSAGE_LOCATOR": "ABSTRACT:p1:s1"})
    return pd.DataFrame([row]), pd.DataFrame([rationale]), pd.DataFrame([claim])


def test_mechanism_route_promotes_with_effect_and_preserves_originals_and_filter_compatibility():
    drugs, rationales, claims = _inputs()
    enriched, assessments = assess_drug_hpo_rows(drugs, rationales, claims, pd.DataFrame())

    assert enriched.loc[0, "ORIGINAL_DECISION"] == "HOLD"
    assert enriched.loc[0, "ORIGINAL_REASON_CODES"] == "ORIGINAL_REASON"
    assert enriched.loc[0, "DECISION"] == "TEST_FIRST"
    assert enriched.loc[0, "MECHANISM_EVIDENCE_IDS"] == "E1"
    assert enriched.loc[0, "PHENOTYPE_EVIDENCE_IDS"] == ""
    assert enriched.loc[0, "PRIORITY_BAND"] == "C"
    assert len(assessments) == 1
    projected = enriched[list(DRUG_REVIEW_COLUMNS)]
    assert len(filter_drug_candidates(projected, "decision_ready")) == 1
    assert enriched.loc[0, "DECISION"] != "ADVANCE_RESEARCH"


def test_unknown_direction_qualifies_only_phenotype_route_and_unknown_safety_is_explicit():
    drugs, rationales, claims = _inputs()
    claims.loc[0, ["TARGET_SYMBOL", "CLAIM_TYPE", "CLAIM_STATUS", "OUTCOME_DIRECTION", "STUDY_DESIGN", "SOURCE_RECORD_ID", "PMID", "OBSERVED_OR_INFERRED"]] = ["", "BENEFIT_OBSERVED", "QUALIFYING", "BENEFIT", "human case report", "LIT-1", "123", "OBSERVED"]
    claims.loc[0, "TARGET_SYMBOL"] = "T1"
    claims.loc[0, "APPLICABILITY"] = "direct"
    claims.loc[0, "INTERVENTION"] = "Drug"
    claims.loc[0, "MODEL_OR_POPULATION"] = "adult patients"
    drugs.loc[0, "ACTION_TYPE"] = "BINDER"
    drugs.loc[0, "ACTION_DIRECTION"] = "UNKNOWN"

    enriched, _ = assess_drug_hpo_rows(drugs, rationales, claims, pd.DataFrame())

    assert enriched.loc[0, "DECISION"] == "TEST_FIRST"
    assert enriched.loc[0, "PROMOTION_ROUTE"] == "PHENOTYPE_DIRECTED"
    assert enriched.loc[0, "DIRECTIONALITY_FIT"] == "UNKNOWN"
    assert "DIRECTION_UNKNOWN" in enriched.loc[0, "REASON_CODES"]
    assert "PHENOTYPE_DIRECTED" in enriched.loc[0, "REASON_CODES"]
    assert "MECHANISM_CORRECTION" not in enriched.loc[0, "REASON_CODES"]
    assert "SAFETY_UNASSESSED" in enriched.loc[0, "REASON_CODES"]
    assert enriched.loc[0, "PRIORITY_BAND"] == "B"


def test_applicable_hazard_and_conflict_block_both_routes_and_reject_relevant_hazard():
    drugs, rationales, claims = _inputs()
    mechanism = claims.iloc[0].copy()
    hazard = claims.iloc[0].copy()
    hazard.update({"EVIDENCE_ID": "H1", "CLAIM_TYPE": "HAZARD", "CLAIM_STATUS": "QUALIFYING", "OUTCOME_DIRECTION": "HAZARD", "SOURCE_TYPE": "FDA", "APPLICABILITY": "direct", "PASSAGE_LOCATOR": "CONTRAINDICATION:p1:s1"})
    claims = pd.DataFrame([mechanism, hazard])
    enriched, _ = assess_drug_hpo_rows(drugs, rationales, claims, pd.DataFrame())

    assert enriched.loc[0, "DECISION"] == "REJECT"
    assert enriched.loc[0, "PRIORITY_BAND"] == ""
    assert "CITED_CONTEXTUAL_HAZARD" in enriched.loc[0, "REASON_CODES"]


def test_conflicting_qualifying_benefit_and_worsening_hold():
    drugs, rationales, claims = _inputs()
    benefit = claims.iloc[0].copy()
    benefit.update({"CLAIM_TYPE": "BENEFIT_OBSERVED", "OUTCOME_DIRECTION": "BENEFIT", "EVIDENCE_ID": "B1"})
    worse = benefit.copy()
    worse.update({"CLAIM_TYPE": "WORSENING_OBSERVED", "OUTCOME_DIRECTION": "WORSENING", "EVIDENCE_ID": "W1"})
    enriched, _ = assess_drug_hpo_rows(drugs, rationales, pd.DataFrame([benefit, worse]), pd.DataFrame())

    assert enriched.loc[0, "DECISION"] == "HOLD"
    assert "EVIDENCE_CONFLICT" in enriched.loc[0, "REASON_CODES"]
    assert enriched.loc[0, "PRIORITY_BAND"] == ""


def test_bands_rank_before_weaker_evidence_and_preexisting_advance_is_preserved():
    drugs, rationales, claims = _inputs()
    second = drugs.iloc[0].copy()
    second.update({"RATIONALE_ID": "R2", "DRUG_ID": "CHEMBL2", "DRUG_NAME": "Drug 2"})
    drugs = pd.concat([drugs, pd.DataFrame([second])], ignore_index=True)
    rationale = rationales.iloc[0].copy()
    rationale["RATIONALE_ID"] = "R2"
    rationales = pd.concat([rationales, pd.DataFrame([rationale])], ignore_index=True)
    claims.loc[0, ["STUDY_DESIGN", "MODEL_OR_POPULATION"]] = ["animal model", "mouse"]
    human = claims.iloc[0].copy()
    human.update({"CANONICAL_DRUG_ID": "CHEMBL2", "EVIDENCE_ID": "H1", "STUDY_DESIGN": "human clinical intervention", "SOURCE_RECORD_ID": "S2", "PMID": "100"})
    human2 = human.copy()
    human2.update({"EVIDENCE_ID": "H2", "SOURCE_RECORD_ID": "S3", "PMID": "101"})
    claims = pd.concat([claims, pd.DataFrame([human, human2])], ignore_index=True)
    drugs.loc[0, "DECISION"] = "ADVANCE_RESEARCH"

    enriched, _ = assess_drug_hpo_rows(drugs, rationales, claims, pd.DataFrame())

    assert enriched.loc[0, "DECISION"] == "ADVANCE_RESEARCH"
    assert enriched.loc[1, "DECISION"] == "TEST_FIRST"
    assert enriched.loc[1, "PRIORITY_BAND"] == "A"
    assert enriched.loc[1, "PRIORITY_RANK"] < enriched.loc[0, "PRIORITY_RANK"]
