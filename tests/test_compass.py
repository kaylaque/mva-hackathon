"""Behavioural tests for Compass: eligibility, decisions, resumption and reporting.

Every public response is synthetic and injected; no test contacts the network, and no
patient, sample or repository table is read.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import re
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).parents[1] / "src"))

import pandas as pd
import pytest

from mva_hackathon import compass
from mva_hackathon import compass_sources as sources
from mva_hackathon.compass import (
    DECISION_ADVANCE,
    DECISION_HOLD,
    DECISION_REJECT,
    DECISION_TEST_FIRST,
    ELIGIBLE,
    STOP_NOT_SIGNIFICANT,
    MECHANISM_PRIOR_PANEL,
    build_eligibility,
    build_rationales,
    build_test_plans,
    load_postanalysis,
    run_compass,
)


class Response:
    def __init__(self, payload: object = None, *, status_code: int = 200, text: str | None = None) -> None:
        self.status_code = status_code
        self.headers: dict[str, str] = {}
        self.text = text if text is not None else ("" if payload is None else json.dumps(payload))
        self._payload = payload

    def json(self) -> object:
        if self._payload is None:
            raise ValueError("no json body")
        return self._payload


class FakeTransport:
    """Injected public requester covering every adapter used by Compass."""

    def __init__(self, failures: set[str] | None = None) -> None:
        self.calls: list[str] = []
        self.failures = failures or set()

    def __call__(self, method: str, url: str, *, params: object = None) -> Response:
        self.calls.append(url)
        for marker in self.failures:
            if marker in url:
                raise TimeoutError(f"synthetic outage for {marker}")
        query = (params or {})
        if "/data/query/" in url:
            return Response({"schemaClass": "Pathway", "stId": url.rsplit("/", 1)[-1], "displayName": "Synthetic event", "speciesName": "Homo sapiens"})
        if "/data/participants/" in url:
            return Response([{"peDbId": 1, "refEntities": [{"schemaClass": "ReferenceGeneProduct", "identifier": "O60566", "displayName": "UniProt:O60566 BUB1B"}]}])
        if "rest.kegg.jp/get/" in url:
            return Response(text="ENTRY       hsa04110                    Pathway\nNAME        Cell cycle - Homo sapiens (human)\nGENE        701  BUB1B; BUB1 mitotic checkpoint kinase B [KO:K06637]\n            7157  TP53; cellular tumor antigen p53 [KO:K04451]\nCOMPOUND    C00002  ATP\n///\n")
        if "rest.kegg.jp/link/" in url:
            return Response(text="path:hsa04110\thsa:701\npath:hsa04110\thsa:7157\n")
        if "ontology/go/terms" in url:
            return Response({"numberOfHits": 1, "results": [{"id": "GO:0007059", "aspect": "biological_process", "name": "chromosome segregation"}]})
        if "annotation/search" in url:
            product = (params or {}).get("geneProductId")
            go_id = (params or {}).get("goId")
            if go_id:
                return Response({"results": [{"id": "a1", "geneProductId": "UniProtKB:O60566", "symbol": "BUB1B", "goEvidence": "IDA", "qualifier": "involved_in"}], "pageInfo": {"total": 1, "resultsPerPage": 100}})
            return Response({"results": [{"id": "a1", "goId": "GO:0007059", "goEvidence": "IDA", "qualifier": "involved_in"}] if product else [], "pageInfo": {"total": 1 if product else 0, "resultsPerPage": 100}})
        if "rest.genenames.org" in url:
            symbol = url.rsplit("/", 1)[-1].upper()
            ensembl = "ENSG" + str(abs(int(hashlib.sha256(symbol.encode()).hexdigest()[:12], 16)) % 10**11).zfill(11)
            return Response({"response": {"numFound": 1, "docs": [{"symbol": symbol, "alias_symbol": [], "prev_symbol": [], "ensembl_gene_id": ensembl, "uniprot_ids": ["O60566"], "hgnc_id": "HGNC:1"}]}})
        if "rest.uniprot.org" in url:
            symbol = re.search(r"gene_exact:([A-Z0-9]+)", str((params or {}).get("query", ""))).group(1)
            return Response({"results": [{"primaryAccession": "O60566", "genes": [{"geneName": {"value": symbol}}], "comments": [{"commentType": "FUNCTION", "texts": [{"value": f"Curated function for {symbol}."}]}]}]})
        if "esearch.fcgi" in url:
            return Response({"esearchresult": {"count": "1", "idlist": ["1"]}})
        if "esummary.fcgi" in url:
            return Response({"result": {"1": {"title": "Synthetic title", "source": "J", "pubdate": "2024", "articleids": []}}})
        if "target.json" in url:
            return Response({"targets": [{"target_chembl_id": "CHEMBL1", "pref_name": "BUB1B", "target_type": "SINGLE PROTEIN", "organism": "Homo sapiens", "target_components": [{"accession": "O60566"}]}], "page_meta": {"next": None}})
        if "mechanism.json" in url:
            return Response({"mechanisms": [{"mec_id": "M1", "molecule_chembl_id": "CHEMBL_OK", "action_type": "ACTIVATOR", "mechanism_of_action": "activator"}, {"mec_id": "M2", "molecule_chembl_id": "CHEMBL_BAD", "action_type": "INHIBITOR", "mechanism_of_action": "inhibitor"}], "page_meta": {"next": None}})
        if "molecule/CHEMBL_OK.json" in url:
            return Response({"molecule_chembl_id": "CHEMBL_OK", "pref_name": "Drug Ok"})
        if "molecule/CHEMBL_BAD.json" in url:
            return Response({"molecule_chembl_id": "CHEMBL_BAD", "pref_name": "Drug Bad"})
        if "activity.json" in url:
            return Response({"activities": [], "page_meta": {"next": None}})
        if "opentargets" in url:
            return Response({"data": {"target": {"id": "ENSG00000156970", "approvedSymbol": "BUB1B", "drugAndClinicalCandidates": {"rows": []}}}})
        raise AssertionError(f"unexpected request: {method} {url}")


def _glue(evidence_id: str, ontology: str, term_id: str, match: str, g_rows: list[dict], d_rows: list[dict], g_state: str, d_state: str) -> dict[str, str]:
    return {
        "EVIDENCE_ID": evidence_id,
        "ONTOLOGY": ontology,
        "TERM_ID": term_id,
        "MATCH_STATUS": match,
        "GPROFILER_ROWS_JSON": json.dumps(g_rows),
        "DAVID_ROWS_JSON": json.dumps(d_rows),
        "GPROFILER_SIGNIFICANCE": g_state,
        "DAVID_SIGNIFICANCE": d_state,
        "EVIDENCE_TIER": "context_only",
        "GPROFILER_CONTRIBUTOR_STATUS": "absent" if g_state == "not_comparable" else "complete",
        "DAVID_CONTRIBUTOR_STATUS": "absent" if d_state == "not_comparable" else "complete",
        "CONTRIBUTOR_STATUS": "complete",
        "MVA_CONTEXT_JSON": "[]",
    }


def write_postanalysis(root: Path) -> Path:
    """A synthetic saved post-analysis with two eligible providers and one context annotation."""
    root.mkdir(parents=True, exist_ok=True)
    evidence = pd.DataFrame([
        _glue("EV-1", "REACTOME", "R-HSA-1", "gprofiler_only", [{"TERM_NAME": "Cell Cycle Checkpoints", "BH_Q": 0.01}], [], "significant", "not_comparable"),
        _glue("EV-2", "KEGG", "hsa04110", "david_only", [], [{"Term": "Cell cycle", "Benjamini": 0.02}], "not_comparable", "significant"),
        _glue("EV-3", "GO", "GO:0007059", "shared", [{"TERM_NAME": "chromosome segregation", "BH_Q": 0.03}], [{"Term": "chromosome segregation", "Benjamini": 0.03}], "significant", "significant"),
        _glue("EV-4", "HP", "HP:0000001", "shared", [{"TERM_NAME": "mosaic variegated aneuploidy", "BH_Q": 0.01}], [{"Term": "mosaic variegated aneuploidy", "Benjamini": 0.01}], "significant", "significant"),
        _glue("EV-5", "REACTOME", "R-HSA-2", "gprofiler_only", [{"TERM_NAME": "Cell cycle", "BH_Q": 0.4}], [], "not_significant", "not_comparable"),
        _glue("EV-6", "REACTOME", "free text keyword", "unmapped_term", [{"TERM_ID": "free text keyword", "TERM_NAME": "unmappable", "BH_Q": 0.01}], [], "significant", "not_comparable"),
    ])
    evidence.to_csv(root / "integrated_pathway_evidence.csv", index=False)
    pd.DataFrame([
        {"EVIDENCE_ID": "EV-1", "ONTOLOGY": "REACTOME", "TERM_ID": "R-HSA-1", "PROVIDER": "g:Profiler", "GENE_SYMBOL": "BUB1B", "MEMBERSHIP_STATUS": "observed"},
        {"EVIDENCE_ID": "EV-2", "ONTOLOGY": "KEGG", "TERM_ID": "hsa04110", "PROVIDER": "DAVID", "GENE_SYMBOL": "BUB1B", "MEMBERSHIP_STATUS": "observed"},
        {"EVIDENCE_ID": "EV-3", "ONTOLOGY": "GO", "TERM_ID": "GO:0007059", "PROVIDER": "g:Profiler", "GENE_SYMBOL": "MAD2L1", "MEMBERSHIP_STATUS": "observed"},
        {"EVIDENCE_ID": "EV-5", "ONTOLOGY": "REACTOME", "TERM_ID": "R-HSA-2", "PROVIDER": "g:Profiler", "GENE_SYMBOL": "TRIP13", "MEMBERSHIP_STATUS": "observed"},
    ]).to_csv(root / "gene_pathway_support.csv", index=False)
    (root / "postanalysis_provenance.json").write_text(json.dumps({"baseline_projections": "abc", "import_version": "def"}), encoding="utf-8")
    return root


def write_network(root: Path) -> Path:
    root.mkdir(parents=True, exist_ok=True)
    pd.DataFrame([
        {"PREFERRED_NAME_A": "BUB1B", "PREFERRED_NAME_B": "MAD2L1"},
        {"PREFERRED_NAME_A": "BUB1B", "PREFERRED_NAME_B": "TRIP13"},
    ]).to_csv(root / "04_phase2_string_edges.csv", index=False)
    pd.DataFrame([{"HGNC_SYMBOL": "BUB1B"}, {"HGNC_SYMBOL": "MAD2L1"}, {"HGNC_SYMBOL": "TRIP13"}, {"HGNC_SYMBOL": "GENE_ISOLATED"}]).to_csv(root / "07_phase2_network_gene_context.csv", index=False)
    return root


def write_reviewed(path: Path, *, hazard_symbol: str | None = "AURKA", conflict_symbol: str | None = None) -> Path:
    columns = list(compass.MECHANISM_EVIDENCE_COLUMNS)

    def row(**values: str) -> dict[str, str]:
        base = {column: "" for column in columns}
        base.update({
            "CLAIM_TYPE": "PROTEIN_ROLE",
            "SOURCE": "manual review",
            "REVIEW_STATUS": "REVIEWED",
            "SOURCE_STATUS": "OK",
            "SOURCE_SECTION": "Abstract results",
            "EVIDENCE_CLASS": "HUMAN_CELL_EXPERIMENT",
            "SUPPORT_DIRECTION": "SUPPORTS",
        })
        base.update(values)
        base["EVIDENCE_ID"] = values.get("EVIDENCE_ID") or f"ME-{abs(hash(json.dumps(values, sort_keys=True))) % 10**10}"
        return base

    rows = [
        row(EVIDENCE_ID="ME-1", TARGET_SYMBOL="BUB1B", CLAIM_TYPE="DISEASE_MECHANISM", SUPPORT_DIRECTION="REDUCED_FUNCTION", CLAIM="Biallelic BUB1B loss reduces checkpoint function.", PMID="1", CITATION_URL="https://pubmed.ncbi.nlm.nih.gov/1/", MODEL_SYSTEM="human patient fibroblasts", SOURCE_SECTION="Case description"),
        row(EVIDENCE_ID="ME-2", TARGET_SYMBOL="BUB1B", CLAIM_TYPE="PERTURBATION_EFFECT", SUPPORT_DIRECTION="SUPPORTS", CLAIM="Restoring BUB1B expression improves checkpoint arrest.", PMID="2", CITATION_URL="https://pubmed.ncbi.nlm.nih.gov/2/", MODEL_SYSTEM="human RPE1 cells", PERTURBATION="BUB1B re-expression", OBSERVED_EFFECT="improved checkpoint arrest"),
    ]
    if hazard_symbol:
        rows.append(row(EVIDENCE_ID="ME-3", TARGET_SYMBOL=hazard_symbol, CLAIM_TYPE="HAZARD", SUPPORT_DIRECTION="WORSENS", CLAIM="Inhibition worsened mitotic fidelity in the cited model.", PMID="3", CITATION_URL="https://pubmed.ncbi.nlm.nih.gov/3/", MODEL_SYSTEM="human cells"))
    if conflict_symbol:
        rows.append(row(EVIDENCE_ID="ME-4", TARGET_SYMBOL=conflict_symbol, CLAIM_TYPE="DISEASE_MECHANISM", SUPPORT_DIRECTION="REDUCED_FUNCTION", CLAIM="Reduced function reported.", PMID="4", CITATION_URL="https://pubmed.ncbi.nlm.nih.gov/4/", MODEL_SYSTEM="human cells"))
        rows.append(row(EVIDENCE_ID="ME-5", TARGET_SYMBOL=conflict_symbol, CLAIM_TYPE="DISEASE_MECHANISM", SUPPORT_DIRECTION="INCREASED_FUNCTION", CLAIM="Increased function reported elsewhere.", PMID="5", CITATION_URL="https://pubmed.ncbi.nlm.nih.gov/5/", MODEL_SYSTEM="human cells"))
    pd.DataFrame(rows, columns=columns).to_csv(path, index=False)
    return path


@pytest.fixture
def artifacts(tmp_path: Path) -> dict[str, Path]:
    return {
        "postanalysis": write_postanalysis(tmp_path / "postanalysis"),
        "network": write_network(tmp_path / "network"),
        "cache": tmp_path / "cache",
        "reviewed": write_reviewed(tmp_path / "reviewed.csv"),
        "output": tmp_path / "compass",
    }


# --------------------------------------------------------------------------- eligibility


def test_eligibility_judges_each_provider_separately(tmp_path: Path) -> None:
    postanalysis = load_postanalysis(write_postanalysis(tmp_path / "post"))
    eligibility = build_eligibility(postanalysis["evidence"], postanalysis["support"])
    by_id = {(row["ONTOLOGY"], row["TERM_ID"], row["PROVIDER"]): row for _, row in eligibility.iterrows()}

    # A single significant provider enters the gate without the other provider validating it.
    assert by_id[("REACTOME", "R-HSA-1", "g:Profiler")]["DECISION"] == ELIGIBLE
    assert by_id[("KEGG", "hsa04110", "DAVID")]["DECISION"] == ELIGIBLE
    assert by_id[("GO", "GO:0007059", "g:Profiler")]["DECISION"] == ELIGIBLE
    assert by_id[("KEGG", "hsa04110", "DAVID")]["SIGNIFICANCE"] == "significant"

    # A nonsignificant provider does not veto the valid provider-specific research.
    assert by_id[("REACTOME", "R-HSA-2", "g:Profiler")]["DECISION"] == STOP_NOT_SIGNIFICANT
    assert by_id[("REACTOME", "R-HSA-1", "g:Profiler")]["DECISION"] == ELIGIBLE

    # An observational phenotype term stays annotation context even when significant.
    hp = by_id[("HP", "HP:0000001", "g:Profiler")]
    assert hp["ENTITY_TYPE"] == "ANNOTATION_CONTEXT"
    assert hp["DECISION"] == STOP_NOT_SIGNIFICANT
    assert "ANNOTATION_CONTEXT_NOT_A_PATHWAY" in hp["REASON_CODES"]

    # A significant unmapped annotation keeps a distinct identity and a ledger entry.
    unmapped = eligibility.loc[eligibility["DECISION"].eq("STOP_UNMAPPED")]
    assert len(unmapped) >= 1
    assert unmapped["EVIDENCE_ID"].is_unique

def test_david_relaxed_admits_only_benjamini_010_observed_panel_pathways() -> None:
    evidence = pd.DataFrame([
        _glue("EV-04", "REACTOME", "R-HSA-04", "david_only", [], [{"Term": "validated", "Benjamini": 0.04}], "not_comparable", "significant"),
        _glue("EV-08", "REACTOME", "R-HSA-08", "david_only", [], [{"Term": "relaxed", "Benjamini": 0.08}], "not_comparable", "significant"),
        _glue("EV-11", "REACTOME", "R-HSA-11", "david_only", [], [{"Term": "stopped", "Benjamini": 0.11}], "not_comparable", "not_significant"),
        _glue("EV-GP", "REACTOME", "R-HSA-GP", "gprofiler_only", [{"TERM_NAME": "g:Profiler", "BH_Q": 0.08}], [], "not_significant", "not_comparable"),
    ])
    support = pd.DataFrame([
        {"EVIDENCE_ID": evidence_id, "ONTOLOGY": "REACTOME", "TERM_ID": term_id, "PROVIDER": provider, "GENE_SYMBOL": "BUB1B", "MEMBERSHIP_STATUS": "observed"}
        for evidence_id, term_id, provider in (
            ("EV-04", "R-HSA-04", "DAVID"),
            ("EV-08", "R-HSA-08", "DAVID"),
            ("EV-11", "R-HSA-11", "DAVID"),
            ("EV-GP", "R-HSA-GP", "g:Profiler"),
        )
    ])

    validated = build_eligibility(evidence, support)
    relaxed = build_eligibility(
        evidence, support, selection_profile="david_relaxed", david_benjamini_threshold=0.10
    ).set_index("EVIDENCE_IDS")

    assert set(validated.loc[validated["DECISION"].eq(ELIGIBLE), "EVIDENCE_IDS"]) == {"EV-04"}
    assert relaxed.loc["EV-04", "DECISION"] == ELIGIBLE
    assert relaxed.loc["EV-08", "DECISION"] == "ELIGIBLE_EXPLORATORY"
    assert relaxed.loc["EV-08", "SIGNIFICANCE"] == "exploratory_relaxed"
    assert "DAVID_BENJAMINI_0_10_EXPLORATORY" in relaxed.loc["EV-08", "REASON_CODES"]
    assert relaxed.loc["EV-11", "DECISION"] == STOP_NOT_SIGNIFICANT
    assert relaxed.loc["EV-GP", "DECISION"] == STOP_NOT_SIGNIFICANT

def test_ninety_percent_david_profile_keeps_provider_and_cutoff_boundaries() -> None:
    evidence = pd.DataFrame([
        _glue("EV-80", "REACTOME", "R-HSA-80", "david_only", [], [{"Term": "exploratory", "Benjamini": 0.80}], "not_comparable", "significant"),
        _glue("EV-92", "REACTOME", "R-HSA-92", "david_only", [], [{"Term": "outside", "Benjamini": 0.92}], "not_comparable", "not_significant"),
        _glue("EV-G", "REACTOME", "R-HSA-G", "gprofiler_only", [{"TERM_NAME": "unvalidated", "BH_Q": 0.80}], [], "not_significant", "not_comparable"),
    ])
    support = pd.DataFrame([
        {"EVIDENCE_ID": evidence_id, "ONTOLOGY": "REACTOME", "TERM_ID": term_id, "PROVIDER": provider, "GENE_SYMBOL": "BUB1B", "MEMBERSHIP_STATUS": "observed"}
        for evidence_id, term_id, provider in (
            ("EV-80", "R-HSA-80", "DAVID"),
            ("EV-92", "R-HSA-92", "DAVID"),
            ("EV-G", "R-HSA-G", "g:Profiler"),
        )
    ])
    rows = build_eligibility(
        evidence, support, selection_profile="david_relaxed", david_benjamini_threshold=0.90
    ).set_index("EVIDENCE_IDS")
    assert rows.loc["EV-80", "DECISION"] == "ELIGIBLE_EXPLORATORY"
    assert "DAVID_BENJAMINI_0_90_EXPLORATORY" in rows.loc["EV-80", "REASON_CODES"]
    assert rows.loc["EV-92", "DECISION"] == STOP_NOT_SIGNIFICANT
    assert rows.loc["EV-G", "DECISION"] == STOP_NOT_SIGNIFICANT



@pytest.mark.parametrize("category,expected", [
    ("GOTERM_CC_DIRECT", STOP_NOT_SIGNIFICANT),
    ("GOTERM_MF_DIRECT", STOP_NOT_SIGNIFICANT),
    ("GOTERM_BP_DIRECT", "ELIGIBLE_EXPLORATORY"),
])
def test_go_annotation_aspect_cannot_nominate_non_process_drugs(category: str, expected: str) -> None:
    evidence = pd.DataFrame([
        _glue("EV-GO", "GO", "GO:0007059", "david_only", [],
              [{"Category": category, "Term": "synthetic annotation", "Benjamini": 0.80}],
              "not_comparable", "significant"),
    ])
    support = pd.DataFrame([{
        "EVIDENCE_ID": "EV-GO", "ONTOLOGY": "GO", "TERM_ID": "GO:0007059",
        "PROVIDER": "DAVID", "GENE_SYMBOL": "BUB1B", "MEMBERSHIP_STATUS": "observed",
    }])
    row = build_eligibility(
        evidence, support, selection_profile="david_relaxed", david_benjamini_threshold=0.90
    ).iloc[0]
    assert row["DECISION"] == expected
    if expected == STOP_NOT_SIGNIFICANT:
        assert row["ENTITY_TYPE"] == compass.ENTITY_CONTEXT
        assert "ANNOTATION_CONTEXT_NOT_A_PATHWAY" in row["REASON_CODES"]


@pytest.mark.parametrize("threshold", (0, -0.01, 1.01, float("inf")))
def test_invalid_david_relaxed_threshold_is_rejected_before_row_processing(threshold: float) -> None:
    with pytest.raises(ValueError, match="threshold"):
        build_eligibility(pd.DataFrame(), pd.DataFrame(), selection_profile="david_relaxed", david_benjamini_threshold=threshold)


def test_drug_candidate_filters_preserve_defined_subsets() -> None:
    frame = pd.DataFrame([
        {"DIRECTIONALITY_FIT": "COMPATIBLE_DIRECTION", "HAZARD_STATUS": "unknown", "DECISION": "TEST_FIRST"},
        {"DIRECTIONALITY_FIT": "MECHANISTIC_HAZARD", "HAZARD_STATUS": "cited_hazard", "DECISION": "REJECT"},
    ])
    assert len(compass.filter_drug_candidates(frame, "all")) == 2
    assert len(compass.filter_drug_candidates(frame, "direction_compatible")) == 1
    assert len(compass.filter_drug_candidates(frame, "no_cited_hazard")) == 1
    assert len(compass.filter_drug_candidates(frame, "decision_ready")) == 1

def test_exploratory_drugs_do_not_override_hazards_or_held_rationales(tmp_path: Path) -> None:
    rationales = pd.DataFrame([
        {"RATIONALE_ID": "RAT-H", "TARGET_SYMBOL": "BUB1B", "DECISION": DECISION_REJECT,
         "REASON_CODES": "CITED_CONTEXTUAL_HAZARD", "COUNTEREVIDENCE_IDS": "ME-H",
         "EXPLORATORY_DAVID_OBSERVED_TERM_IDS": "hsa00080",
         "DESIRED_TARGET_CHANGE": "RESTORE_NORMAL", "DISEASE_STATE_DIRECTION": "REDUCED_FUNCTION"},
        {"RATIONALE_ID": "RAT-U", "TARGET_SYMBOL": "BUB1B", "DECISION": DECISION_HOLD,
         "REASON_CODES": "SAFETY_UNRECORDED", "COUNTEREVIDENCE_IDS": "",
         "EXPLORATORY_DAVID_OBSERVED_TERM_IDS": "hsa00080",
         "DESIRED_TARGET_CHANGE": "RESTORE_NORMAL", "DISEASE_STATE_DIRECTION": "REDUCED_FUNCTION"},
    ])
    drugs, _, _ = compass.fetch_drug_review(
        rationales, tmp_path / "cache", selection_profile="david_relaxed",
        david_benjamini_threshold=0.90, drug_scope="relaxed_pathway", transport=FakeTransport(),
    )
    assert not drugs.empty
    assert set(drugs.loc[drugs.RATIONALE_ID.eq("RAT-H"), "DECISION"]) == {DECISION_REJECT}
    assert set(drugs.loc[drugs.RATIONALE_ID.eq("RAT-H"), "HAZARD_STATUS"]) == {"cited_hazard"}
    assert "RAT-H" not in set(compass.filter_drug_candidates(drugs, "no_cited_hazard").RATIONALE_ID)
    assert set(drugs.loc[drugs.RATIONALE_ID.eq("RAT-U"), "DECISION"]).issubset({DECISION_HOLD, DECISION_REJECT})
    assert not compass.filter_drug_candidates(drugs, "decision_ready").shape[0]
    assert "NONE" not in set(drugs.loc[drugs.RATIONALE_ID.eq("RAT-U"), "EVIDENCE_GAPS"])


def test_incomplete_drug_identity_lookup_marks_compass_source_incomplete(tmp_path: Path, monkeypatch) -> None:
    rationales = pd.DataFrame([{
        "RATIONALE_ID": "R1",
        "TARGET_SYMBOL": "BUB1B",
        "DECISION": DECISION_HOLD,
        "REASON_CODES": "",
        "COUNTEREVIDENCE_IDS": "",
        "EXPLORATORY_DAVID_OBSERVED_TERM_IDS": "hsa05170",
        "DESIRED_TARGET_CHANGE": "RESTORE_NORMAL",
        "DISEASE_STATE_DIRECTION": "REDUCED_FUNCTION",
    }])
    monkeypatch.setattr(
        sources,
        "fetch_drug_evidence",
        lambda *_args, **_kwargs: {
            "status": sources.STATUS_OK,
            "records": [],
            "complete": False,
            "source": "CHEMBL+OPEN_TARGETS",
            "request_id": "identity",
            "response_sha256": "",
            "retrieved_at": "",
            "reason": "molecule name unavailable",
        },
    )

    _, statuses, failures = compass.fetch_drug_review(
        rationales, tmp_path / "cache", selection_profile="david_relaxed",
        david_benjamini_threshold=0.90, drug_scope="relaxed_pathway",
    )

    assert statuses[0]["COMPLETE"] is False
    assert failures == ["Pharmacology:BUB1B"]




def test_eligibility_does_not_require_a_mechanism_panel_member(tmp_path: Path) -> None:
    """An observed contributor outside the MVA panel no longer blocks an eligible term."""
    postanalysis = load_postanalysis(write_postanalysis(tmp_path / "post"))
    support = postanalysis["support"].copy()
    # The significant term now has a non-panel contributor only.
    support.loc[0, "GENE_SYMBOL"] = "RPL13"
    eligibility = build_eligibility(postanalysis["evidence"], support)
    row = eligibility.loc[eligibility["EVIDENCE_ID"].eq("EV-1::gprofiler")].iloc[0]

    assert row["SIGNIFICANCE"] == "significant"
    assert row["DECISION"] == ELIGIBLE
    assert "PROVIDER_SIGNIFICANCE_AND_OBSERVED_MEMBER" in row["REASON_CODES"]
    assert row["MEMBERSHIP_STATUS"] == "observed"


def test_exploratory_term_without_panel_member_is_eligible() -> None:
    """The relaxed profile also admits on statistics plus observed membership alone."""
    evidence = pd.DataFrame([
        _glue("EV-NP", "REACTOME", "R-HSA-NP", "david_only", [],
              [{"Term": "non-panel pathway", "Benjamini": 0.08}],
              "not_comparable", "significant"),
    ])
    support = pd.DataFrame([{
        "EVIDENCE_ID": "EV-NP", "ONTOLOGY": "REACTOME", "TERM_ID": "R-HSA-NP",
        "PROVIDER": "DAVID", "GENE_SYMBOL": "RPL13", "MEMBERSHIP_STATUS": "observed",
    }])
    row = build_eligibility(
        evidence, support, selection_profile="david_relaxed", david_benjamini_threshold=0.10
    ).iloc[0]

    assert row["DECISION"] == "ELIGIBLE_EXPLORATORY"
    assert "DAVID_BENJAMINI_0_10_EXPLORATORY" in row["REASON_CODES"]


# --------------------------------------------------------------------------- rationales


def test_rationales_separate_lanes_and_never_invent_direction() -> None:
    eligibility = pd.DataFrame([{
        "EVIDENCE_ID": "EV-1::gprofiler", "PROVIDER": "g:Profiler", "ONTOLOGY": "REACTOME", "TERM_ID": "R-HSA-1",
        "TERM_NAME": "Cell Cycle Checkpoints", "ENTITY_TYPE": "PATHWAY", "SIGNIFICANCE": "significant",
        "CONTRIBUTOR_STATUS": "complete", "MEMBERSHIP_STATUS": "observed", "DECISION": ELIGIBLE,
        "REASON_CODES": "", "EVIDENCE_IDS": "EV-1",
    }])
    memberships = pd.DataFrame([
        {"EVIDENCE_ID": "EV-1::gprofiler", "PROVIDER": "g:Profiler", "ONTOLOGY": "REACTOME", "TERM_ID": "R-HSA-1", "TARGET_SYMBOL": "BUB1B", "MEMBERSHIP_KIND": "OBSERVED_INTERSECTION", "MEMBERSHIP_STATUS": "observed", "SOURCE_RECORD_ID": "x"},
        {"EVIDENCE_ID": "EV-1::gprofiler", "PROVIDER": "Reactome", "ONTOLOGY": "REACTOME", "TERM_ID": "R-HSA-1", "TARGET_SYMBOL": "TTK", "MEMBERSHIP_KIND": "EXPANDED_DATABASE_MEMBER", "MEMBERSHIP_STATUS": "complete", "SOURCE_RECORD_ID": "y"},
    ])
    protein_roles = pd.DataFrame([{"TARGET_SYMBOL": "TTK", "UNIPROT_ACCESSION": "P33981", "PROTEIN_ROLE": "curated_function", "ROLE_CONTEXT": "kinetochore kinase", "SOURCE_STATUS": "OK", "EVIDENCE_IDS": "ME-x"}])
    network = pd.DataFrame([{"TARGET_SYMBOL": "BUB1B", "MEMBERSHIP_KIND": "NETWORK_COMPONENT", "COMPONENT_ID": "1", "MEMBERSHIP_STATUS": "LEGACY_UNVERIFIED", "SOURCE_RECORD_ID": "c1"}])

    rationales = build_rationales(eligibility, memberships, protein_roles, pd.DataFrame(columns=compass.MECHANISM_EVIDENCE_COLUMNS), network)
    by_symbol = {row["TARGET_SYMBOL"]: row for _, row in rationales.iterrows()}

    supported = by_symbol["BUB1B"]
    assert compass.LABEL_PATHWAY_SUPPORTED in supported["HYPOTHESIS_SOURCE"]
    assert compass.LABEL_MECHANISM_PRIOR in supported["HYPOTHESIS_SOURCE"]
    assert compass.LABEL_NETWORK_CONTEXT in supported["HYPOTHESIS_SOURCE"]
    assert supported["PATHWAY_ENRICHMENT_SUPPORT"] == compass.LABEL_PATHWAY_SUPPORTED

    # A database member is not an enrichment contributor and carries no pathway claim.
    expanded_only = by_symbol["TTK"]
    assert expanded_only["EXPANDED_MEMBERSHIP"] is True or expanded_only["EXPANDED_MEMBERSHIP"] == True  # noqa: E712
    assert expanded_only["PATHWAY_ENRICHMENT_SUPPORT"] == "NOT_IN_OBSERVED_INTERSECTION"

    # Protein role alone leaves the desired change and disease direction unknown.
    assert expanded_only["DESIRED_TARGET_CHANGE"] == "UNKNOWN"
    assert expanded_only["DISEASE_STATE_DIRECTION"] == "UNKNOWN"
    assert expanded_only["DECISION"] == DECISION_HOLD
    assert "AWAITING_ROLE_EVIDENCE" in expanded_only["REASON_CODES"]

    # A mechanism-prior panel member with no significant pathway is still researched.
    panel_only = [row for symbol, row in by_symbol.items() if symbol in MECHANISM_PRIOR_PANEL and symbol not in {"BUB1B"}]
    assert panel_only
    assert all("NO_PATHWAY_ENRICHMENT_SUPPORT" in row["REASON_CODES"] for row in panel_only)


def _verified_identity(*symbols: str) -> pd.DataFrame:
    return pd.DataFrame([
        {"TARGET_SYMBOL": symbol, "UNIPROT_ACCESSION": "O60566", "PROTEIN_ROLE": "curated_function", "ROLE_CONTEXT": "curated role", "SOURCE_STATUS": "OK", "EVIDENCE_IDS": "ME-x"}
        for symbol in symbols
    ])


def test_reviewed_evidence_drives_decisions_and_hazards(tmp_path: Path) -> None:
    postanalysis = load_postanalysis(write_postanalysis(tmp_path / "post"))
    eligibility = build_eligibility(postanalysis["evidence"], postanalysis["support"])
    reviewed, _ = compass.load_reviewed_evidence(write_reviewed(tmp_path / "reviewed.csv", conflict_symbol="TRIP13"))
    identity = _verified_identity("BUB1B", "AURKA", "TRIP13")
    rationales = build_rationales(eligibility, pd.DataFrame(columns=compass.MEMBERSHIP_COLUMNS), identity, reviewed, pd.DataFrame())
    by_symbol = {row["TARGET_SYMBOL"]: row for _, row in rationales.iterrows()}

    bub1b = by_symbol["BUB1B"]
    assert bub1b["DISEASE_STATE_DIRECTION"] == "REDUCED_FUNCTION"
    assert bub1b["DESIRED_TARGET_CHANGE"] == "RESTORE_NORMAL"
    assert bub1b["DECISION"] == DECISION_TEST_FIRST

    assert by_symbol["AURKA"]["DECISION"] == DECISION_REJECT
    assert "CITED_CONTEXTUAL_HAZARD" in by_symbol["AURKA"]["REASON_CODES"]
    assert by_symbol["TRIP13"]["DECISION"] == DECISION_HOLD
    assert "EVIDENCE_CONFLICT" in by_symbol["TRIP13"]["REASON_CODES"]


def test_reviewed_rescue_without_safety_holds_while_complete_rescue_advances(tmp_path: Path) -> None:
    postanalysis = load_postanalysis(write_postanalysis(tmp_path / "post"))
    eligibility = build_eligibility(postanalysis["evidence"], postanalysis["support"])
    reviewed, _ = compass.load_reviewed_evidence(write_reviewed(tmp_path / "reviewed.csv", hazard_symbol=None))
    base = reviewed.copy()

    identity = _verified_identity("BUB1B", "TRIP13")
    without_safety = build_rationales(eligibility, pd.DataFrame(columns=compass.MEMBERSHIP_COLUMNS), identity, base, pd.DataFrame())
    rescue = base.loc[0].copy()
    rescue.update({"EVIDENCE_ID": "ME-9", "TARGET_SYMBOL": "BUB1B", "CLAIM_TYPE": "RESCUE", "SUPPORT_DIRECTION": "SUPPORTS", "CLAIM": "Rescue restored checkpoint function.", "PMID": "9", "CITATION_URL": "https://pubmed.ncbi.nlm.nih.gov/9/", "MODEL_SYSTEM": "human RPE1 cells"})
    with_rescue = pd.concat([base, pd.DataFrame([rescue])], ignore_index=True)

    assert build_rationales(eligibility, pd.DataFrame(columns=compass.MEMBERSHIP_COLUMNS), identity, with_rescue, pd.DataFrame()).set_index("TARGET_SYMBOL").loc["BUB1B", "DECISION"] == DECISION_HOLD

    safety = rescue.copy()
    safety.update({"EVIDENCE_ID": "ME-10", "CLAIM_TYPE": "SAFETY", "SUPPORT_DIRECTION": "SUPPORTS", "CLAIM": "Safety readouts recorded.", "PMID": "10", "CITATION_URL": "https://pubmed.ncbi.nlm.nih.gov/10/"})
    completed = pd.concat([with_rescue, pd.DataFrame([safety])], ignore_index=True)
    final = build_rationales(eligibility, pd.DataFrame(columns=compass.MEMBERSHIP_COLUMNS), identity, completed, pd.DataFrame())
    row = final.set_index("TARGET_SYMBOL").loc["BUB1B"]

    assert row["DECISION"] == DECISION_ADVANCE
    assert "REVIEWED_RESCUE_WITH_MODEL_AND_SAFETY_READOUT" in row["REASON_CODES"]
    assert without_safety.set_index("TARGET_SYMBOL").loc["BUB1B", "DECISION"] == DECISION_TEST_FIRST


def test_unresolved_rationales_get_a_resolution_plan_not_an_experiment() -> None:
    rationales = pd.DataFrame([
        {"RATIONALE_ID": "RAT-1", "TARGET_SYMBOL": "BUB1B", "DECISION": DECISION_HOLD, "REASON_CODES": "AWAITING_ROLE_EVIDENCE", "EVIDENCE_IDS": "", "DESIRED_TARGET_CHANGE": "UNKNOWN", "DISEASE_STATE_DIRECTION": "UNKNOWN"},
        {"RATIONALE_ID": "RAT-2", "TARGET_SYMBOL": "TTK", "DECISION": DECISION_TEST_FIRST, "REASON_CODES": "CONTEXTUAL_HYPOTHESIS_WITH_PERTURBATION", "EVIDENCE_IDS": "ME-1", "DESIRED_TARGET_CHANGE": "RESTORE_NORMAL", "DISEASE_STATE_DIRECTION": "REDUCED_FUNCTION"},
        {"RATIONALE_ID": "RAT-3", "TARGET_SYMBOL": "PLK1", "DECISION": DECISION_TEST_FIRST, "REASON_CODES": "CONTEXTUAL_HYPOTHESIS_WITH_PERTURBATION", "EVIDENCE_IDS": "ME-2", "DESIRED_TARGET_CHANGE": "RESTORE_NORMAL", "DISEASE_STATE_DIRECTION": "REDUCED_FUNCTION"},
    ])
    evidence = pd.DataFrame([{
        **{column: "" for column in compass.MECHANISM_EVIDENCE_COLUMNS},
        "EVIDENCE_ID": "ME-1", "TARGET_SYMBOL": "TTK", "CLAIM_TYPE": "PERTURBATION_EFFECT", "SUPPORT_DIRECTION": "SUPPORTS",
        "REVIEW_STATUS": "REVIEWED", "MODEL_SYSTEM": "human RPE1 cells", "MODEL_SYSTEM_LIMITATION": "",
    }], columns=compass.MECHANISM_EVIDENCE_COLUMNS)

    plans = build_test_plans(rationales, evidence).set_index("RATIONALE_ID")

    assert plans.loc["RAT-1", "PLAN_STATUS"] == "EVIDENCE_RESOLUTION_REQUIRED"
    assert "not applicable" in plans.loc["RAT-1", "PRIMARY_READOUT"]
    assert plans.loc["RAT-2", "PLAN_STATUS"] == "PROPOSED"
    assert "human RPE1 cells" == plans.loc["RAT-2", "MODEL_SYSTEM"]
    assert plans.loc["RAT-3", "PLAN_STATUS"] == "REQUIRES_MODEL_SELECTION"
    assert "experimental calibration" in plans.loc["RAT-2", "PERTURBATION_DIRECTION"]


# --------------------------------------------------------------------------- end to end


def test_run_compass_writes_every_contract_and_source_ledger(artifacts: dict[str, Path]) -> None:
    transport = FakeTransport()
    result = run_compass(
        artifacts["postanalysis"], artifacts["network"], artifacts["output"],
        cache_dir=artifacts["cache"], reviewed_evidence=artifacts["reviewed"], transport=transport,
    )

    assert result["status"] == "completed"
    assert result["advance_research"] == 0
    assert result["test_first"] >= 1
    assert result["reject"] >= 1
    for name in compass._COMPASS_TABLES:
        assert (artifacts["output"] / name).is_file(), name
    provenance = json.loads((artifacts["output"] / "compass_provenance.json").read_text())
    assert provenance["schema_version"] == 1
    assert provenance["identity"]["postanalysis_projection_hash"]
    for name, digest in provenance["output_sha256"].items():
        assert hashlib.sha256((artifacts["output"] / name).read_bytes()).hexdigest() == digest

    # Every cross-file reference resolves.
    rationales = pd.read_csv(artifacts["output"] / "target_mechanism_rationale.csv")
    drug_review = pd.read_csv(artifacts["output"] / "target_drug_review.csv")
    plans = pd.read_csv(artifacts["output"] / "target_test_plan.csv")
    ledger = pd.read_csv(artifacts["output"] / "target_rejection_ledger.csv")
    assert set(plans["RATIONALE_ID"]) == set(rationales["RATIONALE_ID"])
    assert set(drug_review["RATIONALE_ID"]).issubset(set(rationales["RATIONALE_ID"]))
    assert not ledger.empty

    # A compatible action is TEST_FIRST, the opposite action on the same rationale is REJECT.
    decisions = dict(zip(drug_review["DRUG_ID"], drug_review["DECISION"]))
    assert decisions["CHEMBL_OK"] == DECISION_TEST_FIRST
    assert decisions["CHEMBL_BAD"] == DECISION_REJECT

    source_status = pd.read_csv(artifacts["output"] / "source_status.csv")
    assert {"gene_identity", "membership", "protein_role", "literature", "pharmacology"}.issubset(set(source_status["STAGE"]))
    report = (artifacts["output"] / "compass_report.md").read_text()
    assert "PATHWAY_SUPPORTED" in report and "MECHANISM_PRIOR" in report
    assert "not a treatment proposal" in report


def test_ninety_percent_compass_report_and_drug_rows_are_distinct(artifacts: dict[str, Path]) -> None:
    evidence_path = artifacts["postanalysis"] / "integrated_pathway_evidence.csv"
    support_path = artifacts["postanalysis"] / "gene_pathway_support.csv"
    evidence = pd.read_csv(evidence_path)
    evidence = pd.concat([evidence, pd.DataFrame([
        _glue("EV-80", "KEGG", "hsa00080", "david_only", [], [{"Term": "exploratory", "Benjamini": 0.80}], "not_comparable", "significant"),
        _glue("EV-92", "KEGG", "hsa00092", "david_only", [], [{"Term": "outside", "Benjamini": 0.92}], "not_comparable", "not_significant"),
    ])], ignore_index=True)
    evidence.to_csv(evidence_path, index=False)
    support = pd.read_csv(support_path)
    support.loc[len(support)] = {"EVIDENCE_ID": "EV-80", "ONTOLOGY": "KEGG", "TERM_ID": "hsa00080", "PROVIDER": "DAVID", "GENE_SYMBOL": "BUB1B", "MEMBERSHIP_STATUS": "observed"}
    support.to_csv(support_path, index=False)

    result = run_compass(
        artifacts["postanalysis"], artifacts["network"], artifacts["output"],
        cache_dir=artifacts["cache"], reviewed_evidence=artifacts["reviewed"], transport=FakeTransport(),
        selection_profile="david_relaxed", david_benjamini_threshold=0.90, drug_scope="relaxed_pathway",
    )
    assert result["status"] == "completed"
    assert result["david_term_counts"] == {
        "input": 5, "validated_at_005": 3, "at_010": 3,
        "newly_in_cutoff": 1, "above_cutoff": 1, "unknown_statistic": 0,
        "eligible_exploratory": 1, "gated_at_cutoff": 1,
    }
    assert result["pharmacology_targets_queried"] == 1
    report = (artifacts["output"] / "compass_report.md").read_text()
    assert "RELAXED DAVID BENJAMINI 0.90 EXPLORATION" in report
    assert "0.05,0.90" in report
    assert "RELAXED DAVID BENJAMINI 0.10 EXPLORATION" not in report
    rows = pd.read_csv(artifacts["output"] / "target_drug_review.csv")
    assert not rows.empty
    assert set(rows["PATHWAY_EVIDENCE_LEVEL"]) == {"RELAXED_DAVID_BENJAMINI_0_90"}
    assert set(rows["TARGET_SYMBOL"]) == {"BUB1B"}


def test_offline_miss_is_reported_as_source_incomplete(artifacts: dict[str, Path]) -> None:
    result = run_compass(
        artifacts["postanalysis"], artifacts["network"], artifacts["output"],
        cache_dir=artifacts["cache"], reviewed_evidence=artifacts["reviewed"], offline=True,
    )

    assert result["status"] == "source_incomplete"
    assert result["source_failures"]
    source_status = pd.read_csv(artifacts["output"] / "source_status.csv")
    assert set(source_status["STATUS"]) == {sources.STATUS_OFFLINE_MISS}


def test_replaying_the_cache_recomputes_decisions_without_new_requests(artifacts: dict[str, Path]) -> None:
    run_compass(
        artifacts["postanalysis"], artifacts["network"], artifacts["output"],
        cache_dir=artifacts["cache"], reviewed_evidence=artifacts["reviewed"], transport=FakeTransport(),
    )
    before = pd.read_csv(artifacts["output"] / "target_mechanism_rationale.csv")

    changed = write_reviewed(artifacts["reviewed"], hazard_symbol="AURKA", conflict_symbol="BUB1B")
    resent = run_compass(
        artifacts["postanalysis"], artifacts["network"], artifacts["output"],
        cache_dir=artifacts["cache"], reviewed_evidence=changed, resume=True, transport=FakeTransport(failures={"anything"}),
    )

    assert resent["status"] == "completed"
    after = pd.read_csv(artifacts["output"] / "target_mechanism_rationale.csv")
    by_symbol = after.set_index("TARGET_SYMBOL")
    assert by_symbol.loc["BUB1B", "DECISION"] == DECISION_HOLD
    assert "EVIDENCE_CONFLICT" in by_symbol.loc["BUB1B", "REASON_CODES"]
    assert set(before["TARGET_SYMBOL"]) == set(after["TARGET_SYMBOL"])


def test_resume_requires_matching_fixed_identity(artifacts: dict[str, Path]) -> None:
    run_compass(
        artifacts["postanalysis"], artifacts["network"], artifacts["output"],
        cache_dir=artifacts["cache"], reviewed_evidence=artifacts["reviewed"], transport=FakeTransport(),
    )
    with pytest.raises(ValueError, match="already exists"):
        run_compass(
            artifacts["postanalysis"], artifacts["network"], artifacts["output"],
            cache_dir=artifacts["cache"], reviewed_evidence=artifacts["reviewed"], transport=FakeTransport(),
        )

    changed_postanalysis = write_postanalysis(artifacts["postanalysis"].parent / "postanalysis-2")
    pd.DataFrame([{
        "EVIDENCE_ID": "EV-9", "ONTOLOGY": "REACTOME", "TERM_ID": "R-HSA-9", "MATCH_STATUS": "gprofiler_only",
        "GPROFILER_ROWS_JSON": json.dumps([{"TERM_NAME": "Other", "BH_Q": 0.01}]), "DAVID_ROWS_JSON": "[]",
        "GPROFILER_SIGNIFICANCE": "significant", "DAVID_SIGNIFICANCE": "not_comparable",
        "EVIDENCE_TIER": "context_only", "GPROFILER_CONTRIBUTOR_STATUS": "complete",
        "DAVID_CONTRIBUTOR_STATUS": "absent", "CONTRIBUTOR_STATUS": "partial", "MVA_CONTEXT_JSON": "[]",
    }]).to_csv(changed_postanalysis / "integrated_pathway_evidence.csv", index=False)
    pd.DataFrame(columns=["EVIDENCE_ID", "ONTOLOGY", "TERM_ID", "PROVIDER", "GENE_SYMBOL", "MEMBERSHIP_STATUS"]).to_csv(changed_postanalysis / "gene_pathway_support.csv", index=False)

    with pytest.raises(ValueError, match="new output destination"):
        run_compass(
            changed_postanalysis, artifacts["network"], artifacts["output"],
            cache_dir=artifacts["cache"], reviewed_evidence=artifacts["reviewed"], resume=True, transport=FakeTransport(),
        )


def test_same_frozen_inputs_produce_identical_stable_tables(artifacts: dict[str, Path], tmp_path: Path) -> None:
    first = artifacts["output"]
    second = tmp_path / "compass-2"
    run_compass(artifacts["postanalysis"], artifacts["network"], first, cache_dir=artifacts["cache"], reviewed_evidence=artifacts["reviewed"], transport=FakeTransport())
    run_compass(artifacts["postanalysis"], artifacts["network"], second, cache_dir=artifacts["cache"], reviewed_evidence=artifacts["reviewed"], transport=FakeTransport(failures={"rest.uniprot.org"}))

    for name in compass._COMPASS_TABLES:
        if name == "source_status.csv":
            continue
        assert (first / name).read_bytes() == (second / name).read_bytes(), name


def test_cli_refuses_an_existing_destination_and_accepts_resume(artifacts: dict[str, Path], tmp_path: Path) -> None:
    script = Path(__file__).parents[1] / "scripts" / "run_compass.py"

    def command(output: Path, *extra: str) -> list[str]:
        return [
            sys.executable, str(script),
            "--postanalysis-dir", str(artifacts["postanalysis"]),
            "--network-dir", str(artifacts["network"]),
            "--output-dir", str(output),
            "--cache-dir", str(artifacts["cache"]),
            "--reviewed-evidence", str(artifacts["reviewed"]),
            "--offline",
            *extra,
        ]

    cold = subprocess.run(command(tmp_path / "cli-cold"), capture_output=True, text=True, check=False)
    assert cold.returncode == 2, cold.stderr  # a cold cache in offline mode is source-incomplete

    run_compass(
        artifacts["postanalysis"], artifacts["network"], artifacts["output"],
        cache_dir=artifacts["cache"], reviewed_evidence=artifacts["reviewed"], transport=FakeTransport(),
    )
    command = lambda *extra: [  # noqa: E731
        sys.executable, str(script),
        "--postanalysis-dir", str(artifacts["postanalysis"]),
        "--network-dir", str(artifacts["network"]),
        "--output-dir", str(artifacts["output"]),
        "--cache-dir", str(artifacts["cache"]),
        "--reviewed-evidence", str(artifacts["reviewed"]),
        *extra,
    ]
    refused = subprocess.run(command("--offline"), capture_output=True, text=True, check=False)
    assert refused.returncode == 1
    assert "already exists" in refused.stderr

    resumed = subprocess.run(command("--offline", "--resume"), capture_output=True, text=True, check=False)
    assert resumed.returncode == 0, resumed.stderr
    assert json.loads(resumed.stdout)["status"] == "completed"

    with pytest.raises(ValueError, match="missing required post-analysis artifact"):
        compass.load_postanalysis(artifacts["postanalysis"].parent)


def test_canonicalization_resolves_aliases_and_never_chooses_between_them(tmp_path: Path) -> None:
    def routes(method: str, url: str, params: object = None) -> Response:
        if url.endswith("/fetch/symbol/MPS1"):
            return Response({"response": {"numFound": 0, "docs": []}})
        if url.endswith("/fetch/alias_symbol/MPS1"):
            return Response({"response": {"numFound": 2, "docs": [{"symbol": "RPS27"}, {"symbol": "TTK"}]}})
        if url.endswith("/fetch/symbol/BUB1B"):
            return Response({"response": {"numFound": 1, "docs": [{"symbol": "BUB1B", "alias_symbol": ["BUBR1"], "prev_symbol": [], "ensembl_gene_id": "ENSG00000156970", "uniprot_ids": ["O60566"], "hgnc_id": "HGNC:1149"}]}})
        raise AssertionError(url)

    def transport(method: str, url: str, *, params: object = None) -> Response:
        return routes(method, url, params)

    canonical, statuses, ledger, failures = compass.canonicalize_targets(["MPS1", "BUB1B"], tmp_path / "cache", transport=transport)

    assert canonical == {"MPS1": "MPS1", "BUB1B": "BUB1B"}  # an unresolvable identity is never merged
    assert statuses["MPS1"] == sources.STATUS_AMBIGUOUS
    assert failures == []
    assert any(row["STAGE"] == "gene_identity" for row in ledger)
