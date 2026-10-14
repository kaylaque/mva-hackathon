"""Behavioural tests for the Compass public-source adapters.

All responses are synthetic and injected through a fake transport; no test touches
the network, and no patient or repository table is read.
"""

from __future__ import annotations

import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).parents[1] / "src"))

import pytest

from mva_hackathon import compass_sources as sources


class Response:
    """Minimal stand-in for a provider HTTP response."""

    def __init__(self, payload: object = None, *, status_code: int = 200, headers: dict[str, str] | None = None, text: str | None = None) -> None:
        self.status_code = status_code
        self.headers = headers or {}
        if text is None:
            self.text = "" if payload is None else json.dumps(payload)
        else:
            self.text = text
        self._payload = payload

    def json(self) -> object:
        if self._payload is None:
            raise ValueError("no json body")
        return self._payload


class FakeTransport:
    """Route-matched transport that records every call it receives."""

    def __init__(self, routes: list[tuple[str, object]]) -> None:
        self.routes = routes
        self.calls: list[tuple[str, str, object]] = []

    def __call__(self, method: str, url: str, *, params: object = None) -> Response:
        self.calls.append((method, url, params))
        for marker, responder in self.routes:
            if marker in url:
                return responder(method, url, params) if callable(responder) else responder
        raise AssertionError(f"unexpected request: {method} {url}")


def _cache(tmp_path: Path) -> Path:
    return tmp_path / "cache"


# --------------------------------------------------------------------------- requester


def test_verified_cache_response_is_reused_without_a_second_request(tmp_path: Path) -> None:
    transport = FakeTransport(
        [
            ("ontology/go/terms", Response(GO_TERM)),
            ("annotation/search", Response({"results": [], "pageInfo": {"total": 0, "resultsPerPage": 100}})),
        ]
    )
    first = sources.fetch_membership("GO", "GO:0007059", _cache(tmp_path), transport=transport)
    calls = len(transport.calls)
    second = sources.fetch_membership("GO", "GO:0007059", _cache(tmp_path), transport=transport)

    assert sources.REQUEST_TIMEOUT == 30
    assert first["status"] == sources.STATUS_EMPTY
    assert second["status"] == sources.STATUS_EMPTY
    assert calls == 2
    assert len(transport.calls) == 2
    assert sources.cache_file(_cache(tmp_path), "QUICKGO", second["request_id"]).is_file()


def test_offline_cache_miss_is_not_a_negative_result(tmp_path: Path) -> None:
    envelope = sources.fetch_membership("REACTOME", "R-HSA-69620", _cache(tmp_path), offline=True)

    assert envelope["status"] == sources.STATUS_OFFLINE_MISS
    assert envelope["records"] == []
    assert envelope["complete"] is False


def test_corrupt_cache_is_never_trusted(tmp_path: Path) -> None:
    routes = [
        ("ontology/go/terms", Response(GO_TERM)),
        ("annotation/search", Response({"results": [], "pageInfo": {"total": 0, "resultsPerPage": 100}})),
    ]
    transport = FakeTransport(routes)
    assert sources.fetch_membership("GO", "GO:0007059", _cache(tmp_path), transport=transport)["status"] == sources.STATUS_EMPTY
    identity = sources.request_identity("QUICKGO", f"{sources._QUICKGO}/annotation/search", "GET", {"goId": "GO:0007059", "taxonId": 9606, "limit": sources.PAGE_SIZE, "page": 1}, None)
    path = sources.cache_file(_cache(tmp_path), "QUICKGO", identity)
    assert path.is_file()
    path.write_text('{"adapter_version": "compass-sources-1", "response_sha256": "deadbeef", "payload": {"results": []}}', encoding="utf-8")

    # Online mode re-fetches only the corrupt request and stores a trustworthy record again.
    replacement = FakeTransport([("annotation/search", Response({"results": [], "pageInfo": {"total": 0, "resultsPerPage": 100}}))])
    assert sources.fetch_membership("GO", "GO:0007059", _cache(tmp_path), transport=replacement)["status"] == sources.STATUS_EMPTY
    assert len(replacement.calls) == 1

    path.write_text("not json at all", encoding="utf-8")
    offline = sources.fetch_membership("GO", "GO:0007059", _cache(tmp_path), offline=True)
    assert offline["status"] == sources.STATUS_OFFLINE_MISS


def test_transport_failures_retry_classify_and_stay_distinct(tmp_path: Path) -> None:
    def rate_limited(method: str, url: str, params: object = None) -> Response:
        return Response(status_code=429, headers={"Retry-After": "0.01"})

    exhausted = sources.fetch_membership("KEGG", "hsa04110", _cache(tmp_path / "limited"), transport=FakeTransport([("rest.kegg.jp", rate_limited)]))
    assert exhausted["status"] == sources.STATUS_UNAVAILABLE
    assert "429" in exhausted["reason"]

    def broken(method: str, url: str, params: object = None) -> Response:
        raise TimeoutError("synthetic timeout")

    timed_out = sources.fetch_membership("KEGG", "hsa04110", _cache(tmp_path / "timeout"), transport=FakeTransport([("rest.kegg.jp", broken)]))
    assert timed_out["status"] == sources.STATUS_UNAVAILABLE
    assert "TimeoutError" in timed_out["reason"]

    malformed = sources.fetch_membership(
        "KEGG",
        "hsa04110",
        _cache(tmp_path / "malformed"),
        transport=FakeTransport([("rest.kegg.jp", Response(status_code=200, text="<html>nope</html>"))]),
    )
    assert malformed["status"] == sources.STATUS_INVALID_RESPONSE

    server = sources.fetch_membership(
        "KEGG",
        "hsa04110",
        _cache(tmp_path / "server"),
        transport=FakeTransport([("rest.kegg.jp", Response(status_code=503, headers={"Retry-After": "0.01"}))]),
    )
    assert server["status"] == sources.STATUS_UNAVAILABLE


def test_unsupported_ontology_is_context_not_a_pathway(tmp_path: Path) -> None:
    envelope = sources.fetch_membership("HP", "HP:0000001", _cache(tmp_path))

    assert envelope["status"] == sources.STATUS_UNAVAILABLE
    assert "context" in envelope["reason"]


# --------------------------------------------------------------------------- GO


GO_TERM = {"numberOfHits": 1, "pageInfo": {"page": 1, "total": 1, "resultsPerPage": 1}, "results": [{"id": "GO:0007059", "aspect": "biological_process", "name": "chromosome segregation", "isObsolete": False}]}
GO_CC_TERM = {"numberOfHits": 1, "pageInfo": {"page": 1, "total": 1, "resultsPerPage": 1}, "results": [{"id": "GO:0000775", "aspect": "cellular_component", "name": "chromosome", "isObsolete": False}]}


def test_go_pagination_excludes_not_qualifiers_and_keeps_evidence_codes(tmp_path: Path) -> None:
    def annotations(method: str, url: str, params: object = None) -> Response:
        page = (params or {}).get("page")
        if page == 1:
            return Response(
                {
                    "results": [
                        {"id": "r1", "geneProductId": "UniProtKB:O60566", "symbol": "BUB1B", "goId": "GO:0007059", "goEvidence": "IDA", "qualifier": "involved_in", "taxonId": 9606},
                        {"id": "r2", "geneProductId": "UniProtKB:P12345", "symbol": "NOTGENE", "goId": "GO:0007059", "goEvidence": "IEA", "qualifier": "NOT|involved_in", "taxonId": 9606},
                    ],
                    "pageInfo": {"total": 3, "resultsPerPage": 2, "current": 1},
                }
            )
        return Response({"results": [{"id": "r3", "geneProductId": "UniProtKB:Q9H3R5", "symbol": "TTK", "goId": "GO:0007059", "goEvidence": "IMP", "qualifier": "involved_in", "taxonId": 9606}], "pageInfo": {"total": 3, "resultsPerPage": 2, "current": 2}})

    envelope = sources.fetch_membership(
        "GO",
        "GO:0007059",
        _cache(tmp_path),
        transport=FakeTransport([("ontology/go/terms", Response(GO_TERM)), ("annotation/search", annotations)]),
    )

    assert envelope["status"] == sources.STATUS_OK
    assert envelope["complete"] is True
    assert {record["SYMBOL"] for record in envelope["records"]} == {"BUB1B", "TTK"}
    assert [record["EVIDENCE_CODE"] for record in envelope["records"]] == ["IDA", "IMP"]
    assert [record["EVIDENCE_CLASS"] for record in envelope["records"]] == ["EXPERIMENTAL", "EXPERIMENTAL"]
    assert envelope["records"][0]["MEMBERSHIP_ORIGIN"] == "go_annotated_membership"


def test_incomplete_go_pagination_cannot_advance_as_complete(tmp_path: Path) -> None:
    def broken_second_page(method: str, url: str, params: object = None) -> Response:
        if (params or {}).get("page") == 1:
            return Response({"results": [{"id": "r1", "geneProductId": "UniProtKB:O60566", "symbol": "BUB1B", "goEvidence": "IDA"}], "pageInfo": {"total": 500, "resultsPerPage": 1}})
        return Response(status_code=500, headers={"Retry-After": "0.01"})

    envelope = sources.fetch_membership(
        "GO",
        "GO:0007059",
        _cache(tmp_path),
        transport=FakeTransport([("ontology/go/terms", Response(GO_TERM)), ("annotation/search", broken_second_page)]),
    )

    assert envelope["complete"] is False
    assert envelope["status"] == sources.STATUS_OK  # observed rows are kept, coverage is explicitly partial
    assert "paging" in envelope["reason"]


def test_go_cellular_component_is_annotation_context(tmp_path: Path) -> None:
    envelope = sources.fetch_membership("GO", "GO:0000775", _cache(tmp_path), transport=FakeTransport([("ontology/go/terms", Response(GO_CC_TERM))]))

    assert envelope["status"] == sources.STATUS_UNAVAILABLE
    assert "annotation context" in envelope["reason"]


# --------------------------------------------------------------------------- Reactome


EVENT = {"schemaClass": "Pathway", "stId": "R-HSA-69620", "displayName": "Cell Cycle Checkpoints", "speciesName": "Homo sapiens"}


def test_reactome_traverses_nested_complexes_and_excludes_molecules(tmp_path: Path) -> None:
    participants = [
        {"peDbId": 1, "refEntities": [{"schemaClass": "ReferenceGeneProduct", "identifier": "O60566", "displayName": "UniProt:O60566 BUB1B"}]},
        {"peDbId": 2, "refEntities": [{"schemaClass": "ReferenceMolecule", "identifier": "CHEBI:15377", "displayName": "water"}]},
        {"peDbId": 3},
    ]
    resolved = {"schemaClass": "Complex", "dbId": 3, "hasComponent": [{"dbId": 4}]}
    component = {"schemaClass": "EntityWithAccessionedSequence", "dbId": 4, "refEntities": [{"schemaClass": "ReferenceIsoform", "identifier": "Q9H3R5-1", "displayName": "UniProt:Q9H3R5 TTK"}]}

    def query_route(method: str, url: str, params: object = None) -> Response:
        if url.endswith("/data/query/R-HSA-69620"):
            return Response(EVENT)
        if url.endswith("/data/query/3"):
            return Response(resolved)
        if url.endswith("/data/query/4"):
            return Response({"schemaClass": "EntityWithAccessionedSequence", "dbId": 4, "refEntities": component["refEntities"]})
        raise AssertionError(url)

    envelope = sources.fetch_membership(
        "REACTOME",
        "R-HSA-69620",
        _cache(tmp_path),
        transport=FakeTransport([("/data/participants/", Response(participants)), ("/data/query/", query_route)]),
    )

    assert envelope["status"] == sources.STATUS_OK
    assert envelope["complete"] is True
    assert {record["SYMBOL"] for record in envelope["records"]} == {"BUB1B", "TTK"}
    assert {record["UNIPROT_ACCESSION"] for record in envelope["records"]} == {"O60566", "Q9H3R5"}
    assert all(record["MEMBERSHIP_ORIGIN"] == "reactome_participant" for record in envelope["records"])


def test_unresolved_reactome_entity_makes_membership_partial(tmp_path: Path) -> None:
    participants = [{"peDbId": 9}]

    envelope = sources.fetch_membership(
        "REACTOME",
        "R-HSA-69620",
        _cache(tmp_path),
        transport=FakeTransport(
            [
                ("/data/participants/", Response(participants)),
                ("/data/query/R-HSA-69620", Response(EVENT)),
                ("/data/query/9", Response(status_code=404)),
            ]
        ),
    )

    assert envelope["complete"] is False
    assert envelope["status"] == sources.STATUS_EMPTY
    assert "partial" in envelope["reason"]


def test_reactome_requires_a_human_stable_identifier(tmp_path: Path) -> None:
    envelope = sources.fetch_membership("REACTOME", "hsa04110", _cache(tmp_path))

    assert envelope["status"] == sources.STATUS_UNAVAILABLE
    assert "R-HSA-" in envelope["reason"]


# --------------------------------------------------------------------------- KEGG


KEGG_DEFINITION = """ENTRY       hsa04110                    Pathway
NAME        Cell cycle - Homo sapiens (human)
GENE        1017  CDK2; cyclin-dependent kinase 2 [KO:K02206]
            701  BUB1B; BUB1 mitotic checkpoint kinase B [KO:K06637]
            4085  MAD2L1; mitotic spindle assembly checkpoint protein MAD2A [KO:K02537]
COMPOUND    C00002  ATP
REFERENCE   PMID:7908906
///
"""


def test_kegg_maps_human_gene_identifiers_from_the_pathway_definition(tmp_path: Path) -> None:
    links = "path:hsa04110\thsa:1017\npath:hsa04110\thsa:701\npath:hsa04110\thsa:4085\n"

    def routes(method: str, url: str, params: object = None) -> Response:
        if "/get/" in url:
            return Response(text=KEGG_DEFINITION)
        return Response(text=links)

    envelope = sources.fetch_membership("KEGG", "hsa04110", _cache(tmp_path), transport=FakeTransport([("rest.kegg.jp", routes)]))

    assert envelope["status"] == sources.STATUS_OK
    assert envelope["complete"] is True
    assert {record["SYMBOL"] for record in envelope["records"]} == {"CDK2", "BUB1B", "MAD2L1"}
    assert all(record["SOURCE_RECORD_ID"].startswith("hsa:") for record in envelope["records"])


def test_kegg_reports_partial_when_a_linked_gene_has_no_symbol(tmp_path: Path) -> None:
    links = "path:hsa04110\thsa:1017\npath:hsa04110\thsa:999999\n"

    def routes(method: str, url: str, params: object = None) -> Response:
        if "/get/" in url:
            return Response(text=KEGG_DEFINITION)
        return Response(text=links)

    envelope = sources.fetch_membership("KEGG", "hsa04110", _cache(tmp_path), transport=FakeTransport([("rest.kegg.jp", routes)]))

    assert envelope["complete"] is False
    assert envelope["status"] == sources.STATUS_OK
    assert "partial" in envelope["reason"]


# --------------------------------------------------------------------------- WikiPathways


def _wp_document() -> dict[str, object]:
    return {
        "entitiesById": {
            "n1": {"id": "n1", "wpType": "GeneProduct", "gpmlElementName": "DataNode", "textContent": "BUB1B", "type": ["DataNode", "GeneProduct", "HGNC:BUB1B", "UniProt:O60566"], "xrefDataSource": "Ensembl", "xrefIdentifier": "ENSG00000156970"},
            "n2": {"id": "n2", "wpType": "Metabolite", "gpmlElementName": "DataNode", "textContent": "ATP", "type": ["DataNode", "Metabolite"]},
            "n3": {"id": "n3", "wpType": "GeneProduct", "gpmlElementName": "DataNode", "textContent": "unnamed", "type": ["DataNode", "GeneProduct", "Ensembl:ENSG00000000001"]},
        },
        "pathway": {"id": "WP5344", "organism": "Homo sapiens", "name": "Synthetic pathway", "contains": ["n1", "n2", "n3"]},
    }


def test_wikipathways_uses_typed_gene_identifiers_and_reports_partial(tmp_path: Path) -> None:
    envelope = sources.fetch_membership("WIKIPATHWAYS", "WP5344", _cache(tmp_path), transport=FakeTransport([("wikipathways", Response(_wp_document()))]))

    assert envelope["status"] == sources.STATUS_OK
    assert envelope["complete"] is False  # the Ensembl-only node has no resolvable symbol
    assert [record["SYMBOL"] for record in envelope["records"]] == ["BUB1B"]
    assert envelope["records"][0]["UNIPROT_ACCESSION"] == "O60566"


def test_wikipathways_rejects_a_mismatched_or_non_human_document(tmp_path: Path) -> None:
    wrong_id = _wp_document()
    wrong_id["pathway"] = {**wrong_id["pathway"], "id": "WP1111"}
    envelope = sources.fetch_membership("WIKIPATHWAYS", "WP5344", _cache(tmp_path / "id"), transport=FakeTransport([("wikipathways", Response(wrong_id))]))
    assert envelope["status"] == sources.STATUS_UNAVAILABLE
    assert "declared pathway identifier" in envelope["reason"]

    mouse = _wp_document()
    mouse["pathway"] = {**mouse["pathway"], "organism": "Mus musculus"}
    other = sources.fetch_membership("WIKIPATHWAYS", "WP5344", _cache(tmp_path / "mouse"), transport=FakeTransport([("wikipathways", Response(mouse))]))
    assert other["status"] == sources.STATUS_UNAVAILABLE

    changed = sources.fetch_membership("WIKIPATHWAYS", "WP5344", _cache(tmp_path / "shape"), transport=FakeTransport([("wikipathways", Response({"unexpected": True}))]))
    assert changed["status"] == sources.STATUS_INVALID_RESPONSE


# --------------------------------------------------------------------------- UniProt / HGNC


UNIPROT_BUB1B = {
    "results": [
        {
            "primaryAccession": "O60566",
            "genes": [{"geneName": {"value": "BUB1B"}}],
            "comments": [
                {"commentType": "FUNCTION", "texts": [{"value": "Essential component of the mitotic checkpoint."}]},
                {"commentType": "PTM", "texts": [{"value": "Phosphorylated upon DNA damage."}]},
            ],
            "entryAudit": {"lastAnnotationUpdateDate": "2026-01-01"},
        }
    ]
}


def test_uniprot_prefers_the_exact_primary_symbol_and_keeps_evidence_classes(tmp_path: Path) -> None:
    def annotations(method: str, url: str, params: object = None) -> Response:
        return Response({"results": [{"id": "a1", "goId": "GO:0007059", "goEvidence": "IDA", "qualifier": "involved_in"}], "pageInfo": {"total": 1, "resultsPerPage": 100}})

    envelope = sources.fetch_protein_evidence(
        "BUB1B",
        _cache(tmp_path),
        transport=FakeTransport([("rest.uniprot.org", Response(UNIPROT_BUB1B)), ("annotation/search", annotations)]),
    )

    record = envelope["records"][0]
    assert envelope["status"] == sources.STATUS_OK
    assert record["UNIPROT_ACCESSION"] == "O60566"
    assert record["PROTEIN_ROLE"].startswith("Essential component")
    assert record["MODIFICATION_ANNOTATION"].startswith("Phosphorylated")
    assert record["GO_ANNOTATIONS"][0]["EVIDENCE_CLASS"] == "EXPERIMENTAL"
    assert record["SPECIES"] == "Homo sapiens"


def test_uniprot_ambiguous_and_empty_results_are_separate_statuses(tmp_path: Path) -> None:
    ambiguous = {"results": [UNIPROT_BUB1B["results"][0], {**UNIPROT_BUB1B["results"][0], "primaryAccession": "A0A000"}]}
    envelope = sources.fetch_protein_evidence("BUB1B", _cache(tmp_path / "amb"), transport=FakeTransport([("rest.uniprot.org", Response(ambiguous))]))
    assert envelope["status"] == sources.STATUS_AMBIGUOUS

    empty = sources.fetch_protein_evidence("NOSUCHGENE", _cache(tmp_path / "empty"), transport=FakeTransport([("rest.uniprot.org", Response({"results": []}))]))
    assert empty["status"] == sources.STATUS_EMPTY
    assert empty["complete"] is False


def test_hgnc_resolves_aliases_and_reports_multiple_records_as_ambiguous(tmp_path: Path) -> None:
    docs = {"response": {"numFound": 1, "docs": [{"symbol": "TTK", "alias_symbol": ["MPS1", "PYT"], "prev_symbol": [], "ensembl_gene_id": "ENSG00000112742", "uniprot_ids": ["P33981"], "hgnc_id": "HGNC:12401"}]}}
    envelope = sources.fetch_gene_identity("MPS1", _cache(tmp_path), transport=FakeTransport([("rest.genenames.org", Response(docs))]))

    assert envelope["status"] == sources.STATUS_OK
    assert envelope["records"][0]["CANONICAL_SYMBOL"] == "TTK"
    assert envelope["records"][0]["REPORTED_SYMBOL"] == "MPS1"
    assert "MPS1" in envelope["records"][0]["ALIAS_SYMBOLS"]

    multiple = {"response": {"numFound": 2, "docs": docs["response"]["docs"] * 2}}
    other = sources.fetch_gene_identity("TTK", _cache(tmp_path / "multi"), transport=FakeTransport([("rest.genenames.org", Response(multiple))]))
    assert other["status"] == sources.STATUS_AMBIGUOUS


# --------------------------------------------------------------------------- literature


def test_literature_records_total_hits_and_marks_a_capped_search_incomplete(tmp_path: Path) -> None:
    def esearch(method: str, url: str, params: object = None) -> Response:
        term = (params or {}).get("term", "")
        if "mosaic variegated aneuploidy" in term:
            return Response({"esearchresult": {"count": "40", "idlist": ["1", "2", "3", "4", "5"]}})
        return Response({"esearchresult": {"count": "12", "idlist": ["1", "2"]}})

    def esummary(method: str, url: str, params: object = None) -> Response:
        docs = {str(item): {"title": f"T{item}", "source": "J1", "pubdate": "2024 Jan", "articleids": [{"idtype": "doi", "value": "10.1000/x"}]} for item in range(1, 6)}
        return Response({"result": docs})

    envelope = sources.fetch_literature("BUB1B", _cache(tmp_path), transport=FakeTransport([("esearch.fcgi", esearch), ("esummary.fcgi", esummary)]))

    assert envelope["status"] == sources.STATUS_OK
    assert envelope["complete"] is False
    assert {record["PMID"] for record in envelope["records"]} == {"1", "2", "3", "4", "5"}
    assert all(record["TOTAL_HITS"] == 40 for record in envelope["records"])
    assert "not a systematic review" in envelope["reason"]


# --------------------------------------------------------------------------- pharmacology


def test_drug_evidence_requires_authoritative_identity_and_verifies_target_accession(tmp_path: Path) -> None:
    no_ensembl = {"response": {"numFound": 1, "docs": [{"symbol": "BUB1B", "alias_symbol": [], "prev_symbol": [], "ensembl_gene_id": "", "uniprot_ids": ["O60566"], "hgnc_id": "HGNC:1149"}]}}
    unresolved = sources.fetch_drug_evidence("BUB1B", _cache(tmp_path / "identity"), transport=FakeTransport([("rest.genenames.org", Response(no_ensembl))]))
    assert unresolved["status"] == sources.STATUS_AMBIGUOUS

    docs = {"response": {"numFound": 1, "docs": [{"symbol": "BUB1B", "alias_symbol": [], "prev_symbol": [], "ensembl_gene_id": "ENSG00000156970", "uniprot_ids": ["O60566"], "hgnc_id": "HGNC:1149"}]}}
    mismatched = {"targets": [{"target_chembl_id": "CHEMBL1", "target_components": [{"accession": "Q00000"}]}], "page_meta": {"next": None}}

    def routes(method: str, url: str, params: object = None) -> Response:
        if "genenames" in url:
            return Response(docs)
        return Response(mismatched)

    evidence = sources.fetch_drug_evidence("BUB1B", _cache(tmp_path / "verified"), transport=FakeTransport([("rest.genenames.org", Response(docs)), ("chembl", routes)]))
    assert evidence["status"] == sources.STATUS_EMPTY
    assert "accession" in evidence["reason"]


def test_drug_evidence_pages_mechanisms_and_keeps_missing_assays_unavailable(tmp_path: Path) -> None:
    docs = {"response": {"numFound": 1, "docs": [{"symbol": "BUB1B", "alias_symbol": [], "prev_symbol": [], "ensembl_gene_id": "ENSG00000156970", "uniprot_ids": ["O60566"], "hgnc_id": "HGNC:1149"}]}}
    targets = {"targets": [{"target_chembl_id": "CHEMBL1234", "target_components": [{"accession": "O60566"}]}], "page_meta": {"next": None}}
    mechanisms = {"mechanisms": [{"mec_id": "M1", "molecule_chembl_id": "CHEMBL_A", "action_type": "INHIBITOR", "mechanism_of_action": "BUB1B inhibitor", "max_phase": 1}], "page_meta": {"next": None}}
    activities = {"activities": [], "page_meta": {"next": None}}
    ot = {"data": {"target": {"id": "ENSG00000156970", "approvedSymbol": "BUB1B", "drugAndClinicalCandidates": {"rows": []}}}}

    def routes(method: str, url: str, params: object = None) -> Response:
        if "genenames" in url:
            return Response(docs)
        if "target.json" in url:
            return Response(targets)
        if "mechanism.json" in url:
            return Response(mechanisms)
        if "activity.json" in url:
            return Response(activities)
        if "opentargets" in url:
            return Response(ot)
        raise AssertionError(url)

    envelope = sources.fetch_drug_evidence("BUB1B", _cache(tmp_path), transport=FakeTransport([("rest.genenames.org", Response(docs)), ("ebi.ac.uk", routes), ("opentargets", Response(ot))]))

    assert envelope["status"] == sources.STATUS_OK
    record = envelope["records"][0]
    assert record["ACTION_TYPE"] == "INHIBITOR"
    assert record["ASSAY_STATUS"] == "ASSAY_UNAVAILABLE"
    assert record["POTENCY"] == []
    assert record["ENSEMBL_GENE_ID"] == "ENSG00000156970"


def test_replayed_cache_issues_no_http_calls(tmp_path: Path) -> None:
    transport = FakeTransport([("rest.uniprot.org", Response(UNIPROT_BUB1B)), ("annotation/search", Response({"results": [], "pageInfo": {"total": 0, "resultsPerPage": 100}}))])
    first = sources.fetch_protein_evidence("BUB1B", _cache(tmp_path), transport=transport)
    calls = len(transport.calls)

    counting = FakeTransport([])
    second = sources.fetch_protein_evidence("BUB1B", _cache(tmp_path), transport=counting)

    assert calls > 0
    assert counting.calls == []
    assert second["status"] == first["status"] == sources.STATUS_OK
    assert second["response_sha256"] == first["response_sha256"]


def test_request_identity_covers_source_endpoint_method_params_and_species() -> None:
    base = sources.request_identity("GO", "https://example/go", "get", {"a": 1}, "human")
    assert base == sources.request_identity("GO", "https://example/go", "GET", {"a": 1}, "human")
    assert base != sources.request_identity("GO", "https://example/go", "GET", {"a": 2}, "human")
    assert base != sources.request_identity("KEGG", "https://example/go", "GET", {"a": 1}, "human")
    assert base != sources.request_identity("GO", "https://example/go", "GET", {"a": 1}, "mouse")


def test_gene_identity_resolves_a_unique_alias_and_records_alias_ambiguity(tmp_path: Path) -> None:
    def routes(method: str, url: str, params: object = None) -> Response:
        if url.endswith("/fetch/symbol/MPS1"):
            return Response({"response": {"numFound": 0, "docs": []}})
        if url.endswith("/fetch/alias_symbol/MPS1"):
            return Response({"response": {"numFound": 3, "docs": [{"symbol": "IDUA"}, {"symbol": "RPS27"}, {"symbol": "TTK"}]}})
        if url.endswith("/fetch/symbol/TTK"):
            return Response({"response": {"numFound": 0, "docs": []}})
        if url.endswith("/fetch/alias_symbol/TTK"):
            return Response({"response": {"numFound": 1, "docs": [{"symbol": "TTK", "alias_symbol": ["MPS1"], "prev_symbol": [], "ensembl_gene_id": "ENSG00000112742", "uniprot_ids": ["P33981"], "hgnc_id": "HGNC:12401"}]}})
        raise AssertionError(url)

    ambiguous = sources.fetch_gene_identity("MPS1", _cache(tmp_path / "amb"), transport=FakeTransport([("genenames", routes)]))
    assert ambiguous["status"] == sources.STATUS_AMBIGUOUS
    assert ambiguous["records"] == []
    assert "must be resolved manually" in ambiguous["reason"]

    resolved = sources.fetch_gene_identity("TTK", _cache(tmp_path / "ttk"), transport=FakeTransport([("genenames", routes)]))
    assert resolved["status"] == sources.STATUS_OK
    assert resolved["records"][0]["CANONICAL_SYMBOL"] == "TTK"
    assert resolved["records"][0]["RESOLUTION"] == "alias_symbol"
    assert "alias index" in resolved["reason"]
