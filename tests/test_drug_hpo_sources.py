from __future__ import annotations

from mva_hackathon.drug_hpo_sources import (
    fetch_fda_labels,
    fetch_trials,
    resolve_hpo_term,
    search_drug_hpo_literature,
)


class Response:
    def __init__(self, payload, status_code=200):
        self.payload = payload
        self.status_code = status_code
        self.text = payload if isinstance(payload, str) else ""

    def json(self):
        return self.payload


def test_resolve_hpo_term_requires_exact_nonobsolete_single_match(tmp_path):
    calls = []

    def transport(method, url, *, params):
        calls.append((method, url, params))
        return Response({"response": {"docs": [
            {"obo_id": "HP:0000001", "label": "All", "synonyms": ["Everything"], "iri": "http://purl.obolibrary.org/obo/HP_0000001", "is_obsolete": False, "ontology_version": "2025-01"},
            {"obo_id": "HP:0000002", "label": "near miss", "is_obsolete": False},
            {"obo_id": "HP:0000001", "label": "obsolete duplicate", "is_obsolete": True},
        ]}})

    result = resolve_hpo_term("HP:0000001", tmp_path, transport=transport)
    assert calls == [("GET", "https://www.ebi.ac.uk/ols4/api/search", {"q": "HP:0000001", "ontology": "hp", "exact": "true"})]
    assert result["status"] == "OK" and result["complete"] is True
    assert result["records"] == [{"HPO_CODE": "HP:0000001", "HPO_LABEL": "All", "HPO_SYNONYMS": ["Everything"], "IRI": "http://purl.obolibrary.org/obo/HP_0000001", "ONTOLOGY_VERSION": "2025-01", "ONTOLOGY_UPDATED": ""}]


def test_resolve_hpo_term_empty_ambiguous_and_malformed_are_distinct(tmp_path):
    docs = [[], [{"obo_id": "HP:0000001", "label": "one"}, {"obo_id": "HP:0000001", "label": "two"}], None, ["malformed"]]
    statuses = []
    for i, payload in enumerate(docs):
        result = resolve_hpo_term("HP:0000001", tmp_path / str(i), transport=lambda *_args, payload=payload, **_kw: Response({"response": {"docs": payload}}))
        statuses.append(result["status"])
        assert result["complete"] is (i < 2)
    assert statuses == ["EMPTY", "AMBIGUOUS", "INVALID_RESPONSE", "INVALID_RESPONSE"]


def test_literature_full_text_failure_keeps_records_incomplete(tmp_path):
    def transport(method, url, *, params):
        if url.endswith("esearch.fcgi"):
            return Response({"esearchresult": {"idlist": [], "count": "0"}})
        if url.endswith("/search"):
            return Response({"resultList": {"result": [{"pmid": "1", "pmcid": "PMC1", "title": "Article"}]}, "hitCount": 1, "nextCursorMark": "done"})
        if "/articles/" in url:
            return Response({}, 503)
        raise AssertionError(url)

    result = search_drug_hpo_literature("Drug", [], "HP:0000001", "Pain", [], tmp_path, transport=transport)
    assert result["status"] == "OK" and result["complete"] is False


def test_literature_query_pages_pubmed_and_europe_pmc_and_merges_records(tmp_path):
    seen = []

    def transport(method, url, *, params):
        seen.append((url, dict(params)))
        if url.endswith("esearch.fcgi"):
            start = int(params["retstart"])
            ids = ["1", "2"] if start == 0 else ["3"]
            return Response({"esearchresult": {"idlist": ids, "count": "3"}})
        if url.endswith("efetch.fcgi"):
            assert params["retmode"] == "xml"
            articles = "".join(f"<PubmedArticle><MedlineCitation><PMID>{x}</PMID><Article><ArticleTitle>PubMed {x}</ArticleTitle><Abstract><AbstractText>abstract</AbstractText></Abstract><Journal><Title>Journal</Title><JournalIssue><PubDate><Year>2020</Year></PubDate></JournalIssue></Journal><PublicationTypeList><PublicationType>Journal Article</PublicationType></PublicationTypeList></Article></MedlineCitation><PubmedData><ArticleIdList><ArticleId IdType='doi'>10/{x}</ArticleId></ArticleIdList></PubmedData></PubmedArticle>" for x in params["id"].split(","))
            return Response(f"<PubmedArticleSet>{articles}</PubmedArticleSet>")
        if url.endswith("/search"):
            start = params["cursorMark"]
            return Response({"resultList": {"result": [{"pmid": "1", "pmcid": "PMC1", "doi": "10/1", "title": "PubMed 1", "abstractText": "other abstract", "journalTitle": "J", "pubYear": "2021", "pubTypeList": {"pubType": ["Journal Article"]}}] if start == "*" else [], "hitCount": 1, "nextCursorMark": "done"}, "hitCount": 1})
        if "/articles/" in url:
            return Response("<article><body><p>Full text passage.</p></body></article>")
        raise AssertionError(url)

    result = search_drug_hpo_literature("Drug A", ["Alias A", "Drug A"], "HP:0000001", "Pain", ["Ache", "Pain"], tmp_path, transport=transport)
    assert result["source"] == "PUBMED+EUROPE_PMC" and result["status"] == "OK" and result["complete"] is True
    assert len(result["records"]) == 3
    assert result["records"][0]["PMID"] == "1"
    assert result["records"][0]["SOURCE_IDS"] == ["PUBMED", "EUROPE_PMC"]
    assert result["records"][0]["DOI"] == "10/1"
    assert "Drug A" in result["records"][0]["QUERY"] and '"Ache"' in result["records"][0]["QUERY"]
    assert result["records"][0]["PASSAGE_SOURCE"] == "full_text"
    assert any("fullTextXML" in call[0] for call in seen)
    pm = next(params for url, params in seen if url.endswith("esearch.fcgi"))
    assert pm["retmode"] == "json" and pm["retmax"] == 100 and pm["retstart"] == 0
    query = pm["term"]
    assert '"Drug A"' in query and '"Alias A"' in query and '"Pain"' in query and '"Ache"' in query
    assert any(url.endswith("esearch.fcgi") and params["retstart"] == 2 for url, params in seen)
    assert any(url.endswith("/search") and params["format"] == "json" and params["pageSize"] == 100 for url, params in seen)


def test_literature_page_cap_and_failed_page_are_incomplete(tmp_path, monkeypatch):
    import mva_hackathon.drug_hpo_sources as sources

    monkeypatch.setattr(sources, "MAX_PAGES", 1)
    def transport(method, url, **kw):
        if url.endswith("esearch.fcgi"):
            return Response({"esearchresult": {"idlist": ["1"], "count": "4"}})
        if url.endswith("efetch.fcgi"):
            return Response("<PubmedArticleSet><PubmedArticle><MedlineCitation><PMID>1</PMID><Article><ArticleTitle>One</ArticleTitle></Article></MedlineCitation></PubmedArticle></PubmedArticleSet>")
        if url.endswith("/search"):
            return Response({"resultList": {"result": [], "hitCount": 0}, "hitCount": 0})
        raise AssertionError(url)
    result = search_drug_hpo_literature("D", [], "HP:0000001", "Pain", [], tmp_path, transport=transport)
    assert result["complete"] is False
    assert result["status"] == "OK"


def test_fda_label_searches_generic_then_brand_and_paginates(tmp_path):
    seen = []

    def transport(method, url, *, params):
        seen.append((url, dict(params)))
        if "search" not in params:
            return Response({"results": []})
        if "generic_name" in params["search"]:
            return Response({"error": {"code": "NOT_FOUND"}}, 404)
        skip = params.get("skip", 0)
        remaining = max(0, 101 - skip)
        count = min(100, remaining)
        return Response({"meta": {"results": {"total": 101}}, "results": [{"id": str(skip + i), "set_id": "set", "effective_time": "20200101", "openfda": {"brand_name": ["Brand"]}, "indications_and_usage": "use", "boxed_warning": "box", "contraindications": "contra", "warnings": "warn", "adverse_reactions": "adverse", "pediatric_use": "peds"} for i in range(count)]})

    result = fetch_fda_labels("Generic", ["Brand"], tmp_path, transport=transport)
    assert result["source"] == "OPENFDA" and result["status"] == "OK" and result["complete"] is True
    assert len(result["records"]) == 101
    assert {r["SPL_ID"] for r in result["records"]} >= {"0", "100"}
    assert result["records"][0]["INDICATIONS"] == "use"
    assert result["records"][0]["CONTRAINDICATIONS"] == "contra"
    assert [params["skip"] for url, params in seen if params.get("skip") is not None and "brand_name" in params.get("search", "")] == [0, 100]
    assert "openfda.generic_name" in seen[0][1]["search"]
    assert any("openfda.brand_name" in params.get("search", "") for _, params in seen)


def test_fda_empty_is_empty_and_invalid_is_not_negative(tmp_path):
    empty = fetch_fda_labels("Missing", [], tmp_path / "empty", transport=lambda *_a, **_k: Response({"meta": {"results": {"total": 0}}, "results": []}))
    malformed = fetch_fda_labels("Drug", [], tmp_path / "bad", transport=lambda *_a, **_k: Response({"unexpected": []}))
    assert (empty["status"], empty["complete"]) == ("EMPTY", True)
    assert malformed["status"] == "INVALID_RESPONSE" and malformed["complete"] is False


def test_trials_query_and_pages_preserve_registration_context(tmp_path):
    seen = []

    def transport(method, url, *, params):
        seen.append((url, dict(params)))
        if params.get("pageToken"):
            return Response({"studies": []})
        return Response({"studies": [{"protocolSection": {"identificationModule": {"nctId": "NCT1"}, "designModule": {"studyType": "INTERVENTIONAL", "phases": ["PHASE2"]}, "statusModule": {"overallStatus": "COMPLETED", "resultsFirstPostDateStruct": {}}, "armsInterventionsModule": {"interventions": [{"name": "Drug A"}]}, "conditionsModule": {"conditions": ["Pain"]}, "descriptionModule": {"briefSummary": "Summary"}, "outcomesModule": {"primaryOutcomes": [{"measure": "Pain"}]}}}], "nextPageToken": "next"})

    result = fetch_trials("Drug A", ["Alias"], "Pain", ["Ache"], tmp_path, transport=transport)
    assert result["source"] == "CLINICALTRIALS" and result["status"] == "OK" and result["complete"] is True
    record = result["records"][0]
    assert record["NCT_ID"] == "NCT1" and record["STUDY_TYPE"] == "INTERVENTIONAL"
    assert record["OVERALL_STATUS"] == "COMPLETED" and record["PHASES"] == ["PHASE2"]
    assert record["INTERVENTIONS"] == ["Drug A"] and record["CONDITIONS"] == ["Pain"]
    assert record["BRIEF_SUMMARY"] == "Summary" and record["OUTCOME_MEASURES"] == ["Pain"]
    assert record["RESULTS_POSTED"] is False and record["WHY_STOPPED"] == ""
    query = seen[0][1]["filter.advanced"]
    assert "AREA[InterventionName]" in query and "Drug A" in query and "Alias" in query
    assert "AREA[ConditionSearch]" in query and "Pain" in query and "Ache" in query
    assert seen[0][1]["pageSize"] == 100 and seen[1][1]["pageToken"] == "next"


def test_trial_page_failure_is_not_complete(tmp_path):
    def transport(method, url, *, params):
        if not params.get("pageToken"):
            return Response({"studies": [], "nextPageToken": "next"})
        return Response({}, 503)
    result = fetch_trials("Drug", [], "Pain", [], tmp_path, transport=transport)
    assert result["complete"] is False and result["status"] == "UNAVAILABLE"

