"""Cached public evidence adapters for drug--HPO enrichment."""

from __future__ import annotations

import hashlib
import json
import re
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Any, Sequence

from . import compass_sources as cached

Transport = cached.Transport
MAX_PAGES = cached.MAX_PAGES
_PAGE_SIZE = 100
_OLS = "https://www.ebi.ac.uk/ols4/api/search"
_PUBMED = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
_EPMC = "https://www.ebi.ac.uk/europepmc/webservices/rest"
_OPENFDA = "https://api.fda.gov/drug/label.json"
_TRIALS = "https://clinicaltrials.gov/api/v2/studies"


def _envelope(status: str, records: list[dict[str, Any]], source: str, *, complete: bool = False, responses: Sequence[dict[str, Any]] = (), reason: str = "") -> dict[str, Any]:
    responses = list(responses)
    return {
        "status": status,
        "records": records,
        "complete": complete,
        "source": source,
        "request_id": ",".join(dict.fromkeys(r.get("request_id", "") for r in responses if r.get("request_id"))),
        "response_sha256": hashlib.sha256("".join(r.get("response_sha256", "") for r in responses).encode()).hexdigest() if responses else "",
        "retrieved_at": max((r.get("retrieved_at", "") for r in responses), default=""),
        "reason": reason or next((r.get("reason", "") for r in responses if r.get("reason")), ""),
    }


def _request(source: str, endpoint: str, params: Any, cache_dir: Path, offline: bool, transport: Transport | None, *, expect_json: bool = True) -> dict[str, Any]:
    return cached.perform_cached_request(source, endpoint, "GET", params, cache_dir=cache_dir, offline=offline, transport=transport, expect_json=expect_json)


def _failed(response: dict[str, Any], source: str, records: list[dict[str, Any]] | None = None, *, complete: bool = False, responses: Sequence[dict[str, Any]] = ()) -> dict[str, Any]:
    return _envelope(response["status"], records or [], source, complete=complete, responses=[*responses, response], reason=response.get("reason", ""))


def _text(value: Any) -> str:
    return str(value).strip() if value is not None else ""


def _list(value: Any) -> list[str]:
    if isinstance(value, list):
        return list(dict.fromkeys(_text(v) for v in value if _text(v)))
    return [_text(value)] if _text(value) else []


def _stable_terms(canonical: str, aliases: Sequence[str]) -> list[str]:
    return list(dict.fromkeys(v for v in (_text(canonical), *(_text(a) for a in aliases)) if v))


def _quote(value: str) -> str:
    return '"' + value.replace('"', '\\"') + '"'


def resolve_hpo_term(code: str, cache_dir: Path, *, offline: bool = False, transport: Transport | None = None) -> dict[str, Any]:
    source = "OLS4"
    response = _request(source, _OLS, {"q": code, "ontology": "hp", "exact": "true"}, cache_dir, offline, transport)
    if response["status"] != cached.STATUS_OK:
        return _failed(response, source)
    payload = response["payload"]
    docs = payload.get("response", {}).get("docs") if isinstance(payload, dict) and isinstance(payload.get("response"), dict) else None
    if not isinstance(docs, list) or any(not isinstance(doc, dict) or not isinstance(doc.get("obo_id"), str) for doc in docs):
        return _envelope(cached.STATUS_INVALID_RESPONSE, [], source, responses=[response], reason="unrecognized OLS4 search payload")
    exact = [doc for doc in docs if doc.get("obo_id") == code and str(doc.get("is_obsolete", "false")).lower() not in {"true", "1"}]
    if len(exact) > 1:
        return _envelope(cached.STATUS_AMBIGUOUS, [], source, complete=True, responses=[response], reason="multiple non-obsolete exact HPO documents")
    if not exact:
        return _envelope(cached.STATUS_EMPTY, [], source, complete=True, responses=[response], reason="no non-obsolete exact HPO document")
    doc = exact[0]
    if not isinstance(doc.get("label"), str) or not doc["label"].strip():
        return _envelope(cached.STATUS_INVALID_RESPONSE, [], source, responses=[response], reason="OLS4 exact document is missing a label")
    record = {"HPO_CODE": code, "HPO_LABEL": doc["label"].strip(), "HPO_SYNONYMS": _list(doc.get("synonyms")), "IRI": _text(doc.get("iri")), "ONTOLOGY_VERSION": _text(doc.get("ontology_version")), "ONTOLOGY_UPDATED": _text(doc.get("updated"))}
    return _envelope(cached.STATUS_OK, [record], source, complete=True, responses=[response])


def _canonical_query(drugs: Sequence[str], phenotypes: Sequence[str]) -> str:
    d = " OR ".join(_quote(v) for v in drugs)
    h = " OR ".join(_quote(v) for v in phenotypes)
    return f"({d}) AND ({h})"


def _xml_publications(body: str) -> list[dict[str, Any]] | None:
    try:
        root = ET.fromstring(body)
    except ET.ParseError:
        return None
    records = []
    for article in root.findall(".//PubmedArticle"):
        def parts(path: str) -> str:
            return " ".join("".join(node.itertext()).strip() for node in article.findall(path) if "".join(node.itertext()).strip())
        pmid = parts(".//PMID")
        if not pmid:
            continue
        pubtypes = [x.strip() for x in article.findall(".//PublicationType") for x in ["".join(x.itertext())] if x.strip()]
        records.append({"PMID": pmid, "PMCID": parts(".//ArticleId[@IdType='pmc']"), "DOI": parts(".//ArticleId[@IdType='doi']"), "TITLE": parts(".//ArticleTitle"), "ABSTRACT": parts(".//Abstract/AbstractText"), "PUBLICATION_TYPES": pubtypes, "JOURNAL": parts(".//Journal/Title"), "YEAR": parts(".//PubDate/Year") or parts(".//PubDate/MedlineDate")[:4]})
    return records


def _paragraphs(xml_text: str) -> list[str] | None:
    try:
        root = ET.fromstring(xml_text)
    except ET.ParseError:
        return None
    return [" ".join("".join(p.itertext()).split()) for p in root.findall(".//body//p") if " ".join("".join(p.itertext()).split())]


def _record_key(record: dict[str, Any]) -> tuple[str, str]:
    if record.get("PMID"):
        return ("PMID", record["PMID"])
    if record.get("DOI"):
        return ("DOI", record["DOI"].lower())
    return ("TITLE", " ".join(record.get("TITLE", "").casefold().split()))


def search_drug_hpo_literature(drug_name: str, aliases: Sequence[str], hpo_code: str, hpo_label: str, hpo_synonyms: Sequence[str], cache_dir: Path, *, offline: bool = False, transport: Transport | None = None) -> dict[str, Any]:
    drugs = _stable_terms(drug_name, aliases)
    phenotypes = _stable_terms(hpo_label, hpo_synonyms)
    query = _canonical_query(drugs, phenotypes)
    responses: list[dict[str, Any]] = []
    pubs: list[dict[str, Any]] = []
    complete = True
    # PubMed esearch + efetch, paged by the provider's total count.
    start = 0
    total: int | None = None
    for _ in range(MAX_PAGES):
        params = {"db": "pubmed", "term": query, "retmode": "json", "retmax": _PAGE_SIZE, "retstart": start}
        page = _request("PUBMED", f"{_PUBMED}/esearch.fcgi", params, cache_dir, offline, transport)
        responses.append(page)
        if page["status"] != cached.STATUS_OK:
            return _failed(page, "PUBMED+EUROPE_PMC", responses=responses[:-1])
        body = page["payload"].get("esearchresult") if isinstance(page["payload"], dict) else None
        ids = body.get("idlist") if isinstance(body, dict) else None
        try:
            total = int(body["count"])
        except (TypeError, ValueError, KeyError):
            return _envelope(cached.STATUS_INVALID_RESPONSE, [], "PUBMED+EUROPE_PMC", responses=responses, reason="unrecognized PubMed esearch payload")
        if not isinstance(ids, list):
            return _envelope(cached.STATUS_INVALID_RESPONSE, [], "PUBMED+EUROPE_PMC", responses=responses, reason="unrecognized PubMed identifier list")
        if ids:
            fetch = _request("PUBMED", f"{_PUBMED}/efetch.fcgi", {"db": "pubmed", "id": ",".join(map(str, ids)), "retmode": "xml"}, cache_dir, offline, transport, expect_json=False)
            responses.append(fetch)
            if fetch["status"] != cached.STATUS_OK:
                return _failed(fetch, "PUBMED+EUROPE_PMC", responses=responses[:-1])
            parsed = _xml_publications(fetch["payload"].get("text", ""))
            if parsed is None or len(parsed) != len(ids):
                return _envelope(cached.STATUS_INVALID_RESPONSE, [], "PUBMED+EUROPE_PMC", responses=responses, reason="unrecognized or incomplete PubMed efetch XML")
            for record in parsed:
                record["SOURCE_IDS"] = ["PUBMED"]
                record["QUERY"] = query
                record["PASSAGE_SOURCE"] = "abstract" if record["ABSTRACT"] else "title"
                record["RETRACTED"] = any("retract" in x.lower() for x in record["PUBLICATION_TYPES"])
                pubs.append(record)
        start += len(ids)
        if start >= total:
            break
    else:
        complete = False

    cursor = "*"
    seen_cursors: set[str] = set()
    epmc_total: int | None = None
    for _ in range(MAX_PAGES):
        if cursor in seen_cursors:
            complete = False
            break
        seen_cursors.add(cursor)
        page = _request("EUROPE_PMC", f"{_EPMC}/search", {"query": query, "format": "json", "pageSize": _PAGE_SIZE, "cursorMark": cursor}, cache_dir, offline, transport)
        responses.append(page)
        if page["status"] != cached.STATUS_OK:
            return _failed(page, "PUBMED+EUROPE_PMC", responses=responses[:-1])
        payload = page["payload"]
        resultlist = payload.get("resultList") if isinstance(payload, dict) else None
        results = resultlist.get("result") if isinstance(resultlist, dict) else None
        try:
            epmc_total = int(payload["hitCount"])
        except (TypeError, ValueError, KeyError):
            return _envelope(cached.STATUS_INVALID_RESPONSE, [], "PUBMED+EUROPE_PMC", responses=responses, reason="unrecognized Europe PMC search payload")
        if not isinstance(results, list):
            return _envelope(cached.STATUS_INVALID_RESPONSE, [], "PUBMED+EUROPE_PMC", responses=responses, reason="unrecognized Europe PMC result list")
        for entry in results:
            if not isinstance(entry, dict):
                continue
            record = {"PMID": _text(entry.get("pmid")), "PMCID": _text(entry.get("pmcid")), "DOI": _text(entry.get("doi")), "TITLE": _text(entry.get("title")), "ABSTRACT": _text(entry.get("abstractText")), "PUBLICATION_TYPES": _list(entry.get("pubTypeList", {}).get("pubType") if isinstance(entry.get("pubTypeList"), dict) else []), "JOURNAL": _text(entry.get("journalTitle")), "YEAR": _text(entry.get("pubYear")), "SOURCE_IDS": ["EUROPE_PMC"], "QUERY": query, "PASSAGE_SOURCE": "abstract" if entry.get("abstractText") else "title", "RETRACTED": any("retract" in x.lower() for x in _list(entry.get("pubTypeList", {}).get("pubType") if isinstance(entry.get("pubTypeList"), dict) else []))}
            pubs.append(record)
        next_cursor = _text(payload.get("nextCursorMark"))
        retrieved = sum(1 for p in pubs if "EUROPE_PMC" in p["SOURCE_IDS"])
        if retrieved >= epmc_total or not results:
            break
        if not next_cursor or next_cursor == cursor:
            complete = False
            break
        cursor = next_cursor
    else:
        complete = False

    merged: dict[tuple[str, str], dict[str, Any]] = {}
    for record in pubs:
        key = _record_key(record)
        if not key[1]:
            continue
        if key not in merged:
            merged[key] = record
            continue
        prior = merged[key]
        prior["SOURCE_IDS"] = list(dict.fromkeys([*prior["SOURCE_IDS"], *record["SOURCE_IDS"]]))
        for field in ("PMID", "PMCID", "DOI", "TITLE", "ABSTRACT", "JOURNAL", "YEAR"):
            if not prior.get(field) and record.get(field):
                prior[field] = record[field]
        prior["PUBLICATION_TYPES"] = list(dict.fromkeys([*prior["PUBLICATION_TYPES"], *record["PUBLICATION_TYPES"]]))
        prior["RETRACTED"] = bool(prior["RETRACTED"] or record["RETRACTED"])
    for record in merged.values():
        record["SOURCE_RECORD_ID"] = record.get("PMID") or record.get("DOI") or hashlib.sha256(record["TITLE"].casefold().encode()).hexdigest()
        if record.get("PMCID"):
            full = _request("EUROPE_PMC", f"{_EPMC}/articles/{record['PMCID']}/fullTextXML", {"format": "xml"}, cache_dir, offline, transport, expect_json=False)
            responses.append(full)
            if full["status"] != cached.STATUS_OK:
                complete = False
                continue
            paragraphs = _paragraphs(full["payload"].get("text", ""))
            if paragraphs is None:
                complete = False
                continue
            record["FULL_TEXT_PARAGRAPHS"] = paragraphs
            if paragraphs:
                record["PASSAGE_SOURCE"] = "full_text"
        record.setdefault("FULL_TEXT_PARAGRAPHS", [])
    records = sorted(merged.values(), key=lambda r: (r.get("PMID", ""), r.get("DOI", ""), r.get("TITLE", "").casefold()))
    return _envelope(cached.STATUS_OK if records else cached.STATUS_EMPTY, records, "PUBMED+EUROPE_PMC", complete=complete, responses=responses, reason="" if complete else "one or more provider pages or full-text details were incomplete")


def _fda_record(entry: dict[str, Any], query: str) -> dict[str, Any]:
    openfda = entry.get("openfda") if isinstance(entry.get("openfda"), dict) else {}
    return {"SPL_ID": _text(entry.get("id")), "FDA_SET_ID": _text(entry.get("set_id")), "EFFECTIVE_TIME": _text(entry.get("effective_time")), "GENERIC_NAMES": _list(openfda.get("generic_name")), "BRAND_NAMES": _list(openfda.get("brand_name")), "SUBSTANCE_NAMES": _list(openfda.get("substance_name")), "INDICATIONS": _text(" ".join(_list(entry.get("indications_and_usage")))), "BOXED_WARNINGS": _text(" ".join(_list(entry.get("boxed_warning")))), "CONTRAINDICATIONS": _text(" ".join(_list(entry.get("contraindications")))), "WARNINGS": _text(" ".join(_list(entry.get("warnings")))), "ADVERSE_REACTIONS": _text(" ".join(_list(entry.get("adverse_reactions")))), "PEDIATRIC_USE": _text(" ".join(_list(entry.get("pediatric_use")))), "SOURCE_QUERY": query}


def fetch_fda_labels(drug_name: str, aliases: Sequence[str], cache_dir: Path, *, offline: bool = False, transport: Transport | None = None) -> dict[str, Any]:
    responses: list[dict[str, Any]] = []

    def collect(field: str, term: str) -> tuple[list[dict[str, Any]], bool, dict[str, Any] | None]:
        query = f'openfda.{field}:"{term.replace(chr(34), chr(92) + chr(34))}"'
        skip = 0
        total: int | None = None
        records: list[dict[str, Any]] = []
        for _ in range(MAX_PAGES):
            response = _request("OPENFDA", _OPENFDA, {"search": query, "limit": _PAGE_SIZE, "skip": skip}, cache_dir, offline, transport)
            responses.append(response)
            if response["status"] != cached.STATUS_OK:
                if response["status"] == cached.STATUS_UNAVAILABLE and "http 404" in response.get("reason", ""):
                    responses[-1] = {**response, "status": cached.STATUS_EMPTY, "reason": ""}
                    return [], True, None
                return records, False, response
            payload = response["payload"]
            entries = payload.get("results") if isinstance(payload, dict) else None
            if not isinstance(entries, list):
                return records, False, {**response, "status": cached.STATUS_INVALID_RESPONSE, "reason": "unrecognized openFDA label payload"}
            try:
                total = int(payload["meta"]["results"]["total"])
            except (TypeError, ValueError, KeyError):
                return records, False, {**response, "status": cached.STATUS_INVALID_RESPONSE, "reason": "openFDA result total is missing"}
            if any(not isinstance(entry, dict) for entry in entries):
                return records, False, {**response, "status": cached.STATUS_INVALID_RESPONSE, "reason": "openFDA label entry is malformed"}
            records.extend(_fda_record(entry, query) for entry in entries)
            skip += len(entries)
            if skip >= total:
                return records, True, None
            if not entries:
                return records, False, {**response, "status": cached.STATUS_INVALID_RESPONSE, "reason": "openFDA page ended before reported total"}
        return records, False, None

    records, complete, failure = collect("generic_name", drug_name)
    if failure:
        return _failed(failure, "OPENFDA", records, responses=responses[:-1])
    if not records:
        for alias in _stable_terms("", aliases):
            alias_records, alias_complete, failure = collect("brand_name", alias)
            records.extend(alias_records)
            complete = complete and alias_complete
            if failure:
                return _failed(failure, "OPENFDA", records, responses=responses[:-1])
    unique = {r["SPL_ID"] or json.dumps(r, sort_keys=True): r for r in records}
    records = sorted(unique.values(), key=lambda r: (r["FDA_SET_ID"], r["SPL_ID"]))
    return _envelope(cached.STATUS_OK if records else cached.STATUS_EMPTY, records, "OPENFDA", complete=complete, responses=responses, reason="" if complete else "openFDA pagination was incomplete")


def _trial_query(drugs: Sequence[str], phenotypes: Sequence[str]) -> str:
    drug_query = " OR ".join(f'"{x.replace(chr(34), chr(92) + chr(34))}"' for x in drugs)
    hpo_query = " OR ".join(f'"{x.replace(chr(34), chr(92) + chr(34))}"' for x in phenotypes)
    return f"AREA[InterventionName]({drug_query}) AND AREA[ConditionSearch]({hpo_query})"


def _trial_record(study: dict[str, Any]) -> dict[str, Any]:
    protocol = study.get("protocolSection") if isinstance(study.get("protocolSection"), dict) else {}
    def module(key: str) -> dict[str, Any]:
        value = protocol.get(key)
        return value if isinstance(value, dict) else {}
    identification, design, status = module("identificationModule"), module("designModule"), module("statusModule")
    interventions = module("armsInterventionsModule").get("interventions", [])
    conditions = module("conditionsModule").get("conditions", [])
    descriptions = module("descriptionModule")
    outcomes = module("outcomesModule")
    measures = []
    for key in ("primaryOutcomes", "secondaryOutcomes", "otherOutcomes"):
        measures.extend(_text(x.get("measure")) for x in outcomes.get(key, []) if isinstance(x, dict) and _text(x.get("measure")))
    results_module = protocol.get("resultsSection")
    results_posted = bool(results_module or status.get("resultsFirstPostDateStruct") or status.get("resultsFirstPostDate"))
    return {"NCT_ID": _text(identification.get("nctId")), "STUDY_TYPE": _text(design.get("studyType")), "OVERALL_STATUS": _text(status.get("overallStatus")), "PHASES": _list(design.get("phases")), "INTERVENTIONS": [_text(x.get("name")) for x in interventions if isinstance(x, dict) and _text(x.get("name"))] if isinstance(interventions, list) else [], "CONDITIONS": _list(conditions), "BRIEF_SUMMARY": _text(descriptions.get("briefSummary")), "OUTCOME_MEASURES": measures, "RESULTS_POSTED": results_posted, "WHY_STOPPED": _text(status.get("whyStopped"))}


def fetch_trials(drug_name: str, aliases: Sequence[str], hpo_label: str, hpo_synonyms: Sequence[str], cache_dir: Path, *, offline: bool = False, transport: Transport | None = None) -> dict[str, Any]:
    query = _trial_query(_stable_terms(drug_name, aliases), _stable_terms(hpo_label, hpo_synonyms))
    responses: list[dict[str, Any]] = []
    records: list[dict[str, Any]] = []
    token = ""
    seen_tokens: set[str] = set()
    complete = True
    for _ in range(MAX_PAGES):
        if token in seen_tokens:
            complete = False
            break
        seen_tokens.add(token)
        params: dict[str, Any] = {"filter.advanced": query, "pageSize": _PAGE_SIZE}
        if token:
            params["pageToken"] = token
        response = _request("CLINICALTRIALS", _TRIALS, params, cache_dir, offline, transport)
        responses.append(response)
        if response["status"] != cached.STATUS_OK:
            return _failed(response, "CLINICALTRIALS", records, responses=responses[:-1])
        payload = response["payload"]
        studies = payload.get("studies") if isinstance(payload, dict) else None
        if not isinstance(studies, list):
            return _envelope(cached.STATUS_INVALID_RESPONSE, records, "CLINICALTRIALS", responses=responses, reason="unrecognized ClinicalTrials.gov studies payload")
        if any(not isinstance(study, dict) for study in studies):
            return _envelope(cached.STATUS_INVALID_RESPONSE, records, "CLINICALTRIALS", responses=responses, reason="ClinicalTrials.gov study entry is malformed")
        records.extend(_trial_record(s) for s in studies if isinstance(s, dict))
        next_token = _text(payload.get("nextPageToken"))
        if not next_token:
            break
        token = next_token
    else:
        complete = False
    unique = {r["NCT_ID"]: r for r in records if r["NCT_ID"]}
    records = sorted(unique.values(), key=lambda r: r["NCT_ID"])
    return _envelope(cached.STATUS_OK if records else cached.STATUS_EMPTY, records, "CLINICALTRIALS", complete=complete, responses=responses, reason="" if complete else "ClinicalTrials.gov pagination was incomplete")
