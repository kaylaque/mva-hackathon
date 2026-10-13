"""Public-source evidence retrieval for Compass, behind one cached requester.

Every request is addressed by a content hash over source, endpoint, method,
canonical public parameters, species and adapter schema version.  A verified
cache record may be reused; a corrupt record is never trusted.  Failures are
classified, never converted into biological negatives: an unreachable source is
``UNAVAILABLE`` and an offline cache miss is ``OFFLINE_MISS``, both distinct from
a successful ``EMPTY`` result.

Only public identifiers leave this process (gene symbols, pathway and ontology
identifiers).  No patient data, sample identifiers or repository tables are sent.
"""

from __future__ import annotations

from datetime import UTC, datetime
import hashlib
import json
import os
from pathlib import Path
import re
import time
from typing import Any, Callable

import requests

ADAPTER_VERSION = "compass-sources-1"
REQUEST_TIMEOUT = 30
MAX_ATTEMPTS = 3
BACKOFF_SECONDS = (1.0, 2.0, 4.0)
PAGE_SIZE = 100
UNIPROT_PAGE_SIZE = 500
MAX_PAGES = 200
USER_AGENT = "mva-hackathon/0.4.0 (SageBio MVA Hackathon 2026; research exploration; non-diagnostic)"

STATUS_OK = "OK"
STATUS_EMPTY = "EMPTY"
STATUS_UNAVAILABLE = "UNAVAILABLE"
STATUS_INVALID_RESPONSE = "INVALID_RESPONSE"
STATUS_AMBIGUOUS = "AMBIGUOUS"
STATUS_OFFLINE_MISS = "OFFLINE_MISS"

Transport = Callable[..., Any]

_REACTOME_CONTENT = "https://reactome.org/ContentService"
_QUICKGO = "https://www.ebi.ac.uk/QuickGO/services"
_KEGG = "https://rest.kegg.jp"
_WIKIPATHWAYS = "https://www.wikipathways.org/wikipathways-assets/pathways"
_UNIPROT = "https://rest.uniprot.org/uniprotkb/search"
_HGNC = "https://rest.genenames.org/fetch/symbol"
_HGNC_ALIAS = "https://rest.genenames.org/fetch/alias_symbol"
_NCBI = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
_CHEMBL = "https://www.ebi.ac.uk/chembl/api/data"
_OPEN_TARGETS = "https://api.platform.opentargets.org/api/v4/graphql"

_GENE_PRODUCT_CLASSES = {"ReferenceGeneProduct", "ReferenceIsoform"}
_CHEMBL_COLLECTIONS = {"target.json": "targets", "mechanism.json": "mechanisms", "activity.json": "activities"}
_NESTED_ENTITY_KEYS = ("hasComponent", "hasMember", "hasCandidate", "repeatedUnit", "entityOnOtherCell")
_EXPERIMENTAL_CODES = {"EXP", "IDA", "IPI", "IMP", "IGI", "IEP"}
_HIGH_THROUGHPUT_CODES = {"HTP", "HDA", "HMP", "HGI", "HEP"}


def _now() -> str:
    return datetime.now(UTC).isoformat()


def _digest(value: Any) -> str:
    return hashlib.sha256(json.dumps(value, sort_keys=True, default=str).encode()).hexdigest()


def _text(value: Any) -> str:
    if value is None:
        return ""
    return str(value).strip()


def _envelope(
    status: str,
    records: list[dict[str, Any]],
    *,
    source: str,
    complete: bool = False,
    request_id: str = "",
    response_sha256: str = "",
    retrieved_at: str = "",
    reason: str = "",
) -> dict[str, Any]:
    """One source envelope. ``complete`` is true only for a fully retrieved response."""
    return {
        "status": status,
        "records": records,
        "complete": bool(complete),
        "source": source,
        "request_id": request_id,
        "response_sha256": response_sha256,
        "retrieved_at": retrieved_at or _now(),
        "reason": reason,
    }


# --------------------------------------------------------------------------- cache


def request_identity(source: str, endpoint: str, method: str, params: Any, species: str | None = None) -> str:
    """Content address for one public request, including the adapter schema version."""
    return _digest(
        {
            "adapter": ADAPTER_VERSION,
            "source": source,
            "endpoint": endpoint,
            "method": method.upper(),
            "params": params or {},
            "species": species,
        }
    )


def cache_file(cache_dir: Path | str, source: str, request_id: str) -> Path:
    return Path(cache_dir) / "public-responses" / source / f"{request_id}.json"


def read_cached_response(path: Path) -> dict[str, Any] | None:
    """Return one verified cache record, or None when absent, corrupt or schema-stale."""
    if not path.is_file():
        return None
    try:
        record = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError, json.JSONDecodeError):
        return None
    if not isinstance(record, dict) or record.get("adapter_version") != ADAPTER_VERSION:
        return None
    payload = record.get("payload")
    if not isinstance(payload, (dict, list)) or record.get("response_sha256") != _digest(payload):
        return None
    return record


def write_cached_response(
    path: Path,
    *,
    request_id: str,
    source: str,
    endpoint: str,
    payload: Any,
    retrieved_at: str | None = None,
) -> dict[str, Any]:
    """Atomically store one complete public response.

    The record appears only after the temporary file is fully written and renamed,
    so an interrupted run leaves a complete record or no record at all.
    """
    record = {
        "adapter_version": ADAPTER_VERSION,
        "request_id": request_id,
        "source": source,
        "endpoint": endpoint,
        "retrieved_at": retrieved_at or _now(),
        "response_sha256": _digest(payload),
        "payload": payload,
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f"{path.name}.{os.getpid()}.partial")
    temporary.write_text(json.dumps(record, sort_keys=True), encoding="utf-8")
    os.replace(temporary, path)
    return record


# --------------------------------------------------------------------------- request


def default_transport(method: str, url: str, *, params: Any = None) -> Any:
    return requests.request(
        method,
        url,
        params=params,
        timeout=REQUEST_TIMEOUT,
        headers={"User-Agent": USER_AGENT, "Accept": "application/json"},
    )


def _decode(response: Any, expect_json: bool) -> Any | None:
    text = getattr(response, "text", None)
    if not expect_json:
        return {"text": text if isinstance(text, str) else ""}
    loader = getattr(response, "json", None)
    try:
        return loader() if callable(loader) else json.loads(text or "")
    except Exception:
        return None


def _retry_after(response: Any) -> float | None:
    headers = getattr(response, "headers", None)
    if not isinstance(headers, dict):
        return None
    raw = headers.get("Retry-After") or headers.get("retry-after")
    try:
        value = float(str(raw).strip())
    except (TypeError, ValueError):
        return None
    return value if 0 < value <= 60 else None


def _perform(
    source: str,
    endpoint: str,
    method: str,
    params: Any,
    *,
    cache_dir: Path | str,
    offline: bool,
    transport: Transport | None,
    species: str | None = None,
    expect_json: bool = True,
) -> dict[str, Any]:
    """Return a verified payload from cache or one politely retried public request."""
    request_id = request_identity(source, endpoint, method, params, species)
    path = cache_file(cache_dir, source, request_id)
    cached = read_cached_response(path)
    if cached is not None:
        return {
            "status": STATUS_OK,
            "payload": cached["payload"],
            "request_id": request_id,
            "response_sha256": cached["response_sha256"],
            "retrieved_at": cached["retrieved_at"],
            "reason": "reused verified cache snapshot",
        }
    if offline:
        return {
            "status": STATUS_OFFLINE_MISS,
            "payload": None,
            "request_id": request_id,
            "response_sha256": "",
            "retrieved_at": _now(),
            "reason": "offline mode: no verified cached response for this request",
        }
    caller = transport or default_transport
    last_reason = "no attempt was made"
    for attempt in range(MAX_ATTEMPTS):
        try:
            response = caller(method, endpoint, params=params)
        except (OSError, requests.RequestException) as error:  # transport failure, not a negative result
            last_reason = f"transport failure: {type(error).__name__}"
        else:
            code = int(getattr(response, "status_code", 0) or 0)
            if code == 200:
                payload = _decode(response, expect_json)
                if payload is None:
                    return {
                        "status": STATUS_INVALID_RESPONSE,
                        "payload": None,
                        "request_id": request_id,
                        "response_sha256": "",
                        "retrieved_at": _now(),
                        "reason": "response body is not valid JSON for this adapter",
                    }
                record = write_cached_response(path, request_id=request_id, source=source, endpoint=endpoint, payload=payload)
                return {
                    "status": STATUS_OK,
                    "payload": payload,
                    "request_id": request_id,
                    "response_sha256": record["response_sha256"],
                    "retrieved_at": record["retrieved_at"],
                    "reason": "",
                }
            if code == 429 or 500 <= code < 600:
                last_reason = f"http {code}"
                if attempt < MAX_ATTEMPTS - 1:
                    time.sleep(_retry_after(response) or BACKOFF_SECONDS[min(attempt, len(BACKOFF_SECONDS) - 1)])
                continue
            return {
                "status": STATUS_UNAVAILABLE,
                "payload": None,
                "request_id": request_id,
                "response_sha256": "",
                "retrieved_at": _now(),
                "reason": f"http {code}",
            }
        if attempt < MAX_ATTEMPTS - 1:
            time.sleep(BACKOFF_SECONDS[min(attempt, len(BACKOFF_SECONDS) - 1)])
    return {
        "status": STATUS_UNAVAILABLE,
        "payload": None,
        "request_id": request_id,
        "response_sha256": "",
        "retrieved_at": _now(),
        "reason": f"attempts exhausted: {last_reason}",
    }


def _forward(response: dict[str, Any], source: str, records: list[dict[str, Any]] | None = None, *, complete: bool = False) -> dict[str, Any]:
    return _envelope(
        response["status"],
        records or [],
        source=source,
        complete=complete,
        request_id=response["request_id"],
        response_sha256=response["response_sha256"],
        retrieved_at=response["retrieved_at"],
        reason=response["reason"],
    )


# --------------------------------------------------------------------------- memberships


def fetch_membership(
    ontology: str,
    term_id: str,
    cache_dir: Path,
    *,
    offline: bool = False,
    transport: Transport | None = None,
) -> dict[str, Any]:
    """Retrieve complete database membership for one pathway or process identifier."""
    source = _text(ontology).upper()
    term = _text(term_id)
    if not term:
        return _envelope(STATUS_UNAVAILABLE, [], source=source or "UNKNOWN", reason="no term identifier supplied")
    if source in {"REACTOME", "REAC"}:
        return _reactome_membership(term, cache_dir, offline=offline, transport=transport)
    if source in {"GO", "GO:BP"}:
        return _go_membership(term, cache_dir, offline=offline, transport=transport)
    if source == "KEGG":
        return _kegg_membership(term, cache_dir, offline=offline, transport=transport)
    if source in {"WIKIPATHWAYS", "WP"}:
        return _wikipathways_membership(term, cache_dir, offline=offline, transport=transport)
    return _envelope(
        STATUS_UNAVAILABLE,
        [],
        source=source or "UNKNOWN",
        reason="no membership adapter for this ontology; the annotation remains context, not a pathway",
    )


def _reactome_children(entity: dict[str, Any]) -> list[dict[str, Any]]:
    children: list[dict[str, Any]] = []
    for key in _NESTED_ENTITY_KEYS:
        value = entity.get(key)
        if isinstance(value, list):
            children.extend(item for item in value if isinstance(item, dict))
    return children


def _uniprot_accession(identifier: str) -> str:
    """Normalize a UniProt identifier, dropping an isoform suffix but inventing nothing."""
    candidate = _text(identifier)
    base = candidate.split("-")[0]
    return base if re.fullmatch(r"[A-Z0-9]{6,10}", base) else ""


def _reactome_symbol(reference: dict[str, Any]) -> str:
    """Reported gene symbol from a reference entity; never from the event display name."""
    gene = reference.get("geneName")
    if isinstance(gene, list) and gene:
        return _text(gene[0])
    if isinstance(gene, str) and gene.strip():
        return gene.strip()
    tokens = _text(reference.get("displayName")).split()
    tail = tokens[-1] if tokens else ""
    return tail if re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9\-\.]{0,19}", tail) else ""


def _reactome_walk(
    entities: list[Any],
    term: str,
    cache_dir: Path,
    offline: bool,
    transport: Transport | None,
) -> tuple[list[dict[str, Any]], bool]:
    records: list[dict[str, Any]] = []
    seen: set[str] = set()
    visited: set[Any] = set()
    complete = True
    queue = list(entities)
    while queue:
        entity = queue.pop(0)
        if not isinstance(entity, dict):
            complete = False
            continue
        references = entity.get("refEntities")
        if not references:
            db_id = entity.get("dbId") or entity.get("peDbId") or entity.get("id")
            if db_id is None:
                complete = False
                continue
            if db_id in visited:
                # Already expanded: a repeated reference is not missing coverage.
                continue
            visited.add(db_id)
            resolved = _perform("REACTOME", f"{_REACTOME_CONTENT}/data/query/{db_id}", "GET", None, cache_dir=cache_dir, offline=offline, transport=transport)
            if resolved["status"] != STATUS_OK or not isinstance(resolved["payload"], dict):
                complete = False
                continue
            payload = resolved["payload"]
            queue.append(payload)
            if not payload.get("refEntities"):
                queue.extend(_reactome_children(payload))
            continue
        if not isinstance(references, list):
            complete = False
            continue
        for reference in references:
            if not isinstance(reference, dict) or reference.get("schemaClass") not in _GENE_PRODUCT_CLASSES:
                continue
            accession = _text(reference.get("identifier"))
            symbol = _reactome_symbol(reference)
            key = accession or symbol
            if not key:
                complete = False
                continue
            if key in seen:
                continue
            seen.add(key)
            records.append(
                {
                    "SYMBOL": symbol.upper(),
                    "UNIPROT_ACCESSION": _uniprot_accession(accession),
                    "SOURCE_RECORD_ID": key,
                    "MEMBERSHIP_ORIGIN": "reactome_participant",
                }
            )
        queue.extend(_reactome_children(entity))
    return records, complete


def _reactome_membership(term: str, cache_dir: Path, *, offline: bool, transport: Transport | None) -> dict[str, Any]:
    stable = term.upper()
    if re.fullmatch(r"R-HSA-\d+", stable) is None:
        return _envelope(STATUS_UNAVAILABLE, [], source="REACTOME", reason="a human Reactome stable identifier (R-HSA-*) is required")
    event = _perform("REACTOME", f"{_REACTOME_CONTENT}/data/query/{stable}", "GET", None, cache_dir=cache_dir, offline=offline, transport=transport)
    if event["status"] != STATUS_OK:
        return _forward(event, "REACTOME")
    payload = event["payload"]
    if not isinstance(payload, dict) or not isinstance(payload.get("schemaClass"), str):
        return _envelope(STATUS_INVALID_RESPONSE, [], source="REACTOME", reason="unrecognized Reactome event payload")
    if _text(payload.get("speciesName")).lower() != "homo sapiens":
        return _envelope(STATUS_UNAVAILABLE, [], source="REACTOME", reason="only human events are supported")
    participants = _perform("REACTOME", f"{_REACTOME_CONTENT}/data/participants/{stable}", "GET", None, cache_dir=cache_dir, offline=offline, transport=transport)
    if participants["status"] != STATUS_OK:
        return _forward(participants, "REACTOME")
    body = participants["payload"]
    if not isinstance(body, list):
        return _envelope(STATUS_INVALID_RESPONSE, [], source="REACTOME", reason="unrecognized Reactome participant payload")
    records, complete = _reactome_walk(body, stable, cache_dir, offline, transport)
    return _envelope(
        STATUS_OK if records else STATUS_EMPTY,
        records,
        source="REACTOME",
        complete=complete,
        request_id=participants["request_id"],
        response_sha256=participants["response_sha256"],
        retrieved_at=participants["retrieved_at"],
        reason="" if complete else "some referenced entities or gene products could not be resolved; membership is partial",
    )


def _go_pages(
    endpoint: str,
    params: dict[str, Any],
    cache_dir: Path,
    offline: bool,
    transport: Transport | None,
) -> tuple[list[dict[str, Any]], bool, dict[str, Any]]:
    """Follow source end markers with cycle detection; never silently truncate."""
    rows: list[dict[str, Any]] = []
    complete = True
    seen_pages: set[int] = set()
    last: dict[str, Any] = {"status": STATUS_INVALID_RESPONSE, "reason": "", "request_id": "", "response_sha256": "", "retrieved_at": _now()}
    page = 1
    while page <= MAX_PAGES:
        if page in seen_pages:
            complete = False
            break
        seen_pages.add(page)
        response = _perform("QUICKGO", endpoint, "GET", {**params, "page": page}, cache_dir=cache_dir, offline=offline, transport=transport)
        last = response
        if response["status"] != STATUS_OK:
            return rows, False, response
        payload = response["payload"]
        results = payload.get("results") if isinstance(payload, dict) else None
        info = payload.get("pageInfo") if isinstance(payload, dict) else None
        if not isinstance(results, list) or not isinstance(info, dict):
            return rows, False, {**response, "status": STATUS_INVALID_RESPONSE, "reason": "unrecognized QuickGO paged payload"}
        rows.extend(row for row in results if isinstance(row, dict))
        total = info.get("total")
        per_page = info.get("resultsPerPage") or PAGE_SIZE
        try:
            per_page = int(per_page)
        except (TypeError, ValueError):
            per_page = PAGE_SIZE
        projected = page * per_page
        if not results or (isinstance(total, int) and projected >= total) or (not isinstance(total, int) and len(results) < per_page):
            break
        page += 1
    else:
        complete = False
    return rows, complete, last


def _quickgo_annotations(go_id: str, cache_dir: Path, offline: bool, transport: Transport | None) -> tuple[list[dict[str, Any]], bool, dict[str, Any]]:
    rows, complete, last = _go_pages(
        f"{_QUICKGO}/annotation/search",
        {"goId": go_id, "taxonId": 9606, "limit": PAGE_SIZE},
        cache_dir,
        offline,
        transport,
    )
    records: list[dict[str, Any]] = []
    for row in rows:
        qualifier = _text(row.get("qualifier")).upper()
        if qualifier.startswith("NOT") or "|NOT|" in f"|{qualifier}|":
            continue
        product = _text(row.get("geneProductId"))
        records.append(
            {
                "SYMBOL": _text(row.get("symbol")).upper(),
                "UNIPROT_ACCESSION": product.split(":")[-1] if product.startswith("UniProtKB:") else "",
                "EVIDENCE_CODE": _text(row.get("goEvidence")),
                "EVIDENCE_CLASS": _evidence_class(_text(row.get("goEvidence"))),
                "GO_ID": _text(row.get("goId")),
                "TAXON": row.get("taxonId"),
                "SOURCE_RECORD_ID": _text(row.get("id")),
            }
        )
    return records, complete, last


def _evidence_class(code: str) -> str:
    """Map a GO evidence code to its documented class; never promote it to a claim."""
    upper = code.upper()
    if upper in _EXPERIMENTAL_CODES:
        return "EXPERIMENTAL"
    if upper in _HIGH_THROUGHPUT_CODES:
        return "HIGH_THROUGHPUT"
    if upper in {"TAS", "NAS", "IC", "ISS", "ISO", "ISA", "ISM", "IGC", "IBA", "IBD", "IKR", "IRD", "RCA"}:
        return "CURATED_OR_INFERRED"
    if upper == "IEA":
        return "COMPUTATIONAL"
    return "UNCLASSIFIED"


def _go_membership(term: str, cache_dir: Path, *, offline: bool, transport: Transport | None) -> dict[str, Any]:
    go_id = term.upper()
    if re.fullmatch(r"GO:\d{7}", go_id) is None:
        return _envelope(STATUS_UNAVAILABLE, [], source="GO", reason="a GO identifier is required for GO membership")
    identity = _perform("QUICKGO", f"{_QUICKGO}/ontology/go/terms/{go_id}", "GET", None, cache_dir=cache_dir, offline=offline, transport=transport)
    if identity["status"] != STATUS_OK:
        return _forward(identity, "GO")
    payload = identity["payload"]
    results = payload.get("results") if isinstance(payload, dict) else None
    if not isinstance(results, list) or not results or not isinstance(results[0], dict):
        return _envelope(STATUS_INVALID_RESPONSE, [], source="GO", reason="unrecognized QuickGO term payload")
    aspect = _text(results[0].get("aspect"))
    if aspect != "biological_process":
        return _envelope(
            STATUS_UNAVAILABLE,
            [],
            source="GO",
            request_id=identity["request_id"],
            response_sha256=identity["response_sha256"],
            retrieved_at=identity["retrieved_at"],
            reason=f"GO aspect '{aspect or 'unknown'}' is annotation context, not a biological process",
        )
    records, complete, last = _quickgo_annotations(go_id, cache_dir, offline, transport)
    if not records and last["status"] != STATUS_OK:
        # A page-level source failure must stay a source status, never an empty success.
        return _forward(last, "GO")
    for record in records:
        record["MEMBERSHIP_ORIGIN"] = "go_annotated_membership"
    return _envelope(
        STATUS_OK if records else STATUS_EMPTY,
        records,
        source="GO",
        complete=complete,
        request_id=last["request_id"],
        response_sha256=last["response_sha256"],
        retrieved_at=last["retrieved_at"],
        reason="" if complete else "GO annotation paging did not reach a complete end marker",
    )


def _kegg_symbols(definition: str) -> dict[str, str]:
    """Map KEGG gene numbers to symbols from the pathway GENE section only."""
    symbols: dict[str, str] = {}
    inside = False
    for line in definition.splitlines():
        if line[:1].strip() and not line[0].isspace():
            # A new top-level section; only GENE carries gene number/symbol pairs.
            inside = line.startswith("GENE")
            if inside:
                line = line[len("GENE"):]
            else:
                continue
        if not inside:
            continue
        match = re.match(r"\s*(\d+)\s+([A-Za-z0-9_\-\.]+)\s*[;,\[]", line)
        if match:
            symbols[match.group(1)] = match.group(2).upper()
    return symbols


def _kegg_membership(term: str, cache_dir: Path, *, offline: bool, transport: Transport | None) -> dict[str, Any]:
    identifier = term.lower()
    if re.fullmatch(r"hsa\d{5}", identifier) is None:
        return _envelope(STATUS_UNAVAILABLE, [], source="KEGG", reason="a human KEGG pathway identifier (hsa#####) is required")
    definition = _perform("KEGG", f"{_KEGG}/get/{identifier}", "GET", None, cache_dir=cache_dir, offline=offline, transport=transport, expect_json=False)
    if definition["status"] != STATUS_OK:
        return _forward(definition, "KEGG")
    links = _perform("KEGG", f"{_KEGG}/link/hsa/{identifier}", "GET", None, cache_dir=cache_dir, offline=offline, transport=transport, expect_json=False)
    if links["status"] != STATUS_OK:
        return _forward(links, "KEGG")
    definition_text = _text(definition["payload"].get("text"))
    link_text = _text(links["payload"].get("text"))
    if "ENTRY" not in definition_text or "path:" not in link_text:
        return _envelope(STATUS_INVALID_RESPONSE, [], source="KEGG", reason="unrecognized KEGG pathway payload")
    if "Homo sapiens" not in definition_text:
        return _envelope(STATUS_UNAVAILABLE, [], source="KEGG", reason="only human pathways are supported")
    symbols = _kegg_symbols(definition_text)
    gene_ids = sorted({match for match in re.findall(r"hsa:(\d+)", link_text)})
    records: list[dict[str, Any]] = []
    unmapped = 0
    for gene in gene_ids:
        symbol = symbols.get(gene)
        if symbol is None:
            unmapped += 1
            continue
        records.append(
            {
                "SYMBOL": symbol,
                "UNIPROT_ACCESSION": "",
                "SOURCE_RECORD_ID": f"hsa:{gene}",
                "MEMBERSHIP_ORIGIN": "kegg_pathway_link",
            }
        )
    complete = bool(records) and unmapped == 0
    reason = "" if complete else f"{unmapped} linked gene(s) had no symbol in the pathway definition, so membership is partial"
    return _envelope(
        STATUS_OK if records else STATUS_EMPTY,
        records,
        source="KEGG",
        complete=complete,
        request_id=links["request_id"],
        response_sha256=links["response_sha256"],
        retrieved_at=links["retrieved_at"],
        reason=reason,
    )


def _wikipathways_gene(node: dict[str, Any]) -> tuple[str, str]:
    """Gene symbol and accession from a GeneProduct node's database identifiers, not labels."""
    tokens = node.get("type") if isinstance(node.get("type"), list) else []
    symbol = ""
    accession = ""
    for token in tokens:
        text = _text(token)
        if text.startswith("HGNC:") and not symbol:
            symbol = text[5:].strip().upper()
        elif text.startswith("UniProt:") and not accession:
            accession = text[len("UniProt:"):].strip()
    if not symbol:
        source = _text(node.get("xrefDataSource"))
        identifier = _text(node.get("xrefIdentifier"))
        if source in {"HGNC", "HGNC Accession Number"} and identifier and not identifier.isdigit():
            symbol = identifier.upper()
    if not accession and _text(node.get("xrefDataSource")) in {"UniProt", "Uniprot-TrEMBL", "Uniprot-SwissProt"}:
        accession = _text(node.get("xrefIdentifier"))
    return symbol, accession


def _wikipathways_membership(term: str, cache_dir: Path, *, offline: bool, transport: Transport | None) -> dict[str, Any]:
    wp_id = term.upper()
    if re.fullmatch(r"WP\d+", wp_id) is None:
        return _envelope(STATUS_UNAVAILABLE, [], source="WIKIPATHWAYS", reason="a WikiPathways identifier (WP####) is required")
    response = _perform("WIKIPATHWAYS", f"{_WIKIPATHWAYS}/{wp_id}/{wp_id}.json", "GET", None, cache_dir=cache_dir, offline=offline, transport=transport)
    if response["status"] != STATUS_OK:
        return _forward(response, "WIKIPATHWAYS")
    payload = response["payload"]
    pathway = payload.get("pathway") if isinstance(payload, dict) else None
    entities = payload.get("entitiesById") if isinstance(payload, dict) else None
    if not isinstance(pathway, dict) or not isinstance(entities, dict):
        return _envelope(
            STATUS_INVALID_RESPONSE,
            [],
            source="WIKIPATHWAYS",
            request_id=response["request_id"],
            response_sha256=response["response_sha256"],
            retrieved_at=response["retrieved_at"],
            reason="the served WikiPathways document does not match the expected pathway/entities schema",
        )
    if _text(pathway.get("id")).upper() != wp_id:
        return _envelope(STATUS_UNAVAILABLE, [], source="WIKIPATHWAYS", reason="declared pathway identifier does not match the requested pathway")
    if "homo sapiens" not in _text(pathway.get("organism")).lower():
        return _envelope(STATUS_UNAVAILABLE, [], source="WIKIPATHWAYS", reason="only human pathways are supported")
    records: list[dict[str, Any]] = []
    unmapped = 0
    seen: set[str] = set()
    for node_id in pathway.get("contains") or []:
        node = entities.get(node_id)
        if not isinstance(node, dict):
            continue
        node_types = node.get("type") if isinstance(node.get("type"), list) else []
        if node.get("wpType") != "GeneProduct" and "GeneProduct" not in node_types:
            continue
        symbol, accession = _wikipathways_gene(node)
        if not symbol:
            unmapped += 1
            continue
        if symbol in seen:
            continue
        seen.add(symbol)
        records.append(
            {
                "SYMBOL": symbol,
                "UNIPROT_ACCESSION": accession,
                "SOURCE_RECORD_ID": _text(node.get("id")) or node_id,
                "MEMBERSHIP_ORIGIN": "wikipathways_gene_product",
            }
        )
    complete = bool(records) and unmapped == 0
    return _envelope(
        STATUS_OK if records else STATUS_EMPTY,
        records,
        source="WIKIPATHWAYS",
        complete=complete,
        request_id=response["request_id"],
        response_sha256=response["response_sha256"],
        retrieved_at=response["retrieved_at"],
        reason="" if complete else f"{unmapped} gene product node(s) had no resolvable gene identifier, so membership is partial",
    )


# --------------------------------------------------------------------------- protein evidence


def _primary_symbol(entry: dict[str, Any]) -> str:
    genes = entry.get("genes")
    if isinstance(genes, list) and genes and isinstance(genes[0], dict):
        name = genes[0].get("geneName")
        if isinstance(name, dict):
            return _text(name.get("value")).upper()
    return ""


def _comments(entry: dict[str, Any], comment_type: str) -> list[str]:
    texts: list[str] = []
    for comment in entry.get("comments") or []:
        if isinstance(comment, dict) and comment.get("commentType") == comment_type:
            for text in comment.get("texts") or []:
                if isinstance(text, dict):
                    value = _text(text.get("value"))
                    if value:
                        texts.append(value)
    return texts


def _hgnc_docs(payload: Any) -> list[dict[str, Any]]:
    body = payload.get("response") if isinstance(payload, dict) else None
    docs = body.get("docs") if isinstance(body, dict) else None
    return [doc for doc in docs if isinstance(doc, dict)] if isinstance(docs, list) else []


def fetch_gene_identity(
    symbol: str,
    cache_dir: Path,
    *,
    offline: bool = False,
    transport: Transport | None = None,
) -> dict[str, Any]:
    """Resolve a canonical HGNC symbol, recording ambiguity instead of choosing.

    An exact symbol match wins.  When no exact symbol exists the alias index is consulted;
    a unique alias resolves, while several equally valid aliases leave identity unresolved.
    """
    clean = _text(symbol).upper()
    if not clean:
        return _envelope(STATUS_UNAVAILABLE, [], source="HGNC", reason="a gene symbol is required")
    exact = _perform("HGNC", f"{_HGNC}/{clean}", "GET", None, cache_dir=cache_dir, offline=offline, transport=transport)
    if exact["status"] not in {STATUS_OK, STATUS_EMPTY}:
        return _forward(exact, "HGNC")
    resolution = "exact_symbol"
    response = exact
    if not _hgnc_docs(exact["payload"]):
        alias = _perform("HGNC", f"{_HGNC_ALIAS}/{clean}", "GET", None, cache_dir=cache_dir, offline=offline, transport=transport)
        if alias["status"] not in {STATUS_OK, STATUS_EMPTY}:
            return _forward(alias, "HGNC")
        if not _hgnc_docs(alias["payload"]):
            return _envelope(STATUS_EMPTY, [], source="HGNC", reason="no HGNC record for this symbol or alias")
        resolution = "alias_symbol"
        response = alias
    payload = response["payload"]
    body = payload.get("response") if isinstance(payload, dict) else None
    docs = body.get("docs") if isinstance(body, dict) else None
    if not isinstance(docs, list) or not isinstance(body.get("numFound"), int):
        return _envelope(STATUS_INVALID_RESPONSE, [], source="HGNC", reason="unrecognized HGNC payload")
    if not docs:
        return _envelope(STATUS_EMPTY, [], source="HGNC", reason="no HGNC record for this symbol")
    if len(docs) != 1:
        return _envelope(
            STATUS_AMBIGUOUS,
            [],
            source="HGNC",
            reason=f"{len(docs)} HGNC records share this reported symbol or alias; identity must be resolved manually",
        )
    doc = docs[0]
    record = {
        "REPORTED_SYMBOL": clean,
        "CANONICAL_SYMBOL": _text(doc.get("symbol")).upper(),
        "RESOLUTION": resolution,
        "ALIAS_SYMBOLS": [str(item).upper() for item in doc.get("alias_symbol") or []],
        "PREVIOUS_SYMBOLS": [str(item).upper() for item in doc.get("prev_symbol") or []],
        "ENSEMBL_GENE_ID": _text(doc.get("ensembl_gene_id")),
        "UNIPROT_IDS": [_text(item) for item in doc.get("uniprot_ids") or []],
        "SOURCE_RECORD_ID": "hgnc:" + _text(doc.get("hgnc_id")),
    }
    return _envelope(
        STATUS_OK,
        [record],
        source="HGNC",
        complete=True,
        request_id=response["request_id"],
        response_sha256=response["response_sha256"],
        retrieved_at=response["retrieved_at"],
        reason="" if resolution == "exact_symbol" else f"resolved through the alias index: {clean} -> {record['CANONICAL_SYMBOL']}",
    )


def fetch_protein_evidence(
    symbol: str,
    cache_dir: Path,
    *,
    offline: bool = False,
    transport: Transport | None = None,
) -> dict[str, Any]:
    """Reviewed human UniProt entry for one exact gene symbol, plus QuickGO annotation."""
    clean = _text(symbol).upper()
    if re.fullmatch(r"[A-Z0-9][A-Z0-9\-\.]{0,19}", clean) is None:
        return _envelope(STATUS_UNAVAILABLE, [], source="UNIPROT", reason="a gene symbol is required")
    response = _perform(
        "UNIPROT",
        _UNIPROT,
        "GET",
        {"query": f"gene_exact:{clean} AND organism_id:9606 AND reviewed:true", "format": "json", "size": UNIPROT_PAGE_SIZE},
        cache_dir=cache_dir,
        offline=offline,
        transport=transport,
    )
    if response["status"] != STATUS_OK:
        return _forward(response, "UNIPROT")
    payload = response["payload"]
    results = payload.get("results") if isinstance(payload, dict) else None
    if not isinstance(results, list):
        return _envelope(STATUS_INVALID_RESPONSE, [], source="UNIPROT", reason="unrecognized UniProt payload")
    if not results:
        return _envelope(STATUS_EMPTY, [], source="UNIPROT", reason="no reviewed human UniProt entry for this exact gene symbol")
    exact = [entry for entry in results if isinstance(entry, dict) and _primary_symbol(entry) == clean]
    candidates = exact or results
    if len(candidates) != 1:
        return _envelope(
            STATUS_AMBIGUOUS,
            [],
            source="UNIPROT",
            reason=f"{len(candidates)} reviewed human records remain after exact primary-symbol filtering; identity stays unresolved",
        )
    entry = candidates[0]
    accession = _text(entry.get("primaryAccession"))
    annotations, annotations_complete, _ = _quickgo_annotations_by_accession(accession, cache_dir, offline, transport)
    record = {
        "TARGET_SYMBOL": clean,
        "PRIMARY_SYMBOL": _primary_symbol(entry),
        "UNIPROT_ACCESSION": accession,
        "REVIEWED": True,
        "PROTEIN_ROLE": " ".join(_comments(entry, "FUNCTION")),
        "ROLE_CONTEXT": " ".join(_comments(entry, "SUBUNIT") + _comments(entry, "SUBCELLULAR LOCATION")),
        "MODIFICATION_ANNOTATION": " ".join(_comments(entry, "PTM")),
        "DISEASE_COMMENTS": " ".join(_comments(entry, "DISEASE")),
        "GO_ANNOTATIONS": annotations,
        "SPECIES": "Homo sapiens",
        "SOURCE_RELEASE": _text(entry.get("entryAudit", {}).get("lastAnnotationUpdateDate") if isinstance(entry.get("entryAudit"), dict) else ""),
        "SOURCE_RECORD_ID": accession,
    }
    complete = annotations_complete and len(results) < UNIPROT_PAGE_SIZE
    reason = "" if complete else "the query page cap was reached or GO annotation paging was incomplete"
    return _envelope(
        STATUS_OK,
        [record],
        source="UNIPROT",
        complete=complete,
        request_id=response["request_id"],
        response_sha256=response["response_sha256"],
        retrieved_at=response["retrieved_at"],
        reason=reason,
    )


def _quickgo_annotations_by_accession(
    accession: str,
    cache_dir: Path,
    offline: bool,
    transport: Transport | None,
) -> tuple[list[dict[str, Any]], bool, dict[str, Any]]:
    if not accession:
        return [], False, {"status": STATUS_UNAVAILABLE, "reason": "no accession", "request_id": "", "response_sha256": "", "retrieved_at": _now()}
    rows, complete, last = _go_pages(
        f"{_QUICKGO}/annotation/search",
        {"geneProductId": f"UniProtKB:{accession}", "taxonId": 9606, "limit": PAGE_SIZE},
        cache_dir,
        offline,
        transport,
    )
    records = []
    for row in rows:
        code = _text(row.get("goEvidence"))
        records.append(
            {
                "GO_ID": _text(row.get("goId")),
                "EVIDENCE_CODE": code,
                "EVIDENCE_CLASS": _evidence_class(code),
                "QUALIFIER": _text(row.get("qualifier")),
                "SOURCE_RECORD_ID": _text(row.get("id")),
            }
        )
    return records, complete, last


# --------------------------------------------------------------------------- literature


def fetch_literature(
    symbol: str,
    cache_dir: Path,
    *,
    offline: bool = False,
    transport: Transport | None = None,
    recent_limit: int = 5,
) -> dict[str, Any]:
    """Bounded PubMed search per candidate: total hits are recorded, results are capped."""
    clean = _text(symbol).upper()
    if re.fullmatch(r"[A-Z0-9][A-Z0-9\-\.]{0,19}", clean) is None:
        return _envelope(STATUS_UNAVAILABLE, [], source="PUBMED", reason="a gene symbol is required")
    queries = (
        f"{clean}[gene] AND (mosaic variegated aneuploidy OR chromosome segregation OR spindle assembly checkpoint)",
        f"{clean}[gene] AND (loss of function OR gain of function OR inhibition OR depletion OR activation OR rescue)",
    )
    records: list[dict[str, Any]] = []
    seen: set[str] = set()
    totals = 0
    retrieved = 0
    complete = True
    last: dict[str, Any] = {"status": STATUS_OK, "reason": "", "request_id": "", "response_sha256": "", "retrieved_at": _now()}
    for query in queries:
        search = _perform(
            "PUBMED",
            f"{_NCBI}/esearch.fcgi",
            "GET",
            {"db": "pubmed", "term": query, "retmode": "json", "retmax": recent_limit, "sort": "pub+date"},
            cache_dir=cache_dir,
            offline=offline,
            transport=transport,
        )
        last = search
        if search["status"] != STATUS_OK:
            return _forward(search, "PUBMED")
        body = search["payload"].get("esearchresult") if isinstance(search["payload"], dict) else None
        if not isinstance(body, dict) or not isinstance(body.get("idlist"), list):
            return _envelope(STATUS_INVALID_RESPONSE, [], source="PUBMED", reason="unrecognized PubMed esearch payload")
        try:
            totals += int(body.get("count") or 0)
        except (TypeError, ValueError):
            complete = False
        ids = [str(item) for item in body["idlist"] if str(item).isdigit()]
        retrieved += len(ids)
        if len(ids) >= recent_limit and int(body.get("count") or 0) > len(ids):
            complete = False
        fresh = [item for item in ids if item not in seen]
        seen.update(fresh)
        if not fresh:
            continue
        summary = _perform(
            "PUBMED",
            f"{_NCBI}/esummary.fcgi",
            "GET",
            {"db": "pubmed", "id": ",".join(fresh), "retmode": "json"},
            cache_dir=cache_dir,
            offline=offline,
            transport=transport,
        )
        last = summary
        if summary["status"] != STATUS_OK:
            return _forward(summary, "PUBMED")
        result = summary["payload"].get("result") if isinstance(summary["payload"], dict) else None
        if not isinstance(result, dict):
            return _envelope(STATUS_INVALID_RESPONSE, [], source="PUBMED", reason="unrecognized PubMed esummary payload")
        for pmid in fresh:
            doc = result.get(pmid)
            if not isinstance(doc, dict):
                complete = False
                continue
            doi = ""
            for identifier in doc.get("articleids") or []:
                if isinstance(identifier, dict) and identifier.get("idtype") == "doi":
                    doi = _text(identifier.get("value"))
                    break
            year = re.search(r"\b(19|20)\d{2}\b", _text(doc.get("pubdate")))
            records.append(
                {
                    "PMID": pmid,
                    "TITLE": _text(doc.get("title")),
                    "JOURNAL": _text(doc.get("source")),
                    "YEAR": year.group(0) if year else "",
                    "DOI": doi,
                    "QUERY": query,
                    "SEARCH_COMPLETE": complete,
                    "TOTAL_HITS": totals,
                    "SOURCE_RECORD_ID": f"pubmed:{pmid}",
                }
            )
    status = STATUS_OK if records else STATUS_EMPTY
    reason = "" if complete else "bounded search: the retrieved count is capped below total hits, so this is not a systematic review"
    return _envelope(
        status,
        records,
        source="PUBMED",
        complete=complete,
        request_id=last["request_id"],
        response_sha256=last["response_sha256"],
        retrieved_at=last["retrieved_at"],
        reason=reason,
    )


# --------------------------------------------------------------------------- pharmacology


def _chembl_pages(
    endpoint: str,
    params: dict[str, Any],
    cache_dir: Path,
    offline: bool,
    transport: Transport | None,
) -> tuple[list[dict[str, Any]], bool, dict[str, Any]]:
    """Follow ChEMBL's page_meta.next markers until the collection is exhausted."""
    rows: list[dict[str, Any]] = []
    offset = 0
    complete = True
    last: dict[str, Any] = {"status": STATUS_INVALID_RESPONSE, "reason": "", "request_id": "", "response_sha256": "", "retrieved_at": _now()}
    seen_offsets: set[int] = set()
    for _ in range(MAX_PAGES):
        if offset in seen_offsets:
            complete = False
            break
        seen_offsets.add(offset)
        response = _perform("CHEMBL", f"{_CHEMBL}/{endpoint}", "GET", {**params, "limit": PAGE_SIZE, "offset": offset}, cache_dir=cache_dir, offline=offline, transport=transport)
        last = response
        if response["status"] != STATUS_OK:
            return rows, False, response
        payload = response["payload"]
        if not isinstance(payload, dict):
            return rows, False, {**response, "status": STATUS_INVALID_RESPONSE, "reason": "unrecognized ChEMBL payload"}
        collection = payload.get(_CHEMBL_COLLECTIONS.get(endpoint, ""))
        meta = payload.get("page_meta")
        if not isinstance(collection, list) or not isinstance(meta, dict):
            return rows, False, {**response, "status": STATUS_INVALID_RESPONSE, "reason": "unrecognized ChEMBL collection payload"}
        rows.extend(item for item in collection if isinstance(item, dict))
        next_marker = meta.get("next")
        if not next_marker:
            break
        try:
            offset = int(str(next_marker).split("offset=")[-1].split("&")[0])
        except (TypeError, ValueError):
            complete = False
            break
    else:
        complete = False
    return rows, complete, last


def fetch_drug_evidence(
    symbol: str,
    cache_dir: Path,
    *,
    offline: bool = False,
    transport: Transport | None = None,
) -> dict[str, Any]:
    """ChEMBL and Open Targets records for one verified target identity.

    The Ensembl identifier comes from the reviewed HGNC record; fuzzy target-name
    matching is never used.  Potency is reported only from native assay records, and
    absent assays stay ASSAY_UNAVAILABLE rather than implying support.
    """
    clean = _text(symbol).upper()
    identity = fetch_gene_identity(clean, cache_dir, offline=offline, transport=transport)
    if identity["status"] != STATUS_OK:
        return _forward(identity, "CHEMBL")
    gene = identity["records"][0]
    ensembl = gene["ENSEMBL_GENE_ID"]
    accession = (gene["UNIPROT_IDS"] or [""])[0]
    if not ensemble_valid(ensembl):
        return _envelope(STATUS_AMBIGUOUS, [], source="HGNC", reason="no authoritative Ensembl gene identifier for this symbol", complete=False)
    targets, targets_complete, last = _chembl_pages(
        "target.json",
        {"target_components__accession": accession, "organism": "Homo sapiens"} if accession else {"target_components__gene_name": clean, "organism": "Homo sapiens"},
        cache_dir,
        offline,
        transport,
    )
    admitted = [item for item in targets if _target_matches(item, accession)] if accession else []
    if not admitted:
        return _envelope(
            STATUS_EMPTY if targets_complete else STATUS_UNAVAILABLE,
            [],
            source="CHEMBL",
            complete=targets_complete,
            request_id=last["request_id"],
            response_sha256=last["response_sha256"],
            retrieved_at=last["retrieved_at"],
            reason="no human ChEMBL target matched the reviewed accession" if targets_complete else "target lookup was incomplete",
        )
    records: list[dict[str, Any]] = []
    complete = targets_complete
    for target in admitted:
        target_chembl_id = _text(target.get("target_chembl_id"))
        mechanisms, mechanisms_complete, mech_last = _chembl_pages(
            "mechanism.json",
            {"target_chembl_id": target_chembl_id},
            cache_dir,
            offline,
            transport,
        )
        complete = complete and mechanisms_complete
        for mechanism in mechanisms:
            records.append(
                {
                    "TARGET_SYMBOL": clean,
                    "TARGET_CHEMBL_ID": target_chembl_id,
                    "ENSEMBL_GENE_ID": ensembl,
                    "DRUG_ID": _text(mechanism.get("molecule_chembl_id")),
                    "ACTION_TYPE": _text(mechanism.get("action_type")),
                    "MECHANISM": _text(mechanism.get("mechanism_of_action")),
                    "MAX_PHASE": mechanism.get("max_phase"),
                    "SOURCE": "ChEMBL",
                    "SOURCE_RECORD_ID": _text(mechanism.get("mechanism_of_action")) and f"chembl-mechanism:{mechanism.get('mec_id')}",
                    "ASSAY_STATUS": "ASSAY_UNAVAILABLE",
                    "POTENCY": [],
                }
            )
        if mechanisms:
            activities, activities_complete, _ = _chembl_pages(
                "activity.json",
                {"target_chembl_id": target_chembl_id},
                cache_dir,
                offline,
                transport,
            )
            complete = complete and activities_complete
            if not activities:
                continue
            by_molecule: dict[str, list[dict[str, Any]]] = {}
            for activity in activities:
                by_molecule.setdefault(_text(activity.get("molecule_chembl_id")), []).append(
                    {
                        "STANDARD_TYPE": _text(activity.get("standard_type")),
                        "STANDARD_RELATION": _text(activity.get("standard_relation")),
                        "STANDARD_VALUE": activity.get("standard_value"),
                        "STANDARD_UNITS": _text(activity.get("standard_units")),
                        "ASSAY_TYPE": _text(activity.get("assay_type")),
                        "ASSAY_CONFIDENCE": activity.get("confidence_score"),
                        "DOCUMENT_ID": _text(activity.get("document_chembl_id")),
                    }
                )
            for record in records:
                if record["TARGET_CHEMBL_ID"] != target_chembl_id:
                    continue
                measurements = by_molecule.get(record["DRUG_ID"])
                if measurements:
                    record["POTENCY"] = measurements
                    record["ASSAY_STATUS"] = "ASSAY_RECORDED"
    records.extend(_open_targets_records(clean, ensembl, cache_dir, offline, transport))
    return _envelope(
        STATUS_OK if records else STATUS_EMPTY,
        records,
        source="CHEMBL+OPEN_TARGETS",
        complete=complete,
        request_id=last["request_id"],
        response_sha256=last["response_sha256"],
        retrieved_at=last["retrieved_at"],
        reason="" if complete else "one or more source collections were paged incompletely",
    )


def _target_matches(target: dict[str, Any], accession: str) -> bool:
    if not accession:
        return False
    for component in target.get("target_components") or []:
        if isinstance(component, dict) and accession in (component.get("accession") or ""):
            return True
    return False


def ensemble_valid(identifier: str) -> bool:
    return re.fullmatch(r"ENSG\d{11}", _text(identifier)) is not None


def _open_targets_records(
    symbol: str,
    ensembl: str,
    cache_dir: Path,
    offline: bool,
    transport: Transport | None,
) -> list[dict[str, Any]]:
    query = (
        "query Target($id: String!) { target(ensemblId: $id) { id approvedSymbol "
        "drugAndClinicalCandidates { rows { id maxClinicalStage drug { id name drugType "
        "mechanismsOfAction { rows { mechanismOfAction actionType } } drugWarnings { warningType description } } } } } }"
    )
    response = _perform(
        "OPEN_TARGETS",
        _OPEN_TARGETS,
        "POST",
        {"query": query, "variables": {"id": ensembl}},
        cache_dir=cache_dir,
        offline=offline,
        transport=transport,
    )
    if response["status"] != STATUS_OK:
        return []
    payload = response["payload"]
    data = payload.get("data") if isinstance(payload, dict) else None
    target = data.get("target") if isinstance(data, dict) else None
    if not isinstance(target, dict):
        return []
    records: list[dict[str, Any]] = []
    rows = (target.get("drugAndClinicalCandidates") or {}).get("rows") or []
    for row in rows:
        if not isinstance(row, dict):
            continue
        drug = row.get("drug") or {}
        mechanisms = ((drug.get("mechanismsOfAction") or {}).get("rows")) or []
        action = _text(mechanisms[0].get("actionType")) if mechanisms and isinstance(mechanisms[0], dict) else ""
        warnings = [str(item.get("warningType")) for item in drug.get("drugWarnings") or [] if isinstance(item, dict)]
        records.append(
            {
                "TARGET_SYMBOL": symbol,
                "DRUG_ID": _text(drug.get("id")),
                "DRUG_NAME": _text(drug.get("name")),
                "ACTION_TYPE": action,
                "MAX_PHASE": _text(row.get("maxClinicalStage")),
                "WARNINGS": warnings,
                "SOURCE": "Open Targets",
                "SOURCE_RECORD_ID": _text(row.get("id")),
                "ASSAY_STATUS": "ASSAY_UNAVAILABLE",
                "POTENCY": [],
            }
        )
    return records
