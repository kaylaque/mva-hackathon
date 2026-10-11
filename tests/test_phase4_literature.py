"""Unit tests for Phase 4 Literature Adjudication Engine.

Verifies Crossref DOI resolution, NCBI E-utilities PMID verification,
Europe PMC contradiction detection, and in-memory caching.
Mocked offline without live outbound HTTP calls.
"""

from __future__ import annotations

from unittest.mock import MagicMock, patch

import pytest

from mva_hackathon.literature_adjudication import (
    DOI_REGEX,
    PMID_REGEX,
    CrossrefClient,
    LiteratureAdjudicationConfig,
    LiteratureAdjudicationEngine,
    NcbiEUtilsClient,
    VerifiedCitation,
)


def test_doi_and_pmid_regex() -> None:
    """Validate regex conformity for DOIs and PMIDs."""
    assert DOI_REGEX.match("10.1038/ng1449")
    assert DOI_REGEX.match("10.15252/embj.2020106536")
    assert not DOI_REGEX.match("not_a_doi")
    assert not DOI_REGEX.match("http://doi.org/10.1000/182")

    assert PMID_REGEX.match("15475955")
    assert PMID_REGEX.match("34009673")
    assert not PMID_REGEX.match("PMID15475955")
    assert not PMID_REGEX.match("15475955a")


def test_crossref_mock_resolution() -> None:
    """Mock Crossref REST response: 200 matches -> VERIFIED, 404 -> PENDING."""
    config = LiteratureAdjudicationConfig(offline_mode=False, verify_literature=True)
    client = CrossrefClient(config)

    mock_sess = MagicMock()

    # Success case: 200 OK
    mock_resp_200 = MagicMock()
    mock_resp_200.status_code = 200
    mock_resp_200.json.return_value = {
        "message": {
            "title": ["Constitutional aneuploidy and cancer predisposition caused by biallelic mutations in BUB1B."],
            "container-title": ["Nature Genetics"],
            "issued": {"date-parts": [[2004, 10, 10]]},
        }
    }
    mock_sess.get.return_value = mock_resp_200

    cit_success = client.resolve_doi("10.1038/ng1449", session=mock_sess)
    assert cit_success.status == "VERIFIED"
    assert "Nature Genetics" in cit_success.journal
    assert cit_success.year == "2004"
    assert cit_success.verification_source == "CROSSREF"

    # Not found case: 404
    mock_resp_404 = MagicMock()
    mock_resp_404.status_code = 404
    mock_sess.get.return_value = mock_resp_404
    mock_sess.head.return_value = MagicMock(status_code=404)

    cit_fail = client.resolve_doi("10.1038/nonexistent_doi_123", session=mock_sess)
    assert cit_fail.status == "PENDING/UNVERIFIED"
    assert cit_fail.verification_source == "CROSSREF_NOT_FOUND"


def test_ncbi_eutils_mock_verification() -> None:
    """Mock NCBI E-utilities esummary response: 200 matches -> VERIFIED, empty -> PENDING."""
    config = LiteratureAdjudicationConfig(offline_mode=False, verify_literature=True)
    client = NcbiEUtilsClient(config)

    mock_sess = MagicMock()

    # Valid PMID return
    mock_resp_200 = MagicMock()
    mock_resp_200.status_code = 200
    mock_resp_200.json.return_value = {
        "result": {
            "15475955": {
                "title": "Constitutional aneuploidy and cancer predisposition caused by biallelic mutations in BUB1B.",
                "source": "Nat Genet",
                "pubdate": "2004 Nov",
            }
        }
    }
    mock_sess.get.return_value = mock_resp_200

    cit = client.verify_pmid("15475955", session=mock_sess)
    assert cit.status == "VERIFIED"
    assert "Nat Genet" in cit.journal
    assert cit.year == "2004"
    assert cit.verification_source == "NCBI_EUTILS"

    # Empty / nonexistent PMID
    mock_resp_empty = MagicMock()
    mock_resp_empty.status_code = 200
    mock_resp_empty.json.return_value = {"result": {}}
    mock_sess.get.return_value = mock_resp_empty

    cit_empty = client.verify_pmid("999999999", session=mock_sess)
    assert cit_empty.status == "PENDING/UNVERIFIED"


def test_europe_pmc_contradiction_detection() -> None:
    """Mock Europe PMC REST search detecting literature contradiction."""
    from mva_hackathon.literature_adjudication import EuropePmcClient
    config = LiteratureAdjudicationConfig(offline_mode=False, verify_literature=True)
    client = EuropePmcClient(config)

    mock_sess = MagicMock()

    # Hit with contradiction term in title
    mock_resp = MagicMock()
    mock_resp.status_code = 200
    mock_resp.json.return_value = {
        "resultList": {
            "result": [
                {
                    "title": "Compound X accelerates chromosome missegregation and is contraindicated in MVA.",
                    "pmid": "88881234",
                }
            ]
        }
    }
    mock_sess.get.return_value = mock_resp

    is_contra, pmid, note = client.check_contradiction("BUB1B", "COMPOUND_X", session=mock_sess)
    assert is_contra is True
    assert pmid == "88881234"
    assert "contraindicated" in note.lower()


def test_preloaded_registry_cache() -> None:
    """Preloaded registry citations must immediately return VERIFIED from cache."""
    config = LiteratureAdjudicationConfig(offline_mode=True)
    engine = LiteratureAdjudicationEngine(config)

    # BUB1B primary PMID 15475955 is preloaded from literature_registry.json
    cit = engine.verify_pmid("15475955")
    assert cit.status == "VERIFIED"
    assert cit.verification_source == "REGISTRY_PREVERIFIED"

    # BUB1B primary DOI 10.1038/ng1449 is preloaded
    cit_doi = engine.verify_doi("10.1038/ng1449")
    assert cit_doi.status == "VERIFIED"
    assert cit_doi.verification_source == "REGISTRY_PREVERIFIED"
