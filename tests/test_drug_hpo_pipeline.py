from __future__ import annotations

import json
import threading
import time

from pathlib import Path

import pandas as pd
import pytest

from mva_hackathon import compass
from mva_hackathon import drug_hpo_evidence as evidence


def _inputs(root: Path) -> tuple[Path, Path]:
    source = root / "compass"
    source.mkdir()
    rationales = []
    drugs = []
    for index, drug_id in enumerate(("MECH", "PHENO", "HAZARD"), start=1):
        rationale = {column: "" for column in compass.RATIONALE_COLUMNS}
        rationale.update({
            "RATIONALE_ID": f"R{index}",
            "EXPLORATORY_DAVID_OBSERVED_TERM_IDS": "hsa05170",
            "DESIRED_TARGET_CHANGE": "DECREASE" if drug_id == "MECH" else "",
            "DISEASE_STATE_DIRECTION": "INCREASED_FUNCTION" if drug_id == "MECH" else "UNKNOWN",
        })
        rationales.append(rationale)
        drug = {column: "" for column in compass.DRUG_REVIEW_COLUMNS}
        drug.update({
            "RATIONALE_ID": f"R{index}", "DRUG_ID": drug_id, "DRUG_NAME": drug_id.title(),
            "TARGET_SYMBOL": "TP53", "ACTION_TYPE": "INHIBITOR",
            "PATHWAY_EVIDENCE_LEVEL": "RELAXED_DAVID_BENJAMINI_0_90", "DECISION": "HOLD",
        })
        drugs.append(drug)
    pd.DataFrame(rationales, columns=compass.RATIONALE_COLUMNS).to_csv(source / "target_mechanism_rationale.csv", index=False)
    pd.DataFrame(drugs, columns=compass.DRUG_REVIEW_COLUMNS).to_csv(source / "target_drug_review.csv", index=False)
    (source / "compass_provenance.json").write_text(json.dumps({"identity": {"selection_profile": "david_relaxed", "david_benjamini_threshold": 0.9, "drug_scope": "relaxed_pathway"}}))
    hpo = root / "hpo.txt"
    hpo.write_text("HP:0000001\n")
    return source, hpo


def _mock_sources(monkeypatch):
    from mva_hackathon import drug_hpo_sources as sources

    envelope = lambda source_name, records, complete=True: {
        "status": "EMPTY" if not records else "OK", "records": records, "complete": complete,
        "source": source_name, "request_id": "req", "response_sha256": "abc",
        "retrieved_at": "2026-10-25T00:00:00Z", "reason": "",
    }
    monkeypatch.setattr(sources, "resolve_hpo_term", lambda code, *a, **k: envelope("OLS4", [{"HPO_CODE": code, "HPO_LABEL": "Term", "HPO_SYNONYMS": [], "IRI": "iri", "ONTOLOGY_VERSION": "1"}]))
    def literature(drug, *a, **k):
        if drug == "Mech":
            records = [{"PMID": "1", "TITLE": "Mechanism", "ABSTRACT": "Mech inhibits TP53 and reduced activity.", "SOURCE_RECORD_ID": "PMID:1", "PUBLICATION_TYPES": ["Clinical Trial"], "INTERVENTION": "Mech", "MODEL_OR_POPULATION": "human patients"}]
        elif drug == "Pheno":
            records = [{"PMID": "2", "TITLE": "Phenotype", "ABSTRACT": "Pheno improved Term.", "SOURCE_RECORD_ID": "PMID:2", "STUDY_DESIGN": "Case Report", "INTERVENTION": "Pheno", "MODEL_OR_POPULATION": "human patient"}]
        else:
            records = []
        return envelope("LITERATURE", records)
    monkeypatch.setattr(sources, "search_drug_hpo_literature", literature)
    def fda(drug, *a, **k):
        records = [{"SOURCE_RECORD_ID": "FDA:3", "GENERIC_NAMES": ["Hazard"], "BOXED_WARNINGS": "Term is a serious risk."}] if drug == "Hazard" else []
        return envelope("FDA", records)
    monkeypatch.setattr(sources, "fetch_fda_labels", fda)
    monkeypatch.setattr(sources, "fetch_trials", lambda *a, **k: envelope("TRIALS", []))


def test_pipeline_publishes_data_blind_outputs_and_resume_checks_identity(tmp_path, monkeypatch):
    source, hpo = _inputs(tmp_path)
    _mock_sources(monkeypatch)
    output = tmp_path / "result"
    result = evidence.run_drug_hpo_evidence(source, hpo, output, cache_dir=tmp_path / "cache")
    assert result["status"] == "completed"
    assert result["input_drug_rows"] == 3
    assert (result["mechanism_qualified"], result["phenotype_qualified"], result["test_first"]) == (1, 1, 2)
    enriched = pd.read_csv(output / "target_drug_review_enriched.csv", dtype=str, keep_default_na=False).set_index("DRUG_ID")
    assert enriched.loc["MECH", "DECISION"] == "TEST_FIRST"
    assert enriched.loc["MECH", "PROMOTION_ROUTE"] == "MECHANISM_DIRECTED"
    assert enriched.loc["PHENO", "DECISION"] == "TEST_FIRST"
    assert enriched.loc["PHENO", "PROMOTION_ROUTE"] == "PHENOTYPE_DIRECTED"
    assert "DIRECTION_UNKNOWN" in enriched.loc["PHENO", "REASON_CODES"]
    assert enriched.loc["HAZARD", "DECISION"] == "REJECT"
    assert (enriched.loc["MECH", "PRIORITY_BAND"], enriched.loc["PHENO", "PRIORITY_BAND"]) == ("B", "B")
    assert enriched.loc["PHENO", "PRIORITY_RANK"] == "1"
    assert enriched.loc["MECH", "PRIORITY_RANK"] == "2"
    assert len(pd.read_csv(output / "target_drug_review_enriched.csv")) == 3
    assert {p.name for p in output.iterdir()} == {
        "hpo_dictionary.csv", "drug_identity.csv", "search_runs.csv", "public_evidence_claims.csv",
        "drug_hpo_assessments.csv", "target_drug_review_enriched.csv", "target_drug_review_decision_ready.csv",
        "source_status.csv", "drug_hpo_report.md", "drug_hpo_provenance.json",
    }
    report = (output / "drug_hpo_report.md").read_text()
    assert "research-only" in report.lower()
    assert "Term" not in report
    provenance = json.loads((output / "drug_hpo_provenance.json").read_text())
    assert provenance["schema_version"] == "drug-hpo-evidence-1"
    for name, digest in provenance["output_sha256"].items():
        assert evidence.sha256_file(output / name) == digest
    assert evidence.run_drug_hpo_evidence(source, hpo, output, cache_dir=tmp_path / "cache", resume=True, offline=True)["status"] == "completed"
    replay = tmp_path / "replay"
    _mock_sources(monkeypatch)
    replay_result = evidence.run_drug_hpo_evidence(source, hpo, replay, cache_dir=tmp_path / "cache", offline=True)
    assert replay_result["status"] == "completed"
    for name in provenance["output_sha256"]:
        assert (output / name).read_bytes() == (replay / name).read_bytes()


def test_cli_reports_incomplete_sources_and_refuses_mismatched_resume(tmp_path):
    import os
    import subprocess
    import sys

    source, hpo = _inputs(tmp_path)
    output = tmp_path / "cli-output"
    command = [
        sys.executable, "scripts/run_drug_hpo_evidence.py",
        "--compass-dir", str(source), "--hpo-codes", str(hpo),
        "--output-dir", str(output), "--cache-dir", str(tmp_path / "empty-cache"), "--offline",
    ]
    env = {**os.environ, "PYTHONPATH": "src"}
    incomplete = subprocess.run(command, capture_output=True, text=True, env=env)
    assert incomplete.returncode == 2
    assert json.loads(incomplete.stdout)["status"] == "source_incomplete"
    refusal = subprocess.run(command, capture_output=True, text=True, env=env)
    assert refusal.returncode == 1
    changed = tmp_path / "changed-hpo.txt"
    changed.write_text("HP:0000002\n")
    resumed = subprocess.run([*command, "--resume", "--hpo-codes", str(changed)], capture_output=True, text=True, env=env)
    assert resumed.returncode == 1
    assert resumed.stderr.startswith("drug-HPO evidence failed: resumed inputs differ")


def test_cli_subprocess_exits_zero_for_complete_mocked_sources(tmp_path):
    import os
    import subprocess
    import sys

    source, hpo = _inputs(tmp_path)
    bootstrap = tmp_path / "cli_success.py"
    bootstrap.write_text(
        "import runpy, sys\n"
        f"sys.path.insert(0, {str(Path.cwd() / 'src')!r})\n"
        "from mva_hackathon import drug_hpo_sources as sources\n"
        "def envelope(name, records):\n"
        " return {'status':'EMPTY' if not records else 'OK','records':records,'complete':True,'source':name,'request_id':'fixture','response_sha256':'fixture','retrieved_at':'2026-10-25T00:00:00Z','reason':''}\n"
        "sources.resolve_hpo_term=lambda code,*a,**k: envelope('OLS4',[{'HPO_CODE':code,'HPO_LABEL':'Term','HPO_SYNONYMS':[],'IRI':'iri','ONTOLOGY_VERSION':'1'}])\n"
        "sources.search_drug_hpo_literature=lambda *a,**k: envelope('LITERATURE',[])\n"
        "sources.fetch_fda_labels=lambda *a,**k: envelope('FDA',[])\n"
        "sources.fetch_trials=lambda *a,**k: envelope('TRIALS',[])\n"
        f"sys.argv = ['scripts/run_drug_hpo_evidence.py','--compass-dir',{str(source)!r},'--hpo-codes',{str(hpo)!r},'--output-dir',{str(tmp_path / 'cli-success')!r},'--cache-dir',{str(tmp_path / 'cli-cache')!r}]\n"
        f"runpy.run_path({str(Path.cwd() / 'scripts/run_drug_hpo_evidence.py')!r},run_name='__main__')\n"
    )
    env = {**os.environ, "PYTHONPATH": "src"}
    completed = subprocess.run([sys.executable, str(bootstrap)], capture_output=True, text=True, env=env)
    assert completed.returncode == 0, completed.stderr
    summary = json.loads(completed.stdout)
    assert summary["status"] == "completed"
    assert summary["input_drug_rows"] == 3


def test_parallel_retrieval_deduplicates_scopes_preserves_duplicates_and_matches_serial_bytes(tmp_path, monkeypatch):
    source, hpo = _inputs(tmp_path)
    drugs = pd.read_csv(source / "target_drug_review.csv", dtype=str, keep_default_na=False)
    duplicate = drugs.iloc[[0]].assign(TARGET_SYMBOL="BUB1B")
    pd.concat([drugs.iloc[:2], duplicate], ignore_index=True).to_csv(source / "target_drug_review.csv", index=False)
    hpo.write_text("HP:0000001\nHP:0000002\n")

    from mva_hackathon import drug_hpo_sources as sources

    calls: list[tuple[str, str, str]] = []
    active = 0
    peak_active = 0
    lock = threading.Lock()

    def envelope(source_name, records=()):
        return {
            "status": "EMPTY" if not records else "OK",
            "records": list(records),
            "complete": True,
            "source": source_name,
            "request_id": source_name,
            "response_sha256": source_name,
            "retrieved_at": "2026-10-25T00:00:00Z",
            "reason": "",
        }

    def scoped(source_name, key, records=()):
        nonlocal active, peak_active
        with lock:
            calls.append((source_name, *key))
            active += 1
            peak_active = max(peak_active, active)
        time.sleep(0.03)
        with lock:
            active -= 1
        return envelope(source_name, records)

    monkeypatch.setattr(sources, "resolve_hpo_term", lambda code, *a, **k: scoped("OLS4", (code, ""), [{"HPO_CODE": code, "HPO_LABEL": f"Term {code[-1]}", "HPO_SYNONYMS": [], "IRI": "iri", "ONTOLOGY_VERSION": "1"}]))
    monkeypatch.setattr(sources, "fetch_fda_labels", lambda drug, *a, **k: scoped("FDA", (drug, "")))
    monkeypatch.setattr(sources, "search_drug_hpo_literature", lambda drug, _aliases, code, *_a, **_k: scoped("LITERATURE", (drug, code)))
    monkeypatch.setattr(sources, "fetch_trials", lambda drug, _aliases, label, *_a, **_k: scoped("TRIALS", (drug, label[-1])))

    parallel = tmp_path / "parallel"
    evidence.run_drug_hpo_evidence(source, hpo, parallel, cache_dir=tmp_path / "cache", max_workers=4)
    assert peak_active > 1
    assert sorted(calls) == sorted([
        ("OLS4", "HP:0000001", ""), ("OLS4", "HP:0000002", ""),
        ("FDA", "Mech", ""), ("FDA", "Pheno", ""),
        *[("LITERATURE", drug, code) for drug in ("Mech", "Pheno") for code in ("HP:0000001", "HP:0000002")],
        *[("TRIALS", drug, code[-1]) for drug in ("Mech", "Pheno") for code in ("HP:0000001", "HP:0000002")],
    ])
    assert len(pd.read_csv(parallel / "target_drug_review_enriched.csv")) == 3

    calls.clear()
    serial = tmp_path / "serial"
    evidence.run_drug_hpo_evidence(source, hpo, serial, cache_dir=tmp_path / "cache", max_workers=1)
    assert all((serial / name).read_bytes() == (parallel / name).read_bytes() for name in json.loads((parallel / "drug_hpo_provenance.json").read_text())["output_sha256"])
    parallel_provenance = json.loads((parallel / "drug_hpo_provenance.json").read_text())
    serial_provenance = json.loads((serial / "drug_hpo_provenance.json").read_text())
    parallel_provenance.pop("generated_at")
    serial_provenance.pop("generated_at")
    assert parallel_provenance == serial_provenance


def test_cli_rejects_zero_workers_with_validation_message(tmp_path):
    import os
    import subprocess
    import sys

    source, hpo = _inputs(tmp_path)
    completed = subprocess.run(
        [
            sys.executable, "scripts/run_drug_hpo_evidence.py",
            "--compass-dir", str(source), "--hpo-codes", str(hpo),
            "--output-dir", str(tmp_path / "result"), "--cache-dir", str(tmp_path / "cache"),
            "--max-workers", "0",
        ],
        capture_output=True, text=True, env={**os.environ, "PYTHONPATH": "src"},
    )
    assert completed.returncode == 1
    assert completed.stderr == "drug-HPO evidence failed: max_workers must be at least 1\n"


def test_failed_worker_marks_scope_incomplete_without_qualifying_claims(tmp_path, monkeypatch):
    source, hpo = _inputs(tmp_path)
    _mock_sources(monkeypatch)
    from mva_hackathon import drug_hpo_sources as sources

    def fail_literature(drug, *_args, **_kwargs):
        if drug == "Mech":
            raise RuntimeError("network error")
        return {"status": "EMPTY", "records": [], "complete": True, "source": "LITERATURE", "request_id": "", "response_sha256": "", "retrieved_at": "", "reason": ""}

    monkeypatch.setattr(sources, "search_drug_hpo_literature", fail_literature)
    output = tmp_path / "result"
    result = evidence.run_drug_hpo_evidence(source, hpo, output, cache_dir=tmp_path / "cache", max_workers=4)
    assert result["status"] == "source_incomplete"
    statuses = pd.read_csv(output / "source_status.csv", dtype=str, keep_default_na=False)
    failed = statuses.loc[(statuses["SOURCE"] == "LITERATURE") & (statuses["CANONICAL_DRUG_ID"] == "MECH")].iloc[0]
    assert (failed["COMPLETE"], failed["REASON"]) == ("False", "task exception: RuntimeError")
    enriched = pd.read_csv(output / "target_drug_review_enriched.csv", dtype=str, keep_default_na=False).set_index("DRUG_ID")
    assert enriched.loc["MECH", "DECISION"] != "TEST_FIRST"


def test_successful_scope_cannot_supply_claims_for_failed_hpo_scope(tmp_path, monkeypatch):
    source, hpo = _inputs(tmp_path)
    hpo.write_text("HP:0000001\nHP:0000002\n")
    from mva_hackathon import drug_hpo_sources as sources

    def envelope(source_name, records=()):
        return {"status": "EMPTY" if not records else "OK", "records": list(records), "complete": True, "source": source_name, "request_id": source_name, "response_sha256": source_name, "retrieved_at": "", "reason": ""}

    monkeypatch.setattr(sources, "resolve_hpo_term", lambda code, *_a, **_k: envelope("OLS4", [{"HPO_CODE": code, "HPO_LABEL": f"Term {code[-1]}", "HPO_SYNONYMS": [], "IRI": "", "ONTOLOGY_VERSION": "1"}]))
    monkeypatch.setattr(sources, "fetch_fda_labels", lambda *_a, **_k: envelope("FDA"))
    monkeypatch.setattr(sources, "fetch_trials", lambda *_a, **_k: envelope("TRIALS"))

    def literature(drug, _aliases, code, *_args, **_kwargs):
        if drug == "Mech" and code == "HP:0000002":
            raise RuntimeError("failed scope")
        records = [{"PMID": "1", "TITLE": "scope", "ABSTRACT": "Mech improved Term 2.", "SOURCE_RECORD_ID": "PMID:1", "INTERVENTION": "Mech", "MODEL_OR_POPULATION": "human patient"}] if drug == "Mech" else []
        return envelope("LITERATURE", records)

    monkeypatch.setattr(sources, "search_drug_hpo_literature", literature)
    output = tmp_path / "result"
    evidence.run_drug_hpo_evidence(source, hpo, output, cache_dir=tmp_path / "cache", max_workers=4)
    claims = pd.read_csv(output / "public_evidence_claims.csv", dtype=str, keep_default_na=False)
    assert not ((claims["CANONICAL_DRUG_ID"] == "MECH") & (claims["HPO_CODE"] == "HP:0000002")).any()


def test_identifier_only_drug_is_incomplete_without_public_queries(tmp_path, monkeypatch):
    source, hpo = _inputs(tmp_path)
    drugs = pd.read_csv(source / "target_drug_review.csv", dtype=str, keep_default_na=False)
    drugs.loc[0, ["DRUG_ID", "DRUG_NAME"]] = "CHEMBL1"
    drugs.to_csv(source / "target_drug_review.csv", index=False)
    from mva_hackathon import drug_hpo_sources as sources

    calls: list[str] = []
    envelope = lambda source_name, records=(): {"status": "EMPTY", "records": list(records), "complete": True, "source": source_name, "request_id": "", "response_sha256": "", "retrieved_at": "", "reason": ""}
    monkeypatch.setattr(sources, "resolve_hpo_term", lambda code, *_a, **_k: envelope("OLS4", [{"HPO_CODE": code, "HPO_LABEL": "Term", "HPO_SYNONYMS": [], "IRI": "", "ONTOLOGY_VERSION": "1"}]))
    monkeypatch.setattr(sources, "fetch_fda_labels", lambda *_a, **_k: calls.append("FDA") or envelope("FDA"))
    monkeypatch.setattr(sources, "search_drug_hpo_literature", lambda *_a, **_k: calls.append("LITERATURE") or envelope("LITERATURE"))
    monkeypatch.setattr(sources, "fetch_trials", lambda *_a, **_k: calls.append("TRIALS") or envelope("TRIALS"))

    output = tmp_path / "result"
    result = evidence.run_drug_hpo_evidence(source, hpo, output, cache_dir=tmp_path / "cache")

    assert result["status"] == "source_incomplete"
    assert calls.count("FDA") == 2
    statuses = pd.read_csv(output / "source_status.csv", dtype=str, keep_default_na=False)
    unresolved = statuses.loc[statuses["CANONICAL_DRUG_ID"].eq("CHEMBL1")]
    assert unresolved["COMPLETE"].eq("False").all()
    assert set(unresolved["REASON"]) == {"unresolved public drug identity"}
