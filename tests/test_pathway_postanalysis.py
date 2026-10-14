"""Regression tests for offline DAVID and g:Profiler post-analysis."""

import re
import sys
from pathlib import Path


sys.path.insert(0, str(Path(__file__).parents[1] / "src"))



def test_same_title_different_stable_ids_stay_separate() -> None:
    import pandas as pd

    from mva_hackathon.pathway_postanalysis import _build_evidence

    p1 = pd.DataFrame([
        {"SOURCE": "GO:BP", "TERM_ID": "GO:0000001", "TERM_NAME": "one shared title", "BH_Q": 0.01, "CONTRIBUTING_GENES": "BUB1B"},
        {"SOURCE": "GO:BP", "TERM_ID": "GO:0000002", "TERM_NAME": "one shared title", "BH_Q": 0.01, "CONTRIBUTING_GENES": "BUB1B"},
    ])

    evidence, _ = _build_evidence(p1, pd.DataFrame(columns=["Category", "Term", "Benjamini"]), {"BUB1B": "ENSG_A"}, {"ENSG_A"})

    assert set(evidence["TERM_ID"]) == {"GO:0000001", "GO:0000002"}
    assert evidence["EVIDENCE_ID"].is_unique


def _imported_fixture(root: Path, *, empty: bool = False) -> tuple[Path, Path]:
    """Build a synthetic baseline plus a verified empty-or-populated DAVID import."""
    import hashlib
    import json

    import pandas as pd

    from mva_hackathon.david_workflow import run_david_addition

    baseline = root / "baseline"
    root.mkdir(parents=True, exist_ok=True)
    baseline.mkdir()
    pd.DataFrame([
        {"SYMBOL": "BUB1B", "ENSEMBL_ID": "ENSG_A", "STATUS": "mapped", "IN_FOREGROUND": True, "PRIVATE": "never-export"},
        {"SYMBOL": "TRIP13", "ENSEMBL_ID": "ENSG_B", "STATUS": "mapped", "IN_FOREGROUND": True, "PRIVATE": "never-export"},
    ]).to_csv(baseline / "p0_gene_qc.csv", index=False)
    pd.DataFrame([
        {"SOURCE": "GO:BP", "TERM_ID": "GO:0000001", "TERM_NAME": "g:Profiler title", "BH_Q": 0.04, "CONTRIBUTING_GENES": "BUB1B"},
        {"SOURCE": "GO:BP", "TERM_ID": "GO:0000001", "TERM_NAME": "g:Profiler title", "BH_Q": 0.04, "CONTRIBUTING_GENES": "TRIP13"},
    ]).to_csv(baseline / "p1_enrichment_results.csv", index=False)
    david = root / "david"
    assert run_david_addition(baseline, david, None) == 2
    manifest = json.loads((david / "david_manifest.json").read_text())
    exports = root / "exports"
    exports.mkdir()
    chart, clusters = exports / "chart.csv", exports / "clusters.csv"
    if empty:
        chart.write_text("Category,Term,Count,P-Value,Benjamini,User Ids\n", encoding="utf-8")
        clusters.write_text("Cluster,Cluster Enrichment Score,Category,Term,Full Term,Count,P-Value,Benjamini\n", encoding="utf-8")
    else:
        chart.write_text("Category,Term,Count,P-Value,Benjamini,User Ids\nGOTERM_BP_DIRECT,GO:0000001~DAVID title,2,0.01,0.04,ENSG_A\nGOTERM_BP_DIRECT,GO:0000001~DAVID title,2,0.01,0.04,ENSG_B\n", encoding="utf-8")
        clusters.write_text("Cluster,Cluster Enrichment Score,Category,Term,Full Term,Count,P-Value,Benjamini\n1,2.1,GOTERM_BP_DIRECT,DAVID title,GO:0000001~DAVID title,2,0.01,0.04\n1,2.1,GOTERM_BP_DIRECT,DAVID title,GO:0000001~DAVID title,2,0.01,0.04\n", encoding="utf-8")
    digest = lambda path: hashlib.sha256(path.read_bytes()).hexdigest()
    receipt = {
        "schema_version": 1, "run_key": manifest["run_key"], "list_sha256": manifest["list_sha256"], "namespace": "ENSEMBL_GENE_ID", "taxon": 9606,
        "population_confirmed": True, "timestamps": {"exported_at": "2026-10-13T00:00:00Z"}, "selected_population": "custom",
        "foreground_mapping": {"submitted": 2, "mapped": 2, "unique": 2}, "background_mapping": {"mapped": 2, "reason": "tested"},
        "interface": {"selected_categories": ["GOTERM_BP_DIRECT"], "settings": {"chart": {"Count": 2, "EASE": 0.1}, "clustering": {"displayed": "Medium", "actual_categories_and_thresholds": {"overlap": 3}}}},
        "reports": {"chart": {"complete": True, "complete_unfiltered": True, "empty": empty}, "clusters": {"complete": True, "complete_unfiltered": True, "empty": empty}},
        "export_sha256": {"chart.csv": digest(chart), "clusters.csv": digest(clusters)},
    }
    (exports / "receipt.json").write_text(json.dumps(receipt), encoding="utf-8")
    assert run_david_addition(baseline, david, exports) == 0
    return baseline, david


def test_empty_confirmed_exports_publish_header_only_artifacts(tmp_path: Path) -> None:
    import pandas as pd

    from mva_hackathon.pathway_postanalysis import analyze_pathway_results

    baseline, david = _imported_fixture(tmp_path, empty=True)

    result = analyze_pathway_results(baseline, david, tmp_path / "analysis")

    assert result["status"] == "completed"
    assert result["clusters"] == 0
    output = tmp_path / "analysis"
    evidence = pd.read_csv(output / "integrated_pathway_evidence.csv")
    # The g:Profiler term is still reported; its DAVID counterpart is explicitly not comparable.
    assert list(evidence["MATCH_STATUS"]) == ["gprofiler_only"]
    assert list(evidence["DAVID_SIGNIFICANCE"]) == ["not_comparable"]
    cluster_frame = pd.read_csv(output / "david_cluster_summary.csv")
    assert cluster_frame.empty and len(cluster_frame.columns) > 0
    support = pd.read_csv(output / "gene_pathway_support.csv")
    assert set(support["GENE_SYMBOL"]) == {"BUB1B", "TRIP13"}
    assert set(support["PROVIDER"]) == {"g:Profiler"}
    report = (output / "integrated_pathway_report.md").read_text()
    assert "# Integrated pathway report" in report
    assert not re.search(r"\bnan\b", report.lower())
    assert "none: no DAVID row met its own Benjamini threshold" in report


def test_tampered_import_and_mismatched_baseline_are_refused(tmp_path: Path) -> None:
    import json

    import pytest

    from mva_hackathon.pathway_postanalysis import analyze_pathway_results

    baseline, david = _imported_fixture(tmp_path)
    manifest = json.loads((david / "david_manifest.json").read_text())
    version_dir = david / "imports" / manifest["run_key"] / manifest["import_version"]

    (version_dir / "david_chart.csv").write_text("Category,Term,Count,P-Value,Benjamini\nGOTERM_BP_DIRECT,x,1,0.1,0.1\n", encoding="utf-8")
    with pytest.raises(ValueError, match="hash mismatch"):
        analyze_pathway_results(baseline, david, tmp_path / "tampered")

    baseline, david = _imported_fixture(tmp_path / "second")
    (baseline / "p1_enrichment_results.csv").write_text("SOURCE,TERM_ID,BH_Q,CONTRIBUTING_GENES\nGO:BP,GO:0000009,0.02,BUB1B\n", encoding="utf-8")
    with pytest.raises(ValueError, match="baseline projection mismatch"):
        analyze_pathway_results(baseline, david, tmp_path / "second" / "mismatch")


def test_partial_publish_failure_leaves_no_destination(tmp_path: Path) -> None:
    import pytest

    from mva_hackathon import pathway_postanalysis

    baseline, david = _imported_fixture(tmp_path)
    destination = tmp_path / "analysis"
    staging = tmp_path / ".analysis.staging"

    def _explode(*args: object, **kwargs: object) -> None:
        raise RuntimeError("synthetic writer failure")

    original = pathway_postanalysis._write_report
    pathway_postanalysis._write_report = _explode
    try:
        with pytest.raises(RuntimeError, match="synthetic writer failure"):
            pathway_postanalysis.analyze_pathway_results(baseline, david, destination)
    finally:
        pathway_postanalysis._write_report = original

    assert not destination.exists()
    assert not staging.exists()


def test_repository_baseline_and_symlinked_artifacts_are_refused(tmp_path: Path) -> None:
    import subprocess

    import pytest

    from mva_hackathon.pathway_postanalysis import analyze_pathway_results

    repo = tmp_path / "repo"
    repo.mkdir()
    subprocess.run(["git", "init", "-q"], cwd=repo, check=True)
    loose = repo / "loose_baseline"
    loose.mkdir()
    (loose / "p0_gene_qc.csv").write_text("SYMBOL,ENSEMBL_ID,STATUS,IN_FOREGROUND\n", encoding="utf-8")
    (loose / "p1_enrichment_results.csv").write_text("SOURCE,TERM_ID,BH_Q,CONTRIBUTING_GENES\n", encoding="utf-8")

    with pytest.raises(ValueError, match="under data/processed"):
        analyze_pathway_results(loose, tmp_path / "david", tmp_path / "out")

    baseline, david = _imported_fixture(tmp_path / "synthetic")
    escaped = baseline / "p0_gene_qc.csv"
    escaped.unlink()
    escaped.symlink_to(tmp_path / "outside.csv")
    (tmp_path / "outside.csv").write_text("SYMBOL,ENSEMBL_ID,STATUS,IN_FOREGROUND\nBUB1B,ENSG_A,mapped,True\n", encoding="utf-8")

    with pytest.raises(ValueError, match="invalid baseline artifact path"):
        analyze_pathway_results(baseline, david, tmp_path / "synthetic" / "out")

    import hashlib
    import json

    import pandas as pd

    from mva_hackathon.david_workflow import run_david_addition
    from mva_hackathon.pathway_postanalysis import analyze_pathway_results

    baseline = tmp_path / "baseline"
    baseline.mkdir()
    pd.DataFrame([
        {"SYMBOL": "BUB1B", "ENSEMBL_ID": "ENSG_A", "STATUS": "mapped", "IN_FOREGROUND": True, "PRIVATE": "never-export"},
        {"SYMBOL": "TRIP13", "ENSEMBL_ID": "ENSG_B", "STATUS": "mapped", "IN_FOREGROUND": True, "PRIVATE": "never-export"},
    ]).to_csv(baseline / "p0_gene_qc.csv", index=False)
    pd.DataFrame([
        {"SOURCE": "GO:BP", "TERM_ID": "GO:0000001", "TERM_NAME": "g:Profiler title", "BH_Q": 0.04, "PROVIDER_ADJUSTED_P": 0.04, "PROVIDER_CORRECTION": "BH", "RAW_P": 0.01, "FOLD_ENRICHMENT": 2, "CONTRIBUTING_GENES": "BUB1B"},
        {"SOURCE": "GO:BP", "TERM_ID": "GO:0000001", "TERM_NAME": "g:Profiler title", "BH_Q": 0.04, "PROVIDER_ADJUSTED_P": 0.04, "PROVIDER_CORRECTION": "BH", "RAW_P": 0.01, "FOLD_ENRICHMENT": 2, "CONTRIBUTING_GENES": "TRIP13"},
    ]).to_csv(baseline / "p1_enrichment_results.csv", index=False)
    david = tmp_path / "david"
    assert run_david_addition(baseline, david, None) == 2
    manifest = json.loads((david / "david_manifest.json").read_text())
    exports = tmp_path / "exports"
    exports.mkdir()
    chart = exports / "chart.csv"
    clusters = exports / "clusters.csv"
    chart.write_text("Category,Term,Count,P-Value,Benjamini,User Ids\nGOTERM_BP_DIRECT,GO:0000001~DAVID title,2,0.01,0.04,ENSG_A\nGOTERM_BP_DIRECT,GO:0000001~DAVID title,2,0.01,0.04,ENSG_B\n", encoding="utf-8")
    clusters.write_text("Cluster,Cluster Enrichment Score,Category,Term,Full Term,Count,P-Value,Benjamini\n1,2.1,GOTERM_BP_DIRECT,DAVID title,GO:0000001~DAVID title,2,0.01,0.04\n1,2.1,GOTERM_BP_DIRECT,DAVID title,GO:0000001~DAVID title,2,0.01,0.04\n", encoding="utf-8")
    digest = lambda path: hashlib.sha256(path.read_bytes()).hexdigest()
    receipt = {
        "schema_version": 1, "run_key": manifest["run_key"], "list_sha256": manifest["list_sha256"], "namespace": "ENSEMBL_GENE_ID", "taxon": 9606,
        "population_confirmed": True, "timestamps": {"exported_at": "2026-10-13T00:00:00Z"}, "selected_population": "custom",
        "foreground_mapping": {"submitted": 2, "mapped": 2, "unique": 2}, "background_mapping": {"mapped": 2, "reason": "tested"},
        "interface": {"selected_categories": ["GOTERM_BP_DIRECT"], "settings": {"chart": {"Count": 2, "EASE": 0.1}, "clustering": {"displayed": "Medium", "actual_categories_and_thresholds": {"overlap": 3}}}},
        "reports": {"chart": {"complete": True, "complete_unfiltered": True, "empty": False}, "clusters": {"complete": True, "complete_unfiltered": True, "empty": False}},
        "export_sha256": {"chart.csv": digest(chart), "clusters.csv": digest(clusters)},
    }
    (exports / "receipt.json").write_text(json.dumps(receipt), encoding="utf-8")
    assert run_david_addition(baseline, david, exports) == 0

    result = analyze_pathway_results(baseline, david, tmp_path / "analysis")

    evidence = pd.read_csv(tmp_path / "analysis" / "integrated_pathway_evidence.csv")
    clusters_out = pd.read_csv(tmp_path / "analysis" / "david_cluster_summary.csv")
    assert result["status"] == "completed"
    assert evidence.loc[0, "MATCH_STATUS"] == "shared"
    assert evidence.loc[0, "EVIDENCE_TIER"] == "both_significant"
    assert evidence.loc[0, "GENE_JACCARD"] == 1
    support = pd.read_csv(tmp_path / "analysis" / "gene_pathway_support.csv")
    assert set(support["GENE_SYMBOL"]) == {"BUB1B", "TRIP13"}
    assert clusters_out.loc[0, "ANNOTATION_ROW_COUNT"] == 2
    import subprocess

    script = Path(__file__).parents[1] / "scripts" / "analyze_pathway_results.py"
    cli_output = tmp_path / "cli-analysis"
    command = [sys.executable, str(script), "--reuse-pathway-dir", str(baseline), "--david-dir", str(david), "--output-dir", str(cli_output)]
    first = subprocess.run(command, capture_output=True, text=True, check=False)
    second = subprocess.run(command, capture_output=True, text=True, check=False)
    assert first.returncode == 0
    assert second.returncode == 1
    assert pd.read_csv(cli_output / "integrated_pathway_evidence.csv").equals(evidence)
    assert len(clusters_out) == 1
    assert clusters_out.loc[0, "MAPPED_TERM_COUNT"] == 1
    assert "never-export" not in (tmp_path / "analysis" / "integrated_pathway_evidence.csv").read_text()


def test_cli_requires_all_three_paths() -> None:
    import subprocess

    script = Path(__file__).parents[1] / "scripts" / "analyze_pathway_results.py"
    completed = subprocess.run(
        [sys.executable, str(script)],
        capture_output=True,
        text=True,
        check=False,
    )

    assert completed.returncode == 1
    assert "--reuse-pathway-dir" in completed.stderr


def test_repeated_adjusted_values_are_mixed_without_selecting_a_minimum() -> None:
    from mva_hackathon.pathway_postanalysis import _state, _tier

    state = _state([{"BH_Q": 0.01}, {"BH_Q": 0.2}], "BH_Q")

    assert state == "mixed"
    assert _tier("shared", state, "significant") == "review_required"


def test_existing_or_overlapping_output_is_refused_before_writes(tmp_path: Path) -> None:
    from mva_hackathon.pathway_postanalysis import _validate_output

    baseline = tmp_path / "baseline"
    david = tmp_path / "david"
    baseline.mkdir()
    david.mkdir()
    existing = tmp_path / "existing"
    existing.mkdir()

    import pytest

    with pytest.raises(ValueError, match="already exists"):
        _validate_output(existing, baseline, david)
    with pytest.raises(ValueError, match="disjoint"):
        _validate_output(baseline / "child", baseline, david)


def test_unmapped_significant_rows_remain_distinct_and_reach_discordance_review() -> None:
    import json

    import pandas as pd

    from mva_hackathon.pathway_postanalysis import _build_evidence, _discordance

    p1 = pd.DataFrame([
        {"SOURCE": "GO:BP", "TERM_ID": "native keyword A", "BH_Q": 0.01, "CONTRIBUTING_GENES": "BUB1B"},
        {"SOURCE": "GO:BP", "TERM_ID": "native keyword B", "BH_Q": 0.01, "CONTRIBUTING_GENES": "TRIP13"},
    ])
    chart = pd.DataFrame([
        {"Category": "GOTERM_BP_DIRECT", "Term": "native DAVID keyword", "Benjamini": 0.01, "User Ids": "ENSG_A"},
    ])

    evidence, _ = _build_evidence(
        p1,
        chart,
        {"BUB1B": "ENSG_A", "TRIP13": "ENSG_B"},
        {"ENSG_A", "ENSG_B"},
    )
    review = _discordance(evidence, ["GOTERM_BP_DIRECT"])

    assert len(evidence) == 3
    assert evidence["EVIDENCE_ID"].is_unique
    assert set(evidence["EVIDENCE_TIER"]) == {"review_required"}
    assert set(review["MATCH_STATUS"]) == {"unmapped_term"}
    assert review["UNRESOLVED_MATCH_REASON"].str.contains("stable identifier").all()
    assert {json.loads(value)[0]["TERM_ID"] for value in evidence["GPROFILER_ROWS_JSON"] if value != "[]" } == {
        "native keyword A",
        "native keyword B",
    }


def test_contributor_unions_keep_partial_memberships_without_calculating_jaccard() -> None:
    import pandas as pd

    from mva_hackathon.pathway_postanalysis import _build_evidence

    p1 = pd.DataFrame([
        {"SOURCE": "GO:BP", "TERM_ID": "GO:0000001", "BH_Q": 0.01, "CONTRIBUTING_GENES": "BUB1B"},
        {"SOURCE": "GO:BP", "TERM_ID": "GO:0000001", "BH_Q": 0.01, "CONTRIBUTING_GENES": "UNRESOLVED"},
    ])
    chart = pd.DataFrame([
        {"Category": "GOTERM_BP_DIRECT", "Term": "GO:0000001~term", "Benjamini": 0.01, "User Ids": "ENSG_A"},
        {"Category": "GOTERM_BP_DIRECT", "Term": "GO:0000001~term", "Benjamini": 0.01, "User Ids": "ENSG_B"},
    ])

    evidence, support = _build_evidence(
        p1,
        chart,
        {"BUB1B": "ENSG_A", "TRIP13": "ENSG_B"},
        {"ENSG_A", "ENSG_B"},
    )

    row = evidence.iloc[0]
    assert row["GPROFILER_CONTRIBUTOR_STATUS"] == "partial"
    assert row["DAVID_CONTRIBUTOR_STATUS"] == "complete"
    assert row["CONTRIBUTOR_STATUS"] == "partial"
    assert pd.isna(row["GENE_JACCARD"])
    assert set(support.loc[support["PROVIDER"].eq("g:Profiler"), "MEMBERSHIP_STATUS"]) == {"partial"}
    assert set(support.loc[support["PROVIDER"].eq("DAVID"), "MEMBERSHIP_STATUS"]) == {"observed"}


def test_complete_disjoint_contributors_have_zero_jaccard() -> None:
    import pandas as pd

    from mva_hackathon.pathway_postanalysis import _build_evidence

    evidence, _ = _build_evidence(
        pd.DataFrame([{
            "SOURCE": "GO:BP", "TERM_ID": "GO:0000001", "BH_Q": 0.01, "CONTRIBUTING_GENES": "BUB1B",
        }]),
        pd.DataFrame([{
            "Category": "GOTERM_BP_DIRECT", "Term": "GO:0000001~term", "Benjamini": 0.01, "User Ids": "ENSG_B",
        }]),
        {"BUB1B": "ENSG_A", "TRIP13": "ENSG_B"},
        {"ENSG_A", "ENSG_B"},
    )

    row = evidence.iloc[0]
    assert row["CONTRIBUTOR_STATUS"] == "complete"
    assert row["GENE_JACCARD"] == 0


def test_missing_adjusted_value_does_not_hide_later_invalid_value() -> None:
    import pytest

    from mva_hackathon.pathway_postanalysis import _state

    with pytest.raises(ValueError, match="BH_Q"):
        _state([{"BH_Q": None}, {"BH_Q": float("inf")}], "BH_Q")


def test_cluster_summary_uses_full_terms_and_cluster_native_significance() -> None:
    import json

    import pandas as pd

    from mva_hackathon.pathway_postanalysis import _cluster_summary

    evidence = pd.DataFrame([{
        "EVIDENCE_ID": "ev-1",
        "ONTOLOGY": "GO",
        "TERM_ID": "GO:0000001",
        "MATCH_STATUS": "shared",
        "DAVID_SIGNIFICANCE": "significant",
        "MVA_CONTEXT_JSON": json.dumps(["SAC_CORE_GENES"]),
    }])
    clusters = pd.DataFrame([
        {
            "Cluster": "1", "Cluster Enrichment Score": 2.1, "Category": "GOTERM_BP_DIRECT",
            "Term": "display one", "Full Term": "GO:0000001~full one", "Benjamini": 0.01, "User Ids": "ENSG_A",
        },
        {
            "Cluster": "1", "Cluster Enrichment Score": 2.1, "Category": "GOTERM_BP_FAT",
            "Term": "display two", "Full Term": "GO:0000001~full two", "Benjamini": 0.2, "User Ids": "ENSG_A",
        },
        {
            "Cluster": "1", "Cluster Enrichment Score": 2.1, "Category": "GOTERM_BP_DIRECT",
            "Term": "unmapped display", "Full Term": "free-text annotation", "Benjamini": 0.01, "User Ids": "",
        },
    ])

    summary = _cluster_summary(
        clusters,
        evidence,
        {"BUB1B": "ENSG_A"},
        {"ENSG_A"},
    )

    row = summary.iloc[0]
    assert row["ANNOTATION_ROW_COUNT"] == 3
    assert row["MAPPED_TERM_COUNT"] == 2
    assert row["UNMAPPED_TERM_COUNT"] == 1
    assert row["SHARED_TERM_COUNT"] == 1
    assert row["DAVID_SIGNIFICANT_TERM_COUNT"] == 1
    assert row["CONTEXT_MEMBERSHIP_SOURCE"] == "cluster_observed"
    assert json.loads(row["CONTEXT_GENES_JSON"]) == ["BUB1B"]


def test_cluster_summary_rejects_infinite_scores() -> None:
    import pandas as pd
    import pytest

    from mva_hackathon.pathway_postanalysis import _cluster_summary

    clusters = pd.DataFrame([{
        "Cluster": "1", "Cluster Enrichment Score": float("inf"), "Category": "GOTERM_BP_DIRECT",
        "Term": "term", "Full Term": "GO:0000001~term", "Benjamini": 0.01,
    }])
    with pytest.raises(ValueError, match="cluster enrichment score"):
        _cluster_summary(clusters, pd.DataFrame(), {}, set())
