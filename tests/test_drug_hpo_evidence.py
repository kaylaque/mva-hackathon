from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys

import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).parents[1] / "src"))

from mva_hackathon.compass import DRUG_REVIEW_COLUMNS, RATIONALE_COLUMNS
from mva_hackathon.drug_hpo_evidence import load_compass_snapshot, load_hpo_codes


def test_load_hpo_codes_deduplicates_and_hashes_original_bytes(tmp_path: Path) -> None:
    source = tmp_path / "codes.txt"
    source.write_bytes(b"\nHP:0000001\nHP:0000002\nHP:0000001\n")

    codes, digest = load_hpo_codes(source)

    assert codes == ["HP:0000001", "HP:0000002"]
    assert digest == hashlib.sha256(source.read_bytes()).hexdigest()


@pytest.mark.parametrize("contents", ["HP:000001\n", "hp:0000001\n", "HP:00000001\n", " HP:0000001\n"])
def test_load_hpo_codes_rejects_invalid_line_without_echoing_it(tmp_path: Path, contents: str) -> None:
    source = tmp_path / "codes.txt"
    source.write_text(contents, encoding="utf-8")

    with pytest.raises(ValueError, match=r"^HPO code file contains an invalid line: 1$") as error:
        load_hpo_codes(source)

    assert contents.strip() not in str(error.value)


def test_load_hpo_codes_requires_existing_nonempty_input(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="^HPO code file does not exist$"):
        load_hpo_codes(tmp_path / "missing.txt")
    empty = tmp_path / "empty.txt"
    empty.write_text("\n \n", encoding="utf-8")
    with pytest.raises(ValueError, match="^HPO code file contains no valid codes$"):
        load_hpo_codes(empty)


def _snapshot(path: Path, *, wrong_contract: bool = False) -> None:
    path.mkdir()
    rationale = {column: "" for column in RATIONALE_COLUMNS}
    rationale.update({
        "RATIONALE_ID": "R1",
        "TARGET_SYMBOL": "GENE1",
        "EXPLORATORY_DAVID_OBSERVED_TERM_IDS": "hsa05170",
    })
    drug = {column: "" for column in DRUG_REVIEW_COLUMNS}
    drug.update({
        "RATIONALE_ID": "R1",
        "TARGET_SYMBOL": "GENE1",
        "DRUG_ID": "CHEMBL1",
        "DRUG_NAME": "Drug",
        "SOURCE_RECORD_IDS": "record-1",
        "ACTION_TYPE": "INHIBITOR",
        "PATHWAY_EVIDENCE_LEVEL": "WRONG" if wrong_contract else "RELAXED_DAVID_BENJAMINI_0_90",
    })
    pd.DataFrame([drug], columns=DRUG_REVIEW_COLUMNS).to_csv(path / "target_drug_review.csv", index=False)
    pd.DataFrame([rationale], columns=RATIONALE_COLUMNS).to_csv(path / "target_mechanism_rationale.csv", index=False)
    provenance = {"identity": {
        "selection_profile": "david_relaxed",
        "david_benjamini_threshold": 0.9,
        "drug_scope": "relaxed_pathway",
    }}
    (path / "compass_provenance.json").write_text(json.dumps(provenance), encoding="utf-8")


def test_load_compass_snapshot_validates_artifacts_and_preserves_source_row_order(tmp_path: Path) -> None:
    directory = tmp_path / "compass"
    _snapshot(directory)

    drugs, rationales, provenance, hashes = load_compass_snapshot(directory)

    assert list(drugs.columns) == list(DRUG_REVIEW_COLUMNS)
    assert list(rationales.columns) == list(RATIONALE_COLUMNS)
    assert drugs.iloc[0]["DRUG_ID"] == "CHEMBL1"
    assert provenance["identity"]["selection_profile"] == "david_relaxed"
    assert set(hashes) == {
        "target_drug_review.csv",
        "target_mechanism_rationale.csv",
        "compass_provenance.json",
    }
    assert all(len(value) == 64 for value in hashes.values())


def test_load_compass_snapshot_rejects_missing_artifact(tmp_path: Path) -> None:
    directory = tmp_path / "compass"
    directory.mkdir()

    with pytest.raises(ValueError, match="target_drug_review.csv"):
        load_compass_snapshot(directory)


def test_load_compass_snapshot_rejects_wrong_pathway_evidence_level(tmp_path: Path) -> None:
    directory = tmp_path / "compass"
    _snapshot(directory, wrong_contract=True)

    with pytest.raises(ValueError, match="PATHWAY_EVIDENCE_LEVEL"):
        load_compass_snapshot(directory)


def test_row_id_is_deterministic_and_duplicate_inputs_are_separate() -> None:
    from mva_hackathon.drug_hpo_evidence import assign_input_row_ids

    frame = pd.DataFrame([
        {"RATIONALE_ID": "R1", "TARGET_SYMBOL": "GENE1", "DRUG_ID": "D1", "ACTION_TYPE": "INHIBITOR", "SOURCE_RECORD_IDS": "S1"},
        {"RATIONALE_ID": "R1", "TARGET_SYMBOL": "GENE1", "DRUG_ID": "D1", "ACTION_TYPE": "INHIBITOR", "SOURCE_RECORD_IDS": "S1"},
    ])

    result = assign_input_row_ids(frame)

    assert result["INPUT_ROW_NUMBER"].tolist() == [0, 1]
    expected = hashlib.sha256(json.dumps(
        [0, "R1", "GENE1", "D1", "INHIBITOR", "S1"], ensure_ascii=False, separators=(",", ":")
    ).encode("utf-8")).hexdigest()
    assert result.loc[0, "ROW_ID"] == f"DHE-{expected[:16]}"
    assert result["ROW_ID"].str.startswith("DHE-").all()
    assert result["ROW_ID"].nunique() == 2


@pytest.mark.parametrize(
    ("identity_change", "error_field"),
    [
        ({"selection_profile": "validated"}, "selection_profile"),
        ({"david_benjamini_threshold": 0.05}, "david_benjamini_threshold"),
        ({"drug_scope": "decision_ready"}, "drug_scope"),
    ],
)
def test_load_compass_snapshot_rejects_wrong_provenance_identity(
    tmp_path: Path, identity_change: dict[str, object], error_field: str
) -> None:
    directory = tmp_path / "compass"
    _snapshot(directory)
    provenance_path = directory / "compass_provenance.json"
    provenance = json.loads(provenance_path.read_text(encoding="utf-8"))
    provenance["identity"].update(identity_change)
    provenance_path.write_text(json.dumps(provenance), encoding="utf-8")

    with pytest.raises(ValueError, match=error_field):
        load_compass_snapshot(directory)


def test_load_compass_snapshot_rejects_unresolved_rationale_and_missing_pathway_term(tmp_path: Path) -> None:
    directory = tmp_path / "compass"
    _snapshot(directory)
    drugs_path = directory / "target_drug_review.csv"
    drugs = pd.read_csv(drugs_path, dtype=str, keep_default_na=False)
    drugs.loc[0, "RATIONALE_ID"] = "UNKNOWN"
    drugs.to_csv(drugs_path, index=False)
    with pytest.raises(ValueError, match="RATIONALE_ID"):
        load_compass_snapshot(directory)

    drugs.loc[0, "RATIONALE_ID"] = "R1"
    drugs.to_csv(drugs_path, index=False)
    rationale_path = directory / "target_mechanism_rationale.csv"
    rationale = pd.read_csv(rationale_path, dtype=str, keep_default_na=False)
    rationale.loc[0, "EXPLORATORY_DAVID_OBSERVED_TERM_IDS"] = "hsa99999"
    rationale.to_csv(rationale_path, index=False)
    with pytest.raises(ValueError, match="EXPLORATORY_DAVID_OBSERVED_TERM_IDS"):
        load_compass_snapshot(directory)


def test_load_compass_snapshot_rejects_missing_required_csv_column(tmp_path: Path) -> None:
    directory = tmp_path / "compass"
    _snapshot(directory)
    drugs_path = directory / "target_drug_review.csv"
    drugs = pd.read_csv(drugs_path, dtype=str, keep_default_na=False)
    drugs.drop(columns=["DRUG_ID"]).to_csv(drugs_path, index=False)

    with pytest.raises(ValueError, match="DRUG_ID"):
        load_compass_snapshot(directory)

