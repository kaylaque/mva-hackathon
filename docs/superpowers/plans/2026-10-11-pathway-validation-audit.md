# Pathway Validation Audit Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Produce defensible P0-P7 QC, sensitivity, privacy, and comparison reporting without changing branch architecture or making therapeutic claims.

**Architecture:** Extend the existing pathway module with pure table-level audit helpers and make orchestration report zero-result states truthfully. A new audit CLI reads only user-authorized pathway outputs, emits aggregate statistics and Markdown locally, and keeps agent-visible live output data-blind.

**Tech Stack:** Python 3.13, pandas, SciPy, existing `HttpClient` and drug clients, pytest.

**Spec:** `docs/superpowers/specs/2026-10-11-pathway-validation-audit-design.md`

## Global Constraints

- Preserve P0-P7 module and command boundaries.
- Never disclose variants, genotypes, coordinates, phenotype narratives, internal scores, or patient identifiers externally.
- Do not run DAVID.
- Agent data-blind exception: only `data/processed/track2/pathway/`, `data/interim/pathway_cache/`, and `.venv/` may be accessed; all other ignored paths remain immutable.
- Agent-visible live output may contain aggregate counts, statuses, artifact presence, and sanitized failures only.
- Reuse existing HTTP/drug clients; add no transport layer or dependencies.
- Do not claim enrichment, source co-occurrence, or pharmacology proves therapeutic efficacy.

## Review Focus

- No BH-significant terms: report `NO_SIGNIFICANT_ENRICHMENT`, header-only downstream tables, and no drug hypothesis.
- Client failure versus empty response: mark source unavailable, never no hits.
- Sensitivity analyses: threshold-2 and HP-excluded inputs must change actual term eligibility, never merely re-sort baseline rows.
- Comparison input omission: missing network context/drug tables must render unavailable coverage, not completed comparison.
- Privacy inventory: external request classes must list allowed identifier types and never patient fields.

---

### Task 1: Truthful zero-result stage states

**Files:**
- Modify: `src/mva_hackathon/pathway_workflow.py:344-597`
- Modify: `scripts/run_pathway_workflow.py:166-212`
- Test: `tests/test_pathway_workflow.py`

**Interfaces:**
- Produces `stage_summary(frame: pd.DataFrame, empty_reason: str) -> dict[str, object]`.
- Produces `empty_table(columns: list[str]) -> pd.DataFrame`.
- `run_p7` consumes stage status and reason.

- [ ] **Step 1: Write failing tests**

```python
def test_empty_p2_is_no_significant_enrichment():
    summary = stage_summary(pd.DataFrame(), "NO_SIGNIFICANT_ENRICHMENT")
    assert summary == {
        "status": "partial",
        "count": 0,
        "reason": "NO_SIGNIFICANT_ENRICHMENT",
    }

def test_p7_explains_zero_enrichment_chain(tmp_path):
    report = tmp_path / "report.md"
    workflow.run_p7(report, {"p2": {"status": "partial", "count": 0, "reason": "NO_SIGNIFICANT_ENRICHMENT"}})
    assert "NO_SIGNIFICANT_ENRICHMENT" in report.read_text()
```

- [ ] **Step 2: Run tests and verify failure**

Run: `uv run pytest tests/test_pathway_workflow.py -k 'empty_p2 or zero_enrichment' -v`

Expected: fail because `stage_summary` and report reason are absent.

- [ ] **Step 3: Implement minimal behavior**

```python
def stage_summary(frame: pd.DataFrame, empty_reason: str) -> dict[str, object]:
    if frame.empty:
        return {"status": "partial", "count": 0, "reason": empty_reason}
    return {"status": "completed", "count": len(frame)}
```

Use explicit P3/P5 schemas when serializing empty frames. Update P7 stage table and limitations text to render `reason` instead of generic `UNAVAILABLE`.

- [ ] **Step 4: Run focused tests**

Run: `uv run pytest tests/test_pathway_workflow.py -k 'empty_p2 or zero_enrichment' -v`

Expected: pass.

- [ ] **Step 5: Commit**

```bash
git add src/mva_hackathon/pathway_workflow.py scripts/run_pathway_workflow.py tests/test_pathway_workflow.py
git commit -m "fix(pathway): report empty enrichment states"
```

### Task 2: Statistical and sensitivity audit helpers

**Files:**
- Modify: `src/mva_hackathon/pathway_workflow.py`
- Test: `tests/test_pathway_workflow.py`

**Interfaces:**
- Produces `enrichment_qc(enrichment: pd.DataFrame) -> dict[str, object]`.
- Produces `sensitivity_summary(baseline: pd.DataFrame, threshold_two: pd.DataFrame, without_hp: pd.DataFrame) -> dict[str, object]`.

- [ ] **Step 1: Write failing tests**

```python
def test_enrichment_qc_counts_nominal_and_provider_bh_rows():
    qc = enrichment_qc(pd.DataFrame({"RAW_P": [0.01, 0.2], "BH_Q": [0.2, 0.3], "SOURCE": ["REAC", "HP"]}))
    assert qc["nominal_raw_p_lt_0_05"] == 1
    assert qc["provider_bh_le_0_05"] == 0

def test_sensitivity_summary_keeps_hp_exclusion_separate():
    summary = sensitivity_summary(baseline, threshold_two, without_hp)
    assert summary["hp_excluded"] is True
```

- [ ] **Step 2: Run tests and verify failure**

Run: `uv run pytest tests/test_pathway_workflow.py -k 'enrichment_qc or sensitivity_summary' -v`

Expected: fail because helpers do not exist.

- [ ] **Step 3: Implement minimal behavior**

Compute source counts, raw-p availability, nominal raw-p count, provider BH counts, mapping loss, and sensitivity availability. State `GENE_LENGTH_BIAS_UNMEASURED` unless a documented length table is supplied. Do not calculate a gene-length correction.

- [ ] **Step 4: Run focused tests**

Run: `uv run pytest tests/test_pathway_workflow.py -k 'enrichment_qc or sensitivity_summary' -v`

Expected: pass.

- [ ] **Step 5: Commit**

```bash
git add src/mva_hackathon/pathway_workflow.py tests/test_pathway_workflow.py
git commit -m "feat(pathway): add enrichment audit metrics"
```

### Task 3: Drug-source and comparison truthfulness

**Files:**
- Modify: `scripts/run_pathway_workflow.py:58-121, 195-208`
- Modify: `src/mva_hackathon/pathway_workflow.py:482-520`
- Test: `tests/test_pathway_workflow.py`

**Interfaces:**
- `_drug_records` returns `(records, source_statuses)` where every enabled source has `completed`, `empty`, or `unavailable`.
- `run_p6(pathway_targets, network_panel, network_context, network_drugs) -> dict[str, object]`.

- [ ] **Step 1: Write failing tests**

```python
def test_source_failure_is_not_empty_hit_set():
    records, status = _drug_records_with_client(failing_targets, failing_client)
    assert records == []
    assert status["ChEMBL"] == "unavailable"

def test_p6_marks_missing_drug_matrix_unavailable():
    summary = workflow.run_p6(pathway_targets, network_panel, None, None)
    assert summary["drug_comparison_status"] == "UNAVAILABLE"
```

- [ ] **Step 2: Run tests and verify failure**

Run: `uv run pytest tests/test_pathway_workflow.py -k 'source_failure or missing_drug_matrix' -v`

Expected: fail because source status and comparison inputs are ignored.

- [ ] **Step 3: Implement minimal behavior**

Preserve source-native IDs already returned by each client. Preserve DGIdb score/source/directionality and Open Targets clinical/safety metadata from existing response. Mark ChEMBL assay/activity evidence `UNAVAILABLE_NOT_RETRIEVED`; do not fabricate bioactivity. Read network context and drug paths only when supplied. Report unique target/drug counts, source coverage, directional hazards, and descriptive pathway/STRING partitions.

- [ ] **Step 4: Run focused tests**

Run: `uv run pytest tests/test_pathway_workflow.py -k 'source_failure or missing_drug_matrix' -v`

Expected: pass.

- [ ] **Step 5: Commit**

```bash
git add scripts/run_pathway_workflow.py src/mva_hackathon/pathway_workflow.py tests/test_pathway_workflow.py
git commit -m "fix(pathway): distinguish evidence gaps from no hits"
```

### Task 4: Privacy inventory and local audit report

**Files:**
- Create: `src/mva_hackathon/pathway_audit.py`
- Create: `scripts/audit_pathway_workflow.py`
- Modify: `src/mva_hackathon/pathway_workflow.py:522-587`
- Test: `tests/test_pathway_audit.py`

**Interfaces:**
- `audit_pathway_outputs(output_dir: Path) -> tuple[dict[str, object], str]`.
- `write_validation_report(output_dir: Path, summary: dict[str, object]) -> Path`.
- CLI exits `0` for completed/partial report; exits `1` for missing required stage artifacts.

- [ ] **Step 1: Write failing tests**

```python
def test_audit_report_lists_only_allowed_external_identifier_types(tmp_path):
    summary, report = audit_pathway_outputs(fixture_output_dir)
    assert summary["privacy"]["transmitted_identifier_types"] == ["HGNC symbol", "ENSG", "UniProt accession", "source target ID"]
    assert "VARIANT_KEY" not in report

def test_audit_report_explains_no_significant_enrichment(tmp_path):
    path = write_validation_report(tmp_path, {"p2": {"reason": "NO_SIGNIFICANT_ENRICHMENT"}})
    assert "No provider BH-adjusted term passed the configured threshold." in path.read_text()
```

- [ ] **Step 2: Run tests and verify failure**

Run: `uv run pytest tests/test_pathway_audit.py -v`

Expected: fail because audit module and CLI do not exist.

- [ ] **Step 3: Implement minimal behavior**

Read only P0-P6 tables and provenance from provided output directory. Emit aggregate row counts, stage statuses, mapping loss, source coverage, nominal/provider BH counts, P2/P3/P4/P5 zero-result chain, P6 partitions, sensitivity availability, evidence limitations, and prioritized corrections. Never emit patient rows. Write `pathway_validation_report.md` atomically. Include external identifier inventory: HGNC symbols, ENSG, UniProt accessions, ChEMBL/DGIdb/Open Targets source target IDs.

- [ ] **Step 4: Run focused tests**

Run: `uv run pytest tests/test_pathway_audit.py -v`

Expected: pass.

- [ ] **Step 5: Commit**

```bash
git add src/mva_hackathon/pathway_audit.py scripts/audit_pathway_workflow.py src/mva_hackathon/pathway_workflow.py tests/test_pathway_audit.py
git commit -m "feat(pathway): add validation audit report"
```

### Task 5: Verification and data-blind live audit

**Files:**
- Modify: `README.md`
- Modify: `scripts/README.md`
- Modify: `docs/JOURNEY.md`

- [ ] **Step 1: Run synthetic workflow and audit**

Run:

```bash
uv run python scripts/run_pathway_workflow.py --gene-mechanism tests/fixtures/pathway_sample/genes.csv --offline-fixture tests/fixtures/pathway_sample --output-dir /tmp/mva-pathway-sample
uv run python scripts/audit_pathway_workflow.py --output-dir /tmp/mva-pathway-sample
```

Expected: both exit `0`; temporary validation report exists.

- [ ] **Step 2: Run full verification**

Run:

```bash
uv run pytest
uv lock --check
```

Expected: all tests pass; lock validates.

- [ ] **Step 3: Run authorized live audit data-blind**

Run:

```bash
uv run python scripts/audit_pathway_workflow.py --output-dir data/processed/track2/pathway
```

Expected: aggregate-only process output; local `pathway_validation_report.md` exists. Do not print local report rows.

- [ ] **Step 4: Update documentation**

Document audit command, output path, source-failure semantics, no-significant-enrichment interpretation, unmeasured gene-length bias, and data-blind policy. Do not state top hypotheses unless provider BH criteria pass.

- [ ] **Step 5: Commit and publish**

```bash
git add README.md scripts/README.md docs/JOURNEY.md
git commit -m "docs(pathway): document validation audit"
git push
```

## Self-Review

Spec coverage: Tasks 1-2 cover P0-P2 statistics and sensitivity; Task 3 covers P4-P6 source/comparison defects; Task 4 covers privacy and cumulative audit reporting; Task 5 covers synthetic/live verification and documentation. P3 interpretation limits appear in Task 4 limitations because no extra scientific relationship is inferred. DAVID remains excluded. All interfaces are defined before use. Review-focus cases each have an owning test.
