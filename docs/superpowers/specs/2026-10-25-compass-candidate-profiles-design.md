# Design: COMPASS Candidate Evidence Profiles (investigation-status layer)

**Date:** 2026-10-25 · **Status:** Approved design, awaiting implementation plan
**Owner:** David (via omp) · **Scope:** extend COMPASS; no new pipeline, no DuckDB, no LLM

## 1. Goal and research question

COMPASS currently produces per-target decisions (`TEST_FIRST` / `ADVANCE_RESEARCH` / `HOLD` / `REJECT`) across two lanes (`target_mechanism_rationale.csv`, `target_drug_review.csv`) plus a rejection ledger. What it does **not** produce is one consolidated, per-candidate investigation status that answers: *for each drug–target candidate, is the accumulated evidence a supported hypothesis, a conditional hypothesis, insufficient evidence, or a rejection — and what information is still missing?*

This design adds exactly that, as a thin extension of COMPASS, and updates the related report surfaces. It consumes saved evidence (reuse decision confirmed by the user), retrieves only what is not yet covered (openFDA), and exports flat CSVs plus data-blind report/deck updates.

## 2. Non-goals (explicit)

- **No DuckDB / no relational store.** Flat CSV artifacts in the existing COMPASS output directory only. No new Python or system dependency.
- **No LLM layer.** The literature lane stays as cached, keyword-level evidence rows (`mechanism_evidence.csv`). LLM extraction is deferred until those rows prove insufficient; it is not in this slice.
- **No re-retrieval** of sources already covered by COMPASS or the pathway pipeline (Reactome/GO/KEGG/WikiPathways membership, UniProt/HGNC identity, PubMed literature, ChEMBL/Open Targets pharmacology).
- **No change to existing decisions.** `rationales.DECISION`, `drug_review.DECISION`, ledger rows, and test plans are inputs to the new stage, never modified by it.
- **No patient-level data enters COMPASS** (existing contract): no genotype, variant support, phenotype narratives, or patient identifiers. Mechanism direction is used only through saved gene-level columns (`DISEASE_STATE_DIRECTION`, `DESIRED_TARGET_CHANGE`, `ACTION_DIRECTION`, `DIRECTIONALITY_FIT`).

## 3. Changes

### 3.1 New source: `fetch_openfda_evidence()` in `src/mva_hackathon/compass_sources.py`

- Retrieves FDA drug labeling and adverse-event records for the drugs present in `drug_review` (matched by drug name/ID), via the existing cached-transport infrastructure (`perform_cached_request`, request envelopes, response SHA-256).
- Follows the established one-source pattern: `_envelope()` provenance, `SOURCE_STATUS` rows appended to the run's source status, failures recorded in `source_failures`, `offline` mode accepts only declared fixture responses and fails on a missing fixture (repo policy).
- Returns a frame with per-drug records: `DRUG_ID`, `DRUG_NAME`, `LABEL_STATUS` (e.g., labeled indication present / not found), `OPENFDA_FLAG_CODES` (serialized safety flags), `SOURCE_RECORD_IDS`, `SOURCE_STATUS`.
- No new HTTP client code beyond what `compass_sources` already owns.

### 3.2 New stage: `build_candidate_profiles()` in `src/mva_hackathon/compass.py`

Pure function, no I/O, fully unit-testable. Inputs: `rationales`, `drug_review`, `ledger`, `openfda` frame. Output: `candidate_evidence_profile.csv` with columns:

```
RATIONALE_ID, TARGET_SYMBOL, DRUG_ID, DRUG_NAME,
INVESTIGATION_STATUS, REASON_CODES, UNCERTAINTY_CODES,
MISSING_INFO, OPENFDA_FLAGS, EVIDENCE_IDS
```

Status derivation (precedence order, first match wins):

1. **`rejected`** — any `REJECT` in the rationale or drug review for the candidate, or an openFDA safety conflict recorded for the drug. A *safety conflict* means an openFDA flag from the declared conflict set (boxed warning, contraindication, or labeled hazard class) retrieved under a drug-ID match; name-only matches (`UNCERTAINTY_CODES=NAME_MATCH`) never trigger rejection on their own.
2. **`supported_hypothesis`** — rationale decision is `TEST_FIRST` or `ADVANCE_RESEARCH`, drug review `DIRECTIONALITY_FIT` is compatible (no direction conflict), `HAZARD_STATUS` clean, and the merged evidence spans ≥2 independent `EVIDENCE_CLASS`es including at least one non-computational class (curated/clinical/experimental).
3. **`conditional`** — rationale decision is `TEST_FIRST`/`ADVANCE_RESEARCH` but the candidate is single-source (one evidence class) or has coverage gaps (empty openFDA row, missing clinical corroboration). Every gap is enumerated in `MISSING_INFO`.
4. **`insufficient_evidence`** — rationale decision is `HOLD`, or no drug evidence exists for the target and retrieval found none.

`REASON_CODES` reuse the existing reason-code vocabulary wherever possible and add a small new set (`PROFILE_OPENFDA_CONFLICT`, `PROFILE_SINGLE_SOURCE`, `PROFILE_MISSING_OPENFDA`, `PROFILE_MISSING_CLINICAL_CORROBORATION`, `PROFILE_NO_DRUG_EVIDENCE`). `UNCERTAINTY_CODES` records interpretive caveats (e.g., exploratory lane, David-relaxed profile) without changing the status. `MISSING_INFO` is a semicolon-joined enumeration per candidate — the "what information is missing" requirement — never a null or zero.

### 3.3 Orchestration changes in `run_compass()`

After `drug_review` and before the staging write: call `fetch_openfda_evidence()`, then `build_candidate_profiles()`. The profile frame joins the staging frame set and `_COMPASS_TABLES`, gets a SHA-256 in `output_sha256`, and counts of each status are added to the returned summary dict. CLI gains no new flags in this slice (openFDA retrieval is part of the standard run; offline mode already covers fixture-driven execution).

### 3.4 Report updates (data-blind, committed)

- `compass_report.md`: append a section with aggregate status counts (a table of the four statuses with row counts and the top reason codes) and the non-diagnostic scope note. Aggregate counts and status vocabulary only — no patient- or variant-level fields.
- `docs/decks/index.html`: one new page inserted after the existing drug-screening decision page (page 18), in the established deck style (static HTML, no scripts), describing the evidence-integration stage: input tables, status vocabulary, what each status means, and its limits. Page numbering updated (19 pages → 20).
- `docs/JOURNEY.md`: add a short subsection under the COMPASS record documenting the new artifact and status semantics once the run output exists.

Committed content must satisfy the data-blind report policy: statuses, aggregate counts, reason codes, gene/drug decision-table fields are allowed; coordinates, alleles, genotypes, patient identifiers are prohibited.

## 4. Data flow

```text
saved post-analysis tables + network projections
  → (unchanged COMPASS stages: eligibility → memberships → protein roles
      → literature → rationales → drug review → ledger → test plans)
  → fetch_openfda_evidence (new, cached)
  → build_candidate_profiles (new, pure)
  → candidate_evidence_profile.csv + updated compass_report.md + provenance
  → (operator) data-blind deck page + JOURNEY subsection update
```

Outputs land in the operator-specified `--output-dir` (external to the repo, like all COMPASS/DAVID artifacts). Committed surfaces receive only aggregate, data-blind content.

## 5. Testing

- `tests/test_compass.py`: profile-derivation tests covering all four statuses, precedence (safety conflict beats positive decisions), reason/uncertainty code emission, `MISSING_INFO` enumeration on single-source and empty-coverage candidates, and the invariant that existing decision columns are untouched.
- `tests/test_compass_sources.py`: openFDA fixtures — cached hit, offline miss fails loudly, malformed response envelope rejected, source-status rows appended.
- Existing suites stay green: offline runs must produce `source_incomplete` status (exit 2) when openFDA fixtures are absent, not fabricate rows.
- Verification: `uv run pytest tests/test_compass.py tests/test_compass_sources.py`, plus one offline end-to-end run via `scripts/run_compass.py --offline` against fixture inputs to observe the new artifact and report section.

## 6. Risks and mitigations

- **openFDA matching by drug name is lossy.** Mitigation: match on the drug IDs COMPASS already records; when only a name match is possible, mark the row `UNCERTAINTY_CODES=NAME_MATCH` and never treat it as a safety conflict on its own.
- **Status vocabulary drift** (four new statuses vs. existing decision codes). Mitigation: statuses are derived exclusively by the documented precedence from existing decisions; the mapping is one pure function with table-driven tests.
- **Deck renumbering churn.** Mitigation: single insertion point; page-number strings updated mechanically in the same edit.
