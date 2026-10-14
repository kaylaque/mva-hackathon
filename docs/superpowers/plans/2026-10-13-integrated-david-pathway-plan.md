# Integrated DAVID Pathway Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add opt-in DAVID preparation/import to a completed P0–P7 pathway run without rerunning baseline analysis or unsafe DAVID resubmission.

**Architecture:** Keep `david_workflow.py` as the external-artifact boundary. Extend it with a compact, data-blind manifest summary for the main runner. The normal runner owns the new `--with-david` orchestration: P0–P6 create a fresh reusable baseline, DAVID prepares/imports externally, then P7 and provenance record the DAVID state.

**Tech Stack:** Python 3.13, pandas, unittest/pytest, existing `PathwayWorkflow` atomic writers, managed browser for the one native DAVID session.

**Spec:** `docs/superpowers/specs/2026-10-13-integrated-david-pathway-design.md`

## Global Constraints

- `--with-david` is opt-in and valid only with `--stage all`; ordinary P0–P7 invocation remains behaviorally unchanged.
- DAVID outputs require an explicit directory outside the repository; no generated DAVID artifact belongs in `data/`, `results/`, or another ignored repository path.
- Only P0 `SYMBOL`, `ENSEMBL_ID`, `STATUS`, `IN_FOREGROUND` and the specified P1 enrichment columns may cross the DAVID module boundary.
- DAVID never performs HTTP requests, browser automation, credential handling, API use, or automatic resubmission.
- `0` means verified DAVID completion, `2` means awaiting/partial/blocked, and `1` means invalid/stale/failed input; a combined command returns `2` when either P0–P7 or DAVID is partial/awaiting/blocked.
- DAVID P-Value, Benjamini, FDR, and cluster score remain provider-native and are never pooled with g:Profiler values.
- Use synthetic tables and external temporary directories only; never read patient-level data except allowed programmatic P0/P1 projections.

## Review Focus

- Existing full runs without `--with-david` must retain exactly their prior outputs, stage statuses, and exit behavior.
- A DAVID handoff after P0–P6 must not let P7/provenance overwrite its external manifest or rerun P0/P1 during import-only resume.
- A `--with-david` output path under the repository, equal to the pathway baseline, or ignored outside `data/processed` must fail before external output creation.
- A timed-out/unknown browser submission must be representable as `blocked` without changing list identity or triggering an automatic submission retry.
- A successful DAVID import plus an independently partial P0–P7 stage must return `2`, preserving both statuses in provenance.

---

### Task 1: Add data-blind DAVID run summaries

**Files:**
- Modify: `src/mva_hackathon/david_workflow.py:41-50, 808-865`
- Modify: `tests/test_david_workflow.py:19-119`

**Interfaces:**
- Consumes: existing `david_manifest.json` written by `run_david_addition`.
- Produces: `david_run_summary(output_dir: Path, exit_code: int, *, reason: str | None = None) -> dict[str, Any]` for the runner.

- [ ] **Step 1: Write failing summary tests**

```python
def test_summary_exposes_only_status_identity_and_external_location(self) -> None:
    status = run_david_addition(baseline, output, None)
    summary = david_run_summary(output, status)
    self.assertEqual(summary["status"], "awaiting_web_results")
    self.assertEqual(summary["external_output_dir"], str(output.resolve()))
    self.assertIn("run_key", summary)
    self.assertNotIn("foreground", json.dumps(summary))
    self.assertNotIn("PRIVATE_VARIANT", json.dumps(summary))


def test_summary_marks_unknown_submission_as_blocked(self) -> None:
    summary = david_run_summary(output, 2, reason="DAVID_SUBMISSION_UNKNOWN")
    self.assertEqual(summary["status"], "blocked")
    self.assertEqual(summary["reason"], "DAVID_SUBMISSION_UNKNOWN")
```

- [ ] **Step 2: Run the new tests and confirm RED**

Run: `PYTHONDONTWRITEBYTECODE=1 /tmp/mva-david-venv/bin/python -m pytest -p no:cacheprovider tests/test_david_workflow.py -k summary`

Expected: import failure because `david_run_summary` does not exist.

- [ ] **Step 3: Implement the narrow summary helper**

```python
def david_run_summary(output_dir: Path, exit_code: int, *, reason: str | None = None) -> dict[str, Any]:
    manifest = _load_manifest(Path(output_dir) / _MANIFEST_NAME)
    status = "blocked" if reason else str(manifest["status"])
    return {
        "status": status,
        "run_key": str(manifest["run_key"]),
        "external_output_dir": str(Path(output_dir).resolve()),
        **({"reason": reason} if reason else {}),
    }
```

Validate that the manifest has string `run_key` and a known status (`awaiting_web_results` or `verified_complete`) before returning. Raise `DavidWorkflowError` for invalid summary state. Do not expose list hashes, counts, source rows, or receipt contents.

- [ ] **Step 4: Run DAVID workflow tests and confirm GREEN**

Run: `PYTHONDONTWRITEBYTECODE=1 /tmp/mva-david-venv/bin/python -m pytest -p no:cacheprovider tests/test_david_workflow.py`

Expected: all DAVID workflow tests pass.

- [ ] **Step 5: Commit the bounded module change**

```bash
git add src/mva_hackathon/david_workflow.py tests/test_david_workflow.py
git commit -m "feat(pathway): expose DAVID run summary"
```

### Task 2: Orchestrate opt-in DAVID beside P0–P7

**Files:**
- Modify: `scripts/run_pathway_workflow.py:31-51, 127-238`
- Modify: `src/mva_hackathon/pathway_workflow.py:569-660`
- Modify: `tests/test_david_cli.py:22-46`

**Interfaces:**
- Consumes: `run_david_addition(baseline_dir, output_dir, exports_dir, *, force=False) -> int` and `david_run_summary(output_dir, exit_code, *, reason=None) -> dict[str, Any]` from Task 1.
- Produces: a full-run DAVID stage in `workflow.provenance(stages)` and P7 report text; normal mode remains unchanged when `args.with_david` is false.

- [ ] **Step 1: Write failing integrated CLI tests**

```python
def test_with_david_runs_baseline_then_prepares_external_handoff(self) -> None:
    with patch.object(runner, "run_david_addition", return_value=2) as david:
        result = runner.main()  # argv includes --with-david and --david-output-dir
    david.assert_called_once_with(pathway_output, david_output, None, force=False)
    provenance = json.loads((pathway_output / "pathway_provenance.json").read_text())
    self.assertEqual(provenance["stages"]["david"]["status"], "awaiting_web_results")
    self.assertEqual(result, 2)


def test_with_david_rejects_non_all_stage_before_workflow_construction(self) -> None:
    with patch.object(runner, "PathwayWorkflow", side_effect=AssertionError("must not construct")):
        with self.assertRaisesRegex(ValueError, "stage all"):
            runner.main()  # argv includes --with-david --stage p1
```

Use the existing offline fixture path for normal P0–P7 data and an external `TemporaryDirectory` for DAVID output. Add cases that reject missing `--david-output-dir` and `--with-david` combined with `--david`.

- [ ] **Step 2: Run the new CLI tests and confirm RED**

Run: `PYTHONDONTWRITEBYTECODE=1 /tmp/mva-david-venv/bin/python -m pytest -p no:cacheprovider tests/test_david_cli.py -k with_david`

Expected: argument parsing rejects unknown `--with-david`.

- [ ] **Step 3: Add explicit CLI validation and full-run DAVID orchestration**

Add arguments:

```python
parser.add_argument("--with-david", action="store_true", help="Prepare or import opt-in native DAVID evidence after P0–P7.")
parser.add_argument("--david-output-dir", type=Path)
```

Before logger/client construction, reject: `--with-david` with `--david`, `--stage != all`, or missing `--david-output-dir`. Keep the existing DAVID-only early return. After P0–P6 writes fresh P0/P1, call:

```python
david_exit = run_david_addition(output_dir, args.david_output_dir, args.david_exports, force=args.force)
stages["david"] = david_run_summary(args.david_output_dir, david_exit)
```

Call this before P7 so P7 can render the stage. When DAVID returns `1`, set a safe failed summary rather than suppressing the failure. Extend `PathwayWorkflow.run_p7` with a compact DAVID status line that contains only status, run key, external directory, and optional reason. Preserve every existing P7 section.

Return `2` if a DAVID stage has `awaiting_web_results`, `blocked`, `partial`, or `failed`; otherwise preserve existing stage-result behavior.

- [ ] **Step 4: Run integrated CLI tests and confirm GREEN**

Run: `PYTHONDONTWRITEBYTECODE=1 /tmp/mva-david-venv/bin/python -m pytest -p no:cacheprovider tests/test_david_cli.py tests/test_pathway_workflow.py`

Expected: all existing CLI/pathway tests and new integrated tests pass.

- [ ] **Step 5: Commit the orchestration change**

```bash
git add scripts/run_pathway_workflow.py src/mva_hackathon/pathway_workflow.py tests/test_david_cli.py
git commit -m "feat(pathway): run DAVID alongside full workflow"
```

### Task 3: Prove resumption, status precedence, and artifact preservation

**Files:**
- Modify: `tests/test_david_workflow.py:75-116`
- Modify: `tests/test_david_cli.py:22-46`
- Modify: `.gitignore:214-226`

**Interfaces:**
- Consumes: Tasks 1–2 full-run `david` provenance stage and existing DAVID manifest/receipt validation.
- Produces: regression proof that import-only resume never reruns the baseline and that no DAVID artifact is versioned under the repository.

- [ ] **Step 1: Write failing transition tests**

```python
def test_full_run_import_replay_does_not_rerun_baseline(self) -> None:
    first = runner.main()  # offline fixture + --with-david + matching exports
    before = {path.name: path.read_bytes() for path in pathway_output.iterdir()}
    resumed = runner.main()  # --stage david --reuse-pathway-dir pathway_output
    self.assertEqual(resumed, 0)
    self.assertEqual(before, {path.name: path.read_bytes() for path in pathway_output.iterdir()})


def test_partial_baseline_and_verified_david_return_partial(self) -> None:
    result = runner.main()  # fixture induces existing P0–P7 partial stage plus valid DAVID exports
    self.assertEqual(result, 2)
    self.assertEqual(provenance["stages"]["david"]["status"], "verified_complete")
```

Add a test for `DAVID_SUBMISSION_UNKNOWN` summary: no list/export replacement and no second call to any DAVID UI/network abstraction. Assert `/data/processed/track2/pathway_david/` remains the sole future-output ignore rule; do not create it.

- [ ] **Step 2: Run transition tests and confirm RED**

Run: `PYTHONDONTWRITEBYTECODE=1 /tmp/mva-david-venv/bin/python -m pytest -p no:cacheprovider tests/test_david_workflow.py tests/test_david_cli.py -k "replay or partial_baseline or unknown"`

Expected: at least one assertion fails until full-run status precedence and resume behavior are wired.

- [ ] **Step 3: Implement only corrections exposed by the tests**

Keep import-only `--stage david` on the existing early return. Ensure it neither creates the pathway output directory nor writes `pathway_provenance.json`. Ensure the full runner treats `david_exit == 1` as a failed DAVID stage and a final exit `2`, rather than converting it to a successful path. Do not add DAVID browser code, retries, credentials, cookies, or network clients.

- [ ] **Step 4: Run focused regression suite and confirm GREEN**

Run: `PYTHONDONTWRITEBYTECODE=1 /tmp/mva-david-venv/bin/python -m pytest -p no:cacheprovider tests/test_david_workflow.py tests/test_david_cli.py tests/test_pathway_workflow.py`

Expected: all DAVID and pathway tests pass.

- [ ] **Step 5: Commit preservation proof**

```bash
git add tests/test_david_workflow.py tests/test_david_cli.py .gitignore
git commit -m "test(pathway): cover DAVID full-run resumption"
```

### Task 4: Verify end-to-end and perform one controlled native run

**Files:**
- Modify: `README.md` only if the project currently has a runnable-workflow section; otherwise no documentation change.
- External only: a unique DAVID run directory and native DAVID export directory.

**Interfaces:**
- Consumes: combined CLI from Task 2 and verified export importer from Task 3.
- Produces: command evidence for the synthetic full flow and a real-run outcome labeled verified complete or blocked.

- [ ] **Step 1: Run the full automated suite**

Run: `PYTHONDONTWRITEBYTECODE=1 /tmp/mva-david-venv/bin/python -m pytest -p no:cacheprovider`

Expected: every collected test passes.

- [ ] **Step 2: Run the synthetic complete CLI regression**

Run: `PYTHONDONTWRITEBYTECODE=1 /tmp/mva-david-venv/bin/python -m pytest -p no:cacheprovider tests/test_david_cli.py::DavidCliTests::test_full_cli_prepare_import_and_resume`

The test must create its synthetic gene input, offline fixture, pathway output, DAVID output, and native exports under `tempfile.TemporaryDirectory()`. It invokes `scripts/run_pathway_workflow.py` via `subprocess.run` three times: full preparation must exit `2`; full import with the receipt/chart/clusters must exit `0` unless the fixture deliberately has a partial P0–P7 status; DAVID-only resume must exit `0` and preserve pathway-output bytes.

- [ ] **Step 3: Run one native DAVID session after local success**

Use a clean managed browser only. Prepare real lists with the integrated full command through the data-blind targeted notebook procedure. Upload foreground once, upload and confirm custom background, run one chart and clustering report, save native exports externally, write the receipt, and import through the DAVID-only resume command. On UI timeout, auth/CAPTCHA/quota/policy restriction, or missing baseline: preserve the external handoff, report `blocked`, and do not resubmit automatically.

- [ ] **Step 4: Commit verified source/test/documentation changes**

```bash
git add src/mva_hackathon/david_workflow.py src/mva_hackathon/pathway_workflow.py \
  scripts/run_pathway_workflow.py tests/test_david_workflow.py tests/test_david_cli.py \
  .gitignore README.md
git diff --cached --name-only
git commit -m "feat(pathway): integrate resumable DAVID evidence"
```

Stage only paths actually changed; omit `README.md` if unchanged. Never stage generated artifacts, `data/`, `results/`, exports, receipts, manifests, lists, notebooks, or secrets.

## Outcome (executed 2026-10-13)

Steps 1–4 completed. Full suite: 91 tests plus 13 subtests passing.

Native run executed through the managed browser against the real baseline, then
imported through the DAVID-only CLI: import exited `0`, and an identical replay
exited `0` with byte-identical imported artifacts and unchanged baseline size,
mtime, and projection hashes. Measured counts are recorded in `docs/JOURNEY.md`.

Real data exposed three defects in the first implementation, each now fixed with
a regression test:

1. Native exports carry the term identifier in `Full Term`, not `Term`, so the
   provider comparison matched no terms at all (0 shared of 318 chart rows).
2. g:Profiler identifiers are source-prefixed (`REAC:R-HSA-…`, `WP:WP5344`), and
   the KEGG branch compared `"HSA"` against an already-lowercased identifier, so
   those terms never resolved to a stable identity.
3. Resume after a completed import reported stale provenance because the
   top-level manifest kept outdated `result_sha256` values, and `--force` did
   not regenerate prepared lists that had been edited.

Two deliberate limits:

- `--with-david` was not executed against the real baseline, because that would
  re-run the P0–P7 network stages. Orchestration order is covered by the
  synthetic end-to-end run and the CLI tests only.
- The user's main checkout still holds divergent uncommitted copies of the same
  files, so publication used the clean linked worktree. Reconciliation of that
  checkout is left to the user.

Process note for future UI automation: DAVID's upload input is hidden and its
taxonomy/background widgets only initialize when the modal is opened through the
site's own button handler; synthetic change events are ignored. Drive the native
file chooser and dispatch real pointer events on the drop zone instead.
