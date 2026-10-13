# Integrated DAVID Pathway Design

## Goal

Run the existing P0–P7 pathway workflow and an opt-in DAVID evidence addendum from the same fresh P0/P1 outputs, while retaining DAVID's native-web boundary, external artifacts, resumable provenance, and provider-specific interpretation.

## Scope

`--with-david` extends a normal full pathway invocation. It does not replace g:Profiler, rerun P0–P7 during DAVID import, submit to DAVID programmatically, expand targets/drugs, or infer diagnosis/treatment.

## Invocation contract

```text
python scripts/run_pathway_workflow.py \
  --gene-mechanism <gene-level-input> \
  --output-dir <pathway-output> \
  --with-david \
  --david-output-dir <external-david-output> \
  [--david-exports <external-native-export-dir>]
```

- `--with-david` is valid only for the complete P0–P7 workflow (`--stage all`). It requires `--david-output-dir`.
- Without `--david-exports`, the runner completes P0–P7, prepares DAVID handoff artifacts from the resulting P0/P1 tables, records `awaiting_web_results`, and exits `2`.
- With verified native exports, it completes P0–P7 then imports DAVID evidence and returns the ordinary full-workflow status unless DAVID validation fails.
- `--stage david --reuse-pathway-dir <pathway-output> --output-dir <external-david-output> [--david-exports <exports>]` remains the import-only resume route. It must not instantiate network/provider clients or rerun P0–P7.
- `--david` remains an alias for the import-only route for compatibility with the existing CLI tests. It may not be combined with `--with-david`.

## Data flow and ownership

```text
Gene-level input
  └─ P0–P7 ──> pathway output (existing artifacts)
                   ├─ P0/P1 projection ──> DAVID lists + manifest (external)
                   │                           └─ native UI chart/clusters (external)
                   └─ pathway provenance        └─ DAVID import/comparison/report (external)
```

DAVID reads only the existing allowlisted P0/P1 columns. The input lists, manifest, receipt, native downloads, normalized DAVID tables, comparison, report, and telemetry reside in the explicitly supplied external directory, never under the repository. Existing P0–P7 artifacts remain unmodified by DAVID import.

The normal pathway provenance gains a small `david` summary only when `--with-david` is used:

```json
{
  "status": "awaiting_web_results | verified_complete | blocked | failed",
  "run_key": "sha256",
  "external_output_dir": "/absolute/external/path",
  "reason": "optional machine-readable status reason"
}
```

It contains no gene lists, patient/variant fields, browser cookies, credentials, DAVID rows, or export contents. The external manifest remains authoritative for DAVID identity and result hashes.

## Native DAVID execution boundary

The repository prepares/imports files; a managed browser performs the one native UI session. The UI procedure is:

1. Open a clean DAVID workspace and inspect current accessible controls.
2. Upload the immutable foreground list as human Ensembl identifiers; record submitted, mapped, and unique DAVID-gene counts separately.
3. Upload the exact background list and confirm the selected current population is the uploaded custom population.
4. Run one functional annotation chart and one functional annotation clustering evaluation; export their native CSVs.
5. Write a receipt with list/run identity, population confirmation, observed category/settings, report-completion state, and export hashes.
6. Invoke the import-only route with that export directory.

A click or page-observation timeout has an **unknown** submission state. The runner must record `blocked` with `DAVID_SUBMISSION_UNKNOWN`; it must not retry or resubmit automatically. A user may explicitly start a fresh clean session only after confirming whether the prior DAVID list exists. A known pre-submission transient page-load failure permits one normal page reload. Authentication, CAPTCHA, quota, or policy restrictions are terminal blocked states, not conditions to bypass.

DAVID can still be unavailable or slow; the system guarantees safe preservation and resumption, not remote-site uptime.

## Resumption and provenance

`run_david_addition(baseline_dir, output_dir, exports_dir, *, force=False) -> int` remains the DAVID boundary:

- `0`: imports verified exports or returns a hash-verified cached completion.
- `2`: valid lists/manifest prepared and native exports are awaited.
- `1`: invalid path, malformed data/export, stale identity, or failed validation.

The foreground/background list bytes, projection hashes, namespace, taxon, and settings derive a run key. A completed manifest requires its list/result hashes to verify before reuse. A changed identity, missing manifest, changed export hash, unconfirmed population, or partial/unknown receipt fails closed. `--force` may regenerate only external DAVID artifacts; it never reruns P0–P7 or treats old exports as valid for a new identity.

## Evidence semantics

DAVID is secondary functional interpretation alongside g:Profiler. P-Value, Benjamini, FDR, and cluster enrichment score remain separate native fields. Stable identifiers only match provider terms; title matching and pooled p-values are prohibited. Shared annotations are descriptive overlap, not independent replication. DAVID output never causes P2–P7 target/drug expansion.

## Failure behavior

- P0/P1 missing or invalid: no DAVID UI action; DAVID status is `failed`.
- List >3,000 IDs: `blocked` before submission, without truncating/chunking.
- UI uncertainty: `blocked` and preserve external prepared artifacts.
- Native exports absent: `awaiting_web_results` after P0–P7 completed.
- Export/receipt validation failure: `failed`, preserve external source exports, and do not modify baseline P0/P1.
- P0–P7 partial status and DAVID completion coexist in provenance; either condition makes the combined command return `2`.

## Test and verification requirements

Synthetic tests must prove:

- `--with-david` runs P0–P7 exactly once, prepares from fresh P0/P1, and does not construct legacy DAVID-only paths.
- An awaiting DAVID handoff returns `2` while retaining all normal P0–P7 outputs.
- Export import succeeds without rerunning P0–P7 and is byte-identical when replayed.
- DAVID failure/blocked status does not overwrite baseline artifacts or trigger resubmission.
- Native P-Value, Benjamini, FDR, and cluster score remain distinct; stable-ID comparison never title-matches.
- No patient/variant/private columns cross the DAVID boundary.

Verification runs the DAVID/pathway focused suite, the full suite, one synthetic full CLI preparation/import/cache replay, and one native full run only after local checks pass. The native run is reported blocked if browser/site controls fail; it is never replaced by a synthetic success claim.
