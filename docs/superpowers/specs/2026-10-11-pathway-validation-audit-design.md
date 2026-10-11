# Pathway Validation Audit Design

## Purpose

Audit and re-analyze the independent MVA pathway-repurposing workflow (P0-P7) without changing its phase boundaries or existing STRING/drug-evidence architecture. Produce a private, local scientific validation report. Do not infer therapeutic efficacy from enrichment, source co-occurrence, or pharmacology.

## Constraints

- Preserve the existing P0-P7 module and command boundaries.
- Never send variants, genotypes, coordinates, phenotype narratives, internal scores, or patient identifiers externally.
- Do not run DAVID.
- Keep patient-derived tables, API responses, caches, and reports ignored and out of Git.
- Agent inspection remains data-blind: only artifact presence, aggregate counts, statuses, and sanitized failures may enter agent context.
- Use the established HTTP clients and existing drug clients; do not add a second transport layer.

## Audit scope

### P0-P2 statistics

Validate symbol-to-ENSG mapping status, foreground/background membership, GOSt custom-background payload, provider BH handling, effective-domain hypergeometric inputs, and fold enrichment. Re-run sensitivity analyses for passing threshold 2 and HP-term exclusion. Report that gene-length bias is unmeasured unless a documented gene-length source is supplied; do not fabricate a correction.

### P3-P5 evidence

Separate source-documented pathway/function/evidence records from therapeutic hypotheses. Preserve target identity, source-native mechanism and interaction identifiers, source failure status, safety uncertainty, and assay availability. Mark unavailable ChEMBL activity/assay evidence explicitly rather than treating absent mechanism records as no activity.

### P6 comparison

Compare pathway-only, STRING-only, shared, and combined target/drug sets. Report unique targets, unique drugs, available source coverage, directional contradiction flags, and descriptive overlap only. Do not calculate inferential overlap statistics without a prespecified eligible universe.

### Privacy

Emit a provenance inventory of each external request class and allowed identifier type. Include only list names and counts in the report; never copy identifiers from patient-derived input into agent-visible output.

### Reporting

Generate `pathway_validation_report.md` in the ignored pathway output directory. Include per-stage QC, aggregate counts, source coverage/failures, sensitivity status, limitations, corrected interpretation language, and prioritized corrections. The report may include local gene/pathway/drug evidence tables because it stays ignored; the agent will not inspect its rows.

## Implementation boundaries

Add a thin audit module and script. It reads stage outputs locally, computes aggregate QC and comparison tables, writes the ignored report atomically, and returns a sanitized status summary. It does not modify existing P0-P7 result rows. Correct only defects demonstrated by the audit: source failure being conflated with no-hit output, missing required provenance fields, ignored comparison inputs, and inaccurate completed statuses.

## Verification

- Synthetic fixtures prove payload privacy, effective-domain statistics, provider BH preservation, source failure handling, threshold/HP sensitivity, comparison partitioning, and report sections.
- Run the synthetic audit after the workflow fixture run.
- Run the full test suite and `uv lock --check`.
- Run the audit over authorized live local output. Confirm report existence and sanitized aggregate/stage status only.

## Non-goals

- DAVID benchmark or integration.
- Clinical interpretation, treatment recommendation, efficacy ranking, or causal claim.
- New pathway providers, a new database, or architecture rewrite.
