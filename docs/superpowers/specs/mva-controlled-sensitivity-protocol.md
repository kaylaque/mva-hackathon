# Controlled MVA Pathway Sensitivity Protocol

## Purpose

This protocol defines a provenance-bound, sensitivity-only analysis of an existing single-case MVA pathway baseline. It preserves the baseline as a negative-result artifact and does not create target or drug hypotheses. It is research exploration, not diagnosis, treatment advice, or evidence of clinical efficacy.

## Scientific contract

- Baseline enrichment significance governs claims about baseline enrichment. A sensitivity result cannot replace a baseline null.
- Independently supported direct-gene hypotheses do not require pathway enrichment, but must remain on a separately documented mechanism-first route.
- Sensitivity results are exploratory because the multi-hit selection idea and apparent prior result have already been observed; they are not preregistered independent confirmation.
- Returned provider terms are not automatically the complete multiple-testing family.
- Dilution, background composition, and source changes are possible explanations for differences; this workflow does not prove a causal decomposition.
- No analysis may select disease- or pathway-labelled genes as independent enrichment evidence for that same disease or pathway.

## Privacy and input boundary

The tool may project baseline P0 fields `SYMBOL`, `ENSEMBL_ID`, `STATUS`, and `IN_FOREGROUND`, and P1 term/statistical fields required for interpretation. It must not read, hash, write, display, or transmit per-gene variant counts, coordinates, alleles, genotypes, sample identifiers, phenotype narratives, patient identifiers, or profile fields such as `PASSING_VARIANTS` and `TOTAL_VARIANTS`.

The multi-hit input is an operator-prepared CSV with exactly one column, `SYMBOL`. The operator certifies that each listed gene was selected by at least two distinct passing variants. This membership does not establish compound heterozygosity or trans phasing. Symbols are trimmed and uppercased; blanks, unknown symbols, and symbols outside the baseline foreground are invalid. Duplicate symbols collapse. Unmapped eligible symbols remain documented exclusions and are never remapped.

## Immutable artifacts and paths

The baseline is read-only. The output root must be a newly selected, non-ignored directory outside the repository and must not exist at startup. It must be disjoint from the baseline after path resolution; symlink escapes and ignored repository destinations are invalid. There is no overwrite, force, or resume mode.

The sensitivity manifest is a new artifact family. It binds canonical allowlisted P0/P1 projections, input membership, canonical request JSON, response-derived files, and output hashes. It must not imitate or weaken the legacy audit's full-P0 provenance contract.

## Fixed analyses

The tool constructs exactly two requests, each with `organism="hsapiens"`, `domain_scope="custom"`, `significance_threshold_method="fdr"`, `all_results=true`, `no_evidences=false`, and `user_threshold=0.05`.

1. **multihit**: the approved, mapped multi-hit subset as query; the baseline mapped background unchanged; sources `GO:BP`, `REAC`, `WP`, and `HP`.
2. **hp_excluded**: the original mapped baseline foreground and background unchanged; sources `GO:BP`, `REAC`, and `WP`.

The tool freezes and hashes canonical request JSON before any provider call. It records a canonical projected baseline identity, input membership hash, request hash, retrieval date, sources, and any provider source-version metadata. It does not add threshold sweeps, source reselection, a baseline re-query, arbitrary subsets, drug expansion, STRING, DAVID, or ranked/bias-adjusted analysis.

## Execution and interpretation

Without `--execute` or `--offline-fixture`, the status is `PREPARED_NOT_EXECUTED`; this is not a scientific result. `--execute` submits only the two prepared g:Profiler requests using the existing HTTP client. `--offline-fixture` loads exactly `multihit.json` and `hp_excluded.json` and is mutually exclusive with live execution. A failure in one analysis remains unavailable and does not discard the other analysis.

Provider metadata is checked for the requested organism, FDR correction, and custom domain when available. Contradictions fail that analysis. Absent confirmation yields `UNVERIFIED_PROVIDER_METADATA`; it cannot support a claim of verified significance. Missing provider-version data is recorded as unknown. Source-version mismatch or unavailable identity makes baseline comparisons confounded or unverified rather than a foreground-only causal effect.

Responses must have valid finite probability values in [0, 1], non-negative counts, intersections no larger than query or term size, and query/term sizes no larger than their domain. Invalid values are provider-data failures, not negative outcomes. Provider-adjusted values remain separate from locally reconstructed hypergeometric raw p-values; adjusted values are never corrected again. Zero returned terms and zero significant terms are completed negative outcomes with header-only files as needed.

Per analysis, the output contains `request.json`, `manifest.json`, `enrichment.csv`, and `prioritized.csv`. The root report records every attempted analysis and distinguishes nominal p-values, adjusted values, effect sizes, metadata uncertainty, unavailable results, and negative results. It explicitly states that no automated target or drug generation occurred. Console output is limited to fixed run names, statuses, and permitted pathway/FDR summaries.

## Evidence-review protocol

No automated candidate creation is part of this implementation. Human reviewers use one of these mutually explicit routes:

### Mechanism-first direct route

Require credible genetic and inheritance evidence, a variant-specific defect hypothesis, target-specific drug action, pharmacology and safety review, and a falsifiable rescue experiment. Pathway ID and FDR are nullable and explicitly not applicable when no enrichment supports this route.

### Pathway-derived route

Require a significant term with named gene overlap, verified pathway membership, independent functional evidence, an explicit compensatory rationale, and target-specific pharmacology. A sensitivity-only hit remains exploratory until independent work addresses its selection bias and reproducibility. Enrichment alone never promotes a drug.

### Network context

STRING proximity is descriptive functional association. It is not activation, causal rescue, target eligibility, or independent replication.

### Exploratory leads

Biologically motivated non-significant ideas may remain separately labelled exploratory. Quarantine action-unknown compounds. Reject or quarantine mechanistic hazards with documented reasons; retain an audit record instead of silently deleting them.

## Manual candidate record

A manual record contains: route; genetic-anchor reference without patient coordinates; defect hypothesis; nullable pathway ID/FDR; target identity; normalized compound ID; target-specific action; expected rescue direction; source-native evidence IDs or URLs; evidence grade; potency, selectivity, exposure, and safety availability; uncertainty; rescue assay; and verdict. Each claim is labelled observed, curated, inferred, or speculative.

Syndicated ChEMBL evidence through Open Targets or DGIdb is one underlying evidence stream, not multiple independent experiments. Missing assay or exposure evidence means unknown, not safe or ineffective.

Rescue validation requires mechanism-matched endpoints such as segregation errors, micronuclei, or checkpoint timing; appropriate controls; dose response; and separate viability/toxicity measurements. Killing abnormal cells is not automatically normal-tissue rescue. Clinical treatment or efficacy claims require indication-specific clinical evidence, independent specialist review, and appropriate regulatory and ethical oversight; this workflow makes none.

## Deferred analyses and external comparison

Ranked analysis is deferred until there is a defensible continuous score over the eligible background, explicit tie handling, evidence-independent selection, and an appropriate competitive null. Bias-aware analysis is deferred until callable coverage, gene length/mutability, and annotation covariates are documented, followed by matching/permutation or a validated selection-bias model. The current gene count is not a cohort size; cohort inference requires independently ascertained participants and suitable controls.

DAVID remains an optional, complementary native-web annotation route with the same foreground/background, separate native statistics, and stable-ID comparison. This protocol neither launches DAVID nor changes its integration. Absence of a P7 addendum does not prove an external DAVID session never occurred. Provider-shopping to manufacture significance is prohibited.
