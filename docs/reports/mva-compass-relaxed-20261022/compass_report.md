# Compass report — Benjamini-010

Research-exploration evidence integration. Nothing here is a diagnosis, a clinical
recommendation, or experimental proof. Annotation is not efficacy; membership is not
activity; an approved drug is not an approved treatment for this disease.

## Lanes

| Lane | Meaning | Targets |
| --- | --- | --- |
| PATHWAY_SUPPORTED | provider-adjusted significance with an observed panel contributor | 0 |
| MECHANISM_PRIOR | finite MVA / SAC / cohesin panel researched from public gene function | 34 |
| NETWORK_CONTEXT | saved network membership only, never independent confirmation | 5 |

A target may carry several labels; one shared source is never independent replication.

## Pathway evidence

### RELAXED DAVID BENJAMINI 0.10 EXPLORATION
- saved annotations evaluated: 16773 eligibility rows
- stable-ID DAVID input terms: 208
- DAVID terms at validated cutoff (<=0.05): 0
- DAVID terms at earlier exploratory cutoff (<=0.10): 1; g:Profiler remains at its saved 0.05 rule.
- newly in DAVID cutoff (0.05,0.10]: 1; eligible after gates: 0
- DAVID terms above 0.10: 207; unknown statistic: 0; held or context-gated at cutoff: 1
- held rationales: 31
- distinct exploratory targets queried: 0
- drug rows returned by available public sources: 0; source incompleteness is listed below; candidates are not treatment recommendations
- Only DAVID terms with Benjamini <=0.10 could pass the statistical gate; context, observed membership, and MVA relevance gates also apply.
- Exploratory threshold changes do not strengthen evidence or imply treatment benefit.

### Held and stopped annotations

- STOP_NOT_SIGNIFICANT: 16662 annotation(s)
- STOP_UNMAPPED: 111 annotation(s)
- example: g:Profiler GO GO:0000012: STOP_NOT_SIGNIFICANT (ABOVE_PROVIDER_THRESHOLD)
- example: g:Profiler GO GO:0000018: STOP_NOT_SIGNIFICANT (ABOVE_PROVIDER_THRESHOLD)
- example: g:Profiler GO GO:0000027: STOP_NOT_SIGNIFICANT (ABOVE_PROVIDER_THRESHOLD)
- example: g:Profiler GO GO:0000038: STOP_NOT_SIGNIFICANT (ABOVE_PROVIDER_THRESHOLD)
- example: g:Profiler GO GO:0000041: STOP_NOT_SIGNIFICANT (ABOVE_PROVIDER_THRESHOLD)
- example: g:Profiler GO GO:0000045: STOP_NOT_SIGNIFICANT (ABOVE_PROVIDER_THRESHOLD)
- example: g:Profiler GO GO:0000050: STOP_NOT_SIGNIFICANT (ABOVE_PROVIDER_THRESHOLD)
- example: g:Profiler GO GO:0000070: STOP_NOT_SIGNIFICANT (ABOVE_PROVIDER_THRESHOLD)
- example: g:Profiler GO GO:0000075: STOP_NOT_SIGNIFICANT (ABOVE_PROVIDER_THRESHOLD)
- example: g:Profiler GO GO:0000076: STOP_NOT_SIGNIFICANT (ABOVE_PROVIDER_THRESHOLD)
- example: g:Profiler GO GO:0000077: STOP_NOT_SIGNIFICANT (ABOVE_PROVIDER_THRESHOLD)
- example: g:Profiler GO GO:0000079: STOP_NOT_SIGNIFICANT (ABOVE_PROVIDER_THRESHOLD)
- example: g:Profiler GO GO:0000082: STOP_NOT_SIGNIFICANT (ABOVE_PROVIDER_THRESHOLD)
- example: g:Profiler GO GO:0000086: STOP_NOT_SIGNIFICANT (ABOVE_PROVIDER_THRESHOLD)
- example: g:Profiler GO GO:0000096: STOP_NOT_SIGNIFICANT (ABOVE_PROVIDER_THRESHOLD)
- example: g:Profiler GO GO:0000098: STOP_NOT_SIGNIFICANT (ABOVE_PROVIDER_THRESHOLD)
- example: g:Profiler GO GO:0000101: STOP_NOT_SIGNIFICANT (ABOVE_PROVIDER_THRESHOLD)
- example: g:Profiler GO GO:0000122: STOP_NOT_SIGNIFICANT (ABOVE_PROVIDER_THRESHOLD)
- example: g:Profiler GO GO:0000132: STOP_NOT_SIGNIFICANT (ABOVE_PROVIDER_THRESHOLD)
- example: DAVID GO GO:0000132: STOP_NOT_SIGNIFICANT (ABOVE_PROVIDER_THRESHOLD)
- 16753 further stopped annotation(s) are listed in `pathway_eligibility.csv`

## Mechanism-prior research

- curated protein-role rows: 2515
- reviewed mechanism assertions: 9
- unreviewed literature search hits retained for audit: 291
- symbols resolved through the HGNC alias index: none
- identities left unresolved rather than chosen: MPS1

## Decisions

- ADVANCE_RESEARCH (experimental prioritisation only): 0
- TEST_FIRST: 2
- HOLD: 31
- REJECT: 1

- AURKA: HOLD | source=MECHANISM_PRIOR | direction=UNKNOWN | desired=UNKNOWN | reasons=MECHANISM_PRIOR_PANEL;NO_PATHWAY_ENRICHMENT_SUPPORT;AWAITING_ROLE_EVIDENCE
- AURKB: HOLD | source=MECHANISM_PRIOR | direction=UNKNOWN | desired=UNKNOWN | reasons=MECHANISM_PRIOR_PANEL;NO_PATHWAY_ENRICHMENT_SUPPORT;AWAITING_ROLE_EVIDENCE
- BUB1: HOLD | source=MECHANISM_PRIOR | direction=UNKNOWN | desired=UNKNOWN | reasons=MECHANISM_PRIOR_PANEL;NO_PATHWAY_ENRICHMENT_SUPPORT;AWAITING_ROLE_EVIDENCE
- BUB1B: HOLD | source=MECHANISM_PRIOR|NETWORK_CONTEXT | direction=UNKNOWN | desired=UNKNOWN | reasons=MECHANISM_PRIOR_PANEL;NO_PATHWAY_ENRICHMENT_SUPPORT;DIRECTION_UNKNOWN
- BUB3: HOLD | source=MECHANISM_PRIOR | direction=UNKNOWN | desired=UNKNOWN | reasons=MECHANISM_PRIOR_PANEL;NO_PATHWAY_ENRICHMENT_SUPPORT;AWAITING_ROLE_EVIDENCE
- CDC20: HOLD | source=MECHANISM_PRIOR | direction=UNKNOWN | desired=UNKNOWN | reasons=MECHANISM_PRIOR_PANEL;NO_PATHWAY_ENRICHMENT_SUPPORT;AWAITING_ROLE_EVIDENCE
- CENATAC: TEST_FIRST | source=MECHANISM_PRIOR | direction=REDUCED_FUNCTION | desired=RESTORE_NORMAL | reasons=MECHANISM_PRIOR_PANEL;NO_PATHWAY_ENRICHMENT_SUPPORT;CONTEXTUAL_HYPOTHESIS_WITH_PERTURBATION
- CENPE: HOLD | source=MECHANISM_PRIOR | direction=UNKNOWN | desired=UNKNOWN | reasons=MECHANISM_PRIOR_PANEL;NO_PATHWAY_ENRICHMENT_SUPPORT;AWAITING_ROLE_EVIDENCE
- CEP57: HOLD | source=MECHANISM_PRIOR | direction=REDUCED_FUNCTION | desired=RESTORE_NORMAL | reasons=MECHANISM_PRIOR_PANEL;NO_PATHWAY_ENRICHMENT_SUPPORT;MODEL_UNSPECIFIED
- KNL1: HOLD | source=MECHANISM_PRIOR | direction=UNKNOWN | desired=UNKNOWN | reasons=MECHANISM_PRIOR_PANEL;NO_PATHWAY_ENRICHMENT_SUPPORT;AWAITING_ROLE_EVIDENCE
- KNTC1: HOLD | source=MECHANISM_PRIOR | direction=UNKNOWN | desired=UNKNOWN | reasons=MECHANISM_PRIOR_PANEL;NO_PATHWAY_ENRICHMENT_SUPPORT;AWAITING_ROLE_EVIDENCE
- MAD1L1: HOLD | source=MECHANISM_PRIOR|NETWORK_CONTEXT | direction=ALTERED_FUNCTION | desired=UNKNOWN | reasons=MECHANISM_PRIOR_PANEL;NO_PATHWAY_ENRICHMENT_SUPPORT;DIRECTION_UNKNOWN
- MAD2L1: HOLD | source=MECHANISM_PRIOR | direction=UNKNOWN | desired=UNKNOWN | reasons=MECHANISM_PRIOR_PANEL;NO_PATHWAY_ENRICHMENT_SUPPORT;AWAITING_ROLE_EVIDENCE
- MAD2L2: HOLD | source=MECHANISM_PRIOR|NETWORK_CONTEXT | direction=UNKNOWN | desired=UNKNOWN | reasons=MECHANISM_PRIOR_PANEL;NO_PATHWAY_ENRICHMENT_SUPPORT;AWAITING_ROLE_EVIDENCE
- MPS1: HOLD | source=MECHANISM_PRIOR | direction=UNKNOWN | desired=UNKNOWN | reasons=MECHANISM_PRIOR_PANEL;NO_PATHWAY_ENRICHMENT_SUPPORT;SOURCE_INCOMPLETE
- NCAPD2: HOLD | source=MECHANISM_PRIOR|NETWORK_CONTEXT | direction=UNKNOWN | desired=UNKNOWN | reasons=MECHANISM_PRIOR_PANEL;NO_PATHWAY_ENRICHMENT_SUPPORT;AWAITING_ROLE_EVIDENCE
- NCAPD3: HOLD | source=MECHANISM_PRIOR | direction=UNKNOWN | desired=UNKNOWN | reasons=MECHANISM_PRIOR_PANEL;NO_PATHWAY_ENRICHMENT_SUPPORT;AWAITING_ROLE_EVIDENCE
- NCAPG: HOLD | source=MECHANISM_PRIOR | direction=UNKNOWN | desired=UNKNOWN | reasons=MECHANISM_PRIOR_PANEL;NO_PATHWAY_ENRICHMENT_SUPPORT;AWAITING_ROLE_EVIDENCE
- NCAPH: HOLD | source=MECHANISM_PRIOR | direction=UNKNOWN | desired=UNKNOWN | reasons=MECHANISM_PRIOR_PANEL;NO_PATHWAY_ENRICHMENT_SUPPORT;AWAITING_ROLE_EVIDENCE
- NDC80: HOLD | source=MECHANISM_PRIOR | direction=UNKNOWN | desired=UNKNOWN | reasons=MECHANISM_PRIOR_PANEL;NO_PATHWAY_ENRICHMENT_SUPPORT;AWAITING_ROLE_EVIDENCE
- 14 further rationale(s) are listed in `target_mechanism_rationale.csv`

## Counterevidence and safety unknowns

- TTK: cited hazard or opposite action (RV-008)
- Unknown safety is recorded as unknown; it is never reported as safe.

## Test plans

- RAT-3ba05fb06d5f0de4 [EVIDENCE_RESOLUTION_REQUIRED] Resolve MECHANISM_PRIOR_PANEL;NO_PATHWAY_ENRICHMENT_SUPPORT;AWAITING_ROLE_EVIDENCE for AURKA before proposing any experiment
- RAT-d6e04012e738bf77 [EVIDENCE_RESOLUTION_REQUIRED] Resolve MECHANISM_PRIOR_PANEL;NO_PATHWAY_ENRICHMENT_SUPPORT;AWAITING_ROLE_EVIDENCE for AURKB before proposing any experiment
- RAT-d1ce2a6897133fd2 [EVIDENCE_RESOLUTION_REQUIRED] Resolve MECHANISM_PRIOR_PANEL;NO_PATHWAY_ENRICHMENT_SUPPORT;AWAITING_ROLE_EVIDENCE for BUB1 before proposing any experiment
- RAT-927b6231a4071893 [EVIDENCE_RESOLUTION_REQUIRED] Resolve MECHANISM_PRIOR_PANEL;NO_PATHWAY_ENRICHMENT_SUPPORT;DIRECTION_UNKNOWN for BUB1B before proposing any experiment
- RAT-089a7c57dec721e6 [EVIDENCE_RESOLUTION_REQUIRED] Resolve MECHANISM_PRIOR_PANEL;NO_PATHWAY_ENRICHMENT_SUPPORT;AWAITING_ROLE_EVIDENCE for BUB3 before proposing any experiment
- RAT-c285e8c5ebc1bf13 [EVIDENCE_RESOLUTION_REQUIRED] Resolve MECHANISM_PRIOR_PANEL;NO_PATHWAY_ENRICHMENT_SUPPORT;AWAITING_ROLE_EVIDENCE for CDC20 before proposing any experiment
- RAT-73fd735891f71d9d [PROPOSED] RESTORE_NORMAL CENATAC function restores checkpoint or segregation competence in the cited model; this is a research hypothesis, not a treatment proposal
- RAT-8efd2ae1e9dfce84 [EVIDENCE_RESOLUTION_REQUIRED] Resolve MECHANISM_PRIOR_PANEL;NO_PATHWAY_ENRICHMENT_SUPPORT;AWAITING_ROLE_EVIDENCE for CENPE before proposing any experiment
- RAT-6acb2af0662baa96 [EVIDENCE_RESOLUTION_REQUIRED] Resolve MECHANISM_PRIOR_PANEL;NO_PATHWAY_ENRICHMENT_SUPPORT;MODEL_UNSPECIFIED for CEP57 before proposing any experiment
- RAT-fa88af46b6b0bdf0 [EVIDENCE_RESOLUTION_REQUIRED] Resolve MECHANISM_PRIOR_PANEL;NO_PATHWAY_ENRICHMENT_SUPPORT;AWAITING_ROLE_EVIDENCE for KNL1 before proposing any experiment
- RAT-efbee30dc3b63ed0 [EVIDENCE_RESOLUTION_REQUIRED] Resolve MECHANISM_PRIOR_PANEL;NO_PATHWAY_ENRICHMENT_SUPPORT;AWAITING_ROLE_EVIDENCE for KNTC1 before proposing any experiment
- RAT-3a924af92e569c7e [EVIDENCE_RESOLUTION_REQUIRED] Resolve MECHANISM_PRIOR_PANEL;NO_PATHWAY_ENRICHMENT_SUPPORT;DIRECTION_UNKNOWN for MAD1L1 before proposing any experiment
- RAT-a23da88e3b6b40a9 [EVIDENCE_RESOLUTION_REQUIRED] Resolve MECHANISM_PRIOR_PANEL;NO_PATHWAY_ENRICHMENT_SUPPORT;AWAITING_ROLE_EVIDENCE for MAD2L1 before proposing any experiment
- RAT-dc756c9d7f2848d4 [EVIDENCE_RESOLUTION_REQUIRED] Resolve MECHANISM_PRIOR_PANEL;NO_PATHWAY_ENRICHMENT_SUPPORT;AWAITING_ROLE_EVIDENCE for MAD2L2 before proposing any experiment
- RAT-7992837c1a4e5dfd [EVIDENCE_RESOLUTION_REQUIRED] Resolve MECHANISM_PRIOR_PANEL;NO_PATHWAY_ENRICHMENT_SUPPORT;SOURCE_INCOMPLETE for MPS1 before proposing any experiment
- RAT-2479fec050ecfe4c [EVIDENCE_RESOLUTION_REQUIRED] Resolve MECHANISM_PRIOR_PANEL;NO_PATHWAY_ENRICHMENT_SUPPORT;AWAITING_ROLE_EVIDENCE for NCAPD2 before proposing any experiment
- RAT-6642240b68ce6c93 [EVIDENCE_RESOLUTION_REQUIRED] Resolve MECHANISM_PRIOR_PANEL;NO_PATHWAY_ENRICHMENT_SUPPORT;AWAITING_ROLE_EVIDENCE for NCAPD3 before proposing any experiment
- RAT-1b335a4ab1c8f574 [EVIDENCE_RESOLUTION_REQUIRED] Resolve MECHANISM_PRIOR_PANEL;NO_PATHWAY_ENRICHMENT_SUPPORT;AWAITING_ROLE_EVIDENCE for NCAPG before proposing any experiment
- RAT-5cf361198d9ec8b3 [EVIDENCE_RESOLUTION_REQUIRED] Resolve MECHANISM_PRIOR_PANEL;NO_PATHWAY_ENRICHMENT_SUPPORT;AWAITING_ROLE_EVIDENCE for NCAPH before proposing any experiment
- RAT-0b679420d1fcc582 [EVIDENCE_RESOLUTION_REQUIRED] Resolve MECHANISM_PRIOR_PANEL;NO_PATHWAY_ENRICHMENT_SUPPORT;AWAITING_ROLE_EVIDENCE for NDC80 before proposing any experiment
- RAT-40b5e9013ea8add6 [EVIDENCE_RESOLUTION_REQUIRED] Resolve MECHANISM_PRIOR_PANEL;NO_PATHWAY_ENRICHMENT_SUPPORT;AWAITING_ROLE_EVIDENCE for NIPBL before proposing any experiment
- RAT-3075344bd6ec6a16 [EVIDENCE_RESOLUTION_REQUIRED] Resolve MECHANISM_PRIOR_PANEL;NO_PATHWAY_ENRICHMENT_SUPPORT;AWAITING_ROLE_EVIDENCE for PLK1 before proposing any experiment
- RAT-60abf3f4aa042b75 [EVIDENCE_RESOLUTION_REQUIRED] Resolve MECHANISM_PRIOR_PANEL;NO_PATHWAY_ENRICHMENT_SUPPORT;AWAITING_ROLE_EVIDENCE for RAD21 before proposing any experiment
- RAT-76f84e28044c3203 [PROPOSED] RESTORE_NORMAL SLF2 function restores checkpoint or segregation competence in the cited model; this is a research hypothesis, not a treatment proposal
- RAT-5b8a95029977c817 [EVIDENCE_RESOLUTION_REQUIRED] Resolve MECHANISM_PRIOR_PANEL;NO_PATHWAY_ENRICHMENT_SUPPORT;AWAITING_ROLE_EVIDENCE for SMC1A before proposing any experiment
- RAT-4cec2c7fc383c3f4 [EVIDENCE_RESOLUTION_REQUIRED] Resolve MECHANISM_PRIOR_PANEL;NO_PATHWAY_ENRICHMENT_SUPPORT;AWAITING_ROLE_EVIDENCE for SMC2 before proposing any experiment
- RAT-868cceada88efacc [EVIDENCE_RESOLUTION_REQUIRED] Resolve MECHANISM_PRIOR_PANEL;NO_PATHWAY_ENRICHMENT_SUPPORT;AWAITING_ROLE_EVIDENCE for SMC3 before proposing any experiment
- RAT-68170c6f899ffa50 [EVIDENCE_RESOLUTION_REQUIRED] Resolve MECHANISM_PRIOR_PANEL;NO_PATHWAY_ENRICHMENT_SUPPORT;AWAITING_ROLE_EVIDENCE for SMC4 before proposing any experiment
- RAT-c74fadf39a1ace29 [EVIDENCE_RESOLUTION_REQUIRED] Resolve MECHANISM_PRIOR_PANEL;NO_PATHWAY_ENRICHMENT_SUPPORT;MODEL_UNSPECIFIED for SMC5 before proposing any experiment
- RAT-a74e82c1bef3bbfe [EVIDENCE_RESOLUTION_REQUIRED] Resolve MECHANISM_PRIOR_PANEL;NO_PATHWAY_ENRICHMENT_SUPPORT;AWAITING_ROLE_EVIDENCE for STAG2 before proposing any experiment
- RAT-fc45cd9adeb0ec50 [EVIDENCE_RESOLUTION_REQUIRED] Resolve MECHANISM_PRIOR_PANEL;NO_PATHWAY_ENRICHMENT_SUPPORT;AWAITING_ROLE_EVIDENCE for TRIP13 before proposing any experiment
- RAT-5bfe4a3425508f3d [EVIDENCE_RESOLUTION_REQUIRED] Resolve MECHANISM_PRIOR_PANEL;NO_PATHWAY_ENRICHMENT_SUPPORT;CITED_CONTEXTUAL_HAZARD for TTK before proposing any experiment
- RAT-8e8c7aa9dd9b3086 [EVIDENCE_RESOLUTION_REQUIRED] Resolve MECHANISM_PRIOR_PANEL;NO_PATHWAY_ENRICHMENT_SUPPORT;AWAITING_ROLE_EVIDENCE for ZW10 before proposing any experiment
- RAT-9bcfee64c431987c [EVIDENCE_RESOLUTION_REQUIRED] Resolve MECHANISM_PRIOR_PANEL;NO_PATHWAY_ENRICHMENT_SUPPORT;AWAITING_ROLE_EVIDENCE for ZWILCH before proposing any experiment

## Network context

- status: NETWORK_CONTEXT_AVAILABLE (projected membership reconstructed)
- saved network projections are LEGACY_UNVERIFIED and are never independent causal confirmation

## Source completeness

| Stage | Source | Scope | Status | Complete |
| --- | --- | --- | --- | --- |
| gene_identity | HGNC | AURKA | OK | True |
| gene_identity | HGNC | AURKB | OK | True |
| gene_identity | HGNC | BUB1 | OK | True |
| gene_identity | HGNC | BUB1B | OK | True |
| gene_identity | HGNC | BUB3 | OK | True |
| gene_identity | HGNC | CDC20 | OK | True |
| gene_identity | HGNC | CENATAC | OK | True |
| gene_identity | HGNC | CENPE | OK | True |
| gene_identity | HGNC | CEP57 | OK | True |
| gene_identity | HGNC | KNL1 | OK | True |
| gene_identity | HGNC | KNTC1 | OK | True |
| gene_identity | HGNC | MAD1L1 | OK | True |
| gene_identity | HGNC | MAD2L1 | OK | True |
| gene_identity | HGNC | MAD2L2 | OK | True |
| gene_identity | HGNC | MPS1 | AMBIGUOUS | False |
| gene_identity | HGNC | NCAPD2 | OK | True |
| gene_identity | HGNC | NCAPD3 | OK | True |
| gene_identity | HGNC | NCAPG | OK | True |
| gene_identity | HGNC | NCAPH | OK | True |
| gene_identity | HGNC | NDC80 | OK | True |
| gene_identity | HGNC | NIPBL | OK | True |
| gene_identity | HGNC | PLK1 | OK | True |
| gene_identity | HGNC | RAD21 | OK | True |
| gene_identity | HGNC | SLF2 | OK | True |
| gene_identity | HGNC | SMC1A | OK | True |
| gene_identity | HGNC | SMC2 | OK | True |
| gene_identity | HGNC | SMC3 | OK | True |
| gene_identity | HGNC | SMC4 | OK | True |
| gene_identity | HGNC | SMC5 | OK | True |
| gene_identity | HGNC | STAG2 | OK | True |
| gene_identity | HGNC | TRIP13 | OK | True |
| gene_identity | HGNC | TTK | OK | True |
| gene_identity | HGNC | ZW10 | OK | True |
| gene_identity | HGNC | ZWILCH | OK | True |
| protein_role | UniProt | AURKA | OK | True |
| protein_role | UniProt | AURKB | OK | True |
| protein_role | UniProt | BUB1 | OK | True |
| protein_role | UniProt | BUB1B | OK | True |
| protein_role | UniProt | BUB3 | OK | True |
| protein_role | UniProt | CDC20 | OK | True |
| protein_role | UniProt | CENATAC | OK | True |
| protein_role | UniProt | CENPE | OK | True |
| protein_role | UniProt | CEP57 | OK | True |
| protein_role | UniProt | KNL1 | OK | True |
| protein_role | UniProt | KNTC1 | OK | True |
| protein_role | UniProt | MAD1L1 | OK | True |
| protein_role | UniProt | MAD2L1 | OK | True |
| protein_role | UniProt | MAD2L2 | OK | True |
| protein_role | UniProt | MPS1 | AMBIGUOUS | False |
| protein_role | UniProt | NCAPD2 | OK | True |
| protein_role | UniProt | NCAPD3 | OK | True |
| protein_role | UniProt | NCAPG | OK | True |
| protein_role | UniProt | NCAPH | OK | True |
| protein_role | UniProt | NDC80 | OK | True |
| protein_role | UniProt | NIPBL | OK | True |
| protein_role | UniProt | PLK1 | OK | True |
| protein_role | UniProt | RAD21 | OK | True |
| protein_role | UniProt | SLF2 | OK | True |
| protein_role | UniProt | SMC1A | OK | True |
| protein_role | UniProt | SMC2 | OK | True |
| protein_role | UniProt | SMC3 | OK | True |
| protein_role | UniProt | SMC4 | OK | True |
| protein_role | UniProt | SMC5 | OK | True |
| protein_role | UniProt | STAG2 | OK | True |
| protein_role | UniProt | TRIP13 | OK | True |
| protein_role | UniProt | TTK | OK | True |
| protein_role | UniProt | ZW10 | OK | True |
| protein_role | UniProt | ZWILCH | OK | True |
| literature | PubMed | AURKA | OK | False |
| literature | PubMed | AURKB | OK | False |
| literature | PubMed | BUB1 | OK | False |
| literature | PubMed | BUB1B | OK | False |
| literature | PubMed | BUB3 | OK | False |
| literature | PubMed | CDC20 | OK | False |
| literature | PubMed | CENATAC | OK | True |
| literature | PubMed | CENPE | OK | False |
| literature | PubMed | CEP57 | OK | False |
| literature | PubMed | KNL1 | OK | False |
| literature | PubMed | KNTC1 | OK | False |
| literature | PubMed | MAD1L1 | OK | False |
| literature | PubMed | MAD2L1 | OK | False |
| literature | PubMed | MAD2L2 | OK | False |
| literature | PubMed | MPS1 | OK | False |
| literature | PubMed | NCAPD2 | OK | False |
| literature | PubMed | NCAPD3 | OK | False |
| literature | PubMed | NCAPG | OK | False |
| literature | PubMed | NCAPH | OK | False |
| literature | PubMed | NDC80 | OK | False |
| literature | PubMed | NIPBL | OK | False |
| literature | PubMed | PLK1 | OK | False |
| literature | PubMed | RAD21 | OK | False |
| literature | PubMed | SLF2 | OK | False |
| literature | PubMed | SMC1A | OK | False |
| literature | PubMed | SMC2 | OK | False |
| literature | PubMed | SMC3 | OK | False |
| literature | PubMed | SMC4 | OK | False |
| literature | PubMed | SMC5 | OK | False |
| literature | PubMed | STAG2 | OK | False |
| literature | PubMed | TRIP13 | OK | False |
| literature | PubMed | TTK | OK | False |
| literature | PubMed | ZW10 | OK | False |
| literature | PubMed | ZWILCH | OK | False |

## Rejection ledger summary

- STOP_NOT_SIGNIFICANT: 16662
- HOLD: 2015
- STOP_UNMAPPED: 111
- REJECT: 1

## Output links

- [compass_provenance.json](./compass_provenance.json)

## Limitations

- Provider statistics remain provider-specific; a single provider is judged on its own corrected values.
- Term relevance is never established from a title; a held term needs a cited pathway role.
- Public protein-role and modification annotation describe annotation, not measured activity.
- A literature search hit is unreviewed until a reviewed assertion cites its source section.
- Insufficient or failing source access is reported as such and never converted into a negative result.
- No patient direction, variant support or clinical recommendation is produced by this layer.
