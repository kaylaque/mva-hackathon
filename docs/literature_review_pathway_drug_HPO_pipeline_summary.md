# Reproducible Pathway–Drug Literature Review and Drug–HPO Evidence Pipeline

## 1. Purpose and scope

Build a repository pipeline that retrieves pathway/intervention literature, identifies every named drug, gathers drug information, and maps potential conflicts and therapeutic support to eight HPO phenotypes. Preserve evidence and screening decisions so the final matrices can be audited and regenerated.

Inputs: the 27 pathway names and synonyms, eight HPO codes, review eligibility criteria, database selection, publication cutoff, regulatory checking date, and country for affordability assessment. The original RTK/RAS/PI3K/mTOR topic can be configured as an additional pathway; keep it separate from the requested 27-column export unless the scope is expanded.

The earlier export contained 65 agents, but retrieval and review coverage were partial. Reproduction should retain that inventory as a seed, retry failed requests, discover additional agents from eligible papers, and avoid treating the earlier export as a completed PRISMA review.

## 2. End-to-end workflow

| Step | Action | Required saved result |
|---|---|---|
| 1. Configure | Version pathway/HPO dictionaries, eligibility rules, search cutoff and evidence-coding rules. | Configuration and protocol version |
| 2. Resolve terminology | Expand pathway names with validated synonyms and target terms; resolve official HPO labels/synonyms. | Dictionaries with database identifiers and provenance |
| 3. Search literature | Search pathway/synonym plus intervention terms. After drug discovery, search each drug/alias against each HPO independently. | Exact queries, database, date, total hits and pagination |
| 4. Retrieve | Retrieve all result IDs and metadata; fetch abstracts and available full texts. Retry transient failures and resume from checkpoints. | Source responses, retrieval status and checksums |
| 5. Deduplicate | Match PMID/DOI, then review title-based candidates. Link multiple reports of the same study separately from duplicate records. | Article, report and study identifiers |
| 6. Screen | Screen title/abstract and then full text against prespecified criteria. Retain every decision, reviewer and exclusion reason. | Screening ledger and unresolved disagreements |
| 7. Extract | Extract interventions, targets, mechanism, pathway position, effect direction, model/population, regimen, outcomes and limitations. | Structured claims linked to source passages |
| 8. Normalize agents | Resolve active ingredients and aliases, split combination components, preserve formulation differences and classify experimental agents. | Drug registry plus combination membership |
| 9. Map pathways | Distinguish direct targets, downstream effects, activation/inhibition, synthetic lethality and outside-pathway partners/comparators. | Drug–pathway evidence ledger |
| 10. Gather drug status | Check published clinical phase, current trial recruitment/status, termination reasons, FDA indications and commercial status separately. | Dated, indication-specific status records |
| 11. Assess HPO evidence | Extract potential conflicts and beneficial effects for each drug × HPO pair, including negative and contradictory findings. | Drug–HPO evidence ledger and assessment status |
| 12. Assess affordability | Check generic/specialty/investigational access; collect country-specific prices/coverage when available. | Formulation, country, currency/date, source and confidence |
| 13. Appraise | Apply design-appropriate risk-of-bias assessment and assess applicability. Complete independent review and resolve disagreements. | Domain judgments, rationale and review history |
| 14. Validate/export | Derive binary flags and notes from the reviewed evidence, calculate PRISMA counts, and generate reproducible outputs. | CSV/XLSX, dictionaries, logs and review report |

## 3. Skills and information sources

Skills guide implementation and interpretation. Repository code still needs API adapters, schemas, retrieval logic and review interfaces; skill names alone are not executable dependencies.

| Capability | Relevant skills / source |
|---|---|
| PubMed searches and metadata | `life-sciences-literature:ncbi-entrez-skill` |
| PMC access and article metadata | `life-sciences-literature:ncbi-pmc-skill`; retrieve only accessible full text and validate article identity |
| Pathway identifiers and position | `life-sciences-databases:reactome-skill` |
| Drug identity, mechanisms and targets | `life-sciences-databases:chembl-skill`, `pubchem-pug-skill`, `chebi-skill`, `opentargets-skill` |
| Experimental binding evidence | `life-sciences-databases:bindingdb-skill`, when target-binding evidence is needed |
| HPO terms | Official HPO ontology; OLS ontology lookup with the HPO namespace |
| Current clinical trials | `life-sciences-databases:clinicaltrials-skill` |
| FDA/commercial status and safety labels | FDA, Drugs@FDA, Orange Book where applicable, and DailyMed; verify the actual product and indication |
| Review methodology | `scholarresearch:find-sources`, `systematic-review`, `source-evaluation`, `evidence-synthesis` |
| Workbook export | `Spreadsheets` skill; CSV exports can be derived from validated structured tables |

Use primary studies for efficacy claims, official labels for label warnings, and registries for current trial status. Database annotations and review articles can guide discovery but should not replace the underlying evidence when available.

## 4. Search design

Build explicit Boolean queries from versioned term lists:

```text
(pathway name OR validated synonyms OR relevant target terms)
AND (drug OR inhibitor OR agonist OR antagonist OR intervention terms)

(drug name OR validated aliases)
AND (HPO phenotype label OR validated phenotype synonyms)
```

For each discovered drug, separately investigate pathway relationships, HPO conflicts, HPO therapeutic support, adverse effects/contraindications, clinical development, regulatory status, and affordability/access.

Retrieve every result page within the prespecified scope. Save hit counts and retrieved counts; never silently cap a systematic search at ten records. Record per-database coverage and failures. A complete systematic review requires the planned database coverage, citation searching and documented screening; searching PubMed alone does not automatically establish comprehensive coverage.

## 5. Proposed repository structure

```text
evidence_pipeline/
├── config/        # pathways, HPO terms, eligibility, sources, cutoff
├── schemas/       # extraction, screening and evidence contracts
├── src/
│   ├── search/
│   ├── retrieval/
│   ├── deduplication/
│   ├── screening/
│   ├── extraction/
│   ├── drug_normalization/
│   ├── pathway_mapping/
│   ├── regulatory/
│   ├── hpo_assessment/
│   ├── appraisal/
│   └── export/
├── data/
│   ├── raw/       # permitted source captures; immutable per run
│   ├── normalized/
│   └── reviewed/
├── tests/
└── outputs/
```

Use a relational database or structured tables as the authoritative data layer. Generate the wide CSV from those tables. Record the configuration hash, code commit, source checksums, extraction schema, model/prompt version and human edits for each run. Cache requests, honor API limits, retry with backoff, and make stages resumable. Store credentials outside the repository.

## 6. Core data tables

| Table | Minimum fields |
|---|---|
| `search_runs` | Query, database, dates, hit/retrieval counts, pagination, status/errors |
| `articles` | PMID, DOI, title, publication type, access and retraction status, linked study |
| `screening` | Article/report, stage, reviewer, decision, exclusion reason, disagreement resolution |
| `drugs` | Canonical ingredient, aliases, database IDs, agent type, formulations |
| `regimens` | Regimen ID, component drugs, roles, formulation and study context |
| `drug_pathway_evidence` | Drug, pathway, target, position, direction, relationship type, source passage |
| `drug_hpo_evidence` | Drug, HPO, conflict/support, outcome, design/model, context, limitations, source passage |
| `drug_status` | Drug, product, indication, published phase, current trial status, FDA/commercial status, source/date |
| `affordability` | Country, product/formulation, price/access evidence, currency/date, coverage, assessment basis |
| `risk_of_bias` | Study, tool, domains, judgments, rationale, reviewer and resolution |

Give every claim an evidence ID, source URL/identifier, passage locator and checked date. Validate PMID/DOI/title before assigning a full text to a paper. Keep unavailable evidence explicit.

## 7. Binary flags and notes

| Field group | Meaning of `1` | Meaning of `0` |
|---|---|---|
| 27 pathway columns | A documented relationship, with its directness and role specified | No qualifying relationship coded within the searched scope |
| 8 HPO-conflict columns | Evidence supports a potential concern in a stated context | No qualifying positive concern coded; may remain unassessed |
| 8 HPO-support columns | An article reports benefit for the phenotype or a stated related issue | No qualifying positive benefit coded; may remain unassessed |

**Zero is not proof of safety, lack of interaction or lack of efficacy.** Include per-pair statuses such as `not_searched`, `retrieval_failed`, `screening_pending`, `reviewed_no_qualifying_evidence`, and `evidence_coded`.

Keep three grouped notes columns: `pathway_notes`, `HPO_conflict_notes`, and `HPO_support_notes`. Each should identify evidence IDs, directness, population/model, relevant conditions, contradictory results and sources.

Conflict and support may both equal `1`. Distinguish human trials, observational studies, case reports, animal/cell studies, biochemical surrogates and mechanistic inference. A protocol alone cannot create a positive efficacy flag. A combination response cannot establish the independent contribution of each component. Preserve negative clinical findings alongside positive preclinical results. Do not sum binary flags into an unvalidated safety or efficacy ranking.

## 8. HPO dictionary and clinical context

| HPO code | Phenotype |
|---|---|
| HP:0002859 | Rhabdomyosarcoma |
| HP:0000121 | Nephrocalcinosis |
| HP:0004322 | Short stature |
| HP:0001508 | Failure to thrive |
| HP:0003202 | Skeletal muscle atrophy |
| HP:0001622 | Premature birth |
| HP:0001518 | Small for gestational age |
| HP:0200067 | Recurrent spontaneous abortion |

Distinguish the patient's phenotype from a relative's history, active disease from past disease, and maternal pregnancy exposure from treating a person who was born prematurely/SGA. Animal fetal loss or reduced fetal weight is an indirect reproductive concern, not proof of recurrent miscarriage or human SGA. Likewise, muscle injury is not identical to muscle atrophy, and renal toxicity is not identical to calcium deposition.

## 9. Drug information and appraisal rules

- Report published phase separately from the latest verified recruitment status; retain termination reasons and negative outcomes.
- Report FDA status for the actual product, indication and regimen. Approval for another disease does not establish HPO-specific benefit. A substance identifier or ChEMBL phase annotation is not FDA therapeutic approval.
- Identify activation versus inhibition and outside-pathway actions. Do not assume every pathway-associated drug restores a deficient pathway.
- Keep affordability country-, formulation- and date-specific. Separate verified prices/coverage from qualitative generic/specialty access estimates. Research-reagent prices are not clinical treatment costs.
- Use RoB 2 for randomized trials, ROBINS-I for appropriate nonrandomized intervention studies, and SYRCLE for animal studies. Apply an appropriate design-specific method for other evidence. Do not fabricate formal judgments from an abstract alone.
- Derive PRISMA counts from the screening ledger. Independent review, full-text exclusions and unavailable reports must remain visible. Meta-analysis requires sufficiently comparable studies; it is not an automatic pipeline step.

## 10. Validation and deliverables

Before export, check that retrieval counts reconcile with search results or explicitly report gaps; article identities match; aliases and combination components are handled correctly; and every positive flag traces to qualifying evidence. Check all 43 binary columns, required notes, reference integrity and CSV parsing. Calculate review counts from actual decisions and retain unresolved work.

Deliver:

```text
drug_pathway_HPO_matrix.csv
drug_pathway_evidence_ledger.csv
drug_HPO_evidence_ledger.csv
drug_HPO_cell_status.csv
drug_HPO_search_log.csv
drug_matrix_column_dictionary.csv
HPO_code_dictionary.csv
screening_ledger.csv
risk_of_bias.csv
PRISMA_report.md
README.md
```

Use an LLM for proposed screening and structured extraction, with source passages and schema validation. Human reviewers resolve ambiguous identities, clinical applicability, conflicting evidence and risk-of-bias judgments. Record what was automated versus reviewed.

## 11. Implementation priorities

| Priority | Deliverable | Completion gate |
|---|---|---|
| P0 | Configurations, retrieval, deduplication, normalization, evidence ledgers and CSV export | Resumable searches; explicit failures; every positive flag has provenance |
| P1 | Screening interface, HPO extraction, regulatory/development checks and affordability records | Decisions, contexts, dates and contradictions retained |
| P2 | Independent review, formal bias assessment and PRISMA reporting | Planned screening complete, disagreements resolved, exclusions and outstanding gaps disclosed |

This document is an implementation specification. It does not assert that these repository modules exist or that the earlier literature review is complete.
