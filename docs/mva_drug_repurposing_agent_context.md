# MVA Dry-Lab Drug Repurposing — Agent Context

**Updated:** 2026-10-12 · **Scope:** WGS/VCF + HPO → mechanism-aware therapeutic *hypotheses*, not a proven treatment.

## 1. Concepts and decision objective

- **Drug discovery (dry lab):** identify/validate a biological target, find or design active compounds, then computationally assess activity, selectivity, developability, and risk. *New compounds* typically require extensive optimization and experimental validation.
- **Drug repurposing (repositioning):** find a *new indication* for an existing/previously studied compound. Prior knowledge of pharmacology can shorten exploration, but **new-disease efficacy and safety are not established**.
- **Pathway analysis:** ask which annotated biological processes are overrepresented among candidate genes. It proposes *biological context*, **not** directly measured pathway activation or drug sensitivity.
- **Network pharmacology:** connect disease genes, functional partners, drug targets, and compounds; direct or indirect network proximity is **not** itself a mechanism of rescue.

**Goal:** Produce a small, auditable set of hypotheses with the chain **patient genetic evidence → probable gene/protein perturbation → affected process → actionable target → correctly directed drug action → feasibility/safety → proposed rescue experiment**. Distinguish facts, computational predictions, assumptions, and unknowns in every output.

## 2. End-to-end workflow

```text
Patient GRCh38 WGS/VCF + HPO
  → QC/normalization + VEP consequence annotation
  → ClinVar / gnomAD / pathogenicity & splice evidence
  → inheritance + HPO-based variant/gene prioritization
  → gene-level foreground and eligible-gene background
  → g:Convert (symbols → ENSG) → g:GOSt (GO/pathway enrichment; BH-FDR)
  → Reactome pathway members + UniProt / QuickGO functions
  → STRING functional network + mechanism-oriented query panel
  → direct drugs (Open Targets / ChEMBL / DGIdb)
    + indirect/pathway drugs (same sources via supported neighbor targets)
    + optional perturbation-signature and knowledge-graph hypotheses
  → merge/deduplicate → direction-of-effect / evidence / safety gate
  → ranked DRUG–TARGET–MECHANISM hypotheses + validation plan
```

### 2.1. Genetic evidence / handoff (Phase 1)

| Task | Tools and how used | Interpret with care |
|---|---|---|
| Variant QC, normalization | [bcftools](https://samtools.github.io/bcftools/bcftools.html): filter calls, normalize alleles, inspect genotype/read support. | Low-confidence calls and phase ambiguity can mislead downstream results. |
| Functional consequences | [Ensembl VEP](https://mart.ensembl.org/info/docs/tools/vep/index.html): transcript-specific coding, splicing and regulatory consequences. | `missense` / `frameshift` alone does not establish actual protein dysfunction. |
| Rarity and clinical evidence | [gnomAD](https://gnomad.broadinstitute.org/) population allele frequency; [ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/intro/) assertions **and** [review status](https://www.ncbi.nlm.nih.gov/clinvar/docs/review_status/). | Rarity ≠ pathogenicity; ClinVar classifications are not patient-specific diagnoses. |
| In silico evidence | [CADD](https://cadd.gs.washington.edu/), [REVEL](https://sites.google.com/site/revelgenomics/), [SpliceAI paper](https://doi.org/10.1016/j.cell.2018.12.015), [AlphaMissense paper](https://doi.org/10.1126/science.adg7492): consequence-relevant supporting predictions. | Correlated predictors cannot be counted as independent confirmations. |
| Phenotype / inheritance | [Exomiser](https://exomiser.readthedocs.io/en/stable/index.html): VCF + HPO matching, allele frequency, variant pathogenicity, inheritance. | For recessive genes inspect biallelic evidence and phase; single VCF cannot prove aneuploidy or protein function. |

**MVA biological anchors:** [`BUB1B`, `CEP57`, `TRIP13`](https://www.jci.org/articles/view/120316) relate to chromosome segregation/mitotic checkpoint biology, **but different genes can cause different cellular defects**. Do not assume all perturbations have the same direction or respond to the same compound. WGS/VCF + HPO do **not** directly measure protein abundance, pathway activity, chromosome-missegregation rate, or drug rescue.

### 2.2. Pathway / network branch (Phase 2; reported implementation)

| Step | Role / output | Tool / documentation |
|---|---|---|
| P0–P1: Input | Read Phase-1 gene matrix; retain gene-level columns; define foreground (selected genes) and background (all *eligible/evaluated* genes). | Preserve source row counts, filters, and gene IDs. |
| P2: Identifier harmonization | Batch symbols → Ensembl ENSG; save mapped, unmapped, duplicated, and one-to-many cases. | [g:Profiler g:Convert](https://biit.cs.ut.ee/gprofiler/page/docs) / [API](https://biit.cs.ut.ee/gprofiler/page/apis). |
| P3: Enrichment | Query GO biological process/molecular function/cellular component and relevant pathways. Use **custom background** and provider-supplied **Benjamini–Hochberg FDR** (`significance_threshold_method='fdr'`). | [g:GOSt docs](https://biit.cs.ut.ee/gprofiler/page/docs) / [2023 paper](https://doi.org/10.1093/nar/gkad347). |
| P4: Prioritize terms | Retain `term_id`, name, adjusted p/q, foreground overlap, term size, background, matched genes, and version. Collapse redundant parent/child terms for summaries. | Never recompute BH on already adjusted provider values. |
| P5: Functional interpretation | Fetch pathway members/reactions and curated protein functions/GO evidence; build plausible variant → mechanism → process statements. | [Reactome](https://reactome.org/userguide/analysis), [Reactome API](https://reactome.org/dev/analysis), [UniProt](https://www.uniprot.org/), [QuickGO API](https://www.ebi.ac.uk/QuickGO/api/). |
| P6: Therapeutic expansion | Find compounds for causal or functionally supported *alternative* pathway targets; record molecule, action type, target and underlying source. | [Open Targets](https://platform-docs.opentargets.org/getting-started), [ChEMBL API](https://chembl.gitbook.io/chembl-interface-documentation/web-services), [DGIdb](https://dgidb.org/api). |
| P7: Network comparison | Compare pathway-selected targets with the pre-existing STRING component/query panel; label **direct gene**, **pathway member**, or **network neighbor**. | [STRING score/edge semantics](https://string-db.org/help/scores/); [STRING threshold guide](https://en.string-db.org/help/getting_started/). |

**Reported project state (not independently re-audited here):** Phase 1 evaluated **8,876 genes**, **1,985 passed**; the STRING phase reported **477 high-confidence edges/components**, with an MVA-focused panel at `18_phase3_query_gene_panel.csv`. The g:Profiler → Reactome/QuickGO/UniProt → Open Targets/ChEMBL/DGIdb → STRING-comparison branch outputs private P0–P7 tables, Markdown summary, and provenance. **DAVID was not run**; [DAVID](https://davidbioinformatics.nih.gov/documentation/) is an *optional comparison* enrichment engine, not interchangeable results. Confirm exact counts/files from current local artifacts before asserting success.

**Key distinctions:** A significant GO/Reactome result shows *overrepresentation*, not activation/inhibition. A high STRING score is confidence in a **functional association**, not binding affinity, druggability, causal direction or effect magnitude. [Reactome overrepresentation](https://reactome.org/userguide/analysis/) and [g:GOSt enrichment](https://biit.cs.ut.ee/gprofiler/page/docs) are complementary analyses, not direct protein-interaction tests.

## 3. Three drug-candidate generation routes

| Route | How to do it | Evidence and failure mode |
|---|---|---|
| **A. Direct target (primary)** | Input strongest disease gene/protein → [ChEMBL](https://chembl.gitbook.io/chembl-interface-documentation/web-services) (bioactivity, potency, assay conditions, action), [DGIdb](https://dgidb.org/api) (aggregated interactions + source), [Open Targets](https://platform-docs.opentargets.org/data-access/graphql-api) (target–disease, clinical, tractability/safety). | Best provenance for a direct mechanism **if drug action corrects the defect**. A database drug–gene edge alone does not confirm efficacy. |
| **B. Indirect / network compensation** | Mechanism → [Reactome](https://reactome.org/) relevant process → [STRING](https://string-db.org/help/scores/) partners → identify biologically justified upstream/downstream or compensatory targets → query the **same** drug databases. | Record edge evidence, actual functional role, hypothesized compensation and uncertainty; proximity/hub degree ≠ rescue. |
| **C. Perturbation / AI (optional)** | With a **valid disease expression signature**, seek opposite perturbations in [LINCS/CMap (L1000)](https://doi.org/10.1016/j.cell.2017.10.049), or use pathway-oriented [Gene2Drug](https://doi.org/10.1093/bioinformatics/btx800). Explore [PrimeKG](https://doi.org/10.1038/s41597-023-01960-3) for interpretable typed graph paths. | Patient VCF is **not** RNA-seq or a measured disease-expression signature. Predicted network/perturbation links are exploratory and require independent validation. |

**Dry-lab discovery extension (optional, not a substitute for repurposing):** When a specific target and activity assay are credible, inspect active chemical series in ChEMBL; use [RDKit](https://www.rdkit.org/docs/GettingStartedInPython.html) for standardized structures, fingerprints/similarity and basic properties; screen developability with [ADMETlab 3.0](https://doi.org/10.1093/nar/gkae236). Docking/MD is only a **structural plausibility** check for a defined target and site, not proof of target engagement, cellular correction, or safety. [DrugBank API](https://docs.drugbank.com/discovery/v1/) can add mechanism/interactions if licensed/authenticated; do not assume unrestricted access.

## 4. Interpretation and ranking: questions the agent must answer

For **every** candidate, fill these fields explicitly:

1. **Genetic anchor:** Which prioritized variant/gene? How strong are genotype, inheritance, phenotype and read-level evidence?
2. **Defect hypothesis:** Is the effect likely *loss of function*, *gain of function*, altered dosage, mis-splicing, or **unknown**? Cite variant- or function-specific studies, not just gene descriptions.
3. **Pathway evidence:** Which pathway and which overlapping genes? Report background, adjusted FDR, overlap, and whether the annotation supports *membership* or a *directional* process change.
4. **Target opportunity:** Is drug target **direct**, **pathway-based**, or merely **network-associated**? What compensatory rationale connects it to the defect?
5. **Action/direction:** Is drug an inhibitor, agonist, antagonist, activator, etc.? Would its pharmacology *plausibly restore* the dysregulated function—or worsen it? If unknown, flag **DIRECTION UNRESOLVED**.
6. **Pharmacology:** Existing indication/clinical stage, target potency/assay, selectivity, feasible exposure/tissue distribution, known serious risks and interactions; [Open Targets](https://platform-docs.opentargets.org/), [ChEMBL](https://www.ebi.ac.uk/chembl/), and regulatory labeling if relevant.
7. **Evidence grade:** Separate (a) curated clinical/experimental evidence, (b) model-organism/cellular evidence, (c) mechanism analogy, (d) database association, and (e) ML predictions. Do not add overlapping sources as independent votes.
8. **Validation:** Name a *specific falsifiable rescue endpoint*: e.g., chromosome segregation errors/lagging chromosomes, micronuclei, spindle-checkpoint timing/response, viable cell fraction, and dose-response vs toxicity. Match assay to the **specific gene mechanism**.

**Suggested gates (not an opaque numerical score):** `genetic_plausibility`, `mechanistic_direction`, `target_evidence`, `drug_action_evidence`, `exposure_safety`, `independent_validation`. **Reject or quarantine** hypotheses whose pharmacology plausibly aggravates chromosome instability, whose target action is unknown, or whose claimed benefit is supported only by enrichment or graph proximity. Label remaining candidates **high / medium / exploratory mechanistic priority**, never "effective treatment".

## 5. QC checklist and required output

- **Input QC:** genome build, gene universe/background, Ensembl/symbol mappings, duplicate IDs, HPO set, inheritance hypotheses, missing SV/CNV/aneuploidy evidence, variant provenance.
- **Enrichment QC:** foreground ⊆ background, mapped-gene counts, gene-set sizes/overlaps, FDR method, source/version/date, terms with too few supporting genes, correlated GO hierarchy.
- **Network QC:** organism, score cutoff and edge type; avoid interpreting functional STRING associations as proven physical contact or directional regulation.
- **Drug QC:** normalize compound and target IDs; record per-source evidence, target action, molecule status, assay vs prediction, and safety; avoid double-counting ChEMBL records syndicated through Open Targets/DGIdb.
- **Privacy:** Keep patient variant rows, genotype, exact clinical text, coordinates and identifiers local. The reported pathway branch sent gene symbols/ENSG foreground/background to g:Profiler; re-check each external request before execution.

**Deliverables:** `pathway_prioritization.csv`, `target_mechanism_evidence.csv`, `drug_candidates_evidence.csv`, and a short `repurposing_hypotheses.md`. Minimum candidate record:

```text
candidate_id | gene/variant | defect_hypothesis | pathway_id | pathway_FDR |
therapeutic_target | direct_or_indirect | drug_id/name | drug_action |
expected_rescue_direction | evidence_level + URLs | approval/clinical_stage |
exposure_safety_flags | uncertainties | rescue_assay | verdict
```

## 6. Core papers / methodological exemplars

| Reference | Use in this project |
|---|---|
| [g:Profiler update (2023)](https://doi.org/10.1093/nar/gkad347) | Identifier normalization and controlled enrichment with explicit statistical background. |
| [DGIdb 5.0 (2024)](https://doi.org/10.1093/nar/gkad1040) | Evidence-aware retrieval of drug–gene interactions; inspect original sources and direction. |
| [Connectivity Map / L1000 (2017)](https://doi.org/10.1016/j.cell.2017.10.049) | Perturbation-signature strategy **only when expression evidence exists**. |
| [Gene2Drug (2018)](https://doi.org/10.1093/bioinformatics/btx800) | Pathway-oriented drug ranking using compound-induced transcriptomic signatures; not proof of drug rescue. |
| [Sirota et al. (2011)](https://doi.org/10.1126/scitranslmed.3001318) | Expression-signature reversal followed by preclinical testing; illustrates why computational rank needs a validation layer. |
| [Dudley et al. (2011)](https://doi.org/10.1126/scitranslmed.3002648) | Drug repositioning from disease/drug signatures with experimental follow-up. |
| [PrimeKG (2023)](https://doi.org/10.1038/s41597-023-01960-3) | Typed, provenance-oriented disease–gene–pathway–drug graph context (exploratory). |
| [CEP57/MVA mechanism study (2018)](https://www.jci.org/articles/view/120316) | MVA gene-specific cellular defects; caution against treating all MVA genotypes as one phenotype. |

**Implementation / source repositories:** [Exomiser](https://github.com/exomiser/Exomiser) (phenotype-guided variant ranking); [Reactome Analysis Service](https://github.com/reactome/analysis-service) (enrichment API); [DGIdb v5](https://github.com/dgidb/dgidb-v5) (drug–gene integration); [Open Targets Platform API](https://github.com/opentargets/platform-api) (target/drug graph queries); [RDKit](https://github.com/rdkit/rdkit) (chemical representation and fingerprints). These are *reference implementations*, not claims that the project currently runs their code.

## 7. Agent operating instructions

1. Use existing Phase-1 and P0–P7 outputs as **inputs**, not as proof of causality or successful therapy; first audit their actual contents and version.
2. Work **gene/variant → mechanism → pathway → target → compound**, never just “gene list → popular drug”.
3. Cite original experimental papers and database record URLs for every proposed link. Mark results as **observed**, **curated**, **inferred**, or **speculative**.
4. If only WGS/VCF + HPO are available, do **not** claim measured transcription, protein defects, pathway activation, or cell rescue; list which experiment would decide the question.
5. Favor a **few directionally coherent, independently supported, testable candidates** over a long unfiltered compound list. Never suggest a patient treatment regimen.
