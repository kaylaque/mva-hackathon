# Reproduce the HPO–pathway literature workflow in a repository

## What is implemented

The patch adds a repository-local Codex skill, a Python CLI, the actual eight-HPO and 53-pathway synonym rules, eight historical PubMed queries, a frozen scored 538-article CSV, source/checksum manifests, tests, and plotting code. It does not install a global skill, alter existing application code, or publish a website.

The CSV/report stages use Python 3.10+ and its standard library. Optional PNG figures require `matplotlib` and `numpy`. Live retrieval uses the installed Life Sciences Literature Python helpers. No GPU, vector database, model API, or LLM is needed to calculate the counts.

The included full scored CSV has article identifiers, metadata, flags, matched terms, and earlier annotations. Raw abstracts and full article text are not bundled. Retrieve/cache source text locally under the relevant source terms and keep it out of public Git history.

## 1. Import into your existing repo

Run from the repository root:

```bash
git apply --check /path/to/hpo_pathway_reproduction.patch
git apply /path/to/hpo_pathway_reproduction.patch
python -m unittest discover -s tests -p 'test_hpo_pathway*.py'
python scripts/hpo_pathway_pipeline.py report \
  --articles data/hpo_pathway/baseline/articles_scored.csv \
  --out runs/hpo_pathway/baseline \
  --expect-baseline
```

Check for path collisions before applying. If the repo has files with the same names, move this component into a subdirectory and update its ROOT-relative paths; do not overwrite unrelated files. The patch contains new files only.

For figures, install the plotting dependencies listed in `requirements-hpo-pathway.txt` in your repo's normal environment, then add `--plots` to the report command. Do not replace an existing dependency lock file. Pin and lock the resolved environment using your repo's existing dependency manager.

Invoke the local agent workflow with:

```text
Use $hpo-pathway-review to reproduce the saved 538-article corpus,
generate the article filters, all pathway ranks, pathways with >=1 HPO
coverage, and the top-10 HPO heatmap and exact-intersection UpSet plot.
Run the baseline checks and report any deviations.
```

The local skill is in `.agents/skills/hpo-pathway-review/`. A compatible agent discovers that directory; for another agent, load its `SKILL.md` as workflow instructions using that agent's documented configuration. The code itself does not depend on skill discovery.

## 2. Pipeline and stage contracts

| Stage | Inputs | Operation | Persistent outputs |
|---|---|---|---|
| Protocol | Eight HPO IDs; scope; dates | Specify sources, query profiles, eligible record types, screening policies, and extraction definitions | Versioned rules/config; protocol; reviewer forms |
| Search | Exact Boolean queries and supplemental seeds | Run every query; paginate; save reported count, returned IDs, timestamp, and provenance | Raw search responses; `search_log.json`; `registry.csv` |
| Retrieve | Deduplicated registry | Fetch article-level PubMed metadata; retain book records and unresolved placeholders | Raw XML; `records.json`; retrieval manifest |
| Full text | PMCID or DOI plus PMC metadata | Resolve version/license/availability; fetch available XML; extract body paragraphs | XML, body text, availability log, source hashes |
| Annotate | Title, abstract, indexing, available body; frozen rules | Match HPO/clinical/pathway mentions; preserve status, matched term, source, and coverage | `articles_scored.csv`; `article_pathway_evidence.csv` |
| Map | Explicit pathway/process phrase; curated mapping | Assign reference Reactome/KEGG IDs; label exact, broader, related-process, and complex mappings | Versioned dictionary, code/name links, mapping notes |
| Filter | Full scored table | Derive pathway-only and `(any HPO OR clinical mention) AND pathway` views | Filtered article CSVs; no deletion from full registry |
| Aggregate | Article–pathway membership; article HPO flags | Count distinct PMIDs; calculate each HPO count and distinct HPO coverage | All-pathway ranking; minimum-coverage ranking |
| Visualize | Top-N membership matrix | Count exclusive intersections; render ranks, HPO heatmap, UpSet | Region/member CSVs; PNG figures |
| Review | Search records, full texts, protocol | Record independent eligibility decisions, case overlap, and design-appropriate risk of bias | Human-reviewed decisions; exclusions; appraisal; PRISMA accounting |

Do not use the old 48-paper manual subset as an eligibility rule or an input filter. It was a selected narrative set, not an algorithmically derived systematic-review inclusion set.

## 3. Existing skills to connect

| Skill | Role in this pipeline | Boundary |
|---|---|---|
| `life-sciences-literature:ncbi-entrez-skill` | PubMed ESearch/EFetch, metadata and publication links | Use its `scripts/ncbi_entrez.py`; save raw responses for reproducibility; paginate all hits |
| `life-sciences-literature:ncbi-pmc-skill` | Resolve PMC version, availability, license, retraction flag, XML URLs | Its metadata result establishes availability; a PMCID by itself does not establish usable full text |
| `life-sciences-databases:reactome-skill` | Validate stable IDs, official names, reference species and pathway hierarchy | Use `scripts/rest_request.py`; do not infer article pathway activity from a gene name |
| `hpo-pathway-review` (included repo-local skill) | Orchestrate the stage commands, count definitions, regression checks and report generation | Delegate arithmetic to Python; never summarize away unreviewed records |
| A systematic-review protocol/appraisal skill, optionally `scholarresearch:systematic-review` | Help prepare eligibility rules, independent review, PRISMA and appraisal forms | A skill cannot substitute for completed independent screening or manufacture risk-of-bias judgments |
| `work-visualize:visualize` or `answers-charts:answers-charts` | Optional interactive exploration | Use Matplotlib outputs for exportable scientific figures; the CLI needs neither plugin |

Install or connect these through your agent's supported plugin/skill mechanism. Their names are not PyPI package names. Do not copy this session's private plugin-cache paths into your repo. Configure helper paths per machine instead.

## 4. Replay versus refresh

### Replay a frozen annotation snapshot

Use the included `data/hpo_pathway/baseline/articles_scored.csv`. This reproduces aggregation and visualizations exactly from the historical flags, without internet access.

Historical baseline checks:

| Quantity | Expected |
|---|---:|
| Full scored articles | 538 |
| Pathway/process mention articles | 130 |
| `(any HPO OR clinical mention) AND pathway` articles | 86 |
| Named pathway/process groups with >=1 article | 44 |
| Groups covering >=1 of the eight HPO features | 27 |
| Top-10 pathway union | 106 articles |
| Nonempty exact top-10 intersections | 38 |

The attached 130-row pathway-only file also reproduces the 44 groups and 27 HPO-covered groups:

```bash
python scripts/hpo_pathway_pipeline.py report \
  --articles /path/to/Pathway_mentioned_articles.csv \
  --out runs/hpo_pathway/from_pathway_subset --plots
```

Do not use `--expect-baseline` with the 130-row input: its input scope is smaller than 538. Its manifest correctly reports 130 input records. No count is padded to 538.

### Refresh searches and annotations

Configure these local values; the commands below use the installed skill helpers:

```bash
python scripts/hpo_pathway_pipeline.py retrieve \
  --config config/hpo_pathway/rules.json \
  --entrez-helper /path/to/ncbi-entrez-skill/scripts/ncbi_entrez.py \
  --out runs/hpo_pathway/refresh/cache

python scripts/hpo_pathway_pipeline.py fulltext \
  --records runs/hpo_pathway/refresh/cache/records.json \
  --pmc-helper /path/to/ncbi-pmc-skill/scripts/ncbi_pmc.py \
  --out runs/hpo_pathway/refresh/fulltext

python scripts/hpo_pathway_pipeline.py annotate \
  --records runs/hpo_pathway/refresh/cache/records.json \
  --config config/hpo_pathway/rules.json \
  --fulltext-dir runs/hpo_pathway/refresh/fulltext \
  --out runs/hpo_pathway/refresh/articles_scored.csv

python scripts/hpo_pathway_pipeline.py report \
  --articles runs/hpo_pathway/refresh/articles_scored.csv \
  --out runs/hpo_pathway/refresh/reports --plots
```

The shipped queries retain the historical publication cutoff, 2026-10-20. To search through a new date, make a new config version and update every date bound. The baseline snapshot is not modified. The CLI uses a new output folder to keep runs separate; if reusing a folder, cached responses are deliberately reused. Query caches are addressed by the query's identifier, so never change a query and reuse its old cache directory.

For a fixed PMID set with newly fetched metadata, add `--registry data/hpo_pathway/baseline/registry.csv` to `retrieve`. That skips searching and fetches that registry's PMIDs.

The CLI waits 0.4 seconds between uncached Entrez/helper calls and checks page totals. Queries exceeding 9,999 hits stop with an explicit error: partition them by publication date and rerun in a new run folder. For production, implement recursive nonoverlapping date partitions, shared throttling across processes, retries/backoff, and a documented handling policy for exceptionally large one-day result sets. Never accept a truncated first page as the full corpus.

The live retrieval code has offline parser/helper-contract tests and passed one live PubMed metadata fetch for PMID 21552266. It was not validated by rerunning every historical search or downloading every full text during delivery. Service behavior, helper revisions, article indexing and availability can change. Keep `runs/hpo_pathway/`, raw caches, environment files and credentials out of Git using your repo's existing ignore rules.

### Why a refreshed run need not equal 538

The original eight searches yielded 551 hits before deduplication, or 533 unique PubMed IDs. Eleven `PubmedBookArticle` records were skipped by the original parser. Sixteen supplement-only articles were added: `533 - 11 + 16 = 538`. The specific IDs and provenance are preserved in the config and registry.

The refreshed parser retains book records. With otherwise identical search results, it could therefore retain 549 records rather than 538. This is an intentional correction to record accounting, not a new clinical eligibility filter. A formal review must decide eligible document types in its protocol and record any exclusion explicitly.

The historical queries were phenotype-cluster searches with several MVA/gene branches. They are not a complete, disease-neutral search of every literature item associated with any one of the eight HPO terms. Use a separately versioned phenotype-first query profile, controlled vocabulary, additional databases, and citation chasing for broader coverage. Do not silently add those searches to the historical baseline.

## 5. Scoring definitions

HPO columns follow this fixed order:

| Column | ID | Feature |
|---|---|---|
| HP1 | HP:0002859 | Rhabdomyosarcoma |
| HP2 | HP:0000121 | Nephrocalcinosis |
| HP3 | HP:0004322 | Short stature |
| HP4 | HP:0001508 | Failure to thrive |
| HP5 | HP:0003202 | Skeletal muscle atrophy |
| HP6 | HP:0001622 | Premature birth |
| HP7 | HP:0001518 | Small for gestational age |
| HP8 | HP:0200067 | Recurrent spontaneous abortion |

For historical reproduction, an HPO mention flag is 1 if a direct configured term/synonym occurs, including a negated/uncertain mention, or an exact MeSH topic is indexed. Keep the richer status separately: positive, mixed, negated, related term only, not found, or insufficient text. A zero is not proof of absence. `any_nonnegated_HP_topic_0_1` is a separate flag.

The clinical score is a publication-type or patient/case/cohort phrase mention. Reviews and background mentions can score 1. It does not establish an original case, a subject count, or that all mentioned features co-occur in one person.

Pathway inclusion is 1 for a literal `pathway/pathways` word or a configured named pathway/process/complex. Gene-only matches are not added. Generic pathway mentions remain in the article subset with no invented identifier. Only named groups enter the pathway ranking.

Species flags are article-level indexing/text mentions. A human gene name, human cells, human disease background, and a human clinical case are different concepts. The legacy species heuristic and clinical heuristic are preserved for reproduction, not treated as adjudicated study design.

Guard against these common errors in a revised classifier:

- Premature chromatid separation is not premature birth.
- Hypercalciuria is related to, not identical to, nephrocalcinosis.
- IUGR/low birth weight is not automatically small for gestational age.
- `RMS` can name Rabson–Mendenhall syndrome, not rhabdomyosarcoma.
- Maternal pregnancy-loss history is not automatically a phenotype of the child.
- Homologous recombination used to construct a knockout mouse is a methods mention, not necessarily a disease mechanism.

For improved clinical/mechanistic extraction, add `assertion`, `experiencer`, `species`, `section`, and `context` (`background`, `methods`, `result`, `discussion`) to each evidence item. Add entity/offset anchors and manual adjudication. Release this as a new rule version and compare results to the historical version rather than relabeling the baseline.

## 6. Exact aggregation formulas

For pathway label p and HPO feature h:

```text
number_of_articles(p) = count(distinct PMID mentioning p)
HPh(p) = count(distinct PMID mentioning p AND HPh flag == 1)
number_of_articles_with_any_HP(p) = count(distinct PMID mentioning p AND any HPO flag == 1)
number_of_HP_covered(p) = sum(HPh(p) > 0 for h in 1..8)
keep(p) = number_of_HP_covered(p) >= 1
```

An article can contribute to several pathway rows and several HPO columns. Do not add those columns to obtain a unique article count. Rank by distinct article count descending; break ties alphabetically by pathway name. Top ten therefore includes RNA splicing and excludes the tied proteasomal group. No pathway labels are merged solely because they share a broader Reactome/KEGG parent.

The CLI produces an UpSet figure for the ten pathway sets, with all nonempty exact intersections, a membership CSV, and HPO counts per intersection. Use a Venn only for a selected two- or three-set comparison. A ten-circle sketch cannot communicate all 1,023 possible nonempty combinations reliably.

## 7. Pathway IDs and DAVID

`rules.json` contains the actual pathway dictionary and frozen verified ID/name pairs. Validate IDs using the Reactome skill; a compact hierarchy request is:

```json
{"base_url":"https://reactome.org/ContentService","path":"data/eventsHierarchy/9606","headers":{"Accept":"application/json"},"save_raw":true,"raw_output_path":"/absolute/run/path/reactome_hierarchy.json","max_items":1}
```

Validate individual mappings against `data/query/<stableId>` and check the official name, event class, species, and mapping relationship. Store response versions/dates and hashes. Refresh the KEGG reference list from `https://rest.kegg.jp/list/pathway/hsa` in a new dictionary version.

`REACTOME_PATHWAY:R-HSA-...` and `KEGG_PATHWAY:hsa...` are category-compatible reference codes for use with DAVID outputs. They are not separate DAVID pathway ontologies, DAVID gene IDs, evidence that a term is present in the current DAVID release, or enrichment results. Live DAVID membership was not checked. Broader/related references must stay labeled; nonhuman studies are not automatically evidence for the human reference pathway.

The code ships a frozen dictionary; it does not perform gene-list enrichment. If enrichment is later added, require an explicit foreground/background, species, provider/database version, tested-term count and multiple-testing correction. Do not convert article mention frequency to an enrichment p-value.

## 8. Review, PRISMA, duplicate cases and risk of bias

The report command writes an `UNASSESSED` review-decision template. Complete title/abstract and full-text decisions with independent reviewers and reasoned exclusions; preserve disagreements and resolution. Match the appraisal tool to each study design. Leave missing full texts and unassessed risk of bias visible. Do not infer a low-risk judgment from a publication type.

To cross-check repeated cases, compare authors/centres, recruitment windows, patient/family identifiers, variants, age/sex, prior citations and phenotype profiles. Record `SAME_CASE`, `POSSIBLE_OVERLAP`, `DISTINCT_CASE_SUPPORTED`, or `UNRESOLVED`, with evidence. PMID deduplication removes duplicate records, not duplicate patients.

Report search accounting separately from clinical eligibility. The prior 538, 130, 86 and 27 values describe a corpus and automated filters; they are not a completed PRISMA inclusion flow. Do not label a run a full systematic review until screening, appraisal and reporting are actually complete.

## 9. Validation performed and limits

The delivery checks the attached 130-row file, the frozen 538-row snapshot, all filter memberships, 44/27 group counts, each HPO count, top-ten tie handling, and every UpSet marginal/union. The annotation stage is compared against the prior flags using the locally available original metadata and 17 full-text extracts. Those source texts are not included in the patch; their source/version/hash manifest is included.

Replaying the scored CSV reproduces aggregation exactly. Recreating every annotation from the internet requires the same source versions, normalization, rule version and text coverage. Fetching every currently available full text changes coverage and can change the results. Human case/phenotype adjudication and risk-of-bias assessment remain unperformed unless supplied by reviewers.

## 10. Primary references

- PubMed E-utilities parameters: https://www.ncbi.nlm.nih.gov/books/NBK25499/
- NLM E-utilities usage guidance: https://www.nlm.nih.gov/dataguide/eutilities/utilities.html
- Reactome ContentService: https://reactome.org/dev/content-service
- DAVID annotation terminology: https://davidbioinformatics.nih.gov/helps/functional_annotation.html
- KEGG human pathway reference: https://rest.kegg.jp/list/pathway/hsa

## Implementation prompt for your coding agent

```text
Use $hpo-pathway-review. Integrate the supplied component into this repo's
existing conventions without replacing unrelated configuration or lock files.
First run the frozen baseline replay and tests. Keep search, retrieval,
annotation, mapping, aggregation, plotting and human review as separate stages.
Preserve every retrieved record and exact query/source provenance. Make new
queries, rules, references and text-coverage changes create new run versions.
Never force fresh results to 538, treat missing text as phenotype absence,
equate patient mentions with confirmed cases, invent pathway identifiers,
or describe an automated screen as completed PRISMA review.
Improve retries/date-partitioned retrieval and introduce context-aware
evidence extraction as a separately evaluated rule version. Validate against
a manually reviewed stratified sample, including ambiguous and negative cases.
Return reproducible commands, artifact paths, tests and remaining limitations.
```
