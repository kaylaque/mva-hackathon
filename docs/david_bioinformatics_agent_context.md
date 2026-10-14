# DAVID Bioinformatics — Concise Agent Context

**Scope:** Functional annotation, pathway enrichment, clustering, and integration into a gene → pathway → target → drug *hypothesis* pipeline. **Updated:** 2026-10-12. **Status:** Methodology reference, not a report of an actual patient/sample run.

## 1. What DAVID does (and does not do)

**DAVID** (Database for Annotation, Visualization and Integrated Discovery) takes a gene/protein list, resolves identifiers into unified **DAVID Genes**, connects them to biological annotations, tests which terms are **over-represented** against a background, and groups redundant terms or related genes. Annotation sources include GO, pathways (e.g. KEGG/WikiPathways and other supported sources), protein domains, disease associations, and selected drug/gene annotations. It provides **biological context**, **not** gene causality, pathway activity/direction, demonstrated protein dysfunction, or a validated drug treatment.

- [DAVID official website](https://davidbioinformatics.nih.gov/) · [2026 tool documentation](https://davidbioinformatics.nih.gov/documentation/) · [2026 update paper](https://doi.org/10.1093/nar/gkag470) · [2022 knowledgebase update](https://doi.org/10.1093/nar/gkac194).

## 2. Inputs and preparation

| Input | What to supply | Critical check |
|---|---|---|
| **Foreground / gene list** | Meaningfully selected genes (e.g. a prioritized gene matrix); Ensembl, Entrez, gene symbols, protein IDs as supported | Unique genes, species/taxonomy, current symbols, unmapped IDs, duplicated IDs, evidence for selection |
| **Background / universe** | All genes **eligible to have been selected** by the upstream study (where supported), or suitable genome-wide background for genuinely genome-wide selection | Wrong background distorts enrichment; compare mapped/effective universes between tools |
| **Annotation categories** | GO:BP/MF/CC; selected pathways, functional domains, disease annotations | Record database/source, knowledgebase date, ontology version |
| **Analysis configuration** | Species; EASE/count thresholds; adjusted-p-value method; clustering stringency | Keep all parameters, defaults, software/DB versions in provenance |

**For WGS/VCF-derived MVA work:** Do **not** treat every variant-bearing gene as a disease-associated gene. Define a biologically justified foreground from variant quality, rarity/impact, inheritance, HPO, and candidate evidence; preserve the full evaluated/mappable gene universe. Conventional over-representation analysis (ORA) may be biased by gene length, coverage, annotation density and how the candidate list was selected. A pathway from a candidate matrix is a **hypothesis** rather than patient-specific pathway activity.

**Data handling:** DAVID is an external service. Send only necessary de-identified gene identifiers and selected universe—not VCF/genotypes, patient identifiers, or clinical narratives.

## 3. End-to-end DAVID workflow

```text
Input gene matrix + eligible-gene universe
    → select evidence-based foreground; deduplicate; normalize IDs/species
    → DAVID ID mapping + mapping-loss QC; submit foreground/background
    → choose annotation categories
    → Functional Annotation Chart: term-level enrichment (EASE, BH, hits)
    → Functional Annotation Clustering: group redundant terms (Kappa)
    → Gene Functional Classification: group similar genes (Kappa)
    → Functional Annotation Table / Gene Report / pathway viewer
    → export terms + groups + original contributing genes + provenance
    → compare with g:Profiler; investigate mechanisms in Reactome/QuickGO
    → inspect STRING/UniProt and drug-target evidence; form hypotheses
```

### DAVID modules: methods → outputs

| Module | How it works, briefly | Useful output |
|---|---|---|
| **Gene ID conversion / List Manager** | Links heterogeneous identifiers to internal DAVID Gene IDs; tracks organism | Input→DAVID mappings; mapped vs unmapped; deduplicated list |
| **Functional Annotation Table** | Gene-to-term lookup, **no enrichment statistics** | Gene × annotation records (function, domain, pathway, disease) |
| **Functional Annotation Chart** | One-sided over-representation with **EASE** (modified Fisher exact) for each term | Term/source, contributing genes, count, %, EASE p, Benjamini/BH, fold enrichment (if selected) |
| **Functional Annotation Clustering** | Term–term **Kappa similarity** based on shared genes, then **fuzzy heuristic clustering** | Biological-theme clusters, member terms/genes, group enrichment score |
| **Gene Functional Classification** | Gene–gene Kappa based on shared annotations, then fuzzy clustering | Sets of functionally similar genes + consensus terms |
| **Gene Report / Gene Search** | Gene-centric integrated annotation | Functions, identifiers, proteins/domains, annotations, outgoing source links |
| **Pathway Viewer / Related Terms** | Displays mapped list members on pathway maps and related annotation structures | Pathway context and gene overlays, **not** direction of pathway activation |
| **DAVID Ortholog** (optional) | Maps genes to orthologs via OMA/Ensembl Compara | Cross-species functional hints; not necessary for human MVA |

Sources: [Chart](https://davidbioinformatics.nih.gov/documentation/tools/chartReport.html), [Term clustering](https://davidbioinformatics.nih.gov/documentation/tools/term2term.html), [Gene classification](https://davidbioinformatics.nih.gov/documentation/tools/gene2gene.html), [2026 update](https://doi.org/10.1093/nar/gkag470).

## 4. Core algorithms and what the numbers mean

**A. Over-representation analysis (ORA):** For each term, form a 2×2 table of membership in the foreground vs eligible background and membership in the annotation term. Test whether that term is observed in **more** selected genes than expected under a null of no enrichment. [DAVID gene-enrichment documentation](https://davidbioinformatics.nih.gov/documentation/tools/tech-gene-enrichment.html).

**B. EASE score:** DAVID's modified **one-sided Fisher exact p-value**, conservatively recalculating enrichment after subtracting **one** foreground hit from the annotated-term count. Terms supported by only one gene get EASE p=1. A low EASE p is evidence of over-representation, **not** a causal probability or pathway-activity score. The usual Chart default (`EASE ≤ 0.1`, `count ≥ 2`) is a **display filter**, not a sufficient final significance criterion. [Algorithm](https://davidbioinformatics.nih.gov/helps/functional_annotation.html).

**C. Multiple testing:** Many GO and pathway terms are assessed, so inspect **`Benjamini` (BH-adjusted p)** (e.g. ≤0.05 under a prespecified decision rule). DAVID also exposes a separately named **`FDR`** column based on an adaptive method; **do not confuse it with `Benjamini`**. Never pool DAVID p-values with g:Profiler p-values as if from one test family. [DAVID statistics](https://davidbioinformatics.nih.gov/documentation/tools/tech-gene-enrichment.html) · [Benjamini–Hochberg 1995](https://doi.org/10.1111/j.2517-6161.1995.tb02031.x).

**D. Functional Annotation Clustering:** For each pair of terms, compute **Kappa** from common/uncommon assigned genes; fuzzy multiple-linkage heuristics group related terms, potentially with overlapping memberships. A **cluster enrichment score** is roughly `−log10(geometric mean of its member EASE p-values)`. It ranks clusters but is **NOT** a BH-adjusted p-value, independent validation, or a direct drug-target score. [Clustering technical documentation](https://davidbioinformatics.nih.gov/documentation/tools/tech-clustering.html).

**E. Gene Functional Classification:** Applies related Kappa/fuzzy clustering to **genes sharing annotation terms** rather than to terms sharing genes. Useful to find modules within a large candidate list; not a physical protein–protein interaction (PPI) network. [Gene classification](https://davidbioinformatics.nih.gov/documentation/tools/gene2gene.html).

## 5. If you submit 1,000+ genes: what happens?

**Example (illustrative, not experimental data):** Upload 1,200 selected human gene IDs → resolve/deduplicate and record mapping loss → map each unique gene to many GO/pathway/domain annotations → test each eligible annotation term against the selected background → output potentially hundreds/thousands of term rows (depending on list, filters, and knowledgebase) → group highly overlapping terms into interpretable clusters and optionally classify genes into functional groups.

**Expected characteristics:**

1. **Large gene list ≠ 1,000 pathways.** One gene belongs to multiple terms; many genes map to the same term; a term may have few or many hits. Results are **terms enriched in the list**, not automatically causal pathways for each gene.
2. More input genes can yield numerous broad/generic terms (`signal transduction`, `cytoplasm`, etc.). Inspect **specificity, term size, fold enrichment, hit identity, and background**; favor coherent mechanisms rather than ranking solely by tiny p-values.
3. Raw submitted IDs ≠ unique mapped genes. Preserve mapping counts, unmapped genes, species ambiguity, one-to-many mappings, and term coverage; enrichment statistics use internally normalized unique gene entities.
4. Term redundancy is expected. **Functional Annotation Clustering** compresses similar terms but does **not** generate new independent evidence or prove one pathway controls another.
5. **Limits:** DAVID documentation says it is designed/tested for lists of **≤3,000 genes**; **Functional Annotation Clustering and Gene Functional Classification have a 3,000-gene limit**. A 1,000–3,000-gene list is ordinarily in scope, though performance depends on the server. **Do not arbitrarily split a list** merely because it exceeds 1,000: separate tests change significance and biological meaning. Above documented limits, use reasoned modules/subsets or a more scalable alternative, while reporting changed hypotheses/universes. [Official limits](https://davidbioinformatics.nih.gov/documentation/).
6. **Submission method:** DAVID's lightweight URL-based API is intended only for small lists (2026 official tutorial: ~≤400 genes; not suitable for >500); use DAVID UI or its stateful **Web Service** for 1,000+ genes, after checking its current interface/auth and request constraints. Avoid assuming the newer internal REST endpoints are equivalent to an unrestricted public batch API. [2026 tutorial](https://davidbioinformatics.nih.gov/documentation/quick_tutorial.03012026.pdf) · [DAVID Web Service](https://davidbioinformatics.nih.gov/webservice/) · [SOAP operations](https://davidbioinformatics.nih.gov/webservice/services/listServices) · [DAVID-WS paper](https://doi.org/10.1093/bioinformatics/bts251).

**If there are no significant terms:** Do not force a pathway hypothesis. Review foreground selection, unmapped genes, category/term coverage, EASE/count display filters, universe correctness, and statistical power. Report the null result and optionally retain *labeled exploratory* annotation themes.

## 6. How to integrate with your existing g:Profiler branch

**Recommended architecture: DAVID = secondary functional-interpretation branch; g:Profiler = primary ORA branch.**

```text
Candidate gene matrix + same biologically eligible universe
            ├── g:Profiler: GO/Reactome enrichment, corrected statistics
            └── DAVID: EASE/BH enrichment + term & gene clustering
                          ↓
            Normalize source IDs, GO IDs / Reactome IDs,
            contributor genes and mapped universe
                          ↓
            Annotate: shared | g:Profiler-only | DAVID-only
                          ↓
            Reactome/QuickGO/UniProt/STRING mechanism checks
                          ↓
            ChEMBL / Open Targets / DGIdb drug evidence
                          ↓
            Testable drug hypotheses with evidence tier and direction
```

- Match organism, foreground definition, and **intended eligible background**; record if each provider resolves or filters genes differently. DAVID EASE/BH and g:Profiler enrichment tests can legitimately differ due to coverage, universe, test, correction and database versions.
- **g:Profiler's default correction is g:SCS**, not BH. Set its FDR option explicitly if comparing BH-like decisions; even then p-values should remain **provider-specific**. [g:Profiler manual](https://biit.cs.ut.ee/gprofiler/page/docs) · [API docs](https://biit.cs.ut.ee/gprofiler/page/apis).
- Concordance = **interpretive consistency**, not independent statistical replication: GO/pathway annotations overlap between tools. Use DAVID-exclusive terms as *exploratory*, particularly if DAVID `Benjamini > 0.05`.
- **Do not gate direct variant→protein→drug evidence on pathway FDR.** An independently supported, directional gene–drug hypothesis can exist despite weak pathway enrichment. Conversely, significant enrichment alone cannot approve a drug.

### Suggested machine-readable output (minimum)

| Artifact | Essential fields |
|---|---|
| `mapping_qc.csv` | submitted_id, canonical_gene_id, species, mapped, ambiguous, duplicate_group, reason |
| `david_terms.csv` | term_id, term_name, category/source, hit_genes, hit_count, background_hits, EASE_p, Benjamini_p, fold_enrichment, tested_background, source_version |
| `david_term_clusters.csv` | cluster_id, member_term_ids, contributing_genes, enrichment_score, term_count, representative_theme; **not** `cluster_fdr` |
| `david_gene_clusters.csv` | cluster_id, member_genes, shared_terms, stringency |
| `pathway_evidence_integrated.csv` | stable_term_id, provider, source-specific_padj, genes, shared_status, mechanism_note, evidence_tier, caveats |
| `run_provenance.json` | foreground/universe counts & checksums, mapping loss, ID system, organism, date, DAVID/g:Profiler version, annotations, thresholds, exports |

*If some fields are not returned by an endpoint, mark them `null`/unavailable; do not infer or manufacture values.*

## 7. Related tools — when they enter the workflow

| Tool | Main role relative to DAVID | Official reference / code |
|---|---|---|
| **g:Profiler** | Independent-looking but overlapping functional enrichment and ID conversion; primary statistics in this pipeline | [Docs](https://biit.cs.ut.ee/gprofiler/page/docs) · [GitHub](https://github.com/biitgroup/gprofiler-official) |
| **QuickGO / GO** | Resolve GO term definitions, ontology parents/children, annotation evidence | [QuickGO API](https://www.ebi.ac.uk/QuickGO/api/index.html) · [Gene Ontology](https://geneontology.org/) |
| **Reactome** | Explore reaction/event hierarchy, pathway participants and directionally meaningful biochemical mechanism | [Content Service](https://reactome.org/ContentService/) · [Reactome](https://reactome.org/) |
| **STRING** | Functional association network, neighborhood/context; evidence channels ≠ necessarily physical binding | [STRING](https://string-db.org/) · [API](https://string-db.org/cgi/help?subpage=api) |
| **UniProt / InterPro** | Function, protein domain, variant/protein annotations to interpret possible mechanism | [UniProt](https://www.uniprot.org/) · [InterPro](https://www.ebi.ac.uk/interpro/) |
| **Ensembl / g:Convert** | Resolve stable gene IDs, symbols, aliases, species mismatches | [Ensembl REST](https://rest.ensembl.org/) · [g:Convert](https://biit.cs.ut.ee/gprofiler/convert) |
| **GSEA / fgsea** | If all genes have meaningful **continuous ranked scores**, prefer rank-based enrichment over arbitrary cutoff-only ORA, where appropriate | [GSEA](https://www.gsea-msigdb.org/gsea/index.jsp) · [fgsea GitHub](https://github.com/alserglab/fgsea) |
| **clusterProfiler / EnrichmentMap** | Reproducible R enrichment/comparison and network visualization of redundant enriched terms | [clusterProfiler](https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html) · [EnrichmentMap docs](https://enrichmentmap.readthedocs.io/) |
| **Open Targets / ChEMBL / DGIdb** | *Downstream* target–disease/compound activity and mechanism evidence; not DAVID enrichment validation | [Open Targets API](https://platform-docs.opentargets.org/data-access/graphql-api) · [ChEMBL API](https://chembl.gitbook.io/chembl-interface-documentation/web-services/chembl-data-web-services) · [DGIdb](https://www.dgidb.org/) |

## 8. Agent interpretation rules (mandatory)

1. Label each statement **annotation**, **significantly enriched** (provider-adjusted p ≤ prespecified cutoff), **exploratory hypothesis**, or **mechanistically supported**. These are different claims.
2. Preserve input/background, resolved gene counts, annotated hits, raw EASE, **Benjamini BH**, source/version, and overlapping genes; **never** report cluster score as an adjusted p-value.
3. For a pathway with BH ≤0.05, still ask: *Which genes drive it? Is it generic? Which biochemical step is affected? What variant/protein mechanism suggests increase/decrease of function? Is there any independent drug-target/action evidence?*
4. For BH >0.05 or DAVID-only terms, put hypothesis in a **separate exploratory table**, not the significant-results table.
5. **Do not** equate term enrichment, protein interaction, causal variant, pathway activation, target direction, and therapeutic rescue. A drug hit is not an evidence-backed treatment without mechanism, directionality, exposure, safety and rescue validation.
6. Keep DAVID as secondary in this existing workflow unless a test with matched universe shows meaningful incremental insight. Prefer adding **Functional Annotation Clustering** over simply re-running all the same GO ORA terms.

## 9. Primary methods references

- [Sherman et al., **DAVID 2025 update** (published 2026), *Nucleic Acids Research*](https://doi.org/10.1093/nar/gkag470) — newest DAVID architecture, REST interfaces, viewer and tool changes.
- [Sherman et al., **DAVID 2021 update** (2022), *Nucleic Acids Research*](https://doi.org/10.1093/nar/gkac194) — annotation coverage including disease/drug/gene sources.
- [Huang, Sherman & Lempicki (2009), *Nature Protocols*](https://doi.org/10.1038/nprot.2008.211) — canonical DAVID analysis protocol.
- [Huang et al., **Bioinformatics enrichment tools** (2009), *Nature Protocols*](https://doi.org/10.1038/nprot.2008.211) — same methodological overview of annotated gene lists and chart/clustering.
- [Jiao et al. (2012), DAVID-WS, *Bioinformatics*](https://doi.org/10.1093/bioinformatics/bts251) — batch-programmatic access.
- [DAVID EASE and testing rules](https://davidbioinformatics.nih.gov/documentation/tools/tech-gene-enrichment.html) · [DAVID Kappa/fuzzy clustering](https://davidbioinformatics.nih.gov/documentation/tools/tech-clustering.html).

**Agent TL;DR:** `Clean foreground + justified universe → DAVID ID QC → EASE/BH term chart → Kappa term clusters/gene clusters → merge with g:Profiler without pooling statistics → check Reactome/STRING/UniProt mechanism → independent drug-target evidence → evidence-tiered hypotheses`. DAVID is an *interpretation enhancer*, never a causal inference or drug selection engine.
