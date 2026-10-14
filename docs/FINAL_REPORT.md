# Comprehensive Final Report: Track 2 — Computational Pathway Network Context & Pharmacological Evidence Adjudication in Mosaic Variegated Aneuploidy

## 1. Executive Summary

### The Track 2 Scientific Question
This investigation addresses the core translational genomics objective of Track 2 for the **Rare Disease, Real Kid: MVA Hackathon 2026**:
> *Which biological pathways, causal regulatory mechanisms, and pharmacological repurposing candidates are systematically corroborated by whole-genome sequencing, interactome topology, and signed causal pathway networks for an individual presenting with Mosaic Variegated Aneuploidy (MVA) and related chromosome instability syndromes?*

### Headline Analytical Trajectory
To answer this question without introducing unverified assumptions, an end-to-end multi-phase computational discovery architecture was implemented:
1. **Upstream Ingestion & Gating**: From an initial cohort of **18,252** variant calls annotated across **8,876** genes, rigorous rarity ($\text{gnomAD AF} \le 0.001$, homozygous count $< 5$) and functional impact gating distilled **2,215** high-confidence candidate variants across **1,985** seed genes, flagging **168** multi-hit gene candidates and retaining **6** canonical MVA variant calls (**Phase 0**).
2. **Molecular Directionality Classification**: Variants were categorized by molecular mechanism (Loss-of-Function null, hypomorphic missense, and cryptic splice disruption). Genes harboring compound heterozygous archetypes were classified and gated under mandatory physical phasing requirements (`REQUIRES_PHASE_VERIFICATION`) (**Phase 1**).
3. **Interactome & Pathway Topology**: The 1,985 seed genes were mapped to high-confidence physical interactome networks (STRING Database v12.5, score $\ge 0.700$, zero expansion), resolving **477** high-confidence protein-protein interaction edges across **1,600** components and identifying **84** significantly enriched biological pathway terms (Reactome v88 and g:Profiler, $\text{FDR} < 0.05$). The second-largest interactome module (Component 2, 37 genes, 52 edges) formed a dedicated MVA-focus subnetwork anchoring canonical spindle assembly checkpoint regulators `BUB1B` and `MAD1L1` (**Phase 2**).
4. **Pharmacological Space Retrieval**: An 80-gene multi-tiered query panel was constructed spanning core MVA causal genes, graph-derived interactome neighbors, clinical context modifiers, and benchmark controls. Querying Open Targets Platform v4, ChEMBL v2.0, and DGIdb yielded **2,212** drug-target evidence pairs and **2,006** ranked hypotheses, revealing that **38 of 80** targets possessed zero cataloged pharmacological interactions (**Phase 3**).
5. **Signed Causal & Reaction Adjudication**: Shifting the analytical unit from isolated gene targets to signed regulatory cascades and biochemical reactions, **1,088** qualifying candidate pairs entered a sequential 7-gate adjudication pipeline. Reconciling **39,664** SIGNOR 3.0 edges and **1,746** OmniPath 2024 curated edges (resolving **265** cross-database sign conflicts as `AMBIGUOUS`) and mapping **21** Reactome biochemical reactions identified **442 accepted research hypotheses** (`ACCEPTED_RESEARCH_HYPOTHESIS`) led by directionally concordant ALK-inhibitor regulatory cascades (**Phase 4**).

### Primary Findings
- **Leading Molecular Etiology**: The strongest genetic hypothesis remains a potential compound heterozygous disruption in `BUB1B`, comprising an unambiguous nonsense variant (`p.Leu737*`, CADD 36, ClinVar Pathogenic) paired with an ultra-rare missense variant (`p.Asn1002Lys`, CADD 24.5, AlphaMissense 0.9229). However, because short-read sequencing cannot resolve physical chromosome phasing, this finding remains strictly a research candidate requiring orthogonal long-read or parental segregation confirmation.
- **Leading Pharmacological Repurposing Leads**: The best-supported class of research hypotheses belongs to **Class H1 (Disease-Mechanism Correction)**, led by **ALK inhibitors** (`LORLATINIB` rank 1, `PLB1003`, `ALECTINIB`, `BRIGATINIB`, `CRIZOTINIB`, `ENSARTINIB`; Mechanism Confidence = `0.85`). These candidates are corroborated by 6 curated signed regulatory paths from 3 independent repositories demonstrating that ALK activity acts negatively upstream of checkpoint-surveillance endpoints (e.g., downstream p53 phosphorylation and checkpoint surveillance cascades).
- **Core Target Druggability Reality**: Direct small-molecule stimulation of core structural MVA proteins (`BUB1B`, `CEP57`, `TRIP13`, `CENATAC`, `SLF2`, `SMC5`) is completely absent in public pharmacological databases. This reflects fundamental cell biology: structural kinetochore components, centrosomal anchors, and minor spliceosome factors lack druggable enzymatic pockets, establishing that future disease-modifying strategies must focus on upstream/parallel signaling cascades or genetic/RNA therapeutics.

> **CRITICAL ETHICAL & RESEARCH NOTICE**
> All analyses, workflows, models, and outputs documented in this report are conducted strictly for **computational research and candidate prioritization**. They do **not** constitute medical advice, diagnostic certainty, or clinical treatment recommendations. All pharmacological agents discussed are computational research hypotheses. Every output record generated in this project is explicitly labeled with `CLINICAL_RECOMMENDATION_STATUS = 'PROHIBITED / INVESTIGATIONAL ONLY'`.

---

## 2. Scope, Ethics, and Patient Data Governance

### Research-Only Mandate
Mosaic Variegated Aneuploidy is a devastating, ultra-rare pediatric genetic disorder characterized by constitutional mosaic aneuploidies, severe microcephaly, growth deficiency, intellectual disability, and high susceptibility to childhood embryonal malignancies (most notably rhabdomyosarcoma and Wilms tumor). The work presented here was executed under strict computational research boundaries. No findings herein should be construed as clinical trial indications, off-label treatment proposals, or prognostic certainties.

### Data Privacy and Cryptographic Governance
In accordance with international ethical standards for pediatric genomic research and repository-level safety mandates:
1. **Absolute Exclusion of Patient Data from Version Control**: Raw whole-genome sequencing FASTQ files, aligned BAM/CRAM files, primary VCFs, and patient-level variant manifests remain strictly outside version control, isolated within immutable, `.gitignore`-governed storage paths (`input/`, `data/raw/`, `data/interim/`).
2. **Omission of Identifying Coordinates**: No identifiable genomic coordinates, private patient sequence tags, or individual-level medical identifiers are exposed in public report narratives.
3. **Data-Blind Execution Architecture**: Execution scripts and test suites operate deterministically on synthetic mock data or aggregated, sanitized feature matrices. Synthetic regression tests ensure that software operations never depend on direct inspection or leakage of raw patient records.
4. **Non-Diagnostic Framing**: To prevent therapeutic misinterpretation, every drug-target pair evaluated across Phases 3 and 4 carries explicit computer-readable disclaimer columns: `CLINICAL_RECOMMENDATION_STATUS = 'PROHIBITED / INVESTIGATIONAL ONLY'` and `RESEARCH_ONLY_DISCLAIMER`.

---

## 3. Data Intake & Preliminary Upstream Work (Steps 1–3)

Before entering the Track 2 analytical pipeline, the raw patient whole-genome sequencing dataset underwent three foundational preparation steps:

```text
Raw Patient WGS VCF (GRCh38)
       │
       ▼
[Step 1: Quality Control & Normalization] ──► bcftools norm, DP >= 10, GQ >= 20, FILTER == PASS
       │
       ▼
[Step 2: Consequence Annotation]         ──► Ensembl VEP v116 canonical transcripts + snpEff
       │
       ▼
[Step 3: Multimodal Evidence Enrichment] ──► gnomAD v4.1, ClinVar, CADD, REVEL, AlphaMissense, SpliceAI
       │
       ▼
Consolidated Master Table (variant_evidence_master.csv, 18,252 variants)
```

### Step 1: Primary Quality Control & Normalization
- **Objective**: Eliminate sequencing artifacts, decompose multiallelic genomic loci into distinct biallelic records, left-align indels, and enforce strict genotype quality standards against the GRCh38 human reference genome.
- **Execution Thresholds**: Minimum sequencing read depth $\text{DP} \ge 10$; minimum genotype quality $\text{GQ} \ge 20$; variant filter status $\text{FILTER} == \text{"PASS"}$; multiallelic splitting via `bcftools norm -m -any`.
- **Outcome**: Successfully established unambiguous genomic coordinate keys (`CHROM:POS:REF:ALT`) for 18,252 candidate loci.

### Step 2: Functional Consequence Annotation
- **Objective**: Determine transcript-level consequences against standard Ensembl 116 gene models, prioritizing canonical protein-coding transcripts while verifying consequence annotations orthogonally.
- **Tools**: Ensembl Variant Effect Predictor (VEP v116) utilizing flags `--canonical`, `--hgvs`, `--symbol`, and `--biotype protein_coding`, alongside parallel verification via `snpEff` GRCh38 databases.
- **Outcome**: Standardized assignment of sequence ontology consequence terms (nonsense, frameshift, essential splice, missense, inframe indel, synonymous, intronic).

### Step 3: Multimodal Variant Evidence Enrichment
- **Objective**: Augment every candidate variant with orthogonal predictive layers spanning population genetics, clinical precedent, deep learning structural predictions, and cryptic splicing metrics:
  - *Population Allele Frequencies*: `gnomAD v4.1` joint global allele frequency (`AF`), allele count (`AC`), total alleles (`AN`), and homozygous individual counts (`NHOMALT`).
  - *Clinical Evidence*: `ClinVar` clinical assertions (Pathogenic, Likely Pathogenic, VUS, Conflicting).
  - *Deleteriousness Ensembles*: `CADD v1.7 PHRED` score, `REVEL v1.3` ensemble missense score, and `AlphaMissense` structural constraint scores.
  - *Splicing Disruption*: `SpliceAI` neural network predictions, capturing maximum delta scores (`DS_AG`, `DS_AL`, `DS_DG`, `DS_DL`) and nucleotide distance offsets (`DP_*`).
  - *Regulatory & Phenotypic Alignment*: `AlphaGenome` AVI scores and `Exomiser` phenotypic relevance scores incorporating patient Human Phenotype Ontology (HPO) terms.
- **Output Artifact**: Consolidated candidate evidence master table (`variant_evidence_master.csv`), comprising 18,252 enriched variant records across 8,876 unique genes.

---

## 4. Multi-Phase Analysis Workflow (Phases 0–3)

Track 2 was organized into four sequential, modular phases that progressively compressed raw candidate variation into refined biological networks and actionable pharmacological hypotheses.

```mermaid
flowchart LR
    P0["Phase 0: Intake & Triage<br/>(18,252 -> 2,215 variants)"] --> P1["Phase 1: Mechanism & Phase Gating<br/>(8,876 genes, Archetype Fit)"]
    P1 --> P2["Phase 2: PPI & Pathway Networks<br/>(1,985 seed genes, 477 STRING edges)"]
    P2 --> P3["Phase 3: Pharmacology Engine<br/>(80-gene panel, 2,212 drug pairs)"]
    P3 --> P4["Phase 4: Signed Causal Adjudication<br/>(1,088 entry pairs, 442 accepted)"]
```

### Phase 0: Candidate Intake, Gating & Consolidation

$$\text{Upstream Step 3 } (\texttt{variant\_evidence\_master.csv}) \longrightarrow \text{Phase 0 Engine} \longrightarrow \texttt{00\_phase0\_consolidated\_candidates.csv}$$

- **Task & Purpose**: Ingest the 18,252 candidate records from preliminary enrichment, execute transparent triage filter gates, score multi-evidence concordance, flag multi-hit genes, and categorize loci into canonical MVA and Spindle Assembly Checkpoint (SAC) gene sets.
- **Gating Rules & Thresholds**:
  1. *Rarity Gate (`GATE_RARITY`)*: $\text{gnomAD AF} \le 0.001$ (or population frequency missing/unobserved) AND $\text{gnomAD NHOMALT} < 5$. Alleles exceeding 0.1% global frequency or observed homozygously in $\ge 5$ reference individuals fail.
  2. *Impact Gate (`GATE_IMPACT`)*: Variant must satisfy at least one of:
     - High-impact consequence (`stop_gained`, `frameshift_variant`, `splice_acceptor_variant`, `splice_donor_variant`, `start_lost`, `stop_lost`);
     - ClinVar Pathogenic or Likely Pathogenic status;
     - $\text{CADD PHRED} \ge 20.0$ (top 1% most deleterious in the human genome);
     - $\text{SpliceAI MAX\_DS} \ge 0.20$ (cryptic splice alteration threshold);
     - $\text{REVEL} \ge 0.50$ OR $\text{AlphaMissense} \ge 0.50$ (`likely_pathogenic`).
  3. *Compound Overall Gate (`GATE_OVERALL_PASS`)*: Variant must pass coordinate validation AND rarity gating (or unobserved frequency in an MVA gene) AND impact gating.
  4. *Biological Categorization*: Mapped variants to `CANONICAL_MVA` (`BUB1B`, `CEP57`, `TRIP13`), `REPORTED_MVA_CA` (`MAD1L1`, `CENATAC`), `OVERLAPPING_ATELIS` (`SLF2`, `SMC5`), `SAC_CORE` (`AURKB`, `BUB1`, `BUB3`, `CDC20`, `CENPE`, `MAD2L1`, `PLK1`, `TTK`), and `GENOME_WIDE`.
- **Results**:
  - Total variants evaluated: **18,252**
  - Total variants passing all filter gates: **2,215** (~12.1% pass rate)
  - Multi-hit candidate genes ($\ge 2$ passing variants): **168**
  - Canonical MVA variants passing gates: **6**
- **Conclusions**:
  - *Supported by Evidence*: Successfully distilled 18,252 calls down to 2,215 prioritized research variants, isolating 168 candidate multi-hit genes and preserving canonical MVA hits.
  - *Limitations*: Severe somatic mosaicism below standard germline caller sensitivity ($< 15\%$ variant allele fraction) cannot be evaluated without deep targeted sequencing.

### Phase 1: Molecular Mechanism Classification & Biallelic Archetype Gating

$$\texttt{00\_phase0\_consolidated\_candidates.csv} \longrightarrow \text{Phase 1 Engine} \longrightarrow \texttt{01\_phase1\_variant\_mechanism\_matrix.csv} \ \& \ \texttt{02\_phase1\_gene\_mechanism_matrix.csv}$$

- **Task & Purpose**: Classify the precise molecular consequence of each variant into functional mechanism classes (Loss-of-Function null, hypomorphic missense, cryptic splice disruption), evaluate directional effect on protein abundance/activity, and identify recessive compound heterozygous archetypes.
- **Classification Rules & Logic**:
  - *LoF Null*: Stop-gained, frameshift, or essential splice donor/acceptor variants resulting in nonsense-mediated decay or truncated, non-functional polypeptides.
  - *Cryptic Splice Disruption*: Intronic or synonymous variants exhibiting $\text{SpliceAI MAX\_DS} \ge 0.20$.
  - *Hypomorphic Missense*: Missense alterations exhibiting high computational deleteriousness ($\text{REVEL} \ge 0.50$ or $\text{AlphaMissense} \ge 0.50$ or $\text{CADD} \ge 20.0$).
  - *Compound-Het Archetype Fit*: Evaluates whether a gene possesses the classic Mendelian recessive MVA architecture: a severe Loss-of-Function null allele on one parental chromosome complemented by a hypomorphic missense or splice-modulating allele on the other.
  - *Phase Verification Guardrail*: All multi-hit genes are tagged with `PHASE_GATE_STATUS = 'REQUIRES_PHASE_VERIFICATION'`.
- **Results**:
  - Variant mechanism matrix generated across all **18,252** variants (`01_phase1_variant_mechanism_matrix.csv`).
  - Gene-level mechanism matrix generated across **8,876** genes (`02_phase1_gene_mechanism_matrix.csv`).
  - `BUB1B` was flagged as a prime compound heterozygous archetype fit: harboring a nonsense mutation (`c.2210T>G`, `p.Leu737*`, LoF Null) and a rare damaging missense mutation (`c.3006T>G`, `p.Asn1002Lys`, Hypomorphic Missense).
- **Conclusions**:
  - *Supported by Evidence*: Computational evidence confirms `BUB1B` perfectly matches the published biallelic MVA1 disease archetype (Nature Genetics, 2004; DOI: 10.1038/ng1449).
  - *Limitations / Phase Guardrail*: Short-read WGS cannot determine physical phase. Whether these two variants reside in *trans* (compound heterozygous, pathogenic) or in *cis* (same chromosome, carrier state) remains completely unproven without long-read sequencing or parental segregation analysis.

### Phase 2: Protein-Protein Interaction & Pathway Network Context Engine

$$\texttt{02\_phase1\_gene\_mechanism_matrix.csv} \longrightarrow \text{Phase 2 Engine} \longrightarrow \texttt{03\_phase2_seed_genes.csv} \dots \texttt{13\_phase2_pathway_enrichment.html}$$

- **Task & Purpose**: Map the 1,985 seed genes harboring passing candidate variants into physical interactome space using the STRING Database v12.5 and project them onto biological pathways via Reactome and g:Profiler.
- **Topological & API Parameters**:
  - *Seed Definition*: 1,985 genes possessing $\ge 1$ passing candidate variant from Phase 0/1; 100% mapped to standard HGNC nomenclature via g:Profiler `g:Convert`.
  - *STRING DB v12.5 Query*: High-confidence score threshold $\ge 0.700$ ($\ge 700$); zero network expansion (no unseeded interactors added).
  - *Pathway Overrepresentation*: Reactome ContentService v88 and g:Profiler `g:GOSt` across GO:BP, KEGG, Reactome, and WikiPathways, filtered at Benjamini-Hochberg $\text{FDR} < 0.05$.
- **Results**:
  - Total high-confidence PPI edges: **477** (average node degree: 0.48).
  - Network components: **1,600** components (1,503 isolated singletons).
  - Component 1 (largest): **111 genes**, enriched for generic macromolecular complexes.
  - Component 2 (MVA-Focus Component, Rank 2): **37 genes, 52 edges**, densely enriched for spindle assembly checkpoint and mitotic regulators, anchoring `BUB1B`, `MAD1L1`, and `MAD2L2`.
  - Pathway Enrichment: **84 significant pathway terms** ($\text{FDR} < 0.05$), highlighting Spindle Assembly Checkpoint (R-HSA-69620), Mitotic Prometaphase (R-HSA-68877), Resolution of Sister Chromatid Cohesion (R-HSA-2500257), and Kinetochore Organization.
- **Conclusions**:
  - *Supported by Evidence*: Network analysis confirms that candidate genomic variation converges topologically upon mitotic checkpoint machinery, specifically isolating a coherent 37-gene module anchoring canonical MVA drivers.
  - *Limitations*: 75.6% of seed genes (1,503) are isolated singletons at STRING score $\ge 0.700$, indicating that severe filtering without network expansion restricts interactome connectivity.

### Phase 3: Query Gene Panel & Multi-Source Drug Repurposing Evidence Engine

$$\texttt{07\_phase2_network_gene_context.csv} \longrightarrow \text{Phase 3 Panel & Drug Engine} \longrightarrow \texttt{14\_phase3_drug_evidence_matrix.csv} \dots \texttt{19_phase3_query_panel_citations.csv}$$

- **Task & Purpose**: Construct a rigorous, multi-tiered 80-gene pharmacological inquiry panel and query public chemogenomic and drug databases across all 80 targets, evaluating directional alignment and contraindication hazards.
- **80-Gene Query Panel Composition**:
  - `CORE_MVA_CAUSAL` (7 genes): `BUB1B`, `MAD1L1`, `CENATAC`, `SLF2`, `SMC5`, plus injected canonical genes `CEP57` and `TRIP13` (absent from Phase 1 passing variants, injected to ensure comprehensive disease coverage).
  - `GRAPH_DERIVED_NEIGHBOR` (31 genes): High-confidence interactome partners from Component 2 (e.g., `ALK`, `AURKA`, `PLK1`, `CDC20`, `MAD2L1`, `CENPE`, `KNL1`).
  - `CONTEXT_MODIFIER` (29 genes): Modifiers and surveillance hubs, including curated literature controls `TP53`, `ATM`, `NF1`, and `MAD2L2`.
  - `DRUG_TARGET_CONTROL` (13 genes): 9 positive control targets with established pharmacological profiles (`AURKA`, `PLK1`, `BUB1`, `PLK4`, `ATR`, `CHEK1`, `WEE1`, `CDC20`, `MAD2L1`) and 4 negative control olfactory receptors (`OR7G3`, `OR10AD1`, `OR10D3`, `OR10X1`).
- **Pharmacological Search Parameters**:
  - Queried Open Targets Platform v4 GraphQL, ChEMBL v2.0 REST, and DGIdb GraphQL APIs.
  - Triage into 4 Safety Tiers: Tier 1 (Concordant Hypothesis, 5 pairs), Tier 2 (Off-Target Mitotic Compensation, 99 pairs), Tier 3 (Uncertain Directionality, 1,891 pairs), Tier 4 (Contraindicated Hazard, 11 pairs).
- **Results**:
  - Total retrieved drug-target evidence pairs: **2,212** (DGIdb: 1,898; Open Targets: 161; ChEMBL: 153).
  - Total unique drug-target hypotheses ranked: **2,006**.
  - Mechanistic hazards flagged: **15** records (`CONTRADICTION_FLAG = True`).
  - Target Coverage Reality: **38 of 80 genes returned ZERO drug evidence** (`RESEARCH_STATUS = 'NO_DRUG_EVIDENCE'`), including core MVA genes `CEP57`, `TRIP13`, `SLF2`, `SMC5`, and `CENATAC`.
  - Benchmark Controls: Negative controls passed with 4/4 returning zero hits (`PASS_ZERO_HITS`); positive controls yielded 7/9 hits (with `CDC20` and `MAD2L1` yielding zero hits due to public API catalog limitations, logged as `WARNING_ZERO_HITS`).
- **Conclusions**:
  - *Supported by Evidence*: Systematic retrieval established the first comprehensive pharmacological landscape for MVA interactome nodes, identifying 5 concordant Tier 1 leads and rejecting 11 contraindicated hazards.
  - *Key Biological Insight & Phase 4 Motivation*: Direct small-molecule druggability for core MVA causal genes is virtually absent in public chemistry. Structural scaffold proteins cannot be "activated" by standard inhibitors, proving that any therapeutic hypothesis must operate at the **signed pathway and reaction cascade level**.

---

## 5. Phase 4 Methodology: Signed Causal Pathway & Reaction-Level Mechanism Adjudication

### Pre-Registration Principles & Causal Framework
To resolve the direct druggability bottleneck of core MVA proteins, Phase 4 transitions from static target-level drug matching to **signed causal pathway and biochemical reaction network adjudication**. The mathematical and operational architecture was pre-registered prior to live execution:
1. **Directional Sign Conventions**: Every biological interaction and pharmacological action is mapped to a discrete sign:
   - Activation / Agonism / Upregulation: $+1$
   - Inhibition / Antagonism / Degradation / Downregulation: $-1$
   - Ambiguous / Conflicting / Unknown: $0$
2. **Causal Horizon (Max Path Length $\le 2$)**: Causal path traversal is restricted to a maximum path length of 2 steps ($\text{Target} \rightarrow \text{Intermediate} \rightarrow \text{Endpoint}$) to prevent biological signal degradation and spurious cross-pathway drift.
3. **Hierarchy of Evidence**: Curated primary pathway databases take absolute precedence over text-mined or inferred interactions:
   $$\text{Reactome Biochemical Reactions} \equiv \text{SIGNOR 3.0 Curated} \equiv \text{OmniPath Curated} \gg \text{OmniPath Inferred}$$
4. **Frozen Hazard Weights & Thresholds**:
   - Baseline Pediatric Hazard: Tier 4 Contraindicated = 85; Tier 3 Uncertain = 45; Tier 2 Mitotic Compensation = 35; Tier 1 Concordant = 20.
   - Liability Additions: Black-box warning / fatality risk = $+20$; pediatric growth arrest / bone marrow suppression = $+15$; secondary malignancy risk = $+20$; TP53 abrogation in constitutional aneuploidy = 100.
   - Pediatric Hold Threshold: Hazard Score $> 60$ triggers a mandatory pediatric safety hold.

### Knowledge Repositories & Operational Parameters
Phase 4 integrates live, curated pathway and interaction knowledgebases:

| Knowledgebase | Version / Release | Type of Evidence | Operational Role | Status in Phase 4 |
|---|---|---|---|---|
| **SIGNOR** | 3.0 (Homo sapiens 9606) | Curated causal activity flow | Signed directed regulatory edges (39,664 edges) | Live Ingestion |
| **OmniPath** | 2024 Curated Release | Curated signaling meta-network | Multi-source signed causal edges (1,746 edges) | Live Ingestion |
| **Reactome** | ContentService v88 | Curated biochemical reactions | Reaction transitions, catalysts, regulators (21 reactions) | Live Adjudication |
| **Europe PMC** | REST API v2026 | Literature contradiction screening | Text-mining contraindications & adverse reports | Live Screening |
| **NCBI E-utilities / Crossref** | October 2026 Release | Bibliographic verification | Verification of primary discovery DOIs and PMIDs | 100% Live Verified |

### Candidate Pair Entry Criteria (E1–E5)
To ensure computational rigor, exactly 1,088 unique candidate pairs entered Phase 4 satisfying pre-registered criteria E1–E5:
- **E1 (Tier 1 Concordant Leads)**: All 5 concordant hypotheses from `15_phase3_safety_ranking.csv`.
- **E2 (Tier 2 Multi-Source Candidates)**: 25 Tier 2 mitotic compensation hypotheses corroborated by $\ge 2$ independent databases.
- **E3 (Tier 4 Negative Controls)**: 11 contraindicated hazard hypotheses serving as benchmark negative controls.
- **E4 (Malignancy-Indication Targets)**: Pairs targeting MVA-associated pediatric malignancy genes (`ALK`, `KIT`, `PARP1`, `TP53`).
- **E5 (MVA-Focus Interactome Targets)**: Pairs targeting Component 2 interactome genes (`MAD1L1`, `BUB1B`, `ATM`, `NF1`, `HIF1A`), plus 6 core MVA causal genes with zero pharmacology cataloged emitted as explicit `EVIDENCE_GAP` nodes (`BUB1B`, `CEP57`, `TRIP13`, `CENATAC`, `SLF2`, `SMC5`).

### Signed Graph Traversal Mathematics & Causal Net Flux
For any directed causal path $P = (e_1, e_2, \dots, e_k)$ connecting a drug target $T$ to a downstream cellular surveillance endpoint $E$, the net path sign $S(P)$ is the product of constituent edge signs:
$$S(P) = \prod_{i=1}^k s(e_i), \quad \text{where } s(e_i) \in \{+1, -1, 0\}$$
When modulated by a pharmacological compound with action sign $s_{\text{drug}} \in \{+1, -1\}$:
$$F_{\text{net}} = s_{\text{drug}} \times S(P)$$
- **Concordant Flux ($F_{\text{net}} > 0$)**: The compound's action moves the downstream endpoint toward physiological checkpoint restoration or stabilization.
- **Discordant Hazard ($F_{\text{net}} < 0$)**: The compound's action directly exacerbates checkpoint failure, accelerates chromosome missegregation, or abrogates genomic guardians. Flagged as `SIGNED_PATH_HAZARD = True`.
- **Zero / Ambiguous Flux ($F_{\text{net}} = 0$)**: Occurs when any constituent edge is uncharacterized, ambiguous, or contradictory.

### Cross-Database Sign Conflict Resolution Policy
Biological knowledgebases occasionally report conflicting regulatory directions for identical protein pairs (e.g., SIGNOR annotating activation while OmniPath annotates inhibition). Phase 4 enforces a strict, pre-registered **anti-averaging policy**:
- Conflicting curated edges are reconciled as `effect = 'AMBIGUOUS'` with `sign = 0` and `conflict_flag = True`.
- Conflicting edges are **never averaged** into fractional values (e.g., $+0.5$).
- Consequently, paths traversing ambiguous edges evaluate to $F_{\text{net}} = 0$, completely halting causal signal propagation through contested biology.

### Sequential Multi-Gate Architecture (Gates G1–G7)
Candidate pairs traverse seven sequential evaluation gates:
1. **Gate 1 (G1: Identifier Hygiene)**: Validates target symbol, UniProt accession, Ensembl ID, and ChEMBL compound accession. Rejects unmapped or ambiguous entities into quarantine.
2. **Gate 2 (G2: Signed-Tree & Reaction Reachability)**: Searches for directed signed paths ($L \le 2$) and biochemical reaction participations reaching defined endpoints (`BUB1B`, `MAD1L1`, `CDC20`, `AURKA`, `PLK1`, `TP53`, `ATM`). Unreachable pairs are tagged with `EVIDENCE_GAP`.
3. **Gate 3 (G3: Directionality Concordance)**: Evaluates net flux $F_{\text{net}}$ and checks for lethal mechanistic contradictions (e.g., inhibition of tumor suppressor `TP53` or inhibition of genes whose therapeutic requirement is activation). Detects `SIGNED_PATH_HAZARD`.
4. **Gate 4 (G4: Multi-Hypothesis Class Assignment)**: Classifies candidates into H1 (Disease Mechanism), H2 (Cancer Complication), H3 (Supportive Care), or H0 (Unclassified). Pairs may receive multiple valid classifications.
5. **Gate 5 (G5: Pediatric & Regulatory Hazard Gate)**: Computes Pediatric Hazard Score (0–100). Scores $> 60$ trigger a mandatory pediatric hold.
6. **Gate 6 (G6: Literature Validation & Contradiction Screening)**: Validates primary discovery literature via Crossref and NCBI E-utilities, and screens Europe PMC for published contraindications.
7. **Gate 7 (G7: Final Verdict Assignment)**: Issues adjudicated verdicts according to rigorous corroboration rules.

### Multi-Source Corroboration & Verdict Assignment Rules
Verdicts are assigned deterministically under pre-registered standards:
- **`ACCEPTED_RESEARCH_HYPOTHESIS`**: Requires $\ge 1$ curated signed path, corroboration across $\ge 2$ independent curated repositories, verified primary literature, zero literature contradictions, and a Pediatric Hazard Score $\le 60$.
- **`KEEP_WITH_CAVEAT`**: Mechanistically concordant candidates that rely on a single curated repository, have literature pending full verification, present elevated pediatric toxicity liabilities (Hazard Score $> 60$), or represent oncologic agents evaluated in a constitutional aneuploidy background.
- **`EVIDENCE_INSUFFICIENT`**: Assigned to candidate nodes exhibiting reachability gaps (`EVIDENCE_GAP`).
- **`MECHANISTIC_HAZARD_REJECTED`**: Assigned to compounds whose action directly opposes therapeutic requirements or abrogates essential genomic guardians (`SIGNED_PATH_HAZARD = True`).
- **`LITERATURE_CONTRADICTED`**: Assigned if peer-reviewed literature documents direct clinical contraindication or exacerbation in chromosome instability models.
- **`NEGATIVE_NO_SAFE_INTERVENTION`**: Assigned to designated negative controls.

---

## 6. Phase 4 Execution Results: Multi-Gate Adjudication & Best-Supported Options

### Sequential Gate Traversal & Attrition Accounting
Phase 4 executed live (`offline_mode = False`), producing complete traversal accounting across all seven gates:

| Adjudication Gate | Metric Evaluated | Pass Count | Hold / Hazard / Gap Count | Operational Note |
|---|---|---:|---:|---|
| **Gate 1 (G1)** | Identifier Hygiene | **1,088** | 0 failed | 100% valid HGNC and ChEMBL entity mappings. |
| **Gate 2 (G2)** | Reachability & Reactions | **1,058** | 30 evidence gaps | 1,058 reachable to SAC/mitotic endpoints; 30 evidence gaps. |
| **Gate 3 (G3)** | Directional Concordance | **1,077** | 11 mechanistic hazards | 1,077 concordant; 11 flagged for lethal hazard. |
| **Gate 4 (G4)** | Hypothesis Classification | **1,136** | 0 unclassified errors | 1,136 classified instances emitted across 1,088 pairs. |
| **Gate 5 (G5)** | Pediatric Hazard Scoring | **1,077** | 11 pediatric holds | 1,077 score $\le 60$; 11 exceed safety threshold (score $> 60$). |
| **Gate 6 (G6)** | Literature Validation | **942** | 146 pending / 0 contradicted | 942 verified citations; zero published contradictions. |
| **Gate 7 (G7)** | Final Verdict Assignment | **1,136** | — | Final adjudicated verdict distribution established. |

### Final Adjudicated Verdict Distribution
Across the 1,136 classified hypothesis instances, verdicts distributed as follows:

| Adjudicated Verdict | Total Rows | Percentage | Description & Operational Handling |
|---|---:|---:|---|
| **`ACCEPTED_RESEARCH_HYPOTHESIS`** | **442** | 38.9% | Top-tier hypotheses: multi-source curated signed paths, verified literature, acceptable pediatric risk profile. |
| **`KEEP_WITH_CAVEAT`** | **653** | 57.5% | Concordant candidates subject to single-source caveats, elevated pediatric toxicity risks, or oncologic context. |
| **`EVIDENCE_INSUFFICIENT`** | **30** | 2.6% | Target nodes with zero resolvable signed regulatory edges or biochemical reactions in public databases. |
| **`MECHANISTIC_HAZARD_REJECTED`** | **11** | 1.0% | Compounds rejected due to severe directional contradiction or abrogation of genomic guardians (e.g. TP53). |
| **`LITERATURE_CONTRADICTED`** | **0** | 0.0% | Zero candidates contradicted by specific published exacerbation reports. |
| **`NEGATIVE_NO_SAFE_INTERVENTION`**| **0** | 0.0% | Benchmarked within hazard rejection set. |
| **Total Evaluated** | **1,136** | **100.0%** | Comprehensive multi-evidence adjudicated cohort. |

### Hypothesis Class Stratification (H1, H2, H3, H0)
The 1,136 classified instances stratified across the three pre-registered hypothesis classes:

| Hypothesis Class | Total Instances | Description | Status & Biological Reality |
|---|---:|---|---|
| **Class H1: Disease-Mechanism Correction** | **59** | Upstream modulation stabilizing SAC surveillance or KT-MT attachment | Highly enriched for ALK inhibitor cascades converging on p53 surveillance. |
| **Class H2: Cancer-Complication Management** | **799** | Modulating pediatric malignancy indications (RMS, Wilms, leukemia) | Mitotic kinase inhibitors with strict pediatric toxicity caveats. |
| **Class H3: Supportive-Care Symptom Management** | **0** | Secondary phenotypes (nephrocalcinosis, tubular ion transport) | **Coverage Finding**: Zero eligible records in panel / databases. |
| **Class H0: Unclassified / General Baseline** | **278** | Interactors lacking direct signed path to designated endpoints | Documented baseline records. |

### Best-Supported Candidates: Class H1 (Disease-Mechanism Correction) — ALK-Axis Leads
The highest-ranking concordant hypotheses identified in Phase 4 belong to **Class H1**, targeting the anaplastic lymphoma kinase (**ALK**) signaling axis:

| Target Symbol | Drug Candidate | ChEMBL ID | Hypothesis Class | Mechanism Confidence | Pediatric Hazard | Verdict Status | Evidence Nature & Path Rationale |
|---|---|---|---|---:|---:|---|---|
| **ALK** | **LORLATINIB** | CHEMBL3545107 | H1_DISEASE_MECHANISM | 0.85 | 20 | `ACCEPTED_RESEARCH_HYPOTHESIS` | 6 curated signed paths from 3 independent sources; ALK activity acts negatively upstream of p53 surveillance. |
| **ALK** | **PLB1003** | CHEMBL4298135 | H1_DISEASE_MECHANISM | 0.85 | 20 | `ACCEPTED_RESEARCH_HYPOTHESIS` | Multi-source curated signed path; directionally concordant with checkpoint stabilization. |
| **ALK** | **ALECTINIB** | CHEMBL2386561 | H1_DISEASE_MECHANISM | 0.85 | 20 | `ACCEPTED_RESEARCH_HYPOTHESIS` | 6 curated signed paths; approved small-molecule inhibitor; directionally concordant. |
| **ALK** | **BRIGATINIB** | CHEMBL3301610 | H1_DISEASE_MECHANISM | 0.85 | 20 | `ACCEPTED_RESEARCH_HYPOTHESIS` | Multi-source curated signed path; directionally concordant with checkpoint stabilization. |
| **ALK** | **CRIZOTINIB** | CHEMBL584143 | H1_DISEASE_MECHANISM | 0.85 | 20 | `ACCEPTED_RESEARCH_HYPOTHESIS` | Curated signed path; established pediatric oncology dosing data; directionally concordant. |
| **ALK** | **ENSARTINIB** | CHEMBL3989930 | H1_DISEASE_MECHANISM | 0.85 | 20 | `ACCEPTED_RESEARCH_HYPOTHESIS` | Curated signed path; second-generation inhibitor; directionally concordant. |

#### Biological Rationale for ALK Inhibition
In SIGNOR 3.0 and OmniPath 2024 curated signaling networks, ALK receptor tyrosine kinase activity transmits negative regulatory inputs upstream of genomic surveillance endpoints, specifically repressing p53-dependent surveillance and checkpoint arrest cascades. Small-molecule ALK inhibitors exert a negative action sign ($s_{\text{drug}} = -1$) upon an inhibitory regulatory path ($S(P) = -1$). By double-negative logic ($F_{\text{net}} = (-1) \times (-1) = +1$), ALK inhibition results in a net positive, concordant impulse that relieves repression on downstream checkpoint surveillance, thereby promoting cell cycle arrest in the presence of aberrant mitotic attachments.

> **CRITICAL EVIDENTIARY NUANCE: ALK IS A GRAPH NEIGHBOR, NOT AN MVA GENE**
> It is essential to recognize the exact scientific nature of this finding:
> 1. **Verified Literature Supports the Endpoint, Not an ALK-MVA Link**: The verified primary literature associated with these paths supports the downstream **p53 phosphorylation and checkpoint surveillance endpoint mechanism** (BMC Cell Biology, 2001, DOI: 10.1186/1471-2121-2-12; Oncogene, 2001, DOI: 10.1038/sj.onc.1204742; Journal of Cellular Biochemistry, 2003, DOI: 10.1002/jcb.10569). There is **no published literature demonstrating a direct clinical efficacy link between ALK inhibitors and Mosaic Variegated Aneuploidy**.
> 2. **Network Origin**: `ALK` entered the query panel solely as a **graph-derived interactome neighbor (`GRAPH_DERIVED_NEIGHBOR`)** from Component 2 of the Phase 2 protein-protein interaction network.
> 3. **Absence of Patient Mutation**: The patient carries **no pathogenic or damaging variants in `ALK`**. `ALK` is not a causal MVA gene.
> 4. **Translational Status**: ALK inhibition represents a computationally generated, network-level repurposing hypothesis. It must be investigated strictly in cell-line models of aneuploidy before any translational speculation.

### Best-Supported Candidates: Class H2 (Cancer-Complication Management) — Mitotic Modulators & Pediatric Caveats
Children with MVA face an exceptionally high risk of early-onset embryonal malignancies, most notably rhabdomyosarcoma and Wilms tumor. Class H2 evaluated candidates targeting these oncologic vulnerabilities:

| Target Symbol | Drug Candidate | ChEMBL ID | Target Role / Indication | Pediatric Hazard | Verdict Status | Operational Caveat & Pediatric Liability |
|---|---|---|---|---:|---|---|
| **ALK** | **ALECTINIB** | CHEMBL2386561 | Neuroblastoma / Sarcoma | 35 | `KEEP_WITH_CAVEAT` | Oncologic context in constitutional aneuploidy; requires tight pediatric surveillance. |
| **AURKA** | **ALISERTIB** | CHEMBL2029707 | Mitotic Aurora Kinase A | 50 | `KEEP_WITH_CAVEAT` | **Severe Pediatric Toxicity**: Mitotic arrest induces myelosuppression and growth plate damage. |
| **PLK1** | **VOLASERTIB** | CHEMBL1289656 | Polo-like Kinase 1 | 50 | `KEEP_WITH_CAVEAT` | **Severe Pediatric Toxicity**: Bone marrow aplasia and intestinal epithelial cytotoxicity. |
| **PARP1** | **OLAPARIB** | CHEMBL514114 | DNA Damage Repair / Synthetic Lethality | 35 | `KEEP_WITH_CAVEAT` | Secondary malignancy risk in constitutional chromosome instability syndromes. |

#### Severe Pediatric Risk Caveats in Class H2
Mitotic kinase inhibitors (e.g., `AURKA` and `PLK1` inhibitors, supported biologically by discovery literature: Journal of Cell Biology, 1995, DOI: 10.1083/jcb.129.6.1617; EMBO Journal, 1998, DOI: 10.1093/emboj/17.11.3052) are designed to exhaust cancer cells through catastrophic mitotic arrest. In an individual with constitutional chromosome instability and mosaic aneuploidy across all somatic tissues:
- Normal, rapidly dividing pediatric tissues (hematopoietic progenitors, epiphyseal growth plates, neurogenic niches) rely on intact residual mitotic machinery.
- Exposing constitutional aneuploid tissues to potent mitotic inhibitors risks triggering catastrophic bone marrow failure, severe growth arrest, and secondary leukemogenesis.
- Therefore, all Class H2 mitotic compensation candidates are adjudicated under **`KEEP_WITH_CAVEAT`** and assigned elevated Pediatric Hazard Scores.

### Class H3 Coverage Finding: Supportive-Care & Symptom-Management Gaps
- **Execution Metric**: Exactly **0 candidate records** met entry criteria for Class H3 (`H3_SUPPORTIVE_CARE_SYMPTOM_MANAGEMENT`).
- **Nature of the Finding**: This result is a **coverage finding**, **not** a computational or software defect.
- **Underlying Biological Cause**: Secondary phenotypic features of MVA (such as nephrocalcinosis, medullary sponge kidney, or tubular ion transport defects) are driven by specialized transport proteins (e.g., `CLDN16`, `CLDN19`, `CASR`, `SLC12A1`). These targets were not members of the Phase 2 MVA mitotic interactome (Component 2), nor are they currently indexed in public pharmacology databases with causal links to mitotic chromosome segregation. This coverage finding explicitly identifies an unmet need for dedicated supportive-care pharmacology mapping in rare genetic aneuploidy syndromes.

### Negative-Control Benchmark Sanity Check: TP53 Abrogation Rejection
To ensure that the adjudication engine reliably eliminates dangerous mechanistic contradictions, 11 drug-target pairs targeting the tumor suppressor `TP53` with inhibitory mechanisms of action were introduced as negative benchmark controls (Criterion E3):
- **Adjudication Result**: **11 of 11 TP53 inhibitor rows were decisively rejected** with the verdict **`MECHANISTIC_HAZARD_REJECTED`**.
- **Scoring & Flagging**: Every row evaluated with `SIGNED_PATH_HAZARD = True`, triggered adverse event flag `TP53_ABROGATION_LETHAL_HAZARD`, and received a catastrophic Pediatric Hazard Score of **85 to 100**.
- **Biological Validation**: In a patient with constitutional aneuploidy, pharmacological abrogation of p53 surveillance would remove the primary cell cycle checkpoint preventing aneuploid cells from proliferating, drastically accelerating malignant transformation. The immediate, deterministic rejection of all 11 TP53 inhibitor rows demonstrates that the multi-gate adjudication engine possesses absolute fidelity against lethal mechanistic contradictions.

### Literature Adjudication & Citation Verification Audit
- **Verification Methodology**: Every cited publication in the pipeline registry was subjected to live bibliographic resolution via NCBI E-utilities (PubMed API) and Crossref REST API.
- **Audit Outcomes**:
  - Pre-verified Discovery Registry: **38** peer-reviewed discovery publications.
  - Live Adjudication Traversal (Gate 6): **942** verified citations, **146** pending verification, and **0** published contradictions.
  - **DOI Live-Resolution Audit**: All 14 primary references cited in this report were verified live on **2026-10-11** via Crossref and NCBI E-utilities, confirming 100% active DOI resolution and accurate bibliographic metadata.

---

## 7. Integrated Synthesis & Biological Discussion

```text
               Patient Genomic Profile
         BUB1B: c.2210T>G (p.Leu737*, LoF Null)
         BUB1B: c.3006T>G (p.Asn1002Lys, Missense)
                          │
                          ▼ (Recessive Architecture)
     Severe Impairment of Spindle Assembly Checkpoint
                          │
                          ▼
            Persistent Chromosome Missegregation
              & Constitutional Mosaic Aneuploidy
                          │
         ┌────────────────┴────────────────┐
         ▼                                 ▼
   Primary Defect                   Network Consequences
  (BUB1B / SAC)                   (Component 2 Interactome)
  - Zero direct small-molecule    - Upstream ALK / RTK signaling
    drugs in public chemistry       suppresses p53 surveillance
  - Druggability gap = expected   - ALK inhibitors: concordant
    structural biology finding      p53 restoration (H1 Lead)
                                  - Mitotic inhibitors (AURKA/PLK1):
                                    severe pediatric toxicity (H2 Caveat)
```

### Convergence on BUB1B Compound Heterozygosity as Primary Genomic Hypothesis
Across the entire genomic profile, the candidate that demonstrates the highest convergence of biological mechanism, rarity, consequence severity, and clinical precedent is **`BUB1B`**:
1. `c.2210T>G` (`p.Leu737*`): A definitive nonsense mutation creating a premature termination codon upstream of the BUB1B kinase domain. It is classified as Pathogenic in ClinVar, carries a CADD PHRED score of 36, and has a gnomAD population frequency of $9.83 \times 10^{-5}$ with zero observed homozygotes.
2. `c.3006T>G` (`p.Asn1002Lys`): An ultra-rare missense mutation ($\text{gnomAD AF} \approx 8.47 \times 10^{-7}$) altering an evolutionarily conserved residue in the C-terminal kinase domain, supported by high structural constraint scores (AlphaMissense 0.9229, PolyPhen 0.997, CADD 24.5).

This pairing perfectly matches the established Mendelian recessive archetype of Mosaic Variegated Aneuploidy Syndrome 1 (MVA1, OMIM #257300), first described by Hanks et al. (Nature Genetics, 2004; DOI: 10.1038/ng1449). In this archetype, complete biallelic loss of BUB1B is embryonically lethal; affected children survive only by inheriting one severe null allele in combination with a partially functional hypomorphic missense or splice-altering allele that preserves minimal essential mitotic function. 

**The Phasing Caveat**: Because the input sequencing consists of standard short-read WGS, physical chromosome phasing cannot be resolved. These two variants could lie in *trans* (on opposite parental chromosomes, creating compound heterozygosity and causing MVA1) or in *cis* (on the same chromosome, representing carrier status on a single haplotype). Proving compound heterozygosity requires parental sequencing or long-read phasing.

### Network Context & Mechanistic Rationale for ALK Signed-Path Repurposing Leads
When evaluating how small molecules might interface with this checkpoint-deficient system, Phase 4 demonstrated that **ALK inhibitors** represent the most robust computational hypothesis (Class H1):
- `ALK` resides in Component 2 of the Phase 2 interactome, physically and functionally linked to mitotic checkpoint machinery.
- Analysis of SIGNOR 3.0 and OmniPath 2024 networks uncovered 6 curated signed regulatory paths wherein ALK activity suppresses p53 checkpoint surveillance endpoints.
- Relieving this suppression with approved, brain-penetrant ALK inhibitors (`LORLATINIB`, `ALECTINIB`) provides a directionally concordant impulse to restore cell cycle arrest in aneuploid cells.
- However, as emphasized in Section 6, this hypothesis is based on network topology and endpoint biochemistry (DOI: 10.1186/1471-2121-2-12, DOI: 10.1038/sj.onc.1204742, DOI: 10.1002/jcb.10569). It represents a novel, unvalidated research direction, not an established therapy.

### Sparse Direct Druggability of Structural Kinetochore and Centrosomal Checkpoint Genes
A central finding of this investigation is that **38 of the 80 evaluated panel genes returned zero pharmacological evidence** in public drug databases:
- Core causal MVA proteins (`CEP57`, `TRIP13`, `SLF2`, `SMC5`, `CENATAC`) and structural kinetochore anchors lack small-molecule pharmacological modulators.
- This is an **expected biological finding**, not a computational defect. Small molecules typically function by binding defined hydrophobic pockets to inhibit enzymatic activity. Core MVA proteins operate as multi-protein structural scaffolds, microtubule anchors, or components of the minor spliceosome (U12 snRNP). They do not possess druggable active sites that can be activated by conventional small molecules.
- Consequently, future therapeutic efforts aiming at true disease modification must look beyond small molecules toward **AAV-mediated gene replacement therapy, antisense oligonucleotides (ASOs), or base-editing technologies**.

### Data Source Distribution & Literature Mining Dominance
In Phase 3 and Phase 4 pharmacological retrievals:
- `DGIdb` contributed ~86% (1,898 pairs) of all retrieved interactions, reflecting extensive indexing of chemogenomic literature and text-mining screens.
- `Open Targets Platform v4` (161 pairs) and `ChEMBL v2.0` (153 pairs) provided smaller, highly curated subsets featuring verified clinical trial phases, target tractability, and documented mechanism-of-action descriptions.
- The dominance of literature-mined sources underscores the vital importance of the Phase 4 multi-gate adjudication engine, which filtered out thousands of non-directional binding records down to 442 rigorously curated, signed causal hypotheses.

### Coverage Gaps in Supportive-Care Phenotypic Axes
The complete absence of qualifying records for Class H3 emphasizes that current translational oncology and rare disease databases focus heavily on tumor suppression and mitotic kinases, leaving secondary phenotypes (e.g., renal tubular ion transport defects and nephrocalcinosis) unrepresented. Addressing secondary pediatric morbidity in MVA requires building dedicated phenotypic annotation layers for pediatric rare disease complications.

---

## 8. Comprehensive Methodological Limitations

To maintain scientific integrity and transparent research standards, six major methodological limitations must be stated:

1. **Unphased Short-Read Sequencing**:
   - Short-read Illumina sequencing (~150 bp paired-end) cannot span the genomic distance separating the two prioritized `BUB1B` variants (`c.2210T>G` in exon 18 and `c.3006T>G` in exon 23).
   - Physical phase (*cis* vs *trans*) is unproven. If both variants reside in *cis* on a maternally or paternally inherited allele, the child is a carrier and another undetected variant (e.g., a deep intronic splice-disruptor, structural inversion, or non-coding regulatory deletion) must be responsible.
2. **In-Silico Predictive Nature of Computational Ensembles**:
   - Predictors such as CADD, REVEL, AlphaMissense, and SpliceAI are machine-learning algorithms trained on known variation and structural constraints.
   - High computational deleteriousness scores represent statistical probability of pathogenicity, not biochemical or functional proof.
3. **Network Association Does Not Establish Pathogenic Causation**:
   - The identification of Component 2 in STRING v12.5 and the placement of `ALK` within this network reflect known protein-protein interactions and co-expression in published cell models.
   - Network proximity does not prove that ALK activity is aberrant in this specific child's cells, nor does it guarantee that modulating ALK will rescue chromosome missegregation in vivo.
4. **Point-in-Time Knowledgebase Snapshots**:
   - All external queries reflect database releases active as of **October 2026** (STRING v12.5, Reactome v88, SIGNOR 3.0, OmniPath 2024, Open Targets v4, ChEMBL v2.0).
   - As biological knowledgebases evolve, interaction topologies and pharmacological curation will update.
5. **Lack of Independent Wet-Lab Functional Validation**:
   - The pharmacological hypotheses generated in Phase 3 and adjudicated in Phase 4 (including the 442 accepted research hypotheses) are computational predictions.
   - None of these candidates have been tested in patient-derived lymphoblastoid, fibroblast, or organoid models.
6. **Schema-Only State of Synthetic Lethality Ledger (`22_phase4_synthetic_lethality_candidates.csv`)**:
   - Artifact `22_phase4_synthetic_lethality_candidates.csv` contains 0 rows (schema-only) because the external SynLethDB integration was disabled by default in pipeline configuration to preserve strict curation standards.

---

## 9. Reproducibility & Output Artifact Catalog

### Master Flat Artifact Catalog (Artifacts 00–24)
All primary result tables are preserved directly under `data/processed/track2/` using flat, ordered filenames:

| Filename | Phase | Format | Size | Rows / Scope | Primary Contents & Description |
|---|---|---|---:|---:|---|
| `00_phase0_consolidated_candidates.csv` | Phase 0 | CSV | 10.6 MB | 2,215 | Filter-audited candidate table with rarity, impact, and MVA flags. |
| `00_phase0_provenance.json` | Phase 0 | JSON | 1.0 KB | — | Phase 0 parameter provenance, gating thresholds, and disclaimers. |
| `01_phase1_variant_mechanism_matrix.csv` | Phase 1 | CSV | 9.8 MB | 18,252 | Molecular consequence, directionality, and impact scores. |
| `02_phase1_gene_mechanism_matrix.csv` | Phase 1 | CSV | 3.8 MB | 8,876 | Gene-level aggregation, compound-het archetype, phase gating. |
| `02_phase1_provenance.json` | Phase 1 | JSON | 363 B | — | Phase 1 classification audit and candidate totals. |
| `03_phase2_seed_genes.csv` | Phase 2 | CSV | 728 KB | 1,985 | Seed genes passing candidate gating with HGNC mappings. |
| `04_phase2_string_edges.csv` | Phase 2 | CSV | 51 KB | 477 | High-confidence PPI edges (score $\ge 0.700$) from STRING v12.5. |
| `05_phase2_network_summary.csv` | Phase 2 | CSV | 283 B | — | Global topological metrics (degree, components, clustering). |
| `06_phase2_pathway_enrichment.csv` | Phase 2 | CSV | 62 KB | 84 | Reactome and g:Profiler pathway enrichment ($\text{FDR} < 0.05$). |
| `07_phase2_network_gene_context.csv` | Phase 2 | CSV | 823 KB | 1,985 | Unified gene context synthesizing genetics, topology, and pathways. |
| `08_phase2_network_provenance.json` | Phase 2 | JSON | 2.6 KB | — | Complete Phase 2 audit trail capturing API versions and endpoints. |
| `09_phase2_network_nodes.csv` | Phase 2 | CSV | 435 KB | 1,989 | Node-level metrics, component ranks, and isolation statuses. |
| `10_phase2_component_summary.csv` | Phase 2 | CSV | 5.8 KB | 124 | Ranking of connected components, sizes, and gene rosters. |
| `11_phase2_largest_component.html` | Phase 2 | HTML | 47 KB | — | Interactive Plotly visualization of Component 1 (macromolecular). |
| `12_phase2_mva_focus_network.html` | Phase 2 | HTML | 43 KB | — | Interactive Plotly visualization of Component 2 (MVA mitotic module). |
| `13_phase2_pathway_enrichment.html` | Phase 2 | HTML | 54 KB | — | Interactive Plotly bubble chart of top enriched pathways. |
| `14_phase3_drug_evidence_matrix.csv` | Phase 3 | CSV | 1.55 MB | 2,212 | Multi-source drug-target evidence matrix across 80 targets. |
| `15_phase3_safety_ranking.csv` | Phase 3 | CSV | 953 KB | 2,006 | Safety tier ranking (Tier 1: 5, Tier 2: 99, Tier 3: 1,891, Tier 4: 11). |
| `16_phase3_mechanism_action_map.csv` | Phase 3 | CSV | 33 KB | 80 | Target-level pharmacology summary and zero-hit target ledger. |
| `17_phase3_repurposing_provenance.json` | Phase 3 | JSON | 13 KB | — | Complete Phase 3 query audit trail and control outcomes. |
| `18_phase3_query_gene_panel.csv` | Phase 3 | CSV | 122 KB | 80 | Multi-tiered 80-gene panel with biological selection rationale. |
| `18_phase3_query_panel_provenance.json` | Phase 3 | JSON | 1.9 KB | — | Selection gates 0–6 execution audit and panel builder metrics. |
| `19_phase3_query_panel_citations.csv` | Phase 3 | CSV | 29 KB | 80 | Verified literature citations, PMIDs, DOIs, and OMIM entries. |
| `20_phase4_reaction_level_matrix.csv` | Phase 4 | CSV | 3.0 MB | 5,992 | Signed biochemical reactions and causal regulatory paths. |
| `21_phase4_hypothesis_classification.csv` | Phase 4 | CSV | 582 KB | 1,136 | Classification across H1, H2, H3, H0 with pathway justifications. |
| `22_phase4_synthetic_lethality_candidates.csv` | Phase 4 | CSV | 1 KB | 0 | Schema-only synthetic lethality ledger (SynLethDB disabled). |
| `23_phase4_adjudicated_safety_report.csv` | Phase 4 | CSV | 620 KB | 1,136 | Final prioritized safety, pediatric hazard, and verdict report. |
| `24_phase4_adjudication_provenance.json` | Phase 4 | JSON | 9.7 KB | — | Master provenance with SHA256 input fingerprints and gate rules. |

### Machine-Readable Provenance Tracking
Every pipeline phase emits a dedicated JSON provenance file capturing exact execution timestamps, software versions, API endpoint URLs, input SHA256 cryptographic hashes, parameter thresholds, and non-diagnostic disclaimers:
- Phase 0: `00_phase0_provenance.json` (SHA256 fingerprinting of input master tables)
- Phase 1: `02_phase1_provenance.json` (Record of mechanism classification logic)
- Phase 2: `08_phase2_network_provenance.json` (STRING v12.5 and Reactome v88 query parameters)
- Phase 3: `17_phase3_repurposing_provenance.json` and `18_phase3_query_panel_provenance.json` (Panel gates and API response counts)
- Phase 4: `24_phase4_adjudication_provenance.json` (SHA256 input hashes, SIGNOR/OmniPath edge counts, 265 resolved conflicts, and pre-registered scoring rules)

### Automated Verification Suite & Software Governance
The entire Track 2 software stack is managed with Python 3.13 and `uv`. The test suite consists of **53 automated unit and integration tests** executing cleanly under `uv run pytest`:
- `tests/test_migrate_track2_layout.py`: 6 tests verifying file paths and output structure.
- `tests/test_phase0_phase1.py`: 5 tests verifying candidate gating, rarity filters, and mechanism matrices.
- `tests/test_phase2_network.py`: 8 tests verifying STRING DB API integration and graph metrics.
- `tests/test_phase2_visualization.py`: 4 tests validating Plotly HTML rendering and component graphs.
- `tests/test_phase3_drug_evidence.py`: 7 tests verifying Open Targets, ChEMBL, and DGIdb query engines.
- `tests/test_phase3_query_panel.py`: 10 tests validating panel gates 0–6, control benchmarks, and citations.
- `tests/test_phase4_adjudication.py`: 8 tests validating signed path math, reaction adjudication, and gates G1–G7.
- `tests/test_phase4_literature.py`: 5 tests verifying Crossref/NCBI E-utilities verification and contradiction screening.

### Repository Version Pinning
- Workspace Git Commit / Tag Reference: `<commit>`
- Python Environment: Python 3.13 managed via `uv` lockfile (`uv.lock`, `pyproject.toml`).

---

## 10. Verified References & External Database Resources

### Frozen Verified Literature References
All literature citations in this report are strictly restricted to the following verified, peer-reviewed discovery publications, live-resolved and authenticated via Crossref and NCBI E-utilities on **2026-10-11**:

1. Constitutional aneuploidy and cancer predisposition caused by biallelic BUB1B mutation (Nature Genetics, 2004). DOI: 10.1038/ng1449 | PMID: 15475955 — Canonical MVA1 / BUB1B discovery.
2. Mutations in CEP57 cause mosaic variegated aneuploidy syndrome (Nature Genetics, 2011). DOI: 10.1038/ng.822 — MVA2 / CEP57 discovery.
3. Biallelic TRIP13 mutations predispose to Wilms tumor and chromosome missegregation (Nature Genetics, 2017). DOI: 10.1038/ng.3883 | PMID: 28553959 — MVA3 / TRIP13 discovery.
4. Chromosomal instability caused by mutations in the minor spliceosome component CENATAC (The EMBO Journal, 2021). DOI: 10.15252/embj.2020106536 — CENATAC/CCDC84 U12 minor spliceosome discovery.
5. Biallelic germline MAD1L1 mutations induce a syndrome of aneuploidy with high tumor susceptibility (Science Advances, 2022). DOI: 10.1126/sciadv.abq5914 — MAD1L1 mosaic aneuploidy syndrome.
6. Pathogenic variants in SLF2 and SMC5 cause segmented chromosomes and mosaic variegated aneuploidies (Atelis) (Nature Communications, 2022). DOI: 10.1038/s41467-022-34349-8 — SLF2/SMC5 Atelis syndrome discovery.
7. Germ line p53 mutations in a familial syndrome of breast cancer, sarcomas, and other neoplasms (Science, 1990). DOI: 10.1126/science.1978757 | PMID: 1978757 — TP53 Li-Fraumeni context.
8. Type 1 neurofibromatosis gene identification (Science, 1990). DOI: 10.1126/science.2134734 — NF1 genomic context.
9. The role of ATM in DNA damage responses and cancer (Oncogene, 1998). DOI: 10.1038/sj.onc.1202577 | PMID: 9916992 — ATM DNA damage response context.
10. Cell cycle regulation of Plk1 activity and subcellular localization (Journal of Cell Biology, 1995). DOI: 10.1083/jcb.129.6.1617 | PMID: 7790358 — PLK1 mitotic kinase biology.
11. A homologue of Drosophila aurora kinase is oncogenic and amplified in human colorectal cancer (EMBO Journal, 1998). DOI: 10.1093/emboj/17.11.3052 | PMID: 9606188 — AURKA aurora kinase biology.
12. Glycogen synthase kinase-3 beta phosphorylates serine 33 of p53 and activates p53 (BMC Cell Biology, 2001). DOI: 10.1186/1471-2121-2-12 | PMID: 11483158 — p53 phosphorylation endpoint mechanism.
13. Dephosphorylated hypoxia-inducible factor 1alpha as a mediator of p53-dependent apoptosis (Oncogene, 2001). DOI: 10.1038/sj.onc.1204742 | PMID: 11593383 — p53 apoptotic surveillance endpoint.
14. Phosphorylation regulates interaction/complex formation of p53 (Journal of Cellular Biochemistry, 2003). DOI: 10.1002/jcb.10569 | PMID: 12898523 — p53 phosphorylation complex regulation.

### Curated Biological Databases & Computational Knowledgebases
All pathway interactions, pharmacological annotations, and clinical phenotypes lacking a specific primary literature reference are derived directly from the following curated knowledgebase releases:
- **OMIM (Online Mendelian Inheritance in Man)**:
  - MVA1 (`#257300`, Mosaic Variegated Aneuploidy Syndrome 1, `BUB1B`)
  - MVA2 (`#614114`, Mosaic Variegated Aneuploidy Syndrome 2, `CEP57`)
  - MVA3 (`#617598`, Mosaic Variegated Aneuploidy Syndrome 3, `TRIP13`)
  - CENATAC / Microcephaly with Chromosomal Instability (`#620153`, `CENATAC`)
  - MAD1L1 Mosaic Aneuploidy & Tumor Susceptibility (`#620189`, `MAD1L1`)
  - Atelis Syndrome 1 (`#620184`, `SLF2`)
  - Atelis Syndrome 2 (`#620185`, `SMC5`)
- **STRING Database v12.5**: High-confidence protein-protein physical and functional interaction networks (`https://string-db.org`).
- **Reactome Pathway Knowledgebase v88**: Curated human biological pathways and biochemical reaction events (`https://reactome.org`). Pathways referenced include Mitotic Prometaphase (`R-HSA-68877`), Spindle Assembly Checkpoint (`R-HSA-69620`), and Resolution of Sister Chromatid Cohesion (`R-HSA-2500257`).
- **SIGNOR 3.0 (Signaling Network Open Resource)**: Curated human causal signaling network specifying directional mechanism of action and regulatory signs (`https://signor.uniroma2.it`).
- **OmniPath 2024 Curated Database**: Meta-repository of curated human intracellular signaling pathways and activity-flow interactions (`https://omnipathdb.org`).
- **Open Targets Platform v4**: Multi-evidence target-disease associations, target tractability, and clinical phases (`https://platform.opentargets.org`).
- **ChEMBL Database v2.0**: Curated database of bioactive molecules, binding affinities, and drug development phases (`https://www.ebi.ac.uk/chembl`).
- **Drug-Gene Interaction Database (DGIdb)**: GraphQL chemogenomic repository of drug-gene interactions and antineoplastic agents (`https://dgidb.org`).
- **g:Profiler (g:Convert & g:GOSt)**: Functional profiling and gene identifier mapping services (`https://biit.cs.ut.ee/gprofiler`).
