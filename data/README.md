# Private data layout

All patient-derived files and generated analysis tables stay untracked. This repository documents their expected locations only; it does not contain or inspect their contents.

## Directories

| Location | Content | Git policy |
| --- | --- | --- |
| `input/patient/` | Original WGS VCF/index, clinical phenotype document, and any source files received for the case. | Ignored |
| `input/refs/` | Public reference archives, extracted databases, indexes, and tool-specific caches. | Ignored |
| `interim/qc/` | QC summaries from `01_run_qc.sh`. | Ignored |
| `interim/filtered/` | Filtered VCF and index. | Ignored |
| `interim/vep/` | VEP table. Put the latest chunked-processing result at `vep_high_moderate.csv`. | Ignored |
| `interim/cadd/` | `cadd_scores.csv`: candidate-only CADD result table. | Ignored |
| `interim/gnomad/` | `gnomad_candidates.csv`: gnomAD frequency/rarity table. | Ignored |
| `interim/revel/` | `revel.csv`: best REVEL score per `VARIANT_KEY`. | Ignored |
| `interim/alphamissense/` | `alphamissense.csv`: AlphaMissense result table. | Ignored |
| `interim/spliceai/` | `spliceai.csv`: strongest SpliceAI annotation per variant. | Ignored |
| `interim/alphagenome/` | `alphagenome_avi.csv`: summarized AlphaGenome AVI output. | Ignored |
| `interim/exomiser/` | Generated candidate VCF/config and normalized Exomiser results. | Ignored |
| `processed/` | Consolidated `variant_evidence_master.csv` and `prioritized_candidates.csv`. | Ignored |

## Shared join key

Every per-variant table uses `VARIANT_KEY` in the canonical GRCh38 form:

```text
CHROM:POS:REF:ALT
```

Chromosomes are unprefixed (`1`–`22`, `X`, `Y`); alleles are uppercase. One row per genomic variant is retained before evidence is merged. Transcript-level VEP rows are resolved by canonical transcript, then `PRIORITY_SCORE`.

## Stage outputs and order

1. VEP high/moderate table → `interim/vep/vep_high_moderate.csv`
2. CADD and gnomAD candidate annotations
3. Missense branch: REVEL and AlphaMissense
4. Splicing branch: SpliceAI
5. Regulatory branch: AlphaGenome AVI summary
6. Phenotype branch: Exomiser variant results
7. `scripts/run_enrichment.py` produces the processed master evidence table and triage shortlist.

CADD, REVEL, AlphaMissense, SpliceAI, AlphaGenome, and other computational predictors are research-prioritization signals, not independent clinical evidence or diagnoses. Maintain source release/version and GRCh38 provenance alongside every stage output.

## Public references

Large public reference resources remain in `input/refs/`. `scripts/download_references.sh` retrieves only public resources; it never downloads patient data. Reference archives, decompressed files, indexes, and tool-specific caches remain local and must not be committed.

## Independent pathway workflow artifacts

The pathway branch is private and research-only. `.gitignore` prevents ordinary
staging; it does **not** itself authorize an agent or process to inspect data.

| Location | Contents |
| --- | --- |
| `processed/track2/pathway/` | `p0_gene_qc.csv` through `p6_evaluation_summary.csv`, Markdown-only `p7_pathway_report.md`, and stage-keyed `pathway_provenance.json`. |
| `interim/pathway_cache/` | Endpoint-plus-gene-only-payload response cache, when enabled; no credentials. |

P0 emits `SYMBOL, HGNC_ID, ENSEMBL_ID, STATUS, IN_FOREGROUND,
PASSING_VARIANTS, TOTAL_VARIANTS`. P1 preserves provider-adjusted BH values
separately from local raw hypergeometric values. P2 records relevance,
redundancy, and HP observational-annotation flags. P3–P5 retain
mechanism/target/drug source-native evidence; P6 records descriptive STRING
comparison; P7 is a cumulative Markdown report. External services receive only
normalized gene identifiers. Variant keys, genotypes, coordinates, phenotype
text, internal scores, and patient identifiers never leave the local workflow.
