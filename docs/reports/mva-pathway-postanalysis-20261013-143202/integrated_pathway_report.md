# Integrated pathway report

Research-only. Annotation enrichment does not establish altered biology, causality, pathogenicity, or diagnosis.

## Result

- g:Profiler significant terms at BH q <= 0.05: **0**.
- DAVID significant annotations at Benjamini q <= 0.05: **2**.
- Shared supported stable-ID terms: **74**.
- Stable-ID terms significant in both providers: **0**.
- No supported significant MVA-specific pathway theme.

## DAVID significant annotations

- **Acetylation** — DAVID Benjamini q = **3.94e-08**.
  - Protein acetylation is common post-translational modification. It can affect protein interaction, stability, localization, chromatin, and regulation.
  - This label is broad. It does not identify one acetyltransferase pathway or prove altered acetylation activity.
  - No supported stable identifier is retained. Cross-provider comparison is **unavailable**, not g:Profiler absence.
- **Phosphoprotein** — DAVID Benjamini q = **0.00612**.
  - Phosphoprotein means recorded protein phosphorylation. Phosphorylation commonly affects signaling, protein activity, localization, degradation, cell cycle, and mitosis.
  - This label is broad. It does not prove abnormal kinase or phosphatase activity, altered phosphorylation, or spindle-checkpoint failure.
  - No supported stable identifier is retained. Cross-provider comparison is **unavailable**, not g:Profiler absence.

## Meaning of q-values

- DAVID q-values measure annotation overrepresentation against recorded DAVID background after DAVID multiple-testing correction.
- Small q-value does not measure effect size, disease relevance, pathway specificity, patient-level activity, or replication.
- Common, well-studied, multifunctional proteins can have dense modification annotations. Annotation density, background composition, overlapping gene sets, and DAVID keyword-category design can contribute to enrichment.

## MVA interpretation

- Observed MVA/SAC/chromosome-cohesin context memberships: `BUB1B`, `MAD1L1`, `MAD2L2`, `NCAPD2`, `SMC4`.
- These are observed memberships only. They do not show either DAVID label is an MVA mechanism or assign pathogenicity.
- Mechanism-specific support would need a stable mitosis/chromosome-segregation term, corrected significance, coherent contributors, complete contributor data, and independent support.
- Current outputs do not provide this support.

## Clusters and contributors

- DAVID reports 226 cluster IDs and 1,812 annotation rows.
- Cluster enrichment scores are native DAVID scores, not p-values. Clusters can overlap.
- Complete comparable contributors are unavailable for 16,455 integrated rows.
- Null Jaccard means contributor data unavailable. It does not mean disjoint genes or gene absence.

## Method limits

- DAVID and g:Profiler statistics remain provider-specific. No pooled hypothesis test was calculated.
- Different provider coverage, mappings, backgrounds, and export thresholds can produce discordance.
- Significant DAVID labels lack stable ontology mapping and authoritative term definitions in retained artifacts. Biological interpretation remains provisional.
