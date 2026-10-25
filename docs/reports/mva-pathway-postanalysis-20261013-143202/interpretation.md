# Biological interpretation

Research-only. No diagnostic or causal conclusion.

## Conclusion

- No stable-ID term is significant in both DAVID and g:Profiler.
- No significant MVA-specific pathway theme is supported.
- DAVID detects two broad protein-modification annotations. They are hypothesis-generating, not mechanism evidence.

## DAVID annotations

- **Acetylation** — Benjamini q = `3.94e-08`.
  - Common protein modification. Can affect protein interaction, stability, localization, chromatin, and regulation.
  - Does not prove altered acetylation activity or one acetylation pathway.
- **Phosphoprotein** — Benjamini q = `0.00612`.
  - Recorded protein phosphorylation. Common in signaling, cell cycle, and mitosis.
  - Does not prove altered phosphorylation, kinase activity, phosphatase activity, or spindle-checkpoint failure.

## Cross-provider status

- Both DAVID labels lack supported stable identifiers in retained export.
- g:Profiler comparison is unavailable. It is not evidence of g:Profiler absence.
- No pooled hypothesis test was calculated.

## MVA context

- Observed context genes: `BUB1B`, `MAD1L1`, `MAD2L2`, `NCAPD2`, `SMC4`.
- Membership is contextual only. It does not link broad DAVID labels to MVA mechanism or establish pathogenicity.
- Mechanism claim needs stable mitosis/chromosome-segregation term, corrected significance, coherent contributors, and independent support.

## Limits

- DAVID cluster export: 226 cluster IDs; 1,812 annotation rows.
- Cluster scores are not p-values. Clusters can overlap.
- Complete comparable contributors unavailable for 16,455 rows.
- Null Jaccard means unavailable contributor data, not disjoint genes.
- Annotation density, provider coverage, mappings, backgrounds, and export thresholds can affect enrichment.
- Stable ontology mapping and authoritative term definitions are unavailable for significant DAVID labels. Interpretation remains provisional.
