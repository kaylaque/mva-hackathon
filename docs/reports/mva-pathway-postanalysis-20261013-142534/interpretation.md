# Biological interpretation

Research-only; this descriptive post-analysis is not diagnostic and does not establish causality.

## Observed DAVID signal

At DAVID Benjamini q <= 0.05, the saved chart contains two significant annotations:

- **Acetylation** — native DAVID Benjamini q = 3.94e-08; no matched g:Profiler stable identifier, so the g:Profiler state is `absent`.
- **Phosphoprotein** — native DAVID Benjamini q = 0.00612; no matched g:Profiler stable identifier, so the g:Profiler state is `absent`.

Neither label has a supported stable ontology identifier in the retained exports. They are therefore retained as `unmapped_term` / `review_required`, not treated as provider-specific biological absence or cross-provider confirmation. No stable-ID term was significant in both providers.

## Mechanistic context

The observed support table contains memberships in the shared MVA/SAC/chromosome-cohesin context sets for BUB1B, MAD1L1, MAD2L2, NCAPD2, and SMC4. That overlap is contextual observed membership only; it does not show that either broad DAVID annotation is a mosaic-variegated-aneuploidy mechanism or establishes a pathogenicity score.

No significant, stable-ID MVA-consistent pathway theme is established by these saved outputs. The two DAVID-significant broad annotation labels are off-mechanism/non-specific for this purpose unless a separately verified ontology definition and stable term mapping supports a more specific claim.

## Cluster and contributor limits

The DAVID cluster export has 226 distinct cluster IDs across 1,812 annotation rows. Cluster enrichment scores are native DAVID scores, not p-values, and clusters may overlap. The significant annotation labels were not retained as stable-ID mapped cluster terms, so no supporting MVA-specific cluster can be claimed for them.

Complete comparable contributors were unavailable for 16,455 integrated rows. Null Jaccard values indicate incomplete contributor availability, not gene-level disjointness or absence. Provider coverage, mappings, effective backgrounds, and export thresholds may differ; no pooled significance test was performed.
