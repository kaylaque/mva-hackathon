# Committed reports

`docs/reports/<run-name>/` is the repository location for reviewed generated reports that may be committed and pushed.

Before publication, confirm that the run is data-blind: it contains no raw or interim data, source cache, patient or sample identifier, genomic coordinate, allele, HGVS, genotype, read metric, phenotype narrative, or credential. Keep source responses and transient caches outside Git.

Each run directory must retain its command/configuration, provenance and input/output hashes, source-completeness status, and research-only/non-diagnostic scope. Aggregate tables, pathway and drug-decision tables, source-status tables, and static figures are allowed when they meet the data-blind policy.
