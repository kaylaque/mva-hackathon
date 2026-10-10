# MVA Patient Variant-Prioritization Report

## Scope

This report summarizes the saved outputs from `mva_vep_enrichment_analysis.ipynb`. It is a research-prioritization report, not a diagnosis or treatment recommendation. Genomic coordinates and sample identifiers are omitted.

## Main question

Which rare variants and genes best fit a mosaic variegated aneuploidy (MVA) or chromosome-instability mechanism, and which findings warrant follow-up for mechanism and drug-discovery research?

## Patient research profile

The strongest molecular hypothesis is possible recessive `BUB1B`-related MVA. Two rare `BUB1B` variants were prioritized:

| Variant | Consequence | Population frequency | Supporting evidence | Interpretation |
|---|---|---:|---|---|
| `c.2210T>G`, `p.Leu737*` | Stop gained | $9.83 \times 10^{-5}$ | CADD 36; ClinVar pathogenic/likely pathogenic; imported ACMG pathogenic | Strong loss-of-function candidate |
| `c.3006T>G`, `p.Asn1002Lys` | Missense | $8.47 \times 10^{-7}$ | CADD 24.5; SIFT 0.01; PolyPhen 0.997; AlphaMissense 0.9229; REVEL 0.472 | Rare damaging candidate with predictor disagreement |

Together, these variants form a plausible two-variant `BUB1B` candidate. The imported `Compound_Het_Candidate (2_hits)` label does not prove compound heterozygosity. Genotype, zygosity, phase, and parental segregation were unavailable in the processed table.

## Gene priorities

### Priority 1: `BUB1B`

`BUB1B` is a canonical MVA gene. It has the strongest convergence of disease mechanism, rarity, consequence severity, and computational evidence. Immediate follow-up should confirm both calls, determine genotype, and establish whether the variants are in trans.

### Priority 2: `MAD1L1`

`MAD1L1 p.Arg59Cys` was rare, with population frequency approximately 0.0006, and relevant to chromosome-segregation biology. Imported clinical assertions were conflicting, including benign and pathogenic labels. It remains a manual-review candidate, not a result comparable in confidence to the `BUB1B` pair.

### Conditional or low-priority findings

- `CENATAC p.Phe242Leu` was classified as a VUS and had population frequency approximately 0.0072. It requires strong phenotype and inheritance support.
- Observed `SLF2`, `SMC5`, `AURKB`, `CENPE`, and `KNL1` variants were common and predominantly benign.
- `CEP57` and `TRIP13` were not represented in the filtered table. Their absence does not prove that these genes are negative; sequencing coverage and upstream filtering require review.
- Phenotype-matched pathogenic or likely pathogenic variants in non-MVA genes still require review. They may represent an alternative diagnosis, a dual diagnosis, or an incidental finding.

## Dataset summary

| Metric | Result |
|---|---:|
| Annotation rows | 18,252 |
| Distinct variants | 18,252 |
| Represented genes | 8,876 |
| Population frequency greater than 1% | 16,187 (88.69%) |
| Unknown population frequency | 1,312 (7.19%) |
| Population frequency at or below 1% | 753 (4.13%) |
| Population frequency at or below 0.1% | 310 (1.70%) |
| Population frequency at or below 0.01% | 96 (0.53%) |
| Missense annotations | 13,758 |
| Frameshift annotations | 517 |
| Stop-gained annotations | 221 |

Most records were common variants. Rare variants form the relevant review set.

## What processed variant evidence can establish

The processed table supports:

- variant and gene identification;
- transcript and predicted consequence review;
- population-frequency filtering;
- imported ClinVar and ACMG evidence review;
- CADD, SIFT, PolyPhen, REVEL, AlphaMissense, and SpliceAI prioritization;
- selection of variants for manual and experimental follow-up.

These tools support hypotheses. They do not establish causality or drug response.

## Critical missing evidence

The processed table did not provide:

- patient genotype and zygosity;
- sequencing depth, allele balance, genotype quality, and variant filter status;
- cis/trans phase;
- parental segregation;
- confirmed observed HPO terms;
- confirmed cytogenetic mosaic aneuploidy;
- tissue-specific mosaic fraction;
- RNA abundance, splicing, or nonsense-mediated decay;
- BUBR1 protein abundance, localization, stability, or activity;
- spindle-checkpoint functional measurements;
- CNV or structural-variant evidence;
- patient-cell drug response.

`DP_AG`, `DP_AL`, `DP_DG`, and `DP_DL` are SpliceAI delta-position offsets, not sequencing read depth. `PHENOTYPE_PANEL` is an upstream annotation, not confirmation of patient symptoms. `INHERITANCE_MODEL` is an upstream candidate label, not proof of inheritance or phase.

## Importance for drug discovery

Sequence and annotation data indicate what may be altered. Drug selection requires evidence showing how the alteration changes RNA, protein, cells, and disease-relevant biology.

The intervention strategy depends on mechanism:

- absent transcript or protein may favor gene replacement, gene correction, or RNA editing research;
- a surviving premature-stop transcript may permit a readthrough screen;
- an unstable but partly functional missense protein may permit a stabilization screen;
- a stable but inactive protein may require correction or pathway compensation;
- failure to confirm `BUB1B` causality requires reassessment before target-specific drug work.

Because BUBR1 controls chromosome segregation, broad mitotic manipulation carries substantial toxicity and chromosome-instability risk. Functional rescue and safety must be measured together.

## Essential next evidence

1. Confirm both `BUB1B` variants with genotype and sequencing-quality evidence from the original VCF and BAM/CRAM.
2. Establish phase and parental segregation; determine whether the variants are in trans.
3. Confirm relevant transcript and HGVS annotations.
4. Measure allele-specific RNA, splicing, and nonsense-mediated decay.
5. Measure BUBR1 protein abundance, localization, stability, and checkpoint activity.
6. Confirm mosaic aneuploidy using documented cytogenetic methods, tissue, cell count, and mosaic fraction.
7. Test patient-derived cells for chromosome missegregation, micronuclei, and checkpoint failure.
8. Correct the candidate variants in patient cells and test whether the cellular phenotype is rescued.
9. Use the rescued phenotype as the endpoint for compound or modality screening.
10. Measure genotoxicity, chromosome stability, normal-cell toxicity, and transformation risk.

## Data-quality limitations

- Imported rarity flags disagreed with numeric population frequencies for many rows. Numeric frequency requires precedence during review.
- Combined consequence strings were not fully separated because the analysis split `&` while some values used commas.
- Some transcript and SpliceAI gene mismatch counts treated missing string values as real values.
- Some alleles used `-` representations instead of normalized VCF alleles, creating normalization and join risk.
- CADD, REVEL, AlphaMissense, and SpliceAI coverage was sparse. Missing prediction is not benign evidence.
- ClinVar condition, review status, assertion date, and submitter evidence were unavailable.
- Claimed GRCh38 provenance was not independently verified.

## Conclusion

Current evidence supports `BUB1B` as the highest-priority research gene. The combination of `p.Leu737*` and `p.Asn1002Lys` is a plausible rare recessive candidate, but not a confirmed molecular diagnosis.

The decisive next evidence is not another prediction score. It is confirmation of valid calls, appropriate genotype, trans phase, phenotype and cytogenetic fit, altered BUBR1 function, and rescue after genetic correction. Only then can a patient-specific drug mechanism be selected with confidence.
