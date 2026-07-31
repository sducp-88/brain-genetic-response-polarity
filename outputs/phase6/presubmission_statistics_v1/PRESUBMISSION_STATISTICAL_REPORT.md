# Frozen pre-submission statistical strengthening

Run time: 2026-07-31T08:34:59+08:00  
Bootstrap replicates: 10,000  
Seed: 20260731

## RADC co-primary pathology multiplicity

| cohort | model | family | pathology | S | empirical_p_raw | empirical_p_Holm | passes_Holm_FWER05 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| RADC | standard_voom | two_co_primary_global_pathology_axes | CERAD | 0.1064 | 0.3195 | 0.3195 | False |
| RADC | standard_voom | two_co_primary_global_pathology_axes | Braak | 0.1903 | 0.0262 | 0.0524 | False |

Holm correction is applied only to the two standard-voom global co-primary pathology tests. Quality weights remain a sensitivity analysis rather than a second discovery family.

## Direct common-locus RADC–SEA-AD contrasts

| cell_class | pathology | common_anchor_rows | common_loci | S_RADC_common | S_SEAAD_common | Delta_S_RADC_minus_SEAAD | bootstrap_median_Delta | bootstrap_CI95_lower | bootstrap_CI95_upper | bootstrap_centered_two_sided_p | low_locus_count_exploratory | bootstrap_replicates | bootstrap_seed | bootstrap_p_BH_six_contrasts | passes_BH_FDR05 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Immune | CERAD | 24 | 20 | 0.0148 | 0.1110 | -0.0961 | -0.0930 | -0.4990 | 0.2854 | 0.6308 | False | 10000 | 20260731 | 0.7570 | False |
| Immune | Braak | 24 | 20 | 0.1374 | 0.0949 | 0.0426 | 0.0515 | -0.3158 | 0.3737 | 0.8150 | False | 10000 | 20260731 | 0.8150 | False |
| Oligo | CERAD | 9 | 9 | 0.2691 | -0.1292 | 0.3983 | 0.3975 | 0.0852 | 0.7319 | 0.0152 | True | 10000 | 20260731 | 0.0456 | True |
| Oligo | Braak | 9 | 9 | 0.3250 | -0.2093 | 0.5343 | 0.5277 | 0.2019 | 0.8928 | 0.0026 | True | 10000 | 20260731 | 0.0156 | True |
| EN | CERAD | 18 | 13 | 0.1884 | -0.1656 | 0.3540 | 0.3550 | -0.0786 | 0.7957 | 0.1141 | False | 10000 | 20260731 | 0.1711 | False |
| EN | Braak | 18 | 13 | 0.2918 | -0.0829 | 0.3747 | 0.3736 | -0.0523 | 0.8040 | 0.0850 | False | 10000 | 20260731 | 0.1700 | False |

The paired contrast resamples independent LD loci. Its two-sided bootstrap P value compares the observed mean difference with the centered bootstrap null distribution. BH correction is applied across the six fixed contrasts. Comparisons with fewer than 10 loci remain descriptive even after correction.

## Claim boundary

- A confidence interval containing zero is not evidence of equivalence.
- A significant result in one cohort and a nonsignificant result in another is not evidence of a cohort difference.
- No cell-specific, reverse, protective, compensatory, or causal claim is authorized by this analysis.
