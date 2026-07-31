# Presubmission statistical synthesis v2

Generated: 2026-07-31T09:37:28+08:00  
Locus bootstrap: 10,000; seed 20260731

## RADC primary family

| pathology | S | bootstrap_CI95_lower | bootstrap_CI95_upper | empirical_p_raw | empirical_p_Holm_primary_family | passes_primary_Holm05 |
| --- | --- | --- | --- | --- | --- | --- |
| CERAD | 0.1064 | -0.0480 | 0.2595 | 0.3195 | 0.3195 | False |
| Braak | 0.1903 | 0.0308 | 0.3441 | 0.0262 | 0.0524 | False |

No RADC co-primary test passes Holm family-wise correction.

## RADC fixed continuous sensitivities

| model | pathology | S | bootstrap_CI95_lower | bootstrap_CI95_upper | empirical_p_raw | leave_one_locus_min | leave_one_locus_max |
| --- | --- | --- | --- | --- | --- | --- | --- |
| standard | CERAD | 0.1064 | -0.0480 | 0.2595 | 0.3195 | 0.0824 | 0.1353 |
| standard | Braak | 0.1903 | 0.0308 | 0.3441 | 0.0262 | 0.1678 | 0.2225 |
| quality_weights | CERAD | 0.0780 | -0.0889 | 0.2445 | 0.5879 | 0.0524 | 0.1056 |
| quality_weights | Braak | 0.1591 | -0.0144 | 0.3272 | 0.0888 | 0.1355 | 0.1910 |
| min_cells_50 | CERAD | 0.0987 | -0.0576 | 0.2552 | 0.4065 | 0.0744 | 0.1270 |
| min_cells_50 | Braak | 0.1951 | 0.0335 | 0.3488 | 0.0273 | 0.1727 | 0.2276 |
| exclude_age89plus | CERAD | 0.1031 | -0.0313 | 0.2379 | 0.3129 | 0.0790 | 0.1268 |
| exclude_age89plus | Braak | 0.1441 | -0.0088 | 0.2920 | 0.1046 | 0.1204 | 0.1745 |
| omit_nonad | CERAD | 0.1174 | -0.0266 | 0.2564 | 0.2264 | 0.0943 | 0.1463 |
| omit_nonad | Braak | 0.1869 | 0.0421 | 0.3253 | 0.0375 | 0.1649 | 0.2172 |
| omit_log_ncells | CERAD | 0.1082 | -0.0446 | 0.2641 | 0.2847 | 0.0842 | 0.1371 |
| omit_log_ncells | Braak | 0.1938 | 0.0341 | 0.3464 | 0.0334 | 0.1714 | 0.2262 |

## RADC fixed stage contrasts

| pathology | contrast | S | empirical_p_raw | bootstrap_median | bootstrap_CI95_lower | bootstrap_CI95_upper | leave_one_locus_min | leave_one_locus_max | maximum_absolute_locus_contribution_share | loci |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CERAD | middle_vs_low | 0.0333 | 0.9985 | 0.0338 | -0.1289 | 0.1959 | 0.0075 | 0.0614 | 0.0651 | 36 |
| CERAD | high_vs_low | 0.1083 | 0.2773 | 0.1079 | -0.0597 | 0.2767 | 0.0838 | 0.1374 | 0.0601 | 36 |
| Braak | middle_vs_low | 0.1352 | 0.1797 | 0.1350 | -0.0470 | 0.3084 | 0.1115 | 0.1658 | 0.0552 | 36 |
| Braak | high_vs_low | 0.2036 | 0.0501 | 0.2048 | 0.0308 | 0.3718 | 0.1814 | 0.2346 | 0.0571 | 36 |

## Direct RADC minus SEA-AD DFC contrasts

| cell_class | pathology | common_anchor_rows | common_loci | S_RADC_common | S_SEAAD_DFC_common | Delta_S_RADC_minus_SEAAD_DFC | bootstrap_median_Delta | bootstrap_CI95_lower | bootstrap_CI95_upper | bootstrap_centered_two_sided_p | low_locus_count_exploratory | bootstrap_p_BH_six_contrasts | passes_BH05 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Immune | CERAD | 24 | 20 | 0.0148 | 0.2035 | -0.1887 | -0.1819 | -0.5346 | 0.1298 | 0.2652 | False | 0.3182 | False |
| Immune | Braak | 24 | 20 | 0.1374 | 0.0880 | 0.0495 | 0.0500 | -0.2886 | 0.3488 | 0.7604 | False | 0.7604 | False |
| Oligo | CERAD | 9 | 9 | 0.2691 | 0.0611 | 0.2080 | 0.2058 | 0.0176 | 0.4226 | 0.0425 | True | 0.1016 | False |
| Oligo | Braak | 9 | 9 | 0.3250 | -0.0331 | 0.3581 | 0.3488 | 0.0948 | 0.6658 | 0.0156 | True | 0.0936 | False |
| EN | CERAD | 18 | 13 | 0.1884 | -0.1496 | 0.3380 | 0.3327 | -0.1211 | 0.7876 | 0.1530 | False | 0.2295 | False |
| EN | Braak | 18 | 13 | 0.2918 | -0.1256 | 0.4174 | 0.4159 | -0.0040 | 0.8333 | 0.0508 | False | 0.1016 | False |

## Interpretation boundary

- Sensitivities cannot replace the primary model.
- Confidence intervals containing zero do not establish equivalence.
- Fewer than 10 loci is exploratory even after numerical correction.
- Opposed direction is not evidence of protection or compensation.
