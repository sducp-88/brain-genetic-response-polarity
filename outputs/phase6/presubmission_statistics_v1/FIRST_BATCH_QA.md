# First-batch pre-submission statistics: quality-control record

Date: 2026-07-31  
Frozen run: `outputs/phase6/presubmission_statistics_v1`  
Bootstrap: 10,000 independent-LD-locus resamples  
Seed: 20260731

## Reproducibility checks

- The executed script SHA-256 was
  `d39107fddadacc3390e87c1e3f784f7545e62a387a2c1fd83e572942ce14e325`.
- All 21 registered input hashes match the files currently on disk.
- All five registered machine-generated output hashes match the run manifest.
- The standard-voom RADC point estimates reproduced the frozen source tables:
  - CERAD: absolute difference `1.39e-17`;
  - Braak: absolute difference `0`.
- The run completed without missing registered outputs.

## Statistical checks

### Co-primary RADC tests

| Pathology | S | Raw empirical P | Holm P | Confirmatory interpretation |
|---|---:|---:|---:|---|
| CERAD | 0.1064 | 0.3195 | 0.3195 | Not supported |
| Braak | 0.1903 | 0.0262 | 0.0524 | Does not pass the frozen family-wise threshold |

The nominal Braak result is therefore not a multiplicity-controlled positive
finding. Quality-weighted voom remains a sensitivity analysis and does not
create a second discovery family.

### Direct common-locus RADC minus SEA-AD contrasts

| Cell class | Pathology | Common loci | Delta S | Locus-bootstrap 95% CI | BH P | Status |
|---|---|---:|---:|---:|---:|---|
| Immune | CERAD | 20 | -0.0961 | -0.4990 to 0.2854 | 0.7570 | Inconclusive |
| Immune | Braak | 20 | 0.0426 | -0.3158 to 0.3737 | 0.8150 | Inconclusive |
| Oligo | CERAD | 9 | 0.3983 | 0.0852 to 0.7319 | 0.0456 | Exploratory: fewer than 10 loci |
| Oligo | Braak | 9 | 0.5343 | 0.2019 to 0.8928 | 0.0156 | Exploratory: fewer than 10 loci |
| EN | CERAD | 13 | 0.3540 | -0.0786 to 0.7957 | 0.1711 | Inconclusive |
| EN | Braak | 13 | 0.3747 | -0.0523 to 0.8040 | 0.1700 | Inconclusive |

The two Oligo contrasts pass the numerical BH threshold but remain descriptive
under the pre-specified low-locus rule. They cannot support a robust
cell-specific cohort-difference claim.

## Claim boundary after first-batch QA

The defensible result is limited and analysis-context-dependent agreement
between genetically anchored and disease-associated expression directions.
The current results do not authorize claims of protection, compensation,
reverse causation, cell-specific mechanism, causal mediation, equivalence, or
independent replication.

## QA disposition

**PASS WITH INTERPRETIVE RESTRICTION.** The numerical outputs are reproducible
and suitable for downstream figure and manuscript generation. The manuscript
must treat the RADC Braak result as nominal only and the Oligo direct contrasts
as low-locus exploratory observations.
