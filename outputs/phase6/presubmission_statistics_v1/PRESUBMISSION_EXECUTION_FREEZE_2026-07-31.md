# Pre-submission strengthening execution freeze

Freeze date: 2026-07-31  
Status: frozen before executing or inspecting outputs from the new script

## Authorized script

`run_presubmission_statistics.py`

SHA-256:

`D39107FDDADACC3390E87C1E3F784F7545E62A387A2C1FD83E572942CE14E325`

## Fixed parameters

- Bootstrap unit: independent LD locus
- Bootstrap replicates: 10,000
- Seed: 20260731
- RADC primary multiplicity: Holm adjustment across standard-voom CERAD and
  Braak global tests
- Direct comparison family: 3 fixed cell classes × 2 fixed pathology axes
- Direct-comparison multiplicity: BH across six tests
- Low-locus warning: fewer than 10 common loci

## Fixed cells and pathology axes

- Cells: Immune, Oligo, EN
- Pathology: CERAD, Braak
- RADC model for direct comparison: standard voom
- SEA-AD model for this execution: existing three-region-summed standard voom

The region-summed SEA-AD model is retained only as the legacy comparison.
DFC will become the matched-region primary follow-up after the region-specific
expression models are run under a separate frozen execution.

## Expected outputs

1. `radc_primary_multiplicity.csv`
2. `locus_bootstrap_intervals.csv`
3. `radc_seaad_direct_common_locus_contrasts.csv`
4. `radc_seaad_common_locus_components.csv`
5. `PRESUBMISSION_STATISTICAL_REPORT.md`
6. `run_manifest.json`

Every authorized contrast will be reported. No output may be deleted or
replaced because of its direction, confidence interval or P value.

## Interpretation limits

- A confidence interval containing zero is not evidence of equivalence.
- A significant result in one cohort and a nonsignificant result in another
  is not evidence of a direct cohort difference.
- No causal, protective, compensatory, reverse-mechanism or cell-specific
  claim is authorized.

