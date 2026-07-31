# Reproduction and release-boundary QA

QA date: 31 July 2026  
Status: **PASS for privacy-bounded public release v1.0.1**

## Deterministic statistical reproduction

The aggregate-only bundle was used as an independent project root to rerun
`summarize_presubmission_v2.py`.

All five deterministic output tables matched their frozen counterparts by
exact SHA-256:

- `radc_continuous_sensitivities.csv`
- `radc_hc3_sensitivity.csv`
- `radc_seaad_dfc_direct_contrasts.csv`
- `radc_seaad_dfc_locus_components.csv`
- `radc_stage_contrasts.csv`

Result: **5/5 exact matches**.

## Figure reproduction

`make_manuscript_figures_v4.py` independently regenerated five figure
composites and 25 delivery files from the staged bundle.

- All 15 raster deliverables (preview PNG, 600-dpi PNG, and 300-dpi RGB LZW
  TIFF for each figure) matched the frozen V4 files by exact SHA-256.
- PDF and SVG hashes differed as expected because regenerated vector files
  contain time- or session-dependent metadata/identifiers.
- Exact raster identity demonstrates deterministic plotted values, labels,
  geometry, and rendered layout. Original-resolution and final-size visual
  review of all five composites is recorded in `FIGURE_VISUAL_QA_V4.md`.

Result: **15/15 raster files exact; five figure composites passed visual
review**.

## Code QA

- Python scripts parsed: 38.
- R scripts parsed: 14.
- All included Python scripts passed syntax parsing.
- All included R scripts passed `parse()`.
- Aggregate inputs copied from frozen run manifests: 167.
- Aggregate inputs rejected by participant-identifier header gate: 0.

Three local scripts were intentionally not included:

1. `audit_synapse_manifest.py`: credential-handling helper not needed for the
   public CELLxGENE analysis route.
2. `build_phase1_genetic_anchors.py`: contains a machine-specific absolute
   path; the frozen anchor-level inputs and results are included instead.
3. `summarize_radc_metadata.py`: contains a machine-specific absolute path and
   is not required for the frozen presubmission synthesis.

## Data and privacy boundary

The staged bundle contains:

- no H5AD file;
- no publisher XLSX file;
- no donor-level pseudobulk object;
- no Synapse credential or token;
- no aggregate table with an exact `donor_id`, `sample_id`,
  `participant_id`, `individual_id`, `subject_id`, `person_id`, or
  `specimen_id` header.

Third-party source data remain at their official repositories and are not
relicensed by this bundle.

## Release status

- **Completed:** the author-owned code and documentation are covered by the
  MIT License.
- Exact environment lock export.
- **Completed for the public release candidate:** ordered authors,
  affiliations, CRediT roles, funding, correspondence, and available ORCID
  identifiers are recorded in `AUTHORS.md` and `CITATION.cff`.
- The public repository excludes the journal manuscript, street address, and
  non-corresponding-author email addresses.
- The frozen `v1.0.1` release is the version supporting the journal submission
  package;
  it supersedes v1.0.0 only for Figure 2 layout and TIFF color-mode handling;
  analysis inputs, estimates, and inferential decisions are unchanged. The
  optional archival DOI remains pending.
- Manuscript-specific declarations and author approvals are maintained outside
  this public repository.
