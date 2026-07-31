# Aggregate-only reproducibility bundle

Status: **MIT-licensed public reproducibility release candidate (`v0.1.0`)**.

This bundle supports two levels of review:

1. The final presubmission synthesis and publication figures can be regenerated
   from frozen, non-identifying, anchor- or locus-level aggregate inputs.
2. Portable source snapshots for the upstream donor-level pseudobulk and model
   pipeline are provided where they contain no local credentials or absolute
   machine paths. Large third-party data must be downloaded from the official
   sources listed in `provenance/`.

## What is included

- Frozen aggregate inputs and outputs listed by the original run manifests.
- Presubmission statistical synthesis scripts and the V3 figure script.
- Portable upstream Python and R source files.
- Final aggregate tables, figure files, and figure-generation QA.
- Exact public dataset/version identifiers, local verification records, SEA-AD
  object ETags, and a SHA-256 manifest for this bundle.

## What is intentionally excluded

- CELLxGENE H5AD files.
- SingleBrain publisher supplementary workbooks.
- SEA-AD H5AD objects.
- Donor-level pseudobulk matrices, donor metadata, weights, or identifiers.
- Synapse credentials and access helpers.
- Local monitoring/download wrappers containing machine-specific absolute
  paths.

The excluded third-party data remain available from their official sources;
this bundle does not grant new rights to redistribute them.

## Reproduce the frozen statistical synthesis

From the bundle root:

```bash
python scripts/summarize_presubmission_v2.py \
  --project . \
  --output-dir outputs/reproduced/presubmission_synthesis_v2
```

The deterministic CSV outputs should match the corresponding files in
`outputs/phase6/presubmission_synthesis_v2/`. Timestamps and absolute paths in
JSON run manifests are expected to differ.

## Reproduce the figures

```bash
python scripts/make_manuscript_figures_v3.py \
  --project . \
  --output-dir outputs/reproduced/figures_manuscript_v3
```

Raster output can be compared pixel-for-pixel; PDF/SVG metadata may differ
between software versions even when the rendered figure is unchanged.

## Environment

The completed analysis used Python and R in WSL2. `environment.yml` records
the direct packages exported from the analysis environment. Before archival
deposition, an explicit platform lock should also be exported so exact
transitive builds can be reconstructed. Core Python packages are `numpy`,
`pandas`, `scipy`, `matplotlib`, and `Pillow`; upstream work also uses
`anndata`. Core R packages include `limma`, `variancePartition`, and
model-specific dependencies declared in the scripts.

## Citation, authorship, and funding

The ordered authors are Fan Li, Suqin Jin, Feifei Feng, Ping Wang, Hui Yang,
and Peng Cheng (corresponding author). Machine-readable citation metadata are
provided in `CITATION.cff`; affiliations, the corresponding-author contact,
contribution roles, funding, and the AI-use statement are recorded in
`AUTHORS.md`.

This work was supported by the Natural Science Foundation of Shandong Province
(ZR2023MH340; principal investigator: Hui Yang).

## AI-assisted code and language editing

OpenAI Codex was used solely to assist with analysis-code implementation and
debugging and with language editing. The authors determined the scientific
questions, methods, interpretation, and conclusions; reviewed and approved the
code and research outputs; and take responsibility for the work.

## Licence and release boundary

Author-owned software and documentation are released under the MIT License;
see `LICENSE` and `LICENSE_DECISION_RECORD.md`. This licence does not apply to
excluded third-party data.

This public repository intentionally excludes the journal manuscript, street
address, and non-corresponding-author email addresses. Journal submission
declarations are maintained in the separate private submission package. An
archival DOI will be added after the tagged release is deposited in Zenodo.
