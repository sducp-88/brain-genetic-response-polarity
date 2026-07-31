# Release notes

## v1.0.2 - 31 July 2026

This artwork-only patch regularizes Figure 2 panel headings and forest-plot
spacing. Analysis inputs, estimates, confidence intervals, P values,
multiplicity decisions, and manuscript claims are unchanged.

### Figure changes

- Panel letters A, B, and C now use one vertical baseline and the same fixed
  offset from their respective plot areas.
- All three panel titles are centered over their plot areas and shortened
  where necessary to prevent letter-title collisions.
- Panel-A sensitivity labels use concise parallel wording.
- The A-to-B gutter was increased so panel-B labels no longer crowd the
  panel-A frame.
- The V4 vertical padding, wider panel C, and omission of redundant panel-B
  P-value annotations are retained.

### Quality control

- Figure 2 inspected at original 600-dpi resolution and at intended
  manuscript display size;
- all 15 independently regenerated raster deliverables matched V5 by exact
  SHA-256;
- the 12 non-Figure-2 raster files are byte-identical to V4;
- source script: `scripts/make_manuscript_figures_v5.py`;
- full audit: `FIGURE_VISUAL_QA_V5.md`.

## v1.0.1 - 31 July 2026

This patch release revises publication artwork without changing analysis
inputs, estimates, statistical tests, or inferential decisions.

### Figure changes

- Figure 2 now has explicit top and bottom padding for both forest plots.
- Redundant per-row raw P-value labels were removed from panel B because raw
  and Holm-adjusted P values are already consolidated in panel C.
- Panel C was widened and its summary table and decision box were respaced.
- All submission TIFF files are now explicit 24-bit RGB, 300 dpi, and
  losslessly LZW-compressed.
- Figures 1, 3, 4, and Supplementary Figure S1 retain their prior visual
  geometry; only their TIFF color mode changed from RGBA to RGB.

### Quality control

- all five figures inspected at original 600-dpi resolution and at intended
  manuscript display size;
- panel letters, plot areas, margins, gutters, legends, line weights, marker
  weights, clipping, and text overlap: pass;
- all figure files remain below the journal's 10-MB limit;
- source script: `scripts/make_manuscript_figures_v4.py`;
- full audit: `FIGURE_VISUAL_QA_V4.md`.
- privacy-boundary violations: none; Python scripts parsed: 38; R scripts
  parsed: 14; release-manifest mismatches: none.

## v1.0.0 - 31 July 2026

This is the frozen public reproducibility release supporting the journal
submission.

### Included

- analysis and quality-control scripts;
- frozen non-identifying anchor- and locus-level aggregate inputs and outputs;
- presubmission statistical summaries;
- publication figures in raster and vector formats;
- data-provenance manifests, environment metadata, and SHA-256 inventory;
- machine-readable citation and Zenodo metadata.

### Reproducibility checks

- privacy-boundary audit: pass;
- unexpected email addresses or credentials: none;
- Python syntax: 37 files parsed;
- R syntax: 14 files parsed;
- deterministic manifest mismatches: none.

### Excluded

- journal manuscript and cover letter;
- street address and non-corresponding-author email addresses;
- participant-, donor-, or specimen-level data;
- third-party H5AD, XLSX, and other non-redistributable source files;
- credentials and access tokens.

Author-owned code and documentation are released under the MIT License.
