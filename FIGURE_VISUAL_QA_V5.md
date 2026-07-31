# Figure visual and technical QA: layout V5

QA date: 31 July 2026  
Target journal: BMC Genomics  
Status: **PASS**

## Scope and evidence boundary

The V5 revision changes Figure 2 layout and wording only. No analysis input,
estimate, confidence interval, P value, multiplicity decision, or manuscript
claim was changed. Figures 1, 3, 4, and Supplementary Figure S1 are
byte-identical to their V4 raster deliverables.

## Figure 2 changes

- Positioned panel letters A, B, and C at the same vertical baseline and at
  the same fixed offset from the left edge of each plot area.
- Centered each panel title over its own plot area and shortened the titles to
  avoid collision with the panel letters.
- Replaced long, structurally uneven sensitivity labels with concise parallel
  labels: `Primary`, `Q-weighted`, `≥50 nuclei`, `Age <89`,
  `No non-AD adj.`, and `No depth adj.`.
- Increased the A-to-B gutter from 0.40 to 0.58 grid units and reduced the
  panel-B tick-label padding, creating a clear safety gap between the
  panel-B labels and the panel-A frame.
- Retained the V4 vertical forest-plot padding, the wider panel C, and the
  decision to omit redundant P-value annotations from panel B.

## Visual inspection

The Figure 2 600-dpi PNG was inspected at original resolution. Its 180-dpi
preview was inspected at the intended 170-mm full-page display width.

| Check | Original-resolution review | Final-size review | Result |
|---|---|---|---|
| Panel header grid | Letters share one baseline and one left-offset rule; titles are centered over their plot areas | Letter-title pairs remain separated and visually regular | PASS |
| Panel A labels | All six labels are unclipped, right-aligned, and use concise parallel wording | The label column reads as a coherent sensitivity-analysis list | PASS |
| A-to-B gutter | Panel-B labels do not touch or visually merge with the panel-A frame | A visible white-space safety zone remains after reduction | PASS |
| Forest plots | First and last rows retain padding; markers and intervals do not touch the panel edges | Null lines remain behind the data and all intervals remain legible | PASS |
| Panel C | Columns, rules, values, and decision card remain aligned and unclipped | Summary values and the Holm decision remain readable | PASS |

The V5 Figure 2 contains no title collision, label intrusion, clipping,
crowded first or last row, or redundant panel-B P-value annotation.

## Deterministic reproduction

An independent rerun generated five previews, five 600-dpi PNGs, and five
300-dpi TIFFs. All 15 raster files matched the frozen V5 outputs by exact
SHA-256. The 12 non-Figure-2 raster files are also byte-identical to V4.

## Technical audit

BMC Genomics instructions accessed 31 July 2026 specify single composite
files for multipanel figures, accepted TIFF/PDF/PNG formats, 85- or 170-mm
final widths, approximately 300 dpi at final size, LZW compression for TIFF,
and a 10-MB maximum per figure:

https://link.springer.com/journal/12864/submission-guidelines

All final TIFF files are 24-bit RGB, 300 dpi, and LZW-compressed:

| File | Pixels | Bytes | SHA-256 |
|---|---:|---:|---|
| Figure 1 | 2097 × 1384 | 568330 | `7f0f3495b8264c420c7ddaf254748d8df70538edf9a904d44613bb5644400e3a` |
| Figure 2 | 2188 × 1324 | 486622 | `508f62457ab107b38b710fba8bdda8a8d66374245e952d794b153f3643e2e04d` |
| Figure 3 | 2134 × 1342 | 596922 | `bf2229d85f65a70c2d7dd190f207fc263be84129d7316e6c921829e2992d82b8` |
| Figure 4 | 2160 × 1286 | 390870 | `77a1bbe14a187fbd3ba7a1aabaa3b0d2822703903dc8f7d944bc2c93348c430e` |
| Supplementary Figure S1 | 1570 × 1526 | 267324 | `6a25cd28260f15998878fe85a620ecc9537941af830d766d649b36fcfea2ff2c` |

Editable PDF and SVG masters, 600-dpi PNG files, 180-dpi previews, the
generating script, and the frozen aggregate inputs are retained in the public
reproducibility bundle.
