# Figure visual and technical QA: layout V4

QA date: 31 July 2026  
Target journal: BMC Genomics  
Status: **PASS**

## Scope and evidence boundary

The Figure 2 revision changes layout only. No analysis input, estimate,
confidence interval, P value, multiplicity decision, or manuscript claim was
changed. Figures 1, 3, 4, and Supplementary Figure S1 were regenerated from
the same frozen aggregate inputs and retain their previous visual geometry.

## Figure 2 changes

- Increased the canvas height from 4.35 to 4.65 inches.
- Added explicit and more balanced top/bottom limits to panels A and B.
- Removed redundant raw P-value annotations from panel B. Raw and
  Holm-adjusted P values remain in panel C and in the manuscript text.
- Reallocated the three-column grid to widen panel C by approximately one
  third at fixed full-page width.
- Respaced the summary rows and decision card and standardized the expression
  `Holm α = 0.05`.

## Visual inspection

Each 600-dpi PNG was inspected at original resolution. Each 180-dpi preview
was inspected at the intended 170-mm full-page display width.

| Figure | Original-resolution review | Final-size review | Result |
|---|---|---|---|
| Figure 1 | Panel letters and titles align; workflow arrows, boxes, table rows, and interpretation text are unclipped | All text remains legible; whitespace supports the two-level hierarchy | PASS |
| Figure 2 | Forest intervals and markers clear of all edges; null lines remain behind data; no P-value text in panel B; panel C columns and decision card have adequate spacing | Panel headings align; no collisions or wrapping; all numerical values are legible | PASS |
| Figure 3 | Heatmap cells, signed values, daggers, colorbar, and footnote are fully visible | Panel widths and headings are balanced; low-locus warning remains readable | PASS |
| Figure 4 | Markers, paired connectors, confidence intervals, legends, reference lines, and footnotes are unobstructed | First and last forest rows retain adequate padding; panel gutters and axes remain readable | PASS |
| Supplementary Figure S1 | All points, open/filled encodings, identity line, axes, and legends are visible | Shape and fill encodings remain interpretable without reliance on color alone | PASS |

No residual plot titles, clipping, text overlap, neighboring-panel intrusion,
or inconsistent panel-letter placement was observed.

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
| Figure 2 | 2176 × 1324 | 506154 | `a02da6cb3b69485532231093643263b993d48b296d251feb567c02930b377ec2` |
| Figure 3 | 2134 × 1342 | 596922 | `bf2229d85f65a70c2d7dd190f207fc263be84129d7316e6c921829e2992d82b8` |
| Figure 4 | 2160 × 1286 | 390870 | `77a1bbe14a187fbd3ba7a1aabaa3b0d2822703903dc8f7d944bc2c93348c430e` |
| Supplementary Figure S1 | 1570 × 1526 | 267324 | `6a25cd28260f15998878fe85a620ecc9537941af830d766d649b36fcfea2ff2c` |

Editable PDF and SVG masters, 600-dpi PNG files, 180-dpi previews, the
generating script, and the frozen aggregate inputs are retained in the public
reproducibility bundle.
