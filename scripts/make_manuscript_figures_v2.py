#!/usr/bin/env python3
"""Create restrained, publication-ready figures for the v2 synthesis."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
from PIL import Image


PALETTE = {
    "navy": "#1F3A4D",
    "blue": "#4F6D7A",
    "slate": "#6E7781",
    "ochre": "#B38A58",
    "burgundy": "#8A5A5A",
    "pale": "#F5F6F4",
    "ink": "#263238",
    "grid": "#D7DAD8",
    "light_blue": "#A9BAC2",
    "light_ochre": "#D8C4A5",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def set_style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 9,
            "axes.titlesize": 11,
            "axes.labelsize": 9.5,
            "axes.edgecolor": PALETTE["slate"],
            "axes.linewidth": 0.8,
            "axes.facecolor": "white",
            "figure.facecolor": "white",
            "savefig.facecolor": "white",
            "xtick.color": PALETTE["ink"],
            "ytick.color": PALETTE["ink"],
            "text.color": PALETTE["ink"],
            "axes.labelcolor": PALETTE["ink"],
            "axes.titlecolor": PALETTE["ink"],
            "grid.color": PALETTE["grid"],
            "grid.linewidth": 0.6,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
        }
    )


def panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(
        -0.08,
        1.05,
        label,
        transform=ax.transAxes,
        fontsize=13,
        fontweight="bold",
        color=PALETTE["ink"],
        va="top",
    )


def save_figure(
    fig: plt.Figure,
    output: Path,
    name: str,
    registry: list[dict[str, object]],
    adjust_layout: bool = True,
) -> None:
    vector_dir = output / "vector"
    raster_dir = output / "raster"
    preview_dir = output / "preview"
    for directory in (vector_dir, raster_dir, preview_dir):
        directory.mkdir(parents=True, exist_ok=True)
    if adjust_layout:
        fig.subplots_adjust(left=0.07, right=0.97, top=0.91, bottom=0.12)
    pdf = vector_dir / f"{name}.pdf"
    svg = vector_dir / f"{name}.svg"
    png = raster_dir / f"{name}_600dpi.png"
    preview = preview_dir / f"{name}_preview.png"
    tiff = raster_dir / f"{name}_300dpi_lzw.tiff"
    fig.savefig(pdf, bbox_inches="tight")
    fig.savefig(svg, bbox_inches="tight")
    fig.savefig(png, dpi=600, bbox_inches="tight")
    fig.savefig(preview, dpi=180, bbox_inches="tight")
    with Image.open(png) as image:
        target = image.resize(
            (max(1, image.width // 2), max(1, image.height // 2)),
            Image.Resampling.LANCZOS,
        )
        target.save(
            tiff,
            format="TIFF",
            compression="tiff_lzw",
            dpi=(300, 300),
        )
        pixel_size = image.size
    plt.close(fig)
    for path in (pdf, svg, png, preview, tiff):
        registry.append(
            {
                "figure": name,
                "path": str(path),
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
                "pixel_size_600dpi_png": pixel_size if path == png else None,
            }
        )


def add_box(
    ax: plt.Axes,
    xy: tuple[float, float],
    width: float,
    height: float,
    title: str,
    body: str,
    face: str,
) -> None:
    patch = FancyBboxPatch(
        xy,
        width,
        height,
        boxstyle="round,pad=0.012,rounding_size=0.018",
        linewidth=0.9,
        edgecolor=PALETTE["slate"],
        facecolor=face,
    )
    ax.add_patch(patch)
    ax.text(
        xy[0] + width * 0.05,
        xy[1] + height * 0.72,
        title,
        fontsize=9.5,
        fontweight="bold",
        color=PALETTE["navy"],
    )
    ax.text(
        xy[0] + width * 0.05,
        xy[1] + height * 0.46,
        body,
        fontsize=8,
        color=PALETTE["ink"],
        va="top",
        linespacing=1.25,
    )


def figure_design(output: Path, registry: list[dict[str, object]]) -> None:
    fig = plt.figure(figsize=(14.5, 8.2))
    grid = fig.add_gridspec(2, 2, height_ratios=[1.03, 0.97], wspace=0.20)
    ax_a = fig.add_subplot(grid[0, :])
    ax_b = fig.add_subplot(grid[1, 0])
    ax_c = fig.add_subplot(grid[1, 1])
    ax_a.set_axis_off()
    panel_label(ax_a, "A")
    ax_a.set_title(
        "Study logic: compare genetic-expression direction with postmortem expression direction",
        loc="left",
        pad=12,
        fontweight="bold",
    )
    boxes = [
        (
            0.02,
            "Genetic anchors",
            "Cell-class eQTL\ncolocalization + MR\nG1/G2 evidence",
            "#E8EEF0",
        ),
        (
            0.27,
            "Postmortem effects",
            "MSSM case-control\nRADC CERAD/Braak\nSEA-AD DFC/MEC/MTG",
            "#EEF0ED",
        ),
        (
            0.52,
            "Uncertainty propagation",
            "Probability of aligned signs\nAnchor → locus averaging\nEqual weight per LD locus",
            "#F1ECE6",
        ),
        (
            0.77,
            "Falsification tests",
            "Matched-null permutations\nLocus bootstrap + LOO\nRegion/composition sensitivity",
            "#EAECEF",
        ),
    ]
    for x, title, body, face in boxes:
        add_box(ax_a, (x, 0.25), 0.20, 0.50, title, body, face)
    for x in (0.23, 0.48, 0.73):
        ax_a.add_patch(
            FancyArrowPatch(
                (x, 0.50),
                (x + 0.035, 0.50),
                arrowstyle="-|>",
                mutation_scale=12,
                linewidth=1.2,
                color=PALETTE["slate"],
            )
        )

    panel_label(ax_b, "B")
    ax_b.set_title("Evidence hierarchy and statistical unit", loc="left", fontweight="bold")
    rows = [
        ("Genetic anchor", "Independent LD locus", "Direction evidence"),
        ("MSSM", "Donor pseudobulk", "End-stage disease state"),
        ("RADC", "Donor pseudobulk", "Pathology continuum"),
        ("SEA-AD", "Donor × region pseudobulk", "Cross-region sensitivity"),
    ]
    ax_b.set_xlim(0, 3)
    ax_b.set_ylim(-0.5, len(rows) - 0.5)
    for index, (source, unit, role) in enumerate(rows[::-1]):
        y = index
        ax_b.add_patch(
            FancyBboxPatch(
                (0.05, y - 0.32),
                2.82,
                0.64,
                boxstyle="round,pad=0.01,rounding_size=0.025",
                facecolor="white" if index % 2 else PALETTE["pale"],
                edgecolor=PALETTE["grid"],
                linewidth=0.7,
            )
        )
        ax_b.text(0.15, y, source, va="center", fontweight="bold", color=PALETTE["navy"])
        ax_b.text(1.02, y, unit, va="center", fontsize=8.4)
        ax_b.text(2.03, y, role, va="center", fontsize=8.4)
    ax_b.text(0.15, 3.43, "Source", fontweight="bold", fontsize=8)
    ax_b.text(1.02, 3.43, "Biological unit", fontweight="bold", fontsize=8)
    ax_b.text(2.03, 3.43, "Inferential role", fontweight="bold", fontsize=8)
    ax_b.axis("off")

    panel_label(ax_c, "C")
    ax_c.set_title("Interpretation is deliberately bounded", loc="left", fontweight="bold")
    ax_c.set_xlim(0, 1)
    ax_c.set_ylim(0, 1)
    ax_c.set_axis_off()
    ax_c.text(
        0.04,
        0.78,
        r"$S = 2\,P(\mathrm{aligned\ signs}) - 1$",
        fontsize=18,
        color=PALETTE["navy"],
    )
    statements = [
        ("S > 0", "probabilistic sign alignment", PALETTE["navy"]),
        ("S < 0", "probabilistic sign opposition", PALETTE["burgundy"]),
        ("S ≈ 0", "directional uncertainty", PALETTE["slate"]),
    ]
    for index, (symbol, meaning, color) in enumerate(statements):
        y = 0.58 - index * 0.14
        ax_c.scatter([0.08], [y], s=80, color=color, edgecolor="white", linewidth=0.8)
        ax_c.text(0.14, y, symbol, fontweight="bold", color=color, va="center")
        ax_c.text(0.27, y, meaning, va="center")
    ax_c.text(
        0.04,
        0.08,
        "Not authorized: protection, compensation, causal mediation,\n"
        "equivalence, or fully independent genetic replication.",
        fontsize=8.7,
        color=PALETTE["burgundy"],
        linespacing=1.4,
    )
    fig.suptitle(
        "Genetically anchored and postmortem expression directions",
        x=0.07,
        ha="left",
        fontsize=15,
        fontweight="bold",
        color=PALETTE["ink"],
    )
    save_figure(fig, output, "Figure_1_study_design", registry)


def figure_radc(
    continuous: pd.DataFrame,
    stage: pd.DataFrame,
    hc3: pd.DataFrame,
    output: Path,
    registry: list[dict[str, object]],
) -> None:
    labels = {
        "standard": "Primary standard voom",
        "quality_weights": "Quality weights",
        "min_cells_50": "≥50 nuclei",
        "exclude_age89plus": "Exclude 89+ top-code",
        "omit_nonad": "Omit non-AD diagnosis",
        "omit_log_ncells": "Omit nucleus-depth term",
    }
    order = list(labels)
    fig = plt.figure(figsize=(14.8, 7.4))
    grid = fig.add_gridspec(1, 3, width_ratios=[1.65, 1.05, 0.85], wspace=0.30)
    ax_a = fig.add_subplot(grid[0, 0])
    ax_b = fig.add_subplot(grid[0, 1])
    ax_c = fig.add_subplot(grid[0, 2])
    panel_label(ax_a, "A")
    panel_label(ax_b, "B")
    panel_label(ax_c, "C")

    y = np.arange(len(order))[::-1]
    offsets = {"CERAD": -0.13, "Braak": 0.13}
    colors = {"CERAD": PALETTE["blue"], "Braak": PALETTE["ochre"]}
    markers = {"CERAD": "o", "Braak": "s"}
    for axis in ("CERAD", "Braak"):
        subset = (
            continuous.loc[continuous["pathology"].eq(axis)]
            .set_index("model")
            .loc[order]
        )
        yy = y + offsets[axis]
        x = subset["S"].to_numpy()
        low = x - subset["bootstrap_CI95_lower"].to_numpy()
        high = subset["bootstrap_CI95_upper"].to_numpy() - x
        ax_a.errorbar(
            x,
            yy,
            xerr=np.vstack([low, high]),
            fmt=markers[axis],
            color=colors[axis],
            ecolor=colors[axis],
            elinewidth=1.25,
            capsize=2.5,
            markersize=5.8,
            markeredgecolor="white",
            markeredgewidth=0.7,
            label=axis,
            zorder=3,
        )
    ax_a.axvline(0, color=PALETTE["slate"], linewidth=0.9, linestyle="--")
    ax_a.set_yticks(y, [labels[item] for item in order])
    ax_a.set_xlabel("Locus-equal alignment score, S (locus-bootstrap 95% CI)")
    ax_a.set_title("RADC continuous-pathology robustness", loc="left", fontweight="bold")
    ax_a.grid(axis="x", alpha=0.75)
    legend = ax_a.get_legend()
    if legend is not None:
        legend.remove()
    ax_a.text(
        0.72,
        1.01,
        "● CERAD",
        transform=ax_a.transAxes,
        color=PALETTE["blue"],
        fontsize=8.7,
        va="bottom",
    )
    ax_a.text(
        0.87,
        1.01,
        "■ Braak",
        transform=ax_a.transAxes,
        color=PALETTE["ochre"],
        fontsize=8.7,
        va="bottom",
    )

    stage_order = [
        ("CERAD", "middle_vs_low"),
        ("CERAD", "high_vs_low"),
        ("Braak", "middle_vs_low"),
        ("Braak", "high_vs_low"),
    ]
    stage_labels = [
        "CERAD: middle vs low",
        "CERAD: high vs low",
        "Braak: middle vs low",
        "Braak: high vs low",
    ]
    stage_indexed = stage.set_index(["pathology", "contrast"]).loc[stage_order]
    yy = np.arange(4)[::-1]
    for index, ((axis, _), row) in enumerate(stage_indexed.iterrows()):
        ypos = yy[index]
        color = colors[axis]
        ax_b.errorbar(
            row["S"],
            ypos,
            xerr=np.array(
                [
                    [row["S"] - row["bootstrap_CI95_lower"]],
                    [row["bootstrap_CI95_upper"] - row["S"]],
                ]
            ),
            fmt=markers[axis],
            color=color,
            ecolor=color,
            capsize=2.5,
            markersize=6,
            markeredgecolor="white",
            markeredgewidth=0.7,
        )
        ax_b.text(
            0.39,
            ypos,
            f"P={row['empirical_p_raw']:.3f}",
            va="center",
            fontsize=7.8,
            color=PALETTE["slate"],
        )
    ax_b.axvline(0, color=PALETTE["slate"], linewidth=0.9, linestyle="--")
    ax_b.set_yticks(yy, stage_labels)
    ax_b.set_xlim(-0.23, 0.53)
    ax_b.set_xlabel("S (locus-bootstrap 95% CI)")
    ax_b.set_title("Frozen stage contrasts", loc="left", fontweight="bold")
    ax_b.grid(axis="x", alpha=0.75)

    primary = continuous.loc[continuous["model"].eq("standard")].set_index("pathology")
    hc3_index = hc3.set_index("pathology")
    ax_c.set_xlim(0, 1)
    ax_c.set_ylim(0, 1)
    ax_c.axis("off")
    ax_c.set_title("Multiplicity and robust SE", loc="left", fontweight="bold")
    y0 = 0.86
    for index, axis in enumerate(("CERAD", "Braak")):
        yloc = y0 - index * 0.28
        p = primary.loc[axis, "empirical_p_raw"]
        holm = primary.loc[axis, "empirical_p_Holm_primary_family"]
        ax_c.text(0.04, yloc, axis, fontweight="bold", color=colors[axis])
        ax_c.text(0.04, yloc - 0.07, f"Primary S  {primary.loc[axis, 'S']:.3f}", fontsize=8.5)
        ax_c.text(0.04, yloc - 0.13, f"HC3 S      {hc3_index.loc[axis, 'S_HC3']:.3f}", fontsize=8.5)
        ax_c.text(0.04, yloc - 0.19, f"Raw P {p:.4f}  |  Holm {holm:.4f}", fontsize=8.5)
    ax_c.add_patch(
        FancyBboxPatch(
            (0.03, 0.07),
            0.90,
            0.18,
            boxstyle="round,pad=0.012,rounding_size=0.02",
            facecolor=PALETTE["pale"],
            edgecolor=PALETTE["grid"],
        )
    )
    ax_c.text(
        0.07,
        0.20,
        "Result",
        fontweight="bold",
        color=PALETTE["navy"],
    )
    ax_c.text(
        0.07,
        0.14,
        "Positive direction is retained,\n"
        "but neither co-primary test passes Holm 0.05.",
        fontsize=8.2,
        va="top",
        linespacing=1.3,
    )
    fig.suptitle(
        "RADC: a directionally stable but multiplicity-limited Braak signal",
        x=0.07,
        ha="left",
        fontsize=15,
        fontweight="bold",
    )
    save_figure(fig, output, "Figure_2_RADC_robustness", registry)


def figure_seaad(
    seaad: pd.DataFrame, output: Path, registry: list[dict[str, object]]
) -> None:
    row_order = [
        ("region", "DFC", "standard", "DFC · standard"),
        ("region", "DFC", "composition", "DFC · composition"),
        ("region", "MEC", "standard", "MEC · standard"),
        ("region", "MEC", "composition", "MEC · composition"),
        ("region", "MTG", "standard", "MTG · standard"),
        ("region", "MTG", "composition", "MTG · composition"),
        ("repeated", "DFC_MEC_MTG", "standard", "Repeated · donor block"),
    ]
    cmap = LinearSegmentedColormap.from_list(
        "restrained_diverging",
        [PALETTE["burgundy"], "#E6DDDA", PALETTE["pale"], "#D8E1E5", PALETTE["navy"]],
    )
    norm = Normalize(vmin=-0.40, vmax=0.40)
    fig, axes = plt.subplots(1, 3, figsize=(13.8, 7.5), sharey=True)
    for index, cell_class in enumerate(("Immune", "Oligo", "EN")):
        ax = axes[index]
        panel_label(ax, chr(ord("A") + index))
        matrix = np.full((len(row_order), 2), np.nan)
        for row_index, (analysis, region, adjustment, _) in enumerate(row_order):
            for col_index, axis in enumerate(("CERAD", "Braak")):
                selected = seaad.loc[
                    seaad["cell_class"].eq(cell_class)
                    & seaad["pathology"].eq(axis)
                    & seaad["analysis"].eq(analysis)
                    & seaad["region"].eq(region)
                    & seaad["adjustment"].eq(adjustment),
                    "S",
                ]
                if len(selected) != 1:
                    raise RuntimeError("SEA-AD heatmap cell mismatch")
                matrix[row_index, col_index] = float(selected.iloc[0])
        image = ax.imshow(matrix, cmap=cmap, norm=norm, aspect="auto")
        for row_index in range(matrix.shape[0]):
            for col_index in range(matrix.shape[1]):
                value = matrix[row_index, col_index]
                color = "white" if abs(value) > 0.27 else PALETTE["ink"]
                suffix = "†" if cell_class == "Oligo" else ""
                ax.text(
                    col_index,
                    row_index,
                    f"{value:+.2f}{suffix}",
                    ha="center",
                    va="center",
                    fontsize=8.4,
                    color=color,
                    fontweight="bold" if row_index == 0 else "normal",
                )
        ax.set_xticks([0, 1], ["CERAD", "Braak"])
        ax.set_title(
            {"Immune": "Immune", "Oligo": "Oligodendrocyte", "EN": "Excitatory neuron"}[
                cell_class
            ],
            fontweight="bold",
        )
        ax.tick_params(length=0)
        for spine in ax.spines.values():
            spine.set_visible(False)
    axes[0].set_yticks(
        np.arange(len(row_order)), [item[3] for item in row_order], fontsize=8.3
    )
    fig.subplots_adjust(
        left=0.08, right=0.90, top=0.88, bottom=0.14, wspace=0.20
    )
    colorbar_axis = fig.add_axes([0.925, 0.22, 0.014, 0.60])
    colorbar = fig.colorbar(image, cax=colorbar_axis)
    colorbar.set_label("Locus-equal alignment score, S")
    fig.text(
        0.70,
        0.035,
        "† 9 loci; low-locus exploratory evidence",
        fontsize=8,
        color=PALETTE["slate"],
    )
    fig.suptitle(
        "SEA-AD: direction estimates vary by cell class, region, and composition adjustment",
        x=0.06,
        ha="left",
        fontsize=14.5,
        fontweight="bold",
    )
    save_figure(
        fig,
        output,
        "Figure_3_SEAAD_region_composition",
        registry,
        adjust_layout=False,
    )


def figure_direct(
    direct: pd.DataFrame,
    regional: pd.DataFrame,
    output: Path,
    registry: list[dict[str, object]],
) -> None:
    order = [
        ("Immune", "CERAD"),
        ("Immune", "Braak"),
        ("Oligo", "CERAD"),
        ("Oligo", "Braak"),
        ("EN", "CERAD"),
        ("EN", "Braak"),
    ]
    labels = [
        "Immune · CERAD",
        "Immune · Braak",
        "Oligo · CERAD",
        "Oligo · Braak",
        "Excitatory · CERAD",
        "Excitatory · Braak",
    ]
    indexed = direct.set_index(["cell_class", "pathology"]).loc[order]
    fig = plt.figure(figsize=(15.2, 7.3))
    grid = fig.add_gridspec(1, 3, width_ratios=[1.05, 1.10, 1.35], wspace=0.34)
    ax_a = fig.add_subplot(grid[0, 0])
    ax_b = fig.add_subplot(grid[0, 1])
    ax_c = fig.add_subplot(grid[0, 2])
    panel_label(ax_a, "A")
    panel_label(ax_b, "B")
    panel_label(ax_c, "C")
    y = np.arange(len(order))[::-1]
    for index, (_, row) in enumerate(indexed.iterrows()):
        ypos = y[index]
        ax_a.plot(
            [row["S_RADC_common"], row["S_SEAAD_DFC_common"]],
            [ypos, ypos],
            color=PALETTE["grid"],
            linewidth=1.6,
            zorder=1,
        )
        ax_a.scatter(
            row["S_RADC_common"],
            ypos,
            s=48,
            color=PALETTE["navy"],
            edgecolor="white",
            linewidth=0.7,
            zorder=2,
        )
        ax_a.scatter(
            row["S_SEAAD_DFC_common"],
            ypos,
            s=48,
            color=PALETTE["slate"],
            edgecolor="white",
            linewidth=0.7,
            zorder=2,
        )
    ax_a.axvline(0, color=PALETTE["slate"], linestyle="--", linewidth=0.8)
    ax_a.set_yticks(y, labels)
    ax_a.set_xlabel("S among common loci")
    ax_a.set_title("Paired cohort estimates", loc="left", fontweight="bold")
    ax_a.grid(axis="x")
    ax_a.scatter([], [], color=PALETTE["navy"], label="RADC")
    ax_a.scatter([], [], color=PALETTE["slate"], label="SEA-AD DFC")
    ax_a.text(
        0.62,
        1.01,
        "● RADC",
        transform=ax_a.transAxes,
        color=PALETTE["navy"],
        fontsize=8.5,
        va="bottom",
    )
    ax_a.text(
        0.78,
        1.01,
        "● SEA-AD DFC",
        transform=ax_a.transAxes,
        color=PALETTE["slate"],
        fontsize=8.5,
        va="bottom",
    )

    for index, (_, row) in enumerate(indexed.iterrows()):
        ypos = y[index]
        low = row["Delta_S_RADC_minus_SEAAD_DFC"] - row["bootstrap_CI95_lower"]
        high = row["bootstrap_CI95_upper"] - row["Delta_S_RADC_minus_SEAAD_DFC"]
        low_locus = bool(row["low_locus_count_exploratory"])
        ax_b.errorbar(
            row["Delta_S_RADC_minus_SEAAD_DFC"],
            ypos,
            xerr=np.array([[low], [high]]),
            fmt="D" if low_locus else "o",
            color=PALETTE["ochre"] if low_locus else PALETTE["blue"],
            ecolor=PALETTE["ochre"] if low_locus else PALETTE["blue"],
            capsize=2.5,
            markersize=5.5,
            markeredgecolor="white",
            markeredgewidth=0.7,
        )
    ax_b.axvline(0, color=PALETTE["slate"], linestyle="--", linewidth=0.8)
    ax_b.set_yticks(y, [])
    ax_b.set_xlabel("RADC − SEA-AD DFC, ΔS (95% CI)")
    ax_b.set_title("Direct paired differences", loc="left", fontweight="bold")
    ax_b.grid(axis="x")

    standard = regional.loc[regional["adjustment"].eq("standard")].copy()
    regional_order = order
    colors = {"MEC": PALETTE["blue"], "MTG": PALETTE["burgundy"]}
    offsets = {"MEC": 0.12, "MTG": -0.12}
    for region in ("MEC", "MTG"):
        subset = (
            standard.loc[standard["other_region"].eq(region)]
            .set_index(["cell_class", "pathology"])
            .loc[regional_order]
        )
        yy = y + offsets[region]
        x = subset["Delta_S_DFC_minus_other"].to_numpy(float)
        low = x - subset["bootstrap_CI95_lower"].to_numpy(float)
        high = subset["bootstrap_CI95_upper"].to_numpy(float) - x
        ax_c.errorbar(
            x,
            yy,
            xerr=np.vstack([low, high]),
            fmt="o" if region == "MEC" else "s",
            color=colors[region],
            ecolor=colors[region],
            capsize=2.2,
            markersize=5,
            markeredgecolor="white",
            markeredgewidth=0.7,
            label=f"DFC − {region}",
        )
    ax_c.axvline(0, color=PALETTE["slate"], linestyle="--", linewidth=0.8)
    ax_c.set_yticks(y, labels)
    ax_c.yaxis.tick_right()
    ax_c.set_xlabel("Within-SEA-AD regional ΔS (95% CI)")
    ax_c.set_title("Regional context", loc="left", fontweight="bold")
    ax_c.grid(axis="x")
    legend = ax_c.get_legend()
    if legend is not None:
        legend.remove()
    ax_c.text(
        0.52,
        1.01,
        "● DFC−MEC",
        transform=ax_c.transAxes,
        color=PALETTE["blue"],
        fontsize=8.5,
        va="bottom",
    )
    ax_c.text(
        0.75,
        1.01,
        "■ DFC−MTG",
        transform=ax_c.transAxes,
        color=PALETTE["burgundy"],
        fontsize=8.5,
        va="bottom",
    )
    fig.text(
        0.36,
        0.035,
        "◆ fewer than 10 loci; no RADC−SEA-AD contrast passes BH 0.05.",
        fontsize=8,
        color=PALETTE["slate"],
        ha="center",
    )
    fig.text(
        0.78,
        0.035,
        "Oligodendrocyte DFC−MTG differences pass BH only without composition PCs.",
        fontsize=8,
        color=PALETTE["slate"],
        ha="center",
    )
    fig.suptitle(
        "Direct comparisons localize analysis-context dependence without establishing mechanism",
        x=0.065,
        ha="left",
        fontsize=14.5,
        fontweight="bold",
    )
    save_figure(fig, output, "Figure_4_direct_context_contrasts", registry)


def supplementary_composition(
    seaad: pd.DataFrame, output: Path, registry: list[dict[str, object]]
) -> None:
    region = seaad.loc[seaad["analysis"].eq("region")].copy()
    standard = region.loc[region["adjustment"].eq("standard")].set_index(
        ["cell_class", "pathology", "region"]
    )
    composition = region.loc[region["adjustment"].eq("composition")].set_index(
        ["cell_class", "pathology", "region"]
    )
    paired = standard[["S"]].join(
        composition[["S"]], lsuffix="_standard", rsuffix="_composition"
    )
    fig, ax = plt.subplots(figsize=(7.2, 6.5))
    panel_label(ax, "A")
    class_colors = {
        "Immune": PALETTE["navy"],
        "Oligo": PALETTE["ochre"],
        "EN": PALETTE["burgundy"],
    }
    region_markers = {"DFC": "o", "MEC": "s", "MTG": "^"}
    for (cell_class, pathology, brain_region), row in paired.iterrows():
        ax.scatter(
            row["S_standard"],
            row["S_composition"],
            s=58,
            color=class_colors[cell_class],
            marker=region_markers[brain_region],
            edgecolor="white",
            linewidth=0.8,
            alpha=0.92 if pathology == "Braak" else 0.65,
        )
    lim = (-0.48, 0.48)
    ax.plot(lim, lim, color=PALETTE["slate"], linestyle="--", linewidth=0.9)
    ax.axhline(0, color=PALETTE["grid"], linewidth=0.7)
    ax.axvline(0, color=PALETTE["grid"], linewidth=0.7)
    ax.set_xlim(lim)
    ax.set_ylim(lim)
    ax.set_aspect("equal")
    ax.set_xlabel("Standard regional model, S")
    ax.set_ylabel("Composition-adjusted regional model, S")
    ax.set_title(
        "SEA-AD composition sensitivity across fixed class × region × pathology tests",
        loc="left",
        fontweight="bold",
    )
    cell_handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            color="none",
            markerfacecolor=color,
            markeredgecolor="white",
            markersize=7,
            label=cell_class,
        )
        for cell_class, color in class_colors.items()
    ]
    first = ax.legend(
        handles=cell_handles,
        frameon=False,
        loc="upper left",
        title="Cell class",
    )
    ax.add_artist(first)
    region_handles = [
        Line2D(
            [0],
            [0],
            marker=marker,
            color="none",
            markerfacecolor=PALETTE["slate"],
            markeredgecolor="white",
            markersize=7,
            label=brain_region,
        )
        for brain_region, marker in region_markers.items()
    ]
    ax.legend(
        handles=region_handles,
        frameon=False,
        loc="lower right",
        title="Region",
    )
    ax.grid(alpha=0.5)
    save_figure(fig, output, "Supplementary_Figure_S1_composition_shift", registry)


def main() -> None:
    args = parse_args()
    project = args.project.resolve()
    output = args.output_dir.resolve()
    set_style()
    synthesis = project / "outputs/phase6/presubmission_synthesis_v2"
    seaad_summary = (
        project / "outputs/phase6/seaad_region_sensitivity_v2/summary"
    )
    continuous = pd.read_csv(synthesis / "radc_continuous_sensitivities.csv")
    stage = pd.read_csv(synthesis / "radc_stage_contrasts.csv")
    hc3 = pd.read_csv(synthesis / "radc_hc3_sensitivity.csv")
    direct = pd.read_csv(synthesis / "radc_seaad_dfc_direct_contrasts.csv")
    seaad = pd.read_csv(seaad_summary / "seaad_fixed_42_results.csv")
    regional = pd.read_csv(
        seaad_summary / "seaad_regional_paired_contrasts.csv"
    )
    registry: list[dict[str, object]] = []
    figure_design(output, registry)
    figure_radc(continuous, stage, hc3, output, registry)
    figure_seaad(seaad, output, registry)
    figure_direct(direct, regional, output, registry)
    supplementary_composition(seaad, output, registry)
    manifest = {
        "status": "COMPLETE",
        "style": "restrained academic; low saturation; redundant shape encoding",
        "palette": PALETTE,
        "accessibility": [
            "No rainbow scale",
            "Color and shape jointly encode major groups",
            "Zero-reference lines are explicit",
            "Text labels carry all primary numeric conclusions",
        ],
        "outputs": registry,
    }
    (output / "figure_manifest.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    print(
        f"COMPLETE figures={len(set(item['figure'] for item in registry))} "
        f"files={len(registry)} output={output}"
    )


if __name__ == "__main__":
    main()
