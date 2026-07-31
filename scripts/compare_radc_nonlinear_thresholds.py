#!/usr/bin/env python3
"""Compare RADC nonlinear/stage pathology results at 20 and 50 nuclei."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats


CLASSES = ["Astro", "EN", "Endo", "IN", "Immune", "Mural", "OPC", "Oligo"]
PHENOTYPES = ["CERAD", "Braak"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dir20", type=Path, required=True)
    parser.add_argument("--dir50", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    manifest20 = pd.read_csv(args.dir20 / "nonlinear_model_manifest.csv").set_index(
        ["class", "phenotype"]
    )
    manifest50 = pd.read_csv(args.dir50 / "nonlinear_model_manifest.csv").set_index(
        ["class", "phenotype"]
    )
    rows: list[dict] = []
    stable_frames: list[pd.DataFrame] = []
    quadratic_frames: list[pd.DataFrame] = []
    for class_name in CLASSES:
        for phenotype in PHENOTYPES:
            name = f"{class_name}__{phenotype}__nonlinear_results.csv"
            columns = [
                "gene_id",
                "gene_name",
                "high_vs_low_beta",
                "high_vs_low_P",
                "high_vs_low_FDR",
                "stage_omnibus_FDR",
                "quadratic_FDR",
            ]
            a = pd.read_csv(args.dir20 / name, usecols=columns).rename(
                columns={column: f"{column}_20" for column in columns if column != "gene_id"}
            )
            b = pd.read_csv(args.dir50 / name, usecols=columns).rename(
                columns={column: f"{column}_50" for column in columns if column != "gene_id"}
            )
            merged = a.merge(b, on="gene_id", validate="one_to_one")
            rho = stats.spearmanr(
                merged["high_vs_low_beta_20"], merged["high_vs_low_beta_50"]
            )[0]
            prioritized = merged[
                merged[["high_vs_low_P_20", "high_vs_low_P_50"]].min(axis=1) < 0.01
            ]
            stage20 = set(
                merged.loc[merged["stage_omnibus_FDR_20"] < 0.05, "gene_id"]
            )
            stage50 = set(
                merged.loc[merged["stage_omnibus_FDR_50"] < 0.05, "gene_id"]
            )
            high20 = set(merged.loc[merged["high_vs_low_FDR_20"] < 0.05, "gene_id"])
            high50 = set(merged.loc[merged["high_vs_low_FDR_50"] < 0.05, "gene_id"])
            quad20 = set(merged.loc[merged["quadratic_FDR_20"] < 0.05, "gene_id"])
            quad50 = set(merged.loc[merged["quadratic_FDR_50"] < 0.05, "gene_id"])
            rows.append(
                {
                    "class": class_name,
                    "phenotype": phenotype,
                    "n20": int(manifest20.loc[(class_name, phenotype), "n_donors"]),
                    "n50": int(manifest50.loc[(class_name, phenotype), "n_donors"]),
                    "common_genes": len(merged),
                    "high_low_beta_spearman": float(rho),
                    "prioritized_genes_union_p01": len(prioritized),
                    "prioritized_sign_concordance": (
                        float(
                            np.mean(
                                np.sign(prioritized["high_vs_low_beta_20"])
                                == np.sign(prioritized["high_vs_low_beta_50"])
                            )
                        )
                        if len(prioritized)
                        else np.nan
                    ),
                    "high_low_fdr20": len(high20),
                    "high_low_fdr50": len(high50),
                    "high_low_fdr_intersection": len(high20 & high50),
                    "stage_omnibus_fdr20": len(stage20),
                    "stage_omnibus_fdr50": len(stage50),
                    "stage_omnibus_fdr_intersection": len(stage20 & stage50),
                    "quadratic_fdr20": len(quad20),
                    "quadratic_fdr50": len(quad50),
                    "quadratic_fdr_intersection": len(quad20 & quad50),
                }
            )
            if stage20 & stage50:
                frame = merged[merged["gene_id"].isin(stage20 & stage50)].copy()
                frame.insert(0, "phenotype", phenotype)
                frame.insert(0, "class", class_name)
                frame["gene_name"] = frame["gene_name_20"].fillna(
                    frame["gene_name_50"]
                )
                stable_frames.append(frame)
            if quad20 & quad50:
                frame = merged[merged["gene_id"].isin(quad20 & quad50)].copy()
                frame.insert(0, "phenotype", phenotype)
                frame.insert(0, "class", class_name)
                quadratic_frames.append(frame)

    comparison = pd.DataFrame(rows)
    comparison.to_csv(
        args.output_dir / "nonlinear_threshold_stability_by_model.csv", index=False
    )
    stable_stage = (
        pd.concat(stable_frames, ignore_index=True) if stable_frames else pd.DataFrame()
    )
    stable_stage.to_csv(
        args.output_dir / "stable_stage_omnibus_fdr05_hits.csv", index=False
    )
    stable_quadratic = (
        pd.concat(quadratic_frames, ignore_index=True)
        if quadratic_frames
        else pd.DataFrame()
    )
    stable_quadratic.to_csv(
        args.output_dir / "stable_quadratic_fdr05_hits.csv", index=False
    )

    sns.set_theme(style="whitegrid", context="notebook")
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.8))
    specifications = [
        ("high_low_beta_spearman", "High-low effect rho", ".2f", "viridis", 0, 1),
        (
            "prioritized_sign_concordance",
            "High-low sign concordance",
            ".2f",
            "viridis",
            0,
            1,
        ),
        (
            "stage_omnibus_fdr_intersection",
            "Stable stage-omnibus FDR hits",
            ".0f",
            "mako",
            None,
            None,
        ),
    ]
    for axis, (column, title, fmt, cmap, vmin, vmax) in zip(axes, specifications):
        matrix = (
            comparison.pivot(index="class", columns="phenotype", values=column)
            .reindex(index=CLASSES, columns=PHENOTYPES)
        )
        sns.heatmap(
            matrix,
            annot=True,
            fmt=fmt,
            cmap=cmap,
            vmin=vmin,
            vmax=vmax,
            cbar=False,
            ax=axis,
        )
        axis.set_title(title)
    fig.suptitle("RADC nonlinear pathology threshold stability", y=1.02)
    fig.tight_layout()
    fig.savefig(
        args.output_dir / "nonlinear_threshold_stability_heatmaps.png",
        dpi=180,
        bbox_inches="tight",
    )
    plt.close(fig)

    stable_models = comparison[
        (comparison["high_low_beta_spearman"] >= 0.90)
        & (comparison["prioritized_sign_concordance"] >= 0.90)
    ]
    report = [
        "# RADC nonlinear-pathology threshold audit",
        "",
        f"- Models compared: {len(comparison)}/16.",
        (
            f"- Models passing high-low effect stability rule: "
            f"{len(stable_models)}/16."
        ),
        (
            f"- Stage-omnibus FDR rows stable at both thresholds: "
            f"{len(stable_stage)}."
        ),
        (
            f"- Quadratic FDR rows stable at both thresholds: "
            f"{len(stable_quadratic)}."
        ),
        "",
        "A high-versus-low contrast is not itself proof of nonlinearity. Stable",
        "quadratic or omnibus evidence and external stage replication are required",
        "before describing threshold, biphasic, adaptive, or protective responses.",
    ]
    (args.output_dir / "NONLINEAR_THRESHOLD_REPORT.md").write_text(
        "\n".join(report) + "\n", encoding="utf-8"
    )
    manifest = {
        "status": "complete",
        "models_compared": len(comparison),
        "models_passing_effect_stability_rule": len(stable_models),
        "stable_stage_omnibus_fdr05_rows": len(stable_stage),
        "stable_quadratic_fdr05_rows": len(stable_quadratic),
    }
    (args.output_dir / "nonlinear_threshold_manifest.json").write_text(
        json.dumps(manifest, indent=2), encoding="utf-8"
    )


if __name__ == "__main__":
    main()
