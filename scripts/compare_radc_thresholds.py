#!/usr/bin/env python3
"""Compare RADC pathology effects at 20- versus 50-nucleus thresholds."""

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
PHENOTYPES = ["CERAD", "Braak", "Dementia", "AD_status"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dir20", type=Path, required=True)
    parser.add_argument("--dir50", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    manifest20 = pd.read_csv(args.dir20 / "model_manifest.csv").set_index(
        ["class", "phenotype"]
    )
    manifest50 = pd.read_csv(args.dir50 / "model_manifest.csv").set_index(
        ["class", "phenotype"]
    )
    rows: list[dict] = []
    stable_hits: list[pd.DataFrame] = []
    for class_name in CLASSES:
        for phenotype in PHENOTYPES:
            file_name = f"{class_name}__{phenotype}__gene_results.csv"
            a = pd.read_csv(args.dir20 / file_name)[
                ["gene_id", "gene_name", "beta_D", "P.Value", "FDR_within_class"]
            ].rename(
                columns={
                    "gene_name": "gene_name_20",
                    "beta_D": "beta_20",
                    "P.Value": "p_20",
                    "FDR_within_class": "fdr_20",
                }
            )
            b = pd.read_csv(args.dir50 / file_name)[
                ["gene_id", "gene_name", "beta_D", "P.Value", "FDR_within_class"]
            ].rename(
                columns={
                    "gene_name": "gene_name_50",
                    "beta_D": "beta_50",
                    "P.Value": "p_50",
                    "FDR_within_class": "fdr_50",
                }
            )
            merged = a.merge(b, on="gene_id", validate="one_to_one")
            spearman = stats.spearmanr(merged["beta_20"], merged["beta_50"])[0]
            pearson = stats.pearsonr(merged["beta_20"], merged["beta_50"])[0]
            union_p01 = merged[merged[["p_20", "p_50"]].min(axis=1) < 0.01]
            union_fdr = merged[
                (merged["fdr_20"] < 0.05) | (merged["fdr_50"] < 0.05)
            ]
            hits20 = set(merged.loc[merged["fdr_20"] < 0.05, "gene_id"])
            hits50 = set(merged.loc[merged["fdr_50"] < 0.05, "gene_id"])
            intersection = hits20.intersection(hits50)
            union = hits20.union(hits50)
            rows.append(
                {
                    "class": class_name,
                    "phenotype": phenotype,
                    "n20": int(manifest20.loc[(class_name, phenotype), "n_donors"]),
                    "n50": int(manifest50.loc[(class_name, phenotype), "n_donors"]),
                    "common_genes": len(merged),
                    "beta_spearman": float(spearman),
                    "beta_pearson": float(pearson),
                    "union_p01_genes": len(union_p01),
                    "sign_concordance_union_p01": (
                        float(
                            np.mean(
                                np.sign(union_p01["beta_20"])
                                == np.sign(union_p01["beta_50"])
                            )
                        )
                        if len(union_p01)
                        else np.nan
                    ),
                    "fdr20": len(hits20),
                    "fdr50": len(hits50),
                    "fdr_intersection": len(intersection),
                    "fdr_union": len(union),
                    "fdr_jaccard": len(intersection) / len(union) if union else np.nan,
                    "sign_concordance_fdr_union": (
                        float(
                            np.mean(
                                np.sign(union_fdr["beta_20"])
                                == np.sign(union_fdr["beta_50"])
                            )
                        )
                        if len(union_fdr)
                        else np.nan
                    ),
                }
            )
            if intersection:
                stable = merged[merged["gene_id"].isin(intersection)].copy()
                stable.insert(0, "phenotype", phenotype)
                stable.insert(0, "class", class_name)
                stable["gene_name"] = stable["gene_name_20"].fillna(
                    stable["gene_name_50"]
                )
                stable_hits.append(stable)

    comparison = pd.DataFrame(rows)
    comparison.to_csv(args.output_dir / "threshold_stability_by_model.csv", index=False)
    if stable_hits:
        stable = pd.concat(stable_hits, ignore_index=True)
    else:
        stable = pd.DataFrame()
    stable.to_csv(args.output_dir / "fdr05_stable_at_both_thresholds.csv", index=False)

    sns.set_theme(style="whitegrid", context="notebook")
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.8))
    matrices = [
        ("beta_spearman", "Effect Spearman rho", ".2f", "viridis", 0, 1),
        (
            "sign_concordance_union_p01",
            "Sign concordance: union P<0.01",
            ".2f",
            "viridis",
            0,
            1,
        ),
        ("n50", "Eligible donors at >=50 nuclei", ".0f", "mako", None, None),
    ]
    for axis, (column, title, fmt, cmap, vmin, vmax) in zip(axes, matrices):
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
            ax=axis,
            cbar=False,
        )
        axis.set_title(title)
    fig.suptitle("RADC 20- versus 50-nucleus threshold stability", y=1.02)
    fig.tight_layout()
    fig.savefig(
        args.output_dir / "threshold_stability_heatmaps.png",
        dpi=180,
        bbox_inches="tight",
    )
    plt.close(fig)

    high_stability = comparison[
        (comparison["beta_spearman"] >= 0.90)
        & (
            comparison["sign_concordance_union_p01"].isna()
            | (comparison["sign_concordance_union_p01"] >= 0.90)
        )
    ]
    report = [
        "# RADC threshold-stability audit",
        "",
        f"- Models compared: {len(comparison)}/32.",
        (
            f"- Models with beta Spearman >=0.90 and prioritized-sign "
            f"concordance >=0.90: {len(high_stability)}/32."
        ),
        f"- Genes significant at FDR<0.05 at both thresholds: {len(stable)}.",
        "",
        "## Models below the effect-stability rule",
        "",
    ]
    unstable = comparison.drop(high_stability.index).sort_values("beta_spearman")
    if len(unstable):
        for _, row in unstable.iterrows():
            report.append(
                f"- {row['class']} {row['phenotype']}: "
                f"n={int(row['n20'])}->{int(row['n50'])}, "
                f"rho={row['beta_spearman']:.3f}, prioritized sign="
                f"{row['sign_concordance_union_p01']:.3f}."
            )
    else:
        report.append("- None.")
    report += [
        "",
        "FDR overlap is reported but is not the primary stability criterion; "
        "continuous effect and sign stability are less sensitive to threshold-crossing.",
    ]
    (args.output_dir / "THRESHOLD_STABILITY_REPORT.md").write_text(
        "\n".join(report) + "\n", encoding="utf-8"
    )
    manifest = {
        "status": "complete",
        "models_compared": len(comparison),
        "models_passing_effect_stability_rule": len(high_stability),
        "stable_fdr05_rows": len(stable),
        "rule": "beta Spearman >=0.90 and prioritized union P<0.01 sign concordance >=0.90",
    }
    (args.output_dir / "threshold_stability_manifest.json").write_text(
        json.dumps(manifest, indent=2), encoding="utf-8"
    )


if __name__ == "__main__":
    main()
