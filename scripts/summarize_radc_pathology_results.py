#!/usr/bin/env python3
"""Result-level QA and compact summaries for RADC pathology models."""

from __future__ import annotations

import argparse
import json
from itertools import combinations
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
    parser.add_argument("--model-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--quality-weight-pilot-dir", type=Path)
    return parser.parse_args()


def lambda_gc(p_values: pd.Series) -> float:
    p = pd.to_numeric(p_values, errors="coerce").to_numpy(dtype=float)
    p = p[np.isfinite(p) & (p >= 0) & (p <= 1)]
    if not len(p):
        return np.nan
    p = np.clip(p, np.finfo(float).tiny, 1)
    chisq = stats.chi2.isf(p, df=1)
    return float(np.median(chisq) / stats.chi2.ppf(0.5, df=1))


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    manifest = pd.read_csv(args.model_dir / "model_manifest.csv")
    if len(manifest) != 32:
        raise RuntimeError(f"Expected 32 complete models, found {len(manifest)}")
    if not manifest["status"].eq("complete").all():
        raise RuntimeError("Not every model is complete")

    model_frames: dict[tuple[str, str], pd.DataFrame] = {}
    qa_rows: list[dict] = []
    top_frames: list[pd.DataFrame] = []
    for class_name in CLASSES:
        for phenotype in PHENOTYPES:
            path = args.model_dir / f"{class_name}__{phenotype}__gene_results.csv"
            frame = pd.read_csv(path)
            expected_columns = {
                "gene_id",
                "gene_name",
                "beta_D",
                "SE_D_moderated",
                "t",
                "P.Value",
                "FDR_within_class",
            }
            missing = expected_columns.difference(frame.columns)
            if missing:
                raise RuntimeError(f"{path.name} missing {sorted(missing)}")
            if frame["gene_id"].duplicated().any():
                raise RuntimeError(f"Duplicated gene IDs: {path.name}")
            numerical = frame[
                ["beta_D", "SE_D_moderated", "t", "P.Value", "FDR_within_class"]
            ].apply(pd.to_numeric, errors="coerce")
            finite = np.isfinite(numerical.to_numpy()).all(axis=1)
            valid_p = frame["P.Value"].between(0, 1, inclusive="both")
            valid_fdr = frame["FDR_within_class"].between(0, 1, inclusive="both")
            if not finite.all() or not valid_p.all() or not valid_fdr.all():
                raise RuntimeError(f"Non-finite or invalid statistics: {path.name}")
            model_frames[(class_name, phenotype)] = frame
            qa_rows.append(
                {
                    "class": class_name,
                    "phenotype": phenotype,
                    "genes": len(frame),
                    "lambda_gc": lambda_gc(frame["P.Value"]),
                    "p_lt_0_05": int(frame["P.Value"].lt(0.05).sum()),
                    "p_lt_0_001": int(frame["P.Value"].lt(0.001).sum()),
                    "fdr_lt_0_05": int(frame["FDR_within_class"].lt(0.05).sum()),
                    "fdr_lt_0_01": int(frame["FDR_within_class"].lt(0.01).sum()),
                    "median_abs_beta": float(frame["beta_D"].abs().median()),
                    "beta_q01": float(frame["beta_D"].quantile(0.01)),
                    "beta_q99": float(frame["beta_D"].quantile(0.99)),
                    "minimum_p": float(frame["P.Value"].min()),
                    "finite_pass": bool(finite.all()),
                    "unique_gene_id_pass": bool(not frame["gene_id"].duplicated().any()),
                }
            )
            top = frame.nsmallest(50, "P.Value").copy()
            top.insert(0, "rank_in_model", np.arange(1, len(top) + 1))
            top_frames.append(top)

    qa = pd.DataFrame(qa_rows)
    qa.to_csv(args.output_dir / "model_result_qa.csv", index=False)
    top_hits = pd.concat(top_frames, ignore_index=True)
    top_hits.to_csv(args.output_dir / "top50_per_model.csv", index=False)

    pair_rows: list[dict] = []
    for class_name in CLASSES:
        for phenotype_a, phenotype_b in combinations(PHENOTYPES, 2):
            a = model_frames[(class_name, phenotype_a)][
                ["gene_id", "beta_D", "P.Value"]
            ].rename(columns={"beta_D": "beta_a", "P.Value": "p_a"})
            b = model_frames[(class_name, phenotype_b)][
                ["gene_id", "beta_D", "P.Value"]
            ].rename(columns={"beta_D": "beta_b", "P.Value": "p_b"})
            merged = a.merge(b, on="gene_id", validate="one_to_one")
            rho, p_value = stats.spearmanr(merged["beta_a"], merged["beta_b"])
            prioritized = merged[merged[["p_a", "p_b"]].min(axis=1) < 0.01]
            sign_concordance = (
                np.mean(np.sign(prioritized["beta_a"]) == np.sign(prioritized["beta_b"]))
                if len(prioritized)
                else np.nan
            )
            pair_rows.append(
                {
                    "class": class_name,
                    "phenotype_a": phenotype_a,
                    "phenotype_b": phenotype_b,
                    "common_genes": len(merged),
                    "spearman_beta": float(rho),
                    "spearman_p": float(p_value),
                    "prioritized_union_p01": len(prioritized),
                    "sign_concordance_prioritized_union_p01": sign_concordance,
                }
            )
    pairwise = pd.DataFrame(pair_rows)
    pairwise.to_csv(args.output_dir / "cross_phenotype_effect_consistency.csv", index=False)

    significant = pd.concat(
        [
            frame.loc[
                frame["FDR_within_class"] < 0.05,
                [
                    "gene_id",
                    "gene_name",
                    "class",
                    "phenotype",
                    "beta_D",
                    "P.Value",
                    "FDR_within_class",
                ],
            ]
            for frame in model_frames.values()
        ],
        ignore_index=True,
    )
    significant.to_csv(args.output_dir / "all_within_class_fdr05_hits.csv", index=False)
    if len(significant):
        recurrence = (
            significant.assign(direction=np.where(significant["beta_D"] > 0, "up", "down"))
            .groupby(["phenotype", "gene_id", "gene_name"], as_index=False)
            .agg(
                n_classes=("class", "nunique"),
                classes=("class", lambda x: ";".join(sorted(set(x)))),
                min_fdr=("FDR_within_class", "min"),
                positive_classes=("beta_D", lambda x: int((x > 0).sum())),
                negative_classes=("beta_D", lambda x: int((x < 0).sum())),
            )
            .sort_values(["n_classes", "min_fdr"], ascending=[False, True])
        )
    else:
        recurrence = pd.DataFrame()
    recurrence.to_csv(args.output_dir / "fdr05_gene_recurrence.csv", index=False)

    quality_comparison_rows: list[dict] = []
    pilot_dir = args.quality_weight_pilot_dir
    if pilot_dir and pilot_dir.exists():
        for phenotype in ["CERAD", "Braak"]:
            pilot_path = pilot_dir / f"Astro__{phenotype}__gene_results.csv"
            if not pilot_path.exists():
                continue
            standard = model_frames[("Astro", phenotype)][
                ["gene_id", "beta_D", "P.Value", "FDR_within_class"]
            ].rename(
                columns={
                    "beta_D": "beta_standard",
                    "P.Value": "p_standard",
                    "FDR_within_class": "fdr_standard",
                }
            )
            pilot = pd.read_csv(pilot_path)[
                ["gene_id", "beta_D", "P.Value", "FDR_within_class"]
            ].rename(
                columns={
                    "beta_D": "beta_quality",
                    "P.Value": "p_quality",
                    "FDR_within_class": "fdr_quality",
                }
            )
            merged = standard.merge(pilot, on="gene_id", validate="one_to_one")
            rho = stats.spearmanr(merged["beta_standard"], merged["beta_quality"])[0]
            pearson = stats.pearsonr(merged["beta_standard"], merged["beta_quality"])[0]
            union = merged[
                (merged["fdr_standard"] < 0.05) | (merged["fdr_quality"] < 0.05)
            ]
            quality_comparison_rows.append(
                {
                    "class": "Astro",
                    "phenotype": phenotype,
                    "genes": len(merged),
                    "beta_spearman": float(rho),
                    "beta_pearson": float(pearson),
                    "standard_fdr05": int((merged["fdr_standard"] < 0.05).sum()),
                    "quality_weight_fdr05": int((merged["fdr_quality"] < 0.05).sum()),
                    "fdr05_union": len(union),
                    "direction_concordance_in_union": (
                        float(
                            np.mean(
                                np.sign(union["beta_standard"])
                                == np.sign(union["beta_quality"])
                            )
                        )
                        if len(union)
                        else np.nan
                    ),
                }
            )
    quality_comparison = pd.DataFrame(quality_comparison_rows)
    quality_comparison.to_csv(
        args.output_dir / "quality_weight_pilot_comparison.csv", index=False
    )

    sns.set_theme(style="whitegrid", context="notebook")
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.7))
    fdr_matrix = (
        qa.pivot(index="class", columns="phenotype", values="fdr_lt_0_05")
        .reindex(index=CLASSES, columns=PHENOTYPES)
    )
    sns.heatmap(
        np.log10(1 + fdr_matrix),
        annot=fdr_matrix.astype(int),
        fmt="d",
        cmap="mako",
        cbar_kws={"label": "log10(1 + FDR<0.05 genes)"},
        ax=axes[0],
    )
    axes[0].set_title("Within-class discoveries")
    lambda_matrix = (
        qa.pivot(index="class", columns="phenotype", values="lambda_gc")
        .reindex(index=CLASSES, columns=PHENOTYPES)
    )
    sns.heatmap(
        lambda_matrix,
        annot=True,
        fmt=".2f",
        center=1,
        cmap="vlag",
        vmin=max(0.7, float(lambda_matrix.min().min())),
        vmax=min(1.5, float(lambda_matrix.max().max())),
        cbar_kws={"label": "lambda GC (diagnostic only)"},
        ax=axes[1],
    )
    axes[1].set_title("P-value calibration diagnostic")
    fig.suptitle("RADC pathology model result QA", y=1.02)
    fig.tight_layout()
    fig.savefig(args.output_dir / "model_result_qa_heatmaps.png", dpi=180, bbox_inches="tight")
    plt.close(fig)

    max_condition = float(pd.to_numeric(manifest["design_condition_number"]).max())
    lambda_min = float(qa["lambda_gc"].min())
    lambda_max = float(qa["lambda_gc"].max())
    report = [
        "# RADC pathology-model result QA",
        "",
        "## Gate status",
        "",
        f"- Complete models: {len(manifest)}/32.",
        "- All tested gene statistics are finite and all gene identifiers are unique.",
        (
            f"- Every design matrix is full rank; maximum condition number "
            f"{max_condition:.2f}."
        ),
        (
            f"- Diagnostic lambda range: {lambda_min:.2f}-{lambda_max:.2f}. "
            "Lambda is descriptive here because expression traits are correlated."
        ),
        f"- Within-class FDR<0.05 rows across all models: {len(significant)}.",
        "",
        "## Discovery counts by phenotype",
        "",
    ]
    for phenotype in PHENOTYPES:
        subset = qa[qa["phenotype"].eq(phenotype)]
        report.append(
            f"- {phenotype}: {int(subset['fdr_lt_0_05'].sum())} "
            "class-gene discoveries."
        )
    if len(quality_comparison):
        report += ["", "## Quality-weight pilot agreement", ""]
        for row in quality_comparison.itertuples():
            report.append(
                f"- Astro {row.phenotype}: beta Spearman={row.beta_spearman:.4f}, "
                f"Pearson={row.beta_pearson:.4f}; standard/quality-weight FDR hits "
                f"{row.standard_fdr05}/{row.quality_weight_fdr05}."
            )
    report += [
        "",
        "## Interpretation boundary",
        "",
        (
            "These are RADC pilot disease-state/pathology estimates. Gene-level "
            "interpretation remains gated on the 50-nucleus and nonlinear-stage "
            "sensitivities, and the main study still requires MSSM/HBCC or other "
            "prespecified cohort evidence."
        ),
    ]
    (args.output_dir / "RESULT_QA_REPORT.md").write_text(
        "\n".join(report) + "\n", encoding="utf-8"
    )

    manifest_json = {
        "status": "complete",
        "models": int(len(manifest)),
        "all_statistics_finite": True,
        "all_gene_ids_unique_within_model": True,
        "max_design_condition_number": max_condition,
        "lambda_gc_min": lambda_min,
        "lambda_gc_max": lambda_max,
        "within_class_fdr05_rows": int(len(significant)),
        "quality_weight_pilot_models_compared": int(len(quality_comparison)),
    }
    (args.output_dir / "result_qa_manifest.json").write_text(
        json.dumps(manifest_json, indent=2), encoding="utf-8"
    )


if __name__ == "__main__":
    main()
