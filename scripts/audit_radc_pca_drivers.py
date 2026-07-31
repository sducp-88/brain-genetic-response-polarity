#!/usr/bin/env python3
"""Audit technical and biological drivers of RADC pseudobulk PCA structure."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats
from sklearn.mixture import GaussianMixture


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--qc-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def bh_adjust(values: pd.Series) -> pd.Series:
    p = pd.to_numeric(values, errors="coerce").to_numpy(dtype=float)
    out = np.full(p.shape, np.nan)
    valid = np.isfinite(p)
    pv = p[valid]
    if not len(pv):
        return pd.Series(out, index=values.index)
    order = np.argsort(pv)
    ranked = pv[order]
    adjusted = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.clip(adjusted, 0, 1)
    restored = np.empty_like(adjusted)
    restored[order] = adjusted
    out[valid] = restored
    return pd.Series(out, index=values.index)


def eta_squared(groups: list[np.ndarray]) -> float:
    all_values = np.concatenate(groups)
    grand = np.mean(all_values)
    ss_between = sum(len(group) * (np.mean(group) - grand) ** 2 for group in groups)
    ss_total = np.sum((all_values - grand) ** 2)
    return float(ss_between / ss_total) if ss_total > 0 else 0.0


def cramers_v(table: pd.DataFrame) -> float:
    if table.shape[0] < 2 or table.shape[1] < 2:
        return math.nan
    chi2 = stats.chi2_contingency(table, correction=False)[0]
    n = table.to_numpy().sum()
    denom = n * min(table.shape[0] - 1, table.shape[1] - 1)
    return float(np.sqrt(chi2 / denom)) if denom > 0 else math.nan


def safe_numeric(frame: pd.DataFrame, column: str) -> pd.Series:
    return pd.to_numeric(frame[column], errors="coerce")


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    pca = pd.read_csv(args.qc_dir / "pca_scores.csv")
    qc = pd.read_csv(args.qc_dir / "sample_qc.csv", low_memory=False)
    metadata_columns = [
        "sample_id",
        "PMI",
        "CERAD",
        "Braak",
        "MCI",
        "Dementia",
        "ApoE_gt",
        "Diagnosis",
        "Tier1_crossDis",
        "Tier2_AD",
        "Tier3_NPS",
        "qc_outlier_library_or_genes",
    ]
    metadata_columns = [column for column in metadata_columns if column in qc.columns]
    data = pca.merge(qc[metadata_columns], on="sample_id", how="left", validate="one_to_one")
    data["log10_n_cells"] = np.log10(safe_numeric(data, "n_cells"))
    data["log10_library_size"] = np.log10(safe_numeric(data, "library_size"))
    for column in ["age_numeric", "PMI", "CERAD", "Braak", "MCI", "Dementia", "ApoE_gt"]:
        if column in data:
            data[column] = safe_numeric(data, column)

    numeric_variables = [
        "log10_n_cells",
        "log10_library_size",
        "age_numeric",
        "PMI",
        "CERAD",
        "Braak",
        "MCI",
        "Dementia",
        "ApoE_gt",
    ]
    numeric_variables = [column for column in numeric_variables if column in data]
    categorical_variables = [
        "sex",
        "ancestry",
        "AD_status",
        "Diagnosis",
        "Tier1_crossDis",
        "Tier2_AD",
        "Tier3_NPS",
    ]
    categorical_variables = [column for column in categorical_variables if column in data]

    numeric_rows: list[dict] = []
    categorical_rows: list[dict] = []
    for class_name, class_frame in data.groupby("class", sort=True):
        for pc in ["PC1", "PC2"]:
            for variable in numeric_variables:
                subset = class_frame[[pc, variable]].dropna()
                if len(subset) < 20 or subset[variable].nunique() < 3:
                    continue
                rho, p_value = stats.spearmanr(subset[pc], subset[variable])
                numeric_rows.append(
                    {
                        "class": class_name,
                        "PC": pc,
                        "variable": variable,
                        "n": len(subset),
                        "spearman_rho": float(rho),
                        "p_value": float(p_value),
                    }
                )
            for variable in categorical_variables:
                subset = class_frame[[pc, variable]].dropna().copy()
                subset[variable] = subset[variable].astype(str)
                subset = subset[~subset[variable].isin(["nan", "", "Not Evaluated"])]
                counts = subset[variable].value_counts()
                allowed = counts[counts >= 3].index
                subset = subset[subset[variable].isin(allowed)]
                groups = [
                    group[pc].to_numpy()
                    for _, group in subset.groupby(variable, observed=True)
                ]
                if len(groups) < 2 or sum(map(len, groups)) < 20:
                    continue
                f_stat, p_value = stats.f_oneway(*groups)
                categorical_rows.append(
                    {
                        "class": class_name,
                        "PC": pc,
                        "variable": variable,
                        "n": len(subset),
                        "levels": len(groups),
                        "eta_squared": eta_squared(groups),
                        "p_value": float(p_value),
                        "level_counts": ";".join(
                            f"{level}:{count}"
                            for level, count in subset[variable].value_counts().items()
                        ),
                    }
                )

    numeric = pd.DataFrame(numeric_rows)
    categorical = pd.DataFrame(categorical_rows)
    numeric["q_value_global"] = bh_adjust(numeric["p_value"])
    categorical["q_value_global"] = bh_adjust(categorical["p_value"])
    numeric.to_csv(args.output_dir / "pca_numeric_associations.csv", index=False)
    categorical.to_csv(args.output_dir / "pca_categorical_associations.csv", index=False)

    combined_numeric = numeric.rename(columns={"spearman_rho": "signed_effect"}).copy()
    combined_numeric["effect_size"] = combined_numeric["signed_effect"].abs()
    combined_numeric["test"] = "Spearman"
    combined_categorical = categorical.rename(columns={"eta_squared": "effect_size"}).copy()
    combined_categorical["signed_effect"] = np.nan
    combined_categorical["test"] = "ANOVA_eta2"
    keep = [
        "class",
        "PC",
        "variable",
        "n",
        "test",
        "effect_size",
        "signed_effect",
        "p_value",
        "q_value_global",
    ]
    combined = pd.concat(
        [combined_numeric[keep], combined_categorical[keep]], ignore_index=True
    )
    combined = combined.sort_values(
        ["class", "PC", "q_value_global", "effect_size"],
        ascending=[True, True, True, False],
    )
    combined.to_csv(args.output_dir / "pca_all_associations.csv", index=False)
    top = (
        combined.sort_values(
            ["class", "PC", "effect_size", "q_value_global"],
            ascending=[True, True, False, True],
        )
        .groupby(["class", "PC"], as_index=False)
        .head(5)
    )
    top.to_csv(args.output_dir / "pca_top_associations.csv", index=False)

    in_frame = data[data["class"].eq("IN")].copy()
    x = in_frame[["PC2"]].to_numpy()
    bic = {}
    models = {}
    for components in [1, 2, 3]:
        model = GaussianMixture(
            n_components=components,
            covariance_type="full",
            random_state=20260729,
            n_init=50,
        ).fit(x)
        bic[components] = float(model.bic(x))
        models[components] = model
    chosen_components = min(bic, key=bic.get)
    chosen = models[chosen_components]
    raw_labels = chosen.predict(x)
    means = chosen.means_.ravel()
    ordering = {old: rank + 1 for rank, old in enumerate(np.argsort(means))}
    in_frame["PC2_band"] = [
        f"band_{ordering[label]}_of_{chosen_components}" for label in raw_labels
    ]
    in_frame["PC2_band_posterior"] = chosen.predict_proba(x).max(axis=1)
    in_frame.sort_values(["PC2_band", "PC2"]).to_csv(
        args.output_dir / "in_pc2_band_assignments.csv", index=False
    )

    band_rows: list[dict] = []
    if chosen_components == 2:
        levels = sorted(in_frame["PC2_band"].unique())
        low = in_frame[in_frame["PC2_band"].eq(levels[0])]
        high = in_frame[in_frame["PC2_band"].eq(levels[1])]
        for variable in numeric_variables:
            left = low[variable].dropna()
            right = high[variable].dropna()
            if len(left) < 5 or len(right) < 5:
                continue
            u_stat, p_value = stats.mannwhitneyu(left, right, alternative="two-sided")
            rank_biserial = 2 * u_stat / (len(left) * len(right)) - 1
            band_rows.append(
                {
                    "variable": variable,
                    "type": "numeric",
                    "n": len(left) + len(right),
                    "effect_size": abs(float(rank_biserial)),
                    "signed_effect": float(rank_biserial),
                    "p_value": float(p_value),
                    "details": (
                        f"{levels[0]} median={left.median():.4g}; "
                        f"{levels[1]} median={right.median():.4g}"
                    ),
                }
            )
        for variable in categorical_variables:
            subset = in_frame[[variable, "PC2_band"]].dropna().copy()
            subset[variable] = subset[variable].astype(str)
            subset = subset[~subset[variable].isin(["nan", "", "Not Evaluated"])]
            table = pd.crosstab(subset[variable], subset["PC2_band"])
            table = table.loc[table.sum(axis=1) >= 3]
            if table.shape[0] < 2 or table.shape[1] != 2:
                continue
            if table.shape == (2, 2):
                _, p_value = stats.fisher_exact(table.to_numpy())
            else:
                _, p_value, _, _ = stats.chi2_contingency(table)
            band_rows.append(
                {
                    "variable": variable,
                    "type": "categorical",
                    "n": int(table.to_numpy().sum()),
                    "effect_size": cramers_v(table),
                    "signed_effect": np.nan,
                    "p_value": float(p_value),
                    "details": "; ".join(
                        f"{index}=" + "/".join(map(str, row))
                        for index, row in table.iterrows()
                    ),
                }
            )
    band_associations = pd.DataFrame(band_rows)
    if len(band_associations):
        band_associations["q_value"] = bh_adjust(band_associations["p_value"])
        band_associations = band_associations.sort_values(
            ["q_value", "effect_size"], ascending=[True, False]
        )
    band_associations.to_csv(
        args.output_dir / "in_pc2_band_associations.csv", index=False
    )

    sns.set_theme(style="whitegrid", context="notebook")
    fig, axes = plt.subplots(2, 3, figsize=(15, 9))
    sns.scatterplot(
        data=in_frame,
        x="PC1",
        y="PC2",
        hue="PC2_band",
        style="ancestry",
        s=48,
        ax=axes[0, 0],
    )
    axes[0, 0].set_title("IN: PC2 mixture-defined bands")
    sns.scatterplot(
        data=in_frame,
        x="log10_library_size",
        y="PC2",
        hue="PC2_band",
        s=45,
        ax=axes[0, 1],
        legend=False,
    )
    axes[0, 1].set_title("IN PC2 versus library size")
    sns.scatterplot(
        data=in_frame,
        x="log10_n_cells",
        y="PC2",
        hue="PC2_band",
        s=45,
        ax=axes[0, 2],
        legend=False,
    )
    axes[0, 2].set_title("IN PC2 versus nuclei count")
    sns.scatterplot(
        data=in_frame,
        x="PC1",
        y="PC2",
        hue="ancestry",
        s=48,
        ax=axes[1, 0],
    )
    axes[1, 0].set_title("IN: ancestry")
    sns.scatterplot(
        data=in_frame,
        x="PC1",
        y="PC2",
        hue="sex",
        s=48,
        ax=axes[1, 1],
    )
    axes[1, 1].set_title("IN: sex")
    sns.scatterplot(
        data=in_frame,
        x="PC1",
        y="PC2",
        hue="Braak",
        palette="viridis",
        s=48,
        ax=axes[1, 2],
    )
    axes[1, 2].set_title("IN: Braak stage")
    fig.suptitle("RADC inhibitory-neuron PCA driver audit", y=1.01)
    fig.tight_layout()
    fig.savefig(args.output_dir / "in_pca_driver_panels.png", dpi=180, bbox_inches="tight")
    plt.close(fig)

    significant = combined[
        combined["q_value_global"].le(0.05) & combined["effect_size"].ge(0.15)
    ]
    report_lines = [
        "# RADC pseudobulk PCA driver audit",
        "",
        "## Scope",
        "",
        (
            "This audit tests whether donor-level PC1/PC2 structure is associated "
            "with sequencing depth, nuclei count, demographics, ancestry, or "
            "neuropathology. Associations are descriptive QC evidence, not causal results."
        ),
        "",
        "## Mixture audit for inhibitory-neuron PC2",
        "",
        f"- Eligible inhibitory-neuron donors: {len(in_frame)}.",
        "- Gaussian-mixture BIC: "
        + ", ".join(f"{k} components = {v:.1f}" for k, v in sorted(bic.items()))
        + ".",
        f"- BIC-selected component count: {chosen_components}.",
    ]
    for band, group in in_frame.groupby("PC2_band"):
        report_lines.append(
            f"- {band}: n={len(group)}, PC2 median={group['PC2'].median():.2f}, "
            f"minimum posterior={group['PC2_band_posterior'].min():.3f}."
        )
    report_lines += [
        "",
        "## Strongest PC associations",
        "",
    ]
    for (class_name, pc), frame in combined.groupby(["class", "PC"], sort=True):
        best = frame.sort_values(
            ["effect_size", "q_value_global"], ascending=[False, True]
        ).head(3)
        descriptions = []
        for row in best.itertuples():
            direction = (
                f", rho={row.signed_effect:.2f}"
                if np.isfinite(row.signed_effect)
                else ""
            )
            descriptions.append(
                f"{row.variable} ({row.test}, effect={row.effect_size:.2f}"
                f"{direction}, q={row.q_value_global:.3g})"
            )
        report_lines.append(f"- {class_name} {pc}: " + "; ".join(descriptions) + ".")
    report_lines += [
        "",
        "## Prespecified interpretation rule",
        "",
        (
            "- Do not exclude a donor solely because of a PCA or mixture-band label. "
            "Hard exclusion requires a documented sample-integrity failure."
        ),
        (
            "- Technical drivers with reproducible association should enter sensitivity "
            "analysis or be handled by normalization/weights; biological and ancestry "
            "drivers should be represented as covariates where estimable."
        ),
        (
            "- A multimodal PC structure is treated as a modeled source of heterogeneity, "
            "not automatically as an outlier set."
        ),
        "",
        f"Associations meeting q<=0.05 and effect>=0.15: {len(significant)}.",
    ]
    (args.output_dir / "PCA_DRIVER_AUDIT.md").write_text(
        "\n".join(report_lines) + "\n", encoding="utf-8"
    )

    manifest = {
        "status": "complete",
        "input_qc_dir": str(args.qc_dir.resolve()),
        "eligible_sample_class_rows": int(len(data)),
        "classes": sorted(data["class"].unique().tolist()),
        "numeric_tests": int(len(numeric)),
        "categorical_tests": int(len(categorical)),
        "significant_effect_tests_q05_effect015": int(len(significant)),
        "in_pc2_gmm_bic": {str(key): value for key, value in bic.items()},
        "in_pc2_selected_components": int(chosen_components),
        "rule": "No PCA-only exclusions; model or sensitivity-audit identified drivers.",
    }
    (args.output_dir / "pca_driver_manifest.json").write_text(
        json.dumps(manifest, indent=2), encoding="utf-8"
    )


if __name__ == "__main__":
    main()
