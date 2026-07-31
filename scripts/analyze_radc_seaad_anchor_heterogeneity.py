#!/usr/bin/env python3
"""Decompose failed RADC-to-SEA-AD Braak replication at shared anchors."""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats


CLASSES = ("Immune", "Oligo", "EN")
RADC_MODELS = ("standard_voom", "quality_weights")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    return parser.parse_args()


def bh_adjust(values: pd.Series) -> np.ndarray:
    p = values.to_numpy(dtype=float)
    order = np.argsort(p)
    ranked = p[order]
    adjusted_ranked = np.minimum.accumulate(
        (ranked * len(ranked) / np.arange(1, len(ranked) + 1))[::-1]
    )[::-1]
    adjusted = np.empty_like(adjusted_ranked)
    adjusted[order] = np.minimum(adjusted_ranked, 1.0)
    return adjusted


def markdown_table(frame: pd.DataFrame) -> str:
    columns = [str(column) for column in frame.columns]
    rows = [
        [
            str(value).replace("|", r"\|").replace("\n", " ")
            for value in row
        ]
        for row in frame.itertuples(index=False, name=None)
    ]
    return "\n".join(
        [
            "| " + " | ".join(columns) + " |",
            "| " + " | ".join("---" for _ in columns) + " |",
            *["| " + " | ".join(row) + " |" for row in rows],
        ]
    )


def load_anchor(path: Path, prefix: str) -> pd.DataFrame:
    table = pd.read_csv(path, low_memory=False)
    required = {
        "anchor_unit_id",
        "locus_id",
        "gene_symbol",
        "gene_id_clean",
        "cell_class",
        "beta_G",
        "SE_G",
        "evidence_grade",
        "cell_mapping_confidence",
        "beta_D",
        "SE_D",
        "expression_estimable",
    }
    missing = required.difference(table.columns)
    if missing:
        raise RuntimeError(f"Missing columns in {path}: {sorted(missing)}")
    table = table.loc[
        table["expression_estimable"].astype(str).str.lower().eq("true"),
        list(required),
    ].copy()
    if table["anchor_unit_id"].duplicated().any():
        duplicated = table.loc[
            table["anchor_unit_id"].duplicated(False), "anchor_unit_id"
        ].tolist()
        raise RuntimeError(f"Duplicated anchors in {path}: {duplicated[:10]}")
    return table.rename(
        columns={
            "beta_D": f"beta_D_{prefix}",
            "SE_D": f"SE_D_{prefix}",
        }
    )


def safe_correlation(
    x: pd.Series, y: pd.Series, method: str
) -> tuple[float, float]:
    if len(x) < 3 or x.nunique() < 2 or y.nunique() < 2:
        return np.nan, np.nan
    result = (
        stats.pearsonr(x, y)
        if method == "pearson"
        else stats.spearmanr(x, y)
    )
    return float(result.statistic), float(result.pvalue)


def summarize_pairwise(frame: pd.DataFrame) -> dict[str, object]:
    pearson_r, pearson_p = safe_correlation(
        frame["beta_D_radc"], frame["beta_D_seaad"], "pearson"
    )
    spearman_rho, spearman_p = safe_correlation(
        frame["beta_D_radc"], frame["beta_D_seaad"], "spearman"
    )
    sign_match = (
        np.sign(frame["beta_D_radc"]) == np.sign(frame["beta_D_seaad"])
    )
    sign_successes = int(sign_match.sum())
    sign_test = stats.binomtest(
        sign_successes, n=len(frame), p=0.5, alternative="two-sided"
    )
    weights = 1.0 / np.square(frame["SE_delta"])
    weighted_delta = float(
        np.sum(weights * frame["delta_seaad_minus_radc"]) / np.sum(weights)
    )
    weighted_delta_se = float(np.sqrt(1.0 / np.sum(weights)))
    weighted_delta_z = weighted_delta / weighted_delta_se
    weighted_delta_p = float(2 * stats.norm.sf(abs(weighted_delta_z)))
    return {
        "anchors": int(len(frame)),
        "pearson_r": pearson_r,
        "pearson_p": pearson_p,
        "spearman_rho": spearman_rho,
        "spearman_p": spearman_p,
        "beta_sign_agreement_n": sign_successes,
        "beta_sign_agreement_fraction": float(sign_match.mean()),
        "beta_sign_agreement_binomial_p": float(sign_test.pvalue),
        "weighted_mean_delta_seaad_minus_radc": weighted_delta,
        "weighted_mean_delta_SE": weighted_delta_se,
        "weighted_mean_delta_z": float(weighted_delta_z),
        "weighted_mean_delta_p": weighted_delta_p,
        "median_absolute_delta": float(
            frame["delta_seaad_minus_radc"].abs().median()
        ),
        "anchor_heterogeneity_FDR05": int(
            (frame["heterogeneity_FDR_within_class"] < 0.05).sum()
        ),
    }


def build_scatter(pairwise: pd.DataFrame, output_dir: Path) -> None:
    standard = pairwise[pairwise["radc_model"].eq("standard_voom")]
    values = pd.concat(
        [standard["beta_D_radc"], standard["beta_D_seaad"]],
        ignore_index=True,
    )
    maximum = max(float(np.nanmax(np.abs(values))), 0.05) * 1.15
    fig, axes = plt.subplots(1, 3, figsize=(13.2, 4.3), sharex=True, sharey=True)
    colors = {"Immune": "#2A9D8F", "Oligo": "#E9C46A", "EN": "#E76F51"}
    for axis, cell_class in zip(axes, CLASSES, strict=True):
        subset = standard[standard["cell_class"].eq(cell_class)].copy()
        axis.axhline(0, color="#B5B5B5", linewidth=0.8)
        axis.axvline(0, color="#B5B5B5", linewidth=0.8)
        axis.plot(
            [-maximum, maximum],
            [-maximum, maximum],
            linestyle="--",
            linewidth=1,
            color="#666666",
        )
        axis.errorbar(
            subset["beta_D_radc"],
            subset["beta_D_seaad"],
            xerr=subset["SE_D_radc"],
            yerr=subset["SE_D_seaad"],
            fmt="o",
            markersize=4.5,
            alpha=0.78,
            color=colors[cell_class],
            ecolor=colors[cell_class],
            elinewidth=0.7,
            capsize=0,
        )
        labels = subset.nlargest(
            min(3, len(subset)), "absolute_standardized_difference"
        )
        for label_index, row in enumerate(labels.itertuples(index=False)):
            axis.annotate(
                str(row.gene_symbol),
                (row.beta_D_radc, row.beta_D_seaad),
                xytext=(4, 5 + 9 * label_index),
                textcoords="offset points",
                fontsize=6.5,
            )
        axis.set_title(f"{cell_class} (n={len(subset)})")
        axis.set_xlim(-maximum, maximum)
        axis.set_ylim(-maximum, maximum)
        axis.set_aspect("equal", adjustable="box")
        axis.grid(alpha=0.15)
    fig.supxlabel("RADC Braak effect per 1 SD (standard voom)")
    fig.supylabel("SEA-AD Braak effect per 1 SD (standard voom)")
    fig.suptitle(
        "Shared genetic-anchor genes show limited cross-cohort portability",
        y=0.98,
    )
    fig.tight_layout(rect=(0.045, 0.045, 1, 0.93))
    fig.savefig(
        output_dir / "radc_seaad_anchor_effect_scatter.png",
        dpi=300,
        bbox_inches="tight",
    )
    fig.savefig(
        output_dir / "radc_seaad_anchor_effect_scatter.pdf",
        bbox_inches="tight",
    )
    plt.close(fig)


def main() -> None:
    args = parse_args()
    project = args.project.resolve()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    radc_paths = {
        "standard_voom": (
            project
            / "outputs/phase4/radc_stage_coupling_standard_voom/Braak/input"
            / "radc_Braak_G1_G2_anchor_pathology_effects.csv"
        ),
        "quality_weights": (
            project
            / "outputs/phase4/radc_stage_coupling_quality_weights/Braak/input"
            / "radc_Braak_G1_G2_anchor_pathology_effects.csv"
        ),
    }
    radc_coupling_paths = {
        "standard_voom": (
            project
            / "outputs/phase4/radc_stage_coupling_standard_voom/Braak/coupling"
            / "direction_coupling_group_results.csv"
        ),
        "quality_weights": (
            project
            / "outputs/phase4/radc_stage_coupling_quality_weights/Braak/coupling"
            / "direction_coupling_group_results.csv"
        ),
    }
    seaad_paths = {
        cell_class: (
            project
            / f"outputs/phase5/seaad_{cell_class.lower()}_core3/Braak/input"
            / f"seaad_Braak_{cell_class}_G1_G2_anchor_pathology_effects.csv"
        )
        for cell_class in CLASSES
    }

    radc_tables = {
        model: load_anchor(path, "radc")
        for model, path in radc_paths.items()
    }
    seaad_tables = {
        cell_class: load_anchor(path, "seaad")
        for cell_class, path in seaad_paths.items()
    }
    pairwise_parts: list[pd.DataFrame] = []
    summary_records: list[dict[str, object]] = []
    coverage_records: list[dict[str, object]] = []
    unmatched_records: list[dict[str, object]] = []
    for radc_model in RADC_MODELS:
        for cell_class in CLASSES:
            radc = radc_tables[radc_model]
            radc = radc[radc["cell_class"].eq(cell_class)].copy()
            seaad = seaad_tables[cell_class].copy()
            metadata_columns = [
                "anchor_unit_id",
                "locus_id",
                "gene_symbol",
                "gene_id_clean",
                "cell_class",
                "beta_G",
                "SE_G",
                "evidence_grade",
                "cell_mapping_confidence",
            ]
            pair = radc.merge(
                seaad[
                    metadata_columns
                    + ["beta_D_seaad", "SE_D_seaad"]
                ],
                on="anchor_unit_id",
                how="inner",
                validate="one_to_one",
                suffixes=("_radc_meta", "_seaad_meta"),
            )
            if len(pair) < 3:
                raise RuntimeError(
                    f"Fewer than three shared anchors for {radc_model}/"
                    f"{cell_class}: joined={len(pair)}"
                )
            shared_ids = set(pair["anchor_unit_id"])
            radc_only = sorted(
                set(radc["anchor_unit_id"]).difference(shared_ids)
            )
            seaad_only = sorted(
                set(seaad["anchor_unit_id"]).difference(shared_ids)
            )
            coverage_records.append(
                {
                    "radc_model": radc_model,
                    "cell_class": cell_class,
                    "radc_estimable_anchors": len(radc),
                    "seaad_estimable_anchors": len(seaad),
                    "shared_anchors": len(pair),
                    "radc_only_anchors": len(radc_only),
                    "seaad_only_anchors": len(seaad_only),
                }
            )
            unmatched_records.extend(
                {
                    "radc_model": radc_model,
                    "cell_class": cell_class,
                    "availability": "radc_only",
                    "anchor_unit_id": anchor,
                }
                for anchor in radc_only
            )
            unmatched_records.extend(
                {
                    "radc_model": radc_model,
                    "cell_class": cell_class,
                    "availability": "seaad_only",
                    "anchor_unit_id": anchor,
                }
                for anchor in seaad_only
            )
            for column in metadata_columns[1:]:
                left = f"{column}_radc_meta"
                right = f"{column}_seaad_meta"
                if not pair[left].astype(str).equals(pair[right].astype(str)):
                    raise RuntimeError(
                        f"Anchor metadata mismatch for {column}: "
                        f"{radc_model}/{cell_class}"
                    )
                pair[column] = pair[left]
                pair = pair.drop(columns=[left, right])
            pair["radc_model"] = radc_model
            pair["delta_seaad_minus_radc"] = (
                pair["beta_D_seaad"] - pair["beta_D_radc"]
            )
            pair["SE_delta"] = np.sqrt(
                np.square(pair["SE_D_seaad"])
                + np.square(pair["SE_D_radc"])
            )
            pair["heterogeneity_z"] = (
                pair["delta_seaad_minus_radc"] / pair["SE_delta"]
            )
            pair["heterogeneity_p"] = 2 * stats.norm.sf(
                np.abs(pair["heterogeneity_z"])
            )
            pair["heterogeneity_FDR_within_class"] = bh_adjust(
                pair["heterogeneity_p"]
            )
            pair["absolute_standardized_difference"] = pair[
                "heterogeneity_z"
            ].abs()
            pair["beta_sign_agreement"] = (
                np.sign(pair["beta_D_radc"])
                == np.sign(pair["beta_D_seaad"])
            )
            pair["genetic_coupling_sign_radc"] = np.sign(
                pair["beta_G"] * pair["beta_D_radc"]
            )
            pair["genetic_coupling_sign_seaad"] = np.sign(
                pair["beta_G"] * pair["beta_D_seaad"]
            )
            pairwise_parts.append(pair)
            summary = summarize_pairwise(pair)
            summary_records.append(
                {
                    "radc_model": radc_model,
                    "cell_class": cell_class,
                    **summary,
                }
            )
    pairwise = pd.concat(pairwise_parts, ignore_index=True)
    pairwise.to_csv(output_dir / "shared_anchor_pairwise.csv", index=False)
    summary = pd.DataFrame.from_records(summary_records)
    summary.to_csv(output_dir / "class_summary.csv", index=False)
    coverage = pd.DataFrame.from_records(coverage_records)
    coverage.to_csv(output_dir / "anchor_coverage_audit.csv", index=False)
    pd.DataFrame.from_records(
        unmatched_records,
        columns=[
            "radc_model",
            "cell_class",
            "availability",
            "anchor_unit_id",
        ],
    ).to_csv(output_dir / "unmatched_anchors.csv", index=False)

    stability_records: list[dict[str, object]] = []
    standard = radc_tables["standard_voom"]
    quality = radc_tables["quality_weights"]
    for cell_class in CLASSES:
        left = standard[standard["cell_class"].eq(cell_class)][
            ["anchor_unit_id", "beta_D_radc", "SE_D_radc"]
        ]
        right = quality[quality["cell_class"].eq(cell_class)][
            ["anchor_unit_id", "beta_D_radc", "SE_D_radc"]
        ].rename(
            columns={
                "beta_D_radc": "beta_D_quality",
                "SE_D_radc": "SE_D_quality",
            }
        )
        joined = left.merge(
            right, on="anchor_unit_id", validate="one_to_one"
        )
        correlation, correlation_p = safe_correlation(
            joined["beta_D_radc"], joined["beta_D_quality"], "pearson"
        )
        stability_records.append(
            {
                "cell_class": cell_class,
                "anchors": len(joined),
                "standard_quality_pearson_r": correlation,
                "standard_quality_pearson_p": correlation_p,
                "standard_quality_sign_agreement_fraction": float(
                    (
                        np.sign(joined["beta_D_radc"])
                        == np.sign(joined["beta_D_quality"])
                    ).mean()
                ),
                "maximum_absolute_beta_difference": float(
                    (
                        joined["beta_D_radc"] - joined["beta_D_quality"]
                    ).abs().max()
                ),
            }
        )
    stability = pd.DataFrame.from_records(stability_records)
    coupling_by_model: dict[str, pd.DataFrame] = {}
    for model, path in radc_coupling_paths.items():
        table = pd.read_csv(path)
        table = table[
            table["level"].eq("stratum")
            & table["cell_class"].isin(CLASSES)
        ][
            [
                "cell_class",
                "S_locus_equal_moderated_all_anchors",
                "empirical_p_two_sided",
            ]
        ].copy()
        if len(table) != len(CLASSES):
            raise RuntimeError(f"Missing RADC coupling strata in {path}")
        coupling_by_model[model] = table.rename(
            columns={
                "S_locus_equal_moderated_all_anchors": f"S_{model}",
                "empirical_p_two_sided": f"empirical_p_{model}",
            }
        )
    coupling_stability = coupling_by_model["standard_voom"].merge(
        coupling_by_model["quality_weights"],
        on="cell_class",
        validate="one_to_one",
    )
    coupling_stability["absolute_S_difference"] = (
        coupling_stability["S_standard_voom"]
        - coupling_stability["S_quality_weights"]
    ).abs()
    coupling_stability["same_S_sign"] = (
        np.sign(coupling_stability["S_standard_voom"])
        == np.sign(coupling_stability["S_quality_weights"])
    )
    stability = stability.merge(
        coupling_stability, on="cell_class", validate="one_to_one"
    )
    stability.to_csv(output_dir / "radc_model_stability.csv", index=False)
    build_scatter(pairwise, output_dir)

    standard_summary = summary[
        summary["radc_model"].eq("standard_voom")
    ].copy()
    technical_stable = bool(
        (stability["standard_quality_pearson_r"] >= 0.90).all()
        and stability["same_S_sign"].all()
        and (stability["absolute_S_difference"] <= 0.05).all()
    )
    portable = bool(
        (standard_summary["pearson_r"] > 0).all()
        and (
            standard_summary["beta_sign_agreement_fraction"] >= 0.60
        ).all()
    )
    decision = (
        "COHORT_HETEROGENEITY_NOT_MODEL_CHOICE"
        if technical_stable and not portable
        else "MIXED_TECHNICAL_AND_COHORT_HETEROGENEITY"
    )
    decision_json = {
        "status": "COMPLETE",
        "decision": decision,
        "technical_model_stability_gate": technical_stable,
        "cross_cohort_portability_gate": portable,
        "interpretation": (
            "RADC standard and quality-weight anchor effects are internally "
            "stable, while shared-anchor effects have limited portability to "
            "SEA-AD. The failed external replication is therefore not "
            "explained by choosing standard voom versus quality weights."
            if decision == "COHORT_HETEROGENEITY_NOT_MODEL_CHOICE"
            else "Both technical-model sensitivity and cohort heterogeneity "
            "may contribute; no single explanation is established."
        ),
        "boundaries": [
            "This is a post hoc failure decomposition, not a new confirmatory test.",
            "Cohort heterogeneity may include brain-region composition, donor "
            "selection, disease-stage distribution, and technical platform.",
            "The analysis cannot identify which heterogeneity source is causal.",
            "The technical stability gate requires shared-anchor beta Pearson "
            "r >= 0.90, the same coupling-S sign, and absolute coupling-S "
            "difference <= 0.05 for all three tested classes.",
        ],
        "generated_at": datetime.now(timezone.utc).astimezone().isoformat(),
    }
    (output_dir / "decision.json").write_text(
        json.dumps(decision_json, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )

    display_summary = standard_summary[
        [
            "cell_class",
            "anchors",
            "pearson_r",
            "spearman_rho",
            "beta_sign_agreement_fraction",
            "weighted_mean_delta_seaad_minus_radc",
            "weighted_mean_delta_p",
            "anchor_heterogeneity_FDR05",
        ]
    ].copy()
    for column in (
        "pearson_r",
        "spearman_rho",
        "beta_sign_agreement_fraction",
        "weighted_mean_delta_seaad_minus_radc",
        "weighted_mean_delta_p",
    ):
        display_summary[column] = display_summary[column].map(
            lambda value: f"{value:.4f}"
        )
    display_stability = stability.copy()
    for column in (
        "standard_quality_pearson_r",
        "standard_quality_pearson_p",
        "standard_quality_sign_agreement_fraction",
        "maximum_absolute_beta_difference",
        "S_standard_voom",
        "S_quality_weights",
        "absolute_S_difference",
    ):
        display_stability[column] = display_stability[column].map(
            lambda value: f"{value:.4f}"
        )
    lines = [
        "# RADC-to-SEA-AD Braak failure decomposition",
        "",
        f"- Decision: `{decision}`",
        f"- Interpretation: {decision_json['interpretation']}",
        "",
        "## Shared-anchor cross-cohort comparison",
        "",
        markdown_table(display_summary),
        "",
        "## RADC internal model stability",
        "",
        markdown_table(display_stability),
        "",
        "This analysis was performed after external non-replication and is",
        "descriptive failure decomposition, not a new confirmatory result.",
        "",
    ]
    (output_dir / "FAILURE_DECOMPOSITION.md").write_text(
        "\n".join(lines), encoding="utf-8"
    )
    print(json.dumps(decision_json, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
