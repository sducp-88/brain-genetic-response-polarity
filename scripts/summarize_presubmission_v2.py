#!/usr/bin/env python3
"""Build the frozen presubmission v2 statistical synthesis tables."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import norm


BOOTSTRAPS = 10_000
SEED = 20_260_731
AXES = ("CERAD", "Braak")
CLASSES = ("Immune", "Oligo", "EN")
CONTINUOUS_MODELS = {
    "standard": "outputs/phase4/radc_stage_coupling_standard_voom",
    "quality_weights": "outputs/phase4/radc_stage_coupling_quality_weights",
    "min_cells_50": (
        "outputs/phase6/radc_presubmission_sensitivities/min_cells_50"
    ),
    "exclude_age89plus": (
        "outputs/phase6/radc_presubmission_sensitivities/exclude_age89plus"
    ),
    "omit_nonad": (
        "outputs/phase6/radc_presubmission_sensitivities/omit_nonad"
    ),
    "omit_log_ncells": (
        "outputs/phase6/radc_presubmission_sensitivities/omit_log_ncells"
    ),
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


def adjust_holm(values: pd.Series) -> np.ndarray:
    p = values.to_numpy(float)
    order = np.argsort(p)
    ranked = p[order]
    adjusted = np.maximum.accumulate(
        (len(ranked) - np.arange(len(ranked))) * ranked
    )
    restored = np.empty_like(adjusted)
    restored[order] = np.minimum(adjusted, 1.0)
    return restored


def adjust_bh(values: pd.Series) -> np.ndarray:
    p = values.to_numpy(float)
    order = np.argsort(p)
    ranked = p[order]
    adjusted = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    restored = np.empty_like(adjusted)
    restored[order] = np.minimum(adjusted, 1.0)
    return restored


def locus_values(frame: pd.DataFrame, score_column: str) -> pd.Series:
    table = frame[["locus_id", score_column]].copy()
    table[score_column] = pd.to_numeric(table[score_column], errors="coerce")
    table = table.dropna()
    return table.groupby("locus_id", observed=True)[score_column].mean()


def bootstrap_mean(
    values: np.ndarray, rng: np.random.Generator
) -> tuple[float, float, float, np.ndarray]:
    values = np.asarray(values, dtype=float)
    if not len(values) or not np.isfinite(values).all():
        raise RuntimeError("Invalid locus vector")
    draws = rng.choice(values, size=(BOOTSTRAPS, len(values)), replace=True)
    statistics = draws.mean(axis=1)
    lower, median, upper = np.quantile(statistics, [0.025, 0.5, 0.975])
    return float(lower), float(median), float(upper), statistics


def one_global(group_path: Path) -> pd.Series:
    table = pd.read_csv(group_path)
    selected = table.loc[
        table["level"].eq("global")
        & table["disease"].eq("ALL")
        & table["cell_class"].eq("ALL")
    ]
    if len(selected) != 1:
        raise RuntimeError(f"Expected one global row in {group_path}")
    return selected.iloc[0]


def summarize_radc(
    project: Path, rng: np.random.Generator, inputs: set[Path]
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    continuous: list[dict[str, object]] = []
    for model, relative in CONTINUOUS_MODELS.items():
        for axis in AXES:
            root = project / relative / axis / "coupling"
            group_path = root / "direction_coupling_group_results.csv"
            anchor_path = root / "anchor_direction_probabilities.csv.gz"
            inputs.update({group_path, anchor_path})
            group = one_global(group_path)
            anchors = pd.read_csv(anchor_path, low_memory=False)
            loci = locus_values(anchors, "S_moderated_rho_0")
            lower, median, upper, _ = bootstrap_mean(loci.to_numpy(), rng)
            continuous.append(
                {
                    "model": model,
                    "pathology": axis,
                    "S": float(
                        group["S_locus_equal_moderated_all_anchors"]
                    ),
                    "matched_null_mean": float(group["matched_null_mean"]),
                    "empirical_p_raw": float(group["empirical_p_two_sided"]),
                    "bootstrap_median": median,
                    "bootstrap_CI95_lower": lower,
                    "bootstrap_CI95_upper": upper,
                    "leave_one_locus_min": float(group["leave_one_locus_min"]),
                    "leave_one_locus_max": float(group["leave_one_locus_max"]),
                    "maximum_absolute_locus_contribution_share": float(
                        group["maximum_absolute_locus_contribution_share"]
                    ),
                    "anchor_rows": int(group["anchor_rows"]),
                    "loci": int(group["loci"]),
                }
            )
    continuous_table = pd.DataFrame(continuous)
    continuous_table["empirical_p_Holm_primary_family"] = np.nan
    primary = continuous_table["model"].eq("standard")
    continuous_table.loc[
        primary, "empirical_p_Holm_primary_family"
    ] = adjust_holm(
        continuous_table.loc[primary, "empirical_p_raw"]
    )
    continuous_table["passes_primary_Holm05"] = (
        continuous_table["empirical_p_Holm_primary_family"].lt(0.05)
    )

    stage: list[dict[str, object]] = []
    stage_root = (
        project
        / "outputs/phase6/radc_presubmission_sensitivities/stage_contrasts"
    )
    for axis in AXES:
        for contrast in ("middle_vs_low", "high_vs_low"):
            root = stage_root / axis / contrast / "coupling"
            group_path = root / "direction_coupling_group_results.csv"
            anchor_path = root / "anchor_direction_probabilities.csv.gz"
            inputs.update({group_path, anchor_path})
            group = one_global(group_path)
            anchors = pd.read_csv(anchor_path, low_memory=False)
            loci = locus_values(anchors, "S_moderated_rho_0")
            lower, median, upper, _ = bootstrap_mean(loci.to_numpy(), rng)
            stage.append(
                {
                    "pathology": axis,
                    "contrast": contrast,
                    "S": float(
                        group["S_locus_equal_moderated_all_anchors"]
                    ),
                    "empirical_p_raw": float(group["empirical_p_two_sided"]),
                    "bootstrap_median": median,
                    "bootstrap_CI95_lower": lower,
                    "bootstrap_CI95_upper": upper,
                    "leave_one_locus_min": float(group["leave_one_locus_min"]),
                    "leave_one_locus_max": float(group["leave_one_locus_max"]),
                    "maximum_absolute_locus_contribution_share": float(
                        group["maximum_absolute_locus_contribution_share"]
                    ),
                    "loci": int(group["loci"]),
                }
            )
    stage_table = pd.DataFrame(stage)

    hc3: list[dict[str, object]] = []
    hc3_root = project / "outputs/phase6/radc_presubmission_sensitivities/hc3"
    for axis in AXES:
        root = hc3_root / axis / "coupling"
        group_path = root / "direction_coupling_group_results.csv"
        anchor_path = root / "anchor_direction_probabilities.csv.gz"
        inputs.update({group_path, anchor_path})
        group = one_global(group_path)
        anchors = pd.read_csv(anchor_path, low_memory=False)
        loci = locus_values(anchors, "S_HC3_rho_0")
        lower, median, upper, _ = bootstrap_mean(loci.to_numpy(), rng)
        hc3.append(
            {
                "pathology": axis,
                "S_HC3": float(
                    group["S_locus_equal_HC3_matched_anchors"]
                ),
                "bootstrap_median": median,
                "bootstrap_CI95_lower": lower,
                "bootstrap_CI95_upper": upper,
                "anchor_rows": int(group["matched_anchor_rows"]),
                "loci": int(group["loci"]),
            }
        )
    return continuous_table, stage_table, pd.DataFrame(hc3)


def score_from_effects(frame: pd.DataFrame, prefix: str) -> pd.Series:
    beta_g = pd.to_numeric(frame[f"beta_G_{prefix}"], errors="coerce")
    se_g = pd.to_numeric(frame[f"SE_G_{prefix}"], errors="coerce")
    beta_d = pd.to_numeric(frame[f"beta_D_{prefix}"], errors="coerce")
    se_d = pd.to_numeric(frame[f"SE_D_{prefix}"], errors="coerce")
    valid = (
        beta_g.notna()
        & se_g.gt(0)
        & beta_d.notna()
        & se_d.gt(0)
    )
    result = pd.Series(np.nan, index=frame.index, dtype=float)
    pg = norm.cdf((beta_g.loc[valid] / se_g.loc[valid]).to_numpy())
    pdisease = norm.cdf((beta_d.loc[valid] / se_d.loc[valid]).to_numpy())
    result.loc[valid] = 2 * (
        pg * pdisease + (1 - pg) * (1 - pdisease)
    ) - 1
    return result


def read_anchor_input(path: Path, prefix: str) -> pd.DataFrame:
    frame = pd.read_csv(path, low_memory=False)
    estimable = (
        frame["expression_estimable"]
        .astype(str)
        .str.lower()
        .isin({"true", "1", "yes"})
    )
    frame = frame.loc[estimable].copy()
    columns = [
        "anchor_unit_id",
        "locus_id",
        "cell_class",
        "beta_G",
        "SE_G",
        "beta_D",
        "SE_D",
    ]
    return frame[columns].rename(
        columns={
            column: f"{column}_{prefix}"
            for column in columns
            if column != "anchor_unit_id"
        }
    )


def direct_radc_seaad_dfc(
    project: Path, rng: np.random.Generator, inputs: set[Path]
) -> tuple[pd.DataFrame, pd.DataFrame]:
    summaries: list[dict[str, object]] = []
    components: list[pd.DataFrame] = []
    for cell_class in CLASSES:
        for axis in AXES:
            radc_path = (
                project
                / "outputs/phase4/radc_stage_coupling_standard_voom"
                / axis
                / "input"
                / f"radc_{axis}_G1_G2_anchor_pathology_effects.csv"
            )
            seaad_path = (
                project
                / "outputs/phase6/seaad_region_sensitivity_v2"
                / cell_class
                / "coupling/DFC_standard"
                / axis
                / "input"
                / f"seaad_{axis}_{cell_class}_G1_G2_anchor_pathology_effects.csv"
            )
            inputs.update({radc_path, seaad_path})
            radc = read_anchor_input(radc_path, "radc")
            seaad = read_anchor_input(seaad_path, "seaad")
            radc = radc.loc[
                radc["cell_class_radc"].eq(cell_class)
            ].copy()
            seaad = seaad.loc[
                seaad["cell_class_seaad"].eq(cell_class)
            ].copy()
            paired = radc.merge(
                seaad, on="anchor_unit_id", how="inner", validate="one_to_one"
            )
            if not paired["locus_id_radc"].eq(
                paired["locus_id_seaad"]
            ).all():
                raise RuntimeError("Paired locus mismatch")
            paired["locus_id"] = paired["locus_id_radc"]
            paired["S_RADC"] = score_from_effects(paired, "radc")
            paired["S_SEAAD_DFC"] = score_from_effects(paired, "seaad")
            by_locus = (
                paired.groupby("locus_id", observed=True)
                .agg(
                    anchors=("anchor_unit_id", "nunique"),
                    S_RADC=("S_RADC", "mean"),
                    S_SEAAD_DFC=("S_SEAAD_DFC", "mean"),
                )
                .reset_index()
            )
            by_locus["Delta_S_RADC_minus_SEAAD_DFC"] = (
                by_locus["S_RADC"] - by_locus["S_SEAAD_DFC"]
            )
            values = by_locus["Delta_S_RADC_minus_SEAAD_DFC"].to_numpy()
            lower, median, upper, bootstrap = bootstrap_mean(values, rng)
            observed = float(values.mean())
            centered = bootstrap - bootstrap.mean()
            p_value = (1 + np.sum(np.abs(centered) >= abs(observed))) / (
                BOOTSTRAPS + 1
            )
            summaries.append(
                {
                    "cell_class": cell_class,
                    "pathology": axis,
                    "common_anchor_rows": int(len(paired)),
                    "common_loci": int(len(by_locus)),
                    "S_RADC_common": float(by_locus["S_RADC"].mean()),
                    "S_SEAAD_DFC_common": float(
                        by_locus["S_SEAAD_DFC"].mean()
                    ),
                    "Delta_S_RADC_minus_SEAAD_DFC": observed,
                    "bootstrap_median_Delta": median,
                    "bootstrap_CI95_lower": lower,
                    "bootstrap_CI95_upper": upper,
                    "bootstrap_centered_two_sided_p": float(p_value),
                    "low_locus_count_exploratory": bool(len(by_locus) < 10),
                }
            )
            by_locus.insert(0, "pathology", axis)
            by_locus.insert(0, "cell_class", cell_class)
            components.append(by_locus)
    summary = pd.DataFrame(summaries)
    summary["bootstrap_p_BH_six_contrasts"] = adjust_bh(
        summary["bootstrap_centered_two_sided_p"]
    )
    summary["passes_BH05"] = summary[
        "bootstrap_p_BH_six_contrasts"
    ].lt(0.05)
    return summary, pd.concat(components, ignore_index=True)


def markdown(frame: pd.DataFrame) -> str:
    display = frame.copy()
    for column in display.select_dtypes(include=np.number).columns:
        display[column] = display[column].map(
            lambda value: ""
            if pd.isna(value)
            else f"{value:.4f}"
            if isinstance(value, (float, np.floating))
            else str(value)
        )
    header = "| " + " | ".join(display.columns) + " |"
    separator = "| " + " | ".join("---" for _ in display.columns) + " |"
    rows = [
        "| " + " | ".join(map(str, row)) + " |"
        for row in display.itertuples(index=False, name=None)
    ]
    return "\n".join([header, separator, *rows])


def main() -> None:
    args = parse_args()
    project = args.project.resolve()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)
    inputs: set[Path] = set()

    continuous, stage, hc3 = summarize_radc(project, rng, inputs)
    direct, direct_loci = direct_radc_seaad_dfc(project, rng, inputs)
    seaad_path = (
        project
        / "outputs/phase6/seaad_region_sensitivity_v2/summary"
        / "seaad_fixed_42_results.csv"
    )
    seaad_regional_path = (
        project
        / "outputs/phase6/seaad_region_sensitivity_v2/summary"
        / "seaad_regional_paired_contrasts.csv"
    )
    inputs.update({seaad_path, seaad_regional_path})

    paths = {
        "radc_continuous_sensitivities.csv": continuous,
        "radc_stage_contrasts.csv": stage,
        "radc_hc3_sensitivity.csv": hc3,
        "radc_seaad_dfc_direct_contrasts.csv": direct,
        "radc_seaad_dfc_locus_components.csv": direct_loci,
    }
    for name, table in paths.items():
        table.to_csv(output / name, index=False)

    primary = continuous.loc[continuous["model"].eq("standard")].copy()
    report = [
        "# Presubmission statistical synthesis v2",
        "",
        f"Generated: {datetime.now().astimezone().isoformat(timespec='seconds')}  ",
        f"Locus bootstrap: {BOOTSTRAPS:,}; seed {SEED}",
        "",
        "## RADC primary family",
        "",
        markdown(
            primary[
                [
                    "pathology",
                    "S",
                    "bootstrap_CI95_lower",
                    "bootstrap_CI95_upper",
                    "empirical_p_raw",
                    "empirical_p_Holm_primary_family",
                    "passes_primary_Holm05",
                ]
            ]
        ),
        "",
        "No RADC co-primary test passes Holm family-wise correction.",
        "",
        "## RADC fixed continuous sensitivities",
        "",
        markdown(
            continuous[
                [
                    "model",
                    "pathology",
                    "S",
                    "bootstrap_CI95_lower",
                    "bootstrap_CI95_upper",
                    "empirical_p_raw",
                    "leave_one_locus_min",
                    "leave_one_locus_max",
                ]
            ]
        ),
        "",
        "## RADC fixed stage contrasts",
        "",
        markdown(stage),
        "",
        "## Direct RADC minus SEA-AD DFC contrasts",
        "",
        markdown(direct),
        "",
        "## Interpretation boundary",
        "",
        "- Sensitivities cannot replace the primary model.",
        "- Confidence intervals containing zero do not establish equivalence.",
        "- Fewer than 10 loci is exploratory even after numerical correction.",
        "- Opposed direction is not evidence of protection or compensation.",
        "",
    ]
    report_path = output / "PRESUBMISSION_SYNTHESIS_V2_REPORT.md"
    report_path.write_text("\n".join(report), encoding="utf-8")

    outputs = [output / name for name in paths] + [report_path]
    manifest = {
        "status": "COMPLETE",
        "generated_at": datetime.now().astimezone().isoformat(),
        "bootstrap_replicates": BOOTSTRAPS,
        "seed": SEED,
        "script": str(Path(__file__).resolve()),
        "script_sha256": sha256(Path(__file__).resolve()),
        "input_files": [
            {"path": str(path), "sha256": sha256(path)}
            for path in sorted(inputs)
        ],
        "output_files": [
            {"path": str(path), "sha256": sha256(path)} for path in outputs
        ],
        "development_invalid_exclusions": [
            "development_invalid_exclude_age90plus_kept_20260731_0920",
            "development_invalid_numeric_age89_exclusion_kept_20260731_0923",
        ],
    }
    manifest_path = output / "run_manifest.json"
    temporary = manifest_path.with_suffix(".json.tmp")
    temporary.write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    os.replace(temporary, manifest_path)
    print(
        f"COMPLETE continuous={len(continuous)} stage={len(stage)} "
        f"direct={len(direct)} output={output}"
    )


if __name__ == "__main__":
    main()

