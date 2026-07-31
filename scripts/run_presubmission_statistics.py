#!/usr/bin/env python3
"""Frozen pre-submission statistics for the brain-expression manuscript.

This script performs only the analyses authorized before result inspection:

1. Holm correction of the two RADC standard-voom global pathology tests;
2. independent-LD-locus bootstrap confidence intervals for existing coupling
   estimates; and
3. paired common-locus RADC minus SEA-AD contrasts for the six fixed
   cell-class-by-pathology combinations.

It does not fit expression models, add cell classes, change genetic anchors,
or select results according to significance.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
from scipy.stats import norm


BOOTSTRAPS = 10_000
SEED = 20_260_731
CELL_CLASSES = ("Immune", "Oligo", "EN")
PATHOLOGIES = ("CERAD", "Braak")
RADC_MODELS = ("standard_voom", "quality_weights")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    return parser.parse_args()


def timestamp() -> str:
    return datetime.now(timezone.utc).astimezone().isoformat(timespec="seconds")


def sha256(path: Path, block_size: int = 8 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        while block := stream.read(block_size):
            digest.update(block)
    return digest.hexdigest()


def atomic_json(path: Path, payload: object) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(payload, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    temporary.replace(path)


def adjust_bh(values: Iterable[float]) -> np.ndarray:
    series = pd.Series(list(values), dtype=float)
    result = np.full(len(series), np.nan)
    valid = np.flatnonzero(series.notna().to_numpy())
    if len(valid) == 0:
        return result
    p = series.iloc[valid].to_numpy()
    order = np.argsort(p)
    ranked = p[order]
    adjusted = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    restored = np.empty_like(adjusted)
    restored[order] = np.minimum(adjusted, 1.0)
    result[valid] = restored
    return result


def adjust_holm(values: Iterable[float]) -> np.ndarray:
    series = pd.Series(list(values), dtype=float)
    result = np.full(len(series), np.nan)
    valid = np.flatnonzero(series.notna().to_numpy())
    if len(valid) == 0:
        return result
    p = series.iloc[valid].to_numpy()
    order = np.argsort(p)
    ranked = p[order]
    raw = (len(ranked) - np.arange(len(ranked))) * ranked
    adjusted = np.maximum.accumulate(raw)
    restored = np.empty_like(adjusted)
    restored[order] = np.minimum(adjusted, 1.0)
    result[valid] = restored
    return result


def score_from_effects(frame: pd.DataFrame, prefix: str) -> pd.Series:
    beta_g = pd.to_numeric(frame[f"beta_G_{prefix}"], errors="coerce")
    se_g = pd.to_numeric(frame[f"SE_G_{prefix}"], errors="coerce")
    beta_d = pd.to_numeric(frame[f"beta_D_{prefix}"], errors="coerce")
    se_d = pd.to_numeric(frame[f"SE_D_{prefix}"], errors="coerce")
    valid = (
        beta_g.notna()
        & se_g.notna()
        & beta_d.notna()
        & se_d.notna()
        & se_g.gt(0)
        & se_d.gt(0)
    )
    result = pd.Series(np.nan, index=frame.index, dtype=float)
    p_g = norm.cdf((beta_g.loc[valid] / se_g.loc[valid]).to_numpy())
    p_d = norm.cdf((beta_d.loc[valid] / se_d.loc[valid]).to_numpy())
    p_aligned = p_g * p_d + (1.0 - p_g) * (1.0 - p_d)
    result.loc[valid] = 2.0 * p_aligned - 1.0
    return result


def locus_means(frame: pd.DataFrame, score_column: str) -> pd.Series:
    clean = frame.dropna(subset=["locus_id", score_column]).copy()
    clean[score_column] = pd.to_numeric(clean[score_column], errors="coerce")
    clean = clean.dropna(subset=[score_column])
    return clean.groupby("locus_id", observed=True)[score_column].mean()


def bootstrap_mean(
    values: np.ndarray,
    rng: np.random.Generator,
    replicates: int = BOOTSTRAPS,
) -> tuple[float, float, float, np.ndarray]:
    values = np.asarray(values, dtype=float)
    if values.ndim != 1 or len(values) == 0 or not np.isfinite(values).all():
        raise ValueError("Bootstrap values must be a finite non-empty vector")
    draws = rng.choice(values, size=(replicates, len(values)), replace=True)
    statistics = draws.mean(axis=1)
    lower, median, upper = np.quantile(statistics, [0.025, 0.5, 0.975])
    return float(lower), float(median), float(upper), statistics


def summarize_existing_estimate(
    frame: pd.DataFrame,
    score_column: str,
    source: str,
    model: str,
    pathology: str,
    level: str,
    disease: str,
    cell_class: str,
    rng: np.random.Generator,
) -> dict[str, object]:
    loci = locus_means(frame, score_column)
    lower, median, upper, _ = bootstrap_mean(loci.to_numpy(), rng)
    return {
        "source": source,
        "model": model,
        "pathology": pathology,
        "level": level,
        "disease": disease,
        "cell_class": cell_class,
        "anchor_rows": int(frame[score_column].notna().sum()),
        "loci": int(len(loci)),
        "S_observed": float(loci.mean()),
        "bootstrap_median": median,
        "bootstrap_CI95_lower": lower,
        "bootstrap_CI95_upper": upper,
        "low_locus_count_exploratory": bool(len(loci) < 10),
        "bootstrap_replicates": BOOTSTRAPS,
        "bootstrap_seed": SEED,
    }


def read_anchor_probability(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path, low_memory=False)
    required = {"disease", "cell_class", "locus_id", "S_moderated_rho_0"}
    missing = sorted(required - set(frame))
    if missing:
        raise RuntimeError(f"Missing columns in {path}: {missing}")
    return frame


def existing_bootstrap_registry(
    project: Path,
    rng: np.random.Generator,
    input_paths: set[Path],
) -> pd.DataFrame:
    records: list[dict[str, object]] = []

    mssm_path = (
        project
        / "outputs/phase4/mssm_direction_coupling"
        / "anchor_direction_probabilities.csv.gz"
    )
    input_paths.add(mssm_path)
    mssm = read_anchor_probability(mssm_path)
    records.append(
        summarize_existing_estimate(
            mssm,
            "S_moderated_rho_0",
            "MSSM",
            "voomWithQualityWeights",
            "case_control",
            "global",
            "ALL",
            "ALL",
            rng,
        )
    )
    for disease, disease_frame in mssm.groupby("disease", observed=True):
        records.append(
            summarize_existing_estimate(
                disease_frame,
                "S_moderated_rho_0",
                "MSSM",
                "voomWithQualityWeights",
                "case_control",
                "disease",
                str(disease),
                "ALL",
                rng,
            )
        )
    for (disease, cell_class), stratum in mssm.groupby(
        ["disease", "cell_class"], observed=True
    ):
        records.append(
            summarize_existing_estimate(
                stratum,
                "S_moderated_rho_0",
                "MSSM",
                "voomWithQualityWeights",
                "case_control",
                "stratum",
                str(disease),
                str(cell_class),
                rng,
            )
        )

    for model in RADC_MODELS:
        model_folder = (
            "radc_stage_coupling_standard_voom"
            if model == "standard_voom"
            else "radc_stage_coupling_quality_weights"
        )
        for pathology in PATHOLOGIES:
            path = (
                project
                / f"outputs/phase4/{model_folder}/{pathology}/coupling"
                / "anchor_direction_probabilities.csv.gz"
            )
            input_paths.add(path)
            frame = read_anchor_probability(path)
            records.append(
                summarize_existing_estimate(
                    frame,
                    "S_moderated_rho_0",
                    "RADC",
                    model,
                    pathology,
                    "global",
                    "AD",
                    "ALL",
                    rng,
                )
            )
            for cell_class, stratum in frame.groupby(
                "cell_class", observed=True
            ):
                records.append(
                    summarize_existing_estimate(
                        stratum,
                        "S_moderated_rho_0",
                        "RADC",
                        model,
                        pathology,
                        "stratum",
                        "AD",
                        str(cell_class),
                        rng,
                    )
                )

    for cell_class in CELL_CLASSES:
        for pathology in PATHOLOGIES:
            path = (
                project
                / f"outputs/phase5/seaad_{cell_class.lower()}_core3"
                / f"{pathology}/coupling/anchor_direction_probabilities.csv.gz"
            )
            input_paths.add(path)
            frame = read_anchor_probability(path)
            records.append(
                summarize_existing_estimate(
                    frame,
                    "S_moderated_rho_0",
                    "SEA-AD",
                    "three_region_summed_standard_voom",
                    pathology,
                    "stratum",
                    "AD",
                    cell_class,
                    rng,
                )
            )

    return pd.DataFrame.from_records(records)


def primary_multiplicity(
    project: Path, input_paths: set[Path]
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for pathology in PATHOLOGIES:
        path = (
            project
            / "outputs/phase4/radc_stage_coupling_standard_voom"
            / f"{pathology}/coupling/direction_coupling_group_results.csv"
        )
        input_paths.add(path)
        frame = pd.read_csv(path)
        global_row = frame.loc[frame["level"].eq("global")].iloc[0]
        rows.append(
            {
                "cohort": "RADC",
                "model": "standard_voom",
                "family": "two_co_primary_global_pathology_axes",
                "pathology": pathology,
                "S": float(
                    global_row["S_locus_equal_moderated_all_anchors"]
                ),
                "empirical_p_raw": float(global_row["empirical_p_two_sided"]),
            }
        )
    result = pd.DataFrame(rows)
    result["empirical_p_Holm"] = adjust_holm(result["empirical_p_raw"])
    result["passes_Holm_FWER05"] = result["empirical_p_Holm"].lt(0.05)
    return result


def anchor_input_path(
    project: Path, cohort: str, cell_class: str, pathology: str
) -> Path:
    if cohort == "RADC":
        return (
            project
            / "outputs/phase4/radc_stage_coupling_standard_voom"
            / f"{pathology}/input"
            / f"radc_{pathology}_G1_G2_anchor_pathology_effects.csv"
        )
    return (
        project
        / f"outputs/phase5/seaad_{cell_class.lower()}_core3"
        / f"{pathology}/input"
        / f"seaad_{pathology}_{cell_class}_G1_G2_anchor_pathology_effects.csv"
    )


def read_anchor_input(path: Path, prefix: str) -> pd.DataFrame:
    frame = pd.read_csv(path, low_memory=False)
    required = {
        "anchor_unit_id",
        "locus_id",
        "cell_class",
        "beta_G",
        "SE_G",
        "beta_D",
        "SE_D",
        "expression_estimable",
    }
    missing = sorted(required - set(frame))
    if missing:
        raise RuntimeError(f"Missing columns in {path}: {missing}")
    estimable = (
        frame["expression_estimable"]
        .astype(str)
        .str.strip()
        .str.lower()
        .isin({"true", "1", "yes"})
    )
    frame = frame.loc[estimable].copy()
    if frame["anchor_unit_id"].duplicated().any():
        duplicated = frame.loc[
            frame["anchor_unit_id"].duplicated(False), "anchor_unit_id"
        ].tolist()
        raise RuntimeError(f"Duplicated anchors in {path}: {duplicated[:10]}")
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
            "locus_id": f"locus_id_{prefix}",
            "cell_class": f"cell_class_{prefix}",
            "beta_G": f"beta_G_{prefix}",
            "SE_G": f"SE_G_{prefix}",
            "beta_D": f"beta_D_{prefix}",
            "SE_D": f"SE_D_{prefix}",
        }
    )


def direct_common_locus_contrasts(
    project: Path,
    rng: np.random.Generator,
    input_paths: set[Path],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    summary_rows: list[dict[str, object]] = []
    locus_rows: list[pd.DataFrame] = []

    for cell_class in CELL_CLASSES:
        for pathology in PATHOLOGIES:
            radc_path = anchor_input_path(
                project, "RADC", cell_class, pathology
            )
            seaad_path = anchor_input_path(
                project, "SEA-AD", cell_class, pathology
            )
            input_paths.update({radc_path, seaad_path})
            radc = read_anchor_input(radc_path, "radc")
            radc = radc.loc[radc["cell_class_radc"].eq(cell_class)].copy()
            seaad = read_anchor_input(seaad_path, "seaad")
            seaad = seaad.loc[seaad["cell_class_seaad"].eq(cell_class)].copy()

            paired = radc.merge(
                seaad,
                on="anchor_unit_id",
                how="inner",
                validate="one_to_one",
            )
            if paired.empty:
                raise RuntimeError(
                    f"No common anchors for {cell_class}/{pathology}"
                )
            if not paired["locus_id_radc"].eq(
                paired["locus_id_seaad"]
            ).all():
                raise RuntimeError(
                    f"Locus identifiers disagree for {cell_class}/{pathology}"
                )
            paired["locus_id"] = paired["locus_id_radc"]
            paired["S_radc"] = score_from_effects(paired, "radc")
            paired["S_seaad"] = score_from_effects(paired, "seaad")
            paired = paired.dropna(subset=["S_radc", "S_seaad"]).copy()

            by_locus = (
                paired.groupby("locus_id", observed=True)
                .agg(
                    anchors=("anchor_unit_id", "nunique"),
                    S_radc=("S_radc", "mean"),
                    S_seaad=("S_seaad", "mean"),
                )
                .reset_index()
            )
            by_locus["Delta_S_radc_minus_seaad"] = (
                by_locus["S_radc"] - by_locus["S_seaad"]
            )
            deltas = by_locus["Delta_S_radc_minus_seaad"].to_numpy()
            lower, median, upper, bootstrap = bootstrap_mean(deltas, rng)
            observed = float(deltas.mean())
            centered = bootstrap - float(bootstrap.mean())
            p_value = (
                1.0
                + float(np.sum(np.abs(centered) >= abs(observed)))
            ) / (BOOTSTRAPS + 1.0)

            summary_rows.append(
                {
                    "cell_class": cell_class,
                    "pathology": pathology,
                    "common_anchor_rows": int(len(paired)),
                    "common_loci": int(len(by_locus)),
                    "S_RADC_common": float(by_locus["S_radc"].mean()),
                    "S_SEAAD_common": float(by_locus["S_seaad"].mean()),
                    "Delta_S_RADC_minus_SEAAD": observed,
                    "bootstrap_median_Delta": median,
                    "bootstrap_CI95_lower": lower,
                    "bootstrap_CI95_upper": upper,
                    "bootstrap_centered_two_sided_p": p_value,
                    "low_locus_count_exploratory": bool(len(by_locus) < 10),
                    "bootstrap_replicates": BOOTSTRAPS,
                    "bootstrap_seed": SEED,
                }
            )
            by_locus.insert(0, "pathology", pathology)
            by_locus.insert(0, "cell_class", cell_class)
            locus_rows.append(by_locus)

    summary = pd.DataFrame.from_records(summary_rows)
    summary["bootstrap_p_BH_six_contrasts"] = adjust_bh(
        summary["bootstrap_centered_two_sided_p"]
    )
    summary["passes_BH_FDR05"] = summary[
        "bootstrap_p_BH_six_contrasts"
    ].lt(0.05)
    locus_table = pd.concat(locus_rows, ignore_index=True)
    return summary, locus_table


def markdown_table(frame: pd.DataFrame) -> str:
    display = frame.copy()
    for column in display.select_dtypes(include=[np.number]).columns:
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
        "| "
        + " | ".join(str(value).replace("|", "\\|") for value in row)
        + " |"
        for row in display.itertuples(index=False, name=None)
    ]
    return "\n".join([header, separator, *rows])


def main() -> None:
    args = parse_args()
    project = args.project.resolve()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)
    input_paths: set[Path] = set()

    multiplicity = primary_multiplicity(project, input_paths)
    bootstrap = existing_bootstrap_registry(project, rng, input_paths)
    direct, direct_loci = direct_common_locus_contrasts(
        project, rng, input_paths
    )

    multiplicity_path = output_dir / "radc_primary_multiplicity.csv"
    bootstrap_path = output_dir / "locus_bootstrap_intervals.csv"
    direct_path = output_dir / "radc_seaad_direct_common_locus_contrasts.csv"
    direct_loci_path = output_dir / "radc_seaad_common_locus_components.csv"
    multiplicity.to_csv(multiplicity_path, index=False)
    bootstrap.to_csv(bootstrap_path, index=False)
    direct.to_csv(direct_path, index=False)
    direct_loci.to_csv(direct_loci_path, index=False)

    report_lines = [
        "# Frozen pre-submission statistical strengthening",
        "",
        f"Run time: {timestamp()}  ",
        f"Bootstrap replicates: {BOOTSTRAPS:,}  ",
        f"Seed: {SEED}",
        "",
        "## RADC co-primary pathology multiplicity",
        "",
        markdown_table(multiplicity),
        "",
        "Holm correction is applied only to the two standard-voom global "
        "co-primary pathology tests. Quality weights remain a sensitivity "
        "analysis rather than a second discovery family.",
        "",
        "## Direct common-locus RADC–SEA-AD contrasts",
        "",
        markdown_table(direct),
        "",
        "The paired contrast resamples independent LD loci. Its two-sided "
        "bootstrap P value compares the observed mean difference with the "
        "centered bootstrap null distribution. BH correction is applied "
        "across the six fixed contrasts. Comparisons with fewer than 10 loci "
        "remain descriptive even after correction.",
        "",
        "## Claim boundary",
        "",
        "- A confidence interval containing zero is not evidence of "
        "equivalence.",
        "- A significant result in one cohort and a nonsignificant result in "
        "another is not evidence of a cohort difference.",
        "- No cell-specific, reverse, protective, compensatory, or causal "
        "claim is authorized by this analysis.",
        "",
    ]
    report_path = output_dir / "PRESUBMISSION_STATISTICAL_REPORT.md"
    report_path.write_text("\n".join(report_lines), encoding="utf-8")

    output_paths = [
        multiplicity_path,
        bootstrap_path,
        direct_path,
        direct_loci_path,
        report_path,
    ]
    manifest = {
        "created_at": timestamp(),
        "script": str(Path(__file__).resolve()),
        "script_sha256": sha256(Path(__file__).resolve()),
        "project": str(project),
        "bootstrap_replicates": BOOTSTRAPS,
        "seed": SEED,
        "input_files": [
            {
                "path": str(path),
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
            for path in sorted(input_paths)
        ],
        "output_files": [
            {
                "path": str(path),
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
            for path in output_paths
        ],
        "fixed_families": {
            "radc_primary": ["CERAD", "Braak"],
            "direct_contrasts": [
                f"{cell_class}__{pathology}"
                for cell_class in CELL_CLASSES
                for pathology in PATHOLOGIES
            ],
        },
    }
    atomic_json(output_dir / "run_manifest.json", manifest)

    print(
        "COMPLETE "
        f"bootstrap_rows={len(bootstrap)} "
        f"direct_contrasts={len(direct)} "
        f"output={output_dir}"
    )


if __name__ == "__main__":
    main()
