#!/usr/bin/env python3
"""Analyze genetic-risk versus diseased-brain expression direction coupling.

The script implements the pre-result freeze dated 2026-07-30:
continuous direction probabilities, locus-equal aggregation, a matched
disease-expression null, empirical omnibus tests, leave-one-locus-out
diagnostics, HC3 sensitivity, and an error-correlation rho grid.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy.stats import multivariate_normal, norm


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--anchor-results", type=Path, required=True)
    parser.add_argument("--all-gene-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--permutations", type=int, default=10_000)
    parser.add_argument("--candidate-count", type=int, default=50)
    parser.add_argument("--minimum-candidates", type=int, default=20)
    parser.add_argument("--minimum-loci", type=int, default=5)
    parser.add_argument("--seed", type=int, default=20_260_730)
    return parser.parse_args()


def now() -> str:
    return datetime.now(timezone.utc).astimezone().isoformat(timespec="seconds")


def sha256(path: Path, block_size: int = 8 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(block_size):
            digest.update(block)
    return digest.hexdigest()


def atomic_json(path: Path, payload: dict[str, Any]) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    os.replace(temporary, path)


def as_bool(series: pd.Series) -> pd.Series:
    return series.astype(str).str.lower().isin({"true", "1", "yes", "y"})


def bh_adjust(p_values: pd.Series) -> pd.Series:
    result = pd.Series(np.nan, index=p_values.index, dtype=float)
    valid = p_values.notna()
    values = p_values.loc[valid].to_numpy(float)
    if not values.size:
        return result
    order = np.argsort(values)
    ranked = values[order]
    adjusted = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.clip(adjusted, 0, 1)
    inverse = np.empty_like(order)
    inverse[order] = np.arange(len(order))
    result.loc[valid] = adjusted[inverse]
    return result


def independent_score(
    beta_g: np.ndarray,
    se_g: np.ndarray,
    beta_d: np.ndarray,
    se_d: np.ndarray,
) -> np.ndarray:
    genetic_sign = 2 * norm.cdf(beta_g / se_g) - 1
    disease_sign = 2 * norm.cdf(beta_d / se_d) - 1
    return genetic_sign * disease_sign


def correlated_score(
    beta_g: float,
    se_g: float,
    beta_d: float,
    se_d: float,
    rho: float,
) -> float:
    threshold_g = -beta_g / se_g
    threshold_d = -beta_d / se_d
    joint_below = multivariate_normal.cdf(
        [threshold_g, threshold_d],
        mean=[0.0, 0.0],
        cov=[[1.0, rho], [rho, 1.0]],
    )
    aligned = (
        1
        - norm.cdf(threshold_g)
        - norm.cdf(threshold_d)
        + 2 * joint_below
    )
    return float(np.clip(2 * aligned - 1, -1, 1))


def locus_equal_stat(frame: pd.DataFrame, score: str) -> tuple[float, pd.Series]:
    locus_means = frame.groupby("locus_id", observed=True)[score].mean()
    return float(locus_means.mean()), locus_means


def group_definitions(frame: pd.DataFrame) -> list[dict[str, Any]]:
    groups: list[dict[str, Any]] = [
        {
            "level": "global",
            "disease": "ALL",
            "cell_class": "ALL",
            "indices": frame.index.to_numpy(),
        }
    ]
    for disease, disease_frame in frame.groupby("disease", observed=True):
        groups.append(
            {
                "level": "disease",
                "disease": str(disease),
                "cell_class": "ALL",
                "indices": disease_frame.index.to_numpy(),
            }
        )
    for (disease, cell_class), stratum in frame.groupby(
        ["disease", "cell_class"], observed=True
    ):
        groups.append(
            {
                "level": "stratum",
                "disease": str(disease),
                "cell_class": str(cell_class),
                "indices": stratum.index.to_numpy(),
            }
        )
    return groups


def read_all_gene_table(
    directory: Path, disease: str, cell_class: str
) -> pd.DataFrame:
    path = directory / f"{disease}__{cell_class}_all_gene_results.csv.gz"
    if not path.is_file():
        raise FileNotFoundError(path)
    required = [
        "gene_id_clean",
        "beta_D",
        "SE_D",
        "average_log2_expression",
    ]
    table = pd.read_csv(path, usecols=required)
    if table["gene_id_clean"].duplicated().any():
        raise RuntimeError(f"Duplicate gene identifiers in {path}.")
    for column in ["beta_D", "SE_D", "average_log2_expression"]:
        table[column] = pd.to_numeric(table[column], errors="coerce")
    table = table.loc[
        np.isfinite(table["beta_D"])
        & np.isfinite(table["SE_D"])
        & table["SE_D"].gt(0)
        & np.isfinite(table["average_log2_expression"])
    ].copy()
    return table


def robust_coordinates(
    background: pd.DataFrame, queries: pd.DataFrame
) -> tuple[np.ndarray, np.ndarray]:
    background_values = np.column_stack(
        [
            background["average_log2_expression"].to_numpy(float),
            np.log(background["SE_D"].to_numpy(float)),
        ]
    )
    query_values = np.column_stack(
        [
            queries["average_log2_expression"].to_numpy(float),
            np.log(queries["SE_D"].to_numpy(float)),
        ]
    )
    center = np.nanmedian(background_values, axis=0)
    scale = np.nanmedian(np.abs(background_values - center), axis=0) * 1.4826
    fallback = np.nanstd(background_values, axis=0, ddof=1)
    scale = np.where(np.isfinite(scale) & (scale > 0), scale, fallback)
    scale = np.where(np.isfinite(scale) & (scale > 0), scale, 1.0)
    return (
        (background_values - center) / scale,
        (query_values - center) / scale,
    )


def build_matches(
    anchors: pd.DataFrame,
    all_gene_dir: Path,
    candidate_count: int,
    minimum_candidates: int,
) -> tuple[dict[int, dict[str, np.ndarray]], pd.DataFrame]:
    matches: dict[int, dict[str, np.ndarray]] = {}
    diagnostics: list[dict[str, Any]] = []
    for (disease, cell_class), stratum in anchors.groupby(
        ["disease", "cell_class"], observed=True
    ):
        all_genes = read_all_gene_table(all_gene_dir, disease, cell_class)
        anchor_gene_ids = set(stratum["gene_id_clean"].astype(str))
        background = all_genes.loc[
            ~all_genes["gene_id_clean"].astype(str).isin(anchor_gene_ids)
        ].copy()
        anchor_precision = all_genes.set_index("gene_id_clean").reindex(
            stratum["gene_id_clean"].astype(str)
        )
        if anchor_precision[
            ["SE_D", "average_log2_expression"]
        ].isna().any().any():
            raise RuntimeError(
                f"Anchor genes missing from all-gene results: "
                f"{disease} {cell_class}."
            )
        local_queries = stratum.copy()
        local_queries["SE_D"] = anchor_precision["SE_D"].to_numpy()
        local_queries["average_log2_expression"] = anchor_precision[
            "average_log2_expression"
        ].to_numpy()
        if len(background) < minimum_candidates:
            for row_index in stratum.index:
                diagnostics.append(
                    {
                        "row_index": int(row_index),
                        "disease": disease,
                        "cell_class": cell_class,
                        "candidate_count": int(len(background)),
                        "eligible_for_matched_null": False,
                        "nearest_distance": np.nan,
                        "farthest_retained_distance": np.nan,
                    }
                )
            continue
        background_coords, query_coords = robust_coordinates(
            background, local_queries
        )
        retained_count = min(candidate_count, len(background))
        tree = cKDTree(background_coords)
        distances, positions = tree.query(query_coords, k=retained_count)
        if retained_count == 1:
            distances = distances[:, None]
            positions = positions[:, None]
        for local_position, row_index in enumerate(stratum.index):
            candidate_rows = background.iloc[positions[local_position]]
            matches[int(row_index)] = {
                "beta_D": candidate_rows["beta_D"].to_numpy(float),
                "SE_D": candidate_rows["SE_D"].to_numpy(float),
            }
            diagnostics.append(
                {
                    "row_index": int(row_index),
                    "disease": disease,
                    "cell_class": cell_class,
                    "candidate_count": retained_count,
                    "eligible_for_matched_null": (
                        retained_count >= minimum_candidates
                    ),
                    "nearest_distance": float(np.min(distances[local_position])),
                    "farthest_retained_distance": float(
                        np.max(distances[local_position])
                    ),
                }
            )
    return matches, pd.DataFrame(diagnostics)


def null_matrix(
    anchors: pd.DataFrame,
    matches: dict[int, dict[str, np.ndarray]],
    permutations: int,
    seed: int,
) -> tuple[np.ndarray, np.ndarray]:
    eligible_indices = np.asarray(sorted(matches), dtype=int)
    rng = np.random.default_rng(seed)
    matrix = np.empty((len(eligible_indices), permutations), dtype=np.float32)
    for position, row_index in enumerate(eligible_indices):
        candidates = matches[int(row_index)]
        selected = rng.integers(
            0, len(candidates["beta_D"]), size=permutations
        )
        beta_d = candidates["beta_D"][selected]
        se_d = candidates["SE_D"][selected]
        row = anchors.loc[row_index]
        matrix[position] = independent_score(
            np.full(permutations, float(row["beta_G"])),
            np.full(permutations, float(row["SE_G"])),
            beta_d,
            se_d,
        ).astype(np.float32)
    return eligible_indices, matrix


def permutation_group_stat(
    anchors: pd.DataFrame,
    eligible_indices: np.ndarray,
    null_scores: np.ndarray,
    group_indices: np.ndarray,
) -> tuple[np.ndarray, int]:
    selected = np.intersect1d(
        eligible_indices, group_indices, assume_unique=False
    )
    if not len(selected):
        return np.array([], dtype=float), 0
    position = {row_index: i for i, row_index in enumerate(eligible_indices)}
    selected_positions = np.asarray([position[index] for index in selected])
    selected_frame = anchors.loc[selected]
    locus_values: list[np.ndarray] = []
    for _, locus_rows in selected_frame.groupby("locus_id", observed=True):
        local_positions = np.asarray(
            [position[index] for index in locus_rows.index], dtype=int
        )
        locus_values.append(null_scores[local_positions].mean(axis=0))
    return np.vstack(locus_values).mean(axis=0), len(locus_values)


def empirical_two_sided(observed: float, null_values: np.ndarray) -> float:
    center = float(np.mean(null_values))
    distance = abs(observed - center)
    return float(
        (1 + np.count_nonzero(np.abs(null_values - center) >= distance))
        / (len(null_values) + 1)
    )


def main() -> None:
    args = parse_args()
    if args.permutations < 100 or args.candidate_count < 1:
        raise ValueError("Too few permutations or candidates.")
    anchor_path = args.anchor_results.resolve()
    all_gene_dir = args.all_gene_dir.resolve()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    anchors = pd.read_csv(anchor_path, low_memory=False)
    if "disease" not in anchors.columns and "disease.x" in anchors.columns:
        source_disease = anchors["disease.x"]
        if source_disease.isna().any():
            raise RuntimeError(
                "The source anchor disease column disease.x contains missing values."
            )
        if "disease.y" in anchors.columns:
            model_disease = anchors["disease.y"]
            conflicts = (
                model_disease.notna()
                & source_disease.astype(str).ne(model_disease.astype(str))
            )
            if conflicts.any():
                raise RuntimeError(
                    "Conflicting source/model disease labels in anchor results: "
                    f"{int(conflicts.sum())} rows."
                )
        anchors = anchors.rename(columns={"disease.x": "disease"})
    required = {
        "anchor_unit_id",
        "disease",
        "locus_id",
        "gene_id_clean",
        "cell_class",
        "beta_G",
        "SE_G",
        "beta_D",
        "SE_D",
        "expression_estimable",
        "evidence_grade",
        "qtl_resolution",
    }
    missing = sorted(required - set(anchors))
    if missing:
        raise RuntimeError(f"Anchor result columns absent: {missing}")
    anchors = anchors.loc[as_bool(anchors["expression_estimable"])].copy()
    anchors = anchors.loc[anchors["evidence_grade"].isin(["G1", "G2"])].copy()
    anchors["gene_id_clean"] = anchors["gene_id_clean"].astype(str)
    for column in ["beta_G", "SE_G", "beta_D", "SE_D"]:
        anchors[column] = pd.to_numeric(anchors[column], errors="coerce")
    anchors = anchors.loc[
        np.isfinite(anchors["beta_G"])
        & np.isfinite(anchors["SE_G"])
        & anchors["SE_G"].gt(0)
        & np.isfinite(anchors["beta_D"])
        & np.isfinite(anchors["SE_D"])
        & anchors["SE_D"].gt(0)
    ].copy()
    if anchors["anchor_unit_id"].duplicated().any():
        raise RuntimeError("Duplicated anchor_unit_id after model merge.")
    anchors = anchors.reset_index(drop=True)

    anchors["S_moderated_rho_0"] = independent_score(
        anchors["beta_G"].to_numpy(float),
        anchors["SE_G"].to_numpy(float),
        anchors["beta_D"].to_numpy(float),
        anchors["SE_D"].to_numpy(float),
    )
    if {"beta_D_HC3", "SE_D_HC3"} <= set(anchors):
        hc3_beta = pd.to_numeric(anchors["beta_D_HC3"], errors="coerce")
        hc3_se = pd.to_numeric(anchors["SE_D_HC3"], errors="coerce")
        valid_hc3 = np.isfinite(hc3_beta) & np.isfinite(hc3_se) & hc3_se.gt(0)
        anchors["S_HC3_rho_0"] = np.nan
        anchors.loc[valid_hc3, "S_HC3_rho_0"] = independent_score(
            anchors.loc[valid_hc3, "beta_G"].to_numpy(float),
            anchors.loc[valid_hc3, "SE_G"].to_numpy(float),
            hc3_beta.loc[valid_hc3].to_numpy(float),
            hc3_se.loc[valid_hc3].to_numpy(float),
        )
    else:
        anchors["S_HC3_rho_0"] = np.nan

    rho_grid = np.round(np.arange(-0.5, 0.5001, 0.1), 1)
    rho_values = np.concatenate(([-0.75], rho_grid, [0.75]))
    rho_columns: list[str] = []
    for rho in rho_values:
        column = f"S_rho_{rho:+.2f}".replace("+", "p").replace("-", "m")
        rho_columns.append(column)
        anchors[column] = [
            correlated_score(bg, sg, bd, sd, float(rho))
            for bg, sg, bd, sd in anchors[
                ["beta_G", "SE_G", "beta_D", "SE_D"]
            ].itertuples(index=False, name=None)
        ]

    matches, match_diagnostics = build_matches(
        anchors,
        all_gene_dir,
        args.candidate_count,
        args.minimum_candidates,
    )
    eligible_indices, null_scores = null_matrix(
        anchors, matches, args.permutations, args.seed
    )
    match_diagnostics.to_csv(
        output_dir / "matched_null_diagnostics.csv", index=False
    )

    groups = group_definitions(anchors)
    group_rows: list[dict[str, Any]] = []
    group_null: dict[tuple[str, str, str], np.ndarray] = {}
    for group in groups:
        subset = anchors.loc[group["indices"]]
        observed_all, _ = locus_equal_stat(
            subset, "S_moderated_rho_0"
        )
        matched_indices = np.intersect1d(
            eligible_indices, group["indices"], assume_unique=False
        )
        matched_subset = anchors.loc[matched_indices]
        observed = np.nan
        locus_means = pd.Series(dtype=float)
        if len(matched_subset):
            observed, locus_means = locus_equal_stat(
                matched_subset, "S_moderated_rho_0"
            )
        hc3 = np.nan
        if matched_subset["S_HC3_rho_0"].notna().any():
            hc3_subset = matched_subset.dropna(subset=["S_HC3_rho_0"])
            hc3, _ = locus_equal_stat(hc3_subset, "S_HC3_rho_0")
        permutation_values, null_loci = permutation_group_stat(
            anchors,
            eligible_indices,
            null_scores,
            group["indices"],
        )
        p_value = np.nan
        null_mean = np.nan
        null_sd = np.nan
        if null_loci >= args.minimum_loci:
            null_mean = float(np.mean(permutation_values))
            null_sd = float(np.std(permutation_values, ddof=1))
            p_value = empirical_two_sided(observed, permutation_values)
            group_null[
                (group["level"], group["disease"], group["cell_class"])
            ] = permutation_values
        leave_one_values = []
        if len(locus_means) > 1:
            for locus_id in locus_means.index:
                leave_one_values.append(
                    float(locus_means.drop(index=locus_id).mean())
                )
        contribution_denominator = float(np.abs(locus_means).sum())
        max_contribution = (
            float(np.abs(locus_means).max() / contribution_denominator)
            if contribution_denominator > 0
            else np.nan
        )
        group_rows.append(
            {
                "level": group["level"],
                "disease": group["disease"],
                "cell_class": group["cell_class"],
                "anchor_rows": int(len(subset)),
                "matched_anchor_rows": int(len(matched_subset)),
                "loci": int(subset["locus_id"].nunique()),
                "matched_null_loci": int(null_loci),
                "S_locus_equal_moderated_all_anchors": observed_all,
                "S_locus_equal_moderated_matched_anchors": observed,
                "S_locus_equal_HC3_matched_anchors": hc3,
                "matched_null_mean": null_mean,
                "matched_null_sd": null_sd,
                "difference_from_null_mean": observed - null_mean,
                "empirical_p_two_sided": p_value,
                "leave_one_locus_min": (
                    min(leave_one_values) if leave_one_values else np.nan
                ),
                "leave_one_locus_max": (
                    max(leave_one_values) if leave_one_values else np.nan
                ),
                "leave_one_locus_crosses_zero": (
                    min(leave_one_values) <= 0 <= max(leave_one_values)
                    if leave_one_values
                    else np.nan
                ),
                "maximum_absolute_locus_contribution_share": max_contribution,
            }
        )
    group_results = pd.DataFrame(group_rows)
    stratum_mask = group_results["level"].eq("stratum")
    group_results.loc[
        stratum_mask, "empirical_FDR_BH_within_strata"
    ] = bh_adjust(
        group_results.loc[stratum_mask, "empirical_p_two_sided"]
    )
    group_results.to_csv(
        output_dir / "direction_coupling_group_results.csv", index=False
    )

    stratum_groups = [
        group
        for group in groups
        if group["level"] == "stratum"
        and ("stratum", group["disease"], group["cell_class"]) in group_null
    ]
    omnibus_rows: list[dict[str, Any]] = []
    for disease in ["ALL", *sorted(anchors["disease"].unique())]:
        selected_groups = [
            group
            for group in stratum_groups
            if disease == "ALL" or group["disease"] == disease
        ]
        if not selected_groups:
            continue
        observed_values = []
        null_arrays = []
        labels = []
        for group in selected_groups:
            result_row = group_results.loc[
                group_results["level"].eq("stratum")
                & group_results["disease"].eq(group["disease"])
                & group_results["cell_class"].eq(group["cell_class"])
            ].iloc[0]
            observed_values.append(
                result_row["S_locus_equal_moderated_matched_anchors"]
            )
            null_arrays.append(
                group_null[
                    ("stratum", group["disease"], group["cell_class"])
                ]
            )
            labels.append(f"{group['disease']}__{group['cell_class']}")
        null_array = np.vstack(null_arrays)
        centers = null_array.mean(axis=1)
        scales = null_array.std(axis=1, ddof=1)
        valid = np.isfinite(scales) & (scales > 0)
        observed_q = float(
            np.sum(
                (
                    (np.asarray(observed_values)[valid] - centers[valid])
                    / scales[valid]
                )
                ** 2
            )
        )
        null_q = np.sum(
            ((null_array[valid] - centers[valid, None]) / scales[valid, None])
            ** 2,
            axis=0,
        )
        p_value = float(
            (1 + np.count_nonzero(null_q >= observed_q))
            / (args.permutations + 1)
        )
        omnibus_rows.append(
            {
                "scope": "global" if disease == "ALL" else "disease",
                "disease": disease,
                "strata": len(labels),
                "stratum_labels": ";".join(labels),
                "Q_observed": observed_q,
                "empirical_p": p_value,
                "permutations": args.permutations,
            }
        )
    omnibus = pd.DataFrame(omnibus_rows)
    omnibus.to_csv(output_dir / "direction_coupling_omnibus.csv", index=False)

    rho_rows: list[dict[str, Any]] = []
    rho_summary_rows: list[dict[str, Any]] = []
    for group in groups:
        subset = anchors.loc[group["indices"]]
        values = []
        primary_values = []
        for rho, column in zip(rho_values, rho_columns):
            statistic, _ = locus_equal_stat(subset, column)
            values.append(statistic)
            if -0.5 <= rho <= 0.5:
                primary_values.append(statistic)
            rho_rows.append(
                {
                    "level": group["level"],
                    "disease": group["disease"],
                    "cell_class": group["cell_class"],
                    "rho": rho,
                    "grid_type": (
                        "primary" if -0.5 <= rho <= 0.5 else "extreme"
                    ),
                    "S_locus_equal": statistic,
                }
            )
        rho_summary_rows.append(
            {
                "level": group["level"],
                "disease": group["disease"],
                "cell_class": group["cell_class"],
                "primary_grid_min": min(primary_values),
                "primary_grid_max": max(primary_values),
                "primary_grid_crosses_zero": (
                    min(primary_values) <= 0 <= max(primary_values)
                ),
                "extreme_grid_min": min(values),
                "extreme_grid_max": max(values),
                "extreme_grid_crosses_zero": min(values) <= 0 <= max(values),
            }
        )
    pd.DataFrame(rho_rows).to_csv(
        output_dir / "rho_grid_sensitivity.csv", index=False
    )
    pd.DataFrame(rho_summary_rows).to_csv(
        output_dir / "rho_grid_group_summary.csv", index=False
    )

    subset_rows: list[dict[str, Any]] = []
    subset_definitions = {
        "G1_G2_all_mappings": np.ones(len(anchors), dtype=bool),
        "G1_only": anchors["evidence_grade"].eq("G1").to_numpy(),
        "exact_major_class_only": anchors["qtl_resolution"]
        .eq("major_class")
        .to_numpy(),
        "G1_exact_major_class_only": (
            anchors["evidence_grade"].eq("G1")
            & anchors["qtl_resolution"].eq("major_class")
        ).to_numpy(),
    }
    for subset_name, keep in subset_definitions.items():
        subset_anchor_frame = anchors.loc[keep]
        for group in group_definitions(subset_anchor_frame):
            group_frame = subset_anchor_frame.loc[group["indices"]]
            statistic, _ = locus_equal_stat(
                group_frame, "S_moderated_rho_0"
            )
            hc3_statistic = np.nan
            if group_frame["S_HC3_rho_0"].notna().any():
                hc3_statistic, _ = locus_equal_stat(
                    group_frame.dropna(subset=["S_HC3_rho_0"]),
                    "S_HC3_rho_0",
                )
            subset_rows.append(
                {
                    "sensitivity_set": subset_name,
                    "level": group["level"],
                    "disease": group["disease"],
                    "cell_class": group["cell_class"],
                    "anchor_rows": len(group_frame),
                    "loci": group_frame["locus_id"].nunique(),
                    "S_locus_equal_moderated": statistic,
                    "S_locus_equal_HC3": hc3_statistic,
                }
            )
    pd.DataFrame(subset_rows).to_csv(
        output_dir / "direction_coupling_subset_sensitivity.csv",
        index=False,
    )

    anchor_output = output_dir / "anchor_direction_probabilities.csv.gz"
    anchors.to_csv(anchor_output, index=False, compression="gzip")
    outputs = [
        output_dir / "matched_null_diagnostics.csv",
        output_dir / "direction_coupling_group_results.csv",
        output_dir / "direction_coupling_omnibus.csv",
        output_dir / "rho_grid_sensitivity.csv",
        output_dir / "rho_grid_group_summary.csv",
        output_dir / "direction_coupling_subset_sensitivity.csv",
        anchor_output,
    ]
    manifest = {
        "status": "COMPLETE",
        "analysis_freeze": "MSSM_DIRECTION_COUPLING_ANALYSIS_FREEZE_2026-07-30.md",
        "anchor_results": str(anchor_path),
        "anchor_results_sha256": sha256(anchor_path),
        "anchor_rows": int(len(anchors)),
        "loci": int(anchors["locus_id"].nunique()),
        "diseases": sorted(anchors["disease"].unique()),
        "permutations": args.permutations,
        "seed": args.seed,
        "candidate_count": args.candidate_count,
        "minimum_candidates": args.minimum_candidates,
        "matched_rows": len(matches),
        "rho_primary_grid": rho_grid.tolist(),
        "rho_extreme_grid": [-0.75, 0.75],
        "interpretation_boundary": (
            "Opposed direction is not evidence of protection or compensation."
        ),
        "output_sha256": {path.name: sha256(path) for path in outputs},
        "completed_at": now(),
    }
    atomic_json(output_dir / "direction_coupling_manifest.json", manifest)
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
