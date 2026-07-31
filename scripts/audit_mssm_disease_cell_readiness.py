#!/usr/bin/env python3
"""Audit MSSM disease-by-cell-class eligibility and overlap weighting.

This script only uses donor/sample metadata and library summaries. It does not
read gene-level expression effects, so it can be run before disease modeling.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample-metadata", type=Path, required=True)
    parser.add_argument("--eligibility", type=Path, required=True)
    parser.add_argument("--freeze-json", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--primary-min-cells", type=int, default=20)
    parser.add_argument(
        "--thresholds",
        type=int,
        nargs="+",
        default=[10, 20, 50],
    )
    parser.add_argument("--minimum-group-n", type=int, default=30)
    parser.add_argument("--maximum-absolute-smd", type=float, default=0.10)
    parser.add_argument("--minimum-group-ess", type=float, default=30.0)
    return parser.parse_args()


def sha256(path: Path, block_size: int = 16 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(block_size):
            digest.update(block)
    return digest.hexdigest()


def json_safe(value: Any) -> Any:
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, float) and not math.isfinite(value):
        return None
    if isinstance(value, dict):
        return {str(key): json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(item) for item in value]
    return value


def weighted_mean(values: np.ndarray, weights: np.ndarray) -> float:
    return float(np.sum(values * weights) / np.sum(weights))


def effective_sample_size(weights: np.ndarray) -> float:
    return float(np.sum(weights) ** 2 / np.sum(weights**2))


def continuous_smd(
    values: np.ndarray, case: np.ndarray, weights: np.ndarray
) -> float:
    case_values = values[case]
    control_values = values[~case]
    case_weights = weights[case]
    control_weights = weights[~case]
    mean_case = weighted_mean(case_values, case_weights)
    mean_control = weighted_mean(control_values, control_weights)
    unweighted_pooled_sd = np.sqrt(
        (np.var(case_values, ddof=1) + np.var(control_values, ddof=1)) / 2
    )
    if not np.isfinite(unweighted_pooled_sd) or unweighted_pooled_sd == 0:
        return 0.0
    return float((mean_case - mean_control) / unweighted_pooled_sd)


def categorical_smd_rows(
    frame: pd.DataFrame,
    variable: str,
    case: np.ndarray,
    weights: np.ndarray,
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    values = frame[variable].astype(str).to_numpy()
    for level in sorted(np.unique(values)):
        indicator = values == level
        p_case = weighted_mean(indicator[case].astype(float), weights[case])
        p_control = weighted_mean(
            indicator[~case].astype(float), weights[~case]
        )
        pooled = (p_case + p_control) / 2
        denominator = np.sqrt(pooled * (1 - pooled))
        smd = 0.0 if denominator == 0 else (p_case - p_control) / denominator
        rows.append(
            {
                "variable": variable,
                "level": level,
                "weighted_case_mean_or_proportion": p_case,
                "weighted_control_mean_or_proportion": p_control,
                "weighted_smd": float(smd),
            }
        )
    return rows


def fit_overlap_weights(
    frame: pd.DataFrame,
) -> tuple[pd.DataFrame, list[dict[str, Any]], dict[str, Any]]:
    model_frame = frame.copy()
    model_frame["case"] = model_frame["analysis_role"].eq("case").astype(int)
    model_frame["age_numeric"] = pd.to_numeric(
        model_frame["age_numeric"], errors="raise"
    )
    model_frame["PMI_numeric"] = pd.to_numeric(
        model_frame["PMI_numeric"], errors="raise"
    )
    formula = (
        "case ~ bs(age_numeric, df=3, degree=3, include_intercept=False) "
        "+ C(Sex) + C(ancestry_harmonized) + PMI_numeric"
    )
    fitted = smf.glm(
        formula=formula,
        data=model_frame,
        family=sm.families.Binomial(),
    ).fit(maxiter=200)
    propensity = np.clip(
        fitted.predict(model_frame).to_numpy(float), 1e-6, 1 - 1e-6
    )
    case = model_frame["case"].to_numpy(bool)
    overlap_weight = np.where(case, 1 - propensity, propensity)
    model_frame["propensity"] = propensity
    model_frame["overlap_weight"] = overlap_weight

    balance_rows: list[dict[str, Any]] = []
    for variable in ("age_numeric", "PMI_numeric"):
        values = model_frame[variable].to_numpy(float)
        balance_rows.append(
            {
                "variable": variable,
                "level": "",
                "weighted_case_mean_or_proportion": weighted_mean(
                    values[case], overlap_weight[case]
                ),
                "weighted_control_mean_or_proportion": weighted_mean(
                    values[~case], overlap_weight[~case]
                ),
                "weighted_smd": continuous_smd(
                    values, case, overlap_weight
                ),
            }
        )
    balance_rows.extend(
        categorical_smd_rows(model_frame, "Sex", case, overlap_weight)
    )
    balance_rows.extend(
        categorical_smd_rows(
            model_frame, "ancestry_harmonized", case, overlap_weight
        )
    )

    diagnostics = {
        "formula": formula,
        "converged": bool(fitted.converged),
        "iterations": int(fitted.fit_history.get("iteration", -1)),
        "propensity_min": float(propensity.min()),
        "propensity_max": float(propensity.max()),
        "case_weight_sum": float(overlap_weight[case].sum()),
        "control_weight_sum": float(overlap_weight[~case].sum()),
        "case_ess": effective_sample_size(overlap_weight[case]),
        "control_ess": effective_sample_size(overlap_weight[~case]),
    }
    return model_frame, balance_rows, diagnostics


def main() -> None:
    args = parse_args()
    sample_metadata_path = args.sample_metadata.resolve()
    eligibility_path = args.eligibility.resolve()
    freeze_path = args.freeze_json.resolve()
    output_dir = args.output_dir.resolve()
    weights_dir = output_dir / "weights"
    output_dir.mkdir(parents=True, exist_ok=True)
    weights_dir.mkdir(parents=True, exist_ok=True)

    thresholds = sorted(set(args.thresholds))
    if args.primary_min_cells not in thresholds:
        raise ValueError("Primary minimum-cell threshold must be in thresholds.")

    sample_metadata = pd.read_csv(sample_metadata_path, low_memory=False)
    eligibility = pd.read_csv(eligibility_path, low_memory=False)
    freeze = json.loads(freeze_path.read_text(encoding="utf-8"))
    if freeze["status"] != "FROZEN_BEFORE_DISEASE_EXPRESSION_EFFECTS":
        raise RuntimeError("Sample definition is not marked as frozen.")
    if sample_metadata["donor_id"].isna().any():
        raise RuntimeError("Missing donor IDs in pseudobulk sample metadata.")
    if sample_metadata.duplicated(["donor_id", "class"]).any():
        raise RuntimeError("Duplicate donor-by-class sample metadata rows.")

    eligible = eligibility.loc[
        eligibility["primary_donor_eligible_before_cell_filter"]
        .astype(str)
        .eq("True")
    ].copy()
    keep_columns = [
        "target_disease",
        "DonorID",
        "analysis_role",
        "Age",
        "age_topcoded",
        "age_numeric",
        "Sex",
        "Ancestry",
        "ancestry_harmonized",
        "PMI",
        "PMI_numeric",
    ]
    eligible = eligible[keep_columns]
    sample_keep_columns = [
        "donor_id",
        "class",
        "n_cells",
        "total_counts",
        "detected_genes",
        "eligible_min_cells",
    ]
    missing_sample_columns = [
        column
        for column in sample_keep_columns
        if column not in sample_metadata
    ]
    if missing_sample_columns:
        raise RuntimeError(
            f"Required pseudobulk sample columns absent: "
            f"{missing_sample_columns}"
        )
    merged = sample_metadata[sample_keep_columns].merge(
        eligible,
        how="inner",
        left_on="donor_id",
        right_on="DonorID",
        validate="many_to_many",
    )
    expected_rows = len(sample_metadata["class"].unique()) * len(eligible)
    if len(merged) != expected_rows:
        raise RuntimeError(
            f"Unexpected merged sample rows: {len(merged)} != {expected_rows}"
        )

    count_rows: list[dict[str, Any]] = []
    readiness_rows: list[dict[str, Any]] = []
    balance_rows_all: list[dict[str, Any]] = []
    weight_hashes: dict[str, str] = {}

    for disease in sorted(eligible["target_disease"].unique()):
        disease_frame = merged.loc[merged["target_disease"].eq(disease)].copy()
        for class_name in sorted(disease_frame["class"].unique()):
            class_frame = disease_frame.loc[
                disease_frame["class"].eq(class_name)
            ].copy()
            for threshold in thresholds:
                threshold_frame = class_frame.loc[
                    pd.to_numeric(class_frame["n_cells"], errors="raise")
                    >= threshold
                ]
                case_n = int(
                    threshold_frame["analysis_role"].eq("case").sum()
                )
                control_n = int(
                    threshold_frame["analysis_role"].eq("control").sum()
                )
                count_rows.append(
                    {
                        "disease": disease,
                        "class": class_name,
                        "minimum_cells": threshold,
                        "cases": case_n,
                        "controls": control_n,
                        "total": int(len(threshold_frame)),
                        "passes_minimum_group_n": bool(
                            case_n >= args.minimum_group_n
                            and control_n >= args.minimum_group_n
                        ),
                    }
                )

            primary = class_frame.loc[
                pd.to_numeric(class_frame["n_cells"], errors="raise")
                >= args.primary_min_cells
            ].copy()
            case_n = int(primary["analysis_role"].eq("case").sum())
            control_n = int(primary["analysis_role"].eq("control").sum())
            enough_n = (
                case_n >= args.minimum_group_n
                and control_n >= args.minimum_group_n
            )
            readiness: dict[str, Any] = {
                "disease": disease,
                "class": class_name,
                "minimum_cells": args.primary_min_cells,
                "cases": case_n,
                "controls": control_n,
                "minimum_group_n_required": args.minimum_group_n,
                "passes_minimum_group_n": enough_n,
                "propensity_status": "not_run_insufficient_n",
                "max_absolute_weighted_smd": np.nan,
                "case_ess": np.nan,
                "control_ess": np.nan,
                "passes_balance": False,
                "passes_ess": False,
                "confirmatory_ready": False,
            }
            if enough_n:
                weighted, balance_rows, diagnostics = fit_overlap_weights(
                    primary
                )
                for row in balance_rows:
                    row["disease"] = disease
                    row["class"] = class_name
                    balance_rows_all.append(row)
                maximum_smd = max(
                    abs(float(row["weighted_smd"])) for row in balance_rows
                )
                balance_pass = maximum_smd <= args.maximum_absolute_smd
                ess_pass = (
                    diagnostics["case_ess"] >= args.minimum_group_ess
                    and diagnostics["control_ess"] >= args.minimum_group_ess
                )
                readiness.update(
                    {
                        "propensity_status": (
                            "converged"
                            if diagnostics["converged"]
                            else "not_converged"
                        ),
                        "propensity_min": diagnostics["propensity_min"],
                        "propensity_max": diagnostics["propensity_max"],
                        "max_absolute_weighted_smd": maximum_smd,
                        "case_ess": diagnostics["case_ess"],
                        "control_ess": diagnostics["control_ess"],
                        "passes_balance": balance_pass,
                        "passes_ess": ess_pass,
                        "confirmatory_ready": bool(
                            diagnostics["converged"]
                            and balance_pass
                            and ess_pass
                        ),
                    }
                )
                weight_columns = [
                    "target_disease",
                    "class",
                    "donor_id",
                    "analysis_role",
                    "n_cells",
                    "total_counts",
                    "Age",
                    "age_topcoded",
                    "age_numeric",
                    "Sex",
                    "Ancestry",
                    "ancestry_harmonized",
                    "PMI",
                    "PMI_numeric",
                    "propensity",
                    "overlap_weight",
                ]
                weight_path = (
                    weights_dir
                    / f"{disease}_{class_name}_overlap_weights.csv"
                )
                weighted[weight_columns].to_csv(weight_path, index=False)
                weight_hashes[weight_path.name] = sha256(weight_path)
            readiness_rows.append(readiness)

    counts = pd.DataFrame(count_rows)
    readiness = pd.DataFrame(readiness_rows)
    balance = pd.DataFrame(balance_rows_all)
    counts_path = output_dir / "disease_cell_threshold_counts.csv"
    readiness_path = output_dir / "disease_cell_readiness.csv"
    balance_path = output_dir / "overlap_weight_balance.csv"
    counts.to_csv(counts_path, index=False)
    readiness.to_csv(readiness_path, index=False)
    balance.to_csv(balance_path, index=False)

    all_ready = bool(readiness["confirmatory_ready"].all())
    manifest = {
        "status": (
            "ALL_DISEASE_CELL_CLASSES_READY"
            if all_ready
            else "ONE_OR_MORE_DISEASE_CELL_CLASSES_NOT_CONFIRMATORY_READY"
        ),
        "expression_effects_read": False,
        "sample_metadata": str(sample_metadata_path),
        "sample_metadata_sha256": sha256(sample_metadata_path),
        "eligibility": str(eligibility_path),
        "eligibility_sha256": sha256(eligibility_path),
        "freeze_json": str(freeze_path),
        "freeze_json_sha256": sha256(freeze_path),
        "thresholds": thresholds,
        "primary_min_cells": args.primary_min_cells,
        "minimum_group_n": args.minimum_group_n,
        "maximum_absolute_smd": args.maximum_absolute_smd,
        "minimum_group_ess": args.minimum_group_ess,
        "confirmatory_ready_classes": int(
            readiness["confirmatory_ready"].sum()
        ),
        "total_disease_cell_classes": int(len(readiness)),
        "output_sha256": {
            counts_path.name: sha256(counts_path),
            readiness_path.name: sha256(readiness_path),
            balance_path.name: sha256(balance_path),
            **weight_hashes,
        },
        "created_at": pd.Timestamp.now(tz="Asia/Shanghai").isoformat(),
    }
    manifest_path = output_dir / "mssm_disease_cell_readiness_manifest.json"
    manifest_path.write_text(
        json.dumps(json_safe(manifest), ensure_ascii=False, indent=2),
        encoding="utf-8",
    )

    report_lines = [
        "# MSSM 疾病×细胞大类正式建模就绪审计",
        "",
        f"- 状态：`{manifest['status']}`",
        "- 审计只读取供体标签、细胞数和库大小；未读取任何基因疾病效应。",
        f"- 主细胞阈值：≥{args.primary_min_cells} 个细胞核/供体/大类。",
        f"- 最低病例/对照数：各 ≥{args.minimum_group_n}。",
        f"- 加权平衡：最大绝对 SMD ≤{args.maximum_absolute_smd:.2f}。",
        f"- 加权有效样本量：病例和对照各 ≥{args.minimum_group_ess:.0f}。",
        "",
        "| 疾病 | 细胞大类 | 病例 | 对照 | 最大|SMD| | 病例ESS | 对照ESS | 可进入确认性模型 |",
        "|---|---|---:|---:|---:|---:|---:|---|",
    ]
    for row in readiness.to_dict(orient="records"):
        maximum_smd = (
            "NA"
            if pd.isna(row["max_absolute_weighted_smd"])
            else f"{row['max_absolute_weighted_smd']:.3f}"
        )
        case_ess = (
            "NA" if pd.isna(row["case_ess"]) else f"{row['case_ess']:.1f}"
        )
        control_ess = (
            "NA"
            if pd.isna(row["control_ess"])
            else f"{row['control_ess']:.1f}"
        )
        report_lines.append(
            f"| {row['disease']} | {row['class']} | {row['cases']} | "
            f"{row['controls']} | {maximum_smd} | {case_ess} | "
            f"{control_ess} | "
            f"{'是' if row['confirmatory_ready'] else '否'} |"
        )
    report_lines.extend(
        [
            "",
            "未就绪的大类可以保留为探索性或提高细胞阈值/合并层级后的敏感性，"
            "但不得进入最强确认性方向图谱。",
            "",
        ]
    )
    (output_dir / "mssm_disease_cell_readiness_report.md").write_text(
        "\n".join(report_lines), encoding="utf-8"
    )
    print(readiness.to_string(index=False))
    print(f"Outputs: {output_dir}")


if __name__ == "__main__":
    main()
