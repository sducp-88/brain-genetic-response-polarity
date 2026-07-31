#!/usr/bin/env python3
"""Freeze MSSM AD/SCZ donor definitions before expression-effect estimation."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


DISEASES = ("AD", "SCZ")
CONTROL_LABEL = "CTRL"
REQUIRED_COLUMNS = (
    "DonorID",
    "Cohort",
    "Age",
    "Sex",
    "Ancestry",
    "PMI",
    "Diagnosis",
    "Tier1_crossDis",
    "Tier1_crossDis_dx",
)
MODEL_COVARIATES = ("Age", "Sex", "Ancestry_harmonized", "PMI")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("metadata", type=Path)
    parser.add_argument("--cohort", default="MSSM")
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("outputs/phase2/mssm_sample_freeze"),
    )
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
    if isinstance(value, float) and not np.isfinite(value):
        return None
    if isinstance(value, dict):
        return {str(key): json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(item) for item in value]
    return value


def clean_text(series: pd.Series) -> pd.Series:
    return series.astype("string").str.strip().replace("", pd.NA)


def count_topcoded(series: pd.Series) -> int:
    return int(clean_text(series).str.endswith("+", na=False).sum())


def continuous_smd(case: pd.Series, control: pd.Series) -> float:
    case_values = pd.to_numeric(case, errors="coerce").dropna().to_numpy(float)
    control_values = (
        pd.to_numeric(control, errors="coerce").dropna().to_numpy(float)
    )
    pooled_sd = np.sqrt(
        (np.var(case_values, ddof=1) + np.var(control_values, ddof=1)) / 2
    )
    if not np.isfinite(pooled_sd) or pooled_sd == 0:
        return 0.0
    return float((np.mean(case_values) - np.mean(control_values)) / pooled_sd)


def categorical_level_rows(
    frame: pd.DataFrame, disease: str, variable: str
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    levels = sorted(frame[variable].dropna().astype(str).unique())
    for level in levels:
        case = frame.loc[frame["analysis_role"].eq("case"), variable].astype(str)
        control = frame.loc[
            frame["analysis_role"].eq("control"), variable
        ].astype(str)
        case_p = float(case.eq(level).mean())
        control_p = float(control.eq(level).mean())
        pooled = (case_p + control_p) / 2
        denominator = np.sqrt(pooled * (1 - pooled))
        smd = 0.0 if denominator == 0 else (case_p - control_p) / denominator
        rows.append(
            {
                "disease": disease,
                "variable": variable,
                "level": level,
                "case_n": int(case.eq(level).sum()),
                "control_n": int(control.eq(level).sum()),
                "case_proportion": case_p,
                "control_proportion": control_p,
                "standardized_difference": float(smd),
            }
        )
    return rows


def main() -> None:
    args = parse_args()
    metadata_path = args.metadata.resolve()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    metadata = pd.read_csv(metadata_path, dtype="string")
    missing_columns = [field for field in REQUIRED_COLUMNS if field not in metadata]
    if missing_columns:
        raise RuntimeError(f"Required metadata columns absent: {missing_columns}")

    cohort = metadata.loc[metadata["Cohort"].eq(args.cohort)].copy()
    if cohort.empty:
        raise RuntimeError(f"No donors found for cohort {args.cohort}.")
    if cohort["DonorID"].duplicated().any():
        raise RuntimeError("DonorID is not unique within cohort metadata.")

    for field in REQUIRED_COLUMNS:
        cohort[field] = clean_text(cohort[field])
    cohort["age_topcoded"] = cohort["Age"].str.endswith("+", na=False)
    cohort["age_numeric"] = pd.to_numeric(
        cohort["Age"].str.replace("+", "", regex=False), errors="coerce"
    )
    cohort["PMI_numeric"] = pd.to_numeric(cohort["PMI"], errors="coerce")
    cohort["ancestry_harmonized"] = cohort["Ancestry"].replace(
        {"EAS": "EAS_SAS", "SAS": "EAS_SAS"}
    )

    eligibility_frames: list[pd.DataFrame] = []
    summary_rows: list[dict[str, Any]] = []
    covariate_rows: list[dict[str, Any]] = []
    supports: dict[str, dict[str, float]] = {}

    for disease in DISEASES:
        frame = cohort.copy()
        official_tier1 = frame["Tier1_crossDis"].eq("Y").fillna(False)
        target_label = (
            frame["Tier1_crossDis_dx"].eq(disease).fillna(False)
        )
        control_label = (
            frame["Tier1_crossDis_dx"].eq(CONTROL_LABEL).fillna(False)
        )
        frame["analysis_role"] = np.select(
            [official_tier1 & target_label, official_tier1 & control_label],
            ["case", "control"],
            default="excluded",
        )
        candidates = frame.loc[
            frame["analysis_role"].isin(["case", "control"])
        ].copy()
        case_ages = candidates.loc[
            candidates["analysis_role"].eq("case"), "age_numeric"
        ].dropna()
        control_ages = candidates.loc[
            candidates["analysis_role"].eq("control"), "age_numeric"
        ].dropna()
        if case_ages.empty or control_ages.empty:
            raise RuntimeError(f"Missing case or control ages for {disease}.")
        support_lower = float(max(case_ages.min(), control_ages.min()))
        support_upper = float(min(case_ages.max(), control_ages.max()))
        if support_lower > support_upper:
            raise RuntimeError(f"No age common support for {disease}.")
        supports[disease] = {
            "lower": support_lower,
            "upper": support_upper,
        }

        frame["within_fixed_age_support"] = frame["age_numeric"].between(
            support_lower, support_upper, inclusive="both"
        )
        frame["model_covariates_complete"] = (
            frame["age_numeric"].notna()
            & frame["Sex"].notna()
            & frame["ancestry_harmonized"].notna()
            & frame["PMI_numeric"].notna()
        )
        pre_ancestry = frame.loc[
            frame["analysis_role"].isin(["case", "control"])
            & frame["within_fixed_age_support"]
            & frame["model_covariates_complete"]
        ]
        case_ancestries = set(
            pre_ancestry.loc[
                pre_ancestry["analysis_role"].eq("case"),
                "ancestry_harmonized",
            ].astype(str)
        )
        control_ancestries = set(
            pre_ancestry.loc[
                pre_ancestry["analysis_role"].eq("control"),
                "ancestry_harmonized",
            ].astype(str)
        )
        shared_ancestries = sorted(case_ancestries & control_ancestries)
        frame["within_ancestry_common_support"] = (
            frame["ancestry_harmonized"].isin(shared_ancestries)
        )
        frame["primary_donor_eligible_before_cell_filter"] = (
            frame["analysis_role"].isin(["case", "control"])
            & frame["within_fixed_age_support"]
            & frame["model_covariates_complete"]
            & frame["within_ancestry_common_support"]
        )
        supports[disease]["ancestry_levels"] = shared_ancestries

        exclusion_reason = np.full(len(frame), "eligible", dtype=object)
        exclusion_reason[~official_tier1.to_numpy()] = "not_official_tier1"
        other_tier1 = (
            official_tier1
            & ~target_label
            & ~control_label
        )
        exclusion_reason[other_tier1.to_numpy()] = (
            "other_tier1_diagnosis_or_axis"
        )
        candidate_mask = frame["analysis_role"].isin(["case", "control"])
        exclusion_reason[
            (candidate_mask & ~frame["within_fixed_age_support"]).to_numpy()
        ] = "outside_fixed_age_common_support"
        exclusion_reason[
            (
                candidate_mask
                & frame["within_fixed_age_support"]
                & ~frame["model_covariates_complete"]
            ).to_numpy()
        ] = "missing_prespecified_covariate"
        exclusion_reason[
            (
                candidate_mask
                & frame["within_fixed_age_support"]
                & frame["model_covariates_complete"]
                & ~frame["within_ancestry_common_support"]
            ).to_numpy()
        ] = "outside_ancestry_common_support"
        frame["exclusion_reason"] = exclusion_reason
        frame.insert(0, "target_disease", disease)
        eligibility_frames.append(frame)

        analytic = frame.loc[
            frame["primary_donor_eligible_before_cell_filter"]
        ].copy()
        for role in ("case", "control"):
            source = candidates.loc[candidates["analysis_role"].eq(role)]
            selected = analytic.loc[analytic["analysis_role"].eq(role)]
            summary_rows.append(
                {
                    "cohort": args.cohort,
                    "disease": disease,
                    "role": role,
                    "official_tier1_candidates": int(len(source)),
                    "age_support_lower": support_lower,
                    "age_support_upper": support_upper,
                    "eligible_before_cell_filter": int(len(selected)),
                    "retained_percent": (
                        100 * len(selected) / len(source) if len(source) else np.nan
                    ),
                    "topcoded_age_n_eligible": count_topcoded(selected["Age"]),
                    "missing_age_candidates": int(
                        source["age_numeric"].isna().sum()
                    ),
                    "missing_sex_candidates": int(source["Sex"].isna().sum()),
                    "missing_ancestry_candidates": int(
                        source["ancestry_harmonized"].isna().sum()
                    ),
                    "missing_PMI_candidates": int(
                        source["PMI_numeric"].isna().sum()
                    ),
                }
            )

        cases = analytic.loc[analytic["analysis_role"].eq("case")]
        controls = analytic.loc[analytic["analysis_role"].eq("control")]
        for variable, column in (("age", "age_numeric"), ("PMI", "PMI_numeric")):
            covariate_rows.append(
                {
                    "disease": disease,
                    "variable": variable,
                    "level": "",
                    "case_n": int(cases[column].notna().sum()),
                    "control_n": int(controls[column].notna().sum()),
                    "case_proportion": np.nan,
                    "control_proportion": np.nan,
                    "standardized_difference": continuous_smd(
                        cases[column], controls[column]
                    ),
                }
            )
        covariate_rows.extend(categorical_level_rows(analytic, disease, "Sex"))
        covariate_rows.extend(
            categorical_level_rows(analytic, disease, "ancestry_harmonized")
        )

    eligibility = pd.concat(eligibility_frames, ignore_index=True)
    summary = pd.DataFrame(summary_rows)
    covariates = pd.DataFrame(covariate_rows)

    eligibility_path = output_dir / "mssm_donor_analysis_eligibility.csv"
    summary_path = output_dir / "mssm_disease_sample_summary.csv"
    covariates_path = output_dir / "mssm_covariate_balance_unweighted.csv"
    eligibility.to_csv(eligibility_path, index=False)
    summary.to_csv(summary_path, index=False)
    covariates.to_csv(covariates_path, index=False)

    freeze = {
        "status": "FROZEN_BEFORE_DISEASE_EXPRESSION_EFFECTS",
        "version": "1.1",
        "cohort": args.cohort,
        "diseases": list(DISEASES),
        "source_metadata": str(metadata_path),
        "source_metadata_sha256": sha256(metadata_path),
        "case_definition": (
            "Tier1_crossDis == 'Y' and Tier1_crossDis_dx exactly equals "
            "the target disease (AD or SCZ)"
        ),
        "control_definition": (
            "Tier1_crossDis == 'Y' and Tier1_crossDis_dx exactly equals 'CTRL'"
        ),
        "comorbidity_rule": (
            "All non-target Tier-1 labels and all donors outside the official "
            "clean target/control labels are excluded from the primary model"
        ),
        "age_rule": (
            "Restrict each disease comparison to the fixed cohort-level "
            "case-control age-range intersection; public 89+ is coded as 89"
        ),
        "ancestry_rule": (
            "Harmonize EAS, SAS and EAS_SAS as EAS_SAS, then restrict the "
            "primary comparison to ancestry levels represented in both cases "
            "and controls before estimating overlap weights"
        ),
        "topcoded_age_sensitivity": (
            "Repeat key results excluding age-topcoded donors and, separately, "
            "with an age-topcoded indicator when estimable"
        ),
        "model_covariates_required": list(MODEL_COVARIATES),
        "cell_filter_rule": (
            "After pseudobulk construction, require >=20 nuclei for each "
            "donor-by-cell-class sample; report 10/20/50 sensitivity"
        ),
        "disease_cell_minimum": (
            "A confirmatory disease-by-cell-class result requires >=30 eligible "
            "cases and >=30 eligible controls before overlap-weight ESS checks"
        ),
        "common_age_support": supports,
        "propensity_and_balance_rule": (
            "Estimate disease-by-cell-class overlap weights without viewing "
            "gene effects; require absolute SMD <=0.10 for key covariates and "
            "report group-specific effective sample sizes"
        ),
        "primary_effect_model": (
            "donor-level pseudobulk; disease effect adjusted for nonlinear age, "
            "sex, ancestry, PMI, and prespecified library-quality terms; "
            "brain banks analyzed separately"
        ),
        "missing_unavailable_covariates": (
            "Public H5AD/metadata do not provide complete RIN, brain pH, "
            "medication, agonal state, or detailed technical batch; this is a "
            "design limitation and cannot be statistically erased"
        ),
        "output_sha256": {
            eligibility_path.name: sha256(eligibility_path),
            summary_path.name: sha256(summary_path),
            covariates_path.name: sha256(covariates_path),
        },
        "created_at": pd.Timestamp.now(tz="Asia/Shanghai").isoformat(),
    }
    freeze_path = output_dir / "mssm_sample_definition_freeze.json"
    freeze_path.write_text(
        json.dumps(json_safe(freeze), ensure_ascii=False, indent=2),
        encoding="utf-8",
    )

    report_lines = [
        "# MSSM AD/SCZ 正式样本定义冻结报告",
        "",
        f"- 状态：`{freeze['status']}`",
        f"- 版本：`{freeze['version']}`",
        f"- 元数据 SHA-256：`{freeze['source_metadata_sha256']}`",
        "- 本报告在查看 MSSM 疾病表达效应前生成；后续不得按结果修改主要纳排规则。",
        "",
        "## 主要病例与对照",
        "",
        "- AD 病例：官方 `Tier1_crossDis=Y` 且 `Tier1_crossDis_dx=AD`。",
        "- SCZ 病例：官方 `Tier1_crossDis=Y` 且 `Tier1_crossDis_dx=SCZ`。",
        "- 共同对照：官方 `Tier1_crossDis=Y` 且 `Tier1_crossDis_dx=CTRL`。",
        "- 共病、多标签、其他神经退行性或精神疾病均不进入主要病例—对照模型。",
        "- `EAS`、`SAS`、`EAS_SAS` 先统一为 `EAS_SAS`；主分析仅保留病例与对照均出现的祖源层级。",
        "",
        "## 年龄共同支持域与样本量",
        "",
        "| 疾病 | 角色 | Tier-1候选 | 年龄支持域 | 细胞过滤前合格 | 保留率 |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for row in summary.to_dict(orient="records"):
        report_lines.append(
            f"| {row['disease']} | {row['role']} | "
            f"{row['official_tier1_candidates']} | "
            f"{row['age_support_lower']:.0f}–{row['age_support_upper']:.0f} | "
            f"{row['eligible_before_cell_filter']} | "
            f"{row['retained_percent']:.1f}% |"
        )
    report_lines.extend(
        [
            "",
            "## 在查看基因效应前冻结的分析门槛",
            "",
            "- 主阈值：每名供体每个细胞大类至少 20 个细胞核；同时报告 10/20/50 敏感性。",
            "- 确认性疾病×细胞大类：至少 30 例病例和 30 例对照，随后还需通过加权后有效样本量检查。",
            "- 每个疾病×细胞大类独立估计 overlap weights；关键协变量绝对 SMD 必须 ≤0.10。",
            "- 主效应模型使用供体级 pseudobulk，并校正非线性年龄、性别、祖源、PMI 与预设库质量项。",
            "- `89+` 在主分析按公开数据的 89 处理；关键结果另做排除 `89+` 和加入 top-code 指示变量的敏感性。",
            "",
            "## 不能消除的限制",
            "",
            "公开数据缺少完整 RIN、脑 pH、药物暴露、濒死状态和部分技术批次信息。"
            "这些属于残余混杂限制，不能以模型复杂度替代。",
            "",
            "## 机器可读文件",
            "",
            "- `mssm_donor_analysis_eligibility.csv`：每个供体在 AD/SCZ 轴上的角色与排除原因。",
            "- `mssm_disease_sample_summary.csv`：病例/对照、支持域与保留率。",
            "- `mssm_covariate_balance_unweighted.csv`：未加权协变量差异；正式权重尚未估计。",
            "- `mssm_sample_definition_freeze.json`：冻结规则、哈希和版本。",
            "",
        ]
    )
    report_path = output_dir / "mssm_sample_definition_freeze_report.md"
    report_path.write_text("\n".join(report_lines), encoding="utf-8")

    print(summary.to_string(index=False))
    print(f"Outputs: {output_dir}")


if __name__ == "__main__":
    main()
