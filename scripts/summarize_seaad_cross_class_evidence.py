#!/usr/bin/env python3
"""Create an auditable decision summary across SEA-AD validation cell classes."""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd


CLASSES = ("Immune", "Oligo", "EN")
PHENOTYPES = ("CERAD", "Braak")


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


def read_global(path: Path) -> pd.Series:
    table = pd.read_csv(path)
    row = table.loc[table["level"].eq("global")]
    if len(row) != 1:
        raise RuntimeError(f"Expected one global row in {path}; found {len(row)}")
    return row.iloc[0]


def as_bool(value: object) -> bool:
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    return str(value).strip().lower() == "true"


def markdown_table(frame: pd.DataFrame) -> str:
    """Render a small DataFrame without pandas' optional tabulate dependency."""
    columns = [str(column) for column in frame.columns]
    rows = [
        [
            str(value).replace("|", r"\|").replace("\n", " ")
            for value in row
        ]
        for row in frame.itertuples(index=False, name=None)
    ]
    header = "| " + " | ".join(columns) + " |"
    separator = "| " + " | ".join("---" for _ in columns) + " |"
    body = [
        "| " + " | ".join(row) + " |"
        for row in rows
    ]
    return "\n".join([header, separator, *body])


def main() -> None:
    args = parse_args()
    project = args.project.resolve()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    records: list[dict[str, object]] = []
    roots = {
        "Immune": project / "outputs/phase5/seaad_immune_core3",
        "Oligo": project / "outputs/phase5/seaad_oligo_core3",
        "EN": project / "outputs/phase5/seaad_en_core3",
    }
    for cell_class in CLASSES:
        for phenotype in PHENOTYPES:
            path = (
                roots[cell_class]
                / phenotype
                / "coupling"
                / "direction_coupling_group_results.csv"
            )
            row = read_global(path)
            records.append(
                {
                    "cohort": "SEA-AD",
                    "cell_class": cell_class,
                    "phenotype": phenotype,
                    "loci": int(row["loci"]),
                    "S": float(row["S_locus_equal_moderated_all_anchors"]),
                    "matched_null_mean": float(row["matched_null_mean"]),
                    "difference_from_null_mean": float(
                        row["difference_from_null_mean"]
                    ),
                    "empirical_p_two_sided": float(
                        row["empirical_p_two_sided"]
                    ),
                    "leave_one_locus_min": float(row["leave_one_locus_min"]),
                    "leave_one_locus_max": float(row["leave_one_locus_max"]),
                    "leave_one_locus_crosses_zero": as_bool(
                        row["leave_one_locus_crosses_zero"]
                    ),
                    "maximum_absolute_locus_contribution_share": float(
                        row["maximum_absolute_locus_contribution_share"]
                    ),
                    "source_file": str(path),
                }
            )
    evidence = pd.DataFrame.from_records(records)
    evidence["empirical_FDR_BH_six_external_tests"] = bh_adjust(
        evidence["empirical_p_two_sided"]
    )
    evidence.to_csv(output_dir / "seaad_cross_class_evidence.csv", index=False)

    en_braak = evidence[
        evidence["cell_class"].eq("EN")
        & evidence["phenotype"].eq("Braak")
    ].iloc[0]
    positive = en_braak["S"] > 0
    robust_locus = (
        not en_braak["leave_one_locus_crosses_zero"]
        and en_braak["maximum_absolute_locus_contribution_share"] <= 0.20
    )
    raw_p = float(en_braak["empirical_p_two_sided"])
    adjusted_p = float(en_braak["empirical_FDR_BH_six_external_tests"])
    if positive and raw_p <= 0.05 and adjusted_p <= 0.05 and robust_locus:
        decision = "EXPLORATORY_MULTI_TEST_ADJUSTED_EN_SUPPORT"
        interpretation = (
            "SEA-AD EN provides post-selection external cell-class support, "
            "but the RADC quality-weight global result remains only marginal."
        )
    elif positive and raw_p <= 0.05 and robust_locus:
        decision = "NOMINAL_EN_SUPPORT_NOT_MULTI_TEST_ADJUSTED"
        interpretation = (
            "SEA-AD EN is nominally positive but does not survive correction "
            "across the six external class-by-pathology tests."
        )
    elif positive and raw_p < 0.10 and robust_locus:
        decision = "WEAK_EN_TREND_ONLY"
        interpretation = (
            "SEA-AD EN shows only a weak positive trend and cannot be called "
            "external replication."
        )
    elif not positive:
        decision = "EN_DIRECTION_DISCORDANT"
        interpretation = (
            "SEA-AD EN is directionally discordant with the positive RADC EN "
            "candidate and does not externally replicate it."
        )
    else:
        decision = "NO_EN_EXTERNAL_REPLICATION"
        interpretation = (
            "SEA-AD EN does not externally replicate the RADC EN candidate."
        )

    summary = {
        "status": "COMPLETE",
        "decision": decision,
        "interpretation": interpretation,
        "decision_rules": {
            "strong_exploratory_support": (
                "EN Braak S > 0, raw P <= 0.05, BH across six external tests "
                "<= 0.05, leave-one-locus-out does not cross zero, and maximum "
                "absolute locus contribution <= 0.20"
            ),
            "scope_boundary": (
                "Oligo and EN were prioritized after inspecting RADC strata; "
                "therefore even a positive SEA-AD result is post-selection "
                "external support rather than preregistered confirmation."
            ),
            "stop_rule": (
                "After EN, do not screen further SEA-AD cell classes solely "
                "to search for significance."
            ),
        },
        "en_braak": {
            key: (
                bool(en_braak[key])
                if key == "leave_one_locus_crosses_zero"
                else float(en_braak[key])
                if key
                in {
                    "S",
                    "matched_null_mean",
                    "difference_from_null_mean",
                    "empirical_p_two_sided",
                    "empirical_FDR_BH_six_external_tests",
                    "leave_one_locus_min",
                    "leave_one_locus_max",
                    "maximum_absolute_locus_contribution_share",
                }
                else int(en_braak[key])
            )
            for key in (
                "loci",
                "S",
                "matched_null_mean",
                "difference_from_null_mean",
                "empirical_p_two_sided",
                "empirical_FDR_BH_six_external_tests",
                "leave_one_locus_min",
                "leave_one_locus_max",
                "leave_one_locus_crosses_zero",
                "maximum_absolute_locus_contribution_share",
            )
        },
        "generated_at": datetime.now(timezone.utc).astimezone().isoformat(),
    }
    (output_dir / "decision.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )

    table = evidence[
        [
            "cell_class",
            "phenotype",
            "loci",
            "S",
            "matched_null_mean",
            "empirical_p_two_sided",
            "empirical_FDR_BH_six_external_tests",
            "leave_one_locus_crosses_zero",
            "maximum_absolute_locus_contribution_share",
        ]
    ].copy()
    for column in (
        "S",
        "matched_null_mean",
        "empirical_p_two_sided",
        "empirical_FDR_BH_six_external_tests",
        "maximum_absolute_locus_contribution_share",
    ):
        table[column] = table[column].map(lambda value: f"{value:.4f}")
    lines = [
        "# SEA-AD cross-class external validation decision",
        "",
        f"- Decision: `{decision}`",
        f"- Interpretation: {interpretation}",
        "",
        markdown_table(table),
        "",
        "The six class-by-pathology empirical P values are BH-adjusted together.",
        "Oligo and EN were selected after the RADC stratum results and therefore",
        "remain post-selection external analyses. No further SEA-AD cell-class",
        "screening is authorized solely to seek statistical significance.",
        "",
    ]
    (output_dir / "DECISION.md").write_text(
        "\n".join(lines), encoding="utf-8"
    )
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
