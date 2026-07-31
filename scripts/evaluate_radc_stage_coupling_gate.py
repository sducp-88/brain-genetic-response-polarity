#!/usr/bin/env python3
"""Evaluate the frozen RADC stage-coupling escalation gate."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--standard-root", required=True, type=Path)
    parser.add_argument("--quality-root", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    return parser.parse_args()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def one_row(rows: list[dict[str, str]], **filters: str) -> dict[str, str]:
    selected = [
        row for row in rows
        if all(row.get(key) == value for key, value in filters.items())
    ]
    if len(selected) != 1:
        raise ValueError(
            f"Expected one row for {filters}, observed {len(selected)}"
        )
    return selected[0]


def as_float(value: str) -> float:
    return float(value)


def as_bool(value: str) -> bool:
    normalized = value.strip().lower()
    if normalized == "true":
        return True
    if normalized == "false":
        return False
    raise ValueError(f"Not a Boolean value: {value!r}")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_phenotype(root: Path, phenotype: str) -> dict[str, object]:
    coupling = root / phenotype / "coupling"
    group_path = coupling / "direction_coupling_group_results.csv"
    omnibus_path = coupling / "direction_coupling_omnibus.csv"
    subset_path = coupling / "direction_coupling_subset_sensitivity.csv"
    rho_path = coupling / "rho_grid_group_summary.csv"

    group = one_row(
        read_csv(group_path), level="global", disease="ALL", cell_class="ALL"
    )
    omnibus = one_row(read_csv(omnibus_path), scope="global", disease="ALL")
    subsets = read_csv(subset_path)
    g1 = one_row(
        subsets,
        sensitivity_set="G1_only",
        level="global",
        disease="ALL",
        cell_class="ALL",
    )
    exact = one_row(
        subsets,
        sensitivity_set="exact_major_class_only",
        level="global",
        disease="ALL",
        cell_class="ALL",
    )
    rho = one_row(
        read_csv(rho_path), level="global", disease="ALL", cell_class="ALL"
    )

    paths = [group_path, omnibus_path, subset_path, rho_path]
    return {
        "S": as_float(group["S_locus_equal_moderated_all_anchors"]),
        "matched_null_mean": as_float(group["matched_null_mean"]),
        "empirical_p": as_float(group["empirical_p_two_sided"]),
        "loci": int(group["loci"]),
        "anchor_rows": int(group["anchor_rows"]),
        "loo_min": as_float(group["leave_one_locus_min"]),
        "loo_max": as_float(group["leave_one_locus_max"]),
        "loo_crosses_zero": as_bool(group["leave_one_locus_crosses_zero"]),
        "max_locus_share": as_float(
            group["maximum_absolute_locus_contribution_share"]
        ),
        "omnibus_p": as_float(omnibus["empirical_p"]),
        "g1_S": as_float(g1["S_locus_equal_moderated"]),
        "exact_S": as_float(exact["S_locus_equal_moderated"]),
        "rho_min": as_float(rho["primary_grid_min"]),
        "rho_max": as_float(rho["primary_grid_max"]),
        "rho_crosses_zero": as_bool(rho["primary_grid_crosses_zero"]),
        "source_sha256": {str(path): sha256(path) for path in paths},
    }


def same_positive_sign(*values: float) -> bool:
    return all(value > 0 for value in values)


def evaluate(
    standard: dict[str, dict[str, object]],
    quality: dict[str, dict[str, object]],
) -> tuple[str, list[str], dict[str, bool]]:
    braak_s = quality["Braak"]
    checks = {
        "quality_braak_positive": braak_s["S"] > 0,
        "quality_braak_p_le_0_05": braak_s["empirical_p"] <= 0.05,
        "quality_braak_p_lt_0_10": braak_s["empirical_p"] < 0.10,
        "quality_loo_positive": (
            not braak_s["loo_crosses_zero"] and braak_s["loo_min"] > 0
        ),
        "quality_max_locus_share_le_0_20": braak_s["max_locus_share"] <= 0.20,
        "quality_rho_positive": (
            not braak_s["rho_crosses_zero"] and braak_s["rho_min"] > 0
        ),
        "quality_g1_and_exact_positive": same_positive_sign(
            braak_s["g1_S"], braak_s["exact_S"]
        ),
        "standard_quality_same_sign": (
            standard["Braak"]["S"] * braak_s["S"] > 0
        ),
    }
    fatal = (
        not checks["quality_braak_positive"]
        or braak_s["empirical_p"] >= 0.10
        or not checks["quality_max_locus_share_le_0_20"]
    )
    strong_noncomplex = all(
        checks[key]
        for key in [
            "quality_braak_positive",
            "quality_braak_p_le_0_05",
            "quality_loo_positive",
            "quality_max_locus_share_le_0_20",
            "quality_rho_positive",
            "quality_g1_and_exact_positive",
            "standard_quality_same_sign",
        ]
    )
    reasons: list[str] = []
    if fatal:
        decision = "NO_GO"
        reasons.append(
            "Quality-weight Braak failed a frozen no-go criterion."
        )
    elif strong_noncomplex:
        decision = "CONDITIONAL_GO_PENDING_COMPLEX_REGION_AUDIT"
        reasons.append(
            "All computable strong-pass criteria were met; APOE/complex-LD "
            "exclusion remains required before full escalation."
        )
    else:
        decision = "CONDITIONAL_PASS_EXPLORATORY_SEAAD_ONLY"
        reasons.append(
            "The positive direction was retained below the no-go boundary, "
            "but at least one strong-pass robustness criterion was unmet."
        )
    return decision, reasons, checks


def fmt(value: object, digits: int = 4) -> str:
    if isinstance(value, float):
        return f"{value:.{digits}f}"
    return str(value)


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    standard = {
        phenotype: load_phenotype(args.standard_root, phenotype)
        for phenotype in ["CERAD", "Braak"]
    }
    quality = {
        phenotype: load_phenotype(args.quality_root, phenotype)
        for phenotype in ["CERAD", "Braak"]
    }
    decision, reasons, checks = evaluate(standard, quality)
    payload = {
        "status": "COMPLETE",
        "decision": decision,
        "reasons": reasons,
        "checks": checks,
        "standard_voom_primary": standard,
        "quality_weight_sensitivity": quality,
        "complex_region_audit": "PENDING",
        "generated_at": datetime.now(timezone.utc).astimezone().isoformat(),
    }
    json_path = args.output_dir / "radc_stage_coupling_gate.json"
    json_path.write_text(
        json.dumps(payload, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )

    rows = [
        "# RADC stage-coupling gate result",
        "",
        f"Decision: **{decision}**",
        "",
        "The standard-voom analysis remains primary; quality weights are a "
        "sensitivity analysis.",
        "",
        "| Model | Axis | S | Empirical P | Omnibus P | LOO range | Max locus share |",
        "|---|---|---:|---:|---:|---:|---:|",
    ]
    for label, result in [
        ("Standard voom", standard),
        ("Quality weights", quality),
    ]:
        for phenotype in ["CERAD", "Braak"]:
            item = result[phenotype]
            rows.append(
                f"| {label} | {phenotype} | {fmt(item['S'])} | "
                f"{fmt(item['empirical_p'])} | {fmt(item['omnibus_p'])} | "
                f"{fmt(item['loo_min'])} to {fmt(item['loo_max'])} | "
                f"{fmt(item['max_locus_share'])} |"
            )
    rows += [
        "",
        "## Frozen checks",
        "",
    ]
    rows.extend(
        f"- {'PASS' if passed else 'FAIL'}: `{name}`"
        for name, passed in checks.items()
    )
    rows += [
        "",
        "## Remaining limitation",
        "",
        "APOE and other complex-LD exclusion has not yet been computed. "
        "The decision cannot be upgraded beyond conditional go until that "
        "audit is complete.",
        "",
    ]
    md_path = args.output_dir / "RADC_STAGE_COUPLING_GATE_RESULT.md"
    md_path.write_text("\n".join(rows), encoding="utf-8")
    print(json.dumps({"decision": decision, "output": str(json_path)}))


if __name__ == "__main__":
    main()
