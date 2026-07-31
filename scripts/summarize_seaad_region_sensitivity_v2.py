#!/usr/bin/env python3
"""Summarize the frozen SEA-AD region-sensitivity v2 analysis.

This script is intentionally limited to the 42 pre-specified
class-by-model-by-pathology coupling results. It adds locus-bootstrap
uncertainty, fixed-family multiplicity control, and paired regional
contrasts without selecting outcomes.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd


BOOTSTRAPS = 10_000
SEED = 20_260_731
CLASSES = ("Immune", "Oligo", "EN")
PATHOLOGIES = ("CERAD", "Braak")
MODEL_LABELS = (
    "DFC_standard",
    "DFC_composition",
    "MEC_standard",
    "MEC_composition",
    "MTG_standard",
    "MTG_composition",
    "repeated_standard",
)


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


def adjust_bh(values: pd.Series) -> np.ndarray:
    p = pd.to_numeric(values, errors="coerce").to_numpy(dtype=float)
    result = np.full(len(p), np.nan)
    valid = np.flatnonzero(np.isfinite(p))
    if len(valid) == 0:
        return result
    order = np.argsort(p[valid])
    ranked = p[valid][order]
    adjusted = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    restored = np.empty_like(adjusted)
    restored[order] = np.minimum(adjusted, 1.0)
    result[valid] = restored
    return result


def bootstrap_mean(
    values: np.ndarray,
    rng: np.random.Generator,
) -> tuple[float, float, float, np.ndarray]:
    values = np.asarray(values, dtype=float)
    if values.ndim != 1 or len(values) == 0 or not np.isfinite(values).all():
        raise RuntimeError("Bootstrap values must be finite and non-empty")
    draws = rng.choice(
        values,
        size=(BOOTSTRAPS, len(values)),
        replace=True,
    )
    statistics = draws.mean(axis=1)
    lower, median, upper = np.quantile(
        statistics, [0.025, 0.5, 0.975]
    )
    return float(lower), float(median), float(upper), statistics


def model_metadata(model_label: str) -> tuple[str, str, str]:
    if model_label == "repeated_standard":
        return "repeated", "DFC_MEC_MTG", "standard"
    region, adjustment = model_label.split("_", maxsplit=1)
    return "region", region, adjustment


def coupling_root(
    root: Path,
    cell_class: str,
    model_label: str,
    pathology: str,
) -> Path:
    return (
        root
        / cell_class
        / "coupling"
        / model_label
        / pathology
        / "coupling"
    )


def read_anchor_probabilities(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path, low_memory=False)
    required = {
        "anchor_unit_id",
        "locus_id",
        "cell_class",
        "S_moderated_rho_0",
    }
    missing = sorted(required - set(frame))
    if missing:
        raise RuntimeError(f"Missing columns in {path}: {missing}")
    if frame["anchor_unit_id"].duplicated().any():
        raise RuntimeError(f"Duplicated anchor rows in {path}")
    frame["S_moderated_rho_0"] = pd.to_numeric(
        frame["S_moderated_rho_0"], errors="coerce"
    )
    return frame.dropna(
        subset=["locus_id", "S_moderated_rho_0"]
    ).copy()


def summarize_all(
    project: Path,
    rng: np.random.Generator,
    input_paths: set[Path],
) -> tuple[pd.DataFrame, dict[tuple[str, str, str], pd.DataFrame]]:
    root = project / "outputs/phase6/seaad_region_sensitivity_v2"
    records: list[dict[str, object]] = []
    anchors: dict[tuple[str, str, str], pd.DataFrame] = {}

    for cell_class in CLASSES:
        for model_label in MODEL_LABELS:
            analysis, region, adjustment = model_metadata(model_label)
            manifest_path = (
                root
                / cell_class
                / f"models_{model_label}"
                / "model_manifest.csv"
            )
            input_paths.add(manifest_path)
            manifest = pd.read_csv(manifest_path)
            if set(manifest["phenotype"]) != set(PATHOLOGIES):
                raise RuntimeError(
                    f"Incomplete model manifest: {manifest_path}"
                )
            for pathology in PATHOLOGIES:
                path_root = coupling_root(
                    root, cell_class, model_label, pathology
                )
                group_path = (
                    path_root / "direction_coupling_group_results.csv"
                )
                anchor_path = (
                    path_root / "anchor_direction_probabilities.csv.gz"
                )
                input_paths.update({group_path, anchor_path})
                group = pd.read_csv(group_path)
                row = group.loc[
                    group["level"].eq("stratum")
                    & group["cell_class"].eq(cell_class)
                ]
                if len(row) != 1:
                    raise RuntimeError(
                        f"Expected one stratum row in {group_path}"
                    )
                row = row.iloc[0]
                anchor = read_anchor_probabilities(anchor_path)
                anchor = anchor.loc[
                    anchor["cell_class"].eq(cell_class)
                ].copy()
                if anchor.empty:
                    raise RuntimeError(f"No anchors in {anchor_path}")
                by_locus = anchor.groupby(
                    "locus_id", observed=True
                )["S_moderated_rho_0"].mean()
                lower, median, upper, _ = bootstrap_mean(
                    by_locus.to_numpy(), rng
                )
                observed = float(by_locus.mean())
                reported = float(
                    row["S_locus_equal_moderated_all_anchors"]
                )
                if abs(observed - reported) > 1e-12:
                    raise RuntimeError(
                        f"Point-estimate mismatch in {group_path}: "
                        f"{observed} vs {reported}"
                    )
                model_row = manifest.loc[
                    manifest["phenotype"].eq(pathology)
                ].iloc[0]
                records.append(
                    {
                        "cell_class": cell_class,
                        "pathology": pathology,
                        "model_label": model_label,
                        "analysis": analysis,
                        "region": region,
                        "adjustment": adjustment,
                        "primary_DFC_standard": bool(
                            model_label == "DFC_standard"
                        ),
                        "anchor_rows": int(row["anchor_rows"]),
                        "loci": int(row["loci"]),
                        "n_samples": int(model_row["n_samples"]),
                        "n_donors": int(model_row["n_donors"]),
                        "S": observed,
                        "matched_null_mean": float(
                            row["matched_null_mean"]
                        ),
                        "empirical_p_raw": float(
                            row["empirical_p_two_sided"]
                        ),
                        "bootstrap_median": median,
                        "bootstrap_CI95_lower": lower,
                        "bootstrap_CI95_upper": upper,
                        "low_locus_count_exploratory": bool(
                            len(by_locus) < 10
                        ),
                        "bootstrap_replicates": BOOTSTRAPS,
                        "bootstrap_seed": SEED,
                    }
                )
                anchors[
                    (cell_class, model_label, pathology)
                ] = anchor

    result = pd.DataFrame.from_records(records)
    if len(result) != 42:
        raise RuntimeError(f"Expected 42 fixed results; observed {len(result)}")
    result["empirical_p_BH_within_fixed_model_six"] = np.nan
    for model_label in MODEL_LABELS:
        selected = result["model_label"].eq(model_label)
        if int(selected.sum()) != 6:
            raise RuntimeError(f"Incomplete family for {model_label}")
        result.loc[
            selected, "empirical_p_BH_within_fixed_model_six"
        ] = adjust_bh(result.loc[selected, "empirical_p_raw"])
    result["passes_fixed_model_BH05"] = result[
        "empirical_p_BH_within_fixed_model_six"
    ].lt(0.05)
    return result, anchors


def regional_contrasts(
    anchors: dict[tuple[str, str, str], pd.DataFrame],
    rng: np.random.Generator,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    records: list[dict[str, object]] = []
    locus_tables: list[pd.DataFrame] = []
    for adjustment in ("standard", "composition"):
        for cell_class in CLASSES:
            for pathology in PATHOLOGIES:
                dfc = anchors[
                    (cell_class, f"DFC_{adjustment}", pathology)
                ][
                    [
                        "anchor_unit_id",
                        "locus_id",
                        "S_moderated_rho_0",
                    ]
                ].rename(
                    columns={
                        "locus_id": "locus_id_DFC",
                        "S_moderated_rho_0": "S_DFC",
                    }
                )
                for other in ("MEC", "MTG"):
                    other_frame = anchors[
                        (
                            cell_class,
                            f"{other}_{adjustment}",
                            pathology,
                        )
                    ][
                        [
                            "anchor_unit_id",
                            "locus_id",
                            "S_moderated_rho_0",
                        ]
                    ].rename(
                        columns={
                            "locus_id": f"locus_id_{other}",
                            "S_moderated_rho_0": f"S_{other}",
                        }
                    )
                    paired = dfc.merge(
                        other_frame,
                        on="anchor_unit_id",
                        how="inner",
                        validate="one_to_one",
                    )
                    if not paired["locus_id_DFC"].eq(
                        paired[f"locus_id_{other}"]
                    ).all():
                        raise RuntimeError(
                            f"Locus mismatch: {cell_class}/{pathology}/"
                            f"{adjustment}/{other}"
                        )
                    paired["locus_id"] = paired["locus_id_DFC"]
                    by_locus = (
                        paired.groupby("locus_id", observed=True)
                        .agg(
                            anchors=("anchor_unit_id", "nunique"),
                            S_DFC=("S_DFC", "mean"),
                            S_other=(f"S_{other}", "mean"),
                        )
                        .reset_index()
                    )
                    by_locus["Delta_S_DFC_minus_other"] = (
                        by_locus["S_DFC"] - by_locus["S_other"]
                    )
                    values = by_locus[
                        "Delta_S_DFC_minus_other"
                    ].to_numpy()
                    lower, median, upper, bootstrap = bootstrap_mean(
                        values, rng
                    )
                    observed = float(values.mean())
                    centered = bootstrap - float(bootstrap.mean())
                    p_value = (
                        1
                        + float(
                            np.sum(np.abs(centered) >= abs(observed))
                        )
                    ) / (BOOTSTRAPS + 1)
                    records.append(
                        {
                            "adjustment": adjustment,
                            "cell_class": cell_class,
                            "pathology": pathology,
                            "contrast": f"DFC_minus_{other}",
                            "other_region": other,
                            "common_anchor_rows": int(len(paired)),
                            "common_loci": int(len(by_locus)),
                            "S_DFC_common": float(
                                by_locus["S_DFC"].mean()
                            ),
                            "S_other_common": float(
                                by_locus["S_other"].mean()
                            ),
                            "Delta_S_DFC_minus_other": observed,
                            "bootstrap_median_Delta": median,
                            "bootstrap_CI95_lower": lower,
                            "bootstrap_CI95_upper": upper,
                            "bootstrap_centered_two_sided_p": p_value,
                            "low_locus_count_exploratory": bool(
                                len(by_locus) < 10
                            ),
                            "bootstrap_replicates": BOOTSTRAPS,
                            "bootstrap_seed": SEED,
                        }
                    )
                    by_locus.insert(0, "other_region", other)
                    by_locus.insert(0, "pathology", pathology)
                    by_locus.insert(0, "cell_class", cell_class)
                    by_locus.insert(0, "adjustment", adjustment)
                    locus_tables.append(by_locus)
    summary = pd.DataFrame.from_records(records)
    if len(summary) != 24:
        raise RuntimeError(
            f"Expected 24 regional contrasts; observed {len(summary)}"
        )
    summary["bootstrap_p_BH_within_adjustment_12"] = np.nan
    for adjustment in ("standard", "composition"):
        selected = summary["adjustment"].eq(adjustment)
        summary.loc[
            selected, "bootstrap_p_BH_within_adjustment_12"
        ] = adjust_bh(
            summary.loc[
                selected, "bootstrap_centered_two_sided_p"
            ]
        )
    summary["passes_adjustment_BH05"] = summary[
        "bootstrap_p_BH_within_adjustment_12"
    ].lt(0.05)
    loci = pd.concat(locus_tables, ignore_index=True)
    return summary, loci


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

    summary, anchors = summarize_all(project, rng, input_paths)
    contrasts, contrast_loci = regional_contrasts(anchors, rng)

    summary_path = output_dir / "seaad_fixed_42_results.csv"
    contrast_path = output_dir / "seaad_regional_paired_contrasts.csv"
    contrast_loci_path = output_dir / "seaad_regional_locus_components.csv"
    summary.to_csv(summary_path, index=False)
    contrasts.to_csv(contrast_path, index=False)
    contrast_loci.to_csv(contrast_loci_path, index=False)

    primary = summary.loc[
        summary["model_label"].eq("DFC_standard"),
        [
            "cell_class",
            "pathology",
            "n_donors",
            "anchor_rows",
            "loci",
            "S",
            "bootstrap_CI95_lower",
            "bootstrap_CI95_upper",
            "empirical_p_raw",
            "empirical_p_BH_within_fixed_model_six",
            "passes_fixed_model_BH05",
            "low_locus_count_exploratory",
        ],
    ]
    report = [
        "# SEA-AD region-sensitivity v2 summary",
        "",
        f"Generated: {timestamp()}  ",
        f"Fixed results retained: {len(summary)}  ",
        f"Regional paired contrasts retained: {len(contrasts)}  ",
        f"Locus bootstrap: {BOOTSTRAPS:,}; seed {SEED}",
        "",
        "## Frozen DFC primary family",
        "",
        markdown_table(primary),
        "",
        "BH adjustment is across the six fixed DFC-standard "
        "class-by-pathology tests. Composition, MEC, MTG, and repeated "
        "models are mandatory sensitivities and cannot replace DFC.",
        "",
        "## Interpretation boundary",
        "",
        "- Confidence intervals containing zero do not establish equivalence.",
        "- Fewer than 10 contributing loci is low-locus exploratory evidence.",
        "- Region differences require the direct paired regional contrast; "
        "different significance labels alone are not evidence of interaction.",
        "- No protective, compensatory, causal, or replicated-mechanism "
        "claim is authorized by this table.",
    ]
    report_path = output_dir / "SEAAD_REGION_V2_REPORT.md"
    report_path.write_text("\n".join(report) + "\n", encoding="utf-8")

    output_paths = [
        summary_path,
        contrast_path,
        contrast_loci_path,
        report_path,
    ]
    manifest = {
        "created_at": timestamp(),
        "script": str(Path(__file__).resolve()),
        "script_sha256": sha256(Path(__file__).resolve()),
        "bootstrap_replicates": BOOTSTRAPS,
        "seed": SEED,
        "fixed_result_count": int(len(summary)),
        "fixed_regional_contrast_count": int(len(contrasts)),
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
    }
    (output_dir / "run_manifest.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    print(
        f"COMPLETE fixed_results={len(summary)} "
        f"regional_contrasts={len(contrasts)} "
        f"output={output_dir}"
    )


if __name__ == "__main__":
    main()
