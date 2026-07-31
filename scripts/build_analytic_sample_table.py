#!/usr/bin/env python3
"""Build a non-identifying analytic-sample table from frozen aggregate outputs."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from datetime import datetime
from pathlib import Path


SEA_CLASSES = ("Immune", "Oligo", "EN")
RADC_CLASSES = ("Astro", "EN", "IN", "Immune", "OPC", "Oligo")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    return parser.parse_args()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def read_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def int_range(values: list[str]) -> str:
    numbers = sorted({int(float(value)) for value in values})
    if len(numbers) == 1:
        return f"{numbers[0]:,}"
    return f"{numbers[0]:,}-{numbers[-1]:,}"


def float_range(values: list[str]) -> str:
    numbers = sorted(float(value) for value in values)
    if abs(numbers[0] - numbers[-1]) < 0.005:
        return f"{numbers[0]:.1f}"
    return f"{numbers[0]:.1f}-{numbers[-1]:.1f}"


def disease_anchor_counts(
    rows: list[dict[str, str]], disease: str
) -> tuple[int, int]:
    match = next(
        row
        for row in rows
        if row["level"] == "disease" and row["disease"] == disease
    )
    return int(match["anchor_rows"]), int(match["loci"])


def sea_anchor_counts(
    rows: list[dict[str, str]], cell_class: str
) -> tuple[int, int]:
    matches = [
        row
        for row in rows
        if row["cell_class"] == cell_class
        and row["model_label"] == "DFC_standard"
        and row["primary_DFC_standard"] == "True"
    ]
    pairs = {(int(row["anchor_rows"]), int(row["loci"])) for row in matches}
    if len(matches) != 2 or len(pairs) != 1:
        raise ValueError(
            f"Expected two concordant primary DFC rows for {cell_class}; "
            f"found {len(matches)} rows and {len(pairs)} count pairs."
        )
    return pairs.pop()


def markdown_table(rows: list[dict[str, str]]) -> str:
    columns = [
        ("Resource", "resource"),
        ("Role/outcome", "role_outcome"),
        ("Region(s)", "regions"),
        ("Cell-class scope", "cell_scope"),
        ("Eligible processed data", "eligible_data"),
        ("Primary model N", "primary_n"),
        ("Anchors/loci", "anchors_loci"),
    ]
    output = [
        "# Table 1. Cohort roles and analytic sample coverage",
        "",
        "| " + " | ".join(label for label, _ in columns) + " |",
        "|" + "|".join("---" for _ in columns) + "|",
    ]
    for row in rows:
        values = [row[key].replace("|", r"\|") for _, key in columns]
        output.append("| " + " | ".join(values) + " |")
    output.extend(
        [
            "",
            (
                "**Notes.** Counts are non-identifying aggregate values read from "
                "frozen manifests. MSSM case/control and RADC donor counts are ranges "
                "across fitted broad cell classes, not pooled sample sizes. SEA-AD "
                "eligible processed counts precede complete-case covariate filtering; "
                "the repeated oligodendrocyte model therefore used 240 donor-region "
                "samples from 83 donors rather than all 241 samples from 84 donors. "
                "DLPFC, dorsolateral prefrontal cortex; DFC, MEC, and MTG are the "
                "frozen SEA-AD region labels; ESS, effective sample size."
            ),
            "",
            (
                "**Interpretive boundary.** RADC/ROSMAP and SEA-AD overlap the "
                "genetic-anchor ecosystem and are expression-context analyses, not "
                "fully independent genetic replications."
            ),
            "",
        ]
    )
    return "\n".join(output)


def main() -> None:
    args = parse_args()
    root = args.project_root.resolve()
    outdir = args.output_dir.resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    inputs = {
        "mssm_readiness": root
        / "outputs/phase2/mssm_verified_disease_cell_readiness/"
        "disease_cell_readiness.csv",
        "mssm_coupling": root
        / "outputs/phase4/mssm_direction_coupling/"
        "direction_coupling_group_results.csv",
        "radc_model_manifest": root
        / "outputs/phase2/radc_pathology_models/min_cells_20/model_manifest.csv",
        "radc_anchor_manifest": root
        / "outputs/phase4/radc_stage_coupling_standard_voom/CERAD/input/"
        "radc_coupling_input_manifest.json",
        "sea_fixed_results": root
        / "outputs/phase6/seaad_region_sensitivity_v2/summary/"
        "seaad_fixed_42_results.csv",
    }
    for cell_class in SEA_CLASSES:
        inputs[f"sea_{cell_class}_metadata"] = root / (
            "outputs/phase6/seaad_region_sensitivity_v1/"
            f"{cell_class}/metadata_gate.json"
        )
        inputs[f"sea_{cell_class}_dfc_manifest"] = root / (
            "outputs/phase6/seaad_region_sensitivity_v2/"
            f"{cell_class}/models_DFC_standard/model_manifest.csv"
        )
        inputs[f"sea_{cell_class}_repeated_manifest"] = root / (
            "outputs/phase6/seaad_region_sensitivity_v2/"
            f"{cell_class}/models_repeated_standard/model_manifest.csv"
        )

    missing = [str(path) for path in inputs.values() if not path.is_file()]
    if missing:
        raise FileNotFoundError("Missing frozen input(s):\n" + "\n".join(missing))

    mssm = read_csv(inputs["mssm_readiness"])
    mssm_coupling = read_csv(inputs["mssm_coupling"])
    rows: list[dict[str, str]] = []

    for disease in ("AD", "SCZ"):
        fitted = [
            row
            for row in mssm
            if row["disease"] == disease and row["confirmatory_ready"] == "True"
        ]
        anchor_rows, loci = disease_anchor_counts(mssm_coupling, disease)
        rows.append(
            {
                "resource": "MSSM/PsychAD",
                "role_outcome": f"Disease-state test; {disease} vs control",
                "regions": "DLPFC",
                "cell_scope": f"{len(fitted)} broad classes",
                "eligible_data": (
                    f"Cases {int_range([r['cases'] for r in fitted])}; "
                    f"controls {int_range([r['controls'] for r in fitted])} "
                    "per class"
                ),
                "primary_n": (
                    f"Weighted ESS: cases "
                    f"{float_range([r['case_ess'] for r in fitted])}; "
                    f"controls {float_range([r['control_ess'] for r in fitted])}"
                ),
                "anchors_loci": f"{anchor_rows}/{loci}",
                "provenance_note": (
                    "Only confirmatory_ready=True disease-by-class strata included."
                ),
            }
        )

    radc = [
        row
        for row in read_csv(inputs["radc_model_manifest"])
        if row["phenotype"] in {"CERAD", "Braak"}
        and row["class"] in RADC_CLASSES
        and row["status"] == "complete"
    ]
    if len(radc) != 12:
        raise ValueError(f"Expected 12 RADC model rows; found {len(radc)}.")
    radc_anchor = read_json(inputs["radc_anchor_manifest"])
    if int(radc_anchor["estimable_anchor_rows"]) != 72:
        raise ValueError("RADC estimable anchor-row count drifted from 72.")
    rows.append(
        {
            "resource": "RADC/ROSMAP",
            "role_outcome": "Co-primary continuous CERAD and Braak pathology",
            "regions": "DLPFC",
            "cell_scope": "6 anchor-bearing broad classes",
            "eligible_data": (
                f"{int_range([r['n_donors'] for r in radc])} donors per class; "
                f"{int_range([r['genes_tested'] for r in radc])} genes tested"
            ),
            "primary_n": "Same donor range in each pathology model",
            "anchors_loci": "72/36",
            "provenance_note": (
                "Anchor manifest reports 73 source rows and 37 source loci; final "
                "coupling retains 72 estimable rows from 36 loci."
            ),
        }
    )

    sea_fixed = read_csv(inputs["sea_fixed_results"])
    for cell_class in SEA_CLASSES:
        metadata = read_json(inputs[f"sea_{cell_class}_metadata"])
        dfc = [
            row
            for row in read_csv(inputs[f"sea_{cell_class}_dfc_manifest"])
            if row["status"] == "complete"
            and row["phenotype"] in {"CERAD", "Braak"}
        ]
        repeated = [
            row
            for row in read_csv(inputs[f"sea_{cell_class}_repeated_manifest"])
            if row["status"] == "complete"
            and row["phenotype"] in {"CERAD", "Braak"}
        ]
        if len(dfc) != 2 or len(repeated) != 2:
            raise ValueError(f"Unexpected SEA-AD manifest rows for {cell_class}.")
        repeated_pairs = {
            (int(row["n_samples"]), int(row["n_donors"])) for row in repeated
        }
        if len(repeated_pairs) != 1:
            raise ValueError(f"SEA-AD repeated-model N differs by outcome: {cell_class}.")
        repeated_samples, repeated_donors = repeated_pairs.pop()
        anchor_rows, loci = sea_anchor_counts(sea_fixed, cell_class)
        donors_by_region = metadata["donors_by_region"]
        rows.append(
            {
                "resource": f"SEA-AD {cell_class}",
                "role_outcome": "Fixed DFC family; regional/composition sensitivities",
                "regions": "DFC, MEC, MTG",
                "cell_scope": cell_class,
                "eligible_data": (
                    f"{int(metadata['eligible_donor_region_samples'])} donor-region "
                    f"samples/{int(metadata['unique_donors'])} donors "
                    f"(DFC {int(donors_by_region['DFC'])}, "
                    f"MEC {int(donors_by_region['MEC'])}, "
                    f"MTG {int(donors_by_region['MTG'])})"
                ),
                "primary_n": (
                    f"DFC {int(dfc[0]['n_donors'])} donors; repeated model "
                    f"{repeated_samples} samples/{repeated_donors} donors; "
                    f"{int_range([r['genes_tested'] for r in dfc])} DFC genes"
                ),
                "anchors_loci": f"{anchor_rows}/{loci}",
                "provenance_note": (
                    "Eligible counts come from the metadata gate; model N comes "
                    "from complete-case standard manifests."
                ),
            }
        )

    fieldnames = [
        "resource",
        "role_outcome",
        "regions",
        "cell_scope",
        "eligible_data",
        "primary_n",
        "anchors_loci",
        "provenance_note",
    ]
    csv_path = outdir / "TABLE_1_ANALYTIC_SAMPLE_COVERAGE.csv"
    with csv_path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    md_path = outdir / "TABLE_1_ANALYTIC_SAMPLE_COVERAGE.md"
    md_path.write_text(markdown_table(rows), encoding="utf-8")

    manifest = {
        "status": "COMPLETE",
        "generated_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "purpose": (
            "Non-identifying cohort role and analytic-sample coverage table; "
            "no person-level data exported."
        ),
        "row_count": len(rows),
        "input_files": [
            {
                "key": key,
                "relative_path": path.relative_to(root).as_posix(),
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
            for key, path in inputs.items()
        ],
        "output_files": [
            {
                "relative_path": path.relative_to(outdir).as_posix(),
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
            for path in (csv_path, md_path)
        ],
        "validation": {
            "expected_rows": 6,
            "actual_rows": len(rows),
            "contains_person_level_identifiers": False,
            "mssm_confirmatory_ready_filter_applied": True,
            "radc_final_anchor_count_override_verified": "72 rows/36 loci",
            "sea_eligible_and_complete_case_counts_separated": True,
        },
    }
    if len(rows) != 6:
        raise ValueError(f"Expected 6 output rows; found {len(rows)}.")

    manifest_path = outdir / "TABLE_1_PROVENANCE.json"
    manifest_path.write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    print(f"Wrote {csv_path}")
    print(f"Wrote {md_path}")
    print(f"Wrote {manifest_path}")


if __name__ == "__main__":
    main()
