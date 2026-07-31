#!/usr/bin/env python3
"""Adapt frozen RADC stage contrasts to the direction-coupling interface."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
from datetime import datetime
from pathlib import Path

import pandas as pd


CONTRAST_COLUMNS = {
    "middle_vs_low": ("middle_vs_low_beta", "middle_vs_low_SE"),
    "high_vs_low": ("high_vs_low_beta", "high_vs_low_SE"),
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--nonlinear-dir", type=Path, required=True)
    parser.add_argument("--primary-dir", type=Path, required=True)
    parser.add_argument("--anchors", type=Path, required=True)
    parser.add_argument("--phenotype", choices=["CERAD", "Braak"], required=True)
    parser.add_argument(
        "--contrast", choices=sorted(CONTRAST_COLUMNS), required=True
    )
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def clean_gene_id(values: pd.Series) -> pd.Series:
    return values.astype(str).str.replace(r"\..*$", "", regex=True)


def atomic_json(value: dict, path: Path) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(value, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    os.replace(temporary, path)


def main() -> None:
    args = parse_args()
    nonlinear_dir = args.nonlinear_dir.resolve()
    primary_dir = args.primary_dir.resolve()
    anchors_path = args.anchors.resolve()
    output_dir = args.output_dir.resolve()
    all_gene_dir = output_dir / "all_gene_results"
    all_gene_dir.mkdir(parents=True, exist_ok=True)
    beta_column, se_column = CONTRAST_COLUMNS[args.contrast]

    anchors = pd.read_csv(anchors_path, low_memory=False)
    anchors = anchors.loc[
        anchors["disease"].eq("AD")
        & anchors["evidence_grade"].isin(["G1", "G2"])
        & anchors["primary_anchor_eligible"].eq("yes")
    ].copy()
    anchors["gene_id_clean"] = clean_gene_id(anchors["gene_id"])
    if anchors["anchor_unit_id"].duplicated().any():
        raise RuntimeError("Duplicated anchor_unit_id")

    merged_rows: list[pd.DataFrame] = []
    inputs: dict[str, str] = {}
    outputs: dict[str, str] = {}
    diagnostics: list[dict[str, object]] = []

    for cell_class in sorted(anchors["cell_class"].unique()):
        nonlinear_path = (
            nonlinear_dir
            / f"{cell_class}__{args.phenotype}__nonlinear_results.csv"
        )
        primary_path = (
            primary_dir / f"{cell_class}__{args.phenotype}__gene_results.csv"
        )
        nonlinear = pd.read_csv(nonlinear_path, low_memory=False)
        primary = pd.read_csv(primary_path, low_memory=False)
        required_nonlinear = {"gene_id", beta_column, se_column}
        required_primary = {"gene_id", "AveExpr"}
        if not required_nonlinear <= set(nonlinear):
            raise RuntimeError(f"Missing columns in {nonlinear_path}")
        if not required_primary <= set(primary):
            raise RuntimeError(f"Missing columns in {primary_path}")
        nonlinear["gene_id_clean"] = clean_gene_id(nonlinear["gene_id"])
        primary["gene_id_clean"] = clean_gene_id(primary["gene_id"])
        if nonlinear["gene_id_clean"].duplicated().any():
            raise RuntimeError(f"Duplicated genes in {nonlinear_path}")
        if primary["gene_id_clean"].duplicated().any():
            raise RuntimeError(f"Duplicated genes in {primary_path}")

        table = nonlinear[
            ["gene_id_clean", beta_column, se_column]
        ].merge(
            primary[["gene_id_clean", "AveExpr"]],
            on="gene_id_clean",
            how="inner",
            validate="one_to_one",
        )
        table = table.rename(
            columns={
                beta_column: "beta_D",
                se_column: "SE_D",
                "AveExpr": "average_log2_expression",
            }
        )
        for column in ["beta_D", "SE_D", "average_log2_expression"]:
            table[column] = pd.to_numeric(table[column], errors="coerce")
        finite = (
            table["beta_D"].map(math.isfinite)
            & table["SE_D"].map(math.isfinite)
            & table["SE_D"].gt(0)
            & table["average_log2_expression"].map(math.isfinite)
        )
        table = table.loc[finite].copy()
        standardized_path = (
            all_gene_dir / f"AD__{cell_class}_all_gene_results.csv.gz"
        )
        table.to_csv(
            standardized_path,
            index=False,
            compression={"method": "gzip", "compresslevel": 6, "mtime": 0},
        )

        class_anchors = anchors.loc[anchors["cell_class"].eq(cell_class)].copy()
        class_anchors = class_anchors.merge(
            table[["gene_id_clean", "beta_D", "SE_D"]],
            on="gene_id_clean",
            how="left",
            validate="many_to_one",
        )
        class_anchors["expression_estimable"] = (
            class_anchors["beta_D"].notna() & class_anchors["SE_D"].gt(0)
        )
        class_anchors["pathology_phenotype"] = args.phenotype
        class_anchors["pathology_contrast"] = args.contrast
        class_anchors["disease_expression_resource"] = "PsychAD RADC"
        class_anchors["replication_role"] = (
            "stage_mechanism_not_independent_genetic_replication"
        )
        merged_rows.append(class_anchors)
        inputs[str(nonlinear_path)] = sha256(nonlinear_path)
        inputs[str(primary_path)] = sha256(primary_path)
        outputs[str(standardized_path)] = sha256(standardized_path)
        diagnostics.append(
            {
                "cell_class": cell_class,
                "background_genes": int(len(table)),
                "anchor_rows": int(len(class_anchors)),
                "estimable_anchor_rows": int(
                    class_anchors["expression_estimable"].sum()
                ),
            }
        )

    anchor_results = pd.concat(merged_rows, ignore_index=True)
    anchor_path = output_dir / (
        f"radc_{args.phenotype}_{args.contrast}_G1_G2_anchor_effects.csv"
    )
    anchor_results.to_csv(anchor_path, index=False)
    outputs[str(anchor_path)] = sha256(anchor_path)
    manifest = {
        "status": "COMPLETE",
        "phenotype": args.phenotype,
        "contrast": args.contrast,
        "contrast_direction": (
            f"positive beta_D means higher expression for {args.contrast}"
        ),
        "primary_boundary": (
            "Stage contrast is a frozen functional-form sensitivity, not a "
            "replacement for the continuous primary model."
        ),
        "anchor_source": str(anchors_path),
        "anchor_source_sha256": sha256(anchors_path),
        "anchor_rows": int(len(anchor_results)),
        "estimable_anchor_rows": int(
            anchor_results["expression_estimable"].sum()
        ),
        "diagnostics": diagnostics,
        "input_sha256": inputs,
        "output_sha256": outputs,
        "completed_at": datetime.now().astimezone().isoformat(),
    }
    atomic_json(manifest, output_dir / "stage_coupling_input_manifest.json")
    print(
        f"COMPLETE phenotype={args.phenotype} contrast={args.contrast} "
        f"estimable={manifest['estimable_anchor_rows']}"
    )


if __name__ == "__main__":
    main()

