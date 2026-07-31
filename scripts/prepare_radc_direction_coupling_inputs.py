#!/usr/bin/env python3
"""Prepare RADC continuous-pathology results for frozen direction coupling.

The adapter changes column names and file layout only. It does not refit models,
select anchors, or alter effect estimates.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import time
from datetime import datetime
from pathlib import Path

import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--model-dir", type=Path, required=True)
    parser.add_argument("--anchors", type=Path, required=True)
    parser.add_argument("--phenotype", choices=["CERAD", "Braak"], required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def clean_gene_id(values: pd.Series) -> pd.Series:
    return values.astype(str).str.replace(r"\..*$", "", regex=True)


def atomic_json(payload: dict, path: Path) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(payload, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    for attempt in range(30):
        try:
            os.replace(temporary, path)
            return
        except PermissionError:
            if attempt == 29:
                raise
            time.sleep(0.2)


def main() -> None:
    args = parse_args()
    model_dir = args.model_dir.resolve()
    anchor_path = args.anchors.resolve()
    output_dir = args.output_dir.resolve()
    standardized_dir = output_dir / "all_gene_results"
    standardized_dir.mkdir(parents=True, exist_ok=True)

    anchors = pd.read_csv(anchor_path, low_memory=False)
    required_anchor_columns = {
        "anchor_unit_id",
        "disease",
        "locus_id",
        "gene_id",
        "cell_class",
        "beta_G",
        "SE_G",
        "evidence_grade",
        "primary_anchor_eligible",
        "qtl_resolution",
    }
    missing = sorted(required_anchor_columns - set(anchors))
    if missing:
        raise RuntimeError(f"Anchor columns absent: {missing}")
    anchors = anchors.loc[
        anchors["disease"].eq("AD")
        & anchors["evidence_grade"].isin(["G1", "G2"])
        & anchors["primary_anchor_eligible"].eq("yes")
    ].copy()
    anchors["gene_id_clean"] = clean_gene_id(anchors["gene_id"])
    if anchors["anchor_unit_id"].duplicated().any():
        raise RuntimeError("Duplicated AD anchor_unit_id.")

    standardized_hashes: dict[str, str] = {}
    input_hashes: dict[str, str] = {}
    merged_rows: list[pd.DataFrame] = []
    class_diagnostics: list[dict] = []

    for cell_class in sorted(anchors["cell_class"].unique()):
        input_path = model_dir / f"{cell_class}__{args.phenotype}__gene_results.csv"
        if not input_path.is_file():
            raise FileNotFoundError(input_path)
        table = pd.read_csv(input_path, low_memory=False)
        required_model_columns = {
            "gene_id",
            "beta_D",
            "SE_D_moderated",
            "AveExpr",
            "class",
            "phenotype",
        }
        missing = sorted(required_model_columns - set(table))
        if missing:
            raise RuntimeError(f"{input_path.name} columns absent: {missing}")
        if not table["class"].eq(cell_class).all():
            raise RuntimeError(f"Class label mismatch in {input_path.name}.")
        if not table["phenotype"].eq(args.phenotype).all():
            raise RuntimeError(f"Phenotype label mismatch in {input_path.name}.")
        table["gene_id_clean"] = clean_gene_id(table["gene_id"])
        if table["gene_id_clean"].duplicated().any():
            raise RuntimeError(f"Duplicated gene identifiers in {input_path.name}.")
        for column in ["beta_D", "SE_D_moderated", "AveExpr"]:
            table[column] = pd.to_numeric(table[column], errors="coerce")
        finite = (
            table["beta_D"].map(math.isfinite)
            & table["SE_D_moderated"].map(math.isfinite)
            & table["SE_D_moderated"].gt(0)
            & table["AveExpr"].map(math.isfinite)
        )
        if not finite.all():
            raise RuntimeError(
                f"Non-finite model statistics in {input_path.name}: "
                f"{int((~finite).sum())} rows."
            )

        standardized = table[
            ["gene_id_clean", "beta_D", "SE_D_moderated", "AveExpr"]
        ].rename(
            columns={
                "SE_D_moderated": "SE_D",
                "AveExpr": "average_log2_expression",
            }
        )
        standardized_path = (
            standardized_dir / f"AD__{cell_class}_all_gene_results.csv.gz"
        )
        standardized.to_csv(
            standardized_path,
            index=False,
            compression={"method": "gzip", "compresslevel": 6, "mtime": 0},
        )

        class_anchors = anchors.loc[anchors["cell_class"].eq(cell_class)].copy()
        class_anchors = class_anchors.merge(
            standardized[["gene_id_clean", "beta_D", "SE_D"]],
            on="gene_id_clean",
            how="left",
            validate="many_to_one",
        )
        class_anchors["expression_estimable"] = (
            class_anchors["beta_D"].notna() & class_anchors["SE_D"].gt(0)
        )
        class_anchors["pathology_phenotype"] = args.phenotype
        class_anchors["disease_expression_resource"] = "PsychAD RADC"
        class_anchors["replication_role"] = (
            "stage_mechanism_not_independent_genetic_replication"
        )
        merged_rows.append(class_anchors)

        input_hashes[input_path.name] = sha256(input_path)
        standardized_hashes[standardized_path.name] = sha256(standardized_path)
        class_diagnostics.append(
            {
                "cell_class": cell_class,
                "all_gene_rows": int(len(standardized)),
                "anchor_rows": int(len(class_anchors)),
                "estimable_anchor_rows": int(
                    class_anchors["expression_estimable"].sum()
                ),
            }
        )

    anchor_results = pd.concat(merged_rows, ignore_index=True)
    if anchor_results["anchor_unit_id"].duplicated().any():
        raise RuntimeError("Duplicated anchor_unit_id after RADC merge.")
    anchor_output = output_dir / (
        f"radc_{args.phenotype}_G1_G2_anchor_pathology_effects.csv"
    )
    anchor_results.to_csv(anchor_output, index=False)

    manifest = {
        "status": "COMPLETE",
        "phenotype": args.phenotype,
        "effect_direction": (
            "positive beta_D means higher expression with a one-SD increase "
            f"in numeric {args.phenotype}"
        ),
        "adapter_boundary": (
            "Column and file-layout adaptation only; no model refitting, "
            "anchor reselection, or effect modification."
        ),
        "replication_boundary": (
            "RADC/ROSMAP overlaps the genetic-anchor ecosystem and is treated "
            "as stage-mechanism analysis, not independent genetic replication."
        ),
        "model_dir": str(model_dir),
        "anchor_source": str(anchor_path),
        "anchor_source_sha256": sha256(anchor_path),
        "anchor_rows": int(len(anchor_results)),
        "estimable_anchor_rows": int(anchor_results["expression_estimable"].sum()),
        "loci": int(anchor_results["locus_id"].nunique()),
        "class_diagnostics": class_diagnostics,
        "input_sha256": input_hashes,
        "output_sha256": {
            **standardized_hashes,
            anchor_output.name: sha256(anchor_output),
        },
        "completed_at": datetime.now().astimezone().isoformat(),
    }
    atomic_json(manifest, output_dir / "radc_coupling_input_manifest.json")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
