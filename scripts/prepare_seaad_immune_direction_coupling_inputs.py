#!/usr/bin/env python3
"""Adapt SEA-AD Immune pathology models to the frozen coupling engine."""

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
    parser.add_argument("--model-dir", required=True, type=Path)
    parser.add_argument("--anchors", required=True, type=Path)
    parser.add_argument("--phenotype", required=True, choices=["CERAD", "Braak"])
    parser.add_argument(
        "--cell-class",
        default="Immune",
        choices=["Astro", "EN", "Endo", "IN", "Immune", "Mural", "OPC", "Oligo"],
    )
    parser.add_argument("--output-dir", required=True, type=Path)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
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
    output_dir = args.output_dir.resolve()
    standardized_dir = output_dir / "all_gene_results"
    standardized_dir.mkdir(parents=True, exist_ok=True)

    anchors = pd.read_csv(args.anchors, low_memory=False)
    required_anchor = {
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
    missing = sorted(required_anchor - set(anchors))
    if missing:
        raise RuntimeError(f"Anchor columns absent: {missing}")
    anchors = anchors.loc[
        anchors["disease"].eq("AD")
        & anchors["cell_class"].eq(args.cell_class)
        & anchors["evidence_grade"].isin(["G1", "G2"])
        & anchors["primary_anchor_eligible"].eq("yes")
    ].copy()
    if anchors.empty:
        raise RuntimeError(f"No eligible AD {args.cell_class} anchors")
    anchors["gene_id_clean"] = clean_gene_id(anchors["gene_id"])
    if anchors["anchor_unit_id"].duplicated().any():
        raise RuntimeError(f"Duplicated AD {args.cell_class} anchor_unit_id")

    input_path = (
        model_dir / f"{args.cell_class}__{args.phenotype}__gene_results.csv"
    )
    table = pd.read_csv(input_path, low_memory=False)
    required_model = {
        "gene_id", "beta_D", "SE_D_moderated", "AveExpr", "class", "phenotype"
    }
    missing = sorted(required_model - set(table))
    if missing:
        raise RuntimeError(f"Model columns absent: {missing}")
    if not table["class"].eq(args.cell_class).all():
        raise RuntimeError(
            f"Model class is not uniformly {args.cell_class}"
        )
    if not table["phenotype"].eq(args.phenotype).all():
        raise RuntimeError("Model phenotype mismatch")
    table["gene_id_clean"] = clean_gene_id(table["gene_id"])
    if table["gene_id_clean"].duplicated().any():
        raise RuntimeError("Duplicated model gene identifiers")
    for column in ["beta_D", "SE_D_moderated", "AveExpr"]:
        table[column] = pd.to_numeric(table[column], errors="coerce")
    finite = (
        table["beta_D"].map(math.isfinite)
        & table["SE_D_moderated"].map(math.isfinite)
        & table["SE_D_moderated"].gt(0)
        & table["AveExpr"].map(math.isfinite)
    )
    if not finite.all():
        raise RuntimeError(f"Non-finite model rows: {int((~finite).sum())}")

    standardized = table[
        ["gene_id_clean", "beta_D", "SE_D_moderated", "AveExpr"]
    ].rename(
        columns={
            "SE_D_moderated": "SE_D",
            "AveExpr": "average_log2_expression",
        }
    )
    standardized_path = (
        standardized_dir / f"AD__{args.cell_class}_all_gene_results.csv.gz"
    )
    standardized.to_csv(
        standardized_path,
        index=False,
        compression={"method": "gzip", "compresslevel": 6, "mtime": 0},
    )
    merged = anchors.merge(
        standardized[["gene_id_clean", "beta_D", "SE_D"]],
        on="gene_id_clean",
        how="left",
        validate="many_to_one",
    )
    merged["expression_estimable"] = (
        merged["beta_D"].notna() & merged["SE_D"].gt(0)
    )
    merged["pathology_phenotype"] = args.phenotype
    merged["disease_expression_resource"] = "SEA-AD core3 donor pseudobulk"
    merged["replication_role"] = (
        "external_cohort_stage_validation_not_independent_genetic_anchor"
    )
    anchor_path = output_dir / (
        f"seaad_{args.phenotype}_{args.cell_class}_"
        "G1_G2_anchor_pathology_effects.csv"
    )
    merged.to_csv(anchor_path, index=False)
    manifest = {
        "status": "COMPLETE",
        "phenotype": args.phenotype,
        "cell_class": args.cell_class,
        "effect_direction": (
            "positive beta_D means higher expression per one-SD increase "
            f"in numeric {args.phenotype}"
        ),
        "model_dir": str(model_dir),
        "anchor_source": str(args.anchors.resolve()),
        "anchor_source_sha256": sha256(args.anchors),
        "anchor_rows": len(merged),
        "estimable_anchor_rows": int(merged["expression_estimable"].sum()),
        "loci": int(merged["locus_id"].nunique()),
        "replication_boundary": (
            "SEA-AD is an external cohort/stage validation, but upstream "
            "genetic-anchor resource overlap prevents a claim of fully "
            "independent genetic replication."
        ),
        "input_sha256": {input_path.name: sha256(input_path)},
        "output_sha256": {
            standardized_path.name: sha256(standardized_path),
            anchor_path.name: sha256(anchor_path),
        },
        "completed_at": datetime.now().astimezone().isoformat(),
    }
    atomic_json(
        manifest,
        output_dir
        / f"seaad_{args.cell_class.lower()}_coupling_input_manifest.json",
    )
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
