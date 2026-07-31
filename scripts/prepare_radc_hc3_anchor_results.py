#!/usr/bin/env python3
"""Merge RADC anchor-level HC3 estimates with the frozen anchor resource."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from datetime import datetime
from pathlib import Path

import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--hc3-dir", type=Path, required=True)
    parser.add_argument("--anchors", type=Path, required=True)
    parser.add_argument("--phenotype", choices=["CERAD", "Braak"], required=True)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def clean(values: pd.Series) -> pd.Series:
    return values.astype(str).str.replace(r"\..*$", "", regex=True)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    args = parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    anchors = pd.read_csv(args.anchors, low_memory=False)
    anchors = anchors.loc[
        anchors["disease"].eq("AD")
        & anchors["evidence_grade"].isin(["G1", "G2"])
        & anchors["primary_anchor_eligible"].eq("yes")
    ].copy()
    anchors["gene_id_clean"] = clean(anchors["gene_id"])
    inputs: dict[str, str] = {}
    rows: list[pd.DataFrame] = []
    for cell_class in sorted(anchors["cell_class"].unique()):
        path = args.hc3_dir / f"{cell_class}__{args.phenotype}__anchor_hc3.csv"
        table = pd.read_csv(path)
        required = {
            "gene_id_clean",
            "beta_D",
            "SE_D_moderated",
            "beta_D_HC3",
            "SE_D_HC3",
        }
        if not required <= set(table):
            raise RuntimeError(f"Missing HC3 columns in {path}")
        class_anchors = anchors.loc[
            anchors["cell_class"].eq(cell_class)
        ].merge(
            table[list(required)],
            on="gene_id_clean",
            how="left",
            validate="many_to_one",
        )
        class_anchors = class_anchors.rename(
            columns={"SE_D_moderated": "SE_D"}
        )
        class_anchors["expression_estimable"] = (
            class_anchors["beta_D"].notna() & class_anchors["SE_D"].gt(0)
        )
        class_anchors["pathology_phenotype"] = args.phenotype
        class_anchors["disease_expression_resource"] = "PsychAD RADC"
        class_anchors["replication_role"] = (
            "stage_mechanism_not_independent_genetic_replication"
        )
        rows.append(class_anchors)
        inputs[str(path)] = sha256(path)
    result = pd.concat(rows, ignore_index=True)
    result.to_csv(args.output, index=False)
    manifest = {
        "status": "COMPLETE",
        "phenotype": args.phenotype,
        "anchor_rows": int(len(result)),
        "moderated_estimable": int(result["expression_estimable"].sum()),
        "hc3_estimable": int(
            (result["beta_D_HC3"].notna() & result["SE_D_HC3"].gt(0)).sum()
        ),
        "input_sha256": inputs,
        "anchor_sha256": sha256(args.anchors),
        "output_sha256": sha256(args.output),
        "completed_at": datetime.now().astimezone().isoformat(),
    }
    manifest_path = args.output.with_suffix(".manifest.json")
    temporary = manifest_path.with_suffix(".json.tmp")
    temporary.write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    os.replace(temporary, manifest_path)
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()

