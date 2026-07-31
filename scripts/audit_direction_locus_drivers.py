#!/usr/bin/env python3
"""Audit locus influence and named complex regions in coupling output."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path


APOE_CLUSTER = {
    "APOE", "TOMM40", "APOC1", "APOC2", "APOC4", "NECTIN2", "PVRL2"
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--anchor-probabilities", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8-sig", newline="")
    return path.open(encoding="utf-8-sig", newline="")


def named_region(gene: str, locus: str, annotated: str) -> str:
    symbols = {
        item.strip().upper()
        for item in gene.replace("|", ";").split(";")
        if item.strip()
    }
    locus_upper = locus.upper()
    if symbols & APOE_CLUSTER or "APOE" in locus_upper:
        return "APOE_CLUSTER"
    if (
        any(item.startswith("HLA") for item in symbols)
        or "MHC" in locus_upper
    ):
        return "MHC"
    if any(item.startswith("MS4A") for item in symbols) or "MS4A" in locus_upper:
        return "MS4A_CLUSTER"
    if annotated.strip():
        return f"SOURCE_ANNOTATED:{annotated.strip()}"
    return ""


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    with open_text(args.anchor_probabilities) as handle:
        rows = [
            row for row in csv.DictReader(handle)
            if row.get("expression_estimable", "").lower() == "true"
        ]
    if not rows:
        raise ValueError("No estimable anchor rows")

    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        grouped[row["locus_id"]].append(row)

    locus_score = {
        locus: sum(float(row["S_moderated_rho_0"]) for row in members)
        / len(members)
        for locus, members in grouped.items()
    }
    denominator = sum(abs(score) for score in locus_score.values())
    n_loci = len(locus_score)
    overall = sum(locus_score.values()) / n_loci
    output_rows: list[dict[str, object]] = []
    for locus, members in grouped.items():
        genes = sorted({row["gene_symbol"] for row in members})
        cells = sorted({row["cell_class"] for row in members})
        source_annotations = sorted(
            {
                row.get("complex_region", "").strip()
                for row in members
                if row.get("complex_region", "").strip()
            }
        )
        tag = named_region(
            ";".join(genes), locus, ";".join(source_annotations)
        )
        loo = (
            sum(
                score for other, score in locus_score.items()
                if other != locus
            )
            / (n_loci - 1)
        )
        score = locus_score[locus]
        output_rows.append(
            {
                "locus_id": locus,
                "genes": ";".join(genes),
                "cell_classes": ";".join(cells),
                "anchor_rows": len(members),
                "locus_mean_S": score,
                "absolute_contribution_share": (
                    abs(score) / denominator if denominator else 0.0
                ),
                "leave_this_locus_out_S": loo,
                "named_complex_region": tag,
            }
        )
    output_rows.sort(
        key=lambda row: float(row["absolute_contribution_share"]),
        reverse=True,
    )

    table_path = args.output_dir / "locus_driver_audit.csv"
    with table_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(output_rows[0]))
        writer.writeheader()
        writer.writerows(output_rows)

    named = [
        row for row in output_rows if row["named_complex_region"]
    ]
    summary = {
        "status": "COMPLETE",
        "input": str(args.anchor_probabilities),
        "input_sha256": sha256(args.anchor_probabilities),
        "anchor_rows": len(rows),
        "loci": n_loci,
        "overall_S_reconstructed": overall,
        "leave_one_locus_min": min(
            float(row["leave_this_locus_out_S"]) for row in output_rows
        ),
        "leave_one_locus_max": max(
            float(row["leave_this_locus_out_S"]) for row in output_rows
        ),
        "leave_one_locus_crosses_zero": (
            min(float(row["leave_this_locus_out_S"]) for row in output_rows)
            <= 0
            <= max(float(row["leave_this_locus_out_S"]) for row in output_rows)
        ),
        "maximum_absolute_locus_contribution_share": max(
            float(row["absolute_contribution_share"]) for row in output_rows
        ),
        "maximum_contribution_locus": output_rows[0]["locus_id"],
        "named_complex_regions": named,
        "apoe_or_mhc_present": any(
            row["named_complex_region"] in {"APOE_CLUSTER", "MHC"}
            for row in named
        ),
        "coordinate_based_complex_ld_audit": "NOT_AVAILABLE",
        "interpretation": (
            "Named APOE/MHC/MS4A and source-annotated regions were audited. "
            "This does not replace a coordinate-based complex-LD audit."
        ),
        "generated_at": datetime.now(timezone.utc).astimezone().isoformat(),
        "output_sha256": {},
    }
    summary["output_sha256"]["locus_driver_audit.csv"] = sha256(table_path)
    summary_path = args.output_dir / "locus_driver_audit.json"
    summary_path.write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    print(json.dumps({
        "loci": n_loci,
        "overall_S": overall,
        "maximum_share": summary[
            "maximum_absolute_locus_contribution_share"
        ],
        "apoe_or_mhc_present": summary["apoe_or_mhc_present"],
    }))


if __name__ == "__main__":
    main()
