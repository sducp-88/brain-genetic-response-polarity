#!/usr/bin/env python3
"""Aggregate multiple SEA-AD subtype pseudobulks to donor-level major classes."""

from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
from scipy import sparse

from prepare_seaad_immune_core3_donor_pseudobulk import (
    CORE_REGIONS,
    EXPECTED_GENES,
    SINGLEOME_METHODS,
    decode_obs,
    decode_vector,
    prepare_metadata,
)


OBS_COLUMNS = (
    "Donor ID",
    "Brain Region",
    "Supertype",
    "method",
    "Number of nuclei",
    "Age at Death",
    "PMI",
    "Sex",
    "Race (choice=White)",
    "Braak",
    "CERAD score",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", action="append", required=True, type=Path)
    parser.add_argument("--metadata", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--audit", required=True, type=Path)
    parser.add_argument("--cell-class", required=True)
    parser.add_argument("--minimum-region-nuclei", type=int, default=20)
    parser.add_argument("--minimum-donors", type=int, default=60)
    parser.add_argument("--metadata-only", action="store_true")
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_source_metadata(
    path: Path,
    input_index: int,
    reference_gene_ids: list[str] | None,
    reference_gene_names: list[str] | None,
) -> tuple[pd.DataFrame, list[str], list[str], tuple[int, int]]:
    with h5py.File(path, "r") as handle:
        shape = tuple(int(value) for value in handle["X"].shape)
        if len(shape) != 2 or shape[1] != EXPECTED_GENES:
            raise RuntimeError(f"Unexpected X shape for {path.name}: {shape}")
        source = pd.DataFrame(
            {column: decode_obs(handle["obs"][column]) for column in OBS_COLUMNS}
        )
        source["Number of nuclei"] = pd.to_numeric(
            source["Number of nuclei"], errors="raise"
        )
        source["_input_index"] = input_index
        source["_row_in_input"] = np.arange(shape[0], dtype=np.int64)
        var_index = str(handle["var"].attrs.get("_index", "index"))
        gene_names = decode_vector(handle["var"][var_index])
        gene_ids = decode_vector(handle["var"]["gene_ids"])
        if len(set(gene_ids)) != len(gene_ids):
            raise RuntimeError(f"Duplicated gene IDs in {path.name}")
        if reference_gene_ids is not None and gene_ids != reference_gene_ids:
            raise RuntimeError(f"Gene ID order mismatch in {path.name}")
        if reference_gene_names is not None and gene_names != reference_gene_names:
            raise RuntimeError(f"Gene name order mismatch in {path.name}")
    return source, gene_ids, gene_names, shape


def main() -> None:
    args = parse_args()
    inputs = [path.resolve() for path in args.input]
    if len(inputs) < 2:
        raise RuntimeError("Multifile aggregation requires at least two inputs")
    if len(set(inputs)) != len(inputs):
        raise RuntimeError("Duplicate input paths")
    for path in inputs:
        if not path.is_file():
            raise FileNotFoundError(path)

    frames: list[pd.DataFrame] = []
    input_records: list[dict[str, object]] = []
    gene_ids: list[str] | None = None
    gene_names: list[str] | None = None
    for input_index, path in enumerate(inputs):
        source, observed_ids, observed_names, shape = read_source_metadata(
            path, input_index, gene_ids, gene_names
        )
        if gene_ids is None:
            gene_ids = observed_ids
            gene_names = observed_names
        frames.append(source)
        input_records.append(
            {
                "input_index": input_index,
                "path": str(path),
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
                "source_rows": shape[0],
                "genes": shape[1],
                "supertypes": sorted(
                    str(value) for value in source["Supertype"].dropna().unique()
                ),
            }
        )

    if gene_ids is None or gene_names is None:
        raise RuntimeError("No inputs were read")
    source_all = pd.concat(frames, ignore_index=True)
    selected = (
        source_all["method"].isin(SINGLEOME_METHODS)
        & source_all["Brain Region"].isin(CORE_REGIONS)
    )
    selected_source = source_all.loc[selected].copy()
    region_nuclei = (
        selected_source.groupby(
            ["Donor ID", "Brain Region"], observed=True
        )["Number of nuclei"]
        .sum()
        .rename("region_nuclei")
    )
    eligible_region = region_nuclei[
        region_nuclei.ge(args.minimum_region_nuclei)
    ].reset_index()
    donor_region_count = eligible_region.groupby("Donor ID")[
        "Brain Region"
    ].nunique()
    eligible_donors = sorted(
        donor_region_count[donor_region_count.eq(len(CORE_REGIONS))].index
    )
    if len(eligible_donors) < args.minimum_donors:
        raise RuntimeError(
            f"Only {len(eligible_donors)} donors have all core regions; "
            f"minimum is {args.minimum_donors}"
        )
    eligible_pairs = set(
        map(tuple, eligible_region[["Donor ID", "Brain Region"]].values)
    )
    aggregate_mask = np.asarray(
        [
            bool(is_selected)
            and donor in eligible_donors
            and (donor, region) in eligible_pairs
            for is_selected, donor, region in zip(
                selected,
                source_all["Donor ID"],
                source_all["Brain Region"],
                strict=True,
            )
        ],
        dtype=bool,
    )
    aggregate_source = source_all.loc[aggregate_mask].copy()
    donor_nuclei = (
        aggregate_source.groupby("Donor ID", observed=True)["Number of nuclei"]
        .sum()
        .reindex(eligible_donors)
    )
    public_metadata = pd.read_csv(args.metadata.resolve(), low_memory=False)
    donor_metadata = prepare_metadata(
        aggregate_source,
        public_metadata,
        eligible_donors,
        donor_nuclei,
        args.cell_class,
    )
    coverage = (
        eligible_region[eligible_region["Donor ID"].isin(eligible_donors)]
        .groupby("Brain Region", observed=True)
        .agg(donors=("Donor ID", "nunique"), nuclei=("region_nuclei", "sum"))
        .reset_index()
        .to_dict(orient="records")
    )
    per_input_selected = []
    for input_index, path in enumerate(inputs):
        subset = aggregate_source[
            aggregate_source["_input_index"].eq(input_index)
        ]
        per_input_selected.append(
            {
                "input_index": input_index,
                "path": str(path),
                "eligible_source_rows": int(len(subset)),
                "eligible_nuclei": int(subset["Number of nuclei"].sum()),
            }
        )

    audit = {
        "status": "METADATA_COMPLETE" if args.metadata_only else "COMPLETE",
        "cell_class": args.cell_class,
        "identity": f"all_{args.cell_class}_subtype_source_rows",
        "inputs": input_records,
        "source_rows_total": int(len(source_all)),
        "genes": len(gene_ids),
        "core_regions": list(CORE_REGIONS),
        "minimum_region_nuclei": args.minimum_region_nuclei,
        "eligible_donors_all_three_regions": len(eligible_donors),
        "eligible_source_rows": int(aggregate_mask.sum()),
        "per_input_selected": per_input_selected,
        "coverage": coverage,
        "braak_distribution": {
            str(key): int(value)
            for key, value in donor_metadata["Braak"]
            .value_counts()
            .sort_index()
            .items()
        },
        "cerad_distribution": {
            str(key): int(value)
            for key, value in donor_metadata["CERAD"]
            .value_counts()
            .sort_index()
            .items()
        },
        "reported_race_distribution": {
            str(key): int(value)
            for key, value in donor_metadata["reported_race_group"]
            .value_counts()
            .items()
        },
        "generated_at": datetime.now(timezone.utc).astimezone().isoformat(),
    }
    args.audit.parent.mkdir(parents=True, exist_ok=True)
    if args.metadata_only:
        args.audit.write_text(
            json.dumps(audit, indent=2, ensure_ascii=False) + "\n",
            encoding="utf-8",
        )
        print(json.dumps(audit, indent=2))
        return

    donor_index = {donor: index for index, donor in enumerate(eligible_donors)}
    counts = np.zeros((len(eligible_donors), len(gene_ids)), dtype=np.int64)
    input_selected_total = 0
    max_fractional_part = 0.0
    for input_index, path in enumerate(inputs):
        subset = aggregate_source[
            aggregate_source["_input_index"].eq(input_index)
        ]
        selected_indices = subset["_row_in_input"].to_numpy(dtype=np.int64)
        aggregate_codes = np.asarray(
            [donor_index[donor] for donor in subset["Donor ID"]],
            dtype=np.int64,
        )
        if len(selected_indices) == 0:
            continue
        with h5py.File(path, "r") as handle:
            x = handle["X"]
            block_columns = x.chunks[1] if x.chunks else 512
            for start in range(0, len(gene_ids), block_columns):
                end = min(start + block_columns, len(gene_ids))
                block = np.asarray(
                    x[selected_indices, start:end], dtype=np.float64
                )
                if np.any(block < 0):
                    raise RuntimeError(f"Negative raw counts in {path.name}")
                max_fractional_part = max(
                    max_fractional_part,
                    float(np.max(np.abs(block - np.rint(block)))),
                )
                rounded = np.rint(block).astype(np.int64)
                input_selected_total += int(rounded.sum())
                np.add.at(counts[:, start:end], aggregate_codes, rounded)
                print(
                    f"{input_index + 1}/{len(inputs)} {path.name}: "
                    f"genes {start}:{end} of {len(gene_ids)}",
                    flush=True,
                )
    if max_fractional_part > 1e-8:
        raise RuntimeError(
            f"Counts are not integer-like: {max_fractional_part}"
        )
    if int(counts.sum()) != input_selected_total:
        raise RuntimeError("Count-conservation failure")
    if counts.max() > np.iinfo(np.int32).max:
        raise RuntimeError("Count exceeds int32")

    var = pd.DataFrame(
        {"gene_name": gene_names},
        index=pd.Index(gene_ids, name="gene_id"),
    )
    output = ad.AnnData(
        X=sparse.csr_matrix(counts.astype(np.int32)),
        obs=donor_metadata,
        var=var,
    )
    output.uns["aggregation"] = {
        "unit": "donor",
        "cell_class": args.cell_class,
        "identity": f"all_{args.cell_class}_subtype_source_rows",
        "regions": list(CORE_REGIONS),
        "minimum_nuclei_per_donor_region": args.minimum_region_nuclei,
        "inputs": [path.name for path in inputs],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    output.write_h5ad(args.output, compression="gzip")
    audit.update(
        {
            "output": str(args.output.resolve()),
            "output_shape": [output.n_obs, output.n_vars],
            "input_selected_total_counts": input_selected_total,
            "output_total_counts": int(output.X.sum()),
            "count_conservation": int(output.X.sum()) == input_selected_total,
            "max_fractional_part": max_fractional_part,
            "output_sha256": sha256(args.output),
        }
    )
    args.audit.write_text(
        json.dumps(audit, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()
