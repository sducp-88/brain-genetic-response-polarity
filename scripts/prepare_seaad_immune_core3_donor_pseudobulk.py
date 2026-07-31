#!/usr/bin/env python3
"""Build SEA-AD core-three-region resident-immune donor pseudobulks."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from datetime import datetime, timezone
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
from scipy import sparse


EXPECTED_BYTES = 123_541_128
EXPECTED_SHA256 = (
    "c79a9ea30281c51761da0333de90dc6ef0f4cdf525dbdac16d1cd211e74228db"
)
EXPECTED_GENES = 36_601
CORE_REGIONS = ("DFC", "MEC", "MTG")
RESIDENT_STATES = (
    "Micro-PVM_2",
    "Micro-PVM_2_1-SEAAD",
    "Micro-PVM_2_3-SEAAD",
    "Micro-PVM_3-SEAAD",
    "Micro-PVM_4-SEAAD",
)
SINGLEOME_METHODS = ("10Xv3.1", "10xV3.1_HT")
BRAAK_MAP = {
    "Braak 0": 0,
    "Braak I": 1,
    "Braak II": 2,
    "Braak III": 3,
    "Braak IV": 4,
    "Braak V": 5,
    "Braak VI": 6,
}
CERAD_MAP = {"Absent": 1, "Sparse": 2, "Moderate": 3, "Frequent": 4}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--metadata", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--audit", required=True, type=Path)
    parser.add_argument("--cell-class", default="Immune")
    parser.add_argument(
        "--identity-mode", choices=["auto", "resident", "all"], default="auto"
    )
    parser.add_argument("--expected-bytes", type=int)
    parser.add_argument("--expected-sha256")
    parser.add_argument("--minimum-region-nuclei", type=int, default=20)
    parser.add_argument("--metadata-only", action="store_true")
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def decode_vector(dataset: h5py.Dataset) -> list[str]:
    values = dataset.asstr()[...]
    return [str(value) for value in np.asarray(values).tolist()]


def decode_obs(node: h5py.Group | h5py.Dataset) -> np.ndarray:
    if isinstance(node, h5py.Dataset):
        if node.dtype.kind in {"O", "S", "U"}:
            return np.asarray(node.asstr()[...], dtype=object)
        return np.asarray(node[...])
    encoding = node.attrs.get("encoding-type", "")
    if isinstance(encoding, bytes):
        encoding = encoding.decode()
    if encoding != "categorical":
        raise RuntimeError(f"Unsupported obs encoding: {encoding!r}")
    categories = decode_vector(node["categories"])
    codes = np.asarray(node["codes"][...], dtype=int)
    return np.asarray(
        [categories[code] if code >= 0 else None for code in codes],
        dtype=object,
    )


def unique_value(values: pd.Series, field: str, donor: str) -> object:
    unique = values.dropna().astype(str).unique()
    if len(unique) != 1:
        raise RuntimeError(
            f"Expected one {field} value for donor {donor}; observed {unique}"
        )
    return unique[0]


def prepare_metadata(
    source: pd.DataFrame,
    public_metadata: pd.DataFrame,
    eligible_donors: list[str],
    donor_nuclei: pd.Series,
    cell_class: str,
) -> pd.DataFrame:
    records: list[dict[str, object]] = []
    public = public_metadata[
        public_metadata["Donor ID"].isin(eligible_donors)
    ].copy()
    if public.duplicated(["Donor ID", "Brain Region"]).any():
        raise RuntimeError("Duplicated donor-region rows in public metadata")

    for donor in eligible_donors:
        donor_source = source[source["Donor ID"].eq(donor)]
        donor_public = public[public["Donor ID"].eq(donor)]
        age = float(unique_value(donor_source["Age at Death"], "age", donor))
        pmi = float(unique_value(donor_source["PMI"], "PMI", donor))
        sex_raw = unique_value(donor_source["Sex"], "sex", donor)
        sex = {"Female": "F", "Male": "M"}.get(sex_raw)
        if sex is None:
            raise RuntimeError(f"Unexpected sex value for {donor}: {sex_raw}")
        white = unique_value(
            donor_source["Race (choice=White)"], "reported race", donor
        )
        race_group = "White" if white == "Checked" else "Other"
        braak_raw = unique_value(donor_source["Braak"], "Braak", donor)
        cerad_raw = unique_value(
            donor_source["CERAD score"], "CERAD", donor
        )
        if braak_raw not in BRAAK_MAP or cerad_raw not in CERAD_MAP:
            raise RuntimeError(
                f"Unmapped pathology for {donor}: {braak_raw}, {cerad_raw}"
            )
        global_ptau = donor_public["CPS_Global_pTau"].dropna().astype(float)
        if global_ptau.empty:
            raise RuntimeError(f"Missing global pTau CPS for {donor}")
        if global_ptau.max() - global_ptau.min() > 1e-10:
            raise RuntimeError(f"Global pTau CPS varies within donor {donor}")
        records.append(
            {
                "sample_id": f"SEAAD_{donor}_{cell_class}_core3",
                "donor_id": donor,
                "n_cells": int(donor_nuclei.loc[donor]),
                "n_regions": len(CORE_REGIONS),
                "regions": ";".join(CORE_REGIONS),
                "Age": age,
                "Sex": sex,
                "PMI": pmi,
                "reported_race_group": race_group,
                "Braak": BRAAK_MAP[braak_raw],
                "CERAD": CERAD_MAP[cerad_raw],
                "CPS_Global_pTau": float(global_ptau.iloc[0]),
                "Braak_original": braak_raw,
                "CERAD_original": cerad_raw,
            }
        )
    result = pd.DataFrame.from_records(records).set_index("sample_id")
    if result.index.duplicated().any():
        raise RuntimeError("Duplicated output sample identifiers")
    return result


def main() -> None:
    args = parse_args()
    input_path = args.input.resolve()
    metadata_path = args.metadata.resolve()
    expected_bytes = args.expected_bytes
    expected_sha256 = args.expected_sha256
    if args.cell_class == "Immune":
        expected_bytes = expected_bytes or EXPECTED_BYTES
        expected_sha256 = expected_sha256 or EXPECTED_SHA256
    if expected_bytes is not None and input_path.stat().st_size != expected_bytes:
        raise RuntimeError(
            f"Unexpected input size {input_path.stat().st_size}; "
            f"expected {expected_bytes}"
        )
    input_hash = sha256(input_path)
    if expected_sha256 is not None and input_hash != expected_sha256:
        raise RuntimeError(f"Unexpected input SHA256: {input_hash}")

    with h5py.File(input_path, "r") as handle:
        source_shape = tuple(int(value) for value in handle["X"].shape)
        if len(source_shape) != 2 or source_shape[1] != EXPECTED_GENES:
            raise RuntimeError(f"Unexpected X shape: {source_shape}")
        obs = handle["obs"]
        obs_columns = [
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
        ]
        source = pd.DataFrame(
            {column: decode_obs(obs[column]) for column in obs_columns}
        )
        source["Number of nuclei"] = pd.to_numeric(
            source["Number of nuclei"], errors="raise"
        )
        identity_mode = args.identity_mode
        if identity_mode == "auto":
            identity_mode = "resident" if args.cell_class == "Immune" else "all"
        identity_selected = (
            source["Supertype"].isin(RESIDENT_STATES)
            if identity_mode == "resident"
            else pd.Series(True, index=source.index)
        )
        identity_label = (
            "resident_microglia_brain_macrophage"
            if identity_mode == "resident"
            else f"all_{args.cell_class}_source_rows"
        )
        selected = (
            identity_selected
            & source["method"].isin(SINGLEOME_METHODS)
            & source["Brain Region"].isin(CORE_REGIONS)
        )
        selected_source = source.loc[selected].copy()
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
        if len(eligible_donors) < 60:
            raise RuntimeError(
                f"Only {len(eligible_donors)} donors have all core regions"
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
                    source["Donor ID"],
                    source["Brain Region"],
                    strict=True,
                )
            ],
            dtype=bool,
        )
        donor_index = {donor: index for index, donor in enumerate(eligible_donors)}
        aggregate_codes = np.asarray(
            [
                donor_index[donor]
                for donor in source.loc[aggregate_mask, "Donor ID"]
            ],
            dtype=np.int64,
        )
        donor_nuclei = (
            source.loc[aggregate_mask]
            .groupby("Donor ID", observed=True)["Number of nuclei"]
            .sum()
            .reindex(eligible_donors)
        )

        public_metadata = pd.read_csv(metadata_path, low_memory=False)
        donor_metadata = prepare_metadata(
            source.loc[aggregate_mask],
            public_metadata,
            eligible_donors,
            donor_nuclei,
            args.cell_class,
        )

        var_index = str(handle["var"].attrs.get("_index", "index"))
        gene_names = decode_vector(handle["var"][var_index])
        gene_ids = decode_vector(handle["var"]["gene_ids"])
        if len(set(gene_ids)) != len(gene_ids):
            raise RuntimeError("Duplicated Ensembl gene identifiers")

        coverage = (
            eligible_region[
                eligible_region["Donor ID"].isin(eligible_donors)
            ]
            .groupby("Brain Region", observed=True)
            .agg(
                donors=("Donor ID", "nunique"),
                nuclei=("region_nuclei", "sum"),
            )
            .reset_index()
            .to_dict(orient="records")
        )
        base_audit = {
            "status": "METADATA_COMPLETE"
            if args.metadata_only else "COMPLETE",
            "input": str(input_path),
            "input_bytes": input_path.stat().st_size,
            "input_sha256": input_hash,
            "source_shape": list(source_shape),
            "cell_class": args.cell_class,
            "identity": identity_label,
            "core_regions": list(CORE_REGIONS),
            "minimum_region_nuclei": args.minimum_region_nuclei,
            "eligible_donors_all_three_regions": len(eligible_donors),
            "eligible_source_rows": int(aggregate_mask.sum()),
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
                json.dumps(base_audit, indent=2, ensure_ascii=False) + "\n",
                encoding="utf-8",
            )
            print(json.dumps(base_audit, indent=2))
            return

        counts = np.zeros(
            (len(eligible_donors), source_shape[1]), dtype=np.int64
        )
        selected_indices = np.flatnonzero(aggregate_mask)
        input_selected_total = 0
        x = handle["X"]
        block_columns = x.chunks[1] if x.chunks else 512
        max_fractional_part = 0.0
        for start in range(0, source_shape[1], block_columns):
            end = min(start + block_columns, source_shape[1])
            block = np.asarray(
                x[selected_indices, start:end], dtype=np.float64
            )
            if np.any(block < 0):
                raise RuntimeError("Negative raw counts")
            max_fractional_part = max(
                max_fractional_part,
                float(np.max(np.abs(block - np.rint(block)))),
            )
            rounded = np.rint(block).astype(np.int64)
            input_selected_total += int(rounded.sum())
            np.add.at(counts[:, start:end], aggregate_codes, rounded)
            print(
                f"Aggregated genes {start}:{end} of {source_shape[1]}",
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
        "identity": identity_label,
        "regions": list(CORE_REGIONS),
        "minimum_nuclei_per_donor_region": args.minimum_region_nuclei,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    output.write_h5ad(args.output, compression="gzip")

    base_audit.update(
        {
            "output": str(args.output.resolve()),
            "output_shape": [output.n_obs, output.n_vars],
            "input_selected_total_counts": input_selected_total,
            "output_total_counts": int(output.X.sum()),
            "count_conservation": (
                int(output.X.sum()) == input_selected_total
            ),
            "max_fractional_part": max_fractional_part,
            "output_sha256": sha256(args.output),
        }
    )
    args.audit.write_text(
        json.dumps(base_audit, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(base_audit, indent=2))


if __name__ == "__main__":
    main()
