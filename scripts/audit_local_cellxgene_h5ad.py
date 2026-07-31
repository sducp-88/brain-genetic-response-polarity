#!/usr/bin/env python3
"""Audit a downloaded PsychAD CELLxGENE H5AD without loading X into RAM.

The audit verifies the published file size and HDF5 structure, samples the
sparse expression matrix across its full stored range, reconciles H5AD donor
identifiers with the public PsychAD donor table, and writes donor-level cell
counts needed to decide whether pseudobulk analysis is feasible.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
from typing import Any

import h5py
import numpy as np
import pandas as pd


DEFAULT_FIELDS = [
    "donor_id",
    "Source",
    "class",
    "subclass",
    "subtype",
    "development_stage",
    "sex",
    "genetic_ancestry",
    "disease",
    "AD_status",
    "ASCVD_status",
    "Bipolar_Disorder",
    "DLBD_status",
    "FTD_status",
    "Schizophrenia",
    "Tardive_dyskinesia",
    "Tauopathy_status",
    "Vascular_status",
    "Parkinson_disease",
    "n_counts",
    "n_genes",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("h5ad", type=Path)
    parser.add_argument(
        "--metadata",
        type=Path,
        default=Path(
            "external/PsychADxD/1_metadata/PsychAD_SupplementaryTable1.csv"
        ),
    )
    parser.add_argument(
        "--cohort",
        default="RADC",
        help="Cohort label in the public donor metadata.",
    )
    parser.add_argument("--expected-bytes", type=int, default=6_282_475_871)
    parser.add_argument("--expected-cells", type=int, default=693_682)
    parser.add_argument("--expected-genes", type=int, default=34_176)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("outputs/phase0/radc_local_audit"),
    )
    parser.add_argument(
        "--min-cells",
        type=int,
        default=20,
        help="Minimum cells for a donor-cell-type pseudobulk sample.",
    )
    parser.add_argument(
        "--sample-chunks",
        type=int,
        default=12,
        help="Number of evenly spaced X/data chunks to inspect.",
    )
    parser.add_argument("--sample-size", type=int, default=100_000)
    parser.add_argument(
        "--skip-sha256",
        action="store_true",
        help="Skip the full-file SHA-256 pass.",
    )
    return parser.parse_args()


def decode_scalar(value: Any) -> Any:
    if isinstance(value, bytes):
        return value.decode("utf-8", errors="replace")
    if isinstance(value, np.generic):
        return value.item()
    return value


def read_obs_column(obs: h5py.Group, field: str) -> pd.Series:
    node = obs[field]
    if isinstance(node, h5py.Group) and {"categories", "codes"} <= set(node):
        categories = [decode_scalar(value) for value in node["categories"][:]]
        codes = np.asarray(node["codes"][:], dtype=np.int64)
        return pd.Series(pd.Categorical.from_codes(codes, categories=categories))
    values = node[:]
    if values.dtype.kind in {"O", "S", "U"}:
        values = np.asarray([decode_scalar(value) for value in values])
    return pd.Series(values)


def file_sha256(path: Path, chunk_size: int = 16 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def json_safe(value: Any) -> Any:
    if value is pd.NA:
        return None
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, float) and (math.isnan(value) or math.isinf(value)):
        return None
    if isinstance(value, (list, tuple)):
        return [json_safe(item) for item in value]
    if isinstance(value, dict):
        return {str(key): json_safe(item) for key, item in value.items()}
    return value


def matrix_sample_audit(
    data: h5py.Dataset, chunk_count: int, sample_size: int
) -> dict[str, Any]:
    stored_values = int(data.shape[0])
    width = min(sample_size, stored_values)
    if stored_values == 0:
        starts = np.array([], dtype=np.int64)
    else:
        starts = np.linspace(
            0, max(stored_values - width, 0), num=chunk_count, dtype=np.int64
        )
        starts = np.unique(starts)

    observed = 0
    nonfinite = 0
    negative = 0
    noninteger = 0
    sample_min = None
    sample_max = None
    for start in starts:
        values = np.asarray(data[int(start) : int(start) + width])
        observed += int(values.size)
        finite = np.isfinite(values)
        nonfinite += int((~finite).sum())
        finite_values = values[finite]
        if not finite_values.size:
            continue
        negative += int((finite_values < 0).sum())
        noninteger += int(
            (~np.isclose(finite_values, np.rint(finite_values), atol=1e-6)).sum()
        )
        local_min = float(finite_values.min())
        local_max = float(finite_values.max())
        sample_min = local_min if sample_min is None else min(sample_min, local_min)
        sample_max = local_max if sample_max is None else max(sample_max, local_max)

    return {
        "stored_nonzero_values": stored_values,
        "sampled_chunks": int(len(starts)),
        "sampled_values": observed,
        "sample_min": sample_min,
        "sample_max": sample_max,
        "nonfinite_values": nonfinite,
        "negative_values": negative,
        "noninteger_values": noninteger,
        "sample_is_finite_nonnegative_integer_counts": (
            observed > 0 and nonfinite == 0 and negative == 0 and noninteger == 0
        ),
    }


def unique_or_flag(series: pd.Series) -> Any:
    values = series.dropna().astype(str).unique()
    if len(values) == 0:
        return pd.NA
    if len(values) == 1:
        return values[0]
    return "INCONSISTENT: " + " | ".join(sorted(values))


def summarize_level(
    obs: pd.DataFrame, level: str, min_cells: int
) -> tuple[pd.DataFrame, pd.DataFrame]:
    counts = (
        obs.groupby(["donor_id", level], observed=True)
        .size()
        .rename("n_cells")
        .reset_index()
        .sort_values([level, "donor_id"])
    )
    counts["eligible_min_cells"] = counts["n_cells"] >= min_cells
    feasibility = (
        counts.groupby(level, observed=True)
        .agg(
            total_cells=("n_cells", "sum"),
            donors_observed=("donor_id", "nunique"),
            donors_eligible=("eligible_min_cells", "sum"),
            min_cells_per_donor=("n_cells", "min"),
            median_cells_per_donor=("n_cells", "median"),
            max_cells_per_donor=("n_cells", "max"),
        )
        .reset_index()
        .sort_values("total_cells", ascending=False)
    )
    feasibility["min_cells_threshold"] = min_cells
    return counts, feasibility


def main() -> None:
    args = parse_args()
    h5ad = args.h5ad.resolve()
    metadata_path = args.metadata.resolve()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    if not h5ad.is_file():
        raise FileNotFoundError(h5ad)
    actual_bytes = h5ad.stat().st_size
    if actual_bytes != args.expected_bytes:
        raise RuntimeError(
            f"File is incomplete or unexpected: {actual_bytes} bytes; "
            f"expected {args.expected_bytes}."
        )

    audit: dict[str, Any] = {
        "h5ad": str(h5ad),
        "actual_bytes": actual_bytes,
        "expected_bytes": args.expected_bytes,
        "size_matches": actual_bytes == args.expected_bytes,
    }

    with h5py.File(h5ad, "r") as handle:
        if "X" not in handle or not isinstance(handle["X"], h5py.Group):
            raise RuntimeError("Expected sparse X group is absent.")
        x = handle["X"]
        required_sparse_nodes = {"data", "indices", "indptr"}
        if not required_sparse_nodes <= set(x):
            raise RuntimeError("X is not a complete CSR sparse matrix.")
        x_shape = [int(value) for value in x.attrs["shape"]]
        audit["root_keys"] = sorted(handle.keys())
        audit["X"] = {
            "encoding_type": decode_scalar(x.attrs.get("encoding-type")),
            "shape": x_shape,
            "data_dtype": str(x["data"].dtype),
            "indices_dtype": str(x["indices"].dtype),
            "indptr_dtype": str(x["indptr"].dtype),
            "indptr_length": int(x["indptr"].shape[0]),
            "shape_matches_expected": x_shape
            == [args.expected_cells, args.expected_genes],
        }
        audit["X"].update(
            matrix_sample_audit(
                x["data"], args.sample_chunks, args.sample_size
            )
        )
        if x_shape != [args.expected_cells, args.expected_genes]:
            raise RuntimeError(
                f"Unexpected X shape {x_shape}; expected "
                f"{[args.expected_cells, args.expected_genes]}."
            )

        obs_group = handle["obs"]
        missing_fields = [field for field in DEFAULT_FIELDS if field not in obs_group]
        if missing_fields:
            raise RuntimeError(f"Required obs fields absent: {missing_fields}")
        obs = pd.DataFrame(
            {field: read_obs_column(obs_group, field) for field in DEFAULT_FIELDS}
        )

    if len(obs) != args.expected_cells:
        raise RuntimeError(
            f"obs has {len(obs)} rows; expected {args.expected_cells}."
        )

    audit["obs"] = {
        "rows": int(len(obs)),
        "columns_read": DEFAULT_FIELDS,
        "unique_donors": int(obs["donor_id"].nunique()),
        "unique_classes": int(obs["class"].nunique()),
        "unique_subclasses": int(obs["subclass"].nunique()),
        "unique_subtypes": int(obs["subtype"].nunique()),
    }

    invariant_fields = [
        "Source",
        "development_stage",
        "sex",
        "genetic_ancestry",
        "disease",
        "AD_status",
        "ASCVD_status",
        "Bipolar_Disorder",
        "DLBD_status",
        "FTD_status",
        "Schizophrenia",
        "Tardive_dyskinesia",
        "Tauopathy_status",
        "Vascular_status",
        "Parkinson_disease",
    ]
    donor_summary = (
        obs.groupby("donor_id", observed=True)[invariant_fields]
        .agg(unique_or_flag)
        .reset_index()
    )
    donor_summary["total_cells"] = (
        obs.groupby("donor_id", observed=True).size().reindex(
            donor_summary["donor_id"]
        ).to_numpy()
    )
    inconsistent = {
        field: int(
            donor_summary[field].astype("string").str.startswith(
                "INCONSISTENT:", na=False
            ).sum()
        )
        for field in invariant_fields
    }
    audit["donor_invariant_conflicts"] = inconsistent

    metadata = pd.read_csv(metadata_path, dtype={"DonorID": "string"})
    cohort_metadata = metadata.loc[
        metadata["Cohort"].astype(str).eq(args.cohort)
    ].copy()
    h5ad_donors = set(donor_summary["donor_id"].astype(str))
    metadata_donors = set(cohort_metadata["DonorID"].astype(str))
    audit["donor_reconciliation"] = {
        "cohort": args.cohort,
        "h5ad_donors": len(h5ad_donors),
        "metadata_donors": len(metadata_donors),
        "matched_donors": len(h5ad_donors & metadata_donors),
        "h5ad_only": sorted(h5ad_donors - metadata_donors),
        "metadata_only": sorted(metadata_donors - h5ad_donors),
        "exact_set_match": h5ad_donors == metadata_donors,
    }

    donor_summary = donor_summary.merge(
        cohort_metadata,
        how="left",
        left_on="donor_id",
        right_on="DonorID",
        validate="one_to_one",
        suffixes=("_h5ad", "_metadata"),
    )
    donor_summary.to_csv(output_dir / "donor_summary_reconciled.csv", index=False)

    class_counts, class_feasibility = summarize_level(
        obs, "class", args.min_cells
    )
    subclass_counts, subclass_feasibility = summarize_level(
        obs, "subclass", args.min_cells
    )
    class_counts.to_csv(output_dir / "donor_class_cell_counts.csv", index=False)
    subclass_counts.to_csv(
        output_dir / "donor_subclass_cell_counts.csv", index=False
    )
    class_feasibility.to_csv(
        output_dir / "class_pseudobulk_feasibility.csv", index=False
    )
    subclass_feasibility.to_csv(
        output_dir / "subclass_pseudobulk_feasibility.csv", index=False
    )

    disease_flag_fields = [
        "AD_status",
        "DLBD_status",
        "FTD_status",
        "Schizophrenia",
        "Bipolar_Disorder",
        "Tauopathy_status",
        "Vascular_status",
        "ASCVD_status",
        "Parkinson_disease",
        "Tardive_dyskinesia",
    ]
    disease_count_frames = []
    for field in disease_flag_fields:
        counts = (
            donor_summary.groupby(field, dropna=False, observed=True)
            .size()
            .rename("n_donors")
            .reset_index()
            .rename(columns={field: "value"})
        )
        counts.insert(0, "field", field)
        disease_count_frames.append(counts)
    pd.concat(disease_count_frames, ignore_index=True).to_csv(
        output_dir / "donor_disease_flag_counts.csv", index=False
    )

    audit["pseudobulk"] = {
        "minimum_cells_threshold": args.min_cells,
        "classes": class_feasibility.to_dict(orient="records"),
        "subclasses": subclass_feasibility.to_dict(orient="records"),
    }

    if not args.skip_sha256:
        audit["sha256"] = file_sha256(h5ad)
    audit_path = output_dir / "audit.json"
    audit_path.write_text(
        json.dumps(json_safe(audit), ensure_ascii=False, indent=2),
        encoding="utf-8",
    )

    matrix_ok = audit["X"]["sample_is_finite_nonnegative_integer_counts"]
    donors_ok = audit["donor_reconciliation"]["exact_set_match"]
    conflicts_ok = not any(inconsistent.values())
    print(f"H5AD: {h5ad}")
    print(f"Size: {actual_bytes} bytes (published size matched)")
    print(f"X shape: {x_shape}; sampled raw-count-compatible: {matrix_ok}")
    print(
        f"Donors: {len(h5ad_donors)} H5AD / {len(metadata_donors)} metadata; "
        f"exact match: {donors_ok}"
    )
    print(f"Donor-level invariant fields conflict-free: {conflicts_ok}")
    print(f"Outputs: {output_dir}")
    if not (matrix_ok and donors_ok and conflicts_ok):
        raise SystemExit("Audit completed with at least one blocking discrepancy.")


if __name__ == "__main__":
    main()
