#!/usr/bin/env python3
"""Build donor-by-class PsychAD pseudobulk counts from a CSR H5AD.

The full expression matrix is streamed in row blocks. Counts are accumulated
with 64-bit integers, checkpointed atomically, and exported as one H5AD per
cell class. The script never loads the complete cell-by-gene matrix into RAM.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import time
from pathlib import Path
from typing import Any

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp


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
    parser.add_argument("--cohort", default="RADC")
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("data/derived/pseudobulk/RADC"),
    )
    parser.add_argument("--batch-rows", type=int, default=8192)
    parser.add_argument("--checkpoint-every", type=int, default=20)
    parser.add_argument("--min-cells", type=int, default=20)
    parser.add_argument(
        "--keep-checkpoint",
        action="store_true",
        help="Retain the large accumulator checkpoint after successful export.",
    )
    return parser.parse_args()


def decode(value: Any) -> Any:
    if isinstance(value, bytes):
        return value.decode("utf-8", errors="replace")
    if isinstance(value, np.generic):
        return value.item()
    return value


def read_h5ad_column(group: h5py.Group, field: str) -> pd.Series:
    node = group[field]
    if isinstance(node, h5py.Group) and {"categories", "codes"} <= set(node):
        categories = [decode(value) for value in node["categories"][:]]
        codes = np.asarray(node["codes"][:], dtype=np.int64)
        return pd.Series(pd.Categorical.from_codes(codes, categories=categories))
    values = node[:]
    if values.dtype.kind in {"O", "S", "U"}:
        values = np.asarray([decode(value) for value in values])
    return pd.Series(values)


def atomic_json(path: Path, payload: dict[str, Any]) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    replace_with_retry(temporary, path)


def atomic_npy(path: Path, values: np.ndarray) -> None:
    temporary = path.with_suffix(".tmp.npy")
    np.save(temporary, values, allow_pickle=False)
    replace_with_retry(temporary, path)


def replace_with_retry(
    temporary: Path,
    destination: Path,
    attempts: int = 30,
) -> None:
    """Atomically replace a file despite brief Windows/DrvFS read locks."""
    for attempt in range(1, attempts + 1):
        try:
            os.replace(temporary, destination)
            return
        except PermissionError:
            if attempt == attempts:
                raise
            time.sleep(min(0.1 * attempt, 2.0))


def sha256(path: Path, block_size: int = 16 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(block_size):
            digest.update(block)
    return digest.hexdigest()


def load_var(var: h5py.Group) -> pd.DataFrame:
    fields = [
        field
        for field in [
            "_index",
            "feature_name",
            "gene_name",
            "feature_biotype",
            "feature_type",
            "feature_reference",
            "feature_length",
            "n_cells",
        ]
        if field in var
    ]
    frame = pd.DataFrame({field: read_h5ad_column(var, field) for field in fields})
    index_field = decode(var.attrs.get("_index", "_index"))
    if index_field not in frame:
        raise RuntimeError(f"Variable index field is absent: {index_field}")
    frame.index = frame[index_field].astype(str)
    frame.index.name = None
    frame = frame.drop(columns=[index_field])
    return frame


def h5ad_safe_dataframe(frame: pd.DataFrame) -> pd.DataFrame:
    """Convert nullable/object text columns to H5AD-compatible categoricals."""
    result = frame.copy()
    for column in result.columns:
        dtype = result[column].dtype
        if isinstance(dtype, pd.StringDtype) or pd.api.types.is_object_dtype(dtype):
            result[column] = result[column].astype("category")
    return result


def checkpoint_payload(
    h5ad: Path,
    source_bytes: int,
    n_cells: int,
    n_genes: int,
    n_donors: int,
    n_classes: int,
    completed_rows: int,
    total_counts_streamed: int,
) -> dict[str, Any]:
    return {
        "h5ad": str(h5ad),
        "source_bytes": source_bytes,
        "n_cells": n_cells,
        "n_genes": n_genes,
        "n_donors": n_donors,
        "n_classes": n_classes,
        "completed_rows": completed_rows,
        "total_counts_streamed": total_counts_streamed,
        "updated_at": pd.Timestamp.now(tz="Asia/Shanghai").isoformat(),
    }


def main() -> None:
    args = parse_args()
    if args.batch_rows < 1 or args.checkpoint_every < 1:
        raise ValueError("Batch and checkpoint sizes must be positive.")

    h5ad_path = args.h5ad.resolve()
    metadata_path = args.metadata.resolve()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    checkpoint_array = output_dir / "_class_counts_checkpoint.npy"
    checkpoint_state = output_dir / "_class_counts_checkpoint.json"

    started = time.perf_counter()
    source_bytes = h5ad_path.stat().st_size
    with h5py.File(h5ad_path, "r") as handle:
        x_group = handle["X"]
        if decode(x_group.attrs.get("encoding-type")) != "csr_matrix":
            raise RuntimeError("X must use CSR encoding.")
        n_cells, n_genes = (int(value) for value in x_group.attrs["shape"])
        donor_categories = [
            str(decode(value)) for value in handle["obs/donor_id/categories"][:]
        ]
        class_categories = [
            str(decode(value)) for value in handle["obs/class/categories"][:]
        ]
        donor_codes = np.asarray(
            handle["obs/donor_id/codes"][:], dtype=np.int64
        )
        class_codes = np.asarray(handle["obs/class/codes"][:], dtype=np.int64)
        if np.any(donor_codes < 0) or np.any(class_codes < 0):
            raise RuntimeError("Missing donor or class code detected.")
        n_donors = len(donor_categories)
        n_classes = len(class_categories)
        n_groups = n_donors * n_classes
        group_codes = donor_codes * n_classes + class_codes
        cell_counts = np.bincount(group_codes, minlength=n_groups).astype(np.int64)
        indptr = np.asarray(x_group["indptr"][:], dtype=np.int64)
        var = load_var(handle["var"])
        donor_fields = [
            field
            for field in [
                "donor_id",
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
            if field in handle["obs"]
        ]
        cell_donor_metadata = pd.DataFrame(
            {
                field: read_h5ad_column(handle["obs"], field)
                for field in donor_fields
            }
        )
        for field in donor_fields:
            if field == "donor_id":
                continue
            inconsistent = (
                cell_donor_metadata.groupby("donor_id", observed=True)[field]
                .nunique(dropna=False)
                .gt(1)
            )
            if inconsistent.any():
                raise RuntimeError(
                    f"Donor-invariant H5AD field is inconsistent: {field}"
                )
        donor_h5ad_metadata = (
            cell_donor_metadata.drop_duplicates("donor_id")
            .set_index("donor_id")
            .loc[donor_categories]
            .reset_index()
        )

        expected_state = {
            "h5ad": str(h5ad_path),
            "source_bytes": source_bytes,
            "n_cells": n_cells,
            "n_genes": n_genes,
            "n_donors": n_donors,
            "n_classes": n_classes,
        }
        if checkpoint_array.exists() != checkpoint_state.exists():
            raise RuntimeError("Incomplete checkpoint pair; inspect before rerunning.")

        if checkpoint_array.exists():
            state = json.loads(checkpoint_state.read_text(encoding="utf-8"))
            for key, expected in expected_state.items():
                if state.get(key) != expected:
                    raise RuntimeError(
                        f"Checkpoint mismatch for {key}: "
                        f"{state.get(key)!r} != {expected!r}"
                    )
            counts = np.load(checkpoint_array, allow_pickle=False)
            if counts.shape != (n_groups, n_genes) or counts.dtype != np.int64:
                raise RuntimeError("Checkpoint accumulator shape or dtype is invalid.")
            completed_rows = int(state["completed_rows"])
            total_counts_streamed = int(state["total_counts_streamed"])
            print(
                f"Resuming at row {completed_rows:,}/{n_cells:,}",
                flush=True,
            )
        else:
            counts = np.zeros((n_groups, n_genes), dtype=np.int64)
            completed_rows = 0
            total_counts_streamed = 0

        batches_since_checkpoint = 0
        for row_start in range(completed_rows, n_cells, args.batch_rows):
            row_stop = min(row_start + args.batch_rows, n_cells)
            value_start = int(indptr[row_start])
            value_stop = int(indptr[row_stop])
            raw_values = np.asarray(
                x_group["data"][value_start:value_stop], dtype=np.float32
            )
            if (
                not np.isfinite(raw_values).all()
                or np.any(raw_values < 0)
                or not np.array_equal(raw_values, np.rint(raw_values))
            ):
                raise RuntimeError(
                    f"Non-count value detected in rows {row_start}:{row_stop}."
                )
            integer_values = raw_values.astype(np.int64)
            indices = np.asarray(
                x_group["indices"][value_start:value_stop], dtype=np.int64
            )
            local_indptr = indptr[row_start : row_stop + 1] - value_start
            x_batch = sp.csr_matrix(
                (integer_values, indices, local_indptr),
                shape=(row_stop - row_start, n_genes),
            )
            local_groups = group_codes[row_start:row_stop]
            grouping = sp.csr_matrix(
                (
                    np.ones(row_stop - row_start, dtype=np.int8),
                    (local_groups, np.arange(row_stop - row_start)),
                ),
                shape=(n_groups, row_stop - row_start),
            )
            aggregated = (grouping @ x_batch).tocoo()
            np.add.at(
                counts,
                (aggregated.row, aggregated.col),
                aggregated.data.astype(np.int64, copy=False),
            )
            total_counts_streamed += int(integer_values.sum(dtype=np.int64))
            completed_rows = row_stop
            batches_since_checkpoint += 1

            if (
                batches_since_checkpoint >= args.checkpoint_every
                or completed_rows == n_cells
            ):
                atomic_npy(checkpoint_array, counts)
                state = checkpoint_payload(
                    h5ad_path,
                    source_bytes,
                    n_cells,
                    n_genes,
                    n_donors,
                    n_classes,
                    completed_rows,
                    total_counts_streamed,
                )
                atomic_json(checkpoint_state, state)
                elapsed = time.perf_counter() - started
                print(
                    f"Checkpoint: {completed_rows:,}/{n_cells:,} cells "
                    f"({100 * completed_rows / n_cells:.1f}%), "
                    f"{elapsed:.1f}s elapsed",
                    flush=True,
                )
                batches_since_checkpoint = 0

    total_counts_pseudobulk = int(counts.sum(dtype=np.int64))
    if total_counts_pseudobulk != total_counts_streamed:
        raise RuntimeError(
            "Counts conservation failed: "
            f"{total_counts_pseudobulk} != {total_counts_streamed}"
        )

    public_metadata = pd.read_csv(metadata_path, dtype={"DonorID": "string"})
    cohort_metadata = public_metadata.loc[
        public_metadata["Cohort"].astype(str).eq(args.cohort)
    ].copy()
    if set(cohort_metadata["DonorID"].astype(str)) != set(donor_categories):
        raise RuntimeError("Cohort metadata donor set does not match H5AD donors.")
    cohort_metadata = cohort_metadata.set_index("DonorID").loc[donor_categories]

    sample_rows: list[dict[str, Any]] = []
    for donor_code, donor_id in enumerate(donor_categories):
        for class_code, class_name in enumerate(class_categories):
            group_code = donor_code * n_classes + class_code
            row_counts = counts[group_code]
            sample_rows.append(
                {
                    "group_code": group_code,
                    "sample_id": f"{donor_id}__{class_name}",
                    "donor_id": donor_id,
                    "class": class_name,
                    "n_cells": int(cell_counts[group_code]),
                    "total_counts": int(row_counts.sum(dtype=np.int64)),
                    "detected_genes": int(np.count_nonzero(row_counts)),
                    "eligible_min_cells": bool(
                        cell_counts[group_code] >= args.min_cells
                    ),
                }
            )
    sample_metadata = pd.DataFrame(sample_rows)
    sample_metadata = sample_metadata.merge(
        donor_h5ad_metadata,
        how="left",
        on="donor_id",
        validate="many_to_one",
    )
    sample_metadata = sample_metadata.merge(
        cohort_metadata.reset_index(),
        how="left",
        left_on="donor_id",
        right_on="DonorID",
        validate="many_to_one",
    )
    sample_metadata.to_csv(output_dir / "sample_metadata.csv", index=False)
    var.to_csv(output_dir / "gene_metadata.csv", index=True, index_label="gene_id")

    output_hashes: dict[str, str] = {}
    for class_code, class_name in enumerate(class_categories):
        group_rows = (
            np.arange(n_donors, dtype=np.int64) * n_classes + class_code
        )
        class_samples = sample_metadata.loc[
            sample_metadata["class"].eq(class_name)
        ].copy()
        if not np.array_equal(
            class_samples["group_code"].to_numpy(dtype=np.int64), group_rows
        ):
            raise RuntimeError(f"Sample ordering failed for {class_name}.")
        matrix = sp.csr_matrix(counts[group_rows], dtype=np.int64)
        class_samples.index = class_samples["sample_id"].astype(str)
        class_samples.index.name = None
        output = ad.AnnData(
            X=matrix,
            obs=h5ad_safe_dataframe(class_samples),
            var=h5ad_safe_dataframe(var),
        )
        output.uns["source_h5ad"] = str(h5ad_path)
        output.uns["cohort"] = args.cohort
        output.uns["cell_class"] = class_name
        output.uns["minimum_cells_threshold"] = args.min_cells
        output.uns["counts_conservation_total"] = total_counts_pseudobulk
        output_path = (
            output_dir / f"{args.cohort}_{class_name}_pseudobulk_counts.h5ad"
        )
        output.write_h5ad(output_path, compression="gzip")
        output_hashes[output_path.name] = sha256(output_path)
        print(
            f"Wrote {class_name}: {matrix.shape[0]} donors × "
            f"{matrix.shape[1]} genes; {matrix.nnz:,} nonzero entries",
            flush=True,
        )

    manifest = {
        "status": "complete",
        "cohort": args.cohort,
        "source_h5ad": str(h5ad_path),
        "source_bytes": source_bytes,
        "n_cells": n_cells,
        "n_genes": n_genes,
        "n_donors": n_donors,
        "classes": class_categories,
        "observed_donor_class_groups": int(np.count_nonzero(cell_counts)),
        "total_donor_class_groups": int(n_groups),
        "minimum_cells_threshold": args.min_cells,
        "eligible_donor_class_groups": int(np.sum(cell_counts >= args.min_cells)),
        "total_counts_streamed": total_counts_streamed,
        "total_counts_pseudobulk": total_counts_pseudobulk,
        "counts_conservation_passed": True,
        "output_sha256": output_hashes,
        "elapsed_seconds": round(time.perf_counter() - started, 3),
        "completed_at": pd.Timestamp.now(tz="Asia/Shanghai").isoformat(),
    }
    atomic_json(output_dir / "manifest.json", manifest)

    if not args.keep_checkpoint:
        checkpoint_array.unlink(missing_ok=True)
        checkpoint_state.unlink(missing_ok=True)
    print(
        f"Completed pseudobulk export in {manifest['elapsed_seconds']} seconds.",
        flush=True,
    )


if __name__ == "__main__":
    main()
