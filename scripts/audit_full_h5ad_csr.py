#!/usr/bin/env python3
"""Perform a complete sequential readability and count audit of H5AD CSR X.

Unlike the earlier evenly spaced sample audit, this program reads every stored
value and index through h5py. It therefore fails on any unreadable compressed
chunk used by the expression matrix.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import h5py
import numpy as np


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("h5ad", type=Path)
    parser.add_argument("--expected-bytes", type=int, required=True)
    parser.add_argument("--expected-cells", type=int, required=True)
    parser.add_argument("--expected-genes", type=int, required=True)
    parser.add_argument("--batch-rows", type=int, default=4096)
    parser.add_argument("--progress-every", type=int, default=10)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--status", type=Path, required=True)
    return parser.parse_args()


def now() -> str:
    return datetime.now(timezone.utc).astimezone().isoformat(timespec="seconds")


def json_safe(value: Any) -> Any:
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, float) and not math.isfinite(value):
        return None
    if isinstance(value, dict):
        return {str(key): json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(item) for item in value]
    return value


def atomic_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(text, encoding="utf-8")
    os.replace(temporary, path)


def write_report(path: Path, report: dict[str, Any]) -> None:
    atomic_text(
        path,
        json.dumps(json_safe(report), ensure_ascii=False, indent=2) + "\n",
    )


def main() -> None:
    args = parse_args()
    if args.batch_rows < 1 or args.progress_every < 1:
        raise ValueError("Batch and progress intervals must be positive.")
    h5ad = args.h5ad.resolve()
    output = args.output.resolve()
    status = args.status.resolve()
    began = time.perf_counter()
    report: dict[str, Any] = {
        "status": "RUNNING",
        "h5ad": str(h5ad),
        "started_at": now(),
        "expected": {
            "bytes": args.expected_bytes,
            "cells": args.expected_cells,
            "genes": args.expected_genes,
        },
        "batch_rows": args.batch_rows,
        "failure": None,
    }

    def fail(stage: str, error: Exception | str, **details: Any) -> None:
        report["status"] = "FAILED"
        report["failure"] = {
            "stage": stage,
            "error": repr(error) if isinstance(error, Exception) else error,
            **details,
        }
        report["elapsed_seconds"] = round(time.perf_counter() - began, 3)
        report["completed_at"] = now()
        write_report(output, report)
        atomic_text(
            status,
            f"FAILED stage={stage} error={report['failure']['error']}\n",
        )

    try:
        actual_bytes = h5ad.stat().st_size
    except Exception as error:
        fail("file_stat", error)
        raise
    report["actual_bytes"] = actual_bytes
    if actual_bytes != args.expected_bytes:
        error = RuntimeError(
            f"File size {actual_bytes} != expected {args.expected_bytes}."
        )
        fail("file_size", error)
        raise error

    try:
        with h5py.File(h5ad, "r") as handle:
            if "X" not in handle or not isinstance(handle["X"], h5py.Group):
                raise RuntimeError("Sparse X group is absent.")
            x = handle["X"]
            if not {"data", "indices", "indptr"} <= set(x):
                raise RuntimeError("CSR X is missing data, indices, or indptr.")
            encoding = x.attrs.get("encoding-type")
            if isinstance(encoding, bytes):
                encoding = encoding.decode("utf-8", errors="replace")
            if encoding != "csr_matrix":
                raise RuntimeError(f"X encoding is not CSR: {encoding!r}.")
            shape = tuple(int(value) for value in x.attrs["shape"])
            if shape != (args.expected_cells, args.expected_genes):
                raise RuntimeError(
                    f"X shape {shape} != "
                    f"{(args.expected_cells, args.expected_genes)}."
                )
            data = x["data"]
            indices = x["indices"]
            indptr = np.asarray(x["indptr"][:], dtype=np.int64)
            stored_values = int(data.shape[0])
            if indices.shape != data.shape:
                raise RuntimeError(
                    f"data shape {data.shape} != indices shape {indices.shape}."
                )
            if indptr.shape != (args.expected_cells + 1,):
                raise RuntimeError(
                    f"indptr shape {indptr.shape} is unexpected."
                )
            if indptr[0] != 0 or indptr[-1] != stored_values:
                raise RuntimeError(
                    f"CSR pointer endpoints are invalid: "
                    f"{indptr[0]}, {indptr[-1]}, nnz={stored_values}."
                )
            row_nnz = np.diff(indptr)
            if np.any(row_nnz < 0):
                raise RuntimeError("CSR indptr is not monotonic.")
            report["schema"] = {
                "root_keys": sorted(handle.keys()),
                "encoding_type": encoding,
                "shape": list(shape),
                "stored_values": stored_values,
                "data_dtype": str(data.dtype),
                "indices_dtype": str(indices.dtype),
                "indptr_dtype": str(x["indptr"].dtype),
                "row_nnz_min": int(row_nnz.min()),
                "row_nnz_median": float(np.median(row_nnz)),
                "row_nnz_max": int(row_nnz.max()),
                "empty_rows": int(np.count_nonzero(row_nnz == 0)),
            }
            del row_nnz

            total_values_read = 0
            total_count_sum = 0
            global_min = math.inf
            global_max = -math.inf
            total_batches = math.ceil(args.expected_cells / args.batch_rows)
            atomic_text(
                status,
                f"RUNNING stage=csr_scan batches=0/{total_batches} "
                f"rows=0/{args.expected_cells} values=0/{stored_values}\n",
            )

            for batch_index, row_start in enumerate(
                range(0, args.expected_cells, args.batch_rows)
            ):
                row_stop = min(row_start + args.batch_rows, args.expected_cells)
                value_start = int(indptr[row_start])
                value_stop = int(indptr[row_stop])
                try:
                    values = np.asarray(data[value_start:value_stop])
                    gene_indices = np.asarray(
                        indices[value_start:value_stop], dtype=np.int64
                    )
                except Exception as error:
                    fail(
                        "csr_scan_read",
                        error,
                        batch_index=batch_index,
                        row_start=row_start,
                        row_stop=row_stop,
                        value_start=value_start,
                        value_stop=value_stop,
                    )
                    raise
                expected_values = value_stop - value_start
                if values.size != expected_values:
                    error = RuntimeError(
                        f"Read {values.size} values, expected {expected_values}."
                    )
                    fail(
                        "csr_scan_length",
                        error,
                        batch_index=batch_index,
                        row_start=row_start,
                        row_stop=row_stop,
                    )
                    raise error
                if gene_indices.size != expected_values:
                    error = RuntimeError(
                        f"Read {gene_indices.size} indices, "
                        f"expected {expected_values}."
                    )
                    fail(
                        "csr_scan_indices_length",
                        error,
                        batch_index=batch_index,
                        row_start=row_start,
                        row_stop=row_stop,
                    )
                    raise error
                if values.size:
                    if not np.isfinite(values).all():
                        raise RuntimeError(
                            f"Non-finite count in rows {row_start}:{row_stop}."
                        )
                    if np.any(values < 0):
                        raise RuntimeError(
                            f"Negative count in rows {row_start}:{row_stop}."
                        )
                    if not np.equal(values, np.rint(values)).all():
                        raise RuntimeError(
                            f"Non-integer count in rows {row_start}:{row_stop}."
                        )
                    if np.any(gene_indices < 0) or np.any(
                        gene_indices >= args.expected_genes
                    ):
                        raise RuntimeError(
                            f"Out-of-range gene index in rows "
                            f"{row_start}:{row_stop}."
                        )
                    local_min = float(values.min())
                    local_max = float(values.max())
                    global_min = min(global_min, local_min)
                    global_max = max(global_max, local_max)
                    total_count_sum += int(
                        values.sum(dtype=np.float64)
                    )
                total_values_read += int(values.size)

                completed_batches = batch_index + 1
                if (
                    completed_batches % args.progress_every == 0
                    or row_stop == args.expected_cells
                ):
                    elapsed = max(time.perf_counter() - began, 1e-9)
                    percent = 100 * row_stop / args.expected_cells
                    message = (
                        f"RUNNING stage=csr_scan batches={completed_batches}/"
                        f"{total_batches} rows={row_stop}/"
                        f"{args.expected_cells} values={total_values_read}/"
                        f"{stored_values} percent={percent:.2f} "
                        f"elapsed_seconds={elapsed:.1f}"
                    )
                    atomic_text(status, message + "\n")
                    print(message, flush=True)

            if total_values_read != stored_values:
                raise RuntimeError(
                    f"Full scan read {total_values_read} values, "
                    f"expected {stored_values}."
                )
            report["full_csr_scan"] = {
                "batches": total_batches,
                "rows_read": args.expected_cells,
                "values_read": total_values_read,
                "all_values_read": total_values_read == stored_values,
                "all_indices_read": True,
                "counts_are_finite_nonnegative_integers": True,
                "indices_within_gene_bounds": True,
                "count_min": None if global_min == math.inf else global_min,
                "count_max": None if global_max == -math.inf else global_max,
                "total_count_sum": total_count_sum,
            }
    except Exception as error:
        if report["status"] != "FAILED":
            fail("hdf5_structure_or_value_validation", error)
        raise

    report["status"] = "COMPLETE_ALL_CSR_VALUES_AND_INDICES_READABLE"
    report["elapsed_seconds"] = round(time.perf_counter() - began, 3)
    report["completed_at"] = now()
    write_report(output, report)
    atomic_text(
        status,
        "SUCCESS status=COMPLETE_ALL_CSR_VALUES_AND_INDICES_READABLE "
        f"rows={args.expected_cells} values="
        f"{report['full_csr_scan']['values_read']}\n",
    )
    print(report["status"], flush=True)


if __name__ == "__main__":
    main()
