#!/usr/bin/env python3
"""Inspect public CELLxGENE H5AD structure through HTTP range requests.

The script reads HDF5 metadata only. It does not download complete H5AD files
and never loads the expression matrix.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import fsspec
import h5py
import numpy as np


DATASETS = {
    "RADC_Cohort": {
        "url": "https://datasets.cellxgene.cziscience.com/54293783-669c-410e-919d-474960f8761b.h5ad",
        "published_bytes": 6_282_475_871,
    },
    "HBCC_Cohort": {
        "url": "https://datasets.cellxgene.cziscience.com/d27fb144-f105-46c2-b36f-f51421f74e4e.h5ad",
        "published_bytes": 14_150_526_668,
    },
    "Aging_Cohort": {
        "url": "https://datasets.cellxgene.cziscience.com/13e8f1dd-962f-47b7-9cf7-c71d2b21b8a5.h5ad",
        "published_bytes": 12_421_671_635,
    },
    "MSSM_Cohort": {
        "url": "https://datasets.cellxgene.cziscience.com/0e853475-e298-4b09-881a-ed0b60d5a8c9.h5ad",
        "published_bytes": 36_092_176_654,
    },
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--datasets",
        nargs="+",
        choices=sorted(DATASETS),
        default=["RADC_Cohort"],
        help="Datasets to inspect; defaults to the smallest cohort.",
    )
    parser.add_argument(
        "--output",
        default="outputs/phase0/cellxgene_remote_h5ad_schema.json",
    )
    parser.add_argument(
        "--deep",
        action="store_true",
        help="Also read category dictionaries and a small X sample; can be slow.",
    )
    return parser.parse_args()


def node_summary(node: h5py.Group | h5py.Dataset) -> dict[str, object]:
    summary: dict[str, object] = {"kind": type(node).__name__}
    if isinstance(node, h5py.Dataset):
        summary["shape"] = list(node.shape)
        summary["dtype"] = str(node.dtype)
    else:
        summary["keys"] = sorted(node.keys())
    if node.attrs:
        summary["attrs"] = {
            str(key): (
                value.tolist()
                if hasattr(value, "tolist")
                else value.decode()
                if isinstance(value, bytes)
                else str(value)
            )
            for key, value in node.attrs.items()
        }
    return summary


def json_values(values: object, limit: int = 200) -> dict[str, object]:
    array = np.asarray(values)
    flat = array.reshape(-1)
    displayed = []
    for value in flat[:limit]:
        if isinstance(value, bytes):
            displayed.append(value.decode("utf-8", errors="replace"))
        elif hasattr(value, "item"):
            displayed.append(value.item())
        else:
            displayed.append(value)
    return {
        "count": int(flat.size),
        "truncated": bool(flat.size > limit),
        "values": displayed,
    }


def inspect_dataset(name: str, deep: bool = False) -> dict[str, object]:
    spec = DATASETS[name]
    remote = fsspec.open(
        spec["url"],
        mode="rb",
        block_size=8 * 1024 * 1024,
        cache_type="blockcache",
    )
    with remote as handle, h5py.File(handle, "r") as h5:
        result: dict[str, object] = {
            "dataset": name,
            "url": spec["url"],
            "published_bytes": spec["published_bytes"],
            "root_keys": sorted(h5.keys()),
        }
        for key in ["X", "layers", "obs", "var", "obsm", "uns"]:
            if key in h5:
                result[key] = node_summary(h5[key])
        if deep and isinstance(h5.get("X"), h5py.Group) and "data" in h5["X"]:
            x_data = h5["X"]["data"]
            sample = np.asarray(x_data[: min(10_000, x_data.shape[0])])
            result["X_data"] = {
                **node_summary(x_data),
                "sample_count": int(sample.size),
                "sample_min": float(sample.min()) if sample.size else None,
                "sample_max": float(sample.max()) if sample.size else None,
                "sample_all_integer_valued": bool(
                    np.allclose(sample, np.round(sample))
                )
                if sample.size
                else None,
            }
        if deep and "obs" in h5:
            obs_fields: dict[str, object] = {}
            for field in sorted(h5["obs"].keys()):
                node = h5["obs"][field]
                field_summary = node_summary(node)
                if isinstance(node, h5py.Group) and "categories" in node:
                    field_summary["categories"] = json_values(
                        node["categories"][:]
                    )
                obs_fields[field] = field_summary
            result["obs_fields"] = obs_fields
        return result


def main() -> None:
    args = parse_args()
    results = [inspect_dataset(name, deep=args.deep) for name in args.datasets]
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(
        json.dumps(results, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    for result in results:
        obs_keys = result.get("obs", {}).get("keys", [])
        print(
            f"{result['dataset']}: root={len(result['root_keys'])} groups; "
            f"obs_columns={len(obs_keys)}"
        )
        print("  " + ", ".join(obs_keys))
    print(f"Wrote schema audit to {output}")
    print("No complete H5AD file or expression matrix was downloaded.")


if __name__ == "__main__":
    main()
