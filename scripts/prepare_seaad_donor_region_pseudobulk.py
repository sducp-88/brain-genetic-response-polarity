#!/usr/bin/env python3
"""Build auditable SEA-AD donor-by-region pseudobulks.

The analysis unit is one donor, one fixed major cell class, and one fixed brain
region. DFC is the matched-region primary analysis; MEC and MTG are frozen
regional sensitivity analyses. Composition principal components are derived
without using pathology.
"""

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
    BRAAK_MAP,
    CERAD_MAP,
    EXPECTED_GENES,
    RESIDENT_STATES,
    SINGLEOME_METHODS,
    decode_obs,
    decode_vector,
    unique_value,
)


FROZEN_REGIONS = ("DFC", "MEC", "MTG")
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
    parser.add_argument("--composition-table", required=True, type=Path)
    parser.add_argument("--cell-class", required=True)
    parser.add_argument("--minimum-region-nuclei", type=int, default=20)
    parser.add_argument("--minimum-donors-per-region", type=int, default=60)
    parser.add_argument("--metadata-only", action="store_true")
    return parser.parse_args()


def sha256(path: Path, block_size: int = 8 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        while block := stream.read(block_size):
            digest.update(block)
    return digest.hexdigest()


def timestamp() -> str:
    return datetime.now(timezone.utc).astimezone().isoformat(timespec="seconds")


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
        source["_source_file"] = path.name
        source["_composition_label"] = (
            "input"
            + str(input_index).zfill(2)
            + ":"
            + source["Supertype"].fillna("missing").astype(str)
        )

        var_index = str(handle["var"].attrs.get("_index", "index"))
        gene_names = decode_vector(handle["var"][var_index])
        gene_ids = decode_vector(handle["var"]["gene_ids"])
        if len(set(gene_ids)) != len(gene_ids):
            raise RuntimeError(f"Duplicated gene IDs in {path.name}")
        if reference_gene_ids is not None and gene_ids != reference_gene_ids:
            raise RuntimeError(f"Gene ID order mismatch in {path.name}")
        if reference_gene_names is not None and gene_names != reference_gene_names:
            raise RuntimeError(f"Gene-name order mismatch in {path.name}")
    return source, gene_ids, gene_names, shape


def orient_svd(loadings: np.ndarray, scores: np.ndarray) -> None:
    for component in range(loadings.shape[0]):
        pivot = int(np.argmax(np.abs(loadings[component, :])))
        if loadings[component, pivot] < 0:
            loadings[component, :] *= -1
            scores[:, component] *= -1


def composition_pcs(
    eligible_source: pd.DataFrame,
    sample_pairs: pd.DataFrame,
) -> tuple[pd.DataFrame, dict[str, object]]:
    counts = (
        eligible_source.groupby(
            ["Donor ID", "Brain Region", "_composition_label"],
            observed=True,
        )["Number of nuclei"]
        .sum()
        .unstack(fill_value=0)
    )
    expected_index = pd.MultiIndex.from_frame(
        sample_pairs[["Donor ID", "Brain Region"]]
    )
    counts = counts.reindex(expected_index, fill_value=0).sort_index(axis=1)
    if counts.shape[1] == 0:
        raise RuntimeError("No composition labels are available")

    fractions = counts.div(counts.sum(axis=1), axis=0)
    composition = pd.concat(
        [
            sample_pairs.set_index(["Donor ID", "Brain Region"])[
                ["sample_id"]
            ],
            counts.add_prefix("count__"),
            fractions.add_prefix("fraction__"),
        ],
        axis=1,
    ).reset_index()
    for index in range(3):
        composition[f"composition_PC{index + 1}"] = 0.0

    audit: dict[str, object] = {
        "pseudocount": 0.5,
        "transform": "centered_log_ratio",
        "rare_rule": (
            "merge into Other if total-nucleus fraction <0.01 or "
            "prevalence among eligible donors <0.50, within each region"
        ),
        "selection_rule": (
            "smallest K explaining at least 80% CLR variance, capped at 3"
        ),
        "pathology_used": False,
        "regions": {},
    }

    for region in FROZEN_REGIONS:
        region_counts = counts.loc[
            counts.index.get_level_values("Brain Region") == region
        ].copy()
        if region_counts.empty:
            raise RuntimeError(f"No composition rows for {region}")
        abundance = region_counts.sum(axis=0) / float(
            region_counts.to_numpy().sum()
        )
        prevalence = region_counts.gt(0).mean(axis=0)
        retained = abundance.ge(0.01) & prevalence.ge(0.50)
        retained_labels = list(region_counts.columns[retained])
        rare_labels = list(region_counts.columns[~retained])
        collapsed = region_counts.loc[:, retained_labels].copy()
        if rare_labels:
            collapsed["Other"] = region_counts.loc[:, rare_labels].sum(axis=1)
        collapsed = collapsed.loc[:, collapsed.sum(axis=0).gt(0)]

        region_record: dict[str, object] = {
            "eligible_donors": int(len(region_counts)),
            "original_feature_labels": [
                str(value) for value in region_counts.columns
            ],
            "retained_feature_labels": [
                str(value) for value in retained_labels
            ],
            "rare_feature_labels": [
                str(value) for value in rare_labels
            ],
            "collapsed_feature_labels": [
                str(value) for value in collapsed.columns
            ],
            "abundance": {
                str(label): float(value)
                for label, value in abundance.items()
            },
            "prevalence": {
                str(label): float(value)
                for label, value in prevalence.items()
            },
        }
        if collapsed.shape[1] < 2:
            region_record.update(
                {
                    "n_components": 0,
                    "explained_variance_ratio": [],
                    "cumulative_explained_variance": 0.0,
                    "loadings": [],
                }
            )
            audit["regions"][region] = region_record
            continue

        log_counts = np.log(collapsed.to_numpy(dtype=float) + 0.5)
        clr = log_counts - log_counts.mean(axis=1, keepdims=True)
        centered = clr - clr.mean(axis=0, keepdims=True)
        _, singular_values, vt = np.linalg.svd(
            centered, full_matrices=False
        )
        total_variance = float(np.sum(singular_values**2))
        if total_variance <= 0:
            region_record.update(
                {
                    "n_components": 0,
                    "explained_variance_ratio": [],
                    "cumulative_explained_variance": 0.0,
                    "loadings": [],
                }
            )
            audit["regions"][region] = region_record
            continue

        variance_ratio = singular_values**2 / total_variance
        k80 = int(np.searchsorted(np.cumsum(variance_ratio), 0.80) + 1)
        n_components = min(
            3,
            k80,
            collapsed.shape[1] - 1,
            len(collapsed) - 1,
        )
        scores = centered @ vt[:n_components, :].T
        loadings = vt[:n_components, :].copy()
        orient_svd(loadings, scores)
        region_rows = composition["Brain Region"].eq(region)
        if int(region_rows.sum()) != len(scores):
            raise RuntimeError(f"Composition row mismatch for {region}")
        for index in range(n_components):
            composition.loc[
                region_rows, f"composition_PC{index + 1}"
            ] = scores[:, index]
        region_record.update(
            {
                "n_components": int(n_components),
                "explained_variance_ratio": [
                    float(value)
                    for value in variance_ratio[:n_components]
                ],
                "cumulative_explained_variance": float(
                    variance_ratio[:n_components].sum()
                ),
                "loadings": [
                    {
                        "component": f"PC{component + 1}",
                        "values": {
                            str(label): float(value)
                            for label, value in zip(
                                collapsed.columns,
                                loadings[component, :],
                                strict=True,
                            )
                        },
                    }
                    for component in range(n_components)
                ],
            }
        )
        audit["regions"][region] = region_record
    return composition, audit


def prepare_sample_metadata(
    eligible_source: pd.DataFrame,
    public_metadata: pd.DataFrame,
    sample_pairs: pd.DataFrame,
    composition: pd.DataFrame,
    cell_class: str,
) -> pd.DataFrame:
    public = public_metadata[
        public_metadata["Donor ID"].isin(sample_pairs["Donor ID"])
        & public_metadata["Brain Region"].isin(FROZEN_REGIONS)
    ].copy()
    if public.duplicated(["Donor ID", "Brain Region"]).any():
        raise RuntimeError("Duplicated donor-region rows in public metadata")

    records: list[dict[str, object]] = []
    for donor_value, region_value, n_cells_value, sample_id_value in (
        sample_pairs[
            ["Donor ID", "Brain Region", "n_cells", "sample_id"]
        ].itertuples(index=False, name=None)
    ):
        donor = str(donor_value)
        region = str(region_value)
        sample_id = str(sample_id_value)
        n_cells = int(n_cells_value)
        sample_source = eligible_source[
            eligible_source["Donor ID"].eq(donor)
            & eligible_source["Brain Region"].eq(region)
        ]
        donor_source = eligible_source[eligible_source["Donor ID"].eq(donor)]
        donor_public = public[public["Donor ID"].eq(donor)]
        if donor_public.empty:
            raise RuntimeError(f"Missing public metadata for donor {donor}")

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
                "sample_id": sample_id,
                "donor_id": donor,
                "region": region,
                "cell_class": cell_class,
                "n_cells": n_cells,
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

    metadata = pd.DataFrame.from_records(records)
    pc_columns = [
        column for column in composition.columns
        if column.startswith("composition_PC")
    ]
    metadata = metadata.merge(
        composition[["sample_id", *pc_columns]],
        on="sample_id",
        how="left",
        validate="one_to_one",
    )
    if pc_columns and metadata[pc_columns].isna().any().any():
        raise RuntimeError("Missing composition PC values after metadata merge")
    metadata = metadata.set_index("sample_id")
    if metadata.index.duplicated().any():
        raise RuntimeError("Duplicated output sample identifiers")
    return metadata


def main() -> None:
    args = parse_args()
    inputs = [path.resolve() for path in args.input]
    if not inputs or len(set(inputs)) != len(inputs):
        raise RuntimeError("Inputs must be unique and non-empty")
    for path in inputs:
        if not path.is_file():
            raise FileNotFoundError(path)
    metadata_path = args.metadata.resolve()
    if not metadata_path.is_file():
        raise FileNotFoundError(metadata_path)

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
        raise RuntimeError("No source data were read")

    source_all = pd.concat(frames, ignore_index=True)
    identity_selected = (
        source_all["Supertype"].isin(RESIDENT_STATES)
        if args.cell_class == "Immune"
        else pd.Series(True, index=source_all.index)
    )
    identity_label = (
        "resident_microglia_brain_macrophage"
        if args.cell_class == "Immune"
        else f"all_{args.cell_class}_source_rows"
    )
    selected = (
        identity_selected
        & source_all["method"].isin(SINGLEOME_METHODS)
        & source_all["Brain Region"].isin(FROZEN_REGIONS)
    )
    selected_source = source_all.loc[selected].copy()
    region_nuclei = (
        selected_source.groupby(
            ["Donor ID", "Brain Region"], observed=True
        )["Number of nuclei"]
        .sum()
        .rename("n_cells")
        .reset_index()
    )
    sample_pairs = region_nuclei[
        region_nuclei["n_cells"].ge(args.minimum_region_nuclei)
    ].copy()
    sample_pairs = sample_pairs.sort_values(
        ["Brain Region", "Donor ID"]
    ).reset_index(drop=True)
    sample_pairs["sample_id"] = [
        f"SEAAD_{donor}_{args.cell_class}_{region}"
        for donor, region in sample_pairs[
            ["Donor ID", "Brain Region"]
        ].itertuples(index=False, name=None)
    ]

    donors_by_region = (
        sample_pairs.groupby("Brain Region", observed=True)["Donor ID"]
        .nunique()
        .reindex(FROZEN_REGIONS, fill_value=0)
    )
    deficient = donors_by_region[
        donors_by_region.lt(args.minimum_donors_per_region)
    ]
    if not deficient.empty:
        raise RuntimeError(
            "Insufficient donors by region: "
            + ", ".join(
                f"{region}={count}" for region, count in deficient.items()
            )
        )
    eligible_pairs = set(
        map(tuple, sample_pairs[["Donor ID", "Brain Region"]].values)
    )
    aggregate_mask = np.asarray(
        [
            bool(is_selected) and (donor, region) in eligible_pairs
            for is_selected, donor, region in zip(
                selected,
                source_all["Donor ID"],
                source_all["Brain Region"],
                strict=True,
            )
        ],
        dtype=bool,
    )
    eligible_source = source_all.loc[aggregate_mask].copy()
    composition, composition_audit = composition_pcs(
        eligible_source, sample_pairs
    )
    public_metadata = pd.read_csv(metadata_path, low_memory=False)
    sample_metadata = prepare_sample_metadata(
        eligible_source,
        public_metadata,
        sample_pairs,
        composition,
        args.cell_class,
    )

    donor_region_counts = (
        sample_pairs.groupby("Donor ID", observed=True)["Brain Region"]
        .nunique()
        .value_counts()
        .sort_index()
    )
    audit: dict[str, object] = {
        "status": "METADATA_COMPLETE" if args.metadata_only else "COMPLETE",
        "generated_at": timestamp(),
        "cell_class": args.cell_class,
        "identity": identity_label,
        "frozen_regions": list(FROZEN_REGIONS),
        "minimum_region_nuclei": args.minimum_region_nuclei,
        "minimum_donors_per_region": args.minimum_donors_per_region,
        "inputs": input_records,
        "metadata": {
            "path": str(metadata_path),
            "bytes": metadata_path.stat().st_size,
            "sha256": sha256(metadata_path),
        },
        "source_rows_total": int(len(source_all)),
        "eligible_source_rows": int(len(eligible_source)),
        "eligible_donor_region_samples": int(len(sample_pairs)),
        "unique_donors": int(sample_pairs["Donor ID"].nunique()),
        "donors_by_region": {
            str(region): int(value)
            for region, value in donors_by_region.items()
        },
        "donor_region_coverage": {
            str(regions): int(donors)
            for regions, donors in donor_region_counts.items()
        },
        "nuclei_by_region": {
            str(key): int(value)
            for key, value in sample_pairs.groupby(
                "Brain Region", observed=True
            )["n_cells"].sum().items()
        },
        "braak_distribution_by_region": {
            str(region): {
                str(key): int(value)
                for key, value in group["Braak"]
                .value_counts()
                .sort_index()
                .items()
            }
            for region, group in sample_metadata.groupby(
                "region", observed=True
            )
        },
        "cerad_distribution_by_region": {
            str(region): {
                str(key): int(value)
                for key, value in group["CERAD"]
                .value_counts()
                .sort_index()
                .items()
            }
            for region, group in sample_metadata.groupby(
                "region", observed=True
            )
        },
        "composition_pca": composition_audit,
    }

    args.audit.parent.mkdir(parents=True, exist_ok=True)
    args.composition_table.parent.mkdir(parents=True, exist_ok=True)
    composition.to_csv(args.composition_table, index=False)
    if args.metadata_only:
        args.audit.write_text(
            json.dumps(audit, indent=2, ensure_ascii=False) + "\n",
            encoding="utf-8",
        )
        print(json.dumps(audit, indent=2))
        return

    sample_index = {
        (row.donor_id, row.region): index
        for index, row in enumerate(sample_metadata.itertuples())
    }
    counts = np.zeros(
        (len(sample_metadata), len(gene_ids)), dtype=np.int64
    )
    input_selected_total = 0
    max_fractional_part = 0.0
    for input_index, path in enumerate(inputs):
        subset = eligible_source[
            eligible_source["_input_index"].eq(input_index)
        ]
        selected_indices = subset["_row_in_input"].to_numpy(dtype=np.int64)
        aggregate_codes = np.asarray(
            [
                sample_index[(donor, region)]
                for donor, region in subset[
                    ["Donor ID", "Brain Region"]
                ].itertuples(index=False, name=None)
            ],
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
        obs=sample_metadata,
        var=var,
    )
    output.uns["aggregation"] = {
        "unit": "donor_region",
        "cell_class": args.cell_class,
        "identity": identity_label,
        "regions": list(FROZEN_REGIONS),
        "minimum_nuclei_per_donor_region": args.minimum_region_nuclei,
        "composition_pca_json": json.dumps(
            composition_audit,
            ensure_ascii=False,
            sort_keys=True,
        ),
        "inputs": [path.name for path in inputs],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    temporary_output = args.output.with_name(
        args.output.name + ".tmp.h5ad"
    )
    output.write_h5ad(temporary_output, compression="gzip")
    temporary_output.replace(args.output)
    audit.update(
        {
            "output": str(args.output.resolve()),
            "output_shape": [output.n_obs, output.n_vars],
            "input_selected_total_counts": input_selected_total,
            "output_total_counts": int(output.X.sum()),
            "count_conservation": int(output.X.sum())
            == input_selected_total,
            "max_fractional_part": max_fractional_part,
            "output_sha256": sha256(args.output),
            "composition_table": {
                "path": str(args.composition_table.resolve()),
                "sha256": sha256(args.composition_table),
            },
        }
    )
    args.audit.write_text(
        json.dumps(audit, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()
