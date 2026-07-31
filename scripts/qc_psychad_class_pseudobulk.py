#!/usr/bin/env python3
"""QC donor-by-class PsychAD pseudobulk counts before disease modeling."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
from typing import Any

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.sparse as sp
from sklearn.decomposition import PCA


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--pseudobulk-dir",
        type=Path,
        required=True,
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        required=True,
    )
    parser.add_argument(
        "--cohort",
        default=None,
        help="Expected cohort label; defaults to the pseudobulk manifest value.",
    )
    parser.add_argument(
        "--color-field",
        default="Tier1_crossDis_dx",
        help="Donor metadata field used only to color diagnostic QC plots.",
    )
    parser.add_argument("--minimum-cells", type=int, default=20)
    parser.add_argument("--cpm-threshold", type=float, default=1.0)
    parser.add_argument("--minimum-fraction", type=float, default=0.10)
    parser.add_argument("--pca-genes", type=int, default=3000)
    return parser.parse_args()


def sha256(path: Path, block_size: int = 16 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(block_size):
            digest.update(block)
    return digest.hexdigest()


def robust_z(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    median = np.nanmedian(values)
    mad = np.nanmedian(np.abs(values - median))
    if not np.isfinite(mad) or mad == 0:
        return np.zeros_like(values)
    return 0.67448975 * (values - median) / mad


def numeric_age(values: pd.Series) -> pd.Series:
    return pd.to_numeric(
        values.astype(str).str.replace("+", "", regex=False),
        errors="coerce",
    )


def json_safe(value: Any) -> Any:
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, float) and not math.isfinite(value):
        return None
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, dict):
        return {str(key): json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(item) for item in value]
    return value


def main() -> None:
    args = parse_args()
    pseudobulk_dir = args.pseudobulk_dir.resolve()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    manifest_path = pseudobulk_dir / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    expected_hashes = manifest["output_sha256"]
    class_order = manifest["classes"]
    cohort = args.cohort or manifest["cohort"]
    if manifest["cohort"] != cohort:
        raise RuntimeError(
            f"Cohort mismatch: manifest={manifest['cohort']!r}, "
            f"requested={cohort!r}."
        )

    summary_rows: list[dict[str, Any]] = []
    sample_qc_frames: list[pd.DataFrame] = []
    pca_frames: list[pd.DataFrame] = []
    filter_grid_rows: list[dict[str, Any]] = []
    observed_total_counts = 0

    plot_columns = min(4, max(1, len(class_order)))
    plot_rows = int(math.ceil(len(class_order) / plot_columns))
    pca_figure, pca_axes = plt.subplots(
        plot_rows,
        plot_columns,
        figsize=(4.5 * plot_columns, 4.5 * plot_rows),
        squeeze=False,
    )
    library_figure, library_axes = plt.subplots(
        plot_rows,
        plot_columns,
        figsize=(4.5 * plot_columns, 4.5 * plot_rows),
        squeeze=False,
    )
    palette = [
        "#0072B2",
        "#D55E00",
        "#009E73",
        "#CC79A7",
        "#E69F00",
        "#56B4E9",
        "#999999",
        "#000000",
        "#F0E442",
    ]

    for class_index, class_name in enumerate(class_order):
        h5ad_path = (
            pseudobulk_dir
            / f"{cohort}_{class_name}_pseudobulk_counts.h5ad"
        )
        actual_hash = sha256(h5ad_path)
        if actual_hash != expected_hashes[h5ad_path.name]:
            raise RuntimeError(f"SHA-256 mismatch: {h5ad_path}")

        data = ad.read_h5ad(h5ad_path)
        matrix = data.X
        if not sp.issparse(matrix):
            matrix = sp.csr_matrix(matrix)
        else:
            matrix = matrix.tocsr()
        if matrix.shape != (manifest["n_donors"], manifest["n_genes"]):
            raise RuntimeError(f"Unexpected matrix shape for {class_name}.")
        if (
            not np.isfinite(matrix.data).all()
            or np.any(matrix.data < 0)
            or not np.array_equal(matrix.data, np.rint(matrix.data))
        ):
            raise RuntimeError(f"Non-count value in {class_name}.")

        library_size = np.asarray(matrix.sum(axis=1)).ravel().astype(np.int64)
        detected_genes = np.asarray(matrix.getnnz(axis=1), dtype=np.int64)
        observed_total_counts += int(library_size.sum(dtype=np.int64))
        recorded_total = pd.to_numeric(
            data.obs["total_counts"], errors="raise"
        ).to_numpy(dtype=np.int64)
        if not np.array_equal(library_size, recorded_total):
            raise RuntimeError(f"Row total mismatch for {class_name}.")

        n_cells = pd.to_numeric(
            data.obs["n_cells"], errors="raise"
        ).to_numpy(dtype=np.int64)
        eligible = n_cells >= args.minimum_cells
        if not np.array_equal(
            eligible,
            data.obs["eligible_min_cells"].astype(bool).to_numpy(),
        ):
            raise RuntimeError(f"Eligibility mismatch for {class_name}.")
        if np.any((library_size == 0) & eligible):
            raise RuntimeError(f"Eligible zero-library sample in {class_name}.")

        eligible_matrix = matrix[eligible]
        eligible_library = library_size[eligible].astype(float)
        if eligible.sum() < 3:
            raise RuntimeError(
                f"Too few eligible pseudobulk samples for {class_name}."
            )
        minimum_samples = max(
            10, int(math.ceil(args.minimum_fraction * eligible.sum()))
        )
        cpm = eligible_matrix.multiply(
            (1_000_000 / eligible_library)[:, None]
        ).tocsr()
        expressed_samples = np.asarray(
            (cpm >= args.cpm_threshold).sum(axis=0)
        ).ravel()
        retained = expressed_samples >= minimum_samples

        for cpm_threshold in [0.5, 1.0, 2.0]:
            expressed = cpm >= cpm_threshold
            sample_counts = np.asarray(expressed.sum(axis=0)).ravel()
            for fraction in [0.10, 0.20]:
                required = max(10, int(math.ceil(fraction * eligible.sum())))
                filter_grid_rows.append(
                    {
                        "class": class_name,
                        "cpm_threshold": cpm_threshold,
                        "minimum_fraction": fraction,
                        "minimum_samples": required,
                        "retained_genes": int(np.sum(sample_counts >= required)),
                    }
                )

        gene_filter = pd.DataFrame(
            {
                "gene_id": data.var_names.astype(str),
                "samples_cpm_at_least_threshold": expressed_samples,
                "retain_primary_filter": retained,
            }
        )
        gene_filter.to_csv(
            output_dir / f"{class_name}_gene_filter.csv", index=False
        )

        sample_qc = data.obs.copy()
        sample_qc["library_size_recomputed"] = library_size
        sample_qc["detected_genes_recomputed"] = detected_genes
        sample_qc["eligible_recomputed"] = eligible
        sample_qc["library_robust_z"] = robust_z(
            np.log10(np.maximum(library_size, 1))
        )
        sample_qc["detected_genes_robust_z"] = robust_z(detected_genes)
        sample_qc["qc_outlier_library_or_genes"] = (
            np.abs(sample_qc["library_robust_z"]) > 5
        ) | (np.abs(sample_qc["detected_genes_robust_z"]) > 5)
        sample_qc["class"] = class_name
        sample_qc_frames.append(sample_qc.reset_index(names="sample_id_index"))

        if retained.sum() < 2:
            raise RuntimeError(f"Too few retained genes for PCA: {class_name}.")
        log_cpm = np.log2(cpm[:, retained].toarray() + 0.5)
        gene_variance = np.var(log_cpm, axis=0, ddof=1)
        variable_count = min(args.pca_genes, log_cpm.shape[1])
        variable_order = np.argsort(gene_variance)[-variable_count:]
        pca_input = log_cpm[:, variable_order]
        pca_input -= pca_input.mean(axis=0, keepdims=True)
        component_count = min(10, pca_input.shape[0] - 1, pca_input.shape[1])
        pca = PCA(n_components=component_count, svd_solver="full")
        scores = pca.fit_transform(pca_input)

        eligible_obs = data.obs.loc[eligible].copy()
        pca_frame = pd.DataFrame(
            scores,
            index=eligible_obs.index,
            columns=[f"PC{index + 1}" for index in range(component_count)],
        )
        pca_frame["class"] = class_name
        pca_frame["donor_id"] = eligible_obs["donor_id"].astype(str).to_numpy()
        if args.color_field not in eligible_obs:
            raise RuntimeError(
                f"QC color field absent from pseudobulk obs: "
                f"{args.color_field}"
            )
        pca_frame["qc_color_group"] = (
            eligible_obs[args.color_field].astype(str).to_numpy()
        )
        pca_frame["age_numeric"] = numeric_age(eligible_obs["Age"]).to_numpy()
        pca_frame["sex"] = eligible_obs["Sex"].astype(str).to_numpy()
        pca_frame["ancestry"] = eligible_obs["Ancestry"].astype(str).to_numpy()
        pca_frame["n_cells"] = n_cells[eligible]
        pca_frame["library_size"] = library_size[eligible]
        pca_frame["PC1_robust_z"] = robust_z(pca_frame["PC1"].to_numpy())
        pca_frame["PC2_robust_z"] = robust_z(pca_frame["PC2"].to_numpy())
        pca_frame["qc_outlier_pca"] = (
            np.abs(pca_frame["PC1_robust_z"]) > 5
        ) | (np.abs(pca_frame["PC2_robust_z"]) > 5)
        pca_frames.append(pca_frame.reset_index(names="sample_id"))

        pca_axis = pca_axes.flat[class_index]
        plot_groups = sorted(pca_frame["qc_color_group"].unique())
        color_map = {
            group: palette[index % len(palette)]
            for index, group in enumerate(plot_groups)
        }
        for status in plot_groups:
            selected = pca_frame["qc_color_group"].eq(status)
            pca_axis.scatter(
                pca_frame.loc[selected, "PC1"],
                pca_frame.loc[selected, "PC2"],
                s=18,
                alpha=0.75,
                color=color_map[status],
                label=status,
            )
        pca_axis.set_title(class_name)
        pca_axis.set_xlabel(
            f"PC1 ({100 * pca.explained_variance_ratio_[0]:.1f}%)"
        )
        pca_axis.set_ylabel(
            f"PC2 ({100 * pca.explained_variance_ratio_[1]:.1f}%)"
        )

        library_axis = library_axes.flat[class_index]
        statuses = data.obs[args.color_field].astype(str)
        for status in plot_groups:
            selected = eligible & statuses.eq(status).to_numpy()
            library_axis.scatter(
                n_cells[selected],
                library_size[selected],
                s=18,
                alpha=0.75,
                color=color_map[status],
                label=status,
            )
        library_axis.set_xscale("log")
        library_axis.set_yscale("log")
        library_axis.set_title(class_name)
        library_axis.set_xlabel("Nuclei per donor")
        library_axis.set_ylabel("Pseudobulk library size")

        summary_rows.append(
            {
                "class": class_name,
                "donors_total": int(matrix.shape[0]),
                "donors_observed": int(np.sum(n_cells > 0)),
                "donors_eligible_10_cells": int(np.sum(n_cells >= 10)),
                "donors_eligible_20_cells": int(np.sum(n_cells >= 20)),
                "donors_eligible_50_cells": int(np.sum(n_cells >= 50)),
                "cells_total": int(n_cells.sum(dtype=np.int64)),
                "counts_total": int(library_size.sum(dtype=np.int64)),
                "median_cells_eligible": float(np.median(n_cells[eligible])),
                "median_library_eligible": float(
                    np.median(library_size[eligible])
                ),
                "median_detected_genes_eligible": float(
                    np.median(detected_genes[eligible])
                ),
                "retained_genes_primary_filter": int(retained.sum()),
                "minimum_samples_primary_filter": minimum_samples,
                "pca_genes": variable_count,
                "pca_variance_PC1": float(pca.explained_variance_ratio_[0]),
                "pca_variance_PC2": float(pca.explained_variance_ratio_[1]),
                "library_or_gene_outliers": int(
                    sample_qc["qc_outlier_library_or_genes"].sum()
                ),
                "pca_outliers": int(pca_frame["qc_outlier_pca"].sum()),
            }
        )

    for axis_index in range(len(class_order), pca_axes.size):
        pca_axes.flat[axis_index].set_visible(False)
        library_axes.flat[axis_index].set_visible(False)

    if observed_total_counts != manifest["total_counts_pseudobulk"]:
        raise RuntimeError(
            "Cross-file counts conservation failed: "
            f"{observed_total_counts} != {manifest['total_counts_pseudobulk']}"
        )

    summary = pd.DataFrame(summary_rows)
    sample_qc_all = pd.concat(sample_qc_frames, ignore_index=True)
    pca_all = pd.concat(pca_frames, ignore_index=True)
    filter_grid = pd.DataFrame(filter_grid_rows)
    summary.to_csv(output_dir / "class_qc_summary.csv", index=False)
    sample_qc_all.to_csv(output_dir / "sample_qc.csv", index=False)
    pca_all.to_csv(output_dir / "pca_scores.csv", index=False)
    filter_grid.to_csv(output_dir / "gene_filter_sensitivity.csv", index=False)

    handles, labels = pca_axes.flat[0].get_legend_handles_labels()
    pca_figure.legend(
        handles,
        labels,
        loc="upper center",
        ncol=min(6, max(1, len(labels))),
    )
    pca_figure.suptitle(
        f"{cohort} donor-level pseudobulk PCA by cell class "
        f"(color: {args.color_field})",
        y=0.99,
    )
    pca_figure.tight_layout(rect=(0, 0, 1, 0.95))
    pca_figure.savefig(output_dir / "pca_by_class.png", dpi=180)
    plt.close(pca_figure)

    handles, labels = library_axes.flat[0].get_legend_handles_labels()
    library_figure.legend(
        handles,
        labels,
        loc="upper center",
        ncol=min(6, max(1, len(labels))),
    )
    library_figure.suptitle(
        "Pseudobulk library size versus nuclei per donor", y=0.99
    )
    library_figure.tight_layout(rect=(0, 0, 1, 0.95))
    library_figure.savefig(
        output_dir / "library_size_vs_nuclei.png", dpi=180
    )
    plt.close(library_figure)

    result = {
        "status": "complete",
        "cohort": cohort,
        "color_field": args.color_field,
        "pseudobulk_manifest": str(manifest_path),
        "source_hashes_verified": True,
        "counts_conservation_verified": True,
        "total_counts": observed_total_counts,
        "minimum_cells": args.minimum_cells,
        "cpm_threshold": args.cpm_threshold,
        "minimum_fraction": args.minimum_fraction,
        "classes": summary.to_dict(orient="records"),
    }
    (output_dir / "qc_manifest.json").write_text(
        json.dumps(json_safe(result), ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    print(summary.to_string(index=False))
    print(f"QC outputs: {output_dir}")


if __name__ == "__main__":
    main()
