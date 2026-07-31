from __future__ import annotations

import json
import math
import sys
import zipfile
from collections import Counter
from pathlib import Path
from xml.etree import ElementTree as ET

from inspect_xlsx_xml import NS_MAIN, NS_PKG_REL, NS_REL, cell_value, read_shared_strings


MAJOR_QTLS = {"Ast", "End", "Ext", "IN", "MG", "OD", "OPC"}
TARGET_DISEASES = {"AD", "PD", "SCZ", "DLBD"}


def sheet_rows(workbook_path: Path, sheet_name: str) -> list[list[str | None]]:
    with zipfile.ZipFile(workbook_path) as archive:
        shared_strings = read_shared_strings(archive)
        workbook = ET.fromstring(archive.read("xl/workbook.xml"))
        rels = ET.fromstring(archive.read("xl/_rels/workbook.xml.rels"))
        rel_map = {
            rel.attrib["Id"]: rel.attrib["Target"]
            for rel in rels.findall(f"{{{NS_PKG_REL}}}Relationship")
        }
        sheets = workbook.find(f"{{{NS_MAIN}}}sheets")
        if sheets is None:
            raise RuntimeError("Workbook has no sheets")
        target_path = None
        for sheet in sheets:
            if sheet.attrib["name"] == sheet_name:
                rel_id = sheet.attrib[f"{{{NS_REL}}}id"]
                target = rel_map[rel_id].lstrip("/")
                target_path = target if target.startswith("xl/") else f"xl/{target}"
                break
        if target_path is None:
            raise KeyError(f"Sheet not found: {sheet_name}")
        root = ET.fromstring(archive.read(target_path))
        sheet_data = root.find(f"{{{NS_MAIN}}}sheetData")
        rows: list[list[str | None]] = []
        if sheet_data is None:
            return rows
        for row in sheet_data:
            values: list[str | None] = []
            for cell in row.findall(f"{{{NS_MAIN}}}c"):
                ref = cell.attrib.get("r", "")
                letters = "".join(char for char in ref if char.isalpha())
                column = 0
                for char in letters:
                    column = column * 26 + ord(char.upper()) - 64
                while len(values) < column:
                    values.append(None)
                values[column - 1] = cell_value(cell, shared_strings)
            rows.append(values)
        return rows


def records(rows: list[list[str | None]], header_row: int) -> list[dict[str, str | None]]:
    header = [str(value) if value is not None else "" for value in rows[header_row]]
    result: list[dict[str, str | None]] = []
    for row in rows[header_row + 1 :]:
        if not any(value not in (None, "") for value in row):
            continue
        padded = row + [None] * (len(header) - len(row))
        result.append(dict(zip(header, padded)))
    return result


def number(value: str | None) -> float:
    if value in (None, "", "NA"):
        return math.nan
    return float(str(value).replace(",", ""))


def summarize_singlebrain(path: Path) -> dict[str, object]:
    coloc = records(sheet_rows(path, "Table11_COLOC"), 2)
    mr = records(sheet_rows(path, "Table12_MR"), 2)

    coloc_target = [row for row in coloc if row.get("Disease") in TARGET_DISEASES]
    pp4_by_key: dict[tuple[str | None, ...], float] = {}
    for row in coloc_target:
        key = (
            row.get("Disease"),
            row.get("GWAS"),
            row.get("Locus"),
            row.get("GWAS lead SNP"),
            row.get("QTL"),
            row.get("QTL lead SNP"),
            row.get("QTL gene ID"),
        )
        pp4_by_key[key] = max(pp4_by_key.get(key, -math.inf), number(row.get("PP4")))

    strict_rows: list[dict[str, str | None]] = []
    for row in mr:
        disease = row.get("Disease")
        if disease not in TARGET_DISEASES:
            continue
        key = (
            disease,
            row.get("GWAS"),
            row.get("Locus"),
            row.get("GWAS lead SNP"),
            row.get("QTL"),
            row.get("QTL lead SNP"),
            row.get("QTL gene ID"),
        )
        pp4 = pp4_by_key.get(key, math.nan)
        p_value = number(row.get("p-value"))
        egger_p = number(row.get("Egger intercept p-value"))
        nsnp = number(row.get("nsnp"))
        if (
            row.get("Method") == "Inverse variance weighted"
            and pp4 > 0.8
            and p_value < 0.05
            and nsnp >= 3
            and (math.isnan(egger_p) or egger_p > 0.05)
        ):
            enriched = dict(row)
            enriched["PP4"] = str(pp4)
            strict_rows.append(enriched)

    strict_unique: dict[tuple[str | None, ...], dict[str, str | None]] = {}
    for row in strict_rows:
        key = (
            row.get("Disease"),
            row.get("GWAS"),
            row.get("QTL"),
            row.get("QTL gene ID"),
        )
        current = strict_unique.get(key)
        if current is None or number(row.get("p-value")) < number(current.get("p-value")):
            strict_unique[key] = row

    strict = list(strict_unique.values())
    return {
        "coloc_rows_by_disease": dict(Counter(row.get("Disease") for row in coloc_target)),
        "pp4_gt_0_8_unique_gene_qtl_by_disease": dict(
            Counter(
                disease
                for disease, _, _ in {
                    (row.get("Disease"), row.get("QTL"), row.get("QTL gene ID"))
                    for row in coloc_target
                    if number(row.get("PP4")) > 0.8
                }
            )
        ),
        "strict_anchor_definition": (
            "IVW p<0.05, PP4>0.8, nsnp>=3, Egger-intercept p>0.05; "
            "descriptive screen only, not multiplicity-adjusted"
        ),
        "strict_unique_anchors_by_disease": dict(
            Counter(row.get("Disease") for row in strict)
        ),
        "strict_major_class_anchors_by_disease": dict(
            Counter(
                row.get("Disease")
                for row in strict
                if row.get("QTL") in MAJOR_QTLS
            )
        ),
        "strict_direction_by_disease": {
            disease: dict(
                Counter(
                    "positive" if number(row.get("beta")) > 0 else "negative"
                    for row in strict
                    if row.get("Disease") == disease
                )
            )
            for disease in sorted(TARGET_DISEASES)
        },
        "strict_example_rows": strict[:12],
    }


def summarize_sntwas(
    results_path: Path, bulk_path: Path
) -> dict[str, object]:
    table22 = records(sheet_rows(results_path, "TableS22"), 0)
    table23 = records(sheet_rows(results_path, "TableS23"), 0)
    table26 = records(sheet_rows(results_path, "TableS26"), 0)
    bulk = records(sheet_rows(bulk_path, "TableS21"), 0)

    target22 = [row for row in table22 if row.get("Trait") in TARGET_DISEASES]
    target23 = [row for row in table23 if row.get("Trait") in TARGET_DISEASES]
    target26 = [
        row
        for row in table26
        if row.get("Trait") in TARGET_DISEASES and number(row.get("PIP")) >= 0.5
    ]
    target_bulk = [
        row
        for row in bulk
        if row.get("Trait") in TARGET_DISEASES and number(row.get("FDR")) < 0.05
    ]

    return {
        "reported_significant_gta_counts": {
            row.get("Trait"): int(number(row.get("nSig"))) for row in target22
        },
        "novel_gene_celltype_rows_by_trait": dict(
            Counter(row.get("Trait") for row in target23)
        ),
        "focus_pip_ge_0_5_rows_by_trait": dict(
            Counter(row.get("Trait") for row in target26)
        ),
        "bulk_twas_fdr_lt_0_05_direction_by_trait": {
            disease: dict(
                Counter(
                    "positive" if number(row.get("z-score")) > 0 else "negative"
                    for row in target_bulk
                    if row.get("Trait") == disease
                )
            )
            for disease in sorted(TARGET_DISEASES)
        },
        "celltype_direction_field_available_in_downloaded_tables": False,
        "interpretation": (
            "Downloaded snTWAS tables expose cell-type gene lists and fine-mapping, "
            "but not the cell-type-specific z/beta needed as a primary direction anchor. "
            "Bulk S-TWAS z-scores are available but are not cell-type specific."
        ),
    }


def main() -> None:
    if len(sys.argv) not in (5, 6):
        raise SystemExit(
            "Usage: summarize_public_genetic_anchors.py "
            "singlebrain.xlsx sntwas_22_30.xlsx sntwas_bulk.xlsx output.json [--quiet]"
        )
    quiet = len(sys.argv) == 6 and sys.argv[5] == "--quiet"
    if len(sys.argv) == 6 and not quiet:
        raise SystemExit("The only supported optional argument is --quiet")
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")
    output = {
        "singlebrain": summarize_singlebrain(Path(sys.argv[1])),
        "psychad_sntwas": summarize_sntwas(Path(sys.argv[2]), Path(sys.argv[3])),
    }
    output_path = Path(sys.argv[4])
    output_path.write_text(
        json.dumps(output, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    if not quiet:
        print(json.dumps(output, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
