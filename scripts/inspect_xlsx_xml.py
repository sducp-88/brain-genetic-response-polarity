from __future__ import annotations

import json
import re
import sys
import zipfile
from pathlib import Path
from xml.etree import ElementTree as ET


NS_MAIN = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"
NS_REL = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
NS_PKG_REL = "http://schemas.openxmlformats.org/package/2006/relationships"


def column_number(cell_ref: str) -> int:
    letters = re.match(r"[A-Z]+", cell_ref)
    value = 0
    for char in letters.group(0) if letters else "":
        value = value * 26 + ord(char) - 64
    return value


def read_shared_strings(archive: zipfile.ZipFile) -> list[str]:
    path = "xl/sharedStrings.xml"
    if path not in archive.namelist():
        return []
    root = ET.fromstring(archive.read(path))
    strings: list[str] = []
    for item in root.findall(f"{{{NS_MAIN}}}si"):
        strings.append("".join(node.text or "" for node in item.iter(f"{{{NS_MAIN}}}t")))
    return strings


def cell_value(cell: ET.Element, shared_strings: list[str]) -> str | None:
    cell_type = cell.attrib.get("t")
    if cell_type == "inlineStr":
        return "".join(node.text or "" for node in cell.iter(f"{{{NS_MAIN}}}t"))
    value_node = cell.find(f"{{{NS_MAIN}}}v")
    if value_node is None:
        return None
    raw = value_node.text or ""
    if cell_type == "s":
        return shared_strings[int(raw)]
    if cell_type == "b":
        return "TRUE" if raw == "1" else "FALSE"
    return raw


def main() -> None:
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")
    if len(sys.argv) != 2:
        raise SystemExit("Usage: inspect_xlsx_xml.py input.xlsx")
    workbook_path = Path(sys.argv[1])
    result: list[dict[str, object]] = []
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
        for sheet in sheets:
            name = sheet.attrib["name"]
            rel_id = sheet.attrib[f"{{{NS_REL}}}id"]
            target = rel_map[rel_id].lstrip("/")
            xml_path = target if target.startswith("xl/") else f"xl/{target}"
            root = ET.fromstring(archive.read(xml_path))
            dimension = root.find(f"{{{NS_MAIN}}}dimension")
            rows_out: list[list[str | None]] = []
            max_row = 0
            max_column_seen = 0
            sheet_data = root.find(f"{{{NS_MAIN}}}sheetData")
            if sheet_data is not None:
                for row in list(sheet_data):
                    max_row = max(max_row, int(row.attrib.get("r", "0")))
                    values: dict[int, str | None] = {}
                    for cell in row.findall(f"{{{NS_MAIN}}}c"):
                        values[column_number(cell.attrib.get("r", ""))] = cell_value(
                            cell, shared_strings
                        )
                    max_column_seen = max(max_column_seen, max(values, default=0))
                    if any(value not in (None, "") for value in values.values()):
                        row_limit = 30 if name == "Dictionary" else 5
                        if len(rows_out) < row_limit:
                            max_column = min(max(values, default=0), 20)
                            rows_out.append(
                                [values.get(index) for index in range(1, max_column + 1)]
                            )
            result.append(
                {
                    "sheet": name,
                    "dimension": dimension.attrib.get("ref") if dimension is not None else None,
                    "max_row": max_row,
                    "max_column": max_column_seen,
                    "first_nonempty_rows": rows_out,
                }
            )
    print(json.dumps(result, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
