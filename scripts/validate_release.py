"""Validate public-release syntax, metadata, and file-manifest hashes."""

from __future__ import annotations

import ast
import csv
import hashlib
import json
from pathlib import Path

import yaml


ROOT = Path(__file__).resolve().parents[1]


def main() -> None:
    python_files = sorted(ROOT.rglob("*.py"))
    for path in python_files:
        ast.parse(path.read_text(encoding="utf-8"), filename=str(path))

    yaml.safe_load((ROOT / "CITATION.cff").read_text(encoding="utf-8"))
    json.loads((ROOT / ".zenodo.json").read_text(encoding="utf-8"))

    with (ROOT / "BUNDLE_FILE_MANIFEST.csv").open(encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    mismatches = []
    for row in rows:
        path = ROOT / row["path"]
        if not path.is_file():
            mismatches.append({"path": row["path"], "reason": "missing"})
            continue
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        if digest != row["sha256"]:
            mismatches.append({"path": row["path"], "reason": "sha256"})

    result = {
        "python_files_parsed": len(python_files),
        "manifest_rows": len(rows),
        "manifest_mismatches": mismatches,
        "citation_cff": "valid_yaml",
        "zenodo_json": "valid_json",
        "status": "PASS" if not mismatches else "FAIL",
    }
    print(json.dumps(result, indent=2))
    if mismatches:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
