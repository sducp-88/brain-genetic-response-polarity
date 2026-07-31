"""Audit and inventory the privacy-bounded public reproducibility release."""

from __future__ import annotations

import csv
import hashlib
import json
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "BUNDLE_FILE_MANIFEST.csv"
AUDIT = ROOT / "PUBLIC_RELEASE_BOUNDARY_AUDIT.json"
GENERATED = {MANIFEST.name, AUDIT.name}

FORBIDDEN_EXTENSIONS = {
    ".doc",
    ".docx",
    ".h5",
    ".hdf5",
    ".h5ad",
    ".loom",
    ".parquet",
    ".rds",
    ".xls",
    ".xlsx",
}
FORBIDDEN_NAME_PATTERNS = (
    "manuscript_draft",
    "cover_letter",
)
FORBIDDEN_TEXT_PATTERNS = {
    "windows_user_path": re.compile(r"(?i)[A-Z]:\\Users\\"),
    "linux_home_path": re.compile(r"(?<![A-Za-z0-9_])/home/[A-Za-z0-9_.-]+/"),
    "github_pat": re.compile(r"\b(?:ghp|github_pat)_[A-Za-z0-9_]{20,}\b"),
    "synapse_token": re.compile(r"(?i)\b(?:synapse|auth)[_-]?token\s*[:=]\s*\S+"),
    "private_key": re.compile(r"-----BEGIN (?:RSA |EC |OPENSSH )?PRIVATE KEY-----"),
}
EMAIL_RE = re.compile(r"\b[A-Z0-9._%+-]+@[A-Z0-9.-]+\.[A-Z]{2,}\b", re.I)
ALLOWED_EMAILS = {"sducp@email.sdu.edu.cn"}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def text_content(path: Path) -> str | None:
    try:
        return path.read_text(encoding="utf-8")
    except (UnicodeDecodeError, OSError):
        return None


def main() -> None:
    files = sorted(
        path
        for path in ROOT.rglob("*")
        if path.is_file()
        and ".git" not in path.parts
        and path.name not in GENERATED
    )

    rows = []
    violations = []
    found_emails = set()
    for path in files:
        relative = path.relative_to(ROOT).as_posix()
        rows.append(
            {
                "path": relative,
                "size_bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
        )

        suffix = path.suffix.lower()
        if suffix in FORBIDDEN_EXTENSIONS:
            violations.append(
                {"path": relative, "rule": "forbidden_extension", "match": suffix}
            )
        lowered = relative.lower()
        for pattern in FORBIDDEN_NAME_PATTERNS:
            if pattern in lowered:
                violations.append(
                    {"path": relative, "rule": "forbidden_filename", "match": pattern}
                )

        text = text_content(path)
        if text is None:
            continue
        emails = {match.lower() for match in EMAIL_RE.findall(text)}
        found_emails.update(emails)
        for email in sorted(emails - ALLOWED_EMAILS):
            violations.append(
                {"path": relative, "rule": "noncorresponding_email", "match": email}
            )
        for name, pattern in FORBIDDEN_TEXT_PATTERNS.items():
            if pattern.search(text):
                violations.append(
                    {"path": relative, "rule": name, "match": "redacted-in-report"}
                )

    with MANIFEST.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["path", "size_bytes", "sha256"],
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)

    report = {
        "release_root": ".",
        "scope": (
            "Code, environment, non-identifying aggregate inputs/results, figures, "
            "provenance, author names/affiliations/ORCIDs, and corresponding-author "
            "institutional email. The journal manuscript, street address, and "
            "non-corresponding-author emails are excluded."
        ),
        "file_count_excluding_generated_audit_files": len(rows),
        "manifest_rows": len(rows),
        "allowed_emails_found": sorted(found_emails & ALLOWED_EMAILS),
        "unexpected_emails_found": sorted(found_emails - ALLOWED_EMAILS),
        "violations": violations,
        "status": "PASS" if not violations else "FAIL",
    }
    with AUDIT.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    if violations:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
