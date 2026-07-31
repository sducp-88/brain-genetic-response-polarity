#!/usr/bin/env python3
"""Download a large public S3 object as independently verified byte ranges.

The source object's official multipart ETag determines the 8 MiB part layout.
Each range is downloaded to its own temporary file and accepted only when the
HTTP status, Content-Range, and byte count all match. The assembled file is
accepted only when its locally computed multipart ETag equals the official
source ETag.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import hashlib
import json
import math
import os
import re
import subprocess
import threading
import time
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path


@dataclass(frozen=True)
class Part:
    index: int
    start: int
    stop: int
    size: int
    path: Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--url", required=True)
    parser.add_argument("--final", type=Path, required=True)
    parser.add_argument("--expected-size", type=int, required=True)
    parser.add_argument("--official-etag", required=True)
    parser.add_argument("--part-size", type=int, default=8 * 1024 * 1024)
    parser.add_argument("--workers", type=int, default=12)
    parser.add_argument("--retries", type=int, default=20)
    parser.add_argument("--status", type=Path, required=True)
    parser.add_argument("--log", type=Path, required=True)
    return parser.parse_args()


def now() -> str:
    return datetime.now(timezone.utc).astimezone().isoformat(timespec="seconds")


class Reporter:
    def __init__(
        self,
        status: Path,
        log: Path,
        expected_size: int,
        total_parts: int,
    ) -> None:
        self.status = status
        self.log = log
        self.expected_size = expected_size
        self.total_parts = total_parts
        self.lock = threading.Lock()
        self.started = time.monotonic()
        self.completed_parts = 0
        self.completed_bytes = 0
        self.seeded_parts = 0
        self.seeded_bytes = 0

    def append(self, message: str) -> None:
        with self.lock:
            self.log.parent.mkdir(parents=True, exist_ok=True)
            with self.log.open("a", encoding="utf-8") as handle:
                handle.write(f"[{now()}] {message}\n")

    def write_status(self, message: str) -> None:
        with self.lock:
            self.status.parent.mkdir(parents=True, exist_ok=True)
            temporary = self.status.with_suffix(self.status.suffix + ".tmp")
            temporary.write_text(message + "\n", encoding="utf-8")
            os.replace(temporary, self.status)

    def seed(self, parts: int, byte_count: int) -> None:
        with self.lock:
            self.completed_parts = parts
            self.completed_bytes = byte_count
            self.seeded_parts = parts
            self.seeded_bytes = byte_count

    def part_complete(self, part: Part, reused: bool = False) -> None:
        with self.lock:
            self.completed_parts += 1
            self.completed_bytes += part.size
            elapsed = max(time.monotonic() - self.started, 1e-9)
            new_bytes = self.completed_bytes - self.seeded_bytes
            rate = new_bytes / elapsed
            remaining = max(self.expected_size - self.completed_bytes, 0)
            eta_seconds = remaining / rate if rate else math.inf
            percent = 100 * self.completed_bytes / self.expected_size
            message = (
                f"RUNNING stage=parts completed={self.completed_parts}/"
                f"{self.total_parts} bytes={self.completed_bytes}/"
                f"{self.expected_size} percent={percent:.2f} "
                f"reused_parts={self.seeded_parts} "
                f"rate_MiBs={rate / (1024**2):.2f} "
                f"eta_minutes={eta_seconds / 60:.1f}"
            )
            temporary = self.status.with_suffix(self.status.suffix + ".tmp")
            temporary.write_text(message + "\n", encoding="utf-8")
            os.replace(temporary, self.status)
            if (
                reused
                or self.completed_parts % 50 == 0
                or self.completed_parts == self.total_parts
            ):
                with self.log.open("a", encoding="utf-8") as handle:
                    handle.write(f"[{now()}] {message}\n")


def expected_part_count(etag: str) -> int:
    match = re.fullmatch(r"[0-9a-fA-F]{32}-(\d+)", etag)
    if not match:
        raise ValueError(f"Not a multipart ETag: {etag!r}")
    return int(match.group(1))


def make_parts(
    expected_size: int, part_size: int, parts_dir: Path
) -> list[Part]:
    count = math.ceil(expected_size / part_size)
    parts = []
    for index in range(count):
        start = index * part_size
        stop = min(start + part_size, expected_size) - 1
        parts.append(
            Part(
                index=index,
                start=start,
                stop=stop,
                size=stop - start + 1,
                path=parts_dir / f"part_{index:05d}.bin",
            )
        )
    return parts


def valid_existing_part(part: Part) -> bool:
    return part.path.is_file() and part.path.stat().st_size == part.size


def download_part(
    part: Part,
    url: str,
    official_etag: str,
    retries: int,
    reporter: Reporter,
) -> Part:
    if valid_existing_part(part):
        return part

    if part.path.exists():
        part.path.unlink()
    temporary = part.path.with_suffix(part.path.suffix + ".tmp")
    expected_range = (
        f"bytes {part.start}-{part.stop}/{reporter.expected_size}"
    )

    for attempt in range(1, retries + 1):
        if temporary.exists():
            temporary.unlink()
        command = [
            "curl",
            "--http1.1",
            "--silent",
            "--show-error",
            "--location",
            "--fail",
            "--connect-timeout",
            "30",
            "--speed-limit",
            "1024",
            "--speed-time",
            "180",
            "--header",
            f'If-Match: "{official_etag}"',
            "--range",
            f"{part.start}-{part.stop}",
            "--output",
            str(temporary),
            "--write-out",
            "%{http_code}\t%{size_download}\t%header{content-range}",
            url,
        ]
        result = subprocess.run(
            command,
            check=False,
            capture_output=True,
            text=True,
        )
        fields = result.stdout.strip().split("\t")
        http_code = fields[0] if len(fields) >= 1 else ""
        size_download = fields[1] if len(fields) >= 2 else ""
        content_range = fields[2].strip() if len(fields) >= 3 else ""
        actual_size = temporary.stat().st_size if temporary.exists() else -1
        valid = (
            result.returncode == 0
            and http_code == "206"
            and size_download == str(part.size)
            and actual_size == part.size
            and content_range == expected_range
        )
        if valid:
            os.replace(temporary, part.path)
            return part

        reporter.append(
            f"part={part.index} attempt={attempt} failed "
            f"curl_rc={result.returncode} http={http_code!r} "
            f"size_download={size_download!r} actual_size={actual_size} "
            f"content_range={content_range!r} stderr={result.stderr.strip()!r}"
        )
        if temporary.exists():
            temporary.unlink()
        time.sleep(min(2 * attempt, 30))

    raise RuntimeError(
        f"Part {part.index} failed after {retries} attempts "
        f"(bytes {part.start}-{part.stop})."
    )


def assemble_and_verify(
    parts: list[Part],
    final: Path,
    official_etag: str,
    expected_size: int,
    reporter: Reporter,
) -> str:
    assembling = final.with_suffix(final.suffix + ".assembling")
    if assembling.exists():
        assembling.unlink()
    reporter.write_status(
        f"RUNNING stage=assemble completed=0/{len(parts)} bytes=0/"
        f"{expected_size}"
    )
    digests: list[bytes] = []
    written = 0
    with assembling.open("wb") as output:
        for position, part in enumerate(parts, start=1):
            if not valid_existing_part(part):
                raise RuntimeError(f"Missing or invalid part before assembly: {part}")
            payload = part.path.read_bytes()
            digests.append(hashlib.md5(payload).digest())
            output.write(payload)
            written += len(payload)
            if position % 100 == 0 or position == len(parts):
                reporter.write_status(
                    f"RUNNING stage=assemble completed={position}/"
                    f"{len(parts)} bytes={written}/{expected_size}"
                )
        output.flush()
        os.fsync(output.fileno())

    if written != expected_size or assembling.stat().st_size != expected_size:
        raise RuntimeError(
            f"Assembled size mismatch: wrote={written}, "
            f"stat={assembling.stat().st_size}, expected={expected_size}."
        )
    local_etag = (
        f"{hashlib.md5(b''.join(digests)).hexdigest()}-{len(digests)}"
    )
    if local_etag != official_etag:
        mismatch = final.with_suffix(final.suffix + ".etag_mismatch")
        if mismatch.exists():
            mismatch = final.with_suffix(
                final.suffix + f".etag_mismatch.{int(time.time())}"
            )
        os.replace(assembling, mismatch)
        raise RuntimeError(
            f"Multipart ETag mismatch: local={local_etag}, "
            f"official={official_etag}; assembled copy retained at {mismatch}."
        )
    os.replace(assembling, final)
    return local_etag


def main() -> None:
    args = parse_args()
    if args.workers < 1 or args.retries < 1 or args.part_size < 1:
        raise ValueError("workers, retries, and part size must be positive.")
    args.final = args.final.resolve()
    args.status = args.status.resolve()
    args.log = args.log.resolve()
    args.final.parent.mkdir(parents=True, exist_ok=True)
    parts_dir = args.final.with_suffix(args.final.suffix + ".parts")
    parts_dir.mkdir(parents=True, exist_ok=True)
    parts = make_parts(args.expected_size, args.part_size, parts_dir)
    etag_parts = expected_part_count(args.official_etag)
    if len(parts) != etag_parts:
        raise RuntimeError(
            f"Part layout mismatch: computed={len(parts)}, ETag={etag_parts}."
        )
    reporter = Reporter(
        args.status, args.log, args.expected_size, len(parts)
    )
    reporter.append(
        f"START url={args.url} final={args.final} "
        f"size={args.expected_size} etag={args.official_etag} "
        f"part_size={args.part_size} parts={len(parts)} workers={args.workers}"
    )

    if args.final.exists():
        reporter.write_status(
            f"FAILED refusing_to_overwrite_existing_final={args.final}"
        )
        raise FileExistsError(args.final)

    reusable = [part for part in parts if valid_existing_part(part)]
    reusable_bytes = sum(part.size for part in reusable)
    reporter.seed(len(reusable), reusable_bytes)
    if reusable:
        reporter.append(
            f"REUSE valid_size_parts={len(reusable)} bytes={reusable_bytes}"
        )

    pending = [part for part in parts if not valid_existing_part(part)]
    failures: list[str] = []
    with concurrent.futures.ThreadPoolExecutor(
        max_workers=args.workers
    ) as executor:
        futures = {
            executor.submit(
                download_part,
                part,
                args.url,
                args.official_etag,
                args.retries,
                reporter,
            ): part
            for part in pending
        }
        for future in concurrent.futures.as_completed(futures):
            part = futures[future]
            try:
                completed = future.result()
            except Exception as error:
                failures.append(f"part={part.index}: {error!r}")
                reporter.append(failures[-1])
            else:
                reporter.part_complete(completed)

    if failures:
        reporter.write_status(
            f"FAILED stage=parts failed={len(failures)} "
            f"completed={reporter.completed_parts}/{len(parts)}"
        )
        raise RuntimeError("; ".join(failures[:10]))

    reporter.write_status(
        f"RUNNING stage=assemble completed_parts={len(parts)}/{len(parts)}"
    )
    local_etag = assemble_and_verify(
        parts,
        args.final,
        args.official_etag,
        args.expected_size,
        reporter,
    )
    provenance = {
        "status": "SUCCESS",
        "url": args.url,
        "final": str(args.final),
        "file_size": args.final.stat().st_size,
        "part_size": args.part_size,
        "part_count": len(parts),
        "official_etag": args.official_etag,
        "local_etag": local_etag,
        "completed_at": now(),
    }
    provenance_path = args.final.with_suffix(
        args.final.suffix + ".provenance.json"
    )
    provenance_path.write_text(
        json.dumps(provenance, indent=2), encoding="utf-8"
    )
    reporter.append(
        f"SUCCESS final={args.final} size={args.final.stat().st_size} "
        f"etag={local_etag}"
    )
    reporter.write_status(
        f"SUCCESS final={args.final} bytes={args.final.stat().st_size} "
        f"etag={local_etag}"
    )


if __name__ == "__main__":
    main()
