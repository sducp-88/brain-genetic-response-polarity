#!/usr/bin/env python3
"""Compute the Amazon S3 multipart ETag for a local file.

For the CELLxGENE MSSM object, the official ETag suffix reports 4303 parts.
The object length implies an 8 MiB part size:
ceil(36,092,176,654 / 8,388,608) == 4303.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import time
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("path", type=Path)
    parser.add_argument("--part-size", type=int, default=8 * 1024 * 1024)
    parser.add_argument("--progress-every", type=int, default=250)
    args = parser.parse_args()

    size = args.path.stat().st_size
    part_digests: list[bytes] = []
    started = time.time()
    bytes_read = 0

    with args.path.open("rb") as handle:
        part_index = 0
        while True:
            block = handle.read(args.part_size)
            if not block:
                break
            part_index += 1
            bytes_read += len(block)
            part_digests.append(hashlib.md5(block).digest())
            if part_index % args.progress_every == 0:
                elapsed = max(time.time() - started, 1e-9)
                gib = bytes_read / (1024**3)
                rate = bytes_read / elapsed / (1024**2)
                print(
                    f"parts={part_index} read={gib:.2f} GiB "
                    f"rate={rate:.1f} MiB/s",
                    flush=True,
                )

    if not part_digests:
        etag = hashlib.md5(b"").hexdigest()
    elif len(part_digests) == 1:
        etag = part_digests[0].hex()
    else:
        etag = f"{hashlib.md5(b''.join(part_digests)).hexdigest()}-{len(part_digests)}"

    result = {
        "path": str(args.path),
        "file_size": size,
        "part_size": args.part_size,
        "part_count": len(part_digests),
        "etag": etag,
        "elapsed_seconds": round(time.time() - started, 3),
        "basename": os.path.basename(args.path),
    }
    print(json.dumps(result, indent=2), flush=True)


if __name__ == "__main__":
    main()
