#!/usr/bin/env python3
"""Generate or verify complete root package-integrity manifests."""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.audit_public_release_content import public_files  # noqa: E402


PACKAGE_MANIFEST = "PACKAGE_MANIFEST.csv"
CHECKSUM_MANIFEST = "SHA256SUMS.txt"
SELF_MANIFESTS = {PACKAGE_MANIFEST, CHECKSUM_MANIFEST}


def render(root: Path = ROOT) -> tuple[str, str]:
    members = [path for path in public_files(root) if path not in SELF_MANIFESTS]
    csv_buffer = io.StringIO(newline="")
    writer = csv.writer(csv_buffer, lineterminator="\n")
    writer.writerow(["relative_path", "bytes"])
    checksum_lines: list[str] = []
    for relative in members:
        payload = (root / relative).read_bytes()
        writer.writerow([relative, len(payload)])
        checksum_lines.append(f"{hashlib.sha256(payload).hexdigest()}  {relative}")
    return csv_buffer.getvalue(), "\n".join(checksum_lines) + "\n"


def verify(root: Path = ROOT) -> tuple[bool, list[str]]:
    expected_csv, expected_sha = render(root)
    errors: list[str] = []
    manifest_path = root / PACKAGE_MANIFEST
    checksum_path = root / CHECKSUM_MANIFEST
    if not manifest_path.is_file() or manifest_path.read_text(encoding="utf-8") != expected_csv:
        errors.append(PACKAGE_MANIFEST)
    if not checksum_path.is_file() or checksum_path.read_text(encoding="utf-8") != expected_sha:
        errors.append(CHECKSUM_MANIFEST)
    return not errors, errors


def write(root: Path = ROOT) -> None:
    manifest, checksums = render(root)
    (root / PACKAGE_MANIFEST).write_text(manifest, encoding="utf-8", newline="\n")
    (root / CHECKSUM_MANIFEST).write_text(checksums, encoding="utf-8", newline="\n")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=ROOT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    root = args.root.resolve()
    if args.check:
        passed, errors = verify(root)
        print(f"PACKAGE_INTEGRITY={'PASS' if passed else 'FAIL'}")
        if errors:
            print("STALE=" + ",".join(errors))
        return 0 if passed else 1
    write(root)
    print(f"WROTE={root / PACKAGE_MANIFEST}")
    print(f"WROTE={root / CHECKSUM_MANIFEST}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
