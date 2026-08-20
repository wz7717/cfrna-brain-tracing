#!/usr/bin/env python3
"""Generate or verify complete root package-integrity manifests."""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.audit_public_release_content import public_files  # noqa: E402


PACKAGE_MANIFEST = "PACKAGE_MANIFEST.csv"
CHECKSUM_MANIFEST = "SHA256SUMS.txt"
SELF_MANIFESTS = {PACKAGE_MANIFEST, CHECKSUM_MANIFEST}
LFS_POINTER_PREFIX = b"version https://git-lfs.github.com/spec/v1\n"


def _git(root: Path, *args: str, input_bytes: bytes | None = None) -> subprocess.CompletedProcess[bytes]:
    return subprocess.run(
        ["git", *args],
        cwd=root,
        input=input_bytes,
        capture_output=True,
        check=False,
    )


def _index_payloads(root: Path) -> dict[str, bytes] | None:
    """Return canonical stage-0 Git blobs, or ``None`` outside a checkout."""

    listing = _git(root, "ls-files", "--stage", "-z")
    if listing.returncode != 0:
        return None
    oid_by_path: dict[str, str] = {}
    for record in listing.stdout.split(b"\0"):
        if not record:
            continue
        metadata, raw_path = record.split(b"\t", 1)
        _mode, raw_oid, stage = metadata.split()
        if stage == b"0":
            oid_by_path[raw_path.decode("utf-8")] = raw_oid.decode("ascii")

    unique_oids = sorted(set(oid_by_path.values()))
    batch = _git(root, "cat-file", "--batch", input_bytes=("\n".join(unique_oids) + "\n").encode("ascii"))
    if batch.returncode != 0:
        raise RuntimeError(batch.stderr.decode("utf-8", errors="replace"))
    payload_by_oid: dict[str, bytes] = {}
    cursor = 0
    for requested_oid in unique_oids:
        line_end = batch.stdout.index(b"\n", cursor)
        header = batch.stdout[cursor:line_end].decode("ascii").split()
        if len(header) != 3 or header[1] != "blob":
            raise RuntimeError(f"unexpected git cat-file response for {requested_oid}: {header}")
        size = int(header[2])
        start = line_end + 1
        end = start + size
        payload_by_oid[requested_oid] = batch.stdout[start:end]
        if batch.stdout[end : end + 1] != b"\n":
            raise RuntimeError(f"malformed git cat-file payload boundary for {requested_oid}")
        cursor = end + 1
    return {path: payload_by_oid[oid] for path, oid in oid_by_path.items()}


def _changed_paths(root: Path) -> set[str]:
    changed: set[str] = set()
    for args in (
        ("diff", "--name-only", "-z"),
        ("diff", "--cached", "--name-only", "-z"),
        ("ls-files", "--others", "--exclude-standard", "-z"),
    ):
        result = _git(root, *args)
        if result.returncode != 0:
            return set()
        changed.update(item.decode("utf-8") for item in result.stdout.split(b"\0") if item)
    return changed


def _text_attribute(root: Path, relative: str) -> str | None:
    result = _git(root, "check-attr", "-z", "text", "--", relative)
    if result.returncode != 0:
        return None
    fields = result.stdout.rstrip(b"\0").split(b"\0")
    return fields[2].decode("utf-8") if len(fields) == 3 else None


def _looks_text(payload: bytes) -> bool:
    if b"\0" in payload:
        return False
    try:
        payload.decode("utf-8")
    except UnicodeDecodeError:
        return False
    return True


def _canonical_worktree_payload(root: Path, relative: str, payload: bytes) -> bytes:
    """Apply the repository's LF clean rule without mutating the checkout."""

    attribute = _text_attribute(root, relative)
    if attribute == "unset":
        return payload
    if attribute == "set" or (attribute in {"auto", "unspecified", None} and _looks_text(payload)):
        return payload.replace(b"\r\n", b"\n")
    return payload


def canonical_payloads(root: Path, members: list[str]) -> dict[str, bytes]:
    """Resolve package bytes independently of checkout EOL conversion.

    Clean tracked files use their canonical Git blob. Modified or untracked
    text is passed through the repository LF clean rule. Git-LFS members use
    the materialized working-tree payload rather than the pointer blob.
    """

    index = _index_payloads(root)
    changed = _changed_paths(root) if index is not None else set(members)
    resolved: dict[str, bytes] = {}
    for relative in members:
        working = root / relative
        indexed = index.get(relative) if index is not None else None
        if indexed is not None and indexed.startswith(LFS_POINTER_PREFIX):
            resolved[relative] = working.read_bytes()
        elif indexed is not None and relative not in changed:
            resolved[relative] = indexed
        else:
            resolved[relative] = _canonical_worktree_payload(root, relative, working.read_bytes())
    return resolved


def render(root: Path = ROOT) -> tuple[str, str]:
    members = [path for path in public_files(root) if path not in SELF_MANIFESTS]
    payloads = canonical_payloads(root, members)
    csv_buffer = io.StringIO(newline="")
    writer = csv.writer(csv_buffer, lineterminator="\n")
    writer.writerow(["relative_path", "bytes"])
    checksum_lines: list[str] = []
    for relative in members:
        payload = payloads[relative]
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
