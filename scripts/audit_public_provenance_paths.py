#!/usr/bin/env python
"""Fail if public provenance artifacts expose a private machine-local path."""

from __future__ import annotations

import argparse
import json
import re
import subprocess
from pathlib import Path
from typing import Iterable


ROOT = Path(__file__).resolve().parents[1]
PRIVATE_PATTERNS = {
    "windows_drive": re.compile(r"[A-Za-z]:\\"),
    "windows_users": re.compile(r"C:\\Users\\", re.IGNORECASE),
    "macos_users": re.compile(r"/Users/"),
    "unix_home_user": re.compile(r"/home/(?!braintrace(?:/|$))"),
}
PUBLIC_SUFFIXES = {".md", ".json", ".csv", ".txt", ".yml", ".yaml", ".cff"}
PUBLIC_ROOT_FILES = {
    "README.md",
    "DATA_PROVENANCE.md",
    "SCIENTIFIC_REMEDIATION_QA.json",
    "SCIENTIFIC_REMEDIATION_REPORT.md",
    "NONHUANG_SCIENTIFIC_CONFLICT_LEDGER.csv",
    "scientific_claim_ledger.csv",
}


def _fallback_paths() -> Iterable[Path]:
    """Discover the public tree in a source image that deliberately omits .git."""

    for name in PUBLIC_ROOT_FILES:
        yield ROOT / name
    for directory in ("reproducibility", "release", "data", "bo2023_bulk_atlas_buildkit"):
        root = ROOT / directory
        if root.is_dir():
            yield from (path for path in root.rglob("*") if path.is_file())


def tracked_paths() -> Iterable[Path]:
    try:
        result = subprocess.run(["git", "ls-files"], cwd=ROOT, text=True, capture_output=True, check=True)
        candidates = (ROOT / line for line in result.stdout.splitlines() if line.strip())
    except (OSError, subprocess.CalledProcessError):
        candidates = _fallback_paths()
    seen: set[Path] = set()
    for candidate in candidates:
        try:
            path = candidate.relative_to(ROOT)
        except ValueError:
            continue
        absolute = ROOT / path
        if not absolute.is_file():
            continue
        if absolute in seen:
            continue
        seen.add(absolute)
        if path.suffix.lower() not in PUBLIC_SUFFIXES:
            continue
        if path.name in {"Dockerfile", "Dockerfile.repro"}:
            continue
        in_public_tree = path.parts and path.parts[0] in {"reproducibility", "release", "data", "bo2023_bulk_atlas_buildkit"}
        if path.as_posix() in PUBLIC_ROOT_FILES or in_public_tree:
            # Source code is not a generated/public provenance artifact.
            if path.suffix.lower() == ".py":
                continue
            yield absolute


def audit() -> dict[str, object]:
    matches: list[dict[str, object]] = []
    scanned: list[str] = []
    for path in tracked_paths():
        if not path.is_file():
            continue
        try:
            text = path.read_text(encoding="utf-8")
        except UnicodeDecodeError:
            continue
        relative = path.relative_to(ROOT).as_posix()
        scanned.append(relative)
        for label, pattern in PRIVATE_PATTERNS.items():
            for match in pattern.finditer(text):
                line = text.count("\n", 0, match.start()) + 1
                matches.append({"path": relative, "line": line, "pattern": label})
    return {
        "schema": "braintrace.public_provenance_path_audit.v1",
        "status": "PASS" if not matches else "FAIL",
        "CURRENT_PUBLIC_ABSOLUTE_LOCAL_PATH_MATCHES": len(matches),
        "scanned_files": len(scanned),
        "matches": matches,
        "scope": "tracked public provenance/document artifacts; excludes executable source and container-internal paths",
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    payload = audit()
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8", newline="\n")
    print(json.dumps(payload, ensure_ascii=False))
    return 0 if payload["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
