#!/usr/bin/env python
"""Fail closed when a required Git-LFS payload is still a pointer file."""

from __future__ import annotations

import argparse
import json
import subprocess
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
POINTER_PREFIX = b"version https://git-lfs.github.com/spec/v1\n"


def is_lfs_pointer(path: Path) -> bool:
    try:
        with path.open("rb") as handle:
            return handle.read(len(POINTER_PREFIX)) == POINTER_PREFIX
    except OSError:
        return False


def lfs_paths(root: Path) -> list[str]:
    result = subprocess.run(
        ["git", "lfs", "ls-files", "--all", "-n"],
        cwd=root,
        text=True,
        encoding="utf-8",
        errors="replace",
        capture_output=True,
        check=False,
    )
    if result.returncode != 0:
        raise RuntimeError("git lfs ls-files failed")
    return sorted({line.strip().replace("\\", "/") for line in result.stdout.splitlines() if line.strip()})


def audit(root: Path) -> dict[str, Any]:
    paths = lfs_paths(root)
    missing = [relative for relative in paths if not (root / relative).is_file()]
    pointers = [relative for relative in paths if (root / relative).is_file() and is_lfs_pointer(root / relative)]
    errors: list[str] = []
    if missing:
        errors.append("missing required Git-LFS payloads")
    if pointers:
        errors.append("required Git-LFS payloads remain pointer files")
    return {
        "schema": "braintrace.lfs_materialization_audit.v1",
        "status": "PASS" if not errors else "FAIL",
        "required_lfs_files": len(paths),
        "missing_required_files": missing,
        "required_pointer_files": len(pointers),
        "pointer_files": pointers,
        "errors": errors,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, default=ROOT)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    payload = audit(args.repo_root.resolve())
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8", newline="\n")
    print(json.dumps({"status": payload["status"], "required_pointer_files": payload["required_pointer_files"]}))
    return 0 if payload["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
