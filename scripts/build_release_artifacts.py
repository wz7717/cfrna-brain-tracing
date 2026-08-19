#!/usr/bin/env python
"""Build deterministic software and full-reproducibility release ZIP archives."""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
import zipfile
from pathlib import Path
from typing import Any, Callable, Iterable


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
DEFAULT_MANIFEST = ROOT / "release" / "v0.1.17" / "release_manifest.json"
LFS_PAYLOAD = "reproducibility/round5_analysis/p1_fg_permutation_fdr/permutation_fg_max.npz"
ZIP_TIMESTAMP = (1980, 1, 1, 0, 0, 0)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest().upper()


def git_output(root: Path, *args: str) -> str:
    result = subprocess.run(["git", *args], cwd=root, text=True, encoding="utf-8", errors="replace", capture_output=True, check=False)
    if result.returncode != 0:
        raise RuntimeError(f"git {' '.join(args)} failed")
    return result.stdout


def tracked_files(root: Path) -> list[str]:
    result = subprocess.run(["git", "ls-files", "-z"], cwd=root, capture_output=True, check=False)
    if result.returncode != 0:
        raise RuntimeError("git ls-files failed")
    return sorted(item.decode("utf-8") for item in result.stdout.split(b"\0") if item)


def tracked_blob_bytes(root: Path, relative: str) -> bytes:
    """Read a tracked member from Git's canonical blob, not the checkout."""
    result = subprocess.run(["git", "show", f"HEAD:{relative}"], cwd=root, capture_output=True, check=False)
    if result.returncode != 0:
        raise RuntimeError(f"could not read tracked Git blob: {relative}")
    return result.stdout


def member_payload(root: Path, relative: str, materialized_lfs_members: set[str]) -> bytes:
    path = root / relative
    if not path.is_file():
        raise FileNotFoundError(f"tracked release member is absent: {relative}")
    if relative in materialized_lfs_members:
        return path.read_bytes()
    return tracked_blob_bytes(root, relative)


def build_zip(
    root: Path,
    members: Iterable[str],
    output: Path,
    prefix: str,
    *,
    payload_reader: Callable[[str], bytes] | None = None,
) -> dict[str, Any]:
    inventory: list[dict[str, Any]] = []
    output.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(output, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=9, strict_timestamps=True) as archive:
        for relative in sorted(members):
            path = root / relative
            if not path.is_file():
                raise FileNotFoundError(f"tracked release member is absent: {relative}")
            info = zipfile.ZipInfo(f"{prefix}/{relative}", date_time=ZIP_TIMESTAMP)
            info.compress_type = zipfile.ZIP_DEFLATED
            info.external_attr = 0o100644 << 16
            payload = payload_reader(relative) if payload_reader else path.read_bytes()
            archive.writestr(info, payload, compress_type=zipfile.ZIP_DEFLATED, compresslevel=9)
            inventory.append({"path": relative, "bytes": len(payload), "sha256": hashlib.sha256(payload).hexdigest().upper()})
    return {"member_count": len(inventory), "members": inventory, "sha256": sha256_file(output), "bytes": output.stat().st_size}


def build(root: Path, manifest_path: Path, output_dir: Path, *, require_clean: bool = True) -> dict[str, Any]:
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest.get("release_state") != "final":
        raise ValueError("release artifacts may only be built after one-shot manifest finalization")
    if require_clean and git_output(root, "status", "--porcelain").strip():
        raise ValueError("release artifacts require a clean tracked worktree")
    from scripts.audit_lfs_materialization import audit as audit_lfs, lfs_paths

    lfs = audit_lfs(root)
    if lfs["status"] != "PASS":
        raise ValueError("required Git-LFS payloads are not materialized")
    materialized_lfs_members = set(lfs_paths(root))

    def payload_reader(relative: str) -> bytes:
        return member_payload(root, relative, materialized_lfs_members)

    version = str(manifest["version"])
    all_members = tracked_files(root)
    software_members = [relative for relative in all_members if relative != LFS_PAYLOAD]
    source_name = f"braintrace-v{version}-source.zip"
    full_name = f"braintrace-v{version}-full-reproducibility.zip"
    source = build_zip(root, software_members, output_dir / source_name, f"braintrace-v{version}", payload_reader=payload_reader)
    full = build_zip(root, all_members, output_dir / full_name, f"braintrace-v{version}-full-reproducibility", payload_reader=payload_reader)
    checksums = f"{source['sha256']}  {source_name}\n{full['sha256']}  {full_name}\n"
    (output_dir / "SHA256SUMS").write_text(checksums, encoding="utf-8", newline="\n")
    inventory = {
        "schema": "braintrace.release_inventory.v1",
        "version": version,
        "source": source,
        "full_reproducibility": full,
    }
    (output_dir / "RELEASE_INVENTORY.json").write_text(json.dumps(inventory, indent=2) + "\n", encoding="utf-8", newline="\n")
    integrity = {
        "schema": "braintrace.release_integrity.v1",
        "version": version,
        "git_sha": git_output(root, "rev-parse", "HEAD").strip(),
        "software_version_doi": manifest["software_version_doi"],
        "full_repro_version_doi": manifest["full_repro_version_doi"],
        "source_archive": {"name": source_name, "sha256": source["sha256"], "bytes": source["bytes"]},
        "full_repro_archive": {"name": full_name, "sha256": full["sha256"], "bytes": full["bytes"]},
        "lfs_materialization": lfs,
        "determinism": "sorted Git-tracked member order; canonical Git blob bytes for non-LFS members; materialized LFS bytes; fixed ZIP timestamp and permissions; deflate level 9",
    }
    (output_dir / f"RELEASE_INTEGRITY_v{version}.json").write_text(json.dumps(integrity, indent=2) + "\n", encoding="utf-8", newline="\n")
    return integrity


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--allow-dirty", action="store_true", help="Development-only; never use for a release archive.")
    args = parser.parse_args()
    payload = build(ROOT, args.manifest.resolve(), args.output_dir.resolve(), require_clean=not args.allow_dirty)
    print(json.dumps({"status": "PASS", "version": payload["version"]}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
