#!/usr/bin/env python
"""Create a canonical external_data layout from a separately held legacy tree.

This is a local staging helper, never a fallback used by the release runner.
It first verifies that every declared legacy source exists, then creates
hardlinks (no duplicated raw payload) at the manifest's canonical locations.
The destination must be new or empty; existing inputs are never overwritten.
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in os.sys.path:
    os.sys.path.insert(0, str(ROOT))

from core.external_inputs import DEFAULT_MANIFEST, load_manifest, sha256_file  # noqa: E402


def source_target_pairs(manifest: dict[str, Any], source_root: Path, destination_root: Path) -> list[tuple[str, Path, Path]]:
    pairs: list[tuple[str, Path, Path]] = []
    for item in manifest["inputs"]:
        if item["storage"] != "external":
            continue
        legacy_path = item.get("legacy_path")
        if not isinstance(legacy_path, str) or not legacy_path:
            raise ValueError(f"{item['alias']}: staging requires an explicit legacy_path")
        if item["kind"] == "file":
            pairs.append((str(item["alias"]), source_root / legacy_path, destination_root / item["canonical_path"]))
            continue
        for required in item["required_files"]:
            relative = Path(str(required["path"]))
            pairs.append(
                (
                    str(item["alias"]),
                    source_root / legacy_path / relative,
                    destination_root / item["canonical_path"] / relative,
                )
            )
    return pairs


def stage(source_root: Path, destination_root: Path, manifest_path: Path = DEFAULT_MANIFEST) -> dict[str, Any]:
    source_root = source_root.resolve()
    destination_root = destination_root.resolve()
    manifest = load_manifest(manifest_path)
    if destination_root.exists() and any(destination_root.iterdir()):
        raise FileExistsError(f"destination external-data root is not empty: {destination_root}")
    pairs = source_target_pairs(manifest, source_root, destination_root)
    missing = [alias for alias, source, _target in pairs if not source.is_file()]
    if missing:
        raise FileNotFoundError(f"legacy source root is incomplete for aliases: {', '.join(sorted(set(missing)))}")
    source_volumes = {source.drive.lower() for _alias, source, _target in pairs}
    if len(source_volumes | {destination_root.drive.lower()}) != 1:
        raise OSError("hardlink staging requires source and canonical external_data roots on one volume")
    destination_root.mkdir(parents=True, exist_ok=True)
    records: list[dict[str, str]] = []
    try:
        for alias, source, target in pairs:
            target.parent.mkdir(parents=True, exist_ok=True)
            os.link(source, target)
            records.append(
                {
                    "alias": alias,
                    "canonical_locator": f"external_source::{alias}/{target.name}",
                    "sha256": sha256_file(target),
                    "mode": "hardlink",
                }
            )
    except Exception:
        # This helper owns only a newly-created staging root. Removing it after
        # a failed link operation avoids leaving a partially canonical layout.
        for _alias, _source, target in reversed(pairs):
            if target.exists():
                target.unlink()
        for directory in sorted((path for path in destination_root.rglob("*") if path.is_dir()), key=lambda path: len(path.parts), reverse=True):
            try:
                directory.rmdir()
            except OSError:
                pass
        try:
            destination_root.rmdir()
        except OSError:
            pass
        raise
    payload: dict[str, Any] = {
        "schema": "braintrace.external_data_staging.v1",
        "status": "PASS",
        "canonical_root_locator": "external_data/",
        "source_layout": "legacy source tree staged as canonical hardlinks; no release fallback is in use",
        "records": records,
    }
    (destination_root / "EXTERNAL_INPUT_STAGING_AUDIT.json").write_text(
        json.dumps(payload, indent=2) + "\n", encoding="utf-8", newline="\n"
    )
    return payload


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-root", type=Path, required=True, help="Existing separately held legacy-layout source tree.")
    parser.add_argument("--external-data-root", type=Path, required=True, help="New empty canonical external_data root.")
    parser.add_argument("--input-manifest", type=Path, default=DEFAULT_MANIFEST)
    args = parser.parse_args()
    payload = stage(args.source_root, args.external_data_root, args.input_manifest)
    print(json.dumps({"status": payload["status"], "staged_records": len(payload["records"])}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
