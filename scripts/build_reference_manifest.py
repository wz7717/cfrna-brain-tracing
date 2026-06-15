#!/usr/bin/env python
"""Rebuild the combined reference atlas manifest from per-source manifests and files."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
try:
    import yaml
except ModuleNotFoundError:
    yaml = None


ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "data" / "reference_atlases"
MANIFEST_DIR = BASE / "00_manifest"
FIELDS = [
    "atlas_name",
    "publication",
    "doi",
    "source_url",
    "accession",
    "file_name",
    "file_type",
    "expected_format",
    "download_status",
    "md5_or_size_if_available",
    "notes",
    "recommended_use",
]


def read_source_manifests() -> pd.DataFrame:
    frames = []
    for path in sorted(MANIFEST_DIR.glob("*_manifest.tsv")):
        if path.name == "reference_atlas_manifest.tsv":
            continue
        frames.append(pd.read_csv(path, sep="\t", dtype=str).fillna(""))
    if not frames:
        return pd.DataFrame(columns=FIELDS)
    table = pd.concat(frames, ignore_index=True)
    for field in FIELDS:
        if field not in table.columns:
            table[field] = ""
    return table[FIELDS].drop_duplicates(subset=["atlas_name", "source_url", "file_name"], keep="last")


def add_unmanifested_files(table: pd.DataFrame) -> pd.DataFrame:
    known = set(table["file_name"].astype(str)) if not table.empty else set()
    rows = []
    for path in BASE.rglob("*"):
        if not path.is_file() or MANIFEST_DIR in path.parents or path.name.endswith(".part"):
            continue
        if path.name in known:
            continue
        rows.append(
            {
                "atlas_name": "unmanifested_local_file",
                "publication": "",
                "doi": "",
                "source_url": "",
                "accession": "",
                "file_name": path.name,
                "file_type": "local_file",
                "expected_format": path.suffix.lstrip("."),
                "download_status": "local_file_found",
                "md5_or_size_if_available": f"size={path.stat().st_size}",
                "notes": str(path.relative_to(BASE)),
                "recommended_use": "",
            }
        )
    if rows:
        table = pd.concat([table, pd.DataFrame(rows, columns=FIELDS)], ignore_index=True)
    return table[FIELDS]


def main() -> int:
    MANIFEST_DIR.mkdir(parents=True, exist_ok=True)
    table = add_unmanifested_files(read_source_manifests())
    tsv = MANIFEST_DIR / "reference_atlas_manifest.tsv"
    js = MANIFEST_DIR / "reference_atlas_manifest.json"
    yml = MANIFEST_DIR / "reference_atlas_manifest.yaml"
    table.to_csv(tsv, sep="\t", index=False)
    js.write_text(json.dumps(table.to_dict(orient="records"), indent=2, ensure_ascii=False), encoding="utf-8")
    records = table.to_dict(orient="records")
    with yml.open("w", encoding="utf-8") as handle:
        if yaml is not None:
            yaml.safe_dump(records, handle, allow_unicode=True, sort_keys=False)
        else:
            handle.write(json.dumps(records, indent=2, ensure_ascii=False))
    print(f"Wrote {tsv} ({len(table)} rows)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
