#!/usr/bin/env python
from __future__ import annotations

import argparse
import json
import re
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MAP = ROOT / "reports" / "bo2023_publication_label_curation_map_20260704.csv"
DEFAULT_SUMMARY = ROOT / "reports" / "bo2023_publication_result_label_remap_audit_20260704.json"

REGION_COLUMN_HINTS = (
    "region",
    "label",
    "pred",
    "true",
    "truth",
    "group",
    "member",
    "network",
    "anatom",
)
NON_REGION_COLUMN_HINTS = (
    "gene",
    "symbol",
    "ensembl",
    "ensmfag",
    "path",
    "file",
    "evidence",
    "status",
    "query",
    "signature",
    "overlap",
)


def should_remap_column(column: object) -> bool:
    name = str(column).strip().lower()
    if not name:
        return False
    if any(hint in name for hint in NON_REGION_COLUMN_HINTS):
        return False
    return any(hint in name for hint in REGION_COLUMN_HINTS)


def remap_scalar(value: object, mapping: dict[str, str]) -> object:
    if pd.isna(value):
        return value
    text = str(value)
    stripped = text.strip()
    if stripped in mapping:
        return mapping[stripped]
    if "|" in text:
        parts = [part.strip() for part in text.split("|")]
        remapped = [mapping.get(part, part) for part in parts]
        if remapped != parts:
            return " | ".join(remapped)
    if ";" in text:
        parts = [part.strip() for part in text.split(";")]
        remapped = [mapping.get(part, part) for part in parts]
        if remapped != parts:
            return ";".join(remapped)
    remapped_text = text
    for old, new in mapping.items():
        remapped_text = re.sub(
            rf"(?<![A-Za-z0-9_]){re.escape(old)}(?![A-Za-z0-9_])",
            new,
            remapped_text,
        )
    return remapped_text if remapped_text != text else value


def remap_csv(path: Path, mapping: dict[str, str], dry_run: bool) -> dict[str, object]:
    try:
        df = pd.read_csv(path, dtype=str, keep_default_na=False)
    except Exception as exc:
        return {"path": str(path), "status": "read_error", "error": str(exc)}

    old_columns = list(df.columns)
    changed_cells = 0
    changed_columns = []

    new_columns = [remap_scalar(col, mapping) for col in df.columns]
    if new_columns != old_columns:
        df.columns = new_columns
        changed_columns.extend(
            f"header:{old}->{new}" for old, new in zip(old_columns, new_columns) if old != new
        )

    for col in list(df.columns):
        if not should_remap_column(col):
            continue
        before = df[col].copy()
        df[col] = df[col].map(lambda value: remap_scalar(value, mapping))
        changed = int((before != df[col]).sum())
        if changed:
            changed_cells += changed
            changed_columns.append(str(col))

    if changed_cells or changed_columns:
        if not dry_run:
            df.to_csv(path, index=False)
        return {
            "path": str(path),
            "status": "changed" if not dry_run else "would_change",
            "changed_cells": changed_cells,
            "changed_columns": sorted(set(changed_columns)),
        }
    return {"path": str(path), "status": "unchanged", "changed_cells": 0, "changed_columns": []}


def main() -> int:
    parser = argparse.ArgumentParser(description="Remap publication-facing Bo2023 region labels in result CSVs.")
    parser.add_argument("paths", nargs="+", type=Path, help="CSV files or directories to scan recursively.")
    parser.add_argument("--map", type=Path, default=DEFAULT_MAP)
    parser.add_argument("--summary", type=Path, default=DEFAULT_SUMMARY)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    curation = pd.read_csv(args.map, dtype=str).fillna("")
    mapping = {
        str(row.old_region_id).strip(): str(row.new_region_id).strip()
        for row in curation.itertuples(index=False)
        if str(row.old_region_id).strip() and str(row.old_region_id).strip() != str(row.new_region_id).strip()
    }

    csv_paths: list[Path] = []
    for path in args.paths:
        if path.is_dir():
            csv_paths.extend(sorted(path.rglob("*.csv")))
        elif path.suffix.lower() == ".csv":
            csv_paths.append(path)

    csv_paths = [
        path
        for path in csv_paths
        if path.name != args.map.name and "bo2023_signature_genes_humanized" not in path.name
    ]
    reports = [remap_csv(path, mapping, args.dry_run) for path in csv_paths]
    summary = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "map_path": str(args.map),
        "dry_run": bool(args.dry_run),
        "n_csv_scanned": len(csv_paths),
        "n_changed": sum(1 for item in reports if item["status"] in {"changed", "would_change"}),
        "changed_files": [item for item in reports if item["status"] in {"changed", "would_change"}],
    }
    args.summary.parent.mkdir(parents=True, exist_ok=True)
    args.summary.write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
    print(args.summary)
    print(f"changed={summary['n_changed']} scanned={summary['n_csv_scanned']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
