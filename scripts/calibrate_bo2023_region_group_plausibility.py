#!/usr/bin/env python
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MODEL = ROOT / "data" / "models" / "bo2023_region_resolution_groups.json"
DEFAULT_CSV = ROOT / "docs" / "deliverables" / "bo2023_region_group_plausibility_audit.csv"


def group_tier(size: int, members: list[str]) -> tuple[str, list[str]]:
    flags: list[str] = ["same_network_only"]
    if size <= 1:
        return "single_region_low_resolution", flags
    if size <= 4:
        flags.append("compact_group_size")
        return "strong_internal_support", flags
    if size <= 6:
        flags.append("moderate_group_size")
        return "moderate_internal_support", flags
    flags.append("large_group_size")
    return "requires_manual_review", flags


def main() -> int:
    parser = argparse.ArgumentParser(description="Add internal plausibility/risk calibration to Bo2023 Region groups.")
    parser.add_argument("--model", type=Path, default=DEFAULT_MODEL)
    parser.add_argument("--csv", type=Path, default=DEFAULT_CSV)
    args = parser.parse_args()

    model = json.loads(args.model.read_text(encoding="utf-8"))
    entries: dict[str, dict[str, Any]] = model.get("entries", {})
    grouped: dict[tuple[str, str], list[dict[str, Any]]] = {}
    for entry in entries.values():
        grouped.setdefault((str(entry.get("network_id", "")), str(entry.get("resolution_group", ""))), []).append(entry)

    group_calibration: dict[str, dict[str, Any]] = {}
    rows: list[dict[str, Any]] = []
    for (network, group_id), items in sorted(grouped.items()):
        members = sorted({str(member) for item in items for member in item.get("group_members", [])})
        size = len(members)
        tier, flags = group_tier(size, members)
        if any(member.upper().startswith("NA") for member in members):
            flags.append("contains_unresolved_region_name")
            if tier == "strong_internal_support":
                tier = "moderate_internal_support"
        if size >= 7:
            flags.append("do_not_interpret_as_single_functional_area")
        merged = size > 1
        key = f"{network}||{group_id}"
        calibration = {
            "network_id": network,
            "resolution_group": group_id,
            "group_size": size,
            "group_members": members,
            "is_merged_group": merged,
            "plausibility_tier": tier,
            "calibration_flags": flags,
            "interpretation": (
                "Same-network expression/confusion group; use as resolution-aware candidate group, "
                "not as proof of identical neurobiological function."
                if merged
                else "Single Region entry; low_resolution may reflect low sample count or high nearest centroid similarity."
            ),
        }
        group_calibration[key] = calibration
        for item in items:
            item["group_plausibility_tier"] = tier
            item["group_calibration_flags"] = flags
        rows.append(
            {
                "network_id": network,
                "resolution_group": group_id,
                "group_size": size,
                "group_members": " | ".join(members),
                "is_merged_group": "yes" if merged else "no",
                "plausibility_tier": tier,
                "calibration_flags": " | ".join(flags),
            }
        )

    model["group_calibration"] = group_calibration
    model["calibration_summary"] = {
        "n_groups": len(grouped),
        "n_merged_groups": sum(1 for row in rows if row["is_merged_group"] == "yes"),
        "tier_counts": {
            tier: sum(1 for row in rows if row["plausibility_tier"] == tier)
            for tier in sorted({row["plausibility_tier"] for row in rows})
        },
        "scope": "internal calibration only; no claim that merged Regions are functionally identical",
    }
    args.model.write_text(json.dumps(model, ensure_ascii=False, indent=2), encoding="utf-8")

    args.csv.parent.mkdir(parents=True, exist_ok=True)
    with args.csv.open("w", encoding="utf-8-sig", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    print(json.dumps({"model": str(args.model), "csv": str(args.csv), **model["calibration_summary"]}, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
