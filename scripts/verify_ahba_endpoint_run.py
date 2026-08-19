#!/usr/bin/env python
"""Compare a regenerated AHBA endpoint ledger to the frozen scientific contract."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any


ROUTE = "hybrid_projected_network_logcpm_exact"
EXPECTED = {
    "network_top1_hit": (165, 223),
    "network_top3_hit": (211, 223),
    "group_top1_hit": (37, 88),
    "group_top3_hit": (60, 88),
    "region_top1_exact_hit": (24, 88),
    "region_top3_exact_hit": (40, 88),
}


def _boolean(value: object) -> bool:
    text = str(value).strip().lower()
    if text in {"true", "1"}:
        return True
    if text in {"false", "0"}:
        return False
    raise ValueError(f"invalid boolean value: {value!r}")


def validate(path: Path) -> dict[str, Any]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        all_rows = list(csv.DictReader(handle))
    rows = [row for row in all_rows if row.get("route") == ROUTE]
    errors: list[str] = []
    if len(rows) != 231:
        errors.append(f"expected 231 replicate-collapsed hybrid rows, found {len(rows)}")
    results: dict[str, dict[str, int]] = {}
    for column, expected in EXPECTED.items():
        values = [row[column] for row in rows if str(row.get(column, "")).strip()]
        try:
            observed = (sum(_boolean(value) for value in values), len(values))
        except (KeyError, ValueError) as exc:
            errors.append(f"{column}: {exc}")
            continue
        results[column] = {"correct": observed[0], "n": observed[1]}
        if observed != expected:
            errors.append(f"{column} changed: {observed[0]}/{observed[1]}, expected {expected[0]}/{expected[1]}")
    return {
        "schema": "braintrace.ahba_endpoint_run_verification.v1",
        "status": "PASS" if not errors else "BLOCKED: SCIENTIFIC_OUTPUT_DRIFT",
        "route": ROUTE,
        "replicate_collapsed_tissues": len(rows),
        "endpoint_counts": results,
        "errors": errors,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    payload = validate(args.input)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8", newline="\n")
    print(json.dumps({"status": payload["status"]}))
    return 0 if payload["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
