#!/usr/bin/env python
"""Fail closed when a freshly regenerated TCGA/BraTS truth basis drifts."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any


EXPECTED: dict[tuple[str, str], tuple[int, int]] = {
    ("network", "center"): (19, 65),
    ("network", "core"): (12, 65),
    ("network", "edema"): (15, 63),
    ("network", "whole_tumor"): (10, 65),
    ("broad", "center"): (32, 65),
    ("broad", "core"): (45, 65),
    ("broad", "edema"): (52, 63),
    ("broad", "whole_tumor"): (46, 65),
}
EXPECTED_EXCLUDED = {"TCGA-HT-7680", "TCGA-HT-7686"}


def _boolean(value: object) -> bool:
    normalized = str(value).strip().lower()
    if normalized in {"true", "1", "yes"}:
        return True
    if normalized in {"false", "0", "", "nan"}:
        return False
    raise ValueError(f"invalid boolean value: {value!r}")


def validate(path: Path) -> dict[str, Any]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        rows = list(csv.DictReader(handle))
    errors: list[str] = []
    if len(rows) != 65:
        errors.append(f"expected 65 paired cases, found {len(rows)}")
    patients = [str(row.get("patient_barcode", "")).upper() for row in rows]
    if len(patients) != len(set(patients)):
        errors.append("patient_barcode values are not unique")
    edema_eligible = [
        row
        for row in rows
        if float(row["edema_voxels"]) > 0 and str(row["edema_network_dominant"]) != "out_of_scope"
    ]
    excluded = {
        str(row.get("patient_barcode", "")).upper()
        for row in rows
        if row not in edema_eligible
    }
    if excluded != EXPECTED_EXCLUDED:
        errors.append(f"edema exclusions changed: {sorted(excluded)}")
    observed: dict[str, dict[str, dict[str, int]]] = {}
    for (level, truth_basis), expected in EXPECTED.items():
        eligible = edema_eligible if truth_basis == "edema" else rows
        column = f"{truth_basis}_{level}_top3_strict"
        try:
            value = (sum(_boolean(row[column]) for row in eligible), len(eligible))
        except (KeyError, ValueError) as exc:
            errors.append(f"{column}: {exc}")
            continue
        observed.setdefault(level, {})[truth_basis] = {"correct": value[0], "n": value[1]}
        if value != expected:
            errors.append(f"{column} changed: {value[0]}/{value[1]}, expected {expected[0]}/{expected[1]}")
    result = {
        "schema": "braintrace.tcga_brats_truth_basis_run_verification.v1",
        "status": "PASS" if not errors else "BLOCKED: SCIENTIFIC_OUTPUT_DRIFT",
        "total_paired_cases": len(rows),
        "primary_edema_comparator_n": len(edema_eligible),
        "excluded_cases": sorted(excluded),
        "strict_top3": observed,
        "errors": errors,
    }
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    payload = validate(args.input)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8", newline="\n")
    print(json.dumps(payload))
    return 0 if payload["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
