#!/usr/bin/env python
"""Validate a new formal GSE189919 benchmark without comparing hardware timing."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Any


EXPECTED_ROUTE = "projected_vsd_network_top3_logcpm_resolution_local_exact"
EXPECTED_WORKLOAD = 51
EXPECTED_WARM_REPEATS = 3
EXPECTED_BO2023_SOURCE_FEATURE_ROWS = 28_415


def _positive(value: object) -> bool:
    try:
        return math.isfinite(float(value)) and float(value) > 0
    except (TypeError, ValueError):
        return False


def _count_feature_rows(path: Path) -> int:
    with path.open(encoding="utf-8", errors="strict") as handle:
        next(handle, None)
        return sum(1 for _ in handle)


def validate(manifest_path: Path, timing_path: Path, bo2023_counts: Path | None = None) -> dict[str, Any]:
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    errors: list[str] = []
    if manifest.get("status") != "formal_real_input_performance_gate":
        errors.append("manifest is not a formal real-input benchmark")
    if manifest.get("route") != EXPECTED_ROUTE:
        errors.append("frozen route differs from the registered benchmark route")
    input_record = manifest.get("input", {})
    preregistration = manifest.get("preregistration", {})
    warm = manifest.get("warm", {})
    repeats = warm.get("repeats", [])
    aggregate = warm.get("aggregate", {})
    if input_record.get("workload_samples") != EXPECTED_WORKLOAD:
        errors.append("workload is not 51 profiles")
    if input_record.get("total_samples_in_header") != EXPECTED_WORKLOAD:
        errors.append("input header does not contain 51 profiles")
    if preregistration.get("formal_warm_repeats") != EXPECTED_WARM_REPEATS or len(repeats) != EXPECTED_WARM_REPEATS:
        errors.append("benchmark does not contain exactly three warm repeats")
    if aggregate.get("samples") != EXPECTED_WORKLOAD * EXPECTED_WARM_REPEATS:
        errors.append("warm event count is not 153")
    if preregistration.get("formal_execution_authorized") is not True:
        errors.append("formal workload authorization is not recorded")
    for label, timing in [("cold", manifest.get("cold", {}).get("timing", {})), ("warm aggregate", aggregate)]:
        for key in ("wall_total_seconds", "wall_seconds_per_sample", "sample_time_p50_seconds", "sample_time_p95_seconds"):
            if not _positive(timing.get(key)):
                errors.append(f"{label} {key} is not a finite positive value")
    for repeat in repeats:
        timing = repeat.get("timing", {}) if isinstance(repeat, dict) else {}
        if timing.get("samples") != EXPECTED_WORKLOAD:
            errors.append("a warm repeat does not contain 51 profile events")
        if not _positive(timing.get("wall_total_seconds")):
            errors.append("a warm repeat has non-positive wall time")
    with timing_path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    warm_rows = [row for row in rows if row.get("phase") == "warm"]
    cold_rows = [row for row in rows if row.get("phase") == "cold"]
    if len(cold_rows) != EXPECTED_WORKLOAD or len(warm_rows) != EXPECTED_WORKLOAD * EXPECTED_WARM_REPEATS:
        errors.append("timing CSV event counts do not match the manifest")
    source_feature_rows: int | None = None
    if bo2023_counts is not None:
        source_feature_rows = _count_feature_rows(bo2023_counts)
        if source_feature_rows != EXPECTED_BO2023_SOURCE_FEATURE_ROWS:
            errors.append("frozen Bo2023 source feature count is not 28,415")
    return {
        "schema": "braintrace.gse189919_benchmark_run_verification.v1",
        "status": "PASS" if not errors else "FAIL",
        "route": manifest.get("route"),
        "workload_profiles": input_record.get("workload_samples"),
        "warm_repeats": len(repeats),
        "warm_events": aggregate.get("samples"),
        "timing_csv_cold_events": len(cold_rows),
        "timing_csv_warm_events": len(warm_rows),
        "frozen_bo2023_source_feature_rows": source_feature_rows,
        "gene_dimension_semantics": "28,415 is the verified Bo2023 frozen source-feature count; GSE189919 raw source rows are tracked separately and are not relabelled.",
        "hardware_timing_comparison": "not performed; only finite positive timing/memory schema is required on a new environment",
        "errors": errors,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--timing", type=Path, required=True)
    parser.add_argument("--bo2023-counts", type=Path, help="Verified Bo2023 source count matrix; validates its 28,415 source-feature rows.")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    payload = validate(args.manifest, args.timing, args.bo2023_counts)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8", newline="\n")
    print(json.dumps(payload))
    return 0 if payload["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
