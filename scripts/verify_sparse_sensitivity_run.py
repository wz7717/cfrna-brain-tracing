#!/usr/bin/env python
"""Fail closed if the regenerated sparse-sensitivity protocol changes.

This verifier deliberately checks the registered execution contract rather than
trying to reproduce hardware-independent Monte Carlo estimates byte-for-byte.
The baseline scientific endpoint is guarded by ``baseline_gate.json`` produced
by the locked route; every non-baseline condition must retain all 30 frozen-seed
repeats and their deterministic per-sample seeds.
"""

from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path
from typing import Any


EXPECTED_REPLICATES = 30
EXPECTED_PERTURBATION_SEED = 20260711
EXPECTED_BOOTSTRAP_SEED = 20260716
EXPECTED_SAMPLES = 819
EXPECTED_NONBASELINE = ("mild", "moderate", "severe", "extreme")
EXPECTED_SCENARIOS = {
    "baseline": (1.00, 1.00),
    "mild": (0.50, 0.80),
    "moderate": (0.20, 0.60),
    "severe": (0.05, 0.40),
    "extreme": (0.01, 0.20),
}


def _as_int(row: dict[str, str], key: str, errors: list[str]) -> int | None:
    try:
        return int(str(row[key]).strip())
    except (KeyError, TypeError, ValueError):
        errors.append(f"invalid integer {key!r}")
        return None


def _as_float(value: object) -> float | None:
    try:
        return float(str(value).strip())
    except (TypeError, ValueError):
        return None


def _scenario_contract(methods: dict[str, Any], errors: list[str]) -> None:
    scenarios = methods.get("scenarios")
    if not isinstance(scenarios, list):
        errors.append("methods scenarios are absent")
        return
    observed: dict[str, tuple[float, float]] = {}
    for item in scenarios:
        if not isinstance(item, dict):
            errors.append("methods scenario entry is not an object")
            continue
        name = item.get("name")
        depth, retention = _as_float(item.get("depth_fraction")), _as_float(item.get("target_gene_retention"))
        if not isinstance(name, str) or depth is None or retention is None:
            errors.append("methods scenario entry is malformed")
            continue
        observed[name] = (depth, retention)
    if observed != EXPECTED_SCENARIOS:
        errors.append("sparse scenario/depth/retention contract changed")


def validate(methods_path: Path, per_repeat_path: Path, detail_path: Path, baseline_gate_path: Path) -> dict[str, Any]:
    """Validate the frozen sparse protocol and return a portable audit payload."""

    errors: list[str] = []
    try:
        methods = json.loads(methods_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        methods = {}
        errors.append(f"cannot read methods: {exc}")
    if methods.get("replicates_per_nonbaseline_scenario") != EXPECTED_REPLICATES:
        errors.append("nonbaseline replicate count is not 30")
    if methods.get("perturbation_seed") != EXPECTED_PERTURBATION_SEED:
        errors.append("perturbation seed is not 20260711")
    if methods.get("bootstrap_seed") != EXPECTED_BOOTSTRAP_SEED:
        errors.append("bootstrap seed is not 20260716")
    if methods.get("frozen_route", {}).get("enable_pairwise_rescue") is not False:
        errors.append("pairwise rescue must remain disabled")
    _scenario_contract(methods, errors)

    per_repeat_replicates: dict[str, set[int]] = defaultdict(set)
    try:
        with per_repeat_path.open(newline="", encoding="utf-8-sig") as handle:
            per_repeat_rows = list(csv.DictReader(handle))
    except OSError as exc:
        per_repeat_rows = []
        errors.append(f"cannot read per-repeat CSV: {exc}")
    if not per_repeat_rows:
        errors.append("per-repeat CSV has no rows")
    for row in per_repeat_rows:
        scenario = row.get("scenario", "")
        replicate = _as_int(row, "replicate", errors)
        if scenario and replicate is not None:
            per_repeat_replicates[scenario].add(replicate)
    expected_repeat_set = set(range(EXPECTED_REPLICATES))
    for scenario in EXPECTED_NONBASELINE:
        if per_repeat_replicates.get(scenario) != expected_repeat_set:
            errors.append(f"{scenario} does not contain exactly repeats 0..29")
    if per_repeat_replicates.get("baseline") != {0}:
        errors.append("baseline must contain only replicate 0")

    detail_rows_by_scenario: dict[str, int] = defaultdict(int)
    sample_ids_by_scenario_repeat: dict[tuple[str, int], set[str]] = defaultdict(set)
    seeds_by_scenario_repeat: dict[tuple[str, int], set[int]] = defaultdict(set)
    seed_pairs: set[tuple[str, int, str, int]] = set()
    duplicate_seed_rows = 0
    try:
        with detail_path.open(newline="", encoding="utf-8-sig") as handle:
            detail_rows = csv.DictReader(handle)
            for row in detail_rows:
                scenario = row.get("scenario", "")
                replicate = _as_int(row, "replicate", errors)
                seed = _as_int(row, "rng_seed", errors)
                sample_id = str(row.get("sample_id", "")).strip()
                if not scenario or replicate is None or seed is None or not sample_id:
                    errors.append("detail CSV contains an incomplete protocol row")
                    continue
                detail_rows_by_scenario[scenario] += 1
                sample_ids_by_scenario_repeat[(scenario, replicate)].add(sample_id)
                seeds_by_scenario_repeat[(scenario, replicate)].add(seed)
                key = (scenario, replicate, sample_id, seed)
                if key in seed_pairs:
                    duplicate_seed_rows += 1
                seed_pairs.add(key)
    except OSError as exc:
        errors.append(f"cannot read sample-detail CSV: {exc}")
    if duplicate_seed_rows:
        errors.append("sample-detail CSV contains duplicate deterministic seed rows")
    for scenario in EXPECTED_NONBASELINE:
        expected_rows = EXPECTED_SAMPLES * EXPECTED_REPLICATES
        if detail_rows_by_scenario.get(scenario) != expected_rows:
            errors.append(f"{scenario} detail row count is not {expected_rows}")
        for repeat in range(EXPECTED_REPLICATES):
            sample_ids = sample_ids_by_scenario_repeat.get((scenario, repeat), set())
            if len(sample_ids) != EXPECTED_SAMPLES:
                errors.append(f"{scenario} repeat {repeat} does not contain 819 unique samples")
            if len(seeds_by_scenario_repeat.get((scenario, repeat), set())) != EXPECTED_SAMPLES:
                errors.append(f"{scenario} repeat {repeat} does not contain 819 independent deterministic seeds")
    if detail_rows_by_scenario.get("baseline") != EXPECTED_SAMPLES:
        errors.append("baseline detail row count is not 819")
    if len(sample_ids_by_scenario_repeat.get(("baseline", 0), set())) != EXPECTED_SAMPLES:
        errors.append("baseline does not contain 819 unique samples")

    try:
        baseline_gate = json.loads(baseline_gate_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        baseline_gate = {}
        errors.append(f"cannot read baseline gate: {exc}")
    if baseline_gate.get("passed") is not True:
        errors.append("locked baseline scientific endpoint did not pass")
    if baseline_gate.get("n_donors") != 9:
        errors.append("locked baseline donor count is not 9")

    status = "PASS" if not errors else "BLOCKED: SCIENTIFIC_OUTPUT_DRIFT"
    return {
        "schema": "braintrace.sparse_sensitivity_run_verification.v1",
        "status": status,
        "replicates_per_nonbaseline_scenario": EXPECTED_REPLICATES,
        "perturbation_seed": EXPECTED_PERTURBATION_SEED,
        "bootstrap_seed": EXPECTED_BOOTSTRAP_SEED,
        "expected_samples_per_repeat": EXPECTED_SAMPLES,
        "scenario_replicates": {name: sorted(values) for name, values in sorted(per_repeat_replicates.items())},
        "detail_rows_by_scenario": dict(sorted(detail_rows_by_scenario.items())),
        "baseline_gate_passed": baseline_gate.get("passed") is True,
        "errors": errors,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--methods", type=Path, required=True)
    parser.add_argument("--per-repeat", type=Path, required=True)
    parser.add_argument("--detail", type=Path, required=True)
    parser.add_argument("--baseline-gate", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    payload = validate(args.methods, args.per_repeat, args.detail, args.baseline_gate)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8", newline="\n")
    print(json.dumps({"status": payload["status"]}))
    return 0 if payload["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
