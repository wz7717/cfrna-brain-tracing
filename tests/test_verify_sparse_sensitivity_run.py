from __future__ import annotations

import csv
import json
from pathlib import Path

from scripts.verify_sparse_sensitivity_run import (
    EXPECTED_BOOTSTRAP_SEED,
    EXPECTED_NONBASELINE,
    EXPECTED_PERTURBATION_SEED,
    EXPECTED_REPLICATES,
    EXPECTED_SAMPLES,
    EXPECTED_SCENARIOS,
    validate,
)


def _write_contract(tmp_path: Path) -> tuple[Path, Path, Path, Path]:
    methods = tmp_path / "methods.json"
    methods.write_text(json.dumps({
        "replicates_per_nonbaseline_scenario": EXPECTED_REPLICATES,
        "perturbation_seed": EXPECTED_PERTURBATION_SEED,
        "bootstrap_seed": EXPECTED_BOOTSTRAP_SEED,
        "frozen_route": {"enable_pairwise_rescue": False},
        "scenarios": [
            {"name": name, "depth_fraction": value[0], "target_gene_retention": value[1]}
            for name, value in EXPECTED_SCENARIOS.items()
        ],
    }), encoding="utf-8")
    per_repeat = tmp_path / "per_repeat.csv"
    with per_repeat.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["scenario", "replicate", "n_samples"])
        writer.writeheader()
        writer.writerow({"scenario": "baseline", "replicate": 0, "n_samples": EXPECTED_SAMPLES})
        for scenario in EXPECTED_NONBASELINE:
            for repeat in range(EXPECTED_REPLICATES):
                writer.writerow({"scenario": scenario, "replicate": repeat, "n_samples": EXPECTED_SAMPLES})
    detail = tmp_path / "detail.csv"
    with detail.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["scenario", "replicate", "rng_seed", "sample_id"])
        writer.writeheader()
        for sample in range(EXPECTED_SAMPLES):
            writer.writerow({"scenario": "baseline", "replicate": 0, "rng_seed": sample, "sample_id": sample})
        for scenario_index, scenario in enumerate(EXPECTED_NONBASELINE, start=1):
            for repeat in range(EXPECTED_REPLICATES):
                for sample in range(EXPECTED_SAMPLES):
                    writer.writerow({
                        "scenario": scenario,
                        "replicate": repeat,
                        "rng_seed": scenario_index * 1_000_000 + repeat * 10_000 + sample,
                        "sample_id": sample,
                    })
    gate = tmp_path / "baseline_gate.json"
    gate.write_text(json.dumps({"passed": True, "n_donors": 9}), encoding="utf-8")
    return methods, per_repeat, detail, gate


def test_sparse_verifier_rejects_missing_repeat(tmp_path: Path) -> None:
    methods, per_repeat, detail, gate = _write_contract(tmp_path)
    assert validate(methods, per_repeat, detail, gate)["status"] == "PASS"
    rows = list(csv.DictReader(per_repeat.open(encoding="utf-8")))
    rows = [row for row in rows if not (row["scenario"] == "severe" and row["replicate"] == "29")]
    with per_repeat.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["scenario", "replicate", "n_samples"])
        writer.writeheader()
        writer.writerows(rows)
    payload = validate(methods, per_repeat, detail, gate)
    assert payload["status"] == "BLOCKED: SCIENTIFIC_OUTPUT_DRIFT"
    assert any("severe does not contain exactly repeats 0..29" in error for error in payload["errors"])
