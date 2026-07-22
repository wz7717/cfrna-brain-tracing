from __future__ import annotations

import gzip
import json
import os
from pathlib import Path

import pandas as pd
import pytest

from scripts import benchmark_real_input_inference as harness


def test_windows_process_memory_probe_returns_working_set() -> None:
    memory = harness.process_memory()
    if os.name == "nt":
        assert memory.working_set_bytes is not None and memory.working_set_bytes > 0
        assert memory.peak_working_set_bytes is not None and memory.peak_working_set_bytes > 0


def test_harness_checks_current_canonical_contract() -> None:
    contract = harness.assert_canonical_contract()
    assert contract["canonical_regions"] == 110
    assert contract["networks"] == 10
    assert contract["network_top3_beams"] == 120
    assert contract["reference_genes"] == 21668
    assert contract["reference_source"] == harness.PACKAGED_REFERENCE_SOURCE


def test_read_counts_subset_preserves_fixed_header_order(tmp_path: Path) -> None:
    path = tmp_path / harness.COUNT_FILE_NAME
    with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
        handle.write("Geneid,SAMPLE-C,SAMPLE-A,SAMPLE-B,SAMPLE-D\n")
        handle.write("GENE1,1,2,3,4\n")
        handle.write("GENE2,5,6,7,8\n")

    frame, samples, total = harness.read_counts_subset(path, 3)

    assert samples == ["SAMPLE-C", "SAMPLE-A", "SAMPLE-B"]
    assert total == 4
    assert frame.shape == (2, 3)


def test_frozen_route_explicitly_disables_pairwise_rescue(monkeypatch) -> None:
    observed: dict[str, object] = {}

    def fake_network(expression, **kwargs):
        observed.update(kwargs)
        return {
            "results": [{"network_id": f"N{i}"} for i in range(10)],
            "meta": {"pairwise_rescue": {"enabled": False, "switched": False}},
        }

    def fake_region(expression, network_output, db_path, atlas_id, topk):
        observed["db_path"] = db_path
        return {"results": [{"region_id": "R1"}], "meta": {"network_beam": ["N0", "N1", "N2"]}}

    monkeypatch.setattr(harness, "trace_network_expression", fake_network)
    monkeypatch.setattr(harness, "trace_bo2023_secondary_regions", fake_region)

    harness.run_frozen_route(pd.DataFrame({"gene_symbol": ["G1"], "read_count": [1]}))

    assert observed["enable_pairwise_rescue"] is False
    assert observed["project_to_vsd"] is True
    assert observed["db_path"] == "__packaged_reference_only__"


def test_smoke_outputs_are_anonymized_and_formal_run_stays_disabled(tmp_path: Path, monkeypatch) -> None:
    data_dir = tmp_path / "input"
    data_dir.mkdir()
    count_path = data_dir / harness.COUNT_FILE_NAME
    count_path.write_bytes(b"placeholder")
    raw_ids = ["PATIENT-ALPHA", "PATIENT-BETA", "PATIENT-GAMMA"]
    counts = pd.DataFrame(
        {sample: [1.0, 2.0] for sample in raw_ids},
        index=pd.Index(["G1", "G2"], name="Geneid"),
    )

    monkeypatch.setattr(harness, "read_counts_subset", lambda path, n: (counts, raw_ids, 51))
    monkeypatch.setattr(
        harness,
        "assert_canonical_contract",
        lambda: {"canonical_regions": 110, "networks": 10, "network_top3_beams": 120},
    )
    monkeypatch.setattr(
        harness,
        "run_frozen_route",
        lambda expression: (
            {
                "results": [{"network_id": f"N{i}"} for i in range(10)],
                "meta": {"pairwise_rescue": {"enabled": False, "switched": False}},
            },
            {"results": [{"region_id": "R1"}], "meta": {"network_beam": ["N0", "N1", "N2"]}},
        ),
    )
    monkeypatch.setattr(harness, "environment_manifest", lambda: {"git_commit": "test"})
    monkeypatch.setattr(harness, "process_memory", lambda: harness.ProcessMemory(100, 200))

    outdir = tmp_path / "out"
    manifest = harness.run_smoke(data_dir.resolve(), outdir)

    timing_text = (outdir / "timing.csv").read_text(encoding="utf-8")
    manifest_text = (outdir / "manifest.json").read_text(encoding="utf-8")
    for raw_id in raw_ids:
        assert raw_id not in timing_text
        assert raw_id not in manifest_text
    assert manifest["preregistration"]["formal_workload_sizes"] == [1, 8, 51]
    assert manifest["preregistration"]["formal_execution_authorized"] is False
    assert manifest["smoke"]["pairwise_rescue_switches"] == 0
    assert json.loads(manifest_text)["input"]["clinical_metadata_copied"] is False


def test_formal_workload_requires_explicit_authorization(tmp_path: Path) -> None:
    with pytest.raises(PermissionError, match="explicit --authorize-formal"):
        harness.run_formal_workload(
            tmp_path.resolve(),
            tmp_path / "out",
            1,
            authorized=False,
        )


def test_formal_gate_records_fixed_workload_without_identifiers(
    tmp_path: Path,
    monkeypatch,
) -> None:
    data_dir = tmp_path / "input"
    data_dir.mkdir()
    count_path = data_dir / harness.COUNT_FILE_NAME
    count_path.write_bytes(b"placeholder")
    raw_ids = [f"PATIENT-{index}" for index in range(8)]
    counts = pd.DataFrame(
        {sample: [1.0, 2.0] for sample in raw_ids},
        index=pd.Index(["G1", "G2"], name="Geneid"),
    )
    contract = {
        "canonical_regions": 110,
        "networks": 10,
        "network_top3_beams": 120,
        "reference_genes": 21668,
        "reference_source": harness.PACKAGED_REFERENCE_SOURCE,
    }

    monkeypatch.setattr(harness, "read_counts_subset", lambda path, n: (counts.iloc[:, :n], raw_ids[:n], 51))
    monkeypatch.setattr(harness, "assert_canonical_contract", lambda: contract)
    monkeypatch.setattr(
        harness,
        "run_frozen_route",
        lambda expression: (
            {
                "results": [{"network_id": f"N{i}"} for i in range(10)],
                "meta": {"pairwise_rescue": {"enabled": False, "switched": False}},
            },
            {
                "results": [{"region_id": "R1"}],
                "meta": {"network_beam": ["N0", "N1", "N2"]},
            },
        ),
    )
    monkeypatch.setattr(harness, "environment_manifest", lambda: {"git_commit": "test"})
    monkeypatch.setattr(harness, "process_memory", lambda: harness.ProcessMemory(100, 200))

    outdir = tmp_path / "formal"
    manifest = harness.run_formal_workload(
        data_dir.resolve(),
        outdir,
        8,
        authorized=True,
    )

    text = (outdir / "manifest.json").read_text(encoding="utf-8")
    text += (outdir / "timing.csv").read_text(encoding="utf-8")
    assert all(raw_id not in text for raw_id in raw_ids)
    assert manifest["input"]["workload_samples"] == 8
    assert manifest["input"]["clinical_metadata_read"] is False
    assert manifest["input"]["persistent_database_read"] is False
    assert manifest["preregistration"]["formal_warm_repeats"] == 3
    assert manifest["preregistration"]["concurrency"] == 1
    assert len(manifest["warm"]["repeats"]) == 3
    assert manifest["warm"]["aggregate"]["samples"] == 24
    assert manifest["warm"]["aggregate"]["pairwise_rescue_switches"] == 0
