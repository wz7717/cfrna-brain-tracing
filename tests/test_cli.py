from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from app.pages import tracing_page
from core import cli
from core.model_lock import EXPECTED_MODEL_LOCK_ID, ModelLockError
from core.production_route import production_implementation_parameters
from core.query_input import normalize_query_expression, read_expression_file


def _fake_outputs() -> tuple[dict, dict]:
    network = {
        "results": [
            {"network_id": "N1", "rank": 1, "score": 0.9},
            {"network_id": "N2", "rank": 2, "score": 0.8},
            {"network_id": "N3", "rank": 3, "score": 0.7},
        ],
        "meta": {
            "n_overlap_genes": 180,
            "n_model_genes": 200,
            "overlap_fraction": 0.9,
            "model_lock": {"lock_id": EXPECTED_MODEL_LOCK_ID, "status": "frozen"},
        },
    }
    regions = {
        "results": [{"region_id": "R1", "rank": 1}],
        "meta": {
            "reference_expression_source": "packaged_formal_region_assets",
            "region_resolution_annotation": {
                "group_ranking": [{"resolution_group": "G1", "rank": 1}]
            },
            "model_lock": {"lock_id": EXPECTED_MODEL_LOCK_ID, "status": "frozen"},
        },
    }
    return network, regions


def test_help_is_real_cli_help(capsys) -> None:
    try:
        cli.main(["--help"])
    except SystemExit as exc:
        assert exc.code == 0
    output = capsys.readouterr().out
    assert "query" in output
    assert "models" in output
    assert "validate" in output


def test_version_reports_current_software_release(capsys) -> None:
    try:
        cli.main(["--version"])
    except SystemExit as exc:
        assert exc.code == 0
    release = json.loads((Path(__file__).resolve().parents[1] / "release/v0.1.17/release_manifest.json").read_text(encoding="utf-8"))
    assert capsys.readouterr().out.strip() == f"braintrace {release['version']}"


def test_models_reports_frozen_production_inventory(capsys) -> None:
    assert cli.main(["models"]) == 0
    output = capsys.readouterr().out
    assert EXPECTED_MODEL_LOCK_ID in output
    assert "network_count = 10" in output
    assert "region_count = 110" in output
    assert "beam_count = 120" in output
    assert "allow_development_fallback = false" in output
    assert output.count("- data/models/") == 8


def test_validate_verifies_all_eight_artifacts(capsys) -> None:
    assert cli.main(["validate"]) == 0
    output = capsys.readouterr().out
    assert "BrainTrace frozen model bundle: PASS" in output
    assert EXPECTED_MODEL_LOCK_ID in output
    assert "Artifacts verified: 8/8" in output


def test_validate_fails_closed_on_artifact_drift(monkeypatch, capsys) -> None:
    monkeypatch.setattr(
        cli.production_route,
        "verify_production_route",
        lambda: (_ for _ in ()).throw(ModelLockError("locked artifact drift")),
    )
    assert cli.main(["validate"]) == 1
    assert "locked artifact drift" in capsys.readouterr().err


def test_validate_fails_closed_on_production_parameter_drift(monkeypatch, capsys) -> None:
    monkeypatch.setattr(
        cli.production_route,
        "verify_production_route",
        lambda: (_ for _ in ()).throw(ModelLockError("production parameter drift")),
    )
    assert cli.main(["validate"]) == 1
    assert "production parameter drift" in capsys.readouterr().err


def test_query_input_priority_supports_raw_counts_and_logcpm() -> None:
    raw, raw_source = normalize_query_expression(
        pd.DataFrame({"gene": ["A", "B"], "raw_counts": [10, 20], "logCPM": [1.0, 2.0]})
    )
    logcpm, logcpm_source = normalize_query_expression(
        pd.DataFrame({"symbol": ["A", "B"], "logCPM": [1.0, 2.0]})
    )
    assert raw_source == "raw_counts"
    assert list(raw.columns) == ["gene_symbol", "read_count"]
    assert logcpm_source == "logcpm"
    assert list(logcpm.columns) == ["gene_symbol", "log_tpm"]


def test_query_input_reads_csv_tsv_txt_and_xlsx(tmp_path: Path) -> None:
    frame = pd.DataFrame({"symbol": ["A", "B"], "logCPM": [1.0, 2.0]})
    paths = {
        "csv": tmp_path / "sample.csv",
        "tsv": tmp_path / "sample.tsv",
        "txt": tmp_path / "sample.txt",
        "xlsx": tmp_path / "sample.xlsx",
    }
    frame.to_csv(paths["csv"], index=False)
    frame.to_csv(paths["tsv"], sep="\t", index=False)
    frame.to_csv(paths["txt"], sep="\t", index=False)
    frame.to_excel(paths["xlsx"], index=False)

    for path in paths.values():
        parsed, source = read_expression_file(path)
        assert source == "logcpm"
        assert parsed.to_dict("records") == [
            {"gene_symbol": "A", "log_tpm": 1.0},
            {"gene_symbol": "B", "log_tpm": 2.0},
        ]


def test_streamlit_and_cli_call_the_same_locked_route(monkeypatch, tmp_path: Path) -> None:
    calls = []

    def fake_route(expression, **kwargs):
        calls.append({"expression": expression.copy(), "kwargs": kwargs})
        return _fake_outputs()

    monkeypatch.setattr(cli.locked_inference, "run_locked_three_tier_route", fake_route)
    expression = pd.DataFrame({"gene_symbol": ["A"], "read_count": [1.0]})
    tracing_page._run_locked_bo2023_route(expression, atlas_id=None, topk=30)

    input_path = tmp_path / "sample.tsv"
    output_path = tmp_path / "result.json"
    input_path.write_text("gene_symbol\traw_counts\nA\t1\n", encoding="utf-8")
    assert cli.main(["query", "--input", str(input_path), "--output", str(output_path)]) == 0

    assert len(calls) == 2
    assert calls[0]["kwargs"]["topk"] == 30
    assert calls[0]["kwargs"]["atlas_id"] is None
    assert calls[1]["kwargs"] == {}
    payload = json.loads(output_path.read_text(encoding="utf-8"))
    assert payload["network_top3"][0]["network_id"] == "N1"
    assert payload["resolution_group_top3"][0]["resolution_group"] == "G1"
    assert payload["exact_region_exploratory_top3"][0]["region_id"] == "R1"
    assert payload["meta"]["model_lock"]["lock_id"] == EXPECTED_MODEL_LOCK_ID
    assert payload["meta"]["route_name"] == production_implementation_parameters()["route_name"]
    assert payload["meta"]["query_source"] == "raw_counts"


def test_query_inference_error_is_nonzero(monkeypatch, tmp_path: Path, capsys) -> None:
    monkeypatch.setattr(
        cli.locked_inference,
        "run_locked_three_tier_route",
        lambda expression: (_ for _ in ()).throw(ValueError("Insufficient Network model-gene overlap: 1/200.")),
    )
    input_path = tmp_path / "sample.csv"
    output_path = tmp_path / "result.json"
    input_path.write_text("gene_symbol,logCPM\nA,1\n", encoding="utf-8")
    assert cli.main(["query", "--input", str(input_path), "--output", str(output_path)]) == 1
    assert "Insufficient Network model-gene overlap" in capsys.readouterr().err
    assert not output_path.exists()
