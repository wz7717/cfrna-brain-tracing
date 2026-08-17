from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import tomllib

from app import main as app_main
from app.pages import tracing_page
from core import cli
from core.model_lock import EXPECTED_MODEL_LOCK_ID
from examples import verify_examples


ROOT = Path(__file__).resolve().parents[1]
COUNTS = ROOT / "examples/braintrace_example_counts.tsv"
LOGCPM = ROOT / "examples/braintrace_example_logcpm.tsv"


def _fake_outputs() -> tuple[dict, dict]:
    return (
        {
            "results": [{"network_id": "N1", "rank": 1, "score": 1.0}],
            "meta": {
                "n_overlap_genes": 200,
                "n_model_genes": 200,
                "overlap_fraction": 1.0,
                "model_lock": {"lock_id": EXPECTED_MODEL_LOCK_ID, "status": "frozen"},
            },
        },
        {
            "results": [{"region_id": "R1", "rank": 1}],
            "meta": {
                "reference_expression_source": "packaged_formal_region_assets",
                "region_resolution_annotation": {
                    "group_ranking": [{"resolution_group": "G1", "rank": 1}]
                },
                "model_lock": {"lock_id": EXPECTED_MODEL_LOCK_ID, "status": "frozen"},
            },
        },
    )


def test_current_package_version() -> None:
    metadata = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))
    assert metadata["project"]["version"] == "0.1.16"
    assert cli._package_version() == "0.1.16"


def test_public_example_counts_runs() -> None:
    expression, source, network, regions = verify_examples.run_example(COUNTS)
    assert source == "raw_counts"
    assert len(expression) == 21_668
    assert network["results"] and regions["results"]


def test_public_example_logcpm_runs() -> None:
    expression, source, network, regions = verify_examples.run_example(LOGCPM)
    assert source == "logcpm"
    assert len(expression) == 21_668
    assert network["results"] and regions["results"]


def test_public_example_expected_output() -> None:
    for args in verify_examples.EXAMPLES:
        verify_examples.run_and_verify(*args)


def test_public_example_does_not_require_private_data() -> None:
    for path in (COUNTS, LOGCPM):
        _, _, _, region_out = verify_examples.run_example(path)
        verify_examples._assert_public_packaged_reference(region_out)
    source = (ROOT / "examples/verify_examples.py").read_text(encoding="utf-8").lower()
    assert "author-package" not in source
    assert "processed_count" not in source
    assert "processed_vsd" not in source


def test_public_example_network_overlap_gate() -> None:
    for path in (COUNTS, LOGCPM):
        expression, _ = verify_examples.read_expression_file(path)
        network_overlap, regional_overlap = verify_examples._input_overlap(expression)
        assert network_overlap == 200
        assert network_overlap / 200 >= 0.50
        assert regional_overlap >= 20


def test_help_page_exists() -> None:
    assert "help" in app_main.PAGES
    assert app_main.PAGES["help"]["func"] == "app.pages.help_page:display_help_page"
    help_text = (ROOT / "app/pages/help_page.py").read_text(encoding="utf-8")
    for required in ("≥100/200", "at least 20 genes", "exploratory only", "patient-level anatomical truth"):
        assert required in help_text


def test_load_example_uses_locked_route(monkeypatch) -> None:
    calls = []

    def fake_route(expression, **kwargs):
        calls.append((expression.copy(), kwargs))
        return _fake_outputs()

    monkeypatch.setattr(tracing_page.locked_inference, "run_locked_three_tier_route", fake_route)
    expression, source, network, regions = tracing_page._run_public_example_locked_route()
    assert source == "raw_counts"
    assert len(expression) == 21_668
    assert network["results"] and regions["results"]
    assert len(calls) == 1
    assert calls[0][1]["atlas_id"] is None


def test_web_cli_example_same_locked_route(monkeypatch, tmp_path: Path) -> None:
    calls = []

    def fake_route(expression, **kwargs):
        calls.append(kwargs)
        return _fake_outputs()

    monkeypatch.setattr(cli.locked_inference, "run_locked_three_tier_route", fake_route)
    tracing_page._run_public_example_locked_route()

    cli_output = tmp_path / "cli.json"
    assert cli.main(["query", "--input", str(COUNTS), "--output", str(cli_output)]) == 0
    assert json.loads(cli_output.read_text(encoding="utf-8"))["meta"]["model_lock"]["lock_id"] == EXPECTED_MODEL_LOCK_ID

    verify_examples.run_example(COUNTS)
    assert len(calls) == 3
    assert calls[0]["atlas_id"] is None
    assert calls[1] == {}
    assert calls[2] == {}


def test_streamlit_upload_limits_use_defaults() -> None:
    config = tomllib.loads((ROOT / ".streamlit/config.toml").read_text(encoding="utf-8"))
    server = config["server"]
    assert "maxUploadSize" not in server
    assert "maxMessageSize" not in server
