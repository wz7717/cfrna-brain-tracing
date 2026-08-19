from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from scripts import run_gse189919_latest_main_route as runner
from scripts import benchmark_real_input_inference as benchmark


def test_gse_runner_explicitly_disables_pairwise_rescue(monkeypatch) -> None:
    observed: dict[str, object] = {}

    def fake_network(expression, **kwargs):
        observed.update(kwargs)
        return {
            "results": [{"network_id": "N1"}, {"network_id": "N2"}, {"network_id": "N3"}],
            "meta": {"pairwise_rescue": {"enabled": False, "switched": False}},
        }

    def fake_region(expression, network_output, db_path, atlas_id, topk):
        observed.update({"db_path": db_path, "atlas_id": atlas_id, "topk": topk})
        return {"results": [{"region_id": "R1"}], "meta": {}}

    monkeypatch.setattr(runner, "trace_network_expression", fake_network)
    monkeypatch.setattr(runner, "trace_bo2023_secondary_regions", fake_region)

    runner.run_frozen_sample(
        pd.DataFrame({"gene_symbol": ["G1"], "read_count": [1]}),
        min_network_overlap=0.5,
        db_path=Path("unused.db"),
        atlas_id=1,
        topk_regions=30,
    )

    assert observed["enable_pairwise_rescue"] is False
    assert observed["project_to_vsd"] is True
    assert observed["topk"] == 30


def test_gse_runner_fails_if_pairwise_rescue_is_enabled(monkeypatch) -> None:
    monkeypatch.setattr(
        runner,
        "trace_network_expression",
        lambda expression, **kwargs: {
            "results": [],
            "meta": {"pairwise_rescue": {"enabled": True, "switched": False}},
        },
    )

    with pytest.raises(AssertionError, match="pairwise rescue must remain disabled"):
        runner.run_frozen_sample(
            pd.DataFrame({"gene_symbol": ["G1"], "read_count": [1]}),
            min_network_overlap=0.5,
            db_path=Path("unused.db"),
            atlas_id=1,
            topk_regions=30,
        )


def test_overview_uses_tissue_trained_transfer_boundary() -> None:
    source = (Path(__file__).resolve().parents[1] / "app" / "pages" / "overview_page.py").read_text(encoding="utf-8")
    assert "A comprehensive platform for plasma cfRNA tracing" not in source
    assert "macaque plasma cfRNA tracing mode" not in source
    assert "tissue-trained platform for hierarchical brain-origin candidate ranking" in source
    assert "cfRNA use requires separate validation" in source


def test_benchmark_manifest_uses_container_build_provenance_without_git(monkeypatch) -> None:
    def missing_git(*_args, **_kwargs):
        raise FileNotFoundError("git")

    monkeypatch.setenv("BRAINTRACE_GIT_SHA", "0123456789abcdef")
    monkeypatch.setattr(benchmark.subprocess, "run", missing_git)

    assert benchmark.git_commit() == "0123456789abcdef"
