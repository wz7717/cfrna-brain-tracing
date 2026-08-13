from __future__ import annotations

import pandas as pd
import pytest

from app.pages import tracing_page
from core import inference
from core import production_route
from core import bo2023_region_tracing
from core.model_lock import ModelLockError


def _must_not_run(*args, **kwargs):
    raise AssertionError("inference ran after the production lock gate failed")


def test_production_entry_fails_closed_on_artifact_lock_failure(monkeypatch) -> None:
    monkeypatch.setattr(
        production_route,
        "verify_locked_model_bundle",
        lambda: (_ for _ in ()).throw(ModelLockError("locked artifact drift")),
    )
    monkeypatch.setattr(inference.network_tracing, "trace_network_expression", _must_not_run)
    monkeypatch.setattr(inference.region_tracing, "trace_bo2023_secondary_regions", _must_not_run)

    with pytest.raises(ModelLockError, match="locked artifact drift"):
        tracing_page._run_locked_bo2023_route(pd.DataFrame(), atlas_id=None)


def test_production_entry_fails_closed_on_implementation_drift(monkeypatch) -> None:
    monkeypatch.setattr(bo2023_region_tracing, "EXACT_TOP50_GENE_COUNT", 51)
    monkeypatch.setattr(inference.network_tracing, "trace_network_expression", _must_not_run)
    monkeypatch.setattr(inference.region_tracing, "trace_bo2023_secondary_regions", _must_not_run)

    with pytest.raises(ModelLockError, match="exact_top50_gene_count"):
        tracing_page._run_locked_bo2023_route(pd.DataFrame(), atlas_id=None)


def test_production_entry_passes_verified_region_constants(monkeypatch) -> None:
    captured = {"network": {}, "region": {}}

    def fake_network_route(*args, **kwargs):
        captured["network"].update(kwargs)
        return {
            "results": [{"network_id": "N1"}, {"network_id": "N2"}, {"network_id": "N3"}],
            "meta": {"n_networks": 10, "n_model_genes": 200},
        }

    monkeypatch.setattr(
        inference.network_tracing,
        "trace_network_expression",
        fake_network_route,
    )

    def fake_region_route(*args, **kwargs):
        captured["region"].update(kwargs)
        return {
            "results": [{"region_id": "R1"}],
            "meta": {
                "network_top_k": kwargs["network_top_k"],
                "n_local_candidate_genes": kwargs["local_top_n_genes"],
                "n_scoring_genes_top50": kwargs["exact_top50_gene_count"],
                "n_scoring_genes_top100": kwargs["exact_top100_gene_count"],
                "min_required_region_overlap_genes": kwargs["min_region_gene_overlap"],
            },
        }

    monkeypatch.setattr(inference.region_tracing, "trace_bo2023_secondary_regions", fake_region_route)

    tracing_page._run_locked_bo2023_route(pd.DataFrame(), atlas_id=None)

    expected = production_route.production_implementation_parameters()
    assert captured["network"]["min_overlap_fraction"] == expected["network_min_overlap_fraction"] == 0.50
    assert captured["network"]["project_to_vsd"] is expected["project_to_vsd"] is True
    assert captured["network"]["enable_pairwise_rescue"] is expected["enable_pairwise_rescue"] is False
    assert captured["region"]["network_top_k"] == expected["network_top_k"] == 3
    assert captured["region"]["local_top_n_genes"] == expected["region_local_top_n_genes"] == 200
    assert captured["region"]["exact_top50_gene_count"] == expected["exact_top50_gene_count"] == 50
    assert captured["region"]["exact_top100_gene_count"] == expected["exact_top100_gene_count"] == 100
    assert captured["region"]["top50_weight"] == expected["exact_top50_weight"] == 0.25
    assert captured["region"]["min_region_gene_overlap"] == expected["region_min_overlap_genes"] == 20


def test_production_entry_rejects_region_runtime_metadata_drift(monkeypatch) -> None:
    monkeypatch.setattr(
        inference.network_tracing,
        "trace_network_expression",
        lambda *args, **kwargs: {
            "results": [{"network_id": "N1"}],
            "meta": {"n_networks": 10, "n_model_genes": 200},
        },
    )
    monkeypatch.setattr(
        inference.region_tracing,
        "trace_bo2023_secondary_regions",
        lambda *args, **kwargs: {
            "results": [{"region_id": "R1"}],
            "meta": {
                "network_top_k": 3,
                "n_local_candidate_genes": 200,
                "n_scoring_genes_top50": 49,
                "n_scoring_genes_top100": 100,
                "min_required_region_overlap_genes": 20,
            },
        },
    )

    with pytest.raises(ModelLockError, match="n_scoring_genes_top50"):
        tracing_page._run_locked_bo2023_route(pd.DataFrame(), atlas_id=None)
