from __future__ import annotations

import pandas as pd
import pytest

from app.shared import DB_PATH
import core.bo2023_region_tracing as bo2023_region_tracing
from core.bo2023_region_tracing import (
    DEFAULT_BO2023_COUNTS,
    DEFAULT_BO2023_GENE_MAP,
    DEFAULT_BO2023_SAMPLE_INFO,
    ROUTE_NAME,
    trace_bo2023_secondary_regions,
)
from core.network_tracing import load_network_model


ROOT = DEFAULT_BO2023_COUNTS.parents[1]
STAGED_BO2023_INPUTS = (
    ROOT / "tests" / "controlled_data" / "bo2023" / "mfas5_819samples_28415genes_featurecounts_counts.txt",
    ROOT / "tests" / "controlled_data" / "bo2023" / "Information of sequenced samples_update_full878_filter819.xlsx",
    ROOT
    / "tests"
    / "controlled_data"
    / "bo2023"
    / "04_expressed_genes_neocortex_plus_subcortical.cleaned_symbols.csv",
)
DEFAULT_CONTROLLED_BO2023_INPUTS = (
    DEFAULT_BO2023_COUNTS,
    DEFAULT_BO2023_SAMPLE_INFO,
    DEFAULT_BO2023_GENE_MAP,
)


def _controlled_bo2023_inputs_or_skip() -> tuple:
    for paths in (STAGED_BO2023_INPUTS, DEFAULT_CONTROLLED_BO2023_INPUTS):
        if all(path.exists() for path in paths):
            return paths
    missing = [path for path in STAGED_BO2023_INPUTS + DEFAULT_CONTROLLED_BO2023_INPUTS if not path.exists()]
    if missing:
        pytest.skip(
            "Controlled Bo2023 raw expression inputs are not included in the public release: "
            + ", ".join(str(path) for path in missing)
        )
    raise AssertionError("unreachable")


def test_bo2023_secondary_region_tracing_uses_network_top3_beam(monkeypatch):
    counts_path, sample_info_path, gene_map_path = _controlled_bo2023_inputs_or_skip()

    def load_controlled_reference(db_path: str, atlas_id: int):
        return bo2023_region_tracing._load_raw_logcpm_reference_matrix(  # noqa: SLF001
            counts_path,
            sample_info_path,
            gene_map_path,
        )

    monkeypatch.setattr(bo2023_region_tracing, "_load_reference_matrix", load_controlled_reference)

    model = load_network_model()
    expression = pd.DataFrame(
        {
            "gene_symbol": model["genes"],
            "tpm_value": model["reference"][:, 0],
        }
    )
    network_output = {
        "results": [
            {"network_id": str(model["networks"][0]), "rank": 1},
            {"network_id": str(model["networks"][1]), "rank": 2},
            {"network_id": str(model["networks"][2]), "rank": 3},
        ]
    }

    out = trace_bo2023_secondary_regions(expression, network_output, DB_PATH, atlas_id=4, topk=8)

    assert out["meta"]["method"] == ROUTE_NAME
    assert out["meta"]["candidate_region_source"] == "SaleemNetworks Top3 beam"
    assert out["meta"]["n_scoring_genes_top50"] <= 50
    assert out["meta"]["n_scoring_genes_top100"] <= 100
    assert len(out["results"]) == 8
    assert {"region_id", "score", "top50_corr_component", "top100_corr_component"}.issubset(out["results"][0])


def test_bo2023_secondary_region_tracing_rejects_missing_network_beam():
    out = trace_bo2023_secondary_regions(
        pd.DataFrame({"gene_symbol": ["A"], "tpm_value": [1.0]}),
        {"results": []},
        DB_PATH,
        atlas_id=4,
    )

    assert out["results"] == []
    assert out["meta"]["traceability"] == "insufficient"


def test_production_reference_prefers_committed_package_over_local_raw_inputs(monkeypatch):
    packaged = pd.DataFrame({"N1::R1": [1.0]}, index=["GENE1"])
    network_map = {"N1::R1": "N1"}

    monkeypatch.setattr(
        bo2023_region_tracing,
        "_load_packaged_region_reference_matrix",
        lambda *args, **kwargs: (packaged, network_map, "packaged_region_logcpm_reference"),
    )
    monkeypatch.setattr(
        bo2023_region_tracing,
        "_load_raw_logcpm_reference_matrix",
        lambda *args, **kwargs: (_ for _ in ()).throw(
            AssertionError("production inference must not prefer workstation-only raw inputs")
        ),
    )

    matrix, mapping, source = bo2023_region_tracing._load_reference_matrix("unused.db", None)  # noqa: SLF001

    pd.testing.assert_frame_equal(matrix, packaged)
    assert mapping == network_map
    assert source == "packaged_region_logcpm_reference"


def test_formal_group_and_exact_rankings_are_independent():
    rows = [
        {
            "region_id": "R1",
            "resolution_group": "G1",
            "resolution_group_members": "R1",
            "exact_local_score": 3.0,
            "group_local_score": 0.1,
            "resolution_tier": "high_resolution",
            "manual_review_recommended": False,
        },
        {
            "region_id": "R2",
            "resolution_group": "G2",
            "resolution_group_members": "R2",
            "exact_local_score": 1.0,
            "group_local_score": 0.9,
            "resolution_tier": "high_resolution",
            "manual_review_recommended": False,
        },
    ]

    exact = bo2023_region_tracing._rank_exact_regions(rows, topk=2)  # noqa: SLF001
    group_input = sorted(rows, key=lambda row: row["group_local_score"], reverse=True)
    groups = bo2023_region_tracing._rank_resolution_groups(group_input)  # noqa: SLF001

    assert [row["region_id"] for row in exact] == ["R1", "R2"]
    assert [row["resolution_group"] for row in groups] == ["G2", "G1"]
    assert groups[0]["group_score"] == 0.9
