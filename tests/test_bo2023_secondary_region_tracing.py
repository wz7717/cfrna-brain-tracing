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
