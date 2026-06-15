from __future__ import annotations

import pandas as pd

from app.shared import DB_PATH
from core.bo2023_region_tracing import ROUTE_NAME, trace_bo2023_secondary_regions
from core.network_tracing import load_network_model


def test_bo2023_secondary_region_tracing_uses_network_top3_beam():
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
