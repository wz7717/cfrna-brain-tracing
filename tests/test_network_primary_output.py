from __future__ import annotations

import numpy as np
import pandas as pd

from core.network_tracing import (
    DEFAULT_BO2023_NETWORK_METADATA,
    DEFAULT_BO2023_NETWORK_MODEL,
    load_network_model,
    trace_network_expression,
)


def test_network_expression_trace_returns_network_primary_result():
    model = load_network_model(DEFAULT_BO2023_NETWORK_MODEL)
    expression = pd.DataFrame(
        {"gene_symbol": model["genes"], "tpm_value": model["reference"][:, 0]}
    )
    result = trace_network_expression(
        expression, DEFAULT_BO2023_NETWORK_MODEL, DEFAULT_BO2023_NETWORK_METADATA
    )
    assert result["results"][0]["network_id"] == str(model["networks"][0])
    assert result["meta"]["n_model_genes"] == len(model["genes"])


def test_network_expression_trace_rejects_low_gene_overlap():
    expression = pd.DataFrame({"gene_symbol": ["missing"], "tpm_value": [1.0]})
    result = trace_network_expression(expression)
    assert result["results"] == []
    assert result["meta"]["traceability"] == "insufficient"
