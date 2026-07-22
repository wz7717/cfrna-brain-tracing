from __future__ import annotations

from types import MappingProxyType
from typing import Any, Mapping

import core.bo2023_region_tracing as region_tracing
from core.model_lock import (
    assert_locked_production_parameters,
    verify_locked_model_bundle,
)


# This declaration is intentionally independent of the lock manifest. Production
# uses these values, then verifies that they still match the frozen lock before
# loading any model artifact or running inference.
def production_implementation_parameters() -> Mapping[str, Any]:
    """Declare current production values from the constants used by inference."""
    return MappingProxyType(
        {
            "route_name": region_tracing.ROUTE_NAME,
            "region_count": region_tracing.CANONICAL_REGION_COUNT,
            "network_count": region_tracing.CANONICAL_NETWORK_COUNT,
            "beam_count": region_tracing.CANONICAL_BEAM_COUNT,
            "network_top_k": region_tracing.NETWORK_TOP_K,
            "network_gene_count": 200,
            "network_min_overlap_fraction": 0.50,
            "project_to_vsd": True,
            "enable_pairwise_rescue": False,
            "region_local_top_n_genes": region_tracing.DEFAULT_LOCAL_TOP_N_GENES,
            "region_min_overlap_genes": region_tracing.MIN_REGION_GENE_OVERLAP,
            "exact_top50_gene_count": region_tracing.EXACT_TOP50_GENE_COUNT,
            "exact_top100_gene_count": region_tracing.EXACT_TOP100_GENE_COUNT,
            "exact_top50_weight": region_tracing.DEFAULT_TOP50_WEIGHT,
            # Retained only to preserve the frozen manifest contract. The old
            # best-plus-mean group aggregation is inactive: production ranks
            # resolution groups independently with the locked local Top200 panel.
            "resolution_group_mean_weight": 0.10,
            "allow_development_fallback": region_tracing.ALLOW_DEVELOPMENT_FALLBACK,
        }
    )

INACTIVE_LEGACY_PARAMETER_KEYS = frozenset({"resolution_group_mean_weight"})


def verify_production_route() -> tuple[dict[str, Any], dict[str, Any]]:
    """Verify implementation parameters and artifacts before production inference."""
    parameters = dict(production_implementation_parameters())
    assert_locked_production_parameters(parameters)
    manifest = verify_locked_model_bundle()
    return parameters, manifest
