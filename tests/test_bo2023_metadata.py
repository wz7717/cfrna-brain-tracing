from __future__ import annotations

import pandas as pd
import pytest

from core.bo2023_metadata import (
    assert_unique_region_network_mapping,
    normalize_bo2023_network_labels,
)


def test_normalizes_the_two_cross_network_region_labels_without_mutating_source() -> None:
    source = pd.DataFrame(
        {
            "sample_id": ["a", "b", "c", "d"],
            "region_id": ["10m", "10m", "V2", "V2"],
            "network_id": [
                "Lateral Prefrontal Cortex",
                "Orbitomedial Prefrontal Cortex (OMPFC)",
                "Parietal, and Parieto-occipital region",
                "Occipital/Temporal",
            ],
        }
    )

    normalized = normalize_bo2023_network_labels(source)

    assert source.loc[0, "network_id"] == "Lateral Prefrontal Cortex"
    assert normalized.loc[normalized["region_id"].eq("10m"), "network_id"].unique().tolist() == [
        "Orbitomedial Prefrontal Cortex (OMPFC)"
    ]
    assert normalized.loc[normalized["region_id"].eq("V2"), "network_id"].unique().tolist() == [
        "Occipital/Temporal"
    ]
    assert_unique_region_network_mapping(normalized)


def test_unique_mapping_guard_rejects_unresolved_cross_network_region() -> None:
    metadata = pd.DataFrame(
        {
            "region_id": ["X", "X"],
            "network_id": ["Network A", "Network B"],
        }
    )

    with pytest.raises(ValueError, match="X"):
        assert_unique_region_network_mapping(metadata)
