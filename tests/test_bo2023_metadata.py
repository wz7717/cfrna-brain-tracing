from __future__ import annotations

import pandas as pd
import pytest

from core.bo2023_metadata import (
    assert_unique_region_network_mapping,
    canonicalize_bo2023_region_text,
    normalize_bo2023_region_labels,
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


def test_canonicalizes_source_region_artifacts_before_lobe_lookup() -> None:
    assert canonicalize_bo2023_region_text("AI") == "A1"
    assert canonicalize_bo2023_region_text("44563") == "1-2"
    assert canonicalize_bo2023_region_text("44563 | 7op") == "1-2 | 7op"


def test_region_normalization_does_not_touch_non_region_fields() -> None:
    frame = pd.DataFrame({"region_id": ["AI", "44563"], "sample_id": ["AI", "44563"]})
    normalized = normalize_bo2023_region_labels(frame, columns=["region_id"])
    assert normalized["region_id"].tolist() == ["A1", "1-2"]
    assert normalized["sample_id"].tolist() == ["AI", "44563"]
