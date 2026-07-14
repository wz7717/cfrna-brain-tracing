from __future__ import annotations

from scripts.run_ahba_projected_vsd_formal_three_tier_external import group_hits
from scripts.run_bo2023_resolution_tier_validation import normalize_resolution_annotations


def test_normalize_resolution_annotations_removes_repeated_network_prefixes() -> None:
    annotations = {
        "Occipital/Temporal::FST": {
            "region_id": "Occipital/Temporal::FST",
            "network_id": "Occipital/Temporal",
            "resolution_group": (
                "Occipital/Temporal::Occipital/Temporal::FST + Occipital/Temporal::MST"
            ),
            "group_members": ["Occipital/Temporal::FST", "Occipital/Temporal::MST"],
        },
        "Occipital/Temporal::V1": {
            "region_id": "Occipital/Temporal::V1",
            "network_id": "Occipital/Temporal",
            "resolution_group": "Occipital/Temporal::V1",
            "group_members": ["Occipital/Temporal::V1"],
        },
    }

    normalized = normalize_resolution_annotations(annotations)

    assert normalized["Occipital/Temporal::FST"]["resolution_group"] == (
        "Occipital/Temporal::FST + MST"
    )
    assert normalized["Occipital/Temporal::FST"]["group_members"] == ["FST", "MST"]
    assert normalized["Occipital/Temporal::V1"]["resolution_group"] == "V1"
    assert normalized["Occipital/Temporal::V1"]["region_id"] == "V1"


def test_group_hits_keeps_predictions_when_truth_is_outside_candidate_beam() -> None:
    annotations = {
        "N1::R1": {"resolution_group": "N1::R1 + R2"},
        "N1::R2": {"resolution_group": "N1::R1 + R2"},
        "N1::R3": {"resolution_group": "R3"},
    }

    hit1, hit3, ranked_groups, allowed_groups = group_hits(
        ["N1::R1", "N1::R3", "N1::R2"],
        annotations,
        ["N2::TRUTH"],
    )

    assert hit1 is False
    assert hit3 is False
    assert ranked_groups == ["N1::R1 + R2", "R3"]
    assert allowed_groups == []


def test_group_hits_returns_predictions_for_non_exact_mapped_truth() -> None:
    annotations = {"N1::R1": {"resolution_group": "R1"}}

    hit1, hit3, ranked_groups, allowed_groups = group_hits(
        ["N1::R1"],
        annotations,
        [],
    )

    assert hit1 is None
    assert hit3 is None
    assert ranked_groups == ["R1"]
    assert allowed_groups == []
