from __future__ import annotations

import json

from core.region_resolution import annotate_region_candidates


def test_region_resolution_annotation_flags_low_resolution_without_reranking(tmp_path):
    model_path = tmp_path / "resolution.json"
    model_path.write_text(
        json.dumps(
            {
                "entries": {
                    "Temporal||TEa": {
                        "resolution_tier": "low_resolution",
                        "resolution_group": "Temporal::TEa + TEm",
                        "group_members": ["TEa", "TEm"],
                        "resolution_reasons": ["fold_local_merged_confusion_group"],
                    }
                }
            }
        ),
        encoding="utf-8",
    )
    region_output = {
        "results": [
            {"region_id": "TEa", "rank": 1, "score": 0.7},
            {"region_id": "TEm", "rank": 2, "score": 0.6},
        ],
        "meta": {},
    }
    network_output = {"results": [{"network_id": "Temporal", "rank": 1}]}

    out = annotate_region_candidates(region_output, network_output, model_path)

    assert [row["region_id"] for row in out["results"]] == ["TEa", "TEm"]
    assert out["results"][0]["resolution_tier"] == "low_resolution"
    assert out["results"][0]["manual_review_recommended"]
    assert out["meta"]["region_resolution_annotation"]["top1_group_members"] == "TEa | TEm"


def test_region_outside_primary_network_is_low_resolution(tmp_path):
    model_path = tmp_path / "resolution.json"
    model_path.write_text(json.dumps({"entries": {}}), encoding="utf-8")
    out = annotate_region_candidates(
        {"results": [{"region_id": "R1", "rank": 1}], "meta": {}},
        {"results": [{"network_id": "N1", "rank": 1}]},
        model_path,
    )

    assert out["results"][0]["resolution_tier"] == "low_resolution"
    assert out["results"][0]["resolution_reasons"] == "outside_primary_network_resolution_map"


def test_hippocampal_single_region_annotation_survives_missing_model(tmp_path):
    model_path = tmp_path / "missing-resolution-model.json"
    out = annotate_region_candidates(
        {"results": [{"region_id": "Hippocampal formation", "rank": 1}], "meta": {}},
        {"results": [{"network_id": "Hippocampal formation", "rank": 1}]},
        model_path,
    )

    assert out["results"][0]["single_region_network"] is True
    assert out["results"][0]["tier_degeneracy"] == "group=exact"
    assert out["meta"]["single_region_network"]["annotation"] == "single-region Network"
