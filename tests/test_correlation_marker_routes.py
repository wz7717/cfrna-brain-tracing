from __future__ import annotations

import numpy as np

from scripts.run_bo2023_correlation_marker_routes_validation import (
    annotate_correlation_top1,
    build_stable_marker_signature,
    rerank_correlation_topk,
)


def test_stable_markers_use_consistent_region_specific_signal():
    values = np.asarray(
        [
            [8.0, 8.2, 7.8, 8.1, 1.0, 1.1, 0.9, 1.0],
            [9.0, 1.0, 9.0, 1.0, 1.0, 1.0, 1.0, 1.0],
            [1.0, 1.1, 0.9, 1.0, 7.8, 8.0, 8.2, 8.1],
        ]
    )
    result = build_stable_marker_signature(
        values,
        ["A", "unstable", "B"],
        ["R1", "R2"],
        {"R1": np.asarray([0, 1, 2, 3]), "R2": np.asarray([4, 5, 6, 7])},
        heldout_idx=-1,
        topk_per_region=2,
        min_region_train_samples=3,
        min_consistency=0.75,
        min_effect=0.5,
    )
    assert set(result.loc[result["region_id"] == "R1", "gene_symbol"]) == {"A"}
    assert set(result.loc[result["region_id"] == "R2", "gene_symbol"]) == {"B"}


def test_rerank_only_changes_candidates_and_annotation_does_not_predict():
    correlation = np.asarray([1.0, 0.99, 0.98, 0.2])
    marker = np.asarray([0.0, 4.0, -2.0, 100.0])
    reranked, _, _ = rerank_correlation_topk(correlation, marker, candidate_k=3, marker_weight=0.20)
    assert set(reranked[:3]) == {0, 1, 2}
    assert reranked[-1] == 3
    status, marker_rank, _, _ = annotate_correlation_top1(
        np.argsort(correlation)[::-1], marker, np.asarray([10, 10, 10, 10]), candidate_k=3
    )
    assert status == "supported"
    assert marker_rank == 2
