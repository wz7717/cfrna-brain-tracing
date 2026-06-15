from __future__ import annotations

import numpy as np
import pandas as pd

from scripts.run_bo2023_pairwise_marker_annotation_validation import (
    apply_annotation,
    build_pairwise_marker_signature,
    calibrate_annotation_thresholds,
    pairwise_annotation_score,
)


def test_pairwise_markers_keep_directional_discriminators_only():
    values = np.asarray(
        [
            [8.0, 8.2, 7.9, 1.0, 1.2, 0.9],
            [1.0, 1.2, 0.9, 8.0, 8.1, 7.8],
            [3.0, 8.0, 1.0, 3.0, 3.0, 3.0],
        ]
    )
    markers = build_pairwise_marker_signature(
        values,
        ["R1_up", "R2_up", "unstable"],
        {"R1": np.asarray([0, 1, 2]), "R2": np.asarray([3, 4, 5])},
        {"R1": ["R2"], "R2": ["R1"]},
        max_markers_per_pair=2,
        max_markers_per_region=2,
        min_effect=0.5,
        min_consistency=0.75,
    )
    assert set(markers.loc[markers["region_id"] == "R1", "gene_symbol"]) == {"R1_up"}
    assert set(markers.loc[markers["region_id"] == "R2", "gene_symbol"]) == {"R2_up"}
    score, n_pairs, n_markers = pairwise_annotation_score(
        values[:, 0], "R1", ["R2"], markers, min_markers_per_pair=1
    )
    assert score > 0
    assert n_pairs == 1
    assert n_markers == 1


def test_calibrated_annotation_selects_high_and_low_evidence_groups():
    calibration = pd.DataFrame(
        {
            "annotation_score": [-3, -2, -1, -0.5, 0.1, 0.2, 1.0, 2.0, 3.0, 4.0],
            "hit1": [0, 0, 0, 0, 0, 1, 1, 1, 1, 1],
        }
    )
    thresholds = calibrate_annotation_thresholds(calibration, min_group_size=3)
    assert apply_annotation(4.0, thresholds) == "supported"
    assert apply_annotation(-3.0, thresholds) == "contradicted"
    assert apply_annotation(float("nan"), thresholds) == "insufficient_marker_evidence"
