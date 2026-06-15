from __future__ import annotations

import numpy as np

from core.models import apply_value_transform
from scripts.run_bo2023_v2_loso_validation import InMemoryVsdFoldEngine


def test_vsd_values_are_not_transformed_again():
    values = np.asarray([[0.8, 4.2], [2.0, 7.1]], dtype=float)
    np.testing.assert_allclose(apply_value_transform(values, "vsd"), values)


def test_vsd_engine_uses_only_available_adapted_signals_without_marker_weights():
    genes = np.asarray(["G1", "G2", "G3", "G4"], dtype=object)
    regions = ["R1", "R2", "R3"]
    reference = np.asarray(
        [
            [8.0, 2.0, 1.0],
            [7.0, 1.0, 3.0],
            [1.0, 8.0, 4.0],
            [2.0, 7.0, 8.0],
        ]
    )
    engine = InMemoryVsdFoldEngine(genes, regions, reference, "S1", reference[:, 0])

    result = engine.trace(
        sample_id="S1",
        method="ensemble",
        atlas_id=4,
        use_value="vsd",
        vsd_compatible=True,
        return_all=True,
        persist=False,
        bootstrap_n=0,
        min_overlap_genes=1,
        min_overlap_fraction=0.0,
    )

    weights = result["meta"]["signal_component_weights"]
    assert set(weights) == {"corr", "nnls", "rank"}
    assert np.isclose(weights["corr"], 0.3125)
    assert np.isclose(weights["nnls"], 0.3125)
    assert np.isclose(weights["rank"], 0.375)
    assert result["meta"]["use_value"] == "vsd"
    assert result["results"][0]["region_id"] == "R1"


def test_vsd_fold_marker_weights_enable_marker_and_support_but_not_detection():
    genes = np.asarray(["G1", "G2", "G3", "G4"], dtype=object)
    regions = ["R1", "R2", "R3"]
    reference = np.asarray(
        [
            [8.0, 2.0, 1.0],
            [7.0, 1.0, 3.0],
            [1.0, 8.0, 4.0],
            [2.0, 7.0, 8.0],
        ]
    )
    marker_weights = np.ones_like(reference)
    marker_weights[0, 0] = 2.5
    marker_weights[2, 1] = 2.5
    marker_weights[3, 2] = 2.5
    engine = InMemoryVsdFoldEngine(
        genes, regions, reference, "S1", reference[:, 0], marker_weights=marker_weights
    )

    result = engine.trace(
        sample_id="S1",
        method="ensemble",
        atlas_id=4,
        use_value="vsd",
        vsd_compatible=True,
        return_all=True,
        persist=False,
        bootstrap_n=0,
        min_overlap_genes=1,
        min_overlap_fraction=0.0,
    )

    weights = result["meta"]["signal_component_weights"]
    assert "marker" in weights
    assert "support" in weights
    assert "detect" not in weights
    assert result["meta"]["reference_weighting"]["weighting_active"]
