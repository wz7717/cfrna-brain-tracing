from __future__ import annotations

import math

import numpy as np
import pytest

from scripts.run_huang2025_external_candidate import (
    benjamini_hochberg,
    cliffs_delta,
    jaccard,
    parse_sample_id,
    source_log2_rpm_to_log1p_rpm,
)


def test_parse_sample_id_locks_pair_key() -> None:
    glioma = parse_sample_id("GLI_CSF16")
    plasma = parse_sample_id("GLI_plasma16")
    assert glioma["disease"] == "glioma"
    assert glioma["specimen"] == "CSF"
    assert plasma["specimen"] == "plasma"
    assert glioma["patient_key"] == plasma["patient_key"] == "GLI_16"


def test_parse_sample_id_rejects_unknown_contract() -> None:
    with pytest.raises(ValueError, match="Unexpected Huang"):
        parse_sample_id("GBM_CSF1")


def test_log2_rpm_conversion_is_exact_scale_change() -> None:
    source = np.array([0.0, 1.0, 2.0, 3.5])
    converted = source_log2_rpm_to_log1p_rpm(source)
    np.testing.assert_allclose(converted, source * math.log(2.0), rtol=1e-7)


def test_log2_rpm_conversion_rejects_negative_values() -> None:
    with pytest.raises(ValueError, match="negative"):
        source_log2_rpm_to_log1p_rpm(np.array([0.0, -0.1]))


def test_jaccard_and_cliffs_delta() -> None:
    assert jaccard(["A", "B", "C"], ["B", "C", "D"]) == pytest.approx(0.5)
    assert cliffs_delta(np.array([3.0, 4.0]), np.array([1.0, 2.0])) == pytest.approx(1.0)


def test_benjamini_hochberg_is_monotone_by_p_value() -> None:
    p_values = [0.04, 0.001, 0.02]
    adjusted = benjamini_hochberg(p_values)
    ordered = [adjusted[i] for i in np.argsort(p_values)]
    assert ordered == sorted(ordered)
    assert all(p <= q <= 1.0 for p, q in zip(p_values, adjusted))
