from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from scripts.run_p0_4_sparse_domain_shift_sensitivity import (
    EXPECTED_BASELINE,
    METRICS,
    _cluster_bootstrap_draws,
    across_repeat_summary,
    read_detail_csv,
    repeat_summary,
    validate_baseline,
)


def _synthetic_detail() -> pd.DataFrame:
    rows = []
    for scenario in ("baseline", "mild"):
        for repeat in range(2):
            for donor, values in (("A", (1.0, 0.0)), ("B", (0.0,))):
                for sample_i, value in enumerate(values):
                    row = {
                        "scenario": scenario,
                        "replicate": repeat,
                        "sample_id": f"{scenario}-{repeat}-{donor}-{sample_i}",
                        "monkey_id": donor,
                        "depth_fraction": 1.0 if scenario == "baseline" else 0.5,
                        "target_gene_retention": 1.0 if scenario == "baseline" else 0.8,
                    }
                    row.update({metric: value for metric in METRICS})
                    rows.append(row)
    return pd.DataFrame(rows)


def test_repeat_summary_reports_both_estimators() -> None:
    summary = repeat_summary(_synthetic_detail())
    rows = summary[(summary["scenario"] == "baseline") & (summary["replicate"] == 0) & (summary["metric"] == "network_hit1")]
    observed = rows.set_index("estimator")["value"].to_dict()
    assert observed["sample_weighted"] == pytest.approx(1 / 3)
    assert observed["donor_macro"] == pytest.approx(0.25)


def test_vectorized_bootstrap_is_deterministic_and_finite() -> None:
    detail = _synthetic_detail()
    first = _cluster_bootstrap_draws(detail[detail["scenario"] == "mild"], "network_hit1", 1000, 20260716)
    second = _cluster_bootstrap_draws(detail[detail["scenario"] == "mild"], "network_hit1", 1000, 20260716)
    assert set(first) == {"sample_weighted", "donor_macro"}
    assert all(draw.shape == (1000,) for draw in first.values())
    assert all(np.isfinite(draw).all() for draw in first.values())
    assert all(np.array_equal(first[name], second[name]) for name in first)


def test_across_repeat_summary_has_registered_statistics_without_nan() -> None:
    detail = _synthetic_detail()
    result = across_repeat_summary(detail, repeat_summary(detail), 500, 20260716)
    expected = {"mean", "sd", "min", "max", "mc_q025", "mc_q975", "donor_bootstrap_ci_low", "donor_bootstrap_ci_high"}
    assert expected.issubset(result.columns)
    assert not result[list(expected)].isna().any().any()
    assert set(result["estimator"]) == {"sample_weighted", "donor_macro"}


def test_canonical_baseline_gate_is_fail_closed() -> None:
    rows = []
    for sample_idx in range(819):
        row = {"scenario": "baseline", "monkey_id": f"M{sample_idx % 9}"}
        for metric, (hits, denominator) in EXPECTED_BASELINE.items():
            row[metric] = int(sample_idx < hits) if sample_idx < denominator else np.nan
        rows.append(row)
    detail = pd.DataFrame(rows)
    ontology = {"regions": 110, "networks": 10, "total_labels": 120}
    assert validate_baseline(detail, ontology, True)["passed"] is True
    with pytest.raises(RuntimeError, match="requires all 819"):
        validate_baseline(detail, ontology, False)
    detail.loc[0, "network_hit1"] = 0
    with pytest.raises(RuntimeError, match="network_hit1"):
        validate_baseline(detail, ontology, True)


def test_detail_reader_preserves_numeric_looking_donor_ids_as_strings(tmp_path) -> None:
    detail_path = tmp_path / "detail.csv"
    detail_path.write_text("sample_id,monkey_id,scenario,replicate\na,2,baseline,0\nb,qbt,baseline,0\nc,2,mild,0\n", encoding="utf-8")
    detail = read_detail_csv(detail_path)
    assert str(detail["monkey_id"].dtype) == "string"
    assert detail["monkey_id"].tolist() == ["2", "qbt", "2"]
    assert detail["monkey_id"].nunique() == 2
