from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd
import pytest


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = (
    ROOT
    / "validation_runs"
    / "r08_rf_fair_comparator_20260717"
    / "run_rf_fair_comparator.py"
)
SPEC = importlib.util.spec_from_file_location("r08_rf_fair_comparator", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def test_rf_contract_is_frozen_historical_comparator() -> None:
    contract = MODULE.rf_contract()
    assert contract["feature_selection"] == {
        "class": "SelectKBest",
        "score_func": "f_classif",
        "k": 500,
        "fit_scope": "training fold only",
    }
    assert contract["random_forest"] == {
        "n_estimators": 300,
        "class_weight": "balanced_subsample",
        "min_samples_leaf": 2,
        "random_state": 20260629,
        "n_jobs": -1,
        "remaining_parameters": "sklearn constructor defaults",
    }


def test_ranked_top3_uses_fitted_rf_classes() -> None:
    pipe = MODULE.rf_pipeline()
    x = np.asarray(
        [[float(i + j) for j in range(500)] for i in range(12)],
        dtype=np.float32,
    )
    y = np.asarray(["a", "b", "c"] * 4)
    pipe.fit(x, y)
    ranked = MODULE.ranked_top3(pipe, x[0])
    assert len(ranked) == 3
    assert set(ranked) == {"a", "b", "c"}


def test_group_parameters_are_prespecified() -> None:
    assert MODULE.GROUP_PARAMETERS == {
        "top_n_genes": 200,
        "min_resolution_samples": 8,
        "min_merge_samples": 3,
        "min_pair_errors": 2,
        "min_confusion_rate": 0.20,
        "similarity_threshold": 0.95,
        "merge_similarity_threshold": 0.90,
        "max_group_size": 4,
    }


def test_full_inference_constants_are_frozen() -> None:
    assert MODULE.BOOTSTRAP_REPS == 50_000
    assert MODULE.BOOTSTRAP_SEED == 20260717
    assert MODULE.EXPECTED_DONORS == 9
    assert MODULE.EXPECTED_INPUT_SHA256["vsd"] == (
        "286aeab66b21b7fa012fac8ceaa24497894327e0736f9f6b200334c57089a1b3"
    )
    assert set(MODULE.EXPECTED_INPUT_SHA256) == {
        "vsd",
        "network_detail",
        "exact_detail",
        "group_detail",
        "reference_matrix",
        "historical_rf_code",
        "group_code",
        "group_gene_code",
    }


def test_exact_sign_flip_enumerates_all_512_combinations() -> None:
    differences = np.asarray([0.01 * i for i in range(1, 10)], dtype=float)
    raw_p, extreme, observed = MODULE.exact_sign_flip(differences)
    assert observed == pytest.approx(differences.mean())
    assert raw_p == pytest.approx(extreme / 512)
    assert 1 <= extreme <= 512


def test_bh_is_applied_to_three_prespecified_values() -> None:
    adjusted = MODULE.bh_adjust([0.01, 0.04, 0.03])
    assert adjusted == pytest.approx([0.03, 0.04, 0.04])


def test_bootstrap_summary_requires_all_nine_donors() -> None:
    donor = pd.DataFrame(
        {
            "numerator": np.arange(1, 10),
            "denominator": np.repeat(10, 9),
        }
    )
    summary = MODULE.bootstrap_summary(
        donor, np.random.default_rng(20260717), 100
    )
    assert summary["sample_weighted"] == pytest.approx(0.5)
    assert summary["donor_macro"] == pytest.approx(0.5)
    with pytest.raises(ValueError, match="nine donors"):
        MODULE.bootstrap_summary(
            donor.iloc[:-1], np.random.default_rng(20260717), 100
        )


def test_group_label_audit_marks_nonidentical_class_universes() -> None:
    sample_ids = [f"s{i}" for i in range(812)]
    donors = [str(i % 9) for i in range(812)]
    predictions = pd.DataFrame(
        {
            "endpoint": ["Group"] * 812,
            "sample_id": sample_ids,
            "heldout_monkey_id": donors,
            "truth": ["fold_group"] * 812,
        }
    )
    formal = pd.DataFrame(
        {
            "sample_id": sample_ids,
            "monkey_id": donors,
            "true_resolution_group": ["fold_group"] * 132
            + ["beam_group"] * 680,
        }
    )
    table, payload = MODULE.audit_group_label_universe(predictions, formal)
    assert len(table) == 10
    assert payload["status"] == "class_universes_not_identical"
    assert payload["matching_group_label_strings"] == 132
    assert payload["mismatching_group_label_strings"] == 680
    assert "cannot support superiority" in payload["interpretation"]
