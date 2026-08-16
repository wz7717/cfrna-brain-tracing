from __future__ import annotations

import csv
import math
from pathlib import Path

from core.lomo_f1 import (
    CANONICAL_FORMAL_PATH,
    FORMAL_ROUTE_FAMILY,
    compute_lomo_network_metrics,
    load_formal_predictions,
)


ROOT = Path(__file__).resolve().parents[1]


def test_formal_prediction_source_and_integer_accounting() -> None:
    rows = load_formal_predictions(CANONICAL_FORMAL_PATH)
    metrics = compute_lomo_network_metrics(rows)
    summary = metrics["summary"]

    assert len(rows) == 819
    assert {row["route_family"] for row in rows} == {FORMAL_ROUTE_FAMILY}
    assert summary["top1_correct"] == 455
    assert summary["top3_correct"] == 750
    assert sum(int(row["tp"]) for row in metrics["classes"]) == 455
    assert math.isclose(summary["micro_f1"], 455 / 819, rel_tol=0, abs_tol=1e-15)
    assert summary["n_zero_f1_classes"] == 0


def test_formal_class_counts_match_frozen_top1() -> None:
    metrics = compute_lomo_network_metrics(load_formal_predictions(CANONICAL_FORMAL_PATH))
    by_class = {row["class"]: row for row in metrics["classes"]}

    assert by_class["Subcortical"]["support"] == 54
    assert by_class["Subcortical"]["tp"] == 41
    assert by_class["Subcortical"]["fp"] == 0
    assert by_class["Subcortical"]["fn"] == 13
    assert by_class["Hippocampal formation"]["support"] == 8
    assert by_class["Hippocampal formation"]["tp"] == 8
    assert math.isclose(
        by_class["Temporal"]["f1"], 0.5443425076452599, rel_tol=0, abs_tol=1e-12
    )


def test_rf_comparator_is_not_overwritten_by_formal_source() -> None:
    comparator = ROOT / "reproducibility" / "p2_publication_completeness" / "P2_RF200_lomo_detail.csv"
    with comparator.open(newline="", encoding="utf-8-sig") as handle:
        rows = list(csv.DictReader(handle))
    assert len(rows) == 819
    assert sum(str(row["hit1"]) == "1" for row in rows) == 389
    assert sum(str(row["hit3"]) == "1" for row in rows) == 680


def test_derived_csv_summaries_match_prediction_level_metrics() -> None:
    metrics = compute_lomo_network_metrics(load_formal_predictions(CANONICAL_FORMAL_PATH))
    summary = metrics["summary"]

    with (ROOT / "reproducibility" / "p1_cross1_5" / "cross3_f1_distribution_summary.csv").open(
        newline="", encoding="utf-8-sig"
    ) as handle:
        cross3 = next(row for row in csv.DictReader(handle) if row["endpoint"] == "LOMO_Network")
    assert math.isclose(float(cross3["macro_f1"]), summary["macro_f1"], abs_tol=1e-15)
    assert math.isclose(float(cross3["weighted_f1"]), summary["weighted_f1"], abs_tol=1e-15)
    assert math.isclose(float(cross3["micro_f1"]), 455 / 819, abs_tol=1e-15)

    with (ROOT / "reproducibility" / "v4_p0_13_macro_f1.csv").open(
        newline="", encoding="utf-8-sig"
    ) as handle:
        class_rows = {
            row["class"]: row
            for row in csv.DictReader(handle)
            if row["endpoint"] == "LOMO_Network"
        }
    for expected in metrics["classes"]:
        actual = class_rows[expected["class"]]
        assert int(actual["n"]) == expected["support"]
        assert math.isclose(float(actual["f1"]), expected["f1"], abs_tol=1e-10)
