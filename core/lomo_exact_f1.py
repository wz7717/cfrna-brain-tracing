"""Prediction-level provenance helpers for the formal LOMO Exact endpoint.

The class universe for this endpoint is the set of truth labels in the frozen
formal prediction table.  It deliberately does not expand to prediction-only
labels: those predictions are accounted for as false positives against the
truth-label universe, which is the manuscript's stated evaluability endpoint.
"""

from __future__ import annotations

import csv
import math
from pathlib import Path
from typing import Iterable, Mapping


FORMAL_ROUTE_FAMILY = "hybrid_projected_network_logcpm_exact"
FORMAL_ROUTE = "top3_beam_local_top50_top100_zfusion_w0p25"
FORMAL_N = 812
FORMAL_TOP1 = 177
FORMAL_TOP3 = 346
FORMAL_N_CLASSES = 104

CANONICAL_FORMAL_PATH = (
    Path(__file__).resolve().parents[1]
    / "reproducibility"
    / "p2_publication_completeness"
    / "formal_lomo_exact_region_detail.csv"
)


def _first_present(row: Mapping[str, str], *names: str) -> str:
    for name in names:
        value = row.get(name)
        if value is not None and str(value).strip() != "":
            return str(value).strip()
    raise ValueError(f"None of the required columns are present: {names}")


def load_formal_predictions(
    path: Path,
    *,
    expected_route_family: str = FORMAL_ROUTE_FAMILY,
    expected_route: str = FORMAL_ROUTE,
) -> list[dict[str, str]]:
    """Load the locked route and recompute its Top1/Top3 hit indicators."""

    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames:
            raise ValueError(f"CSV has no header: {path}")
        raw_rows = list(reader)

    selected: list[dict[str, str]] = []
    for raw in raw_rows:
        route_family = str(raw.get("route_family", "")).strip()
        if route_family != expected_route_family:
            continue
        route = _first_present(raw, "route")
        if route != expected_route:
            raise ValueError(
                f"Unexpected route for {expected_route_family}: {route!r} in {path}"
            )
        truth = _first_present(raw, "truth", "label")
        pred_top1 = _first_present(raw, "pred_top1")
        pred_top2 = _first_present(raw, "pred_top2")
        pred_top3 = _first_present(raw, "pred_top3")
        selected.append(
            {
                "route": route,
                "sample_id": _first_present(raw, "sample_id"),
                "monkey_id": _first_present(raw, "monkey_id"),
                "truth": truth,
                "pred_top1": pred_top1,
                "pred_top2": pred_top2,
                "pred_top3": pred_top3,
                "hit1": "1" if pred_top1 == truth else "0",
                "hit3": "1" if truth in (pred_top1, pred_top2, pred_top3) else "0",
                "route_family": route_family,
            }
        )

    sample_ids = [row["sample_id"] for row in selected]
    if len(selected) != FORMAL_N or len(set(sample_ids)) != FORMAL_N:
        raise ValueError("Formal LOMO Exact source must contain 812 unique sample_id values")
    if sum(row["hit1"] == "1" for row in selected) != FORMAL_TOP1:
        raise ValueError("Formal LOMO Exact Top1 does not equal 177/812")
    if sum(row["hit3"] == "1" for row in selected) != FORMAL_TOP3:
        raise ValueError("Formal LOMO Exact Top3 does not equal 346/812")
    return selected


def _sample_sd(values: list[float]) -> float:
    if len(values) < 2:
        return 0.0
    mean = sum(values) / len(values)
    return math.sqrt(sum((value - mean) ** 2 for value in values) / (len(values) - 1))


def _linear_percentile(values: Iterable[float], quantile: float) -> float:
    """NumPy-compatible linear percentile without requiring NumPy."""

    ordered = sorted(float(value) for value in values)
    if not ordered:
        return float("nan")
    if len(ordered) == 1:
        return ordered[0]
    position = (len(ordered) - 1) * quantile
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return ordered[lower]
    return ordered[lower] + (position - lower) * (ordered[upper] - ordered[lower])


def _class_sort_key(row: Mapping[str, object]) -> tuple[int, str]:
    return (-int(row["support"]), str(row["class"]))


def compute_lomo_exact_metrics(rows: list[Mapping[str, str]]) -> dict[str, object]:
    """Compute integer class accounting and all F1 summaries from predictions."""

    if len(rows) != FORMAL_N:
        raise ValueError(f"Expected {FORMAL_N} rows, got {len(rows)}")
    classes = sorted({str(row["truth"]) for row in rows})
    if len(classes) != FORMAL_N_CLASSES:
        raise ValueError(f"Expected {FORMAL_N_CLASSES} truth classes, got {len(classes)}")

    class_rows: list[dict[str, object]] = []
    for cls in classes:
        tp = sum(row["truth"] == cls and row["pred_top1"] == cls for row in rows)
        fp = sum(row["truth"] != cls and row["pred_top1"] == cls for row in rows)
        fn = sum(row["truth"] == cls and row["pred_top1"] != cls for row in rows)
        support = sum(row["truth"] == cls for row in rows)
        predicted_positive = tp + fp
        precision = tp / predicted_positive if predicted_positive else 0.0
        recall = tp / support if support else 0.0
        f1 = 2 * tp / (2 * tp + fp + fn) if 2 * tp + fp + fn else 0.0
        class_rows.append(
            {
                "class": cls,
                "support": int(support),
                "tp": int(tp),
                "fp": int(fp),
                "fn": int(fn),
                "predicted_positive": int(predicted_positive),
                "precision": float(precision),
                "recall": float(recall),
                "f1": float(f1),
            }
        )
    class_rows.sort(key=_class_sort_key)

    total_tp = sum(int(row["tp"]) for row in class_rows)
    total_support = sum(int(row["support"]) for row in class_rows)
    top1_correct = sum(row["truth"] == row["pred_top1"] for row in rows)
    top3_correct = sum(
        row["truth"] in (row["pred_top1"], row["pred_top2"], row["pred_top3"])
        for row in rows
    )
    if total_tp != top1_correct or total_support != len(rows):
        raise AssertionError("Class accounting does not equal prediction-level accounting")

    f1_values = [float(row["f1"]) for row in class_rows]
    nonzero = [value for value in f1_values if value > 0]
    truth_universe = set(classes)
    predicted_outside_truth = sorted(
        {str(row["pred_top1"]) for row in rows} - truth_universe
    )
    summary = {
        "endpoint": "LOMO_Exact",
        "n_classes": len(class_rows),
        "n_samples": total_support,
        "top1_correct": int(top1_correct),
        "top3_correct": int(top3_correct),
        "top1_accuracy": top1_correct / total_support,
        "top3_accuracy": top3_correct / total_support,
        "macro_f1": sum(f1_values) / len(f1_values),
        "sd_class_f1": _sample_sd(f1_values),
        "median_class_f1": _linear_percentile(f1_values, 0.5),
        "q1_class_f1": _linear_percentile(f1_values, 0.25),
        "q3_class_f1": _linear_percentile(f1_values, 0.75),
        "iqr_class_f1": _linear_percentile(f1_values, 0.75)
        - _linear_percentile(f1_values, 0.25),
        "weighted_f1": sum(
            int(row["support"]) * float(row["f1"]) for row in class_rows
        )
        / total_support,
        "micro_f1": total_tp / total_support,
        "n_zero_f1_classes": sum(value == 0 for value in f1_values),
        "fraction_zero_f1_classes": sum(value == 0 for value in f1_values)
        / len(f1_values),
        "conditional_macro_f1_nonzero": sum(nonzero) / len(nonzero),
        "conditional_median_f1_nonzero": _linear_percentile(nonzero, 0.5),
        "truth_label_universe": "truth labels only; prediction-only labels remain FP",
        "predicted_top1_labels_outside_truth_universe": predicted_outside_truth,
    }
    return {"classes": class_rows, "summary": summary}


def macro_class_rows(metrics: Mapping[str, object]) -> list[dict[str, str]]:
    return [
        {
            "endpoint": "LOMO_Exact",
            "class": str(row["class"]),
            "n": str(int(row["support"])),
            "precision": f"{float(row['precision']):.10f}",
            "recall": f"{float(row['recall']):.10f}",
            "f1": f"{float(row['f1']):.10f}",
        }
        for row in metrics["classes"]  # type: ignore[index]
    ]


def cross3_summary_row(metrics: Mapping[str, object]) -> dict[str, object]:
    summary = metrics["summary"]  # type: ignore[index]
    return {
        "endpoint": "LOMO_Exact",
        "n_classes": int(summary["n_classes"]),
        "n_samples": int(summary["n_samples"]),
        "macro_f1": float(summary["macro_f1"]),
        "sd_class_f1": float(summary["sd_class_f1"]),
        "median_class_f1": float(summary["median_class_f1"]),
        "weighted_f1": float(summary["weighted_f1"]),
        "micro_f1": float(summary["micro_f1"]),
        "n_zero_f1_classes": int(summary["n_zero_f1_classes"]),
        "fraction_zero_f1_classes": float(summary["fraction_zero_f1_classes"]),
        "conditional_macro_f1_nonzero": float(
            summary["conditional_macro_f1_nonzero"]
        ),
        "conditional_median_f1_nonzero": float(
            summary["conditional_median_f1_nonzero"]
        ),
    }
