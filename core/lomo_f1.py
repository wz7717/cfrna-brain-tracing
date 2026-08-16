"""Prediction-level LOMO Network F1 provenance and reporting helpers.

The formal LOMO Network endpoint is defined by the frozen, prediction-level
source rather than by rounded class-level summaries.  This module keeps the
source selection and all integer-count metric calculations in one place so
the evidence CSVs, audit scripts, and tests cannot silently drift apart.
"""

from __future__ import annotations

import csv
import hashlib
import math
from pathlib import Path
from typing import Iterable, Mapping


FORMAL_ROUTE_FAMILY = "hybrid_projected_network_logcpm_exact"
FORMAL_ROUTE = "network_discriminative_correlation_top200"
FORMAL_N = 819
FORMAL_TOP1 = 455
FORMAL_TOP3 = 750

CANONICAL_FORMAL_PATH = (
    Path(__file__).resolve().parents[1]
    / "reproducibility"
    / "p2_publication_completeness"
    / "formal_lomo_network_detail.csv"
)


def sha256_file(path: Path) -> str:
    """Return the uppercase SHA256 digest of *path*."""

    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest().upper()


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
    """Load and validate the prediction-level formal LOMO Network source.

    The historical full source contains several route families.  Only the
    exact frozen formal route is selected.  The returned records use a stable
    schema and recompute hit indicators from truth/prediction labels.
    """

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
        sample_id = _first_present(raw, "sample_id")
        monkey_id = _first_present(raw, "monkey_id")
        selected.append(
            {
                "route": route,
                "sample_id": sample_id,
                "monkey_id": monkey_id,
                "truth": truth,
                "pred_top1": pred_top1,
                "pred_top2": pred_top2,
                "pred_top3": pred_top3,
                "hit1": "1" if pred_top1 == truth else "0",
                "hit3": "1"
                if truth in (pred_top1, pred_top2, pred_top3)
                else "0",
                "route_family": route_family,
            }
        )

    if len(selected) != FORMAL_N:
        raise ValueError(
            f"Formal route selection must contain {FORMAL_N} rows, got {len(selected)}"
        )
    sample_ids = [row["sample_id"] for row in selected]
    if len(set(sample_ids)) != FORMAL_N:
        raise ValueError("Formal prediction source contains duplicate sample_id values")
    if sum(row["hit1"] == "1" for row in selected) != FORMAL_TOP1:
        raise ValueError("Formal LOMO Network Top1 does not equal 455/819")
    if sum(row["hit3"] == "1" for row in selected) != FORMAL_TOP3:
        raise ValueError("Formal LOMO Network Top3 does not equal 750/819")
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
    fraction = position - lower
    return ordered[lower] + fraction * (ordered[upper] - ordered[lower])


def _class_sort_key(row: Mapping[str, object]) -> tuple[int, str]:
    return (-int(row["support"]), str(row["class"]))


def compute_lomo_network_metrics(rows: list[Mapping[str, str]]) -> dict[str, object]:
    """Compute class and summary metrics from prediction-level rows.

    All confusion counts are integer TP/FP/FN values.  For a single-label
    multiclass endpoint, micro-F1 is therefore the exact Top1 accuracy.
    """

    if len(rows) != FORMAL_N:
        raise ValueError(f"Expected {FORMAL_N} rows, got {len(rows)}")
    classes = sorted({str(row["truth"]) for row in rows})
    class_rows: list[dict[str, object]] = []
    for cls in classes:
        tp = sum(row["truth"] == cls and row["pred_top1"] == cls for row in rows)
        fp = sum(row["truth"] != cls and row["pred_top1"] == cls for row in rows)
        fn = sum(row["truth"] == cls and row["pred_top1"] != cls for row in rows)
        support = sum(row["truth"] == cls for row in rows)
        predicted_positive = tp + fp
        precision = tp / predicted_positive if predicted_positive else 0.0
        recall = tp / support if support else 0.0
        f1 = 2 * tp / (2 * tp + fp + fn) if (2 * tp + fp + fn) else 0.0
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
    if total_tp != top1_correct:
        raise AssertionError("Sum of class TP must equal prediction-level Top1 hits")
    if total_support != len(rows):
        raise AssertionError("Class supports must sum to prediction-level n")

    f1_values = [float(row["f1"]) for row in class_rows]
    macro = sum(f1_values) / len(f1_values)
    weighted = sum(
        int(row["support"]) * float(row["f1"]) for row in class_rows
    ) / total_support
    nonzero = [value for value in f1_values if value > 0]
    summary = {
        "endpoint": "LOMO_Network",
        "n_classes": len(class_rows),
        "n_samples": total_support,
        "top1_correct": int(top1_correct),
        "top3_correct": int(top3_correct),
        "top1_accuracy": top1_correct / total_support,
        "top3_accuracy": top3_correct / total_support,
        "macro_f1": macro,
        "sd_class_f1": _sample_sd(f1_values),
        "median_class_f1": _linear_percentile(f1_values, 0.5),
        "q1_class_f1": _linear_percentile(f1_values, 0.25),
        "q3_class_f1": _linear_percentile(f1_values, 0.75),
        "iqr_class_f1": _linear_percentile(f1_values, 0.75)
        - _linear_percentile(f1_values, 0.25),
        "weighted_f1": weighted,
        "micro_f1": total_tp / total_support,
        "n_zero_f1_classes": sum(value == 0 for value in f1_values),
        "fraction_zero_f1_classes": sum(value == 0 for value in f1_values)
        / len(f1_values),
        "conditional_macro_f1_nonzero": sum(nonzero) / len(nonzero)
        if nonzero
        else 0.0,
        "conditional_median_f1_nonzero": _linear_percentile(nonzero, 0.5)
        if nonzero
        else 0.0,
    }
    return {"classes": class_rows, "summary": summary}


def macro_class_rows(metrics: Mapping[str, object]) -> list[dict[str, str]]:
    """Convert computed class metrics to the repository macro-F1 schema."""

    return [
        {
            "endpoint": "LOMO_Network",
            "class": str(row["class"]),
            "n": str(int(row["support"])),
            "precision": f"{float(row['precision']):.10f}",
            "recall": f"{float(row['recall']):.10f}",
            "f1": f"{float(row['f1']):.10f}",
        }
        for row in metrics["classes"]  # type: ignore[index]
    ]


def cross3_summary_row(metrics: Mapping[str, object]) -> dict[str, object]:
    """Convert computed metrics to the CROSS1--5 summary row schema."""

    summary = metrics["summary"]  # type: ignore[index]
    return {
        "endpoint": "LOMO_Network",
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

