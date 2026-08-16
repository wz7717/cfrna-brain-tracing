#!/usr/bin/env python
"""Stage and regenerate the formal prediction-level LOMO Network evidence."""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from core.lomo_f1 import (  # noqa: E402
    CANONICAL_FORMAL_PATH,
    FORMAL_ROUTE,
    FORMAL_ROUTE_FAMILY,
    compute_lomo_network_metrics,
    cross3_summary_row,
    load_formal_predictions,
    macro_class_rows,
    sha256_file,
)


MACRO_JSON = ROOT / "reproducibility" / "macro_f1_class_data.json"
MACRO_CSV = ROOT / "reproducibility" / "v4_p0_13_macro_f1.csv"
CROSS3_CSV = (
    ROOT / "reproducibility" / "p1_cross1_5" / "cross3_f1_distribution_summary.csv"
)
FORMAL_DETAIL = ROOT / "reproducibility" / "formal_lomo_network_f1.csv"
PROVENANCE = ROOT / "reproducibility" / "lomo_network_f1_provenance.json"
PROVENANCE_MD = ROOT / "reproducibility" / "LOMO_NETWORK_F1_PROVENANCE.md"

DETAIL_FIELDS = [
    "class",
    "support",
    "tp",
    "fp",
    "fn",
    "predicted_positive",
    "precision",
    "recall",
    "f1",
]
PREDICTION_FIELDS = [
    "sample_id",
    "monkey_id",
    "truth",
    "pred_top1",
    "pred_top2",
    "pred_top3",
    "hit1",
    "hit3",
    "route",
    "route_family",
]
CROSS3_FIELDS = [
    "endpoint",
    "n_classes",
    "n_samples",
    "macro_f1",
    "sd_class_f1",
    "median_class_f1",
    "weighted_f1",
    "micro_f1",
    "n_zero_f1_classes",
    "fraction_zero_f1_classes",
    "conditional_macro_f1_nonzero",
    "conditional_median_f1_nonzero",
]


def write_csv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def stage_predictions(source: Path) -> list[dict[str, str]]:
    rows = load_formal_predictions(source)
    write_csv(CANONICAL_FORMAL_PATH, rows, PREDICTION_FIELDS)
    return rows


def update_macro_json(metrics: dict[str, object]) -> None:
    payload = json.loads(MACRO_JSON.read_text(encoding="utf-8"))
    class_rows = macro_class_rows(metrics)
    data = [
        row
        for row in payload["data"]
        if row.get("endpoint") not in {"LOMO_Network", "SUMMARY"}
        or not str(row.get("class", "")).startswith("LOMO_Network_")
    ]
    data = [row for row in data if row.get("endpoint") != "LOMO_Network"]
    data.extend(class_rows)
    summary = metrics["summary"]
    data.extend(
        [
            {
                "endpoint": "SUMMARY",
                "class": "LOMO_Network_macro",
                "n": str(summary["n_classes"]),
                "precision": "",
                "recall": "",
                "f1": f"{float(summary['macro_f1']):.10f}",
            },
            {
                "endpoint": "SUMMARY",
                "class": "LOMO_Network_weighted",
                "n": str(summary["n_samples"]),
                "precision": "",
                "recall": "",
                "f1": f"{float(summary['weighted_f1']):.10f}",
            },
        ]
    )
    payload["data"] = data
    payload["provenance"] = {
        "formal_lomo_network_source": str(CANONICAL_FORMAL_PATH.relative_to(ROOT)),
        "route": FORMAL_ROUTE,
        "route_family": FORMAL_ROUTE_FAMILY,
        "metrics_are_prediction_level": True,
        "micro_f1_definition": "sum integer TP / sum support = Top1 accuracy",
    }
    MACRO_JSON.write_text(
        json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )


def update_macro_csv(metrics: dict[str, object]) -> None:
    with MACRO_CSV.open(newline="", encoding="utf-8-sig") as handle:
        existing = list(csv.DictReader(handle))
    class_rows = macro_class_rows(metrics)
    retained = [row for row in existing if row.get("endpoint") != "LOMO_Network"]
    retained = [
        row
        for row in retained
        if not (
            row.get("endpoint") == "SUMMARY"
            and str(row.get("class", "")).startswith("LOMO_Network_")
        )
    ]
    summary = metrics["summary"]
    retained.extend(class_rows)
    retained.extend(
        [
            {
                "endpoint": "SUMMARY",
                "class": "LOMO_Network_macro",
                "n": str(summary["n_classes"]),
                "precision": "",
                "recall": "",
                "f1": f"{float(summary['macro_f1']):.10f}",
            },
            {
                "endpoint": "SUMMARY",
                "class": "LOMO_Network_weighted",
                "n": str(summary["n_samples"]),
                "precision": "",
                "recall": "",
                "f1": f"{float(summary['weighted_f1']):.10f}",
            },
        ]
    )
    write_csv(MACRO_CSV, retained, ["endpoint", "class", "n", "precision", "recall", "f1"])


def update_cross3(metrics: dict[str, object]) -> None:
    with CROSS3_CSV.open(newline="", encoding="utf-8-sig") as handle:
        existing = list(csv.DictReader(handle))
    retained = [row for row in existing if row.get("endpoint") != "LOMO_Network"]
    row = cross3_summary_row(metrics)
    write_csv(CROSS3_CSV, retained + [row], CROSS3_FIELDS)


def write_derived(metrics: dict[str, object]) -> None:
    write_csv(FORMAL_DETAIL, metrics["classes"], DETAIL_FIELDS)  # type: ignore[arg-type]
    update_macro_json(metrics)
    update_macro_csv(metrics)
    update_cross3(metrics)


def verify_derived(metrics: dict[str, object]) -> None:
    """Check that all generated evidence files match the prediction source."""

    with FORMAL_DETAIL.open(newline="", encoding="utf-8-sig") as handle:
        detail_rows = list(csv.DictReader(handle))
    expected_detail = metrics["classes"]
    if len(detail_rows) != len(expected_detail):
        raise ValueError("formal_lomo_network_f1.csv has the wrong row count")
    for actual, expected in zip(detail_rows, expected_detail):
        for key in DETAIL_FIELDS:
            if key in {"precision", "recall", "f1"}:
                if abs(float(actual[key]) - float(expected[key])) > 1e-12:
                    raise ValueError(f"Formal detail mismatch in {key}: {actual[key]}")
            elif str(actual[key]) != str(expected[key]):
                raise ValueError(f"Formal detail mismatch in {key}: {actual[key]}")

    with CROSS3_CSV.open(newline="", encoding="utf-8-sig") as handle:
        cross3_rows = list(csv.DictReader(handle))
    cross3 = next(row for row in cross3_rows if row["endpoint"] == "LOMO_Network")
    expected_cross3 = cross3_summary_row(metrics)
    for key in CROSS3_FIELDS:
        if key == "endpoint":
            continue
        if abs(float(cross3[key]) - float(expected_cross3[key])) > 1e-12:
            raise ValueError(f"CROSS3 mismatch in {key}: {cross3[key]}")


def write_provenance(source: Path, metrics: dict[str, object]) -> None:
    summary = metrics["summary"]
    try:
        staged_from = str(source.resolve().relative_to(ROOT)).replace("\\", "/")
    except ValueError:
        staged_from = f"external_source::{source.name}"
    payload = {
        "schema": "braintrace.lomo_network_f1_provenance.v1",
        "formal_source": {
            "path": str(CANONICAL_FORMAL_PATH.relative_to(ROOT)).replace("\\", "/"),
            "sha256": sha256_file(CANONICAL_FORMAL_PATH),
            "staged_from": staged_from,
            "staged_from_sha256": sha256_file(source),
            "route": FORMAL_ROUTE,
            "route_family": FORMAL_ROUTE_FAMILY,
            "n": summary["n_samples"],
        },
        "prediction_level_metrics": {
            "top1_correct": summary["top1_correct"],
            "top1_accuracy": summary["top1_accuracy"],
            "top3_correct": summary["top3_correct"],
            "top3_accuracy": summary["top3_accuracy"],
            "micro_f1": summary["micro_f1"],
            "macro_f1": summary["macro_f1"],
            "sd_class_f1_sample": summary["sd_class_f1"],
            "median_class_f1": summary["median_class_f1"],
            "q1_class_f1": summary["q1_class_f1"],
            "q3_class_f1": summary["q3_class_f1"],
            "iqr_class_f1": summary["iqr_class_f1"],
            "weighted_f1": summary["weighted_f1"],
            "n_zero_f1_classes": summary["n_zero_f1_classes"],
            "fraction_zero_f1_classes": summary["fraction_zero_f1_classes"],
            "conditional_macro_f1_nonzero": summary["conditional_macro_f1_nonzero"],
            "conditional_median_f1_nonzero": summary["conditional_median_f1_nonzero"],
        },
        "integer_confusion_accounting": {
            "sum_tp": sum(int(row["tp"]) for row in metrics["classes"]),  # type: ignore[index]
            "sum_support": sum(int(row["support"]) for row in metrics["classes"]),  # type: ignore[index]
            "micro_f1_equals_top1_accuracy": True,
        },
        "root_cause": (
            "Historical macro_f1_class_data.json used rounded class recall values and "
            "a separate/stale LOMO route; the old script estimated micro-F1 as "
            "sum(recall*n)/819 rather than counting prediction-level Top1 hits."
        ),
        "superseded_route": {
            "route": "network_pairwise_correlation_rescue_top3",
            "status": "historical; not the formal v0.1.15 reporting endpoint",
        },
    }
    PROVENANCE.write_text(
        json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )

    lines = [
        "# Formal LOMO Network F1 provenance",
        "",
        "This report is generated from the prediction-level formal LOMO Network source.",
        "It is a reporting/provenance correction; the frozen model, ontology, route and",
        "prediction set are unchanged.",
        "",
        f"- Route: `{FORMAL_ROUTE}`",
        f"- Route family: `{FORMAL_ROUTE_FAMILY}`",
        f"- Canonical source: `{CANONICAL_FORMAL_PATH.relative_to(ROOT).as_posix()}`",
        f"- Canonical source SHA-256: `{sha256_file(CANONICAL_FORMAL_PATH)}`",
        f"- Prediction rows: `{summary['n_samples']}`",
        f"- Top1: `{summary['top1_correct']}/{summary['n_samples']} = {summary['top1_accuracy']:.10f}`",
        f"- Top3: `{summary['top3_correct']}/{summary['n_samples']} = {summary['top3_accuracy']:.10f}`",
        "",
        "## Derived metrics",
        "",
        "| class | support | TP | FP | FN | precision | recall | F1 |",
        "|---|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in metrics["classes"]:  # type: ignore[union-attr]
        lines.append(
            f"| {row['class']} | {row['support']} | {row['tp']} | {row['fp']} | "
            f"{row['fn']} | {row['precision']:.10f} | {row['recall']:.10f} | {row['f1']:.10f} |"
        )
    lines.extend(
        [
            "",
            f"- Macro-F1: `{summary['macro_f1']:.10f}`",
            f"- Class-level sample SD: `{summary['sd_class_f1']:.10f}`",
            f"- Median: `{summary['median_class_f1']:.10f}`",
            f"- IQR: `{summary['iqr_class_f1']:.10f}`",
            f"- Weighted-F1: `{summary['weighted_f1']:.10f}`",
            f"- Micro-F1: `{summary['micro_f1']:.10f}` (= Top1 accuracy)",
            f"- Zero-F1 classes: `{summary['n_zero_f1_classes']}` (`{summary['fraction_zero_f1_classes']:.10f}`)",
            "",
            "## Root cause of the superseded values",
            "",
            "The historical `0.61845` macro-F1 and `0.5812749695` weighted-F1 came from "
            "the stale rounded `macro_f1_class_data.json` LOMO rows. The historical "
            "`0.5802962149` micro-F1 was then estimated as `sum(recall*n)/819`; because "
            "the recalls were rounded, this was not an integer prediction-level TP count. "
            "The historical validation report also described a separate pairwise-rescue "
            "route (`network_pairwise_correlation_rescue_top3`), which is retained as "
            "historical evidence but is not the current formal endpoint.",
            "",
            "Integer accounting check: `sum(TP) = "
            f"{sum(int(row['tp']) for row in metrics['classes'])}` and `sum(support) = "
            f"{sum(int(row['support']) for row in metrics['classes'])}`.",
            "",
        ]
    )
    PROVENANCE_MD.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--source",
        type=Path,
        default=CANONICAL_FORMAL_PATH,
        help="Full prediction-level source containing route families, or the staged canonical source",
    )
    parser.add_argument(
        "--verify-only",
        action="store_true",
        help="Validate the canonical formal source and print metrics without writing",
    )
    args = parser.parse_args()

    if args.verify_only:
        rows = load_formal_predictions(CANONICAL_FORMAL_PATH)
        source = CANONICAL_FORMAL_PATH
    else:
        rows = stage_predictions(args.source)
        source = args.source
    metrics = compute_lomo_network_metrics(rows)
    if args.verify_only:
        verify_derived(metrics)
    else:
        write_derived(metrics)
        write_provenance(source, metrics)
    print(json.dumps(metrics["summary"], indent=2, ensure_ascii=False))
    for row in metrics["classes"]:  # type: ignore[union-attr]
        print(
            f"{row['class']}: support={row['support']} TP={row['tp']} "
            f"FP={row['fp']} FN={row['fn']} P={row['precision']:.10f} "
            f"R={row['recall']:.10f} F1={row['f1']:.10f}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
