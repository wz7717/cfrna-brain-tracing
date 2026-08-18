#!/usr/bin/env python
"""Stage and regenerate formal prediction-level LOMO Exact F1 evidence."""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from core.lomo_exact_f1 import (  # noqa: E402
    CANONICAL_FORMAL_PATH,
    FORMAL_ROUTE,
    FORMAL_ROUTE_FAMILY,
    compute_lomo_exact_metrics,
    cross3_summary_row,
    load_formal_predictions,
    macro_class_rows,
)
from core.lomo_f1 import sha256_file  # noqa: E402


MACRO_JSON = ROOT / "reproducibility" / "macro_f1_class_data.json"
MACRO_CSV = ROOT / "reproducibility" / "v4_p0_13_macro_f1.csv"
CROSS3_CSV = (
    ROOT / "reproducibility" / "p1_cross1_5" / "cross3_f1_distribution_summary.csv"
)
FORMAL_DETAIL = ROOT / "reproducibility" / "formal_lomo_exact_region_f1.csv"
PROVENANCE = ROOT / "reproducibility" / "lomo_exact_region_f1_provenance.json"
PROVENANCE_MD = ROOT / "reproducibility" / "LOMO_EXACT_REGION_F1_PROVENANCE.md"

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


def _replace_endpoint_rows(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    return [
        row
        for row in rows
        if row.get("endpoint") != "LOMO_Exact"
        and not (
            row.get("endpoint") == "SUMMARY"
            and str(row.get("class", "")).startswith("LOMO_Exact_")
        )
    ]


def _summary_rows(summary: dict[str, object]) -> list[dict[str, str]]:
    return [
        {
            "endpoint": "SUMMARY",
            "class": "LOMO_Exact_macro",
            "n": str(summary["n_classes"]),
            "precision": "",
            "recall": "",
            "f1": f"{float(summary['macro_f1']):.10f}",
        },
        {
            "endpoint": "SUMMARY",
            "class": "LOMO_Exact_weighted",
            "n": str(summary["n_samples"]),
            "precision": "",
            "recall": "",
            "f1": f"{float(summary['weighted_f1']):.10f}",
        },
    ]


def update_macro_json(metrics: dict[str, object], origin_sha256: str) -> None:
    payload = json.loads(MACRO_JSON.read_text(encoding="utf-8"))
    data = _replace_endpoint_rows(list(payload["data"]))
    data.extend(macro_class_rows(metrics))
    data.extend(_summary_rows(metrics["summary"]))  # type: ignore[arg-type]
    payload["data"] = data
    provenance = dict(payload.get("provenance", {}))
    provenance.update(
        {
            "formal_lomo_exact_source": CANONICAL_FORMAL_PATH.relative_to(ROOT).as_posix(),
            "formal_lomo_exact_source_sha256": origin_sha256,
            "formal_lomo_exact_staged_sha256": sha256_file(CANONICAL_FORMAL_PATH),
            "formal_lomo_exact_route": FORMAL_ROUTE,
            "formal_lomo_exact_route_family": FORMAL_ROUTE_FAMILY,
            "formal_lomo_exact_truth_class_universe": True,
            "formal_lomo_exact_micro_f1_definition": "sum integer TP / sum support = Top1 accuracy",
        }
    )
    payload["provenance"] = provenance
    MACRO_JSON.write_text(
        json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )


def update_macro_csv(metrics: dict[str, object]) -> None:
    with MACRO_CSV.open(newline="", encoding="utf-8-sig") as handle:
        existing = list(csv.DictReader(handle))
    retained = _replace_endpoint_rows(existing)
    retained.extend(macro_class_rows(metrics))
    retained.extend(_summary_rows(metrics["summary"]))  # type: ignore[arg-type]
    write_csv(
        MACRO_CSV,
        retained,
        ["endpoint", "class", "n", "precision", "recall", "f1"],
    )


def update_cross3(metrics: dict[str, object]) -> None:
    with CROSS3_CSV.open(newline="", encoding="utf-8-sig") as handle:
        existing = list(csv.DictReader(handle))
    retained = [row for row in existing if row.get("endpoint") != "LOMO_Exact"]
    write_csv(CROSS3_CSV, retained + [cross3_summary_row(metrics)], CROSS3_FIELDS)


def write_derived(metrics: dict[str, object], origin_sha256: str) -> None:
    write_csv(FORMAL_DETAIL, metrics["classes"], DETAIL_FIELDS)  # type: ignore[arg-type]
    update_macro_json(metrics, origin_sha256)
    update_macro_csv(metrics)
    update_cross3(metrics)


def write_provenance(source: Path, metrics: dict[str, object], origin_sha256: str) -> None:
    summary = metrics["summary"]  # type: ignore[assignment]
    try:
        staged_from = source.resolve().relative_to(ROOT).as_posix()
    except ValueError:
        staged_from = f"external_source::{source.name}"
    payload = {
        "schema": "braintrace.lomo_exact_f1_provenance.v1",
        "formal_source": {
            "path": CANONICAL_FORMAL_PATH.relative_to(ROOT).as_posix(),
            "sha256": origin_sha256,
            "staged_sha256": sha256_file(CANONICAL_FORMAL_PATH),
            "staged_from": staged_from,
            "route": FORMAL_ROUTE,
            "route_family": FORMAL_ROUTE_FAMILY,
            "n": summary["n_samples"],
        },
        "prediction_level_metrics": {
            key: summary[key]
            for key in (
                "top1_correct",
                "top1_accuracy",
                "top3_correct",
                "top3_accuracy",
                "macro_f1",
                "sd_class_f1",
                "median_class_f1",
                "q1_class_f1",
                "q3_class_f1",
                "iqr_class_f1",
                "weighted_f1",
                "micro_f1",
                "n_zero_f1_classes",
                "fraction_zero_f1_classes",
                "conditional_macro_f1_nonzero",
                "conditional_median_f1_nonzero",
                "predicted_top1_labels_outside_truth_universe",
            )
        },
        "integer_confusion_accounting": {
            "sum_tp": sum(int(row["tp"]) for row in metrics["classes"]),  # type: ignore[index]
            "sum_support": sum(int(row["support"]) for row in metrics["classes"]),  # type: ignore[index]
            "micro_f1_equals_top1_accuracy": True,
            "class_universe": "104 truth labels; prediction-only labels are FP and do not expand the macro denominator",
        },
        "root_cause": "Superseded LOMO Exact class rows were derived from a stale route and rounded summaries rather than the current frozen prediction-level formal route.",
        "superseded_artifact_status": "historical/superseded; current generators must read the staged formal prediction detail",
    }
    PROVENANCE.write_text(
        json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )

    lines = [
        "# Formal LOMO Exact F1 provenance",
        "",
        "The current LOMO Exact endpoint is regenerated from the frozen prediction-level formal route.",
        "Its macro denominator is the 104-label truth universe; Top1 labels outside that universe",
        "remain false positives and do not create extra macro classes.",
        "",
        f"- Route: `{FORMAL_ROUTE}`",
        f"- Route family: `{FORMAL_ROUTE_FAMILY}`",
        f"- Source SHA-256: `{origin_sha256}`",
        f"- Staged prediction SHA-256: `{sha256_file(CANONICAL_FORMAL_PATH)}`",
        f"- Top1 / Top3: `{summary['top1_correct']}/{summary['n_samples']}`; `{summary['top3_correct']}/{summary['n_samples']}`",
        f"- Prediction-only Top1 labels: `{', '.join(summary['predicted_top1_labels_outside_truth_universe'])}`",
        "",
        "## Derived summary",
        "",
    ]
    for key in (
        "macro_f1",
        "sd_class_f1",
        "median_class_f1",
        "iqr_class_f1",
        "weighted_f1",
        "micro_f1",
        "n_zero_f1_classes",
        "conditional_macro_f1_nonzero",
    ):
        lines.append(f"- {key}: `{summary[key]}`")
    lines.extend(
        [
            "",
            "## Integer accounting",
            "",
            f"`sum(TP) = {sum(int(row['tp']) for row in metrics['classes'])}`; "
            f"`sum(support) = {sum(int(row['support']) for row in metrics['classes'])}`.",
            "",
        ]
    )
    PROVENANCE_MD.write_text("\n".join(lines), encoding="utf-8")


def verify_derived(metrics: dict[str, object]) -> None:
    provenance = json.loads(PROVENANCE.read_text(encoding="utf-8"))
    source = provenance["formal_source"]
    if source["staged_sha256"] != sha256_file(CANONICAL_FORMAL_PATH):
        raise ValueError("LOMO Exact source changed without regeneration")
    summary = metrics["summary"]  # type: ignore[assignment]
    if provenance["prediction_level_metrics"]["micro_f1"] != summary["micro_f1"]:
        raise ValueError("LOMO Exact provenance micro-F1 is stale")
    with CROSS3_CSV.open(newline="", encoding="utf-8-sig") as handle:
        row = next(row for row in csv.DictReader(handle) if row["endpoint"] == "LOMO_Exact")
    expected = cross3_summary_row(metrics)
    for key, value in expected.items():
        if key != "endpoint" and abs(float(row[key]) - float(value)) > 1e-12:
            raise ValueError(f"LOMO Exact CROSS3 mismatch in {key}")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, help="Full current prediction-level source")
    parser.add_argument("--verify-only", action="store_true")
    args = parser.parse_args()

    if args.verify_only:
        rows = load_formal_predictions(CANONICAL_FORMAL_PATH)
        metrics = compute_lomo_exact_metrics(rows)
        verify_derived(metrics)
    else:
        source = args.source or CANONICAL_FORMAL_PATH
        origin_sha256 = sha256_file(source)
        rows = stage_predictions(source)
        metrics = compute_lomo_exact_metrics(rows)
        write_derived(metrics, origin_sha256)
        write_provenance(source, metrics, origin_sha256)
    print(json.dumps(metrics["summary"], indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
