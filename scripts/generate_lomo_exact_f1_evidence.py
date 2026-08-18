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
GENERATOR = "scripts/generate_lomo_exact_f1_evidence.py"
GENERATOR_INPUT_BINDING = "core.lomo_exact_f1.CANONICAL_FORMAL_PATH"

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


def _canonical_path(path: Path) -> str:
    """Record a supplied external origin verbatim as an absolute path."""

    return str(path.resolve())


def source_and_origin(source: Path | None) -> tuple[Path, str, str]:
    """Return the staging input and retain its distinct origin provenance.

    A verify/regenerate run without ``--source`` consumes the already staged
    table.  It must not replace the external-origin hash with the staged hash.
    """

    if source is not None:
        return source, _canonical_path(source), sha256_file(source)

    previous: dict[str, object] = {}
    if PROVENANCE.exists():
        previous = json.loads(PROVENANCE.read_text(encoding="utf-8"))
    chain = previous.get("source_chain", {})
    origin = chain.get("origin", {}) if isinstance(chain, dict) else {}
    legacy = previous.get("formal_source", {})
    if not isinstance(legacy, dict):
        legacy = {}
    origin_path = str(
        origin.get("path")
        or legacy.get("origin_path")
        or legacy.get("staged_from")
        or CANONICAL_FORMAL_PATH.relative_to(ROOT).as_posix()
    )
    origin_sha256 = str(
        origin.get("sha256")
        or legacy.get("sha256")
        or sha256_file(CANONICAL_FORMAL_PATH)
    )
    return CANONICAL_FORMAL_PATH, origin_path, origin_sha256


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


def update_macro_json(
    metrics: dict[str, object], origin_path: str, origin_sha256: str
) -> None:
    payload = json.loads(MACRO_JSON.read_text(encoding="utf-8"))
    data = _replace_endpoint_rows(list(payload["data"]))
    data.extend(macro_class_rows(metrics))
    data.extend(_summary_rows(metrics["summary"]))  # type: ignore[arg-type]
    payload["data"] = data
    provenance = dict(payload.get("provenance", {}))
    provenance.update(
        {
            "formal_lomo_exact_origin_path": origin_path,
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

    # Preserve all endpoint rows and use a canonical summary block at EOF.
    # This makes source staging idempotent instead of moving unrelated endpoint
    # summaries when LOMO Exact is regenerated.
    exact_rows = macro_class_rows(metrics)
    summary = metrics["summary"]  # type: ignore[assignment]
    exact_summaries = {
        "LOMO_Exact_macro": {
            "endpoint": "SUMMARY",
            "class": "LOMO_Exact_macro",
            "n": str(summary["n_classes"]),
            "precision": "",
            "recall": "",
            "f1": f"{float(summary['macro_f1']):.4f}",
        },
        "LOMO_Exact_weighted": {
            "endpoint": "SUMMARY",
            "class": "LOMO_Exact_weighted",
            "n": str(summary["n_samples"]),
            "precision": "",
            "recall": "",
            "f1": f"{float(summary['weighted_f1']):.4f}",
        },
    }
    retained: list[dict[str, object]] = []
    existing_summaries: dict[str, dict[str, str]] = {}
    existing_summary_order: list[str] = []
    exact_rows_written = False
    for row in existing:
        if row.get("endpoint") == "SUMMARY":
            class_name = str(row.get("class", ""))
            existing_summaries[class_name] = row
            if class_name not in existing_summary_order:
                existing_summary_order.append(class_name)
            continue
        if row.get("endpoint") == "LOMO_Exact":
            if not exact_rows_written:
                retained.extend(exact_rows)
                exact_rows_written = True
            continue
        retained.append(row)
    if not exact_rows_written:
        retained.extend(exact_rows)
    existing_summaries.update(exact_summaries)
    canonical_summary_order = [
        "LOMO_Exact_macro",
        "LOMO_Exact_weighted",
        "LOMO_Network_macro",
        "LOMO_Network_weighted",
        "LOSO_Exact_macro",
        "LOSO_Exact_weighted",
        "LOSO_Network_macro",
        "LOSO_Network_weighted",
    ]
    for class_name in canonical_summary_order + existing_summary_order:
        if class_name in existing_summaries:
            retained.append(existing_summaries.pop(class_name))
    retained.extend(existing_summaries.values())
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


def write_derived(
    metrics: dict[str, object], origin_path: str, origin_sha256: str
) -> None:
    write_csv(FORMAL_DETAIL, metrics["classes"], DETAIL_FIELDS)  # type: ignore[arg-type]
    update_macro_json(metrics, origin_path, origin_sha256)
    update_macro_csv(metrics)
    update_cross3(metrics)


def write_provenance(
    origin_path: str, metrics: dict[str, object], origin_sha256: str
) -> None:
    summary = metrics["summary"]  # type: ignore[assignment]
    staged_path = CANONICAL_FORMAL_PATH.relative_to(ROOT).as_posix()
    staged_sha256 = sha256_file(CANONICAL_FORMAL_PATH)
    source_chain = {
        "origin": {
            "path": origin_path,
            "sha256": origin_sha256,
        },
        "staged": {
            "path": staged_path,
            "sha256": staged_sha256,
        },
        "generator_input": {
            "path": staged_path,
            "sha256": staged_sha256,
            "consumer": GENERATOR,
            "binding": GENERATOR_INPUT_BINDING,
            "equals_staged": True,
        },
        "staging_transform": (
            "select the 812 frozen route rows and serialize canonical prediction fields "
            "with recomputed hit1/hit3 indicators"
        ),
    }
    payload = {
        "schema": "braintrace.lomo_exact_f1_provenance.v2",
        "source_chain": source_chain,
        "formal_source": {
            "path": staged_path,
            "sha256": origin_sha256,
            "origin_path": origin_path,
            "staged_sha256": staged_sha256,
            "staged_from": origin_path,
            "generator_input_path": staged_path,
            "generator_input_sha256": staged_sha256,
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
        "## Origin / staged / generator-input pairing",
        "",
        "| Role | Path | SHA-256 | Binding |",
        "| --- | --- | --- | --- |",
        f"| Origin | `{origin_path}` | `{origin_sha256}` | frozen external formal prediction detail |",
        f"| Staged | `{staged_path}` | `{staged_sha256}` | repository-staged canonical table |",
        f"| Generator input | `{staged_path}` | `{staged_sha256}` | `{GENERATOR}` via `{GENERATOR_INPUT_BINDING}` |",
        "",
        "The staged and generator-input path/SHA pairs must be identical; only the origin may differ because staging normalizes the frozen route rows.",
        "",
        f"- Route: `{FORMAL_ROUTE}`",
        f"- Route family: `{FORMAL_ROUTE_FAMILY}`",
        f"- Origin SHA-256: `{origin_sha256}`",
        f"- Staged / generator-input SHA-256: `{staged_sha256}`",
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
    chain = provenance.get("source_chain", {})
    if not isinstance(chain, dict):
        raise ValueError("LOMO Exact provenance is missing its source chain")
    staged = chain.get("staged", {})
    generator_input = chain.get("generator_input", {})
    if not isinstance(staged, dict) or not isinstance(generator_input, dict):
        raise ValueError("LOMO Exact source chain is malformed")
    if staged.get("path") != generator_input.get("path"):
        raise ValueError("LOMO Exact staged and generator-input paths differ")
    if staged.get("sha256") != generator_input.get("sha256"):
        raise ValueError("LOMO Exact staged and generator-input SHA-256 values differ")
    if source["staged_sha256"] != sha256_file(CANONICAL_FORMAL_PATH):
        raise ValueError("LOMO Exact source changed without regeneration")
    if staged.get("sha256") != sha256_file(CANONICAL_FORMAL_PATH):
        raise ValueError("LOMO Exact source-chain staged SHA-256 is stale")
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
        source, origin_path, origin_sha256 = source_and_origin(args.source)
        rows = stage_predictions(source)
        metrics = compute_lomo_exact_metrics(rows)
        write_derived(metrics, origin_path, origin_sha256)
        write_provenance(origin_path, metrics, origin_sha256)
    print(json.dumps(metrics["summary"], indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
