from __future__ import annotations

import csv
import json
import math
from pathlib import Path

from core.lomo_exact_f1 import (
    CANONICAL_FORMAL_PATH,
    FORMAL_N,
    FORMAL_N_CLASSES,
    FORMAL_ROUTE,
    FORMAL_ROUTE_FAMILY,
    FORMAL_TOP1,
    compute_lomo_exact_metrics,
    load_formal_predictions,
)
from core.lomo_f1 import sha256_file


ROOT = Path(__file__).resolve().parents[1]
EXACT_ORIGIN_SHA256 = "EB5F10F01B122F68D09256EA6866DEAE2B439AABAD27E076181EC8760E7AAF36"


def test_lomo_exact_prediction_level_integer_accounting() -> None:
    rows = load_formal_predictions(CANONICAL_FORMAL_PATH)
    metrics = compute_lomo_exact_metrics(rows)
    summary = metrics["summary"]

    assert len(rows) == FORMAL_N
    assert {row["route"] for row in rows} == {FORMAL_ROUTE}
    assert {row["route_family"] for row in rows} == {FORMAL_ROUTE_FAMILY}
    assert summary["n_classes"] == FORMAL_N_CLASSES
    assert summary["top1_correct"] == FORMAL_TOP1
    assert sum(int(row["tp"]) for row in metrics["classes"]) == FORMAL_TOP1
    assert sum(int(row["support"]) for row in metrics["classes"]) == FORMAL_N
    assert math.isclose(summary["micro_f1"], FORMAL_TOP1 / FORMAL_N, abs_tol=1e-15)
    assert summary["predicted_top1_labels_outside_truth_universe"] == ["36c", "Cla", "MT"]


def test_lomo_exact_derived_artifacts_share_current_source() -> None:
    metrics = compute_lomo_exact_metrics(load_formal_predictions(CANONICAL_FORMAL_PATH))
    summary = metrics["summary"]
    provenance = json.loads(
        (ROOT / "reproducibility" / "lomo_exact_region_f1_provenance.json").read_text(
            encoding="utf-8"
        )
    )
    assert provenance["formal_source"]["staged_sha256"] == sha256_file(
        CANONICAL_FORMAL_PATH
    )
    assert provenance["prediction_level_metrics"]["micro_f1"] == summary["micro_f1"]

    macro = json.loads(
        (ROOT / "reproducibility" / "macro_f1_class_data.json").read_text(
            encoding="utf-8"
        )
    )
    assert macro["provenance"]["formal_lomo_exact_source_sha256"] == provenance[
        "formal_source"
    ]["sha256"]
    qa = json.loads((ROOT / "SCIENTIFIC_REMEDIATION_QA.json").read_text(encoding="utf-8"))
    assert qa["lomo_exact_f1"]["source_sha256"] == provenance["formal_source"]["sha256"]
    lomo_rows = [row for row in macro["data"] if row.get("endpoint") == "LOMO_Exact"]
    assert len(lomo_rows) == FORMAL_N_CLASSES
    assert sum(int(row["n"]) for row in lomo_rows) == FORMAL_N

    with (
        ROOT / "reproducibility" / "p1_cross1_5" / "cross3_f1_distribution_summary.csv"
    ).open(newline="", encoding="utf-8-sig") as handle:
        cross3 = next(row for row in csv.DictReader(handle) if row["endpoint"] == "LOMO_Exact")
    assert math.isclose(float(cross3["macro_f1"]), summary["macro_f1"], abs_tol=1e-15)
    assert math.isclose(float(cross3["micro_f1"]), FORMAL_TOP1 / FORMAL_N, abs_tol=1e-15)

    with (ROOT / "reproducibility" / "formal_lomo_exact_region_f1.csv").open(
        newline="", encoding="utf-8-sig"
    ) as handle:
        class_rows = list(csv.DictReader(handle))
    assert len(class_rows) == FORMAL_N_CLASSES
    assert sum(int(row["tp"]) for row in class_rows) == FORMAL_TOP1
    assert sum(int(row["support"]) for row in class_rows) == FORMAL_N


def test_lomo_exact_origin_staged_and_generator_input_are_explicitly_paired() -> None:
    provenance = json.loads(
        (ROOT / "reproducibility" / "lomo_exact_region_f1_provenance.json").read_text(
            encoding="utf-8"
        )
    )
    chain = provenance["source_chain"]
    origin = chain["origin"]
    staged = chain["staged"]
    generator_input = chain["generator_input"]

    assert origin["sha256"] == EXACT_ORIGIN_SHA256
    assert origin["path"].replace("\\", "/").endswith(
        "/lomo/formal_lomo_exact_region_detail.csv"
    )
    assert staged["path"] == CANONICAL_FORMAL_PATH.relative_to(ROOT).as_posix()
    assert staged["sha256"] == sha256_file(CANONICAL_FORMAL_PATH)
    assert generator_input["path"] == staged["path"]
    assert generator_input["sha256"] == staged["sha256"]
    assert generator_input["equals_staged"] is True
    assert generator_input["consumer"] == "scripts/generate_lomo_exact_f1_evidence.py"
    assert generator_input["binding"] == "core.lomo_exact_f1.CANONICAL_FORMAL_PATH"
