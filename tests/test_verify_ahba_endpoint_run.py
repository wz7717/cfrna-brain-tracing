from __future__ import annotations

import csv
from pathlib import Path

from scripts.verify_ahba_endpoint_run import EXPECTED, ROUTE, validate


def test_ahba_endpoint_verifier_rejects_changed_fraction(tmp_path: Path) -> None:
    path = tmp_path / "detail.csv"
    fields = ["route", *EXPECTED]
    rows = []
    for index in range(231):
        row = {"route": ROUTE}
        for column, (correct, total) in EXPECTED.items():
            row[column] = "" if index >= total else str(index < correct)
        rows.append(row)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    assert validate(path)["status"] == "PASS"
    rows[0]["network_top1_hit"] = "False"
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    assert validate(path)["status"] == "BLOCKED: SCIENTIFIC_OUTPUT_DRIFT"
