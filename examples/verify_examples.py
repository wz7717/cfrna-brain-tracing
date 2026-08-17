#!/usr/bin/env python3
"""Verify BrainTrace's public synthetic smoke-test examples."""

from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Any

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import core.inference as locked_inference  # noqa: E402
from core.cli import _json_default, _query_payload  # noqa: E402
from core.model_lock import EXPECTED_MODEL_LOCK_ID  # noqa: E402
from core.query_input import read_expression_file  # noqa: E402


EXAMPLES = (
    ("braintrace_example_counts.tsv", "expected_output_counts.json", "raw_counts"),
    ("braintrace_example_logcpm.tsv", "expected_output_logcpm.json", "logcpm"),
)
NETWORK_GENES = ROOT / "data/models/bo2023_saleem_network_top200_model_genes.csv"
REGION_REFERENCE = ROOT / "data/models/bo2023_formal_region_logcpm_reference_matrix.npz"
ALLOWED_ARTIFACT_ROOT = (ROOT / "data/models").resolve()


def _json_ready(value: Any) -> Any:
    return json.loads(json.dumps(value, ensure_ascii=False, default=_json_default))


def _assert_public_packaged_reference(region_out: dict[str, Any]) -> None:
    source = str(region_out.get("meta", {}).get("reference_expression_source", ""))
    if "packaged" not in source.lower() and "formal" not in source.lower():
        raise AssertionError(f"Unexpected non-packaged region reference source: {source!r}")
    if not REGION_REFERENCE.resolve().is_relative_to(ALLOWED_ARTIFACT_ROOT):
        raise AssertionError("Regional reference escaped the public packaged-artifact directory.")


def _input_overlap(expression: pd.DataFrame) -> tuple[int, int]:
    query_genes = set(expression["gene_symbol"].astype(str))
    network_genes = set(pd.read_csv(NETWORK_GENES)["gene_symbol"].astype(str))
    network_overlap = len(query_genes & network_genes)

    import numpy as np

    with np.load(REGION_REFERENCE, allow_pickle=False) as archive:
        regional_genes = set(archive["genes"].astype(str))
    return network_overlap, len(query_genes & regional_genes)


def run_example(input_path: Path) -> tuple[pd.DataFrame, str, dict[str, Any], dict[str, Any]]:
    """Run one example through the public locked production entry point."""
    expression, query_source = read_expression_file(input_path)
    network_out, region_out = locked_inference.run_locked_three_tier_route(expression)
    return expression, query_source, network_out, region_out


def run_and_verify(input_name: str, expected_name: str, expected_source: str) -> None:
    input_path = Path(__file__).with_name(input_name)
    expected_path = Path(__file__).with_name(expected_name)
    expression, query_source, network_out, region_out = run_example(input_path)
    if query_source != expected_source:
        raise AssertionError(f"{input_name}: expected source {expected_source}, got {query_source}")

    network_overlap, regional_overlap = _input_overlap(expression)
    if network_overlap < 100 or network_overlap / 200 < 0.50:
        raise AssertionError(f"{input_name}: Network overlap gate failed ({network_overlap}/200)")
    if regional_overlap < 20:
        raise AssertionError(f"{input_name}: regional overlap gate failed ({regional_overlap})")

    if not network_out.get("results"):
        raise AssertionError(f"{input_name}: empty Network output")
    group_rows = region_out.get("meta", {}).get("region_resolution_annotation", {}).get("group_ranking", [])
    if not group_rows:
        raise AssertionError(f"{input_name}: empty Resolution Group output")
    if not region_out.get("results"):
        raise AssertionError(f"{input_name}: empty Exact-Region output")

    lock = network_out.get("meta", {}).get("model_lock", {})
    if lock.get("lock_id") != EXPECTED_MODEL_LOCK_ID:
        raise AssertionError(f"{input_name}: wrong model lock {lock!r}")
    _assert_public_packaged_reference(region_out)

    actual = _json_ready(_query_payload(network_out, region_out, query_source))
    expected = json.loads(expected_path.read_text(encoding="utf-8"))
    if actual != expected:
        raise AssertionError(f"{input_name}: output differs from {expected_name}")

    print(
        f"PASS {input_name}: Network {network_overlap}/200; "
        f"regional {regional_overlap}; lock {EXPECTED_MODEL_LOCK_ID}"
    )


def main() -> int:
    for item in EXAMPLES:
        run_and_verify(*item)
    print("PASS: both public examples match committed locked-route outputs")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
