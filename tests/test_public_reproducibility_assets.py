from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_saleem_crosswalk_is_public_and_complete() -> None:
    path = ROOT / "reproducibility" / "crosswalks" / "P2_Bo2023_Saleem_crosswalk.csv"
    with path.open(encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle))
    assert len(rows) == 110
    assert len({row["bo2023_region_id"] for row in rows}) == 110
    assert all(row["locked_network"] for row in rows)


def test_independent_enrichment_frozen_counts_and_panel_hash() -> None:
    base = ROOT / "reproducibility" / "independent_enrichment"
    manifest = json.loads((base / "independent_enrichment_manifest.json").read_text(encoding="utf-8"))
    assert manifest["panel_n"] == 200
    assert manifest["background_n"] == 21668
    assert manifest["gprofiler_primary_significant"] == {"GO:BP": 446, "KEGG": 11}
    assert manifest["gprofiler_sensitivity_significant"] == {"GO:BP": 410, "KEGG": 8}
    assert len(manifest["celltype_results"]) == 14

    panel = ROOT / "data" / "models" / "bo2023_saleem_network_top200_model_genes.csv"
    expected = manifest["input_sha256"]["code\\data\\models\\bo2023_saleem_network_top200_model_genes.csv"]
    assert _sha256(panel) == expected


def test_s18_s19_public_inputs_are_present() -> None:
    required = [
        ROOT / "reproducibility" / "v4_p0_10_lambda_sensitivity.csv",
        ROOT / "reproducibility" / "v4_p0_10_lambda_friedman.csv",
        ROOT / "reproducibility" / "v4_p0_5_ahba_trace.csv",
        ROOT / "supplementary" / "p1_sensitivity_analysis.csv",
        ROOT / "reproducibility" / "s18_s19" / "run_s18_s19_sensitivity.py",
    ]
    assert all(path.is_file() and path.stat().st_size > 0 for path in required)
