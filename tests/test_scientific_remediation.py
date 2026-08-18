from __future__ import annotations

import json
import math
from pathlib import Path

from core.lomo_f1 import sha256_file
from core.resolution_group_baselines import (
    CANONICAL_PATHS,
    compute_baseline_record,
    load_formal_rows,
)


ROOT = Path(__file__).resolve().parents[1]
GROUP_ORIGIN_SHA256 = "B9A17D20BA434F52BD812FAE361E1A3F51C55B705AEFA745B8853544276390F1"


def _qa() -> dict:
    return json.loads((ROOT / "SCIENTIFIC_REMEDIATION_QA.json").read_text(encoding="utf-8"))


def test_tcga_broad_range_is_dynamic_max_minus_min() -> None:
    qa = _qa()
    range_record = qa["tcga_broad_strict_top3_range"]
    rates = list(range_record["source_rates_percentage_points"].values())
    assert math.isclose(
        range_record["reported_range_percentage_points"],
        max(rates) - min(rates),
        abs_tol=1e-12,
    )
    assert range_record["formula"] == "max(source_rates) - min(source_rates)"


def test_benchmark_has_profiles_events_and_separate_memory_stages() -> None:
    benchmark = _qa()["benchmark"]
    assert benchmark["warm_events"] == benchmark["n_profiles"] * benchmark["n_warm_repeats"]
    assert "warm_samples" not in benchmark
    assert benchmark["cold"]["peak_working_set_mib"] != benchmark["warm"][
        "maximum_working_set_mib"
    ]
    assert benchmark["source"]["staged_sha256"] == sha256_file(
        ROOT / benchmark["source"]["staged_path"]
    )


def test_resolution_group_baselines_are_current_source_derived() -> None:
    qa = _qa()
    payload = qa["resolution_group_random_baselines"]
    actual = {record["endpoint"]: record for record in payload["records"]}
    for endpoint in ("LOSO", "LOMO"):
        source = payload["sources"][endpoint]
        assert source["staged_sha256"] == sha256_file(ROOT / source["staged_path"])
        rows, _ = load_formal_rows(CANONICAL_PATHS[endpoint], endpoint)
        expected = compute_baseline_record(rows, endpoint)
        for field in (
            "n_profiles",
            "observed_hits",
            "uniform_random_rate",
            "weighted_random_rate",
        ):
            if isinstance(expected[field], float):
                assert math.isclose(actual[endpoint][field], expected[field], abs_tol=1e-15)
            else:
                assert actual[endpoint][field] == expected[field]


def test_lomo_exact_and_group_path_sha_pairs_are_explicit_and_current() -> None:
    qa = _qa()
    pairing = qa["lomo_input_path_sha_pairing"]
    on_disk = json.loads(
        (ROOT / "reproducibility" / "lomo_input_chain_provenance.json").read_text(
            encoding="utf-8"
        )
    )
    assert on_disk == pairing

    exact = pairing["chains"]["LOMO Exact"]
    group = pairing["chains"]["LOMO Group"]
    for label, chain in (("LOMO Exact", exact), ("LOMO Group", group)):
        origin = chain["origin"]
        staged = chain["staged"]
        generator_input = chain["generator_input"]
        assert origin["path"]
        assert len(origin["sha256"]) == 64
        assert staged["path"] == generator_input["path"], label
        assert staged["sha256"] == generator_input["sha256"], label
        assert generator_input["equals_staged"] is True
        assert staged["sha256"] == sha256_file(ROOT / staged["path"]), label

    assert group["origin"]["sha256"] == GROUP_ORIGIN_SHA256
    assert group["origin"]["path"].replace("\\", "/").endswith(
        "/lomo/formal_lomo_resolution_group_detail.csv"
    )
    assert group["generator_input"]["consumer"] == (
        "scripts/generate_scientific_remediation_artifacts.py"
    )
    assert group["generator_input"]["binding"] == (
        "core.resolution_group_baselines.CANONICAL_PATHS['LOMO']"
    )

    baseline_chain = qa["resolution_group_random_baselines"]["sources"]["LOMO"][
        "source_chain"
    ]
    assert baseline_chain == group

    pairing_markdown = (
        ROOT / "reproducibility" / "LOMO_INPUT_CHAIN_PROVENANCE.md"
    ).read_text(encoding="utf-8")
    assert "LOMO Exact" in pairing_markdown
    assert "LOMO Group" in pairing_markdown
    assert GROUP_ORIGIN_SHA256 in pairing_markdown


def test_only_supported_friedman_claim_is_retained() -> None:
    qa = _qa()
    friedman = qa["friedman_exact_enumeration"]
    assert friedman["exact_enumeration"] == "REMOVED"
    assert friedman["chi2"] == 0.5385
    assert friedman["df"] == 2
    assert friedman["p_value"] == 0.764
    assert qa["stale_current_value_scan"]["status"] == "PASS"
