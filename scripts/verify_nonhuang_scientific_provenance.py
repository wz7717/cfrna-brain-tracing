#!/usr/bin/env python3
"""Arithmetic QA for the non-Huang scientific-provenance patch candidate."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[1]


def rounded_percent(correct: int, denominator: int) -> float:
    return round(correct / denominator * 100, 2)


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def load_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def check_fraction(
    checks: list[dict[str, str]],
    name: str,
    correct: int,
    denominator: int,
    expected_correct: int,
    expected_denominator: int,
    expected_percent: float,
) -> None:
    require((correct, denominator) == (expected_correct, expected_denominator), f"{name}: unexpected fraction {correct}/{denominator}")
    require(rounded_percent(correct, denominator) == expected_percent, f"{name}: unexpected percent")
    checks.append({"check": name, "status": "PASS", "value": f"{correct}/{denominator}={expected_percent:.2f}%"})


def check_ahba(root: Path, checks: list[dict[str, str]]) -> None:
    summary = load_json(root / "reproducibility" / "ahba" / "ahba_endpoint_evaluability_summary.json")
    endpoints = summary["endpoint_evaluability"]
    require(endpoints["replicate_collapsed_tissues"] == 231, "AHBA tissue count changed")
    require(endpoints["network"]["evaluable"] == 223, "AHBA Network denominator changed")
    require(endpoints["resolution_group"]["evaluable"] == 88, "AHBA group denominator changed")
    require(endpoints["exact_region"]["evaluable"] == 88, "AHBA exact denominator changed")
    checks.append({"check": "AHBA endpoint counts", "status": "PASS", "value": "231 tissues; Network n=223; group/exact n=88"})
    for name, endpoint, topk, correct, denominator, percent in (
        ("AHBA Network Top1", "network", "top1", 165, 223, 73.99),
        ("AHBA Network Top3", "network", "top3", 211, 223, 94.62),
        ("AHBA group Top1", "resolution_group", "top1", 37, 88, 42.05),
        ("AHBA group Top3", "resolution_group", "top3", 60, 88, 68.18),
        ("AHBA exact Top1", "exact_region", "top1", 24, 88, 27.27),
        ("AHBA exact Top3", "exact_region", "top3", 40, 88, 45.45),
        ("AHBA Network single-label Top1", "network_unique_single_label_sensitivity", "top1", 39, 56, 69.64),
        ("AHBA Network single-label Top3", "network_unique_single_label_sensitivity", "top3", 44, 56, 78.57),
        ("AHBA group single-label Top1", "group_single_label_sensitivity", "top1", 19, 40, 47.50),
        ("AHBA group single-label Top3", "group_single_label_sensitivity", "top3", 29, 40, 72.50),
        ("AHBA exact single-label Top1", "exact_single_label_sensitivity", "top1", 11, 40, 27.50),
        ("AHBA exact single-label Top3", "exact_single_label_sensitivity", "top3", 20, 40, 50.00),
    ):
        value = endpoints[endpoint][topk]
        check_fraction(checks, name, value["correct"], value["n"], correct, denominator, percent)
    ledger = load_csv(root / "reproducibility" / "ahba" / "ahba_endpoint_evaluability_ledger.csv")
    require(len(ledger) == 231, "AHBA ledger row count changed")
    require(all(not row["group_truth_label_count"] for row in ledger), "AHBA ledger invented group truth-label counts")
    checks.append({"check": "AHBA ledger semantics", "status": "PASS", "value": "group_truth_label_count intentionally blank"})


def row_for(rows: list[dict[str, str]], truth_basis: str, level: str) -> dict[str, str]:
    matches = [row for row in rows if row["truth_basis"] == truth_basis and row["level"] == level]
    require(len(matches) == 1, f"expected exactly one TCGA row for {truth_basis}/{level}")
    return matches[0]


def check_tcga(root: Path, checks: list[dict[str, str]]) -> None:
    summary = load_json(root / "reproducibility" / "tcga_brats_truth_basis_top3_summary.json")
    require(summary["total_paired_cases"] == 65, "TCGA paired total changed")
    require(summary["primary_edema_comparator"]["n"] == 63, "TCGA primary edema denominator changed")
    exclusions = summary["primary_edema_comparator"]["excluded_cases"]
    require(len(exclusions) == 2, "TCGA edema exclusions changed")
    checks.append({"check": "TCGA primary edema comparator", "status": "PASS", "value": "63 of 65 paired cases; 2 source-derived exclusions"})
    rows = load_csv(root / "reproducibility" / "tcga_brats_truth_basis_top3_summary.csv")
    for level, expected in {
        "network": {"center": (19, 65, 29.23), "core": (12, 65, 18.46), "edema": (15, 63, 23.81), "whole_tumor": (10, 65, 15.38)},
        "broad": {"center": (32, 65, 49.23), "core": (45, 65, 69.23), "edema": (52, 63, 82.54), "whole_tumor": (46, 65, 70.77)},
    }.items():
        for truth_basis, (correct, denominator, percent) in expected.items():
            row = row_for(rows, truth_basis, level)
            check_fraction(checks, f"TCGA {level} {truth_basis} strict Top3", int(row["correct"]), int(row["n"]), correct, denominator, percent)
    require(round(summary["range_across_truth_bases_percentage_points"]["network"], 2) == 13.85, "TCGA Network range changed")
    require(round(summary["range_across_truth_bases_percentage_points"]["broad"], 2) == 33.31, "TCGA broad range changed")
    checks.append({"check": "TCGA strict Top3 ranges", "status": "PASS", "value": "Network 13.85 pp; broad 33.31 pp"})


def check_orthology(root: Path, checks: list[dict[str, str]]) -> None:
    summary = load_json(root / "reproducibility" / "orthology_humanization_summary.json")
    rows = summary["region_signature_rows"]
    require(rows["unit"] == "gene-by-region row occurrence", "orthology denominator unit changed")
    check_fraction(checks, "Orthology humanized row occurrences", rows["humanized"]["correct"], rows["humanized"]["n"], 5324, 8800, 60.50)
    check_fraction(checks, "Orthology unmapped row occurrences", rows["unmapped"]["correct"], rows["unmapped"]["n"], 3476, 8800, 39.50)
    require(rows["humanized"]["correct"] + rows["unmapped"]["correct"] == rows["total"], "orthology row partition changed")
    check_fraction(
        checks,
        "Top200 orthology humanizable",
        summary["network_top200_orthology_humanizable"]["correct"],
        summary["network_top200_orthology_humanizable"]["n"],
        188,
        200,
        94.00,
    )
    require(summary["gprofiler_mapped"]["mapped"] == 179 and summary["gprofiler_mapped"]["input_panel_n"] == 200, "g:Profiler mapping count changed")
    checks.append({"check": "g:Profiler mapping universe", "status": "PASS", "value": "179/200, distinct from frozen orthology universe"})


def check_tier(root: Path, checks: list[dict[str, str]]) -> None:
    summary = load_json(root / "reproducibility" / "tier_cascade_loso_summary.json")
    exact = summary["exact_evaluable"]
    network = summary["network_candidate_set_within_same_universe"]
    require(exact["top3_hit"] + exact["top3_miss"] == exact["n"] == 814, "tier exact-evaluable universe changed")
    require(network["truth_retained"] + network["truth_missed"] == exact["n"], "tier Network status universe changed")
    checks.append({"check": "Tier exact-evaluable universe", "status": "PASS", "value": "368+446=814 and 750+64=814"})
    check_fraction(
        checks,
        "Tier Network candidate misses among exact Top3 misses",
        summary["network_candidate_miss_share_of_exact_top3_misses"]["correct"],
        summary["network_candidate_miss_share_of_exact_top3_misses"]["n"],
        64,
        446,
        14.35,
    )
    check_fraction(
        checks,
        "Tier Exact Top3 conditional on Network truth retained",
        summary["exact_top3_given_network_truth_retained"]["correct"],
        summary["exact_top3_given_network_truth_retained"]["n"],
        368,
        750,
        49.07,
    )
    check_fraction(
        checks,
        "Tier group Top3 conditional on Network truth retained",
        summary["group_top3_given_network_truth_retained"]["correct"],
        summary["group_top3_given_network_truth_retained"]["n"],
        590,
        750,
        78.67,
    )
    require(summary["recovery_after_network_candidate_miss"]["exact_top3"]["correct"] == 0, "tier exact recovery changed")
    require(summary["recovery_after_network_candidate_miss"]["group_top3"]["correct"] == 0, "tier group recovery changed")
    checks.append({"check": "Tier recovery after Network candidate miss", "status": "PASS", "value": "0% exact and group"})


def check_sign_flip(root: Path, checks: list[dict[str, str]]) -> None:
    rows = load_csv(root / "reproducibility" / "sign_flip_current_family.csv")
    expected = [
        ("Network Top1", "0.031250", "0.125000"),
        ("Network Top3", "0.375000", "0.500000"),
        ("resolution-group Top3", "0.593750", "0.593750"),
        ("exact-region Top3", "0.324219", "0.500000"),
    ]
    observed = [(row["endpoint"], row["raw_p_display"], row["bh_p_m4_display"]) for row in rows]
    require(observed == expected, "current sign-flip family changed")
    require(all(row["significant_bh_0_05"].lower() == "false" for row in rows), "a sign-flip endpoint is BH significant")
    checks.append({"check": "Current four-test sign-flip family", "status": "PASS", "value": "raw/BH values verified; none significant"})


def run_checks(root: Path) -> dict[str, Any]:
    checks: list[dict[str, str]] = []
    check_ahba(root, checks)
    check_tcga(root, checks)
    check_orthology(root, checks)
    check_tier(root, checks)
    check_sign_flip(root, checks)
    return {"status": "PASS", "checks": checks, "n_checks": len(checks)}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, default=REPO_ROOT)
    parser.add_argument("--write-json", type=Path)
    args = parser.parse_args()
    payload = run_checks(args.repo_root.resolve())
    if args.write_json:
        args.write_json.parent.mkdir(parents=True, exist_ok=True)
        args.write_json.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
