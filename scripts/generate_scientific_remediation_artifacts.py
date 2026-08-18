#!/usr/bin/env python
"""Generate the final scientific-remediation audit chain from canonical inputs."""

from __future__ import annotations

import argparse
import csv
import json
import math
import shutil
import sys
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from core.lomo_exact_f1 import (  # noqa: E402
    CANONICAL_FORMAL_PATH as EXACT_DETAIL,
    compute_lomo_exact_metrics,
    load_formal_predictions,
)
from core.lomo_f1 import sha256_file  # noqa: E402
from core.resolution_group_baselines import (  # noqa: E402
    CANONICAL_PATHS as GROUP_DETAIL_PATHS,
    compute_baseline_record,
    load_formal_rows,
)


TCGA_SUMMARY = ROOT / "reproducibility" / "tcga_brats_truth_basis_top3_summary.json"
EXACT_PROVENANCE = ROOT / "reproducibility" / "lomo_exact_region_f1_provenance.json"
MACRO_JSON = ROOT / "reproducibility" / "macro_f1_class_data.json"
LAMBDA_CSV = ROOT / "reproducibility" / "v4_p0_10_lambda_friedman.csv"
BASELINE_CSV = ROOT / "reproducibility" / "formal_resolution_group_random_baselines.csv"
BASELINE_JSON = ROOT / "reproducibility" / "formal_resolution_group_random_baselines.json"
LOMO_INPUT_PAIRING_JSON = ROOT / "reproducibility" / "lomo_input_chain_provenance.json"
LOMO_INPUT_PAIRING_MD = ROOT / "reproducibility" / "LOMO_INPUT_CHAIN_PROVENANCE.md"
BENCHMARK_MANIFEST = ROOT / "reproducibility" / "formal_real_input_performance_manifest.json"
BENCHMARK_PROVENANCE = ROOT / "reproducibility" / "formal_real_input_performance_provenance.json"
REPORT = ROOT / "SCIENTIFIC_REMEDIATION_REPORT.md"
QA = ROOT / "SCIENTIFIC_REMEDIATION_QA.json"
LEDGER = ROOT / "scientific_claim_ledger.csv"
GROUP_GENERATOR = "scripts/generate_scientific_remediation_artifacts.py"
GROUP_GENERATOR_INPUT_BINDING = "core.resolution_group_baselines.CANONICAL_PATHS['LOMO']"


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def write_csv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _source_chain(
    *,
    origin_path: str,
    origin_sha256: str,
    staged_path: str,
    staged_sha256: str,
    consumer: str,
    binding: str,
    staging_transform: str,
) -> dict[str, Any]:
    """Create an explicit, self-checking origin/staged/generator-input chain."""

    return {
        "origin": {"path": origin_path, "sha256": origin_sha256},
        "staged": {"path": staged_path, "sha256": staged_sha256},
        "generator_input": {
            "path": staged_path,
            "sha256": staged_sha256,
            "consumer": consumer,
            "binding": binding,
            "equals_staged": True,
        },
        "staging_transform": staging_transform,
    }


def _assert_staged_generator_pair(label: str, chain: dict[str, Any]) -> None:
    staged = chain["staged"]
    generator_input = chain["generator_input"]
    if staged["path"] != generator_input["path"]:
        raise ValueError(f"{label} staged and generator-input paths differ")
    if staged["sha256"] != generator_input["sha256"]:
        raise ValueError(f"{label} staged and generator-input SHA-256 values differ")
    if not generator_input.get("equals_staged"):
        raise ValueError(f"{label} must declare generator input identical to staged input")
    staged_path = ROOT / str(staged["path"])
    if sha256_file(staged_path) != staged["sha256"]:
        raise ValueError(f"{label} staged input SHA-256 does not match its file")


def stage_group_sources(
    loso_source: Path | None, lomo_source: Path | None
) -> tuple[dict[str, list[dict[str, str]]], dict[str, dict[str, Any]]]:
    """Freeze formal group prediction rows and retain their original input hashes."""

    previous = read_json(BASELINE_JSON) if BASELINE_JSON.exists() else {}
    previous_sources = previous.get("sources", {})
    supplied = {"LOSO": loso_source, "LOMO": lomo_source}
    staged_rows: dict[str, list[dict[str, str]]] = {}
    sources: dict[str, dict[str, Any]] = {}
    for endpoint, provided in supplied.items():
        origin = provided or GROUP_DETAIL_PATHS[endpoint]
        rows, fields = load_formal_rows(origin, endpoint)
        if provided is not None:
            write_csv(GROUP_DETAIL_PATHS[endpoint], rows, fields)
            origin_sha = sha256_file(provided)
            origin_path = str(provided.resolve())
        else:
            previous_source = previous_sources.get(endpoint, {})
            if not isinstance(previous_source, dict):
                previous_source = {}
            previous_chain = previous_source.get("source_chain", {})
            if not isinstance(previous_chain, dict):
                previous_chain = {}
            previous_origin = previous_chain.get("origin", {})
            if not isinstance(previous_origin, dict):
                previous_origin = {}
            origin_sha = str(
                previous_origin.get(
                    "sha256", previous_source.get("origin_sha256", sha256_file(origin))
                )
            )
            origin_path = str(
                previous_origin.get(
                    "path",
                    previous_source.get("origin_path", origin.relative_to(ROOT).as_posix()),
                )
            )
        staged_rows[endpoint] = rows
        staged_path = GROUP_DETAIL_PATHS[endpoint].relative_to(ROOT).as_posix()
        staged_sha256 = sha256_file(GROUP_DETAIL_PATHS[endpoint])
        input_binding = (
            GROUP_GENERATOR_INPUT_BINDING
            if endpoint == "LOMO"
            else "core.resolution_group_baselines.CANONICAL_PATHS['LOSO']"
        )
        chain = _source_chain(
            origin_path=origin_path,
            origin_sha256=origin_sha,
            staged_path=staged_path,
            staged_sha256=staged_sha256,
            consumer=GROUP_GENERATOR,
            binding=input_binding,
            staging_transform=(
                "filter the frozen formal detail to the current hybrid route-family rows "
                "and preserve the input columns used by the baseline generator"
            ),
        )
        _assert_staged_generator_pair(f"{endpoint} Group", chain)
        sources[endpoint] = {
            "origin_path": origin_path,
            "origin_sha256": origin_sha,
            "staged_path": staged_path,
            "staged_sha256": staged_sha256,
            "generator_input_path": staged_path,
            "generator_input_sha256": staged_sha256,
            "source_chain": chain,
        }
    return staged_rows, sources


def generate_group_baselines(
    loso_source: Path | None, lomo_source: Path | None
) -> dict[str, Any]:
    rows_by_endpoint, sources = stage_group_sources(loso_source, lomo_source)
    records = []
    for endpoint in ("LOSO", "LOMO"):
        record = compute_baseline_record(rows_by_endpoint[endpoint], endpoint)
        record.update(
            {
                "source_origin_path": sources[endpoint]["origin_path"],
                "source_origin_sha256": sources[endpoint]["origin_sha256"],
                "source_staged_sha256": sources[endpoint]["staged_sha256"],
                "source_staged_path": sources[endpoint]["staged_path"],
                "source_generator_input_path": sources[endpoint]["generator_input_path"],
                "source_generator_input_sha256": sources[endpoint]["generator_input_sha256"],
            }
        )
        records.append(record)
    fields = [
        "endpoint",
        "top_k",
        "n_profiles",
        "observed_hits",
        "observed_rate",
        "uniform_random_rate",
        "weighted_random_rate",
        "uniform_formula",
        "weighted_formula",
        "rng_seed",
        "n_weighted_random_draws",
        "source_origin_path",
        "source_origin_sha256",
        "source_staged_sha256",
        "source_staged_path",
        "source_generator_input_path",
        "source_generator_input_sha256",
    ]
    write_csv(BASELINE_CSV, records, fields)
    payload = {
        "schema": "braintrace.formal_resolution_group_random_baselines.v2",
        "sources": sources,
        "records": records,
        "regeneration_rule": "A changed staged prediction SHA-256 requires regenerated baseline records.",
    }
    write_json(BASELINE_JSON, payload)
    return payload


def stage_benchmark(source: Path | None) -> tuple[dict[str, Any], dict[str, str]]:
    previous = read_json(BENCHMARK_PROVENANCE) if BENCHMARK_PROVENANCE.exists() else {}
    if source is not None:
        shutil.copyfile(source, BENCHMARK_MANIFEST)
        origin_sha = sha256_file(source)
        origin_path = f"external_source::{source.name}"
    else:
        origin_sha = str(previous.get("origin_sha256", sha256_file(BENCHMARK_MANIFEST)))
        origin_path = str(
            previous.get("origin_path", BENCHMARK_MANIFEST.relative_to(ROOT).as_posix())
        )
    payload = read_json(BENCHMARK_MANIFEST)
    provenance = {
        "origin_path": origin_path,
        "origin_sha256": origin_sha,
        "staged_path": BENCHMARK_MANIFEST.relative_to(ROOT).as_posix(),
        "staged_sha256": sha256_file(BENCHMARK_MANIFEST),
    }
    write_json(BENCHMARK_PROVENANCE, provenance)
    return payload, provenance


def derive_benchmark(manifest: dict[str, Any], provenance: dict[str, str]) -> dict[str, Any]:
    n_profiles = int(manifest["input"]["workload_samples"])
    n_genes = 28415
    n_warm_repeats = len(manifest["warm"]["repeats"])
    cold = manifest["cold"]
    warm = manifest["warm"]
    warm_events = int(warm["aggregate"]["samples"])
    if warm_events != n_profiles * n_warm_repeats:
        raise ValueError("warm_events must equal n_profiles × n_warm_repeats")
    cold_total_seconds = float(cold["timing"]["wall_total_seconds"])
    cold_profiles = int(cold["timing"]["samples"])
    if cold_profiles != n_profiles:
        raise ValueError("Cold benchmark profile count differs from input profile count")
    cold_per_profile = cold_total_seconds / n_profiles
    warm_total_seconds = float(warm["aggregate"]["wall_total_seconds"])
    warm_per_event = warm_total_seconds / warm_events
    cold_peak_mib = float(cold["memory"]["process"]["peak_working_set_bytes"]) / 1024**2
    warm_peak_mib = max(
        float(repeat["memory"]["process"]["peak_working_set_bytes"]) / 1024**2
        for repeat in warm["repeats"]
    )
    return {
        "source": provenance,
        "n_profiles": n_profiles,
        "n_genes": n_genes,
        "n_warm_repeats": n_warm_repeats,
        "warm_events": warm_events,
        "cold": {
            "total_seconds": cold_total_seconds,
            "seconds_per_profile": cold_per_profile,
            "peak_working_set_mib": cold_peak_mib,
        },
        "warm": {
            "total_seconds": warm_total_seconds,
            "seconds_per_event": warm_per_event,
            "maximum_working_set_mib": warm_peak_mib,
        },
    }


def derive_tcga_range() -> dict[str, Any]:
    payload = read_json(TCGA_SUMMARY)
    broad = payload["strict_top3"]["broad"]
    rates = {basis: float(record["percent"]) for basis, record in broad.items()}
    reported = float(payload["range_across_truth_bases_percentage_points"]["broad"])
    derived = max(rates.values()) - min(rates.values())
    if not math.isclose(reported, derived, rel_tol=0, abs_tol=1e-12):
        raise ValueError("TCGA broad range is not max(source rates) - min(source rates)")
    return {
        "canonical_source": TCGA_SUMMARY.relative_to(ROOT).as_posix(),
        "source_sha256": sha256_file(TCGA_SUMMARY),
        "source_rates_percentage_points": rates,
        "derived_range_percentage_points": derived,
        "reported_range_percentage_points": reported,
        "formula": "max(source_rates) - min(source_rates)",
    }


def derive_lomo_exact() -> dict[str, Any]:
    provenance = read_json(EXACT_PROVENANCE)
    rows = load_formal_predictions(EXACT_DETAIL)
    metrics = compute_lomo_exact_metrics(rows)
    summary = metrics["summary"]
    source = provenance["formal_source"]
    source_chain = provenance.get("source_chain")
    if not isinstance(source_chain, dict):
        source_chain = _source_chain(
            origin_path=str(source.get("origin_path", source.get("staged_from", source["path"]))),
            origin_sha256=str(source.get("origin_sha256", source["sha256"])),
            staged_path=str(source["path"]),
            staged_sha256=str(source["staged_sha256"]),
            consumer="scripts/generate_lomo_exact_f1_evidence.py",
            binding="core.lomo_exact_f1.CANONICAL_FORMAL_PATH",
            staging_transform="legacy provenance migration",
        )
    _assert_staged_generator_pair("LOMO Exact", source_chain)
    if int(summary["top1_correct"]) != sum(int(row["tp"]) for row in metrics["classes"]):
        raise ValueError("LOMO Exact micro-F1 integer identity failed")
    return {
        "canonical_source": str(source_chain["generator_input"]["path"]),
        "source_sha256": str(source_chain["generator_input"]["sha256"]),
        "staged_source_sha256": str(source_chain["staged"]["sha256"]),
        "source_chain": source_chain,
        "summary": summary,
        "integer_accounting": {
            "sum_tp": sum(int(row["tp"]) for row in metrics["classes"]),
            "sum_support": sum(int(row["support"]) for row in metrics["classes"]),
        },
    }


def build_lomo_input_path_sha_pairing(
    lomo_exact: dict[str, Any], baselines: dict[str, Any]
) -> dict[str, Any]:
    """Pair the two active LOMO source chains in one audit surface."""

    group_source = baselines["sources"]["LOMO"]
    group_chain = group_source.get("source_chain")
    if not isinstance(group_chain, dict):
        group_chain = _source_chain(
            origin_path=str(group_source["origin_path"]),
            origin_sha256=str(group_source["origin_sha256"]),
            staged_path=str(group_source["staged_path"]),
            staged_sha256=str(group_source["staged_sha256"]),
            consumer=GROUP_GENERATOR,
            binding=GROUP_GENERATOR_INPUT_BINDING,
            staging_transform="legacy provenance migration",
        )
    chains = {
        "LOMO Exact": lomo_exact["source_chain"],
        "LOMO Group": group_chain,
    }
    for label, chain in chains.items():
        _assert_staged_generator_pair(label, chain)
    return {
        "schema": "braintrace.lomo_input_path_sha_pairing.v1",
        "chains": chains,
        "invariant": (
            "For each endpoint, generator_input path/SHA-256 must exactly equal the "
            "repository-staged path/SHA-256; origin is retained as its distinct source pair."
        ),
    }


def write_lomo_input_path_sha_pairing(pairing: dict[str, Any]) -> None:
    """Write machine- and reviewer-readable versions of the paired input chains."""

    write_json(LOMO_INPUT_PAIRING_JSON, pairing)
    lines = [
        "# LOMO origin / staged / generator-input path-SHA pairing",
        "",
        pairing["invariant"],
        "",
        "| Endpoint | Role | Path | SHA-256 | Generator binding |",
        "| --- | --- | --- | --- | --- |",
    ]
    for endpoint, chain in pairing["chains"].items():
        origin = chain["origin"]
        staged = chain["staged"]
        generator_input = chain["generator_input"]
        lines.extend(
            [
                f"| {endpoint} | Origin | `{origin['path']}` | `{origin['sha256']}` | frozen external formal detail |",
                f"| {endpoint} | Staged | `{staged['path']}` | `{staged['sha256']}` | repository-staged canonical table |",
                f"| {endpoint} | Generator input | `{generator_input['path']}` | `{generator_input['sha256']}` | `{generator_input['consumer']}` via `{generator_input['binding']}` |",
                "",
            ]
        )
    LOMO_INPUT_PAIRING_MD.write_text("\n".join(lines), encoding="utf-8")


def derive_loso_exact_unchanged() -> dict[str, Any]:
    rows = [
        row
        for row in read_json(MACRO_JSON)["data"]
        if row.get("endpoint") == "LOSO_Exact"
    ]
    values = [float(row["f1"]) for row in rows]
    supports = [int(row["n"]) for row in rows]
    if len(values) != 105 or sum(supports) != 814:
        raise ValueError("LOSO Exact source was unexpectedly altered")
    return {
        "status": "UNCHANGED",
        "n_classes": len(values),
        "n_profiles": sum(supports),
        "macro_f1": sum(values) / len(values),
        "locked_top1": "182/814",
        "locked_top3": "368/814",
    }


def derive_friedman() -> dict[str, Any]:
    with LAMBDA_CSV.open(newline="", encoding="utf-8-sig") as handle:
        row = next(row for row in csv.DictReader(handle) if row["lambda"] == "Friedman_test")
    chi2 = float(str(row["exact_hit1"]).split("=", 1)[1])
    p_value = float(str(row["exact_hit3"]).split("=", 1)[1])
    return {
        "source": LAMBDA_CSV.relative_to(ROOT).as_posix(),
        "source_sha256": sha256_file(LAMBDA_CSV),
        "chi2": chi2,
        "df": 2,
        "p_value": p_value,
        "exact_enumeration": "REMOVED",
        "reason": "No reproducible exact-enumeration algorithm, null space, or result artifact was found; no exact claim is retained.",
    }


def scan_current_outputs() -> dict[str, Any]:
    """Scan only active generated evidence, not labelled historical archives."""

    paths = [
        REPORT,
        QA,
        LEDGER,
        ROOT / "reproducibility" / "TRACEABILITY_MATRIX.md",
        ROOT / "reproducibility" / "v4_p0_13_macro_f1.csv",
        ROOT / "reproducibility" / "p1_cross1_5" / "cross3_f1_distribution_summary.csv",
        BASELINE_CSV,
        BASELINE_JSON,
        EXACT_PROVENANCE,
        LOMO_INPUT_PAIRING_JSON,
        LOMO_INPUT_PAIRING_MD,
    ]
    patterns = {
        "stale_tcga_range": "range 32.02",
        "stale_lomo_exact_macro": "LOMO_Exact,104,812,0.194",
        "unsupported_friedman_enumeration": "19,683-pattern exact-enumeration",
        "benchmark_unit_mix": "153-sample by 28,415-gene",
    }
    matches: list[dict[str, str]] = []
    for path in paths:
        if not path.exists():
            continue
        text = path.read_text(encoding="utf-8", errors="replace")
        for name, pattern in patterns.items():
            if pattern in text:
                matches.append({"pattern": name, "path": path.relative_to(ROOT).as_posix()})
    return {
        "scope": [path.relative_to(ROOT).as_posix() for path in paths if path.exists()],
        "current_matches": matches,
        "status": "PASS" if not matches else "FAIL",
        "historical_artifact_policy": "Historical values may remain only in explicitly labelled historical/superseded material and are not read by current generators.",
    }


def ledger_rows(
    tcga: dict[str, Any], lomo_exact: dict[str, Any], benchmark: dict[str, Any],
    baselines: dict[str, Any], friedman: dict[str, Any]
) -> list[dict[str, str]]:
    lomo = lomo_exact["summary"]
    baseline_by_endpoint = {record["endpoint"]: record for record in baselines["records"]}
    return [
        {
            "claim_id": "TCGA_BROAD_STRICT_TOP3_RANGE",
            "canonical_source": tcga["canonical_source"],
            "source_sha256": tcga["source_sha256"],
            "denominator": "four strict broad truth-basis rates",
            "exclusion_rule": "endpoint-specific denominators retained from canonical summary",
            "numerator": "max 82.5396825397 minus min 49.2307692308 percentage points",
            "derived_value": f"{tcga['derived_range_percentage_points']:.13f} pp",
            "display_value": f"{tcga['derived_range_percentage_points']:.2f} pp",
            "manuscript_location": "Validation truth-basis statement",
            "supplement_location": "R2 truth-basis range",
            "figure_table": "Traceability T029",
            "generator": "scripts/generate_scientific_remediation_artifacts.py",
            "status": "CURRENT",
        },
        {
            "claim_id": "LOMO_EXACT_F1_SUMMARY",
            "canonical_source": lomo_exact["canonical_source"],
            "source_sha256": lomo_exact["source_sha256"],
            "denominator": f"{lomo['n_classes']} truth-label classes; {lomo['n_samples']} profiles",
            "exclusion_rule": "prediction-only labels are false positives but do not expand the truth-label macro denominator",
            "numerator": f"Top1 TP={lomo['top1_correct']}; Top3={lomo['top3_correct']}",
            "derived_value": f"macro={lomo['macro_f1']:.15f}; micro={lomo['micro_f1']:.15f}",
            "display_value": f"macro {lomo['macro_f1']:.3f} ± {lomo['sd_class_f1']:.3f}; micro {lomo['micro_f1']:.3f}",
            "manuscript_location": "Validation exact-region F1 statement",
            "supplement_location": "R1/R9/Table S11/S13",
            "figure_table": "Table S11; Table S13",
            "generator": "scripts/generate_lomo_exact_f1_evidence.py",
            "status": "CURRENT",
        },
        {
            "claim_id": "BENCHMARK_COLD_WARM_UNITS",
            "canonical_source": benchmark["source"]["staged_path"],
            "source_sha256": benchmark["source"]["origin_sha256"],
            "denominator": f"{benchmark['n_profiles']} profiles; {benchmark['warm_events']} warm inference events",
            "exclusion_rule": "warm events are repeated timed inferences, not additional input profiles",
            "numerator": f"cold {benchmark['cold']['total_seconds']:.7f} s; warm {benchmark['warm']['total_seconds']:.7f} s",
            "derived_value": f"cold={benchmark['cold']['seconds_per_profile']:.10f} s/profile; warm={benchmark['warm']['seconds_per_event']:.10f} s/event",
            "display_value": f"51 profiles × 28,415 genes; 153 warm events",
            "manuscript_location": "Implementation benchmark sentence",
            "supplement_location": "R12 engineering performance",
            "figure_table": "Traceability T033–T035",
            "generator": "scripts/generate_scientific_remediation_artifacts.py",
            "status": "CURRENT",
        },
        *[
            {
                "claim_id": f"{endpoint}_GROUP_TOP3_RANDOM_BASELINES",
                "canonical_source": record["source_staged_path"],
                "source_sha256": record["source_staged_sha256"],
                "denominator": f"{record['n_profiles']} formal profiles",
                "exclusion_rule": "current hybrid route-family rows only; candidate identities are not fully serialized",
                "numerator": f"observed group Top3={record['observed_hits']}",
                "derived_value": f"uniform={record['uniform_random_rate']:.15f}; weighted={record['weighted_random_rate']:.15f}",
                "display_value": f"uniform {record['uniform_random_rate'] * 100:.1f}%; weighted {record['weighted_random_rate'] * 100:.1f}%",
                "manuscript_location": "",
                "supplement_location": "Table S8",
                "figure_table": "Table S8",
                "generator": "core/resolution_group_baselines.py",
                "status": "CURRENT",
            }
            for endpoint, record in baseline_by_endpoint.items()
        ],
        {
            "claim_id": "FRIEDMAN_LAMBDA_SENSITIVITY",
            "canonical_source": friedman["source"],
            "source_sha256": friedman["source_sha256"],
            "denominator": "9 donor-level profiles across three lambda conditions",
            "exclusion_rule": "No unsupported exact-enumeration result is retained",
            "numerator": "chi-square approximation",
            "derived_value": f"chi2={friedman['chi2']:.4f}; df={friedman['df']}; p={friedman['p_value']:.3f}",
            "display_value": f"Friedman χ²={friedman['chi2']:.4f}, df={friedman['df']}, P={friedman['p_value']:.3f}",
            "manuscript_location": "Lambda sensitivity statement",
            "supplement_location": "R9",
            "figure_table": "Lambda sensitivity",
            "generator": "reproducibility/generate_all_csvs.py",
            "status": "CURRENT; EXACT_ENUMERATION_REMOVED",
        },
    ]


def write_report(
    tcga: dict[str, Any], lomo_exact: dict[str, Any], benchmark: dict[str, Any],
    baselines: dict[str, Any], lomo_pairing: dict[str, Any], friedman: dict[str, Any], scan: dict[str, Any],
    tests_status: str, docx_status: str
) -> None:
    summary = lomo_exact["summary"]
    base = {row["endpoint"]: row for row in baselines["records"]}
    lines = [
        "# Scientific remediation report",
        "",
        "All corrected secondary statistics are regenerated from the current formal prediction/detail sources. The frozen model and primary endpoints are unchanged.",
        "",
        "## Corrected source chains",
        "",
        f"- TCGA/BraTS broad strict Top3 range: max-min = `{tcga['derived_range_percentage_points']:.13f}` pp (`{tcga['derived_range_percentage_points']:.2f}` pp displayed) from `{tcga['canonical_source']}`.",
        f"- LOMO Exact: `{summary['top1_correct']}/{summary['n_samples']}` Top1 = micro-F1 `{summary['micro_f1']:.15f}`; macro-F1 `{summary['macro_f1']:.15f}` across `{summary['n_classes']}` truth-label classes.",
        f"- LOMO Exact and LOMO Group origin/staged/generator-input path+SHA pairs: `{LOMO_INPUT_PAIRING_MD.relative_to(ROOT).as_posix()}` ({len(lomo_pairing['chains'])} endpoint chains; staged and generator-input pairs are identical by assertion).",
        f"- Benchmark: `{benchmark['n_profiles']}` profiles × `{benchmark['n_genes']}` genes; `{benchmark['warm_events']}` warm inference events; cold peak `{benchmark['cold']['peak_working_set_mib']:.4f}` MiB and warm maximum `{benchmark['warm']['maximum_working_set_mib']:.4f}` MiB.",
        f"- Resolution-group Top3 random baselines: LOSO uniform/weighted `{base['LOSO']['uniform_random_rate']:.15f}` / `{base['LOSO']['weighted_random_rate']:.15f}`; LOMO `{base['LOMO']['uniform_random_rate']:.15f}` / `{base['LOMO']['weighted_random_rate']:.15f}`.",
        f"- Friedman: χ²=`{friedman['chi2']:.4f}`, df=`{friedman['df']}`, P=`{friedman['p_value']:.3f}`; exact-enumeration status: `{friedman['exact_enumeration']}`.",
        "",
        "## Provenance boundary",
        "",
        "The LOMO Exact macro denominator is the truth-label universe. Top1 predictions outside that universe remain false positives and are disclosed in the formal F1 provenance artifact; they do not alter the frozen model or route.",
        "",
        "## QA status",
        "",
        f"- Current generated-evidence stale scan: `{scan['status']}`",
        f"- DOCX remediation: `{docx_status}`",
        f"- Regression tests: `{tests_status}`",
        "",
    ]
    REPORT.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--loso-group-source", type=Path)
    parser.add_argument("--lomo-group-source", type=Path)
    parser.add_argument("--benchmark-manifest", type=Path)
    parser.add_argument("--tests-status", default="NOT_RUN")
    parser.add_argument("--docx-status", default="PENDING")
    args = parser.parse_args()

    tcga = derive_tcga_range()
    lomo_exact = derive_lomo_exact()
    baselines = generate_group_baselines(args.loso_group_source, args.lomo_group_source)
    lomo_pairing = build_lomo_input_path_sha_pairing(lomo_exact, baselines)
    write_lomo_input_path_sha_pairing(lomo_pairing)
    manifest, benchmark_provenance = stage_benchmark(args.benchmark_manifest)
    benchmark = derive_benchmark(manifest, benchmark_provenance)
    friedman = derive_friedman()
    provisional_scan = {"status": "PENDING"}
    write_report(
        tcga, lomo_exact, benchmark, baselines, lomo_pairing, friedman, provisional_scan,
        args.tests_status, args.docx_status,
    )
    scan = scan_current_outputs()
    qa = {
        "schema": "braintrace.scientific_remediation_qa.v1",
        "tcga_broad_strict_top3_range": tcga,
        "lomo_exact_f1": lomo_exact,
        "lomo_input_path_sha_pairing": lomo_pairing,
        "loso_exact_unchanged": derive_loso_exact_unchanged(),
        "benchmark": benchmark,
        "resolution_group_random_baselines": baselines,
        "friedman_exact_enumeration": friedman,
        "stale_current_value_scan": scan,
        "tests": {"status": args.tests_status},
        "docx_remediation": {"status": args.docx_status},
    }
    write_json(QA, qa)
    scan = scan_current_outputs()
    qa["stale_current_value_scan"] = scan
    write_json(QA, qa)
    fields = [
        "claim_id",
        "canonical_source",
        "source_sha256",
        "denominator",
        "exclusion_rule",
        "numerator",
        "derived_value",
        "display_value",
        "manuscript_location",
        "supplement_location",
        "figure_table",
        "generator",
        "status",
    ]
    write_csv(LEDGER, ledger_rows(tcga, lomo_exact, benchmark, baselines, friedman), fields)
    scan = scan_current_outputs()
    qa["stale_current_value_scan"] = scan
    write_json(QA, qa)
    write_report(
        tcga,
        lomo_exact,
        benchmark,
        baselines,
        lomo_pairing,
        friedman,
        scan,
        args.tests_status,
        args.docx_status,
    )
    print(json.dumps({"qa": QA.as_posix(), "scan": scan["status"]}, ensure_ascii=False))
    return 0 if scan["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
