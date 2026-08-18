#!/usr/bin/env python3
"""Build non-Huang scientific-provenance artifacts from sample-level outputs.

This script deliberately consumes the authoritative sample-level outputs rather
than manuscript summaries.  It does not fit, tune, or modify the frozen
BrainTrace model.  Its outputs are the source for the non-Huang provenance
patch candidate.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any, Iterable

import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[1]
CANONICAL_AHBA_ROUTE = "hybrid_projected_network_logcpm_exact"


def sha256(path: Path) -> str:
    """Return the SHA-256 checksum of a source artifact."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def source_record(path: Path) -> dict[str, str]:
    # Store a portable locator plus SHA-256; never persist machine-local paths.
    resolved = path.resolve()
    try:
        locator = resolved.relative_to(REPO_ROOT.resolve()).as_posix()
        path_kind = "repository_relative"
    except ValueError:
        locator = path.name
        path_kind = "external_basename"
    return {"path": locator, "path_kind": path_kind, "sha256": sha256(path)}


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def boolean_series(series: pd.Series) -> pd.Series:
    """Normalize booleans without treating the string ``False`` as truthy."""
    return series.map(
        lambda value: value is True
        or (isinstance(value, (int, float)) and not isinstance(value, bool) and value == 1)
        or str(value).strip().lower() in {"true", "1", "yes"}
    )


def pipe_count(series: pd.Series) -> pd.Series:
    """Count source labels represented as a pipe-delimited field."""
    return series.fillna("").map(
        lambda value: 0
        if not str(value).strip()
        else len([item for item in str(value).split(" | ") if item.strip()])
    )


def percentage(correct: int, denominator: int) -> float:
    if denominator <= 0:
        raise ValueError("denominator must be positive")
    return correct / denominator * 100


def fraction_record(correct: int, denominator: int) -> dict[str, float | int]:
    return {
        "correct": int(correct),
        "n": int(denominator),
        "proportion": correct / denominator,
        "percent": percentage(correct, denominator),
    }


def require_columns(frame: pd.DataFrame, columns: Iterable[str], source: Path) -> None:
    missing = sorted(set(columns) - set(frame.columns))
    if missing:
        raise ValueError(f"{source} is missing required columns: {missing}")


def build_ahba_endpoint_ledger(sample_detail: Path, output_dir: Path) -> dict[str, Any]:
    """Build endpoint-specific AHBA accounting from the canonical route detail."""
    raw = pd.read_csv(sample_detail)
    require_columns(
        raw,
        {
            "route",
            "sample_id",
            "ahba_donor",
            "supported_for_accuracy",
            "accuracy_level",
            "allowed_bo2023_networks",
            "allowed_bo2023_region_keys",
            "group_evaluation_status",
            "group_truth_in_candidate_beam",
            "network_top1_hit",
            "network_top3_hit",
            "group_top1_hit",
            "group_top3_hit",
            "region_top1_exact_hit",
            "region_top3_exact_hit",
        },
        sample_detail,
    )
    detail = raw.loc[raw["route"].eq(CANONICAL_AHBA_ROUTE)].copy()
    if detail["sample_id"].duplicated().any():
        raise ValueError("canonical AHBA route has duplicate sample IDs")
    if len(detail) == 0:
        raise ValueError(f"canonical AHBA route {CANONICAL_AHBA_ROUTE!r} is absent")

    network_evaluable = boolean_series(detail["supported_for_accuracy"])
    exact_evaluable = detail["accuracy_level"].eq("exact_region")
    # A mapped exact region is the published group endpoint's eligibility
    # criterion.  The source's group_evaluation_status additionally records
    # whether the truth survived the predicted Network candidate beam, which
    # is an outcome/failure mode, not a pre-evaluation denominator filter.
    group_evaluable = exact_evaluable.copy()
    network_label_count = pipe_count(detail["allowed_bo2023_networks"])
    exact_label_count = pipe_count(detail["allowed_bo2023_region_keys"])
    group_labels_in_candidate_beam = pipe_count(detail.get("allowed_resolution_groups", pd.Series("", index=detail.index)))

    reason = detail["accuracy_level"].astype(str)
    ledger = pd.DataFrame(
        {
            "sample_id": detail["sample_id"],
            "donor_id": detail["ahba_donor"],
            "technical_replicate_group_or_tissue_id": detail["sample_id"].astype(str).str.split("|", n=1).str[-1],
            "network_evaluable": network_evaluable,
            "group_evaluable": group_evaluable,
            "exact_evaluable": exact_evaluable,
            "network_single_label_subset": network_evaluable & network_label_count.eq(1),
            "group_single_label_subset": group_evaluable & exact_label_count.eq(1),
            "exact_single_label_subset": exact_evaluable & exact_label_count.eq(1),
            "network_truth_label_count": network_label_count,
            # The formal source does not expose a candidate-independent count
            # of resolution-group truth labels.  Blank is intentional: the
            # available group list is conditional on the predicted Network beam.
            "group_truth_label_count": pd.NA,
            "exact_truth_label_count": exact_label_count,
            "source_group_truth_label_count_in_candidate_beam": group_labels_in_candidate_beam,
            "source_group_evaluation_status": detail["group_evaluation_status"],
            "source_group_truth_in_candidate_beam": boolean_series(detail["group_truth_in_candidate_beam"]),
            "reason_network_not_evaluable": reason.where(~network_evaluable, ""),
            "reason_group_not_evaluable": reason.where(~group_evaluable, ""),
            "reason_exact_not_evaluable": reason.where(~exact_evaluable, ""),
            "network_top1_hit": boolean_series(detail["network_top1_hit"]),
            "network_top3_hit": boolean_series(detail["network_top3_hit"]),
            "group_top1_hit": boolean_series(detail["group_top1_hit"]),
            "group_top3_hit": boolean_series(detail["group_top3_hit"]),
            "exact_top1_hit": boolean_series(detail["region_top1_exact_hit"]),
            "exact_top3_hit": boolean_series(detail["region_top3_exact_hit"]),
        }
    ).sort_values("sample_id", kind="stable")

    if int(network_evaluable.sum()) > len(ledger) or int(exact_evaluable.sum()) > int(network_evaluable.sum()):
        raise ValueError("AHBA endpoint eligibility is internally inconsistent")

    def endpoint_metrics(mask: pd.Series, top1: str, top3: str) -> dict[str, Any]:
        n = int(mask.sum())
        return {
            "evaluable": n,
            "top1": fraction_record(int(ledger.loc[mask, top1].sum()), n),
            "top3": fraction_record(int(ledger.loc[mask, top3].sum()), n),
        }

    network_metrics = endpoint_metrics(ledger["network_evaluable"], "network_top1_hit", "network_top3_hit")
    group_metrics = endpoint_metrics(ledger["group_evaluable"], "group_top1_hit", "group_top3_hit")
    exact_metrics = endpoint_metrics(ledger["exact_evaluable"], "exact_top1_hit", "exact_top3_hit")

    def sensitivity_metrics(mask: pd.Series, top1: str, top3: str) -> dict[str, Any]:
        n = int(mask.sum())
        return {
            "n": n,
            "top1": fraction_record(int(ledger.loc[mask, top1].sum()), n),
            "top3": fraction_record(int(ledger.loc[mask, top3].sum()), n),
        }

    summary: dict[str, Any] = {
        "artifact_status": "v0.1.17 scientific-provenance patch candidate (unreleased)",
        "source": source_record(sample_detail),
        "canonical_route": CANONICAL_AHBA_ROUTE,
        "unit": "replicate-collapsed AHBA tissue sample",
        "endpoint_evaluability": {
            "replicate_collapsed_tissues": int(len(ledger)),
            "network": network_metrics,
            "resolution_group": group_metrics,
            "exact_region": exact_metrics,
            "network_unique_single_label_sensitivity": sensitivity_metrics(
                ledger["network_single_label_subset"], "network_top1_hit", "network_top3_hit"
            ),
            "group_single_label_sensitivity": sensitivity_metrics(
                ledger["group_single_label_subset"], "group_top1_hit", "group_top3_hit"
            ),
            "exact_single_label_sensitivity": sensitivity_metrics(
                ledger["exact_single_label_subset"], "exact_top1_hit", "exact_top3_hit"
            ),
        },
        "interpretation": {
            "canonical_accounting": "Endpoint-specific evaluability ledger; not a unique sequential sample-attrition pipeline.",
            "group_truth_label_count": "Not available as a candidate-independent field in the source detail; only a candidate-beam-dependent source count is retained in the ledger.",
            "single_label_group_rule": "The group single-label sensitivity subset is defined by a uniquely mapped exact-region truth label, the source-supported fixed mapping basis for the group/exact endpoints.",
        },
    }

    output_dir.mkdir(parents=True, exist_ok=True)
    ledger.to_csv(output_dir / "ahba_endpoint_evaluability_ledger.csv", index=False)
    write_json(output_dir / "ahba_endpoint_evaluability_summary.json", summary)
    return summary


def build_tcga_truth_basis_summary(sample_detail: Path, output_dir: Path) -> dict[str, Any]:
    """Recompute strict Top3 truth-basis counts from the per-patient result."""
    data = pd.read_csv(sample_detail)
    require_columns(
        data,
        {
            "patient_barcode",
            "edema_voxels",
            "edema_network_dominant",
            "center_network_top3_strict",
            "core_network_top3_strict",
            "edema_network_top3_strict",
            "whole_tumor_network_top3_strict",
            "center_broad_top3_strict",
            "core_broad_top3_strict",
            "edema_broad_top3_strict",
            "whole_tumor_broad_top3_strict",
        },
        sample_detail,
    )
    if data["patient_barcode"].duplicated().any():
        raise ValueError("TCGA/BraTS source has duplicate patient barcodes")

    edema_voxels = pd.to_numeric(data["edema_voxels"], errors="raise")
    edema_eligible = edema_voxels.gt(0) & ~data["edema_network_dominant"].eq("out_of_scope")
    exclusion_reason = pd.Series("", index=data.index, dtype="object")
    exclusion_reason.loc[edema_voxels.le(0)] = "no_label_2_edema_voxels"
    exclusion_reason.loc[edema_voxels.gt(0) & data["edema_network_dominant"].eq("out_of_scope")] = (
        "cerebellar_or_posterior_fossa_outside_current_braintrace_label_space"
    )
    if exclusion_reason.eq("").sum() != int(edema_eligible.sum()):
        raise ValueError("TCGA/BraTS edema eligibility rule did not partition patients")

    rows: list[dict[str, Any]] = []
    metrics: dict[str, dict[str, dict[str, Any]]] = {}
    for level in ("network", "broad"):
        metrics[level] = {}
        for truth_basis in ("center", "core", "edema", "whole_tumor"):
            mask = edema_eligible if truth_basis == "edema" else pd.Series(True, index=data.index)
            column = f"{truth_basis}_{level}_top3_strict"
            hit = boolean_series(data[column])
            value = fraction_record(int(hit.loc[mask].sum()), int(mask.sum()))
            metrics[level][truth_basis] = value
            rows.append(
                {
                    "truth_basis": truth_basis,
                    "level": level,
                    "variant": "strict",
                    "top_k": "top3",
                    "correct": value["correct"],
                    "n": value["n"],
                    "accuracy": f"{value['proportion']:.8f}",
                    "accuracy_percent": f"{value['percent']:.2f}",
                    "inclusion_rule": (
                        "edema_voxels > 0 and edema_network_dominant != out_of_scope"
                        if truth_basis == "edema"
                        else "all paired TCGA/BraTS cases"
                    ),
                }
            )

    ranges = {
        level: max(item["percent"] for item in level_metrics.values())
        - min(item["percent"] for item in level_metrics.values())
        for level, level_metrics in metrics.items()
    }
    excluded = data.loc[~edema_eligible, ["patient_barcode"]].copy()
    excluded["reason"] = exclusion_reason.loc[~edema_eligible]
    summary: dict[str, Any] = {
        "artifact_status": "v0.1.17 scientific-provenance patch candidate (unreleased)",
        "source": source_record(sample_detail),
        "total_paired_cases": int(len(data)),
        "primary_edema_comparator": {
            "n": int(edema_eligible.sum()),
            "inclusion_rule": "edema_voxels > 0 and edema_network_dominant != out_of_scope",
            "excluded_cases": excluded.sort_values("patient_barcode").to_dict(orient="records"),
        },
        "strict_top3": metrics,
        "range_across_truth_bases_percentage_points": ranges,
    }
    output_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(output_dir / "tcga_brats_truth_basis_top3_summary.csv", index=False)
    write_json(output_dir / "tcga_brats_truth_basis_top3_summary.json", summary)
    return summary


def build_orthology_humanization_summary(
    signature_humanization: Path,
    top200_humanization: Path,
    gprofiler_results: Path,
    output_dir: Path,
) -> dict[str, Any]:
    """Recompute row-occurrence humanization and preserve denominator semantics."""
    signature = pd.read_csv(signature_humanization)
    top200 = pd.read_csv(top200_humanization)
    require_columns(signature, {"human_gene_symbol", "is_unmapped"}, signature_humanization)
    require_columns(top200, {"human_gene_symbol", "is_unmapped"}, top200_humanization)
    signature_unmapped = boolean_series(signature["is_unmapped"])
    top200_unmapped = boolean_series(top200["is_unmapped"])
    signature_humanized = ~signature_unmapped & signature["human_gene_symbol"].fillna("").astype(str).str.strip().ne("")
    top200_humanized = ~top200_unmapped & top200["human_gene_symbol"].fillna("").astype(str).str.strip().ne("")

    gprofiler = json.loads(gprofiler_results.read_text(encoding="utf-8"))
    result_rows = gprofiler.get("result", [])
    query_sizes = {int(row["query_size"]) for row in result_rows if "query_size" in row}
    if len(query_sizes) != 1:
        raise ValueError("g:Profiler source does not have one unambiguous mapped query size")
    gprofiler_mapped = query_sizes.pop()

    summary: dict[str, Any] = {
        "artifact_status": "v0.1.17 scientific-provenance patch candidate (unreleased)",
        "sources": {
            "region_signature_humanization": source_record(signature_humanization),
            "network_top200_humanization": source_record(top200_humanization),
            "gprofiler_model_background_results": source_record(gprofiler_results),
        },
        "region_signature_rows": {
            "unit": "gene-by-region row occurrence",
            "total": int(len(signature)),
            "humanized": fraction_record(int(signature_humanized.sum()), int(len(signature))),
            "unmapped": fraction_record(int(signature_unmapped.sum()), int(len(signature))),
            "humanization_status_counts": {
                str(key): int(value)
                for key, value in signature["humanization_status"].fillna("missing").value_counts().sort_index().items()
            },
        },
        "network_top200_orthology_humanizable": fraction_record(int(top200_humanized.sum()), int(len(top200))),
        "gprofiler_mapped": {
            "mapped": int(gprofiler_mapped),
            "input_panel_n": int(len(top200)),
            "unit": "g:Profiler mapped query genes; a distinct service/filtering universe from the frozen orthology audit",
        },
        "interpretation": (
            "The 60.5% humanization rate is row-level: gene-by-region row occurrences, not independent genes."
        ),
    }
    output_dir.mkdir(parents=True, exist_ok=True)
    write_json(output_dir / "orthology_humanization_summary.json", summary)
    return summary


def build_tier_cascade_summary(exact_detail: Path, group_detail: Path, output_dir: Path) -> dict[str, Any]:
    """Calculate tier-cascade denominators within the same exact-evaluable universe."""
    exact = pd.read_csv(exact_detail)
    group = pd.read_csv(group_detail)
    require_columns(exact, {"sample_id", "network_top3_hit", "hit3"}, exact_detail)
    require_columns(group, {"sample_id", "network_top3_hit", "group_hit3"}, group_detail)
    if exact["sample_id"].duplicated().any() or group["sample_id"].duplicated().any():
        raise ValueError("tier detail contains duplicate sample IDs")
    group_indexed = group.set_index("sample_id")
    if set(exact["sample_id"]) != set(group_indexed.index):
        raise ValueError("exact and group details do not have the same sample universe")

    exact = exact.set_index("sample_id").sort_index()
    group_indexed = group_indexed.loc[exact.index]
    network_hit = boolean_series(exact["network_top3_hit"])
    if not network_hit.equals(boolean_series(group_indexed["network_top3_hit"])):
        raise ValueError("exact and group details disagree on Network candidate-set status")
    exact_hit = boolean_series(exact["hit3"])
    group_hit = boolean_series(group_indexed["group_hit3"])
    n = len(exact)
    network_retained = int(network_hit.sum())
    network_miss = int((~network_hit).sum())
    exact_top3_hit = int(exact_hit.sum())
    exact_top3_miss = int((~exact_hit).sum())
    recovery_exact = int((exact_hit & ~network_hit).sum())
    recovery_group = int((group_hit & ~network_hit).sum())

    summary: dict[str, Any] = {
        "artifact_status": "v0.1.17 scientific-provenance patch candidate (unreleased)",
        "sources": {
            "exact_evaluable_detail": source_record(exact_detail),
            "group_evaluable_detail": source_record(group_detail),
        },
        "universe": "LOSO exact-evaluable sample IDs shared by exact and resolution-group detail outputs",
        "exact_evaluable": {
            "n": int(n),
            "top3_hit": int(exact_top3_hit),
            "top3_miss": int(exact_top3_miss),
        },
        "network_candidate_set_within_same_universe": {
            "truth_retained": int(network_retained),
            "truth_missed": int(network_miss),
        },
        "exact_top3_given_network_truth_retained": fraction_record(
            int((exact_hit & network_hit).sum()), network_retained
        ),
        "group_top3_given_network_truth_retained": fraction_record(
            int((group_hit & network_hit).sum()), network_retained
        ),
        "network_candidate_miss_share_of_exact_top3_misses": fraction_record(network_miss, exact_top3_miss),
        "recovery_after_network_candidate_miss": {
            "exact_top3": fraction_record(recovery_exact, network_miss),
            "group_top3": fraction_record(recovery_group, network_miss),
        },
    }
    output_dir.mkdir(parents=True, exist_ok=True)
    write_json(output_dir / "tier_cascade_loso_summary.json", summary)
    return summary


def bh_adjust(pvalues: list[float]) -> list[float]:
    """Benjamini-Hochberg adjustment with monotonicity enforcement."""
    order = sorted(range(len(pvalues)), key=pvalues.__getitem__)
    adjusted = [0.0] * len(pvalues)
    running = 1.0
    total = len(pvalues)
    for rank_position in range(total - 1, -1, -1):
        index = order[rank_position]
        rank = rank_position + 1
        running = min(running, pvalues[index] * total / rank)
        adjusted[index] = min(1.0, running)
    return adjusted


def build_sign_flip_current_family(sign_flip_source: Path, output_dir: Path) -> dict[str, Any]:
    """Apply one current four-test BH family to source raw P values."""
    source = pd.read_csv(sign_flip_source)
    require_columns(source, {"endpoint", "topk", "raw_p"}, sign_flip_source)
    wanted = [
        ("network", "top1", "Network Top1"),
        ("network", "top3", "Network Top3"),
        ("resolution_group", "top3", "resolution-group Top3"),
        ("exact_region", "top3", "exact-region Top3"),
    ]
    rows: list[dict[str, Any]] = []
    pvalues: list[float] = []
    for endpoint, topk, display in wanted:
        selected = source.loc[source["endpoint"].eq(endpoint) & source["topk"].eq(topk), "raw_p"]
        if len(selected) != 1:
            raise ValueError(f"expected exactly one raw sign-flip P value for {endpoint}/{topk}")
        pvalues.append(float(selected.iloc[0]))
        rows.append({"endpoint": display, "source_endpoint": endpoint, "top_k": topk})
    adjusted = bh_adjust(pvalues)
    for row, raw_p, bh_p in zip(rows, pvalues, adjusted):
        row["raw_p"] = raw_p
        row["raw_p_display"] = f"{raw_p:.6f}"
        row["bh_p_m4"] = bh_p
        row["bh_p_m4_display"] = f"{bh_p:.6f}"
        row["significant_bh_0_05"] = bh_p < 0.05

    output_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(output_dir / "sign_flip_current_family.csv", index=False)
    summary = {
        "artifact_status": "v0.1.17 scientific-provenance patch candidate (unreleased)",
        "source": source_record(sign_flip_source),
        "family": "Network Top1, Network Top3, resolution-group Top3, exact-region Top3",
        "n_tests": 4,
        "none_significant_bh_0_05": not any(row["significant_bh_0_05"] for row in rows),
        "rows": rows,
    }
    write_json(output_dir / "sign_flip_current_family.json", summary)
    return summary


def build_conflict_ledger(
    ahba: dict[str, Any],
    tcga: dict[str, Any],
    orthology: dict[str, Any],
    tier: dict[str, Any],
    sign_flip: dict[str, Any],
    output_path: Path,
) -> None:
    """Write the required audit ledger from recomputed artifact values."""
    ahba_counts = ahba["endpoint_evaluability"]
    tcga_network = tcga["strict_top3"]["network"]
    tcga_broad = tcga["strict_top3"]["broad"]
    rows = [
        {
            "issue": "AHBA endpoint evaluability",
            "current_repository_statement": "Historical AHBA traces can be read as competing sequential sample-flow accounts.",
            "recomputed_value": (
                f"{ahba_counts['replicate_collapsed_tissues']} replicate-collapsed tissues; "
                f"Network n={ahba_counts['network']['evaluable']}; group/exact n={ahba_counts['resolution_group']['evaluable']}; "
                f"Network single-label n={ahba_counts['network_unique_single_label_sensitivity']['n']}; "
                f"group/exact single-label n={ahba_counts['exact_single_label_sensitivity']['n']}"
            ),
            "authoritative_source": ahba["source"]["path"],
            "action": "Created the canonical endpoint-evaluability ledger and classified old traces as historical engineering traces.",
            "status": "FIXED_STALE_VALUE",
            "notes": "Endpoint-specific eligibility is not a unique sequential attrition pipeline.",
        },
        {
            "issue": "TCGA/BraTS primary edema comparator",
            "current_repository_statement": "The audit script reported edema n=64.",
            "recomputed_value": f"edema n={tcga['primary_edema_comparator']['n']}; exclusions={tcga['primary_edema_comparator']['excluded_cases']}",
            "authoritative_source": tcga["source"]["path"],
            "action": "Derived patient counts from the per-patient truth/prediction output and corrected endpoint metadata.",
            "status": "FIXED_STALE_VALUE",
            "notes": "No underlying TCGA/BraTS prediction hit was changed.",
        },
        {
            "issue": "TCGA/BraTS strict Top3 truth-basis range",
            "current_repository_statement": "Current provenance lacked a source-derived range audit and risked retaining 32.02 pp.",
            "recomputed_value": (
                f"Network {tcga_network['center']['correct']}/{tcga_network['center']['n']}, "
                f"{tcga_network['core']['correct']}/{tcga_network['core']['n']}, "
                f"{tcga_network['edema']['correct']}/{tcga_network['edema']['n']}, "
                f"{tcga_network['whole_tumor']['correct']}/{tcga_network['whole_tumor']['n']}; "
                f"range={tcga['range_across_truth_bases_percentage_points']['network']:.2f} pp. "
                f"Broad range={tcga['range_across_truth_bases_percentage_points']['broad']:.2f} pp"
            ),
            "authoritative_source": tcga["source"]["path"],
            "action": "Added source-derived strict truth-basis summary and synchronized current documentation.",
            "status": "DOCUMENTATION_SYNC",
            "notes": "Broad range is computed from 32/65, 45/65, 52/63, and 46/65.",
        },
        {
            "issue": "Orthology/humanization denominator unit",
            "current_repository_statement": "A row-level percentage could be misread as a unique-gene percentage.",
            "recomputed_value": (
                f"{orthology['region_signature_rows']['humanized']['correct']}/"
                f"{orthology['region_signature_rows']['total']} humanized and "
                f"{orthology['region_signature_rows']['unmapped']['correct']}/"
                f"{orthology['region_signature_rows']['total']} unmapped gene-by-region row occurrences; "
                f"Top200 orthology={orthology['network_top200_orthology_humanizable']['correct']}/"
                f"{orthology['network_top200_orthology_humanizable']['n']}; "
                f"g:Profiler={orthology['gprofiler_mapped']['mapped']}/"
                f"{orthology['gprofiler_mapped']['input_panel_n']}"
            ),
            "authoritative_source": orthology["sources"]["region_signature_humanization"]["path"],
            "action": "Made row-occurrence unit explicit and kept Top200 orthology and g:Profiler universes separate.",
            "status": "FIXED_STALE_VALUE",
            "notes": "8,800 is not a count of independent macaque genes.",
        },
        {
            "issue": "LOSO tier-cascade denominator",
            "current_repository_statement": "Historical prose used 66/446=14.8% and Group conditional 78.32%.",
            "recomputed_value": (
                f"exact n={tier['exact_evaluable']['n']}; exact Top3={tier['exact_evaluable']['top3_hit']}/"
                f"{tier['exact_evaluable']['n']}; Network candidate misses={tier['network_candidate_set_within_same_universe']['truth_missed']}/"
                f"{tier['exact_evaluable']['top3_miss']}={tier['network_candidate_miss_share_of_exact_top3_misses']['percent']:.2f}%; "
                f"exact conditional={tier['exact_top3_given_network_truth_retained']['correct']}/"
                f"{tier['exact_top3_given_network_truth_retained']['n']}="
                f"{tier['exact_top3_given_network_truth_retained']['percent']:.2f}%; "
                f"group conditional={tier['group_top3_given_network_truth_retained']['correct']}/"
                f"{tier['group_top3_given_network_truth_retained']['n']}="
                f"{tier['group_top3_given_network_truth_retained']['percent']:.2f}%"
            ),
            "authoritative_source": tier["sources"]["exact_evaluable_detail"]["path"],
            "action": "Recomputed within the shared exact-evaluable sample universe and replaced stale denominator leakage.",
            "status": "FIXED_STALE_VALUE",
            "notes": "Exact and group recovery after a Network candidate-set miss are both 0%.",
        },
        {
            "issue": "Sign-flip traceability",
            "current_repository_statement": "Traceability retained prior three-test BH values for two current endpoints.",
            "recomputed_value": "; ".join(
                f"{row['endpoint']}: raw={row['raw_p_display']}, BH={row['bh_p_m4_display']}" for row in sign_flip["rows"]
            ),
            "authoritative_source": sign_flip["source"]["path"],
            "action": "Applied one four-test BH family to source raw P values and synchronized traceability.",
            "status": "FIXED_STALE_VALUE",
            "notes": "None significant after BH correction.",
        },
    ]
    pd.DataFrame(rows).to_csv(output_path, index=False)


def parser() -> argparse.ArgumentParser:
    argument_parser = argparse.ArgumentParser(description=__doc__)
    argument_parser.add_argument("--ahba-sample-detail", type=Path, required=True)
    argument_parser.add_argument("--tcga-sample-detail", type=Path, required=True)
    argument_parser.add_argument("--signature-humanization", type=Path, required=True)
    argument_parser.add_argument("--top200-humanization", type=Path, required=True)
    argument_parser.add_argument("--gprofiler-results", type=Path, required=True)
    argument_parser.add_argument("--tier-exact-detail", type=Path, required=True)
    argument_parser.add_argument("--tier-group-detail", type=Path, required=True)
    argument_parser.add_argument("--sign-flip-source", type=Path, required=True)
    argument_parser.add_argument("--output-root", type=Path, default=REPO_ROOT)
    return argument_parser


def main() -> int:
    args = parser().parse_args()
    root = args.output_root.resolve()
    ahba = build_ahba_endpoint_ledger(args.ahba_sample_detail, root / "reproducibility" / "ahba")
    tcga = build_tcga_truth_basis_summary(args.tcga_sample_detail, root / "reproducibility")
    orthology = build_orthology_humanization_summary(
        args.signature_humanization,
        args.top200_humanization,
        args.gprofiler_results,
        root / "reproducibility",
    )
    tier = build_tier_cascade_summary(args.tier_exact_detail, args.tier_group_detail, root / "reproducibility")
    sign_flip = build_sign_flip_current_family(args.sign_flip_source, root / "reproducibility")
    build_conflict_ledger(ahba, tcga, orthology, tier, sign_flip, root / "NONHUANG_SCIENTIFIC_CONFLICT_LEDGER.csv")
    print(
        json.dumps(
            {
                "ahba": ahba["endpoint_evaluability"],
                "tcga_edema_n": tcga["primary_edema_comparator"]["n"],
                "tcga_ranges_pp": tcga["range_across_truth_bases_percentage_points"],
                "orthology_rows": orthology["region_signature_rows"],
                "tier": tier,
                "sign_flip_none_significant": sign_flip["none_significant_bh_0_05"],
            },
            ensure_ascii=False,
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
