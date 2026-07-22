#!/usr/bin/env python
from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, spearmanr


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from core.bo2023_region_tracing import ROUTE_NAME, trace_bo2023_secondary_regions  # noqa: E402
from core.network_tracing import (  # noqa: E402
    DEFAULT_BO2023_NETWORK_MODEL,
    DEFAULT_BO2023_REFERENCE_PROJECTOR,
    trace_network_expression,
)


DEFAULT_INPUT = (
    ROOT
    / "external_inputs"
    / "huang2025_pmc12041490"
    / "41698_2025_909_MOESM2_ESM.csv"
)
DEFAULT_SOURCE_XLSB = DEFAULT_INPUT.with_suffix(".xlsb")
DEFAULT_OUTDIR = ROOT / "validation_runs" / "huang2025_external_candidate_20260719"
DEFAULT_DB_PATH = ROOT / "braintrace_source_tracing.db"
REGION_REFERENCE = ROOT / "data" / "models" / "bo2023_formal_region_logcpm_reference_matrix.npz"
BEAM_PANELS = ROOT / "data" / "models" / "bo2023_formal_region_beam_gene_panels.json"
SAMPLE_RE = re.compile(r"^(GLI|MEN|NOR)_(CSF|plasma)(\d+)$")
DISEASE_LABEL = {"GLI": "glioma", "MEN": "meningioma", "NOR": "control"}
SPECIMEN_LABEL = {"CSF": "CSF", "plasma": "plasma"}
EXPECTED_COUNTS = {
    ("glioma", "CSF"): 16,
    ("glioma", "plasma"): 18,
    ("meningioma", "CSF"): 43,
    ("meningioma", "plasma"): 46,
    ("control", "CSF"): 18,
    ("control", "plasma"): 18,
}
N_BOOTSTRAP = 10_000
N_PERMUTATIONS = 10_000
RANDOM_SEED = 20250719


def json_default(value: Any) -> Any:
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating,)):
        return float(value)
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, Path):
        return str(value)
    return str(value)


def write_json(path: Path, payload: Any) -> None:
    path.write_text(
        json.dumps(payload, indent=2, ensure_ascii=False, default=json_default) + "\n",
        encoding="utf-8",
    )


def markdown_table(frame: pd.DataFrame) -> str:
    if frame.empty:
        return "_No rows._"
    text = frame.astype(object).where(pd.notna(frame), "")
    columns = [str(column) for column in text.columns]
    lines = [
        "| " + " | ".join(columns) + " |",
        "| " + " | ".join("---" for _ in columns) + " |",
    ]
    for row in text.itertuples(index=False):
        values = [str(value).replace("\n", " ").replace("|", "\\|") for value in row]
        lines.append("| " + " | ".join(values) + " |")
    return "\n".join(lines)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def parse_sample_id(sample_id: str) -> dict[str, Any]:
    match = SAMPLE_RE.fullmatch(str(sample_id))
    if not match:
        raise ValueError(f"Unexpected Huang 2025 sample ID: {sample_id!r}")
    prefix, specimen_raw, number = match.groups()
    disease = DISEASE_LABEL[prefix]
    specimen = SPECIMEN_LABEL[specimen_raw]
    patient_key = f"{prefix}_{int(number):02d}"
    return {
        "sample_id": str(sample_id),
        "disease": disease,
        "tumor_status": "control" if disease == "control" else "tumor",
        "specimen": specimen,
        "patient_key": patient_key,
        "patient_number": int(number),
    }


def source_log2_rpm_to_log1p_rpm(values: np.ndarray) -> np.ndarray:
    numeric = np.asarray(values, dtype=np.float64)
    if not np.isfinite(numeric).all() or (numeric < 0).any():
        raise ValueError("Huang matrix contains negative or non-finite values")
    return (numeric * math.log(2.0)).astype(np.float32)


def normalized_entropy(probabilities: list[float]) -> float:
    values = np.asarray(probabilities, dtype=float)
    values = np.clip(values, 0.0, None)
    total = float(values.sum())
    if total <= 0 or len(values) <= 1:
        return float("nan")
    values = values / total
    nonzero = values[values > 0]
    return float(-(nonzero * np.log(nonzero)).sum() / np.log(len(values)))


def jaccard(left: list[str], right: list[str]) -> float:
    a, b = set(map(str, left)), set(map(str, right))
    if not a and not b:
        return float("nan")
    return float(len(a & b) / len(a | b))


def cliffs_delta(left: np.ndarray, right: np.ndarray) -> float:
    x = np.asarray(left, dtype=float)
    y = np.asarray(right, dtype=float)
    return float((np.greater.outer(x, y).sum() - np.less.outer(x, y).sum()) / (len(x) * len(y)))


def bootstrap_median_difference(
    left: np.ndarray,
    right: np.ndarray,
    *,
    n_resamples: int,
    seed: int,
) -> tuple[float, float, float]:
    x = np.asarray(left, dtype=float)
    y = np.asarray(right, dtype=float)
    rng = np.random.default_rng(seed)
    differences = np.empty(n_resamples, dtype=float)
    for i in range(n_resamples):
        differences[i] = float(
            np.median(rng.choice(x, size=len(x), replace=True))
            - np.median(rng.choice(y, size=len(y), replace=True))
        )
    return (
        float(np.median(x) - np.median(y)),
        float(np.quantile(differences, 0.025)),
        float(np.quantile(differences, 0.975)),
    )


def benjamini_hochberg(p_values: list[float]) -> list[float]:
    values = np.asarray(p_values, dtype=float)
    order = np.argsort(values)
    adjusted = np.empty(len(values), dtype=float)
    running = 1.0
    for reverse_rank, index in enumerate(order[::-1], start=1):
        rank = len(values) - reverse_rank + 1
        running = min(running, float(values[index]) * len(values) / rank)
        adjusted[index] = min(running, 1.0)
    return adjusted.tolist()


def top_ids(rows: list[dict[str, Any]], key: str, k: int = 3) -> list[str]:
    return [str(row.get(key, "")) for row in rows[:k] if str(row.get(key, "")).strip()]


def candidate_flags(network_top1: str) -> dict[str, bool]:
    anterior = {
        "Frontal (agranular frontal motor areas)",
        "Lateral Prefrontal Cortex",
        "Orbitomedial Prefrontal Cortex (OMPFC)",
    }
    return {
        "anterior_fossa_candidate": network_top1 in anterior,
        "middle_fossa_candidate": network_top1 == "Temporal",
        "orbital_candidate": network_top1 == "Orbitomedial Prefrontal Cortex (OMPFC)",
        "posterior_fossa_out_of_atlas": True,
    }


def read_expression(path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    header = pd.read_csv(path, nrows=0).columns.astype(str).tolist()
    if len(header) != 160:
        raise ValueError(f"Expected 160 columns, found {len(header)}")
    sample_ids = header[1:]
    metadata = pd.DataFrame([parse_sample_id(sample_id) for sample_id in sample_ids])
    observed = metadata.groupby(["disease", "specimen"]).size().to_dict()
    if observed != EXPECTED_COUNTS:
        raise ValueError(f"Unexpected sample counts: {observed}")

    with np.load(DEFAULT_BO2023_REFERENCE_PROJECTOR, allow_pickle=False) as payload:
        projector_genes = set(payload["genes"].astype(str))
    with np.load(DEFAULT_BO2023_NETWORK_MODEL, allow_pickle=False) as payload:
        network_genes = set(payload["genes"].astype(str))
    required_genes = projector_genes | network_genes

    raw = pd.read_csv(path, dtype={header[0]: str}, low_memory=False)
    raw = raw.rename(columns={header[0]: "gene_symbol"})
    raw["gene_symbol"] = raw["gene_symbol"].astype(str).str.strip()
    if raw.shape != (83_929, 160):
        raise ValueError(f"Expected 83,929 x 160 matrix, found {raw.shape}")
    if raw["gene_symbol"].eq("").any() or raw["gene_symbol"].duplicated().any():
        raise ValueError("Gene symbols must be nonblank and unique")
    selected = raw.loc[raw["gene_symbol"].isin(required_genes)].copy()
    values = selected[sample_ids].to_numpy(dtype=np.float64)
    selected.loc[:, sample_ids] = source_log2_rpm_to_log1p_rpm(values)
    selected = selected.set_index("gene_symbol")
    return selected.astype(np.float32), metadata


def compare_tumor_control(sample_df: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    metrics = ["atlas_fit_score", "network_margin", "network_entropy"]
    counter = 0
    for specimen in ["CSF", "plasma"]:
        subset = sample_df[sample_df["specimen"].eq(specimen)]
        tumor = subset[subset["tumor_status"].eq("tumor")]
        control = subset[subset["tumor_status"].eq("control")]
        for metric in metrics:
            counter += 1
            x = tumor[metric].to_numpy(dtype=float)
            y = control[metric].to_numpy(dtype=float)
            statistic, p_value = mannwhitneyu(x, y, alternative="two-sided")
            difference, ci_low, ci_high = bootstrap_median_difference(
                x,
                y,
                n_resamples=N_BOOTSTRAP,
                seed=RANDOM_SEED + counter,
            )
            rows.append(
                {
                    "specimen": specimen,
                    "comparison": "tumor_vs_control",
                    "metric": metric,
                    "n_tumor": len(x),
                    "n_control": len(y),
                    "tumor_median": float(np.median(x)),
                    "control_median": float(np.median(y)),
                    "median_difference_tumor_minus_control": difference,
                    "bootstrap_ci95_low": ci_low,
                    "bootstrap_ci95_high": ci_high,
                    "mann_whitney_u": float(statistic),
                    "p_value": float(p_value),
                    "cliffs_delta": cliffs_delta(x, y),
                }
            )
    frame = pd.DataFrame(rows)
    frame["bh_fdr"] = benjamini_hochberg(frame["p_value"].tolist())
    return frame


def permutation_pvalue(
    left_records: list[dict[str, Any]],
    right_records: list[dict[str, Any]],
    metric: Callable[[dict[str, Any], dict[str, Any]], float],
    *,
    n_permutations: int,
    seed: int,
) -> tuple[float, float]:
    n = len(left_records)
    metric_matrix = np.asarray(
        [
            [metric(left_records[i], right_records[j]) for j in range(n)]
            for i in range(n)
        ],
        dtype=float,
    )
    rows = np.arange(n)
    observed = float(np.nanmean(metric_matrix[rows, rows]))
    rng = np.random.default_rng(seed)
    null = np.empty(n_permutations, dtype=float)
    for i in range(n_permutations):
        order = rng.permutation(n)
        null[i] = float(np.nanmean(metric_matrix[rows, order]))
    return observed, float((1 + np.sum(null >= observed)) / (n_permutations + 1))


def paired_stability(records: list[dict[str, Any]]) -> tuple[pd.DataFrame, pd.DataFrame]:
    by_sample = {row["sample_id"]: row for row in records}
    pair_rows: list[dict[str, Any]] = []
    metric_functions: dict[str, Callable[[dict[str, Any], dict[str, Any]], float]] = {
        "network_top1_agreement": lambda a, b: float(a["network_top1"] == b["network_top1"]),
        "network_top3_jaccard": lambda a, b: jaccard(a["network_top3"], b["network_top3"]),
        "network_score_spearman": lambda a, b: float(
            spearmanr(a["network_scores"], b["network_scores"]).statistic
        ),
        "region_top1_agreement": lambda a, b: float(a["region_top1"] == b["region_top1"]),
        "region_top3_jaccard": lambda a, b: jaccard(a["region_top3"], b["region_top3"]),
    }
    for disease, prefix, n_pairs in [
        ("glioma", "GLI", 16),
        ("meningioma", "MEN", 43),
        ("control", "NOR", 18),
    ]:
        for number in range(1, n_pairs + 1):
            csf = by_sample[f"{prefix}_CSF{number}"]
            plasma = by_sample[f"{prefix}_plasma{number}"]
            row = {
                "disease": disease,
                "tumor_status": "control" if disease == "control" else "tumor",
                "patient_key": f"{prefix}_{number:02d}",
                "csf_sample_id": csf["sample_id"],
                "plasma_sample_id": plasma["sample_id"],
            }
            row.update({name: function(csf, plasma) for name, function in metric_functions.items()})
            pair_rows.append(row)

    pair_frame = pd.DataFrame(pair_rows)
    summary_rows: list[dict[str, Any]] = []
    for group in ["tumor", "control"]:
        group_pairs = pair_frame[pair_frame["tumor_status"].eq(group)]
        left_records = [
            by_sample[sample_id] for sample_id in group_pairs["csf_sample_id"].astype(str)
        ]
        right_records = [
            by_sample[sample_id] for sample_id in group_pairs["plasma_sample_id"].astype(str)
        ]
        for offset, (name, function) in enumerate(metric_functions.items(), start=1):
            observed, p_value = permutation_pvalue(
                left_records,
                right_records,
                function,
                n_permutations=N_PERMUTATIONS,
                seed=RANDOM_SEED + 100 + offset + (1000 if group == "control" else 0),
            )
            values = group_pairs[name].to_numpy(dtype=float)
            summary_rows.append(
                {
                    "pair_group": group,
                    "n_pairs": len(group_pairs),
                    "metric": name,
                    "mean": float(np.nanmean(values)),
                    "median": float(np.nanmedian(values)),
                    "within_group_permutation_p_value": p_value,
                }
            )
    return pair_frame, pd.DataFrame(summary_rows)


def run(args: argparse.Namespace) -> int:
    outdir = args.outdir.resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    expression, metadata = read_expression(args.input_csv)
    sample_ids = metadata["sample_id"].astype(str).tolist()

    sample_rows: list[dict[str, Any]] = []
    network_rows: list[dict[str, Any]] = []
    region_rows: list[dict[str, Any]] = []
    records: list[dict[str, Any]] = []

    for index, meta in enumerate(metadata.to_dict("records"), start=1):
        sample_id = str(meta["sample_id"])
        frame = pd.DataFrame(
            {
                "gene_symbol": expression.index.astype(str),
                "log_tpm": expression[sample_id].to_numpy(dtype=np.float32),
            }
        )
        network_out = trace_network_expression(
            frame,
            min_overlap_fraction=0.50,
            project_to_vsd=True,
            enable_pairwise_rescue=False,
        )
        region_out = trace_bo2023_secondary_regions(
            frame,
            network_out,
            str(args.db_path),
            int(args.atlas_id),
            topk=30,
        )
        network_ranked = network_out.get("results", [])
        region_ranked = region_out.get("results", [])
        if len(network_ranked) != 10 or not region_ranked:
            raise RuntimeError(f"Untraceable formal route output for {sample_id}")

        network_top3 = top_ids(network_ranked, "network_id", 3)
        region_top3 = top_ids(region_ranked, "region_id", 3)
        network_top1 = network_top3[0]
        network_scores_by_id = {
            str(row["network_id"]): float(row["score"]) for row in network_ranked
        }
        network_score_order = sorted(network_scores_by_id)
        network_scores = [network_scores_by_id[name] for name in network_score_order]
        confidences = [float(row["confidence"]) for row in network_ranked]
        atlas_fit_score = float(network_ranked[0]["score"])
        network_margin = float(network_ranked[0]["score"] - network_ranked[1]["score"])
        flags = candidate_flags(network_top1)
        sample_row = {
            **meta,
            "route": ROUTE_NAME,
            "input_scale": "source_log2_rpm_plus1_times_ln2",
            "network_traceability": network_out["meta"].get("traceability"),
            "region_traceability": region_out["meta"].get("traceability"),
            "n_input_genes": int(expression.shape[0]),
            "n_projector_genes": network_out["meta"]["reference_projection"].get("n_projector_genes"),
            "n_projector_overlap_genes": network_out["meta"]["reference_projection"].get(
                "n_input_projector_overlap_genes"
            ),
            "n_network_model_genes": network_out["meta"].get("n_model_genes"),
            "n_network_overlap_genes": network_out["meta"].get("n_overlap_genes"),
            "network_overlap_fraction": network_out["meta"].get("overlap_fraction"),
            "atlas_fit_score": atlas_fit_score,
            "network_margin": network_margin,
            "network_entropy": normalized_entropy(confidences),
            "network_top1": network_top1,
            "network_top2": network_top3[1],
            "network_top3": " | ".join(network_top3),
            "region_top1": region_top3[0],
            "region_top2": region_top3[1],
            "region_top3": " | ".join(region_top3),
            "n_candidate_regions": region_out["meta"].get("n_candidate_regions"),
            "n_region_overlap_genes": region_out["meta"].get("n_overlap_genes"),
            **flags,
        }
        sample_rows.append(sample_row)
        for row in network_ranked:
            network_rows.append(
                {
                    "sample_id": sample_id,
                    "disease": meta["disease"],
                    "specimen": meta["specimen"],
                    **dict(row),
                }
            )
        for row in region_ranked:
            region_rows.append(
                {
                    "sample_id": sample_id,
                    "disease": meta["disease"],
                    "specimen": meta["specimen"],
                    **dict(row),
                }
            )
        records.append(
            {
                **meta,
                "network_top1": network_top1,
                "network_top3": network_top3,
                "network_score_order": network_score_order,
                "network_scores": network_scores,
                "region_top1": region_top3[0],
                "region_top3": region_top3,
            }
        )
        print(f"[{index:03d}/{len(sample_ids):03d}] {sample_id}: {network_top1} -> {region_top3[0]}")

    sample_df = pd.DataFrame(sample_rows)
    network_df = pd.DataFrame(network_rows)
    region_df = pd.DataFrame(region_rows)
    comparisons = compare_tumor_control(sample_df)
    pairs, pair_summary = paired_stability(records)
    network_distribution = (
        sample_df.groupby(["disease", "specimen", "network_top1"], dropna=False)
        .size()
        .reset_index(name="n")
    )
    coarse_candidates = (
        sample_df.groupby(["disease", "specimen"], dropna=False)[
            [
                "anterior_fossa_candidate",
                "middle_fossa_candidate",
                "orbital_candidate",
            ]
        ]
        .sum()
        .reset_index()
    )
    published_locations = pd.DataFrame(
        [
            {"disease": "glioma", "location": "anterior_fossa", "n": 6},
            {"disease": "glioma", "location": "calvarium", "n": 2},
            {"disease": "glioma", "location": "middle_fossa", "n": 8},
            {"disease": "glioma", "location": "orbital", "n": 0},
            {"disease": "glioma", "location": "posterior_fossa", "n": 1},
            {"disease": "glioma", "location": "others", "n": 1},
            {"disease": "meningioma", "location": "anterior_fossa", "n": 5},
            {"disease": "meningioma", "location": "calvarium", "n": 9},
            {"disease": "meningioma", "location": "middle_fossa", "n": 10},
            {"disease": "meningioma", "location": "orbital", "n": 3},
            {"disease": "meningioma", "location": "posterior_fossa", "n": 10},
            {"disease": "meningioma", "location": "others", "n": 9},
        ]
    )

    sample_df.to_csv(outdir / "sample_summary.csv", index=False)
    network_df.to_csv(outdir / "network_rankings.csv", index=False)
    region_df.to_csv(outdir / "exact_region_rankings.csv", index=False)
    comparisons.to_csv(outdir / "tumor_control_comparisons.csv", index=False)
    pairs.to_csv(outdir / "paired_csf_plasma_concordance.csv", index=False)
    pair_summary.to_csv(outdir / "paired_csf_plasma_summary.csv", index=False)
    network_distribution.to_csv(outdir / "network_top1_distribution.csv", index=False)
    coarse_candidates.to_csv(outdir / "coarse_candidate_flags.csv", index=False)
    published_locations.to_csv(outdir / "published_cohort_location_counts.csv", index=False)

    asset_paths = [
        DEFAULT_BO2023_NETWORK_MODEL,
        DEFAULT_BO2023_REFERENCE_PROJECTOR,
        REGION_REFERENCE,
        BEAM_PANELS,
    ]
    manifest = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "protocol_status": "frozen_before_model_outputs",
        "claim": "external_cfRNA_transfer_and_candidate_set_plausibility",
        "not_a_claim": "sample_level_anatomical_localization_validation",
        "route": ROUTE_NAME,
        "pairwise_rescue": False,
        "source_xlsb": {
            "path": str(args.source_xlsb.resolve()),
            "sha256": sha256_file(args.source_xlsb),
            "bytes": args.source_xlsb.stat().st_size,
        },
        "input_csv": {
            "path": str(args.input_csv.resolve()),
            "sha256": sha256_file(args.input_csv),
            "bytes": args.input_csv.stat().st_size,
        },
        "model_assets": {
            str(path.relative_to(ROOT)): sha256_file(path) for path in asset_paths
        },
        "n_features_public_matrix": 83_929,
        "n_samples": int(len(sample_df)),
        "n_traceable_network": int(sample_df["network_traceability"].eq("high").sum()),
        "n_traceable_region": int(sample_df["region_traceability"].eq("high").sum()),
        "n_input_model_union_genes": int(expression.shape[0]),
        "n_projector_overlap_genes": int(sample_df["n_projector_overlap_genes"].iloc[0]),
        "n_projector_genes": int(sample_df["n_projector_genes"].iloc[0]),
        "n_network_overlap_genes": int(sample_df["n_network_overlap_genes"].iloc[0]),
        "n_network_model_genes": int(sample_df["n_network_model_genes"].iloc[0]),
        "bootstrap_resamples": N_BOOTSTRAP,
        "paired_permutations": N_PERMUTATIONS,
        "random_seed": RANDOM_SEED,
        "limitations": [
            "No public patient-to-tumor-location mapping.",
            "Atlas fit is not a calibrated brain-origin probability.",
            "Huang log2(RPM+1) and Bo2023 count-derived logCPM remain cross-study domain shifted.",
            "Exact-region rankings are conditional on the locked Network Top3 beam.",
            "Coarse fossa candidate flags are descriptive and cannot produce localization accuracy.",
        ],
    }
    write_json(outdir / "manifest.json", manifest)

    primary = comparisons[comparisons["specimen"].eq("CSF")].copy()
    sensitivity = comparisons[comparisons["specimen"].eq("plasma")].copy()
    tumor_pairs = pair_summary[pair_summary["pair_group"].eq("tumor")]
    ompfc = "Orbitomedial Prefrontal Cortex (OMPFC)"
    csf_tumor = sample_df[
        sample_df["specimen"].eq("CSF") & sample_df["tumor_status"].eq("tumor")
    ]
    csf_control = sample_df[
        sample_df["specimen"].eq("CSF") & sample_df["tumor_status"].eq("control")
    ]
    plasma_tumor = sample_df[
        sample_df["specimen"].eq("plasma") & sample_df["tumor_status"].eq("tumor")
    ]
    plasma_control = sample_df[
        sample_df["specimen"].eq("plasma") & sample_df["tumor_status"].eq("control")
    ]
    min_fdr = float(comparisons["bh_fdr"].min())
    min_tumor_pair_p = float(tumor_pairs["within_group_permutation_p_value"].min())
    report = [
        "# Huang 2025 external cfRNA candidate-set analysis",
        "",
        "## Claim boundary",
        "",
        "This is an external cfRNA transfer and cohort-level candidate-set plausibility analysis. "
        "It is not sample-level anatomical localization validation.",
        "",
        "## Technical transfer",
        "",
        f"- Traceable Network outputs: {manifest['n_traceable_network']}/{manifest['n_samples']}.",
        f"- Traceable exact-region outputs: {manifest['n_traceable_region']}/{manifest['n_samples']}.",
        f"- Projector overlap: {manifest['n_projector_overlap_genes']}/{manifest['n_projector_genes']} "
        f"({manifest['n_projector_overlap_genes'] / manifest['n_projector_genes']:.1%}).",
        f"- Network-model overlap: {manifest['n_network_overlap_genes']}/{manifest['n_network_model_genes']} "
        f"({manifest['n_network_overlap_genes'] / manifest['n_network_model_genes']:.1%}).",
        "",
        "## Primary: CSF tumor versus control",
        "",
        markdown_table(primary),
        "",
        "## Sensitivity: plasma tumor versus control",
        "",
        markdown_table(sensitivity),
        "",
        "## Paired tumor CSF-plasma stability",
        "",
        markdown_table(tumor_pairs),
        "",
        "## Result assessment",
        "",
        f"- None of the six tumor-control comparisons passed BH-FDR < 0.05 "
        f"(minimum FDR {min_fdr:.3f}).",
        f"- CSF Network Top1 collapsed to OMPFC in "
        f"{int(csf_tumor['network_top1'].eq(ompfc).sum())}/{len(csf_tumor)} tumors and "
        f"{int(csf_control['network_top1'].eq(ompfc).sum())}/{len(csf_control)} controls.",
        f"- Plasma Network Top1 was OMPFC in "
        f"{int(plasma_tumor['network_top1'].eq(ompfc).sum())}/{len(plasma_tumor)} tumors and "
        f"{int(plasma_control['network_top1'].eq(ompfc).sum())}/{len(plasma_control)} controls.",
        f"- None of the five paired-tumor stability metrics exceeded the within-cohort "
        f"permutation null at p < 0.05 (minimum p {min_tumor_pair_p:.3f}).",
        "- The dataset therefore supports technical portability only. It does not provide "
        "positive evidence for tumor-associated brain-source discrimination or anatomical "
        "localization with the frozen route.",
        "",
        "## Interpretation",
        "",
        "- A significant atlas-fit difference would support transfer plausibility, not tissue specificity.",
        "- Pair concordance tests modality stability; it does not test anatomical correctness.",
        "- Published fossa counts and model candidate flags remain separate because no public patient-level join exists.",
        "",
    ]
    (outdir / "RESULTS.md").write_text("\n".join(report), encoding="utf-8")
    manifest["files"] = sorted(path.name for path in outdir.iterdir() if path.is_file())
    write_json(outdir / "manifest.json", manifest)
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description="Run the frozen Huang 2025 external cfRNA candidate-set analysis.")
    parser.add_argument("--input-csv", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--source-xlsb", type=Path, default=DEFAULT_SOURCE_XLSB)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument("--db-path", type=Path, default=DEFAULT_DB_PATH)
    parser.add_argument("--atlas-id", type=int, default=1)
    return run(parser.parse_args())


if __name__ == "__main__":
    raise SystemExit(main())
