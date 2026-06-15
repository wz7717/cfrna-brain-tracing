#!/usr/bin/env python
from __future__ import annotations

import argparse
import gzip
import json
import re
import sys
from collections import Counter
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from core.network_tracing import trace_network_expression  # noqa: E402
from core.tumor_adapted_network_tracing import (  # noqa: E402
    trace_tumor_adapted_network_expression,
)
from scripts.analyze_gse228512_ev_rna_tracing import (  # noqa: E402
    ADAPTED_MODEL,
    BASELINE_MODEL,
    NETWORK_TO_PRIMARY_BROAD,
    add_route,
    model_scores,
)


DATA_DIR = ROOT / "data" / "external_validation" / "GSE189919"
OUTDIR = ROOT / "results" / "gse189919_csf_tracing_validation_20260613"


def sample_key(value: str) -> str:
    return re.sub(r"[^a-z0-9]", "", str(value).lower())


def read_soft_metadata(path: Path) -> pd.DataFrame:
    raw = gzip.open(path, "rt", encoding="utf-8", errors="replace").read()
    rows = []
    for block in raw.split("^SAMPLE = ")[1:]:
        accession = re.search(r"!Sample_geo_accession = (.*)", block)
        title = re.search(r"!Sample_title = (.*)", block)
        source = re.search(r"!Sample_source_name_ch1 = (.*)", block)
        characteristics: dict[str, str] = {}
        for item in re.findall(r"!Sample_characteristics_ch1 = (.*)", block):
            if ":" in item:
                key, value = item.split(":", 1)
                characteristics[key.strip().lower()] = value.strip()
        disease = characteristics.get("disease", "")
        rows.append(
            {
                "geo_accession": accession.group(1) if accession else "",
                "geo_title": title.group(1) if title else "",
                "sample_key": sample_key(title.group(1) if title else ""),
                "group": "MB" if disease == "Medulloblastoma" else "Normal",
                "molecular_subgroup": characteristics.get("sub group", "not available"),
                "biofluid": source.group(1) if source else "",
            }
        )
    return pd.DataFrame(rows)


def read_expression(
    data_dir: Path,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    counts = pd.read_csv(data_dir / "GSE189919_count.csv.gz", index_col="Geneid")
    tpm = pd.read_csv(data_dir / "GSE189919_tpm_count.csv.gz", index_col="Geneid")
    for frame in (counts, tpm):
        frame.index = frame.index.astype(str).str.strip()
        frame.columns = frame.columns.astype(str)
    counts = counts.apply(pd.to_numeric, errors="raise")
    tpm = tpm.apply(pd.to_numeric, errors="raise")

    if not counts.index.equals(tpm.index):
        raise ValueError("Count and TPM gene orders differ")
    if not counts.columns.equals(tpm.columns):
        raise ValueError("Count and TPM sample orders differ")
    if counts.index.duplicated().any() or tpm.index.duplicated().any():
        raise ValueError("Duplicated gene symbols in official matrix")
    if counts.isna().any().any() or tpm.isna().any().any():
        raise ValueError("Missing expression value")
    if (counts < 0).any().any() or (tpm < 0).any().any():
        raise ValueError("Negative expression value")

    soft = read_soft_metadata(data_dir / "GSE189919_family.soft.gz")
    by_key = soft.set_index("sample_key")
    matrix_keys = pd.Index([sample_key(column) for column in counts.columns])
    missing = matrix_keys.difference(by_key.index)
    extra = by_key.index.difference(matrix_keys)
    if len(missing) or len(extra):
        raise ValueError(
            f"Matrix/SOFT mismatch: missing metadata={missing.tolist()}, "
            f"extra metadata={extra.tolist()}"
        )
    metadata = by_key.loc[matrix_keys].reset_index(drop=True)
    metadata.insert(0, "sample_id", counts.columns)

    library_sizes = counts.sum(axis=0)
    if (library_sizes <= 0).any():
        raise ValueError("Zero-sized count library")
    cpm = counts.divide(library_sizes, axis=1) * 1_000_000.0
    return counts.astype(float), tpm.astype(float), cpm.astype(float), metadata


def production_predictions(
    expression: pd.DataFrame,
    output: pd.DataFrame,
    route: str,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    transformed = np.log1p(expression)
    top1, top3, switched = [], [], []
    score_by_network = []
    for sample in transformed.columns:
        frame = pd.DataFrame(
            {
                "gene_symbol": transformed.index,
                "tpm_value": transformed[sample].to_numpy(dtype=float),
            }
        )
        traced = trace_network_expression(frame, min_overlap_fraction=0.50)
        ranked = traced.get("results", [])
        if len(ranked) != 10:
            raise ValueError(f"{route}/{sample}: expected 10 Network results")
        top1.append(str(ranked[0]["network_id"]))
        top3.append(" | ".join(str(row["network_id"]) for row in ranked[:3]))
        switched.append(
            bool(traced.get("meta", {}).get("pairwise_rescue", {}).get("switched", False))
        )
        score_by_network.append(
            {str(row["network_id"]): float(row["score"]) for row in ranked}
        )
    output[f"{route}_top1_network"] = top1
    output[f"{route}_top3_network"] = top3
    output[f"{route}_top1_broad"] = output[f"{route}_top1_network"].map(
        NETWORK_TO_PRIMARY_BROAD
    )
    output[f"{route}_top3_broad"] = output[f"{route}_top3_network"].map(
        lambda value: " | ".join(
            dict.fromkeys(
                NETWORK_TO_PRIMARY_BROAD[item.strip()]
                for item in str(value).split("|")
                if item.strip()
            )
        )
    )
    output[f"{route}_pairwise_switched"] = switched
    return output, {
        "score_by_network": score_by_network,
        "n_switched": int(sum(switched)),
        "fraction_switched": float(np.mean(switched)),
    }


def adapted_production_predictions(
    expression: pd.DataFrame,
    output: pd.DataFrame,
    route: str,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    top1, top3, accepted = [], [], []
    maximum_correlation, score_by_network = [], []
    for sample in expression.columns:
        frame = pd.DataFrame(
            {
                "gene_symbol": expression.index,
                "tpm_value": expression[sample].to_numpy(dtype=float),
            }
        )
        traced = trace_tumor_adapted_network_expression(frame)
        ranked = traced.get("results", [])
        if len(ranked) != 10:
            raise ValueError(f"{route}/{sample}: expected 10 Network results")
        top1.append(str(ranked[0]["network_id"]))
        top3.append(" | ".join(str(row["network_id"]) for row in ranked[:3]))
        accepted.append(bool(traced["meta"]["ood_accepted"]))
        maximum_correlation.append(float(traced["meta"]["maximum_correlation"]))
        score_by_network.append(
            {str(row["network_id"]): float(row["score"]) for row in ranked}
        )
    output[f"{route}_top1_network"] = top1
    output[f"{route}_top3_network"] = top3
    output[f"{route}_top1_broad"] = output[f"{route}_top1_network"].map(
        NETWORK_TO_PRIMARY_BROAD
    )
    output[f"{route}_top3_broad"] = output[f"{route}_top3_network"].map(
        lambda value: " | ".join(
            dict.fromkeys(
                NETWORK_TO_PRIMARY_BROAD[item.strip()]
                for item in str(value).split("|")
                if item.strip()
            )
        )
    )
    output[f"{route}_maximum_raw_correlation"] = maximum_correlation
    output[f"{route}_ood_accepted"] = accepted
    return output, {"score_by_network": score_by_network}


def sample_detection(
    expression: pd.DataFrame,
    genes: np.ndarray,
    metadata: pd.DataFrame,
    model_name: str,
) -> pd.DataFrame:
    values = expression.reindex(genes).fillna(0.0)
    rows = []
    for column, row in zip(expression.columns, metadata.itertuples(index=False)):
        detected = int((values[column] > 0).sum())
        rows.append(
            {
                "sample_id": column,
                "group": row.group,
                "molecular_subgroup": row.molecular_subgroup,
                "model": model_name,
                "detected_markers": detected,
                "total_markers": int(len(genes)),
                "detected_fraction": float(detected / len(genes)),
            }
        )
    return pd.DataFrame(rows)


def prediction_distribution(
    output: pd.DataFrame, routes: list[str]
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    top1_rows, top3_rows, collapse_rows = [], [], []
    for route in routes:
        for group, sub in output.groupby("group"):
            counts = sub[f"{route}_top1_network"].value_counts()
            proportions = counts / len(sub)
            hhi = float(np.square(proportions).sum())
            collapse_rows.append(
                {
                    "route": route,
                    "group": group,
                    "n": int(len(sub)),
                    "active_networks": int(len(counts)),
                    "dominant_network": str(counts.index[0]),
                    "dominant_fraction": float(proportions.iloc[0]),
                    "hhi": hhi,
                    "effective_network_count": float(1.0 / hhi),
                }
            )
            for network, count in counts.items():
                top1_rows.append(
                    {
                        "route": route,
                        "group": group,
                        "network": network,
                        "broad_anatomy": NETWORK_TO_PRIMARY_BROAD[network],
                        "count": int(count),
                        "fraction": float(count / len(sub)),
                    }
                )
            occurrences = Counter(
                item.strip()
                for value in sub[f"{route}_top3_network"]
                for item in str(value).split("|")
                if item.strip()
            )
            for network, count in occurrences.items():
                top3_rows.append(
                    {
                        "route": route,
                        "group": group,
                        "network": network,
                        "broad_anatomy": NETWORK_TO_PRIMARY_BROAD[network],
                        "samples_with_label_in_top3": int(count),
                        "fraction": float(count / len(sub)),
                    }
                )
    return (
        pd.DataFrame(top1_rows),
        pd.DataFrame(top3_rows),
        pd.DataFrame(collapse_rows),
    )


def route_agreement(output: pd.DataFrame, routes: list[str]) -> pd.DataFrame:
    rows = []
    for left_index, left in enumerate(routes):
        for right in routes[left_index + 1 :]:
            for group, sub in output.groupby("group"):
                top1 = (
                    sub[f"{left}_top1_network"].to_numpy()
                    == sub[f"{right}_top1_network"].to_numpy()
                )
                left_top3 = sub[f"{left}_top3_network"].map(
                    lambda value: set(item.strip() for item in str(value).split("|"))
                )
                right_top3 = sub[f"{right}_top3_network"].map(
                    lambda value: set(item.strip() for item in str(value).split("|"))
                )
                overlap = [
                    len(a & b) / len(a | b) if a | b else 1.0
                    for a, b in zip(left_top3, right_top3)
                ]
                rows.append(
                    {
                        "left_route": left,
                        "right_route": right,
                        "group": group,
                        "n": int(len(sub)),
                        "top1_agreement": float(top1.mean()),
                        "mean_top3_jaccard": float(np.mean(overlap)),
                    }
                )
    return pd.DataFrame(rows)


def audit_algorithms(
    output: pd.DataFrame,
    expression_by_source: dict[str, pd.DataFrame],
    direct_routes: dict[str, tuple[str, Path, str]],
    production_audits: dict[str, dict[str, Any]],
    adapted_production_audits: dict[str, dict[str, Any]],
    adapted_threshold: float,
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []

    def record(check: str, passed: bool, detail: str) -> None:
        rows.append({"check": check, "passed": bool(passed), "detail": detail})

    for route, (source, model_path, transform) in direct_routes.items():
        scores, networks, _, _ = model_scores(
            expression_by_source[source], model_path, transform
        )
        probabilities = np.exp(scores - scores.max(axis=0, keepdims=True))
        probabilities /= probabilities.sum(axis=0, keepdims=True)
        stored_scores = output[
            [f"{route}_score_{index + 1:02d}" for index in range(len(networks))]
        ].to_numpy(dtype=float).T
        record(
            f"{route}: stored scores equal recomputed scores",
            np.allclose(scores, stored_scores, atol=1e-12),
            f"max_abs_error={np.max(np.abs(scores - stored_scores)):.3e}",
        )
        manual = np.empty_like(scores)
        model = np.load(model_path, allow_pickle=False)
        genes = model["genes"].astype(str)
        values = expression_by_source[source].reindex(genes).fillna(0.0).to_numpy(float)
        if transform == "logcpm":
            values = np.log1p(values)
        elif transform == "harmonized":
            target = np.sort(model["target_quantiles"].astype(float))
            mapped = np.empty_like(values)
            for column in range(values.shape[1]):
                order = np.argsort(values[:, column], kind="mergesort")
                mapped[order, column] = target
            values = mapped
        reference = model["reference"].astype(float)
        for i in range(reference.shape[1]):
            for j in range(values.shape[1]):
                ref = reference[:, i]
                val = values[:, j]
                if np.std(ref) <= 1e-12 or np.std(val) <= 1e-12:
                    manual[i, j] = 0.0
                else:
                    manual[i, j] = np.corrcoef(ref, val)[0, 1]
        record(
            f"{route}: vectorized Pearson equals independent np.corrcoef",
            np.allclose(scores, manual, atol=1e-10),
            f"max_abs_error={np.max(np.abs(scores - manual)):.3e}",
        )
        record(
            f"{route}: softmax probabilities sum to one",
            np.allclose(probabilities.sum(axis=0), 1.0, atol=1e-12),
            f"max_sum_error={np.max(np.abs(probabilities.sum(axis=0) - 1)):.3e}",
        )
        expected_top1 = networks[np.argmax(scores, axis=0)]
        record(
            f"{route}: Top1 equals score argmax",
            np.array_equal(expected_top1, output[f"{route}_top1_network"].to_numpy()),
            f"n={len(output)}",
        )
        record(
            f"{route}: margins are nonnegative",
            bool((output[f"{route}_raw_margin"] >= -1e-12).all()),
            f"min_margin={output[f'{route}_raw_margin'].min():.6g}",
        )
        entropy = output[f"{route}_normalized_entropy"]
        record(
            f"{route}: normalized entropy lies in [0,1]",
            bool(((entropy >= -1e-12) & (entropy <= 1 + 1e-12)).all()),
            f"range={entropy.min():.6g}..{entropy.max():.6g}",
        )
        if route.startswith("adapted_"):
            expected_ood = output[f"{route}_max_score"] >= adapted_threshold
            record(
                f"{route}: OOD decision equals threshold comparison",
                expected_ood.equals(output[f"{route}_ood_accepted"]),
                f"threshold={adapted_threshold:.6g}",
            )

    for source_name, production_route, direct_route in (
        ("author_tpm", "baseline_production_tpm", "baseline_tpm_log"),
        ("count_cpm", "baseline_production_count_cpm", "baseline_count_logcpm"),
    ):
        model = np.load(BASELINE_MODEL, allow_pickle=False)
        networks = model["networks"].astype(str)
        direct_scores = output[
            [f"{direct_route}_score_{index + 1:02d}" for index in range(len(networks))]
        ].to_numpy(dtype=float)
        production_scores = np.asarray(
            [
                [sample_scores[network] for network in networks]
                for sample_scores in production_audits[production_route]["score_by_network"]
            ]
        )
        record(
            f"{production_route}: production scores equal direct {source_name} scores",
            np.allclose(direct_scores, production_scores, atol=1e-12),
            f"max_abs_error={np.max(np.abs(direct_scores - production_scores)):.3e}",
        )
        nonswitched = ~output[f"{production_route}_pairwise_switched"].to_numpy(bool)
        record(
            f"{production_route}: non-rescued Top1 equals raw argmax",
            np.array_equal(
                output.loc[nonswitched, f"{production_route}_top1_network"].to_numpy(),
                output.loc[nonswitched, f"{direct_route}_top1_network"].to_numpy(),
            ),
            f"non_rescued_n={int(nonswitched.sum())}",
        )

    adapted_model = np.load(ADAPTED_MODEL, allow_pickle=False)
    adapted_genes = adapted_model["genes"].astype(str)
    adapted_networks = adapted_model["networks"].astype(str)
    target = adapted_model["target_quantiles"].astype(float)
    offsets = adapted_model["calibration_offsets"].astype(float)
    reference = adapted_model["reference"].astype(float)
    for source_name, route in (
        ("author_tpm", "adapted_production_tpm"),
        ("count_cpm", "adapted_production_count_cpm"),
    ):
        values = expression_by_source[source_name].reindex(adapted_genes).fillna(0.0)
        expected_scores = []
        expected_raw_max = []
        for sample in values.columns:
            vector = np.log1p(np.clip(values[sample].to_numpy(dtype=float), 0, None))
            order = np.argsort(vector, kind="mergesort")
            harmonized = np.empty_like(vector)
            harmonized[order] = target
            raw_scores = np.asarray(
                [
                    0.0
                    if np.std(reference[:, index]) <= 1e-12
                    or np.std(harmonized) <= 1e-12
                    else np.corrcoef(reference[:, index], harmonized)[0, 1]
                    for index in range(reference.shape[1])
                ]
            )
            expected_scores.append(raw_scores + offsets)
            expected_raw_max.append(float(raw_scores.max()))
        expected_scores = np.asarray(expected_scores)
        stored_scores = np.asarray(
            [
                [sample_scores[network] for network in adapted_networks]
                for sample_scores in adapted_production_audits[route][
                    "score_by_network"
                ]
            ]
        )
        record(
            f"{route}: production calibrated scores equal independent calculation",
            np.allclose(expected_scores, stored_scores, atol=1e-10),
            f"max_abs_error={np.max(np.abs(expected_scores - stored_scores)):.3e}",
        )
        record(
            f"{route}: OOD uses maximum raw correlation before class offsets",
            np.array_equal(
                np.asarray(expected_raw_max) >= adapted_threshold,
                output[f"{route}_ood_accepted"].to_numpy(bool),
            ),
            f"threshold={adapted_threshold:.6g}",
        )
        record(
            f"{route}: Top1 equals calibrated-score argmax",
            np.array_equal(
                adapted_networks[np.argmax(expected_scores, axis=1)],
                output[f"{route}_top1_network"].to_numpy(),
            ),
            f"n={len(output)}",
        )

    prediction_columns = [
        column for column in output.columns if column.endswith("_top3_network")
    ]
    top3_valid = all(
        output[column].map(
            lambda value: len([x for x in str(value).split("|") if x.strip()]) == 3
            and len(set(x.strip() for x in str(value).split("|"))) == 3
        ).all()
        for column in prediction_columns
    )
    record("All Top3 outputs contain three unique Networks", top3_valid, f"routes={len(prediction_columns)}")
    broad_columns = [
        column for column in output.columns if column.endswith("_top1_broad")
    ]
    record(
        "All Top1 Broad anatomy labels are mapped",
        bool(output[broad_columns].notna().all().all()),
        f"routes={len(broad_columns)}",
    )
    numeric = output.select_dtypes(include=[np.number])
    record(
        "Prediction table contains no non-finite numeric values",
        bool(np.isfinite(numeric.to_numpy()).all()),
        f"numeric_cells={numeric.size}",
    )
    return pd.DataFrame(rows)


def methodological_review(
    expression_by_source: dict[str, pd.DataFrame],
) -> pd.DataFrame:
    model = np.load(ADAPTED_MODEL, allow_pickle=False)
    genes = model["genes"].astype(str)
    target = model["target_quantiles"].astype(float)
    rows = []
    for source, expression in expression_by_source.items():
        values = expression.reindex(genes).fillna(0.0)
        affected = 0
        maximum_tie_spread = 0.0
        for sample in values.columns:
            vector = np.log1p(np.clip(values[sample].to_numpy(dtype=float), 0, None))
            order = np.argsort(vector, kind="mergesort")
            mapped = np.empty_like(vector)
            mapped[order] = target
            for value in np.unique(vector):
                indices = np.flatnonzero(vector == value)
                if len(indices) > 1:
                    spread = float(np.ptp(mapped[indices]))
                    maximum_tie_spread = max(maximum_tie_spread, spread)
                    if spread > 1e-12:
                        affected += 1
                        break
        rows.append(
            {
                "component": "adapted production quantile mapping",
                "input_source": source,
                "status": "warning",
                "affected_samples": int(affected),
                "n_samples": int(values.shape[1]),
                "maximum_mapped_spread_within_equal_input_tie": maximum_tie_spread,
                "finding": (
                    "Equal input values receive different target quantiles according "
                    "to stable gene order. This is deterministic and matches the trained "
                    "model, but can create order-dependent signal in sparse liquid-biopsy data."
                ),
                "required_action": (
                    "Develop a tie-aware mapping, then retrain calibration offsets and "
                    "the OOD threshold before replacing the current production route."
                ),
            }
        )
    return pd.DataFrame(rows)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", type=Path, default=DATA_DIR)
    parser.add_argument("--outdir", type=Path, default=OUTDIR)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    counts, tpm, cpm, metadata = read_expression(args.data_dir)
    output = metadata.copy()
    expression_by_source = {"author_tpm": tpm, "count_cpm": cpm}

    direct_specs: dict[str, tuple[str, Path, str]] = {
        "baseline_tpm_log": ("author_tpm", BASELINE_MODEL, "logcpm"),
        "baseline_count_logcpm": ("count_cpm", BASELINE_MODEL, "logcpm"),
        "adapted_tpm": ("author_tpm", ADAPTED_MODEL, "cpm"),
        "adapted_logtpm": ("author_tpm", ADAPTED_MODEL, "logcpm"),
        "adapted_harmonized_tpm": ("author_tpm", ADAPTED_MODEL, "harmonized"),
        "adapted_count_cpm": ("count_cpm", ADAPTED_MODEL, "cpm"),
        "adapted_count_logcpm": ("count_cpm", ADAPTED_MODEL, "logcpm"),
        "adapted_harmonized_count": ("count_cpm", ADAPTED_MODEL, "harmonized"),
    }
    adapted = np.load(ADAPTED_MODEL, allow_pickle=False)
    adapted_threshold = float(adapted["ood_max_correlation_threshold"][0])
    model_coverage: dict[str, Any] = {}
    for route, (source, model_path, transform) in direct_specs.items():
        scores, networks, genes, present = model_scores(
            expression_by_source[source], model_path, transform
        )
        threshold = adapted_threshold if route.startswith("adapted_") else None
        output = add_route(output, route, scores, networks, threshold)
        model_coverage[route] = {
            "present": int(present.sum()),
            "total": int(len(genes)),
            "fraction": float(present.mean()),
        }

    production_audits: dict[str, dict[str, Any]] = {}
    for source, route in (
        ("author_tpm", "baseline_production_tpm"),
        ("count_cpm", "baseline_production_count_cpm"),
    ):
        output, production_audits[route] = production_predictions(
            expression_by_source[source], output, route
        )

    adapted_production_audits: dict[str, dict[str, Any]] = {}
    for source, route in (
        ("author_tpm", "adapted_production_tpm"),
        ("count_cpm", "adapted_production_count_cpm"),
    ):
        output, adapted_production_audits[route] = adapted_production_predictions(
            expression_by_source[source], output, route
        )

    prediction_routes = [
        "baseline_production_tpm",
        "baseline_production_count_cpm",
        "adapted_production_tpm",
        "adapted_production_count_cpm",
        *direct_specs.keys(),
    ]
    top1, top3, collapse = prediction_distribution(output, prediction_routes)
    agreement = route_agreement(output, prediction_routes)

    baseline = np.load(BASELINE_MODEL, allow_pickle=False)
    detection = pd.concat(
        [
            sample_detection(tpm, baseline["genes"].astype(str), metadata, "baseline_tpm"),
            sample_detection(tpm, adapted["genes"].astype(str), metadata, "adapted_tpm"),
            sample_detection(cpm, baseline["genes"].astype(str), metadata, "baseline_count"),
            sample_detection(cpm, adapted["genes"].astype(str), metadata, "adapted_count"),
        ],
        ignore_index=True,
    )
    audit = audit_algorithms(
        output,
        expression_by_source,
        direct_specs,
        production_audits,
        adapted_production_audits,
        adapted_threshold,
    )
    method_review = methodological_review(expression_by_source)
    if not audit["passed"].all():
        failed = audit.loc[~audit["passed"], "check"].tolist()
        raise AssertionError(f"Algorithm audit failed: {failed}")

    sample_qc = metadata.copy()
    sample_qc["count_library_size"] = counts.sum(axis=0).to_numpy(dtype=float)
    sample_qc["count_detected_genes"] = (counts > 0).sum(axis=0).to_numpy(dtype=int)
    sample_qc["tpm_detected_genes"] = (tpm > 0).sum(axis=0).to_numpy(dtype=int)
    sample_qc["tpm_sum"] = tpm.sum(axis=0).to_numpy(dtype=float)

    output.to_csv(
        args.outdir / "gse189919_sample_predictions.csv",
        index=False,
        encoding="utf-8-sig",
    )
    sample_qc.to_csv(args.outdir / "sample_input_qc.csv", index=False, encoding="utf-8-sig")
    detection.to_csv(
        args.outdir / "sample_marker_detection.csv", index=False, encoding="utf-8-sig"
    )
    top1.to_csv(
        args.outdir / "network_broad_top1_distribution.csv",
        index=False,
        encoding="utf-8-sig",
    )
    top3.to_csv(
        args.outdir / "network_broad_top3_occurrence.csv",
        index=False,
        encoding="utf-8-sig",
    )
    collapse.to_csv(
        args.outdir / "prediction_collapse_metrics.csv",
        index=False,
        encoding="utf-8-sig",
    )
    agreement.to_csv(
        args.outdir / "route_agreement.csv", index=False, encoding="utf-8-sig"
    )
    audit.to_csv(
        args.outdir / "algorithm_step_audit.csv", index=False, encoding="utf-8-sig"
    )
    method_review.to_csv(
        args.outdir / "methodological_review.csv", index=False, encoding="utf-8-sig"
    )

    summary = {
        "dataset": "GSE189919",
        "validation_scope": (
            "External technical validation of Network/Broad anatomy tracing. "
            "No cerebellum assumption and no localization accuracy calculation."
        ),
        "n_samples": int(len(output)),
        "group_counts": output["group"].value_counts().to_dict(),
        "molecular_subgroup_counts": output["molecular_subgroup"].value_counts().to_dict(),
        "n_gene_rows": int(len(tpm)),
        "input_qc": {
            "count_tpm_gene_order_identical": bool(counts.index.equals(tpm.index)),
            "count_tpm_sample_order_identical": bool(counts.columns.equals(tpm.columns)),
            "median_count_library_size": float(counts.sum(axis=0).median()),
            "median_detected_genes": float((counts > 0).sum(axis=0).median()),
            "median_tpm_sum": float(tpm.sum(axis=0).median()),
        },
        "model_gene_coverage": model_coverage,
        "ood_threshold": adapted_threshold,
        "ood_acceptance": {
            route: {
                group: {
                    "accepted": int(
                        output.loc[output.group.eq(group), f"{route}_ood_accepted"].sum()
                    ),
                    "n": int(output.group.eq(group).sum()),
                    "fraction": float(
                        output.loc[output.group.eq(group), f"{route}_ood_accepted"].mean()
                    ),
                }
                for group in ("MB", "Normal")
            }
            for route in direct_specs
            if route.startswith("adapted_")
        },
        "pairwise_rescue": {
            route: {
                "n_switched": item["n_switched"],
                "fraction_switched": item["fraction_switched"],
            }
            for route, item in production_audits.items()
        },
        "algorithm_audit": {
            "checks": int(len(audit)),
            "passed": int(audit["passed"].sum()),
            "failed": int((~audit["passed"]).sum()),
        },
        "methodological_review": {
            "warnings": int(method_review.status.eq("warning").sum()),
            "components": method_review.component.unique().tolist(),
        },
    }
    (args.outdir / "analysis_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
