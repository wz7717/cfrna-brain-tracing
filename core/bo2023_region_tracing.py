from __future__ import annotations

import json
import sqlite3
from functools import lru_cache
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from core.models import softmax_confidence, trace_corr
from core.reference_projection import (
    compute_logcpm,
    map_index_to_symbols,
    read_bo2023_gene_matrix,
    read_gene_map,
)
from core.region_resolution import load_region_resolution_model


DEFAULT_TOP50_WEIGHT = 0.25
DEFAULT_LOCAL_TOP_N_GENES = 200
ROOT = Path(__file__).resolve().parents[1]
DEFAULT_BO2023_COUNTS = ROOT / "bo2023 data" / "mfas5_819samples_28415genes_featurecounts_counts.txt"
DEFAULT_BO2023_SAMPLE_INFO = ROOT / "bo2023 data" / "Information of sequenced samples_update_full878_filter819.xlsx"
DEFAULT_BO2023_GENE_MAP = ROOT / "bo2023_bulk_atlas_buildkit" / "04_expressed_genes_neocortex_plus_subcortical.cleaned_symbols.csv"
DEFAULT_VALIDATION_SUMMARY = (
    ROOT
    / "data"
    / "models"
    / "bo2023_exact_region_validation_summary.json"
)
DEFAULT_FORMAL_THREE_TIER_SUMMARY = (
    ROOT
    / "results"
    / "bo2023_reference_projection_20260616_cleaned_symbols"
    / "formal_three_tier_loso_hybrid"
    / "hybrid_formal_loso_summary.json"
)
ROUTE_NAME = "projected_vsd_network_top3_logcpm_resolution_local_exact"
LEGACY_EXACT_ROUTE_NAME = "top3_beam_local_top50_top100_zfusion_w0p25"


def _table_exists(conn: sqlite3.Connection, table: str) -> bool:
    return conn.execute(
        "SELECT 1 FROM sqlite_master WHERE type='table' AND name=?",
        (table,),
    ).fetchone() is not None


def _read_bo2023_sample_metadata(path: Path) -> pd.DataFrame:
    info = pd.read_excel(path, sheet_name="mfas5_819samples_phenSet4", usecols=["No.", "Region", "SaleemNetworks"])
    info["sample_id"] = info["No."].astype(str).str.strip()
    info["region_id"] = info["Region"].astype(str).str.strip()
    info["network_id"] = info["SaleemNetworks"].astype(str).str.strip()
    info = info.drop_duplicates("sample_id").set_index("sample_id")
    return info[info["region_id"].ne("") & info["network_id"].ne("")].copy()


@lru_cache(maxsize=2)
def _load_raw_logcpm_reference_matrix(
    counts_path: Path = DEFAULT_BO2023_COUNTS,
    sample_info_path: Path = DEFAULT_BO2023_SAMPLE_INFO,
    gene_map_path: Path = DEFAULT_BO2023_GENE_MAP,
) -> tuple[pd.DataFrame, dict[str, str], str]:
    if not counts_path.exists() or not sample_info_path.exists() or not gene_map_path.exists():
        return pd.DataFrame(), {}, "raw_logcpm_files_missing"

    gene_map = read_gene_map(gene_map_path)
    counts, _ = map_index_to_symbols(read_bo2023_gene_matrix(counts_path, dtype="float32"), gene_map)
    metadata = _read_bo2023_sample_metadata(sample_info_path)
    samples = [sample for sample in counts.columns.astype(str) if sample in metadata.index]
    if not samples:
        return pd.DataFrame(), {}, "raw_logcpm_no_metadata_sample_overlap"

    logcpm = compute_logcpm(counts.loc[:, samples])
    region_series = metadata.loc[samples, "region_id"].astype(str)
    region_network: dict[str, str] = {}
    for region, rows in metadata.loc[samples].groupby("region_id"):
        networks = sorted(rows["network_id"].astype(str).dropna().unique().tolist())
        if len(networks) == 1:
            region_network[str(region)] = networks[0]

    region_matrix = logcpm.T.assign(region_id=region_series.to_numpy()).groupby("region_id").mean().T
    region_matrix.index = region_matrix.index.astype(str)
    region_matrix.columns = region_matrix.columns.astype(str)
    region_matrix = region_matrix.loc[region_matrix.abs().sum(axis=1) > 0].sort_index()
    return region_matrix.astype("float32"), region_network, "raw_featurecounts_logcpm"


@lru_cache(maxsize=4)
def _load_db_reference_matrix(db_path: str, atlas_id: int) -> tuple[pd.DataFrame, dict[str, str], str]:
    conn = sqlite3.connect(db_path)
    try:
        if not _table_exists(conn, "reference_expression") or not _table_exists(conn, "macaque_brain_atlas"):
            return pd.DataFrame(), {}, "db_reference_missing_tables"
        ref = pd.read_sql_query(
            """
            SELECT gene_symbol, region_id, avg_tpm
            FROM reference_expression
            WHERE atlas_id = ?
            """,
            conn,
            params=(int(atlas_id),),
        )
        atlas = pd.read_sql_query(
            """
            SELECT region_id, coordinates
            FROM macaque_brain_atlas
            WHERE atlas_id = ?
            """,
            conn,
            params=(int(atlas_id),),
        )
    finally:
        conn.close()

    if ref.empty or atlas.empty:
        return pd.DataFrame(), {}, "db_reference_empty"
    ref = ref.dropna(subset=["gene_symbol", "region_id", "avg_tpm"]).copy()
    ref["gene_symbol"] = ref["gene_symbol"].astype(str)
    ref["region_id"] = ref["region_id"].astype(str)
    ref["avg_tpm"] = pd.to_numeric(ref["avg_tpm"], errors="coerce").fillna(0.0)
    matrix = (
        ref.pivot_table(index="gene_symbol", columns="region_id", values="avg_tpm", aggfunc="mean")
        .fillna(0.0)
        .sort_index()
    )
    matrix = matrix.loc[matrix.abs().sum(axis=1) > 0]

    region_network: dict[str, str] = {}
    for row in atlas.itertuples(index=False):
        coordinates = {}
        try:
            coordinates = json.loads(row.coordinates) if row.coordinates else {}
        except Exception:
            coordinates = {}
        network = str(coordinates.get("saleem_network", "") or "").strip()
        if network:
            region_network[str(row.region_id)] = network
    return matrix, region_network, "db_reference_expression_avg_tpm_fallback"


def _load_reference_matrix(db_path: str, atlas_id: int) -> tuple[pd.DataFrame, dict[str, str], str]:
    matrix, region_network, source = _load_raw_logcpm_reference_matrix()
    if not matrix.empty and region_network:
        return matrix, region_network, source
    return _load_db_reference_matrix(db_path, atlas_id)


def _candidate_gene_order(reference: np.ndarray, max_genes: int) -> np.ndarray:
    if reference.size == 0:
        return np.array([], dtype=int)
    row_std = reference.std(axis=1)
    row_range = reference.max(axis=1) - reference.min(axis=1)
    score = row_std + 0.1 * row_range
    return np.argsort(score)[::-1][: min(int(max_genes), reference.shape[0])]


def _zscore(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    std = float(values.std())
    if std <= 1e-12:
        return np.zeros_like(values, dtype=float)
    return (values - float(values.mean())) / std


def _validation_metrics(
    summary_path: Path = DEFAULT_VALIDATION_SUMMARY,
    formal_summary_path: Path = DEFAULT_FORMAL_THREE_TIER_SUMMARY,
) -> dict[str, Any]:
    formal: dict[str, Any] = {}
    if formal_summary_path.exists():
        try:
            formal = json.loads(formal_summary_path.read_text(encoding="utf-8"))
        except Exception:
            formal = {}
    if not summary_path.exists():
        return {"route": ROUTE_NAME, "formal_three_tier_loso": formal} if formal else {}
    try:
        summary = json.loads(summary_path.read_text(encoding="utf-8"))
    except Exception:
        return {"route": ROUTE_NAME, "formal_three_tier_loso": formal} if formal else {}
    route = summary.get("routes", {}).get(LEGACY_EXACT_ROUTE_NAME, {})
    baseline = summary.get("routes", {}).get(summary.get("baseline_route", ""), {})
    return {
        "route": ROUTE_NAME,
        "legacy_exact_route": LEGACY_EXACT_ROUTE_NAME,
        "n_test_samples": summary.get("n_test_samples"),
        "legacy_exact_top1_accuracy": route.get("top1_accuracy"),
        "legacy_exact_top3_accuracy": route.get("top3_accuracy"),
        "legacy_baseline_top1_accuracy": baseline.get("top1_accuracy"),
        "legacy_baseline_top3_accuracy": baseline.get("top3_accuracy"),
        "legacy_validation_design": summary.get("validation_design"),
        "formal_three_tier_loso": formal,
    }


def _sample_logcpm_series(expression: pd.DataFrame) -> tuple[pd.Series, str]:
    sample = expression.dropna(subset=["gene_symbol"]).copy()
    sample["gene_symbol"] = sample["gene_symbol"].astype(str).str.strip()
    if "read_count" in sample.columns:
        read_count = pd.to_numeric(sample["read_count"], errors="coerce").fillna(0.0).clip(lower=0.0)
        if float(read_count.sum()) > 0:
            cpm = read_count / float(read_count.sum()) * 1_000_000.0
            sample["_score_value"] = np.log1p(cpm)
            return sample.groupby("gene_symbol")["_score_value"].mean(), "read_count_logcpm"
    if "log_tpm" in sample.columns:
        sample["_score_value"] = pd.to_numeric(sample["log_tpm"], errors="coerce")
        if sample["_score_value"].notna().any():
            return sample.groupby("gene_symbol")["_score_value"].mean().fillna(0.0), "stored_log_tpm_fallback"
    sample["tpm_value"] = pd.to_numeric(sample.get("tpm_value", 0.0), errors="coerce").fillna(0.0)
    sample["_score_value"] = np.log1p(sample["tpm_value"].clip(lower=0.0))
    return sample.groupby("gene_symbol")["_score_value"].mean(), "log1p_tpm_fallback"


def _resolution_entry(
    model: dict[str, Any],
    region: str,
    region_network: dict[str, str],
) -> dict[str, Any]:
    network = str(region_network.get(region, ""))
    entry = model.get("entries", {}).get(f"{network}||{region}") if model else None
    if not entry:
        return {
            "resolution_tier": "low_resolution",
            "resolution_group": region,
            "group_members": [region],
            "resolution_reasons": ["outside_region_resolution_model"],
            "group_plausibility_tier": "not_calibrated",
            "group_calibration_flags": [],
        }
    return entry


def _rank_resolution_groups(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    groups: dict[str, dict[str, Any]] = {}
    for row in rows:
        group = str(row.get("resolution_group", row["region_id"]))
        item = groups.setdefault(
            group,
            {
                "resolution_group": group,
                "best_region_id": row["region_id"],
                "group_members": row.get("resolution_group_members", row["region_id"]),
                "member_region_ids": [],
                "member_scores": [],
                "resolution_tier": row.get("resolution_tier", "low_resolution"),
                "manual_review_recommended": bool(row.get("manual_review_recommended", True)),
                "group_plausibility_tier": row.get("group_plausibility_tier", "not_calibrated"),
            },
        )
        item["member_region_ids"].append(row["region_id"])
        item["member_scores"].append(float(row["exact_local_score"]))
        item["manual_review_recommended"] = bool(item["manual_review_recommended"]) or bool(row.get("manual_review_recommended", True))
        if float(row["exact_local_score"]) > max(item["member_scores"][:-1] or [float("-inf")]):
            item["best_region_id"] = row["region_id"]
            item["resolution_tier"] = row.get("resolution_tier", "low_resolution")
            item["group_plausibility_tier"] = row.get("group_plausibility_tier", "not_calibrated")

    group_rows = []
    for item in groups.values():
        scores = item.pop("member_scores")
        group_score = float(max(scores) + 0.10 * np.mean(scores))
        item["group_score"] = group_score
        item["mean_member_score"] = float(np.mean(scores))
        item["n_returned_group_members"] = int(len(scores))
        group_rows.append(item)
    group_rows.sort(key=lambda row: float(row["group_score"]), reverse=True)
    for rank, row in enumerate(group_rows, start=1):
        row["rank"] = int(rank)
    return group_rows


def _apply_group_first_order(rows: list[dict[str, Any]], group_rows: list[dict[str, Any]], topk: int) -> list[dict[str, Any]]:
    group_rank = {row["resolution_group"]: int(row["rank"]) for row in group_rows}
    group_score = {row["resolution_group"]: float(row["group_score"]) for row in group_rows}
    ordered = sorted(
        rows,
        key=lambda row: (group_rank.get(str(row.get("resolution_group")), 9999), -float(row["exact_local_score"])),
    )
    final_scores = np.asarray([float(row["exact_local_score"]) for row in ordered], dtype=float)
    confidence = softmax_confidence(final_scores) if len(final_scores) else np.array([], dtype=float)
    for rank, row in enumerate(ordered[: max(1, int(topk))], start=1):
        group = str(row.get("resolution_group"))
        row["rank"] = int(rank)
        row["score"] = float(row["exact_local_score"])
        row["confidence"] = float(confidence[rank - 1]) if rank - 1 < len(confidence) else 0.0
        row["resolution_group_rank"] = int(group_rank.get(group, 9999))
        row["resolution_group_score"] = float(group_score.get(group, float("nan")))
    return ordered[: max(1, int(topk))]


def trace_bo2023_secondary_regions(
    expression: pd.DataFrame,
    network_output: dict[str, Any],
    db_path: str,
    atlas_id: int,
    topk: int = 30,
    top50_weight: float = DEFAULT_TOP50_WEIGHT,
    local_top_n_genes: int = DEFAULT_LOCAL_TOP_N_GENES,
) -> dict[str, Any]:
    """Formal three-tier Bo2023 route.

    Stage 1 is supplied by ``network_output``: projected VSD SaleemNetworks Top3 beam.
    Stages 2-3 use logCPM-compatible local evidence: resolution group rerank, then
    exact-region rerank inside the ordered groups.
    """
    network_rows = network_output.get("results", [])
    top_networks = [str(row.get("network_id", "")) for row in network_rows[:3] if str(row.get("network_id", "")).strip()]
    if not top_networks:
        return {
            "results": [],
            "meta": {
                "endpoint": "Bo2023 exact Region",
                "method": ROUTE_NAME,
                "traceability": "insufficient",
                "error": "missing Network Top3 beam",
            },
        }

    matrix, region_network, reference_source = _load_reference_matrix(db_path, atlas_id)
    if matrix.empty or not region_network:
        return {
            "results": [],
            "meta": {
                "endpoint": "Bo2023 exact Region",
                "method": ROUTE_NAME,
                "traceability": "insufficient",
                "error": "missing Bo2023 logCPM reference matrix or region-network mapping",
            },
        }

    candidate_regions = [
        region for region in matrix.columns.astype(str).tolist()
        if region_network.get(region) in set(top_networks)
    ]
    if len(candidate_regions) < 2:
        return {
            "results": [],
            "meta": {
                "endpoint": "Bo2023 exact Region",
                "method": ROUTE_NAME,
                "traceability": "insufficient",
                "error": "fewer than two candidate regions in Network Top3 beam",
                "network_beam": top_networks,
            },
        }

    series, query_scale = _sample_logcpm_series(expression)
    overlap_genes = matrix.index.intersection(series.index)
    if len(overlap_genes) < 20:
        return {
            "results": [],
            "meta": {
                "endpoint": "Bo2023 exact Region",
                "method": ROUTE_NAME,
                "traceability": "insufficient",
                "error": "insufficient Bo2023 Region gene overlap",
                "n_overlap_genes": int(len(overlap_genes)),
                "network_beam": top_networks,
            },
        }

    candidate_matrix = matrix.loc[overlap_genes, candidate_regions].to_numpy(dtype=float)
    vector = series.reindex(overlap_genes).fillna(0.0).to_numpy(dtype=float)
    gene_order = _candidate_gene_order(candidate_matrix, max_genes=local_top_n_genes)
    rows50 = gene_order[: min(50, len(gene_order))]
    rows100 = gene_order[: min(100, len(gene_order))]
    scores50 = trace_corr(candidate_matrix[rows50, :], vector[rows50])
    scores100 = trace_corr(candidate_matrix[rows100, :], vector[rows100])
    fused = float(top50_weight) * _zscore(scores50) + (1.0 - float(top50_weight)) * _zscore(scores100)
    model = load_region_resolution_model()
    rows = []
    for idx in np.argsort(fused)[::-1].tolist():
        region = candidate_regions[int(idx)]
        entry = _resolution_entry(model, region, region_network)
        group_members = [str(x) for x in entry.get("group_members", [region])]
        rows.append(
            {
                "region_id": region,
                "network_id": region_network.get(region),
                "exact_local_score": float(fused[int(idx)]),
                "top50_corr_component": float(scores50[int(idx)]),
                "top100_corr_component": float(scores100[int(idx)]),
                "resolution_tier": str(entry.get("resolution_tier", "low_resolution")),
                "resolution_group": str(entry.get("resolution_group", region)),
                "resolution_group_members": " | ".join(group_members),
                "resolution_reasons": ";".join(map(str, entry.get("resolution_reasons", []))),
                "group_plausibility_tier": str(entry.get("group_plausibility_tier", "not_calibrated")),
                "group_calibration_flags": ";".join(map(str, entry.get("group_calibration_flags", []))),
                "manual_review_recommended": str(entry.get("resolution_tier", "low_resolution")) == "low_resolution",
            }
        )
    group_rows = _rank_resolution_groups(rows)
    rows = _apply_group_first_order(rows, group_rows, topk)
    top = rows[0] if rows else {}
    return {
        "results": rows,
        "meta": {
            "endpoint": "Bo2023 formal three-tier Region",
            "method": ROUTE_NAME,
            "route_stages": [
                "projected_vsd_network_top3_beam",
                "logcpm_resolution_group_rerank",
                "logcpm_local_exact_rerank",
            ],
            "network_beam": top_networks,
            "n_candidate_regions": int(len(candidate_regions)),
            "candidate_region_source": "SaleemNetworks Top3 beam",
            "reference_expression_source": reference_source,
            "query_expression_source": query_scale,
            "n_reference_genes": int(matrix.shape[0]),
            "n_overlap_genes": int(len(overlap_genes)),
            "n_local_candidate_genes": int(len(gene_order)),
            "n_scoring_genes_top50": int(len(rows50)),
            "n_scoring_genes_top100": int(len(rows100)),
            "top50_weight": float(top50_weight),
            "traceability": "high",
            "result_interpretation": (
                "Formal three-tier source tracing: projected VSD is used only for the validated Network Top3 "
                "candidate beam; logCPM-compatible local evidence reranks resolution groups and exact regions."
            ),
            "region_resolution_annotation": {
                "enabled": True,
                "primary_network": str(top_networks[0]) if top_networks else None,
                "top1_resolution_tier": top.get("resolution_tier"),
                "top1_resolution_group": top.get("resolution_group"),
                "top1_group_members": top.get("resolution_group_members"),
                "top1_group_plausibility_tier": top.get("group_plausibility_tier"),
                "manual_review_recommended": bool(top.get("manual_review_recommended", False)),
                "interpretation": "Resolution group ranking is an active stage before exact-region reranking.",
                "group_ranking_method": "best_exact_local_score_plus_0p10_mean_returned_member_score",
                "group_ranking": group_rows[:10],
            },
            "full_loso_validation": _validation_metrics(),
        },
    }
