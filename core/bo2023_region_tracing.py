from __future__ import annotations

import json
import sqlite3
from functools import lru_cache
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from core.models import softmax_confidence, trace_corr


DEFAULT_TOP50_WEIGHT = 0.25
DEFAULT_VALIDATION_SUMMARY = (
    Path(__file__).resolve().parents[1]
    / "data"
    / "models"
    / "bo2023_exact_region_validation_summary.json"
)
ROUTE_NAME = "top3_beam_local_top50_top100_zfusion_w0p25"


def _table_exists(conn: sqlite3.Connection, table: str) -> bool:
    return conn.execute(
        "SELECT 1 FROM sqlite_master WHERE type='table' AND name=?",
        (table,),
    ).fetchone() is not None


@lru_cache(maxsize=4)
def _load_reference_matrix(db_path: str, atlas_id: int) -> tuple[pd.DataFrame, dict[str, str]]:
    conn = sqlite3.connect(db_path)
    try:
        if not _table_exists(conn, "reference_expression") or not _table_exists(conn, "macaque_brain_atlas"):
            return pd.DataFrame(), {}
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
        return pd.DataFrame(), {}
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
    return matrix, region_network


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


def _validation_metrics(summary_path: Path = DEFAULT_VALIDATION_SUMMARY) -> dict[str, Any]:
    if not summary_path.exists():
        return {}
    try:
        summary = json.loads(summary_path.read_text(encoding="utf-8"))
    except Exception:
        return {}
    route = summary.get("routes", {}).get(ROUTE_NAME, {})
    baseline = summary.get("routes", {}).get(summary.get("baseline_route", ""), {})
    return {
        "route": ROUTE_NAME,
        "n_test_samples": summary.get("n_test_samples"),
        "top1_accuracy": route.get("top1_accuracy"),
        "top3_accuracy": route.get("top3_accuracy"),
        "baseline_top1_accuracy": baseline.get("top1_accuracy"),
        "baseline_top3_accuracy": baseline.get("top3_accuracy"),
        "validation_design": summary.get("validation_design"),
    }


def trace_bo2023_secondary_regions(
    expression: pd.DataFrame,
    network_output: dict[str, Any],
    db_path: str,
    atlas_id: int,
    topk: int = 30,
    top50_weight: float = DEFAULT_TOP50_WEIGHT,
) -> dict[str, Any]:
    """Score Bo2023 exact Region candidates inside the validated SaleemNetworks Top3 beam."""
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

    matrix, region_network = _load_reference_matrix(db_path, atlas_id)
    if matrix.empty or not region_network:
        return {
            "results": [],
            "meta": {
                "endpoint": "Bo2023 exact Region",
                "method": ROUTE_NAME,
                "traceability": "insufficient",
                "error": "missing Bo2023 reference matrix or region-network mapping",
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

    sample = expression[["gene_symbol", "tpm_value"]].dropna().copy()
    sample["gene_symbol"] = sample["gene_symbol"].astype(str)
    sample["tpm_value"] = pd.to_numeric(sample["tpm_value"], errors="coerce").fillna(0.0)
    series = sample.groupby("gene_symbol")["tpm_value"].mean()
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
    gene_order = _candidate_gene_order(candidate_matrix, max_genes=100)
    rows50 = gene_order[: min(50, len(gene_order))]
    rows100 = gene_order[: min(100, len(gene_order))]
    scores50 = trace_corr(candidate_matrix[rows50, :], vector[rows50])
    scores100 = trace_corr(candidate_matrix[rows100, :], vector[rows100])
    fused = float(top50_weight) * _zscore(scores50) + (1.0 - float(top50_weight)) * _zscore(scores100)
    confidence = softmax_confidence(fused)
    order = np.argsort(fused)[::-1]
    rows = []
    for rank, idx in enumerate(order[: max(1, int(topk))], start=1):
        region = candidate_regions[int(idx)]
        rows.append(
            {
                "region_id": region,
                "rank": int(rank),
                "score": float(fused[int(idx)]),
                "confidence": float(confidence[int(idx)]),
                "network_id": region_network.get(region),
                "top50_corr_component": float(scores50[int(idx)]),
                "top100_corr_component": float(scores100[int(idx)]),
            }
        )
    return {
        "results": rows,
        "meta": {
            "endpoint": "Bo2023 exact Region",
            "method": ROUTE_NAME,
            "network_beam": top_networks,
            "n_candidate_regions": int(len(candidate_regions)),
            "candidate_region_source": "SaleemNetworks Top3 beam",
            "n_reference_genes": int(matrix.shape[0]),
            "n_overlap_genes": int(len(overlap_genes)),
            "n_scoring_genes_top50": int(len(rows50)),
            "n_scoring_genes_top100": int(len(rows100)),
            "top50_weight": float(top50_weight),
            "traceability": "high",
            "result_interpretation": (
                "Secondary exact Region candidates scored within the validated SaleemNetworks Top3 beam "
                "using Top50/Top100 discriminative-gene correlation z-fusion."
            ),
            "full_loso_validation": _validation_metrics(),
        },
    }
