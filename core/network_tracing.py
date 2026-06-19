from __future__ import annotations

import json
from functools import lru_cache
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from core.models import softmax_confidence, trace_corr
from core.reference_projection import apply_projector, load_projector_npz


DEFAULT_BO2023_NETWORK_MODEL = (
    Path(__file__).resolve().parents[1]
    / "data"
    / "models"
    / "bo2023_saleem_network_top200_model.npz"
)
DEFAULT_BO2023_NETWORK_METADATA = DEFAULT_BO2023_NETWORK_MODEL.with_suffix(".json")
DEFAULT_BO2023_NETWORK_PAIRWISE_RESCUE = (
    Path(__file__).resolve().parents[1]
    / "data"
    / "models"
    / "bo2023_saleem_network_pairwise_rescue_model.json"
)
DEFAULT_BO2023_REFERENCE_PROJECTOR = (
    Path(__file__).resolve().parents[1]
    / "data"
    / "models"
    / "bo2023_reference_projector_linear_full.npz"
)


@lru_cache(maxsize=4)
def load_network_model(model_path: Path = DEFAULT_BO2023_NETWORK_MODEL) -> dict[str, Any]:
    with np.load(model_path, allow_pickle=False) as data:
        return {
            "genes": data["genes"].astype(str),
            "networks": data["networks"].astype(str),
            "reference": data["reference"].astype(float),
            "fisher_scores": data["fisher_scores"].astype(float),
        }


def _pair_key(left: str, right: str) -> str:
    return "||".join(sorted((str(left), str(right))))


@lru_cache(maxsize=4)
def _load_pairwise_rescue_model(model_path: Path = DEFAULT_BO2023_NETWORK_PAIRWISE_RESCUE) -> dict[str, Any]:
    if not model_path.exists():
        return {}
    try:
        model = json.loads(model_path.read_text(encoding="utf-8"))
    except Exception:
        return {}
    pairs = {}
    for item in model.get("pairs", []):
        left = str(item.get("left_network", ""))
        right = str(item.get("right_network", ""))
        genes = [str(gene) for gene in item.get("genes", [])]
        reference = item.get("reference", {})
        if not left or not right or not genes or left not in reference or right not in reference:
            continue
        pairs[_pair_key(left, right)] = {
            "left": left,
            "right": right,
            "genes": genes,
            "reference": np.column_stack(
                [
                    np.asarray(reference[left], dtype=float),
                    np.asarray(reference[right], dtype=float),
                ]
            ),
            "training_error_count": int(item.get("training_error_count", 0)),
        }
    model["pair_lookup"] = pairs
    return model


def _apply_pairwise_rescue(
    scores: np.ndarray,
    networks: np.ndarray,
    series: pd.Series,
    rescue_model: dict[str, Any],
) -> tuple[list[int], dict[str, Any]]:
    order = np.argsort(scores)[::-1].tolist()
    if not rescue_model:
        return order, {"enabled": False, "switched": False}
    pairs = rescue_model.get("pair_lookup", {})
    params = rescue_model.get("parameters", {})
    candidate_k = int(params.get("candidate_k", 3))
    min_margin = float(params.get("pair_min_margin", 0.0))
    anchor = order[0]
    anchor_network = str(networks[anchor])
    best_position = -1
    best_margin = min_margin
    best_info: dict[str, Any] = {}
    for position in range(1, min(candidate_k, len(order))):
        challenger = order[position]
        challenger_network = str(networks[challenger])
        pair = pairs.get(_pair_key(anchor_network, challenger_network))
        if not pair:
            continue
        pair_genes = pair["genes"]
        present = [gene for gene in pair_genes if gene in series.index]
        if len(present) < max(20, int(0.5 * len(pair_genes))):
            continue
        gene_pos = [pair_genes.index(gene) for gene in present]
        reference = pair["reference"][gene_pos, :]
        vector = series.reindex(present).fillna(0.0).to_numpy(dtype=float)
        pair_scores = trace_corr(reference, vector)
        left = pair["left"]
        right = pair["right"]
        pair_index = {left: 0, right: 1}
        margin = float(pair_scores[pair_index[challenger_network]] - pair_scores[pair_index[anchor_network]])
        if margin > best_margin:
            best_position = position
            best_margin = margin
            best_info = {
                "pair": f"{left} <> {right}",
                "challenger": challenger_network,
                "anchor": anchor_network,
                "margin": margin,
                "n_pair_genes": int(len(pair_genes)),
                "n_overlap_pair_genes": int(len(present)),
                "training_error_count": int(pair.get("training_error_count", 0)),
            }
    if best_position >= 1:
        order[0], order[best_position] = order[best_position], order[0]
        return order, {"enabled": True, "switched": True, **best_info}
    return order, {"enabled": True, "switched": False}


@lru_cache(maxsize=2)
def _load_reference_projector(projector_path: Path = DEFAULT_BO2023_REFERENCE_PROJECTOR):
    if not projector_path.exists():
        return None
    return load_projector_npz(projector_path)


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


def trace_network_expression(
    expression: pd.DataFrame,
    model_path: Path = DEFAULT_BO2023_NETWORK_MODEL,
    metadata_path: Path = DEFAULT_BO2023_NETWORK_METADATA,
    pairwise_rescue_path: Path = DEFAULT_BO2023_NETWORK_PAIRWISE_RESCUE,
    projector_path: Path = DEFAULT_BO2023_REFERENCE_PROJECTOR,
    min_overlap_fraction: float = 0.50,
    project_to_vsd: bool = True,
) -> dict[str, Any]:
    model = load_network_model(model_path)
    input_series, input_scale = _sample_logcpm_series(expression)
    projector = _load_reference_projector(projector_path) if project_to_vsd else None
    projection_meta: dict[str, Any] = {
        "enabled": bool(projector is not None and project_to_vsd),
        "input_scale": input_scale,
        "projector_path": str(projector_path) if projector_path else None,
    }
    if projector is not None and project_to_vsd:
        projected = apply_projector(projector, pd.DataFrame({"query": input_series}))
        series = projected["query"]
        projection_meta.update(
            {
                "output_scale": "projected_vsd",
                "n_projector_genes": int(len(projector.genes)),
                "n_input_projector_overlap_genes": int(pd.Index(projector.genes).intersection(input_series.index).size),
            }
        )
    else:
        series = input_series
        projection_meta["output_scale"] = input_scale
    present = np.asarray([gene in input_series.index for gene in model["genes"]], dtype=bool)
    overlap_fraction = float(present.mean()) if len(present) else 0.0
    if overlap_fraction < float(min_overlap_fraction):
        return {
            "results": [],
            "meta": {
                "endpoint": "SaleemNetworks",
                "method": "discriminative_gene_pearson_correlation",
                "n_networks": int(len(model["networks"])),
                "n_model_genes": int(len(model["genes"])),
                "n_overlap_genes": int(present.sum()),
                "overlap_fraction": overlap_fraction,
                "traceability": "insufficient",
                "error": "insufficient network-model gene overlap",
                "reference_projection": projection_meta,
            },
        }
    vector = series.reindex(model["genes"]).fillna(0.0).to_numpy(dtype=float)
    scores = trace_corr(model["reference"], vector)
    confidence = softmax_confidence(scores)
    pairwise_model = _load_pairwise_rescue_model(pairwise_rescue_path)
    order, rescue_meta = _apply_pairwise_rescue(scores, model["networks"], series, pairwise_model)
    rows = [
        {
            "network_id": str(model["networks"][int(j)]),
            "rank": int(rank),
            "score": float(scores[int(j)]),
            "confidence": float(confidence[int(j)]),
        }
        for rank, j in enumerate(order, start=1)
    ]
    metadata: dict[str, Any] = {}
    if metadata_path.exists():
        metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    return {
        "results": rows,
        "meta": {
            "endpoint": "SaleemNetworks",
            "method": (
                "discriminative_gene_pearson_correlation"
                if not rescue_meta.get("enabled")
                else "discriminative_gene_pearson_correlation_with_top3_pairwise_rescue"
            ),
            "n_networks": int(len(model["networks"])),
            "n_model_genes": int(len(model["genes"])),
            "n_overlap_genes": int(present.sum()),
            "overlap_fraction": overlap_fraction,
            "traceability": "high",
            "pairwise_rescue": rescue_meta,
            "reference_projection": projection_meta,
            "model_metadata": metadata,
            "pairwise_rescue_validation": pairwise_model.get("validation", {}) if pairwise_model else {},
        },
    }
