from __future__ import annotations

import json
from functools import lru_cache
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from core.models import softmax_confidence, trace_corr


DEFAULT_MODEL = (
    Path(__file__).resolve().parents[1]
    / "data"
    / "models"
    / "bo2023_tcga_tumor_adapted_network_model.npz"
)


@lru_cache(maxsize=2)
def load_tumor_adapted_model(path: Path = DEFAULT_MODEL) -> dict[str, Any]:
    with np.load(path, allow_pickle=False) as data:
        model = {
            "genes": data["genes"].astype(str),
            "networks": data["networks"].astype(str),
            "reference": data["reference"].astype(float),
            "target_quantiles": data["target_quantiles"].astype(float),
            "calibration_offsets": data["calibration_offsets"].astype(float),
            "ood_threshold": float(data["ood_max_correlation_threshold"][0]),
        }
    metadata_path = path.with_suffix(".json")
    model["metadata"] = (
        json.loads(metadata_path.read_text(encoding="utf-8")) if metadata_path.exists() else {}
    )
    return model


def quantile_map(values: np.ndarray, target_quantiles: np.ndarray) -> np.ndarray:
    order = np.argsort(values, kind="mergesort")
    mapped = np.empty_like(values, dtype=float)
    mapped[order] = target_quantiles
    return mapped


def trace_tumor_adapted_network_expression(
    expression: pd.DataFrame,
    model_path: Path = DEFAULT_MODEL,
    min_overlap_fraction: float = 0.80,
) -> dict[str, Any]:
    model = load_tumor_adapted_model(model_path)
    sample = expression[["gene_symbol", "tpm_value"]].dropna().copy()
    sample["gene_symbol"] = sample["gene_symbol"].astype(str)
    sample["tpm_value"] = pd.to_numeric(sample["tpm_value"], errors="coerce").fillna(0.0)
    series = sample.groupby("gene_symbol")["tpm_value"].mean()
    present = np.asarray([gene in series.index for gene in model["genes"]], dtype=bool)
    overlap_fraction = float(present.mean()) if len(present) else 0.0
    if overlap_fraction < min_overlap_fraction:
        return {
            "results": [],
            "meta": {
                "method": "tumor_adapted_quantile_mapped_correlation",
                "traceability": "insufficient",
                "n_model_genes": int(len(model["genes"])),
                "n_overlap_genes": int(present.sum()),
                "overlap_fraction": overlap_fraction,
                "error": "insufficient tumor-adapted model gene overlap",
            },
        }

    vector = series.reindex(model["genes"]).fillna(0.0).to_numpy(dtype=float)
    harmonized = quantile_map(np.log1p(np.clip(vector, 0, None)), model["target_quantiles"])
    raw_scores = trace_corr(model["reference"], harmonized)
    scores = raw_scores + model["calibration_offsets"]
    order = np.argsort(scores)[::-1]
    confidence = softmax_confidence(scores)
    maximum_correlation = float(np.max(raw_scores))
    accepted = maximum_correlation >= model["ood_threshold"]
    rows = [
        {
            "network_id": str(model["networks"][int(index)]),
            "rank": rank,
            "score": float(scores[int(index)]),
            "confidence": float(confidence[int(index)]),
            "exploratory_only": not accepted,
        }
        for rank, index in enumerate(order, start=1)
    ]
    return {
        "results": rows,
        "meta": {
            "method": "tumor_adapted_quantile_mapped_correlation_with_ood_rejection",
            "traceability": "exploratory" if not accepted else "accepted",
            "n_model_genes": int(len(model["genes"])),
            "n_overlap_genes": int(present.sum()),
            "overlap_fraction": overlap_fraction,
            "ood_accepted": accepted,
            "maximum_correlation": maximum_correlation,
            "ood_threshold": model["ood_threshold"],
            "class_prior_offsets": model["calibration_offsets"].tolist(),
            "model_metadata": model["metadata"],
        },
    }
