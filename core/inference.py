from __future__ import annotations

from pathlib import Path

import pandas as pd

import core.bo2023_region_tracing as region_tracing
import core.network_tracing as network_tracing
from core.model_lock import ModelLockError
from core.production_route import verify_production_route


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_DB_PATH = str((ROOT / "braintrace_source_tracing.db").resolve())


def run_locked_three_tier_route(
    expression: pd.DataFrame,
    *,
    atlas_id: int | None = None,
    db_path: str = DEFAULT_DB_PATH,
    topk: int = 30,
) -> tuple[dict, dict]:
    """Run the single locked production route used by both Streamlit and CLI."""
    parameters, model_lock = verify_production_route()
    network_out = network_tracing.trace_network_expression(
        expression,
        min_overlap_fraction=float(parameters["network_min_overlap_fraction"]),
        project_to_vsd=bool(parameters["project_to_vsd"]),
        enable_pairwise_rescue=bool(parameters["enable_pairwise_rescue"]),
    )
    network_meta = network_out.get("meta", {})
    if (
        int(network_meta.get("n_networks", -1)) != int(parameters["network_count"])
        or int(network_meta.get("n_model_genes", -1)) != int(parameters["network_gene_count"])
    ):
        raise ModelLockError("production Network runtime metadata differs from the frozen route")
    if not network_out.get("results"):
        raise ValueError(
            f"Insufficient Network model-gene overlap: {network_meta.get('n_overlap_genes', 0)}/"
            f"{network_meta.get('n_model_genes', 0)}."
        )
    lock_meta = {"lock_id": model_lock["lock_id"], "status": model_lock["status"]}
    network_out.setdefault("meta", {})["model_lock"] = lock_meta

    out = region_tracing.trace_bo2023_secondary_regions(
        expression,
        network_out,
        db_path,
        atlas_id,
        topk=max(int(topk), int(parameters["network_top_k"])),
        network_top_k=int(parameters["network_top_k"]),
        top50_weight=float(parameters["exact_top50_weight"]),
        exact_top50_gene_count=int(parameters["exact_top50_gene_count"]),
        exact_top100_gene_count=int(parameters["exact_top100_gene_count"]),
        local_top_n_genes=int(parameters["region_local_top_n_genes"]),
        min_region_gene_overlap=int(parameters["region_min_overlap_genes"]),
    )
    if not out.get("results"):
        meta = out.get("meta", {})
        raise ValueError(str(meta.get("error") or "Bo2023 three-tier route returned no region candidates."))

    region_meta = out.get("meta", {})
    expected_region_runtime = {
        "network_top_k": int(parameters["network_top_k"]),
        "n_local_candidate_genes": int(parameters["region_local_top_n_genes"]),
        "n_scoring_genes_top50": int(parameters["exact_top50_gene_count"]),
        "n_scoring_genes_top100": int(parameters["exact_top100_gene_count"]),
        "min_required_region_overlap_genes": int(parameters["region_min_overlap_genes"]),
    }
    changed_runtime = sorted(
        key
        for key, expected in expected_region_runtime.items()
        if int(region_meta.get(key, -1)) != expected
    )
    if changed_runtime:
        raise ModelLockError(
            "production Region runtime metadata differs from the frozen route: "
            + ", ".join(changed_runtime)
        )
    out.setdefault("meta", {})["model_lock"] = lock_meta
    return network_out, out
