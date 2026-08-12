from __future__ import annotations

import itertools
import json
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "reproducibility" / "round5_analysis"
ABLATION = ROOT / "reproducibility" / "round5_projection_ablation" / "bo2023_projected_vsd_loso_detail.csv"
PARAMS = ROOT / "reproducibility" / "p0_bio3_projector" / "projector_gene_parameters.csv"
MODEL_GENES = ROOT / "data" / "models" / "bo2023_saleem_network_top200_model_genes.csv"
META = ROOT / "bo2023 data" / "Information of sequenced samples_update_full878_filter819.xlsx"
CANON = ROOT / "reports" / "validation_recheck_20260713_canonical110"
LOMO_LOCKED = ROOT / "reports" / "model_validation" / "model_validation_stage2_20260716" / "lomo_frozen_no_pairwise"


def exact_sign_flip(values: np.ndarray) -> float:
    values = np.asarray(values, dtype=float)
    observed = abs(float(values.mean()))
    stats = []
    for signs in itertools.product((-1.0, 1.0), repeat=len(values)):
        stats.append(abs(float(np.mean(values * np.asarray(signs)))))
    return float(np.mean(np.asarray(stats) >= observed - 1e-15))


def projection_ablation() -> tuple[pd.DataFrame, pd.DataFrame]:
    detail = pd.read_csv(ABLATION)
    metadata = pd.read_excel(META, sheet_name="mfas5_819samples_phenSet4", usecols=["No.", "MonkeyID"])
    metadata["sample_id"] = metadata["No."].astype(str).str.strip()
    monkey = metadata.drop_duplicates("sample_id").set_index("sample_id")["MonkeyID"]
    detail["monkey_id"] = detail["sample_id"].map(monkey)
    detail["reciprocal_rank"] = 1.0 / detail["true_rank"].astype(float)
    detail["ndcg_at_3"] = np.where(
        detail["true_rank"].le(3), 1.0 / np.log2(detail["true_rank"].astype(float) + 1.0), 0.0
    )
    summary = detail.groupby("route", as_index=False).agg(
        n=("sample_id", "size"),
        top1=("hit1", "mean"),
        top3=("hit3", "mean"),
        mrr=("reciprocal_rank", "mean"),
        ndcg_at_3=("ndcg_at_3", "mean"),
        median_true_rank=("true_rank", "median"),
        mean_decision_margin=("decision_margin", "mean"),
    )
    donor = detail.groupby(["route", "monkey_id"], as_index=False).agg(
        n=("sample_id", "size"), top1=("hit1", "mean"), top3=("hit3", "mean"),
        mrr=("reciprocal_rank", "mean"), ndcg_at_3=("ndcg_at_3", "mean")
    )
    rows = []
    wide = donor.pivot(index="monkey_id", columns="route", values=["top1", "top3", "mrr", "ndcg_at_3"])
    for metric in ["top1", "top3", "mrr", "ndcg_at_3"]:
        for comparator in ["logcpm_baseline", "native_vsd"]:
            diffs = (wide[(metric, "projected_vsd")] - wide[(metric, comparator)]).to_numpy()
            rows.append({
                "metric": metric,
                "contrast": f"projected_vsd - {comparator}",
                "n_donors": len(diffs),
                "mean_donor_difference": float(diffs.mean()),
                "exact_sign_flip_p_two_sided": exact_sign_flip(diffs),
            })
    return summary, donor, pd.DataFrame(rows)


def ols_quality() -> tuple[pd.DataFrame, dict]:
    params = pd.read_csv(PARAMS)
    model = pd.read_csv(MODEL_GENES)
    gene_col = "gene_symbol" if "gene_symbol" in model.columns else model.columns[0]
    locked = set(model[gene_col].astype(str))
    params["panel"] = np.where(params["gene_symbol"].astype(str).isin(locked), "locked_network_200", "other")
    metrics = ["r2", "spearman_r", "residual_sd", "slope", "intercept"]
    rows = []
    for panel_name, frame in [("all_21668", params), ("locked_network_200", params[params["panel"].eq("locked_network_200")])]:
        for metric in metrics:
            x = pd.to_numeric(frame[metric], errors="coerce").dropna()
            rows.append({
                "panel": panel_name, "metric": metric, "n": len(x), "mean": float(x.mean()),
                "sd": float(x.std(ddof=1)), "p10": float(x.quantile(.10)), "q1": float(x.quantile(.25)),
                "median": float(x.median()), "q3": float(x.quantile(.75)), "p90": float(x.quantile(.90)),
            })
    fallback = params["fallback_reason"].fillna("").astype(str).ne("")
    summary = {
        "n_all_genes": int(len(params)),
        "n_locked_network_genes": int(params["panel"].eq("locked_network_200").sum()),
        "n_fallback_genes": int(fallback.sum()),
        "fraction_r2_below_0p10_all": float(params["r2"].lt(.10).mean()),
        "fraction_r2_below_0p10_locked": float(params.loc[params["panel"].eq("locked_network_200"), "r2"].lt(.10).mean()),
        "interpretation": "In-sample paired-reference engineering fit; not a proof of DESeq2 VST equivalence or external calibration.",
    }
    return pd.DataFrame(rows), summary


def ranking_metrics() -> pd.DataFrame:
    paths = {
        ("LOSO", "Network"): CANON / "loso" / "hybrid_formal_loso_network_detail.csv",
        ("LOSO", "ResolutionGroup"): CANON / "loso" / "hybrid_formal_loso_resolution_group_detail.csv",
        ("LOSO", "ExactRegion"): CANON / "loso" / "hybrid_formal_loso_exact_region_detail.csv",
        ("LOMO", "Network"): LOMO_LOCKED / "formal_lomo_network_detail.csv",
        ("LOMO", "ResolutionGroup"): LOMO_LOCKED / "formal_lomo_resolution_group_detail.csv",
        ("LOMO", "ExactRegion"): LOMO_LOCKED / "formal_lomo_exact_region_detail.csv",
    }
    rows = []
    for (scheme, level), path in paths.items():
        frame = pd.read_csv(path)
        if "route_family" in frame.columns and frame["sample_id"].duplicated().any():
            frame = frame.loc[frame["route_family"].eq("hybrid_projected_network_logcpm_exact")].copy()
        if frame["sample_id"].duplicated().any():
            raise ValueError(f"Duplicate sample IDs remain after route selection: {path}")
        rank_col = "group_true_rank" if level == "ResolutionGroup" else "true_rank"
        candidate_col = "n_candidate_groups" if level == "ResolutionGroup" else ("n_candidate_regions" if level == "ExactRegion" else None)
        rank = pd.to_numeric(frame[rank_col], errors="coerce")
        retrieved = rank.notna()
        if candidate_col:
            candidates = pd.to_numeric(frame[candidate_col], errors="coerce")
            retrieved &= rank.le(candidates)
        rr = np.where(retrieved, 1.0 / rank, 0.0)
        ndcg3 = np.where(retrieved & rank.le(3), 1.0 / np.log2(rank + 1.0), 0.0)
        rows.append({
            "scheme": scheme, "level": level, "n": int(len(frame)),
            "mrr": float(np.mean(rr)), "ndcg_at_3": float(np.mean(ndcg3)),
            "top1": float(np.mean(rank.eq(1))), "top3": float(np.mean(rank.le(3) & retrieved)),
            "unretrieved_truth_n": int((~retrieved).sum()),
            "definition": "MRR assigns 0 when the true fine label is outside the retained candidate set; NDCG@3 uses one binary relevant label and ideal DCG=1.",
        })
    return pd.DataFrame(rows)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    ablation, donor, contrasts = projection_ablation()
    ols, ols_summary = ols_quality()
    ranking = ranking_metrics()
    ablation.to_csv(OUT / "network_projection_ablation_summary.csv", index=False)
    donor.to_csv(OUT / "network_projection_ablation_by_donor.csv", index=False)
    contrasts.to_csv(OUT / "network_projection_ablation_signflip.csv", index=False)
    ols.to_csv(OUT / "projector_ols_quality_summary.csv", index=False)
    (OUT / "projector_ols_quality_manifest.json").write_text(json.dumps(ols_summary, indent=2), encoding="utf-8")
    ranking.to_csv(OUT / "ranking_metrics_mrr_ndcg.csv", index=False)
    print(ablation.to_string(index=False))
    print(contrasts.to_string(index=False))
    print(ols.to_string(index=False))
    print(ranking.to_string(index=False))


if __name__ == "__main__":
    main()
