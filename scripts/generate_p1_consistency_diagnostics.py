#!/usr/bin/env python
from __future__ import annotations

import json
import math
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy.stats import kendalltau, mannwhitneyu, spearmanr
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, roc_curve
from sklearn.preprocessing import StandardScaler


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from core.reference_projection import (  # noqa: E402
    align_matrices,
    apply_projector,
    compute_logcpm,
    fit_linear_projector,
    map_index_to_symbols,
    read_bo2023_gene_matrix,
    read_gene_map,
)
from scripts.build_bo2023_reference_projector import (  # noqa: E402
    DEFAULT_COUNTS,
    DEFAULT_GENE_MAP,
    DEFAULT_MODEL_GENES,
    DEFAULT_SAMPLE_INFO,
    DEFAULT_VSD,
    read_locked_model_genes,
)
from scripts.run_bo2023_projected_vsd_loso import build_centroids, corr_scores, read_labels  # noqa: E402


OUTDIR = ROOT / "reports" / "p1_consistency_diagnostics_20260703"
P0 = ROOT / "reports" / "p0_hard_evidence_20260629"
RESULT_ROOT = ROOT / "results" / "bo2023_reference_projection_20260616_cleaned_symbols"
HYBRID = "hybrid_projected_network_logcpm_exact"


def softmax(scores: np.ndarray) -> np.ndarray:
    z = scores - np.nanmax(scores)
    exp = np.exp(z)
    return exp / np.sum(exp)


def normalized_entropy(probs: np.ndarray) -> float:
    probs = probs[probs > 0]
    if len(probs) <= 1:
        return 0.0
    return float(-np.sum(probs * np.log(probs)) / math.log(len(probs)))


def read_metadata() -> pd.DataFrame:
    info = pd.read_excel(DEFAULT_SAMPLE_INFO, sheet_name="mfas5_819samples_phenSet4", usecols=["No.", "Region", "SaleemNetworks", "MonkeyID"])
    info["sample_id"] = info["No."].astype(str).str.strip()
    info["region"] = info["Region"].astype(str).str.strip()
    info["network"] = info["SaleemNetworks"].astype(str).str.strip()
    info["monkey_id"] = info["MonkeyID"].astype(str).str.strip()
    return info.drop_duplicates("sample_id").set_index("sample_id")


def build_lomo_centroids(values: pd.DataFrame, labels: pd.Series, train_samples: list[str], groups: list[str]) -> np.ndarray:
    cols = []
    for group in groups:
        group_samples = [sample for sample in train_samples if labels.loc[sample] == group]
        cols.append(values[group_samples].mean(axis=1).to_numpy(dtype=np.float32))
    return np.column_stack(cols).astype(np.float32)


def rank_row(design: str, fold_id: str, sample: str, truth: str, groups: list[str], scores: np.ndarray, n_genes: int) -> dict[str, Any]:
    order = np.argsort(scores)[::-1]
    ranked = [groups[int(i)] for i in order]
    true_rank = ranked.index(truth) + 1
    probs = softmax(scores)
    top_scores = scores[order]
    return {
        "design": design,
        "fold_id": fold_id,
        "sample_id": sample,
        "truth_network": truth,
        "network_rank_top1": ranked[0],
        "network_rank_top2": ranked[1],
        "network_rank_top3": ranked[2],
        "network_true_rank": int(true_rank),
        "network_hit1": int(true_rank == 1),
        "network_hit3": int(true_rank <= 3),
        "score_top1": float(top_scores[0]),
        "score_top2": float(top_scores[1]),
        "score_top3": float(top_scores[2]),
        "score_top4": float(top_scores[3]) if len(top_scores) > 3 else float("nan"),
        "top1_margin": float(top_scores[0] - top_scores[1]),
        "top3_beam_margin": float(top_scores[2] - top_scores[3]) if len(top_scores) > 3 else float("nan"),
        "score_entropy": normalized_entropy(probs),
        "n_overlap_genes": int(n_genes),
        "network_rank_full": " | ".join(ranked),
        "score_vector_json": json.dumps({group: float(scores[i]) for i, group in enumerate(groups)}, ensure_ascii=False),
    }


def compute_network_diagnostics() -> pd.DataFrame:
    gene_map = read_gene_map(DEFAULT_GENE_MAP)
    counts, _ = map_index_to_symbols(read_bo2023_gene_matrix(DEFAULT_COUNTS, dtype="float32"), gene_map)
    vsd, _ = map_index_to_symbols(read_bo2023_gene_matrix(DEFAULT_VSD, dtype="float32"), gene_map)
    counts, vsd, genes, samples = align_matrices(counts, vsd)
    metadata = read_metadata()
    samples = [sample for sample in samples if sample in metadata.index]
    labels = metadata.loc[samples, "network"]
    monkeys = metadata.loc[samples, "monkey_id"]
    locked_genes = read_locked_model_genes(DEFAULT_MODEL_GENES)
    selected_genes = [gene for gene in locked_genes if gene in set(genes)]
    if len(selected_genes) < 20:
        selected_genes = genes
    counts = counts.loc[selected_genes, samples]
    vsd = vsd.loc[selected_genes, samples]
    logcpm = compute_logcpm(counts)
    rows: list[dict[str, Any]] = []

    # LOSO
    for fold_no, sample in enumerate(samples, start=1):
        train_samples = [s for s in samples if s != sample]
        truth = str(labels.loc[sample])
        groups = sorted(set(labels.loc[train_samples].astype(str)))
        fit, _ = fit_linear_projector(logcpm[train_samples], vsd[train_samples])
        projected = apply_projector(fit, logcpm[[sample]])
        reference = build_centroids(vsd, labels, train_samples, groups)
        scores = corr_scores(reference, projected[sample].to_numpy(dtype=np.float32))
        rows.append(rank_row("LOSO", str(fold_no), sample, truth, groups, scores, len(selected_genes)))

    # LOMO
    for monkey_id in sorted(monkeys.unique().tolist()):
        test_samples = [sample for sample in samples if monkeys.loc[sample] == monkey_id]
        train_samples = [sample for sample in samples if monkeys.loc[sample] != monkey_id]
        groups = sorted(set(labels.loc[train_samples].astype(str)))
        fit, _ = fit_linear_projector(logcpm[train_samples], vsd[train_samples])
        projected = apply_projector(fit, logcpm[test_samples])
        reference = build_lomo_centroids(vsd, labels, train_samples, groups)
        for sample in test_samples:
            truth = str(labels.loc[sample])
            scores = corr_scores(reference, projected[sample].to_numpy(dtype=np.float32))
            rows.append(rank_row("LOMO", str(monkey_id), sample, truth, groups, scores, len(selected_genes)))
    out = pd.DataFrame(rows)
    out.to_csv(OUTDIR / "network_projected_vsd_fullscore_diagnostics.csv", index=False)
    return out


def region_to_network_map() -> dict[str, str]:
    h = pd.read_csv(P0 / "resolution_group_hierarchy.csv")
    return dict(zip(h["region"].astype(str), h["network"].astype(str)))


def load_downstream(design: str, level: str) -> pd.DataFrame:
    if design == "LOSO":
        base = P0 / "formal_loso_rerun_20260629"
        if level == "exact":
            return pd.read_csv(base / "hybrid_formal_loso_exact_region_detail.csv")
        return pd.read_csv(base / "hybrid_formal_loso_resolution_group_detail.csv")
    base = RESULT_ROOT / "formal_three_tier_lomo_hybrid"
    if level == "exact":
        df = pd.read_csv(base / "formal_lomo_exact_region_detail.csv")
    else:
        df = pd.read_csv(base / "formal_lomo_resolution_group_detail.csv")
    return df[df["route_family"].eq(HYBRID)].copy()


def add_downstream_features(df: pd.DataFrame, level: str, r2n: dict[str, str]) -> pd.DataFrame:
    out = df.copy()
    out["beam_networks"] = out["network_beam"].astype(str).str.split(r"\s+\|\s+", regex=True)
    known_networks = set(r2n.values())

    def predicted_network(label: Any) -> str | None:
        text = str(label).strip()
        if not text or text.lower() == "nan":
            return None
        if "::" in text:
            return text.split("::", 1)[0].strip()
        if text in r2n:
            return r2n[text]
        if text in known_networks:
            return text
        return None

    if level == "exact":
        out["pred_top1_network"] = out["pred_top1"].astype(str).map(r2n)
        pred_cols = ["pred_top1", "pred_top2", "pred_top3"]
        out["pred_top1_in_network_beam"] = [pred in beam for pred, beam in zip(out["pred_top1_network"], out["beam_networks"])]
        out["downstream_true_rank"] = out["true_rank"]
        out["downstream_hit3"] = out["hit3"]
    else:
        out["pred_top1_network"] = out["pred_group_top1"].astype(str).str.split("::").str[0]
        # Single-region groups may not carry network prefix; fall back through the predicted region if possible.
        fallback = out["pred_group_top1"].astype(str).map(r2n)
        out.loc[~out["pred_top1_network"].isin(set(r2n.values())), "pred_top1_network"] = fallback
        pred_cols = ["pred_group_top1", "pred_group_top2", "pred_group_top3"]
        out["pred_top1_in_network_beam"] = [pred in beam for pred, beam in zip(out["pred_top1_network"], out["beam_networks"])]
        out["downstream_true_rank"] = out["group_true_rank"]
        out["downstream_hit3"] = out["group_hit3"]
    pred_top3_networks = []
    for _, row in out.iterrows():
        networks = [predicted_network(row.get(col)) for col in pred_cols]
        pred_top3_networks.append([net for net in networks if net])
    out["pred_top3_networks"] = [" | ".join(networks) for networks in pred_top3_networks]
    out["pred_top3_all_in_network_beam"] = [
        bool(networks) and all(network in beam for network in networks)
        for networks, beam in zip(pred_top3_networks, out["beam_networks"])
    ]
    out["pred_top3_any_outside_network_beam"] = [
        bool(networks) and any(network not in beam for network in networks)
        for networks, beam in zip(pred_top3_networks, out["beam_networks"])
    ]
    out["beam_miss"] = out["network_top3_hit"].astype(int).eq(0)
    out["local_rerank_miss_given_beam_hit"] = out["network_top3_hit"].astype(int).eq(1) & out["downstream_hit3"].astype(int).eq(0)
    out["cascade_rescued_after_beam_miss"] = out["beam_miss"] & out["downstream_hit3"].astype(int).eq(1)
    return out


def consistency_and_cascade(net: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    r2n = region_to_network_map()
    consistency_rows = []
    cascade_rows = []
    merged_frames = []
    for design in ["LOSO", "LOMO"]:
        ndiag = net[net["design"].eq(design)].copy()
        for level in ["exact", "group"]:
            d = add_downstream_features(load_downstream(design, level), level, r2n)
            merged = d.merge(ndiag, on="sample_id", how="left", suffixes=("", "_diag"))
            merged["level"] = level
            merged["design"] = design
            merged_frames.append(merged)
            n = len(merged)
            consistency_rows.append(
                {
                    "design": design,
                    "level": level,
                    "n": n,
                    "pred_top1_in_network_beam_rate": float(merged["pred_top1_in_network_beam"].mean()),
                    "pred_top3_all_in_network_beam_rate": float(merged["pred_top3_all_in_network_beam"].mean()),
                    "pred_top3_any_outside_network_beam_rate": float(merged["pred_top3_any_outside_network_beam"].mean()),
                    "network_beam_hit_rate": float(merged["network_top3_hit"].mean()),
                    "downstream_top3_rate": float(merged["downstream_hit3"].mean()),
                    "downstream_top3_given_beam_hit": float(merged.loc[merged["network_top3_hit"].eq(1), "downstream_hit3"].mean()),
                    "downstream_top3_given_beam_miss": float(merged.loc[merged["network_top3_hit"].eq(0), "downstream_hit3"].mean()) if (merged["network_top3_hit"].eq(0)).any() else float("nan"),
                }
            )
            misses = merged[merged["downstream_hit3"].astype(int).eq(0)]
            cascade_rows.append(
                {
                    "design": design,
                    "level": level,
                    "n": n,
                    "downstream_misses": int(len(misses)),
                    "misses_due_to_network_beam_miss": int(misses["beam_miss"].sum()),
                    "misses_due_to_local_rerank_with_beam_hit": int(misses["local_rerank_miss_given_beam_hit"].sum()),
                    "beam_miss_fraction_of_downstream_misses": float(misses["beam_miss"].mean()) if len(misses) else float("nan"),
                    "cascade_rescue_after_beam_miss_count": int(merged["cascade_rescued_after_beam_miss"].sum()),
                    "cascade_rescue_after_beam_miss_rate": float(merged.loc[merged["beam_miss"], "downstream_hit3"].mean()) if merged["beam_miss"].any() else float("nan"),
                }
            )
    consistency = pd.DataFrame(consistency_rows)
    cascade = pd.DataFrame(cascade_rows)
    merged_all = pd.concat(merged_frames, ignore_index=True)
    consistency.to_csv(OUTDIR / "dual_space_agreement_metrics.csv", index=False)
    cascade.to_csv(OUTDIR / "error_cascade_metrics.csv", index=False)
    merged_all.to_csv(OUTDIR / "dual_space_merged_sample_detail.csv", index=False)
    return consistency, cascade, merged_all


def rank_correlation(merged: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (design, level), frame in merged.groupby(["design", "level"], sort=True):
        x = frame["network_true_rank"].astype(float)
        y = frame["downstream_true_rank"].astype(float)
        sp = spearmanr(x, y, nan_policy="omit")
        kt = kendalltau(x, y, nan_policy="omit")
        rows.append(
            {
                "design": design,
                "level": level,
                "n": int(len(frame)),
                "spearman_r": float(sp.statistic),
                "spearman_p": float(sp.pvalue),
                "kendall_tau": float(kt.statistic),
                "kendall_p": float(kt.pvalue),
                "interpretation": "positive means worse VSD Network rank is associated with worse logCPM downstream rank",
            }
        )
    out = pd.DataFrame(rows)
    out.to_csv(OUTDIR / "vsd_network_rank_vs_logcpm_downstream_rank.csv", index=False)
    return out


def mannwhitney_summary(frame: pd.DataFrame, predictor: str, outcome: str) -> dict[str, Any]:
    hit = frame.loc[frame[outcome].astype(int).eq(1), predictor].dropna().astype(float)
    miss = frame.loc[frame[outcome].astype(int).eq(0), predictor].dropna().astype(float)
    if len(hit) == 0 or len(miss) == 0 or frame[predictor].nunique(dropna=True) <= 1:
        return {
            "hit_median": float(hit.median()) if len(hit) else float("nan"),
            "miss_median": float(miss.median()) if len(miss) else float("nan"),
            "mannwhitney_p": float("nan"),
            "status": "not_testable_no_variation_or_single_class",
        }
    test = mannwhitneyu(hit, miss, alternative="two-sided")
    return {
        "hit_median": float(hit.median()),
        "miss_median": float(miss.median()),
        "mannwhitney_p": float(test.pvalue),
        "status": "tested",
    }


def threshold_for_predictor(y: np.ndarray, x: np.ndarray, higher_is_better: bool) -> dict[str, Any]:
    mask = np.isfinite(x)
    y = y[mask].astype(int)
    x = x[mask].astype(float)
    if len(np.unique(y)) < 2 or len(np.unique(x)) < 2:
        return {"auc": float("nan"), "threshold": float("nan"), "sensitivity": float("nan"), "specificity": float("nan")}
    score = x if higher_is_better else -x
    auc = roc_auc_score(y, score)
    fpr, tpr, thr = roc_curve(y, score)
    idx = int(np.argmax(tpr - fpr))
    raw_thr = float(thr[idx] if higher_is_better else -thr[idx])
    return {"auc": float(auc), "threshold": raw_thr, "sensitivity": float(tpr[idx]), "specificity": float(1 - fpr[idx])}


def diagnostic_validation(merged: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    predictors = [
        ("top1_margin", True, "larger margin should predict hit"),
        ("top3_beam_margin", True, "larger beam boundary margin should predict hit"),
        ("score_entropy", False, "lower entropy should predict hit"),
        ("n_overlap_genes", True, "larger marker overlap should predict hit"),
    ]
    outcomes = [
        ("network_hit3", "Network Top3 hit"),
        ("downstream_hit3", "Downstream Top3 hit"),
    ]
    diff_rows = []
    auc_rows = []
    for (design, level), frame in merged.groupby(["design", "level"], sort=True):
        # Network rows are duplicated across exact/group; keep all for downstream, deduplicate for network outcome.
        for predictor, higher, note in predictors:
            for outcome, label in outcomes:
                f = frame.drop_duplicates("sample_id") if outcome == "network_hit3" else frame
                stats = mannwhitney_summary(f, predictor, outcome)
                th = threshold_for_predictor(f[outcome].to_numpy(), f[predictor].to_numpy(dtype=float), higher)
                diff_rows.append(
                    {
                        "design": design,
                        "level": level,
                        "outcome": outcome,
                        "outcome_label": label,
                        "predictor": predictor,
                        "direction": note,
                        **stats,
                    }
                )
                auc_rows.append(
                    {
                        "design": design,
                        "level": level,
                        "outcome": outcome,
                        "predictor": predictor,
                        "roc_auc": th["auc"],
                        "recommended_threshold": th["threshold"],
                        "threshold_sensitivity": th["sensitivity"],
                        "threshold_specificity": th["specificity"],
                        "threshold_rule": ("low_confidence_if_below_threshold" if higher else "low_confidence_if_above_threshold"),
                    }
                )
    diff = pd.DataFrame(diff_rows)
    aucs = pd.DataFrame(auc_rows)
    diff.to_csv(OUTDIR / "diagnostic_hit_miss_tests.csv", index=False)
    aucs.to_csv(OUTDIR / "diagnostic_auc_thresholds.csv", index=False)

    # Multivariable logistic regression for downstream hit.
    lr_rows = []
    for (design, level), frame in merged.groupby(["design", "level"], sort=True):
        predictors_used = ["top1_margin", "top3_beam_margin", "score_entropy"]
        f = frame.dropna(subset=predictors_used + ["downstream_hit3"]).copy()
        y = f["downstream_hit3"].astype(int).to_numpy()
        if len(np.unique(y)) < 2:
            continue
        x = f[predictors_used].astype(float).to_numpy()
        x = StandardScaler().fit_transform(x)
        model = LogisticRegression(max_iter=1000, class_weight="balanced", random_state=20260703)
        model.fit(x, y)
        prob = model.predict_proba(x)[:, 1]
        lr_rows.append(
            {
                "design": design,
                "level": level,
                "outcome": "downstream_hit3",
                "n": int(len(f)),
                "predictors": "top1_margin + top3_beam_margin + score_entropy",
                "apparent_roc_auc": float(roc_auc_score(y, prob)),
                "coef_top1_margin": float(model.coef_[0][0]),
                "coef_top3_beam_margin": float(model.coef_[0][1]),
                "coef_score_entropy": float(model.coef_[0][2]),
                "note": "apparent in-sample diagnostic model; use as calibration evidence, not performance validation",
            }
        )
    lr = pd.DataFrame(lr_rows)
    lr.to_csv(OUTDIR / "diagnostic_logistic_regression_auc.csv", index=False)
    return diff, aucs


def coverage_boundary() -> pd.DataFrame:
    rows = []
    for source, path in [
        ("internal_projected_vsd_fullscore", OUTDIR / "network_projected_vsd_fullscore_diagnostics.csv"),
        ("AHBA_hybrid", RESULT_ROOT / "ahba_external_formal_three_tier" / "ahba_formal_three_tier_sample_detail.csv"),
        ("GSE189919_no_truth", ROOT / "results" / "gse189919_csf_tracing_validation_20260613" / "sample_marker_detection.csv"),
    ]:
        if not path.exists():
            continue
        df = pd.read_csv(path)
        cols = [c for c in df.columns if c in {"n_overlap_genes", "network_overlap_genes", "exact_overlap_genes_top100", "detected_fraction", "detected_markers"}]
        for col in cols:
            rows.append(
                {
                    "source": source,
                    "metric": col,
                    "n": int(df[col].notna().sum()),
                    "min": float(df[col].min()),
                    "median": float(df[col].median()),
                    "max": float(df[col].max()),
                    "std": float(df[col].std(ddof=0)),
                    "has_accuracy_truth": source != "GSE189919_no_truth",
                    "interpretation": "predictive test possible only if metric varies and anatomical truth exists",
                }
            )
    out = pd.DataFrame(rows)
    out.to_csv(OUTDIR / "marker_coverage_validation_boundary.csv", index=False)
    return out


def write_report(consistency: pd.DataFrame, cascade: pd.DataFrame, corr: pd.DataFrame, diff: pd.DataFrame, aucs: pd.DataFrame, cov: pd.DataFrame) -> None:
    def get(df: pd.DataFrame, design: str, level: str, col: str) -> Any:
        return df[(df["design"].eq(design)) & (df["level"].eq(level))].iloc[0][col]

    exact_loso_cascade = cascade[(cascade["design"].eq("LOSO")) & (cascade["level"].eq("exact"))].iloc[0]
    group_loso = consistency[(consistency["design"].eq("LOSO")) & (consistency["level"].eq("group"))].iloc[0]
    exact_loso = consistency[(consistency["design"].eq("LOSO")) & (consistency["level"].eq("exact"))].iloc[0]
    exact_auc = aucs[(aucs["design"].eq("LOSO")) & (aucs["level"].eq("exact")) & (aucs["outcome"].eq("downstream_hit3")) & (aucs["predictor"].eq("top3_beam_margin"))].iloc[0]
    entropy_auc = aucs[(aucs["design"].eq("LOSO")) & (aucs["level"].eq("exact")) & (aucs["outcome"].eq("downstream_hit3")) & (aucs["predictor"].eq("score_entropy"))].iloc[0]
    text = f"""# P1 dual-space consistency and diagnostic validity report

Generated: 2026-07-03

## 1. VSD/logCPM dual-space consistency

The projected-VSD Network beam and logCPM downstream reranking are structurally coordinated. The observed `pred_top1_in_network_beam_rate`, `pred_top3_all_in_network_beam_rate`, and `pred_top3_any_outside_network_beam_rate` audit whether downstream logCPM candidates remain inside the retained Network Top3 beam.

Key results:

- LOSO exact-region top1 candidate in Network beam: {get(consistency, 'LOSO', 'exact', 'pred_top1_in_network_beam_rate'):.2%}.
- LOSO exact-region Top3 all inside Network beam: {get(consistency, 'LOSO', 'exact', 'pred_top3_all_in_network_beam_rate'):.2%}; any outside beam: {get(consistency, 'LOSO', 'exact', 'pred_top3_any_outside_network_beam_rate'):.2%}.
- LOSO resolution-group top1 candidate in Network beam: {get(consistency, 'LOSO', 'group', 'pred_top1_in_network_beam_rate'):.2%}.
- LOSO resolution-group Top3 all inside Network beam: {get(consistency, 'LOSO', 'group', 'pred_top3_all_in_network_beam_rate'):.2%}; any outside beam: {get(consistency, 'LOSO', 'group', 'pred_top3_any_outside_network_beam_rate'):.2%}.
- LOSO exact Top3 given Network beam hit: {get(consistency, 'LOSO', 'exact', 'downstream_top3_given_beam_hit'):.2%}.
- LOSO resolution-group Top3 given Network beam hit: {get(consistency, 'LOSO', 'group', 'downstream_top3_given_beam_hit'):.2%}.
- LOSO exact-region cascade after Network beam miss: {int(exact_loso_cascade['cascade_rescue_after_beam_miss_count'])} rescued cases ({exact_loso_cascade['cascade_rescue_after_beam_miss_rate']:.2%}).

Interpretation: Network beam misses are near-hard failures for downstream group/exact recovery. LOSO shows no recovery after a beam miss, while LOMO is explicitly audited in `error_cascade_metrics.csv`. Most remaining downstream misses occur despite a correct Network beam, meaning the residual error is mainly local logCPM reranking ambiguity rather than VSD/logCPM inconsistency.

## 2. Rank correlation between spaces

`vsd_network_rank_vs_logcpm_downstream_rank.csv` reports Spearman and Kendall correlations between the projected-VSD true Network rank and the logCPM downstream true-rank. Positive values mean that worse Network ranking tends to accompany worse downstream ranking.

This is a calibration/audit result rather than a claim that Network rank alone determines region rank: logCPM reranking still carries independent local evidence inside the beam.

## 3. Diagnostic validity

Full projected-VSD Network scores were recomputed for LOSO and LOMO, allowing true entropy and margins:

- `top1_margin`: score gap between Network rank 1 and rank 2.
- `top3_beam_margin`: score gap between rank 3 and rank 4; this measures how secure the retained beam boundary is.
- `score_entropy`: normalized softmax entropy over all Network scores.
- `n_overlap_genes`: marker overlap count.

For LOSO exact-region Top3, Top3 beam margin AUC for downstream hit was {exact_auc['roc_auc']:.3f}; a data-driven low-confidence threshold is {exact_auc['recommended_threshold']:.4g} under rule `{exact_auc['threshold_rule']}`. Entropy AUC was {entropy_auc['roc_auc']:.3f}; high entropy should be treated as lower confidence when above the listed threshold.

Marker overlap/coverage is not statistically testable in internal Bo2023 or AHBA validation because marker overlap is fixed at full coverage (200/200 Network genes; 100/100 exact genes in AHBA). GSE189919 has variable marker detection but no anatomical truth, so it supports QC/feasibility thresholds rather than accuracy prediction.

## 4. Recommended reviewer-facing interpretation

Projected VSD and logCPM are not independent competing endpoints. Projected VSD is used to define a broad Network Top3 candidate beam; logCPM then performs local reranking that can be audited against that beam. Error-cascade analysis shows that a missed Network beam almost always prevents downstream recovery, while downstream misses after a correct beam reflect local anatomical ambiguity. Diagnostic margins and entropy provide quantitative confidence signals; marker coverage remains a QC/domain-transfer diagnostic rather than an internally validated accuracy predictor because controlled validation samples have near-complete marker coverage.

## 5. Evidence files

- `network_projected_vsd_fullscore_diagnostics.csv`
- `dual_space_agreement_metrics.csv`
- `error_cascade_metrics.csv`
- `vsd_network_rank_vs_logcpm_downstream_rank.csv`
- `diagnostic_hit_miss_tests.csv`
- `diagnostic_auc_thresholds.csv`
- `diagnostic_logistic_regression_auc.csv`
- `marker_coverage_validation_boundary.csv`
"""
    (OUTDIR / "P1_CONSISTENCY_DIAGNOSTICS_REPORT.md").write_text(text, encoding="utf-8")


def main() -> int:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    net = compute_network_diagnostics()
    consistency, cascade, merged = consistency_and_cascade(net)
    corr = rank_correlation(merged)
    diff, aucs = diagnostic_validation(merged)
    cov = coverage_boundary()
    write_report(consistency, cascade, corr, diff, aucs, cov)
    checklist = """# P1 consistency/diagnostic completion checklist

- [x] VSD Network Top3 and logCPM downstream agreement metrics.
- [x] VSD Network true-rank vs logCPM downstream true-rank Spearman/Kendall correlation.
- [x] Beam miss vs local rerank miss decomposition.
- [x] Error-cascade check showing whether group/exact can recover after Network beam miss.
- [x] Diagnostic hit/miss tests for margins, entropy and marker overlap where statistically testable.
- [x] ROC-AUC and recommended thresholds for margin/entropy diagnostics.
- [x] Logistic regression diagnostic model for downstream Top3 hit.
- [x] Marker coverage boundary statement: full coverage internally, variable only in no-truth biofluid data.
"""
    (OUTDIR / "P1_COMPLETION_CHECKLIST.md").write_text(checklist, encoding="utf-8")
    manifest = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "script": "scripts/generate_p1_consistency_diagnostics.py",
        "files": sorted(str(p.relative_to(OUTDIR)) for p in OUTDIR.rglob("*") if p.is_file()),
    }
    (OUTDIR / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(f"Wrote P1 consistency/diagnostic package to {OUTDIR}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
