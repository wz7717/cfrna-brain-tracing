#!/usr/bin/env python
from __future__ import annotations

import argparse
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

from scripts.run_bo2023_loso_validation import read_vsd_matrix  # noqa: E402
from scripts.run_bo2023_network_correlation_validation import (  # noqa: E402
    build_group_reference,
    select_group_discriminative_genes,
)
from scripts.run_bo2023_v2_loso_validation import (  # noqa: E402
    DEFAULT_GENE_MAP,
    DEFAULT_MATRIX,
    DEFAULT_SAMPLE_INFO,
    map_matrix_to_symbols,
)


TRUTH_FILE = (
    ROOT
    / "results"
    / "brats_tcga_lgg_65_mri_truth_evaluation_20260609"
    / "brats_tcga_lgg_65_mri_truth_and_predictions.csv"
)
TCGA_MATRIX = (
    ROOT
    / "data"
    / "tcga_brain_tumor_expression"
    / "tcga_gbm_lgg_primary_tumor_tpm_unstranded_sample_mean.tsv"
)
TCGA_SUMMARY = (
    ROOT
    / "results"
    / "tcga_gbm_lgg_sample_mri_label_tracing_20260605"
    / "tcga_gbm_lgg_sample_mri_label_tracing_summary.csv"
)
OUTDIR = ROOT / "results" / "tcga_lgg_65_tumor_adapted_network_20260612"
HALLMARK_FILE = ROOT / "data" / "gene_sets" / "MSigDB_Hallmark_2020.txt"
MODEL_OUT = ROOT / "data" / "models" / "bo2023_tcga_tumor_adapted_network_model.npz"
ORTHOLOGY_FILE = ROOT / "data" / "orthology" / "ensembl_mfascicularis_hsapiens_homology.tsv"
SYMBOL_RE = re.compile(r"^[A-Z][A-Z0-9.-]{1,30}$")
TUMOR_STATE_HALLMARKS = {
    "TNF-alpha Signaling via NF-kB",
    "Hypoxia",
    "Cholesterol Homeostasis",
    "Mitotic Spindle",
    "TGF-beta Signaling",
    "IL-6/JAK/STAT3 Signaling",
    "DNA Repair",
    "G2-M Checkpoint",
    "Apoptosis",
    "Notch Signaling",
    "Interferon Alpha Response",
    "Interferon Gamma Response",
    "Complement",
    "Unfolded Protein Response",
    "PI3K/AKT/mTOR  Signaling",
    "mTORC1 Signaling",
    "E2F Targets",
    "Myc Targets V1",
    "Myc Targets V2",
    "Epithelial Mesenchymal Transition",
    "Inflammatory Response",
    "Glycolysis",
    "Reactive Oxygen Species Pathway",
    "p53 Pathway",
    "Angiogenesis",
    "Coagulation",
    "IL-2/STAT5 Signaling",
    "Allograft Rejection",
    "KRAS Signaling Up",
}


def corr_scores(reference: np.ndarray, samples: np.ndarray) -> np.ndarray:
    ref = np.asarray(reference, dtype=float)
    x = np.asarray(samples, dtype=float)
    ref0 = ref - ref.mean(axis=0, keepdims=True)
    x0 = x - x.mean(axis=0, keepdims=True)
    numerator = ref0.T @ x0
    denominator = np.sqrt(
        np.square(ref0).sum(axis=0)[:, None] * np.square(x0).sum(axis=0)[None, :] + 1e-12
    )
    return np.nan_to_num(numerator / denominator)


def percentile_ranks(values: np.ndarray) -> np.ndarray:
    order = np.argsort(values, axis=0, kind="mergesort")
    ranks = np.empty_like(order, dtype=float)
    base = np.linspace(0.0, 1.0, values.shape[0], dtype=float)
    for column in range(values.shape[1]):
        ranks[order[:, column], column] = base
    return ranks


def quantile_map_to_vsd(tumor_values: np.ndarray, monkey_values: np.ndarray) -> np.ndarray:
    target = np.median(np.sort(monkey_values, axis=0), axis=1)
    mapped = np.empty_like(tumor_values, dtype=float)
    for column in range(tumor_values.shape[1]):
        order = np.argsort(tumor_values[:, column], kind="mergesort")
        mapped[order, column] = target
    return mapped


def valid_human_symbol(symbol: str) -> bool:
    text = str(symbol).strip()
    return bool(SYMBOL_RE.fullmatch(text)) and not text.startswith(("ENSMFAG", "RF0"))


def load_tumor_state_blacklist(path: Path) -> set[str]:
    genes: set[str] = set()
    for line in path.read_text(encoding="utf-8").splitlines():
        parts = line.rstrip().split("\t")
        if parts and parts[0] in TUMOR_STATE_HALLMARKS:
            genes.update(item.strip() for item in parts[2:] if item.strip())
    return genes


def map_monkey_to_one2one_human(matrix: pd.DataFrame, orthology_path: Path) -> pd.DataFrame:
    orthology = pd.read_csv(orthology_path, sep="\t", dtype=str)
    confidence_column = "Human orthology confidence [0 low, 1 high]"
    orthology = orthology[
        orthology["Human homology type"].eq("ortholog_one2one")
        & orthology[confidence_column].eq("1")
        & orthology["Human gene name"].notna()
    ].copy()
    orthology["Gene stable ID"] = orthology["Gene stable ID"].astype(str).str.strip()
    orthology["Human gene name"] = orthology["Human gene name"].astype(str).str.strip()
    symbol_by_id = (
        orthology.drop_duplicates("Gene stable ID")
        .set_index("Gene stable ID")["Human gene name"]
    )
    symbols = matrix.index.to_series().map(symbol_by_id)
    keep = symbols.notna() & symbols.map(valid_human_symbol)
    conserved = matrix.loc[keep].groupby(symbols.loc[keep], sort=True).mean()
    conserved.index.name = "gene_symbol"
    return conserved


def load_tcga_subset(path: Path, sample_ids: list[str], genes: set[str]) -> pd.DataFrame:
    frames: list[pd.DataFrame] = []
    usecols = ["gene_symbol", *sample_ids]
    for chunk in pd.read_csv(path, sep="\t", usecols=usecols, chunksize=5000):
        keep = chunk["gene_symbol"].astype(str).isin(genes)
        if keep.any():
            frames.append(chunk.loc[keep])
    if not frames:
        raise ValueError("No candidate genes found in TCGA matrix")
    result = pd.concat(frames, ignore_index=True)
    result["gene_symbol"] = result["gene_symbol"].astype(str)
    return result.groupby("gene_symbol", sort=False)[sample_ids].mean()


def tumor_pc1_loadings(values: np.ndarray) -> np.ndarray:
    x = np.log1p(np.clip(values, 0, None))
    x = x - x.mean(axis=1, keepdims=True)
    scale = x.std(axis=1, keepdims=True)
    scale[scale < 1e-8] = 1.0
    z = x / scale
    u, _, _ = np.linalg.svd(z, full_matrices=False)
    return np.abs(u[:, 0])


def candidate_truths(value: Any) -> set[str]:
    return {item.strip() for item in str(value or "").split("|") if item.strip()}


def rank_predictions(scores: np.ndarray, groups: list[str]) -> tuple[list[str], list[str]]:
    order = np.argsort(scores, axis=0)[::-1]
    top1 = [groups[int(order[0, i])] for i in range(scores.shape[1])]
    top3 = [" | ".join(groups[int(j)] for j in order[:3, i]) for i in range(scores.shape[1])]
    return top1, top3


def nested_bias_prior_calibration(
    scores: np.ndarray,
    labels: np.ndarray,
    groups: list[str],
) -> tuple[np.ndarray, pd.DataFrame]:
    n = scores.shape[1]
    calibrated = np.empty_like(scores)
    choices: list[dict[str, Any]] = []
    bias_strengths = (0.0, 0.5, 0.75, 1.0, 1.25)
    prior_strengths = (0.0, 0.1, 0.25, 0.5)

    def offsets(train_idx: np.ndarray, bias_strength: float, prior_strength: float) -> np.ndarray:
        bias = scores[:, train_idx].mean(axis=1)
        counts = np.asarray([(labels[train_idx] == group).sum() for group in groups], dtype=float)
        prior = (counts + 1.0) / (counts.sum() + len(groups))
        return -bias_strength * bias + prior_strength * np.log((1.0 / len(groups)) / prior)

    def macro_recall(predictions: list[str], truths: np.ndarray) -> float:
        recalls = []
        for group in sorted(set(truths)):
            mask = truths == group
            recalls.append(np.mean(np.asarray(predictions, dtype=str)[mask] == group))
        return float(np.mean(recalls)) if recalls else 0.0

    for heldout in range(n):
        outer_train = np.asarray([i for i in range(n) if i != heldout], dtype=int)
        best = None
        for bias_strength in bias_strengths:
            for prior_strength in prior_strengths:
                hits = []
                predictions = []
                truths = []
                for inner in outer_train:
                    inner_train = outer_train[outer_train != inner]
                    adjusted = scores[:, inner] + offsets(inner_train, bias_strength, prior_strength)
                    prediction = groups[int(np.argmax(adjusted))]
                    predictions.append(prediction)
                    truths.append(labels[inner])
                    hits.append(prediction == labels[inner])
                balanced = macro_recall(predictions, np.asarray(truths, dtype=str))
                result = (balanced, float(np.mean(hits)), -(bias_strength + prior_strength))
                if best is None or result > best[0]:
                    best = (result, bias_strength, prior_strength)
        assert best is not None
        _, bias_strength, prior_strength = best
        calibrated[:, heldout] = scores[:, heldout] + offsets(outer_train, bias_strength, prior_strength)
        choices.append(
            {
                "heldout_index": heldout,
                "bias_strength": bias_strength,
                "prior_strength": prior_strength,
                "inner_macro_recall": best[0][0],
                "inner_top1_accuracy": best[0][1],
            }
        )
    return calibrated, pd.DataFrame(choices)


def monkey_loo_scores(values: np.ndarray, labels: np.ndarray, groups: list[str]) -> np.ndarray:
    scores = np.empty((len(groups), values.shape[1]), dtype=float)
    sums = {group: values[:, labels == group].sum(axis=1) for group in groups}
    counts = Counter(labels)
    for index in range(values.shape[1]):
        columns = []
        for group in groups:
            total = sums[group] - (values[:, index] if labels[index] == group else 0.0)
            count = counts[group] - int(labels[index] == group)
            columns.append(total / count)
        scores[:, index] = corr_scores(np.column_stack(columns), values[:, index : index + 1])[:, 0]
    return scores


def summarize_route(
    frame: pd.DataFrame,
    route: str,
    top1_col: str,
    top3_col: str,
    accepted_col: str | None = None,
) -> dict[str, Any]:
    accepted = frame[accepted_col].astype(bool) if accepted_col else pd.Series(True, index=frame.index)
    strict_truth = frame["whole_tumor_network_dominant"].astype(str)
    tolerant_truth = frame["whole_tumor_network_candidates"].map(candidate_truths)
    top1 = frame[top1_col].astype(str)
    top3 = frame[top3_col].map(candidate_truths)
    strict1 = top1.eq(strict_truth)
    tolerant1 = [prediction in truth for prediction, truth in zip(top1, tolerant_truth)]
    strict3 = [truth in predictions for truth, predictions in zip(strict_truth, top3)]
    tolerant3 = [bool(predictions & truth) for predictions, truth in zip(top3, tolerant_truth)]
    metrics = {
        "route": route,
        "n_total": int(len(frame)),
        "n_accepted": int(accepted.sum()),
        "coverage": float(accepted.mean()),
    }
    for name, values in (
        ("strict_top1", strict1),
        ("tolerant_top1", tolerant1),
        ("strict_top3", strict3),
        ("tolerant_top3", tolerant3),
    ):
        series = pd.Series(values, index=frame.index).loc[accepted]
        metrics[name] = float(series.mean()) if len(series) else None
        metrics[f"{name}_hits"] = int(series.sum()) if len(series) else 0
    return metrics


def main() -> int:
    parser = argparse.ArgumentParser(description="Build and evaluate a tumor-adapted TCGA-LGG Network route.")
    parser.add_argument("--truth", type=Path, default=TRUTH_FILE)
    parser.add_argument("--tcga-matrix", type=Path, default=TCGA_MATRIX)
    parser.add_argument("--tcga-summary", type=Path, default=TCGA_SUMMARY)
    parser.add_argument("--monkey-matrix", type=Path, default=DEFAULT_MATRIX)
    parser.add_argument("--sample-info", type=Path, default=DEFAULT_SAMPLE_INFO)
    parser.add_argument("--gene-map", type=Path, default=DEFAULT_GENE_MAP)
    parser.add_argument("--orthology", type=Path, default=ORTHOLOGY_FILE)
    parser.add_argument("--candidate-genes", type=int, default=1200)
    parser.add_argument("--final-genes", type=int, default=200)
    parser.add_argument("--rank-shift-quantile", type=float, default=0.75)
    parser.add_argument("--pc-loading-quantile", type=float, default=0.90)
    parser.add_argument("--ood-quantile", type=float, default=0.05)
    parser.add_argument("--outdir", type=Path, default=OUTDIR)
    parser.add_argument("--hallmark-file", type=Path, default=HALLMARK_FILE)
    parser.add_argument("--disable-hallmark-filter", action="store_true")
    parser.add_argument("--model-out", type=Path, default=MODEL_OUT)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    truth = pd.read_csv(args.truth)
    truth["patient_barcode"] = truth["patient_barcode"].astype(str).str.upper()
    sample_ids = truth["sample_id"].astype(str).tolist()
    tcga_summary = pd.read_csv(args.tcga_summary)
    lgg_ids = (
        tcga_summary.loc[tcga_summary["project_id"].eq("TCGA-LGG"), "sample_id"]
        .astype(str)
        .drop_duplicates()
        .tolist()
    )

    raw_monkey = read_vsd_matrix(args.monkey_matrix)
    monkey = map_monkey_to_one2one_human(raw_monkey, args.orthology)
    annotations = pd.read_excel(args.sample_info, sheet_name="mfas5_819samples_phenSet4")
    annotations["sample_id"] = annotations["No."].astype(str).str.strip()
    annotations["label"] = annotations["SaleemNetworks"].fillna("NA").astype(str).str.strip()
    annotations = annotations[annotations["sample_id"].isin(monkey.columns)].copy()
    monkey = monkey.loc[:, annotations["sample_id"].tolist()]
    labels = annotations.set_index("sample_id").reindex(monkey.columns)["label"].to_numpy(dtype=str)
    groups = sorted(set(labels))
    values = monkey.to_numpy(dtype=np.float32)
    reference, training = build_group_reference(values, labels, groups, heldout_idx=-1)
    candidate_idx, candidate_audit = select_group_discriminative_genes(
        values, groups, training, args.candidate_genes
    )
    candidate_symbols = monkey.index.to_numpy(dtype=str)[candidate_idx]
    valid_mask = np.asarray([valid_human_symbol(gene) for gene in candidate_symbols])
    candidate_idx = candidate_idx[valid_mask]
    candidate_symbols = candidate_symbols[valid_mask]
    candidate_audit = candidate_audit.loc[valid_mask].reset_index(drop=True)

    tumor_all = load_tcga_subset(args.tcga_matrix, lgg_ids, set(candidate_symbols))
    present = np.asarray([gene in tumor_all.index for gene in candidate_symbols])
    candidate_idx = candidate_idx[present]
    candidate_symbols = candidate_symbols[present]
    candidate_audit = candidate_audit.loc[present].reset_index(drop=True)
    tumor_all = tumor_all.reindex(candidate_symbols)
    monkey_candidates = values[candidate_idx, :]
    tumor_values = tumor_all.to_numpy(dtype=float)

    monkey_rank_mean = percentile_ranks(monkey_candidates).mean(axis=1)
    tumor_rank_mean = percentile_ranks(np.log1p(np.clip(tumor_values, 0, None))).mean(axis=1)
    rank_shift = np.abs(tumor_rank_mean - monkey_rank_mean)
    pc_loading = tumor_pc1_loadings(tumor_values)
    tumor_state_blacklist = set() if args.disable_hallmark_filter else load_tumor_state_blacklist(args.hallmark_file)
    hallmark_blacklisted = np.asarray([gene in tumor_state_blacklist for gene in candidate_symbols])
    rank_cut = float(np.quantile(rank_shift, args.rank_shift_quantile))
    pc_cut = float(np.quantile(pc_loading, args.pc_loading_quantile))
    keep = (rank_shift <= rank_cut) & (pc_loading <= pc_cut) & ~hallmark_blacklisted

    audit = candidate_audit.copy()
    audit["gene_symbol"] = candidate_symbols
    audit["monkey_mean_rank"] = monkey_rank_mean
    audit["tumor_mean_rank"] = tumor_rank_mean
    audit["absolute_rank_shift"] = rank_shift
    audit["tumor_pc1_abs_loading"] = pc_loading
    audit["hallmark_tumor_state_blacklist"] = hallmark_blacklisted
    audit["passes_domain_filter"] = keep
    audit = audit.sort_values(["passes_domain_filter", "fisher_score"], ascending=[False, False])
    passing = audit[audit["passes_domain_filter"]].copy()
    passing_positions = np.asarray([monkey.index.get_loc(gene) for gene in passing["gene_symbol"]], dtype=int)
    passing_values = values[passing_positions, :].astype(float)
    owner_rows: list[dict[str, Any]] = []
    per_network = max(1, args.final_genes // len(groups))
    selected_symbols_set: set[str] = set()
    for group in groups:
        in_group = labels == group
        out_group = ~in_group
        mean_in = passing_values[:, in_group].mean(axis=1)
        mean_out = passing_values[:, out_group].mean(axis=1)
        var_in = passing_values[:, in_group].var(axis=1)
        var_out = passing_values[:, out_group].var(axis=1)
        effect = np.abs(mean_in - mean_out) / np.sqrt(var_in + var_out + 1e-8)
        order = np.argsort(effect)[::-1][:per_network]
        for index in order:
            gene = str(passing.iloc[index]["gene_symbol"])
            selected_symbols_set.add(gene)
            owner_rows.append({"gene_symbol": gene, "balanced_owner_network": group, "one_vs_rest_effect": effect[index]})
    if len(selected_symbols_set) < args.final_genes:
        for gene in passing.sort_values("fisher_score", ascending=False)["gene_symbol"].astype(str):
            selected_symbols_set.add(gene)
            if len(selected_symbols_set) >= args.final_genes:
                break
    selected = passing[passing["gene_symbol"].astype(str).isin(selected_symbols_set)].copy()
    owner = pd.DataFrame(owner_rows)
    if len(owner):
        owner_summary = (
            owner.sort_values("one_vs_rest_effect", ascending=False)
            .groupby("gene_symbol", as_index=False)
            .agg(
                balanced_owner_network=("balanced_owner_network", lambda values: " | ".join(sorted(set(values)))),
                one_vs_rest_effect=("one_vs_rest_effect", "max"),
            )
        )
        selected = selected.merge(owner_summary, on="gene_symbol", how="left")
    selected = selected.sort_values("fisher_score", ascending=False).head(args.final_genes)
    selected_symbols = selected["gene_symbol"].astype(str).to_numpy()
    selected_idx = np.asarray([monkey.index.get_loc(gene) for gene in selected_symbols], dtype=int)
    selected.to_csv(args.outdir / "tumor_adapted_network_genes.csv", index=False, encoding="utf-8-sig")
    audit.to_csv(args.outdir / "candidate_gene_domain_filter_audit.csv", index=False, encoding="utf-8-sig")

    monkey_selected = values[selected_idx, :].astype(float)
    adapted_reference, _ = build_group_reference(monkey_selected, labels, groups, heldout_idx=-1)
    tumor_65 = tumor_all.reindex(selected_symbols)[sample_ids].to_numpy(dtype=float)
    raw_tpm_scores = corr_scores(adapted_reference, tumor_65)
    log_scores = corr_scores(adapted_reference, np.log1p(np.clip(tumor_65, 0, None)))
    harmonized = quantile_map_to_vsd(np.log1p(np.clip(tumor_65, 0, None)), monkey_selected)
    harmonized_scores = corr_scores(adapted_reference, harmonized)

    strict_labels = truth["whole_tumor_network_dominant"].astype(str).to_numpy()
    calibrated_scores, calibration = nested_bias_prior_calibration(harmonized_scores, strict_labels, groups)
    calibration["patient_barcode"] = truth["patient_barcode"]
    calibration.to_csv(args.outdir / "nested_calibration_choices.csv", index=False, encoding="utf-8-sig")

    monkey_scores = monkey_loo_scores(monkey_selected, labels, groups)
    monkey_max = monkey_scores.max(axis=0)
    monkey_order = np.argsort(monkey_scores, axis=0)[::-1]
    monkey_top1 = np.asarray([groups[int(monkey_order[0, i])] for i in range(monkey_order.shape[1])])
    monkey_top3_hit = np.asarray(
        [labels[i] in {groups[int(j)] for j in monkey_order[:3, i]} for i in range(monkey_order.shape[1])]
    )
    ood_threshold = float(np.quantile(monkey_max, args.ood_quantile))
    tumor_max = harmonized_scores.max(axis=0)
    accepted = tumor_max >= ood_threshold

    output = truth[
        [
            "patient_barcode",
            "sample_id",
            "whole_tumor_network_dominant",
            "whole_tumor_network_candidates",
            "network_prediction_top1",
            "network_prediction_top3",
        ]
    ].copy()
    output = output.rename(
        columns={
            "network_prediction_top1": "baseline_top1",
            "network_prediction_top3": "baseline_top3",
        }
    )
    for route, scores in (
        ("adapted_raw_tpm", raw_tpm_scores),
        ("adapted_log1p", log_scores),
        ("adapted_harmonized", harmonized_scores),
        ("adapted_harmonized_calibrated", calibrated_scores),
    ):
        top1, top3 = rank_predictions(scores, groups)
        output[f"{route}_top1"] = top1
        output[f"{route}_top3"] = top3
        output[f"{route}_max_score"] = scores.max(axis=0)
        sorted_scores = np.sort(scores, axis=0)
        output[f"{route}_margin"] = sorted_scores[-1] - sorted_scores[-2]
    output["ood_threshold"] = ood_threshold
    output["ood_accepted"] = accepted
    output.to_csv(args.outdir / "tcga_lgg_65_tumor_adapted_predictions.csv", index=False, encoding="utf-8-sig")

    metrics = [
        summarize_route(output, "baseline_production", "baseline_top1", "baseline_top3"),
        summarize_route(output, "adapted_raw_tpm", "adapted_raw_tpm_top1", "adapted_raw_tpm_top3"),
        summarize_route(output, "adapted_log1p", "adapted_log1p_top1", "adapted_log1p_top3"),
        summarize_route(
            output,
            "adapted_harmonized",
            "adapted_harmonized_top1",
            "adapted_harmonized_top3",
        ),
        summarize_route(
            output,
            "adapted_harmonized_calibrated",
            "adapted_harmonized_calibrated_top1",
            "adapted_harmonized_calibrated_top3",
        ),
        summarize_route(
            output,
            "adapted_harmonized_calibrated_ood",
            "adapted_harmonized_calibrated_top1",
            "adapted_harmonized_calibrated_top3",
            "ood_accepted",
        ),
    ]
    metrics_df = pd.DataFrame(metrics)
    metrics_df.to_csv(args.outdir / "route_ablation_metrics.csv", index=False, encoding="utf-8-sig")

    args.model_out.parent.mkdir(parents=True, exist_ok=True)
    target_quantiles = np.median(np.sort(monkey_selected, axis=0), axis=1)
    np.savez_compressed(
        args.model_out,
        genes=selected_symbols.astype(str),
        networks=np.asarray(groups, dtype=str),
        reference=adapted_reference.astype(np.float32),
        target_quantiles=target_quantiles.astype(np.float32),
        calibration_offsets=np.zeros(len(groups), dtype=np.float32),
        ood_max_correlation_threshold=np.asarray([ood_threshold], dtype=np.float32),
    )
    model_metadata = {
        "created_at": "2026-06-12",
        "purpose": "tumor-domain exploratory Network tracing with strict OOD rejection",
        "input_scale": "TPM; inference applies log1p and within-sample quantile mapping to monkey VSD target",
        "n_genes": int(len(selected_symbols)),
        "n_networks": int(len(groups)),
        "gene_filters": {
            "ensembl_one2one_high_confidence_human_macaque_ortholog": True,
            "valid_human_symbol_and_present_in_tcga_lgg": True,
            "rank_shift_quantile": args.rank_shift_quantile,
            "tumor_pc1_loading_quantile": args.pc_loading_quantile,
            "hallmark_tumor_state_filter": not args.disable_hallmark_filter,
            "balanced_markers_per_network": max(1, args.final_genes // len(groups)),
        },
        "class_prior_calibration": {
            "implemented": True,
            "nested_gate_result": "prior_strength_zero_in_all_outer_folds",
            "score_bias_result": "exploratory_fold_specific_bias_not_deployed",
            "offsets": [0.0] * len(groups),
        },
        "ood_rejection": {
            "metric": "maximum Pearson correlation",
            "calibration": "monkey leave-one-out",
            "quantile": args.ood_quantile,
            "threshold": ood_threshold,
        },
        "source_domain_fixed_panel_diagnostic": {
            "n": int(len(labels)),
            "top1_accuracy": float(np.mean(monkey_top1 == labels)),
            "top3_accuracy": float(monkey_top3_hit.mean()),
            "note": "Fixed-panel LOO diagnostic; gene selection used full source labels and is not formal leakage-free LOSO.",
        },
    }
    args.model_out.with_suffix(".json").write_text(
        json.dumps(model_metadata, ensure_ascii=False, indent=2), encoding="utf-8"
    )

    summary = {
        "n_patients": int(len(output)),
        "n_candidate_valid_human_genes": int(len(candidate_symbols)),
        "n_domain_filter_pass": int(keep.sum()),
        "n_hallmark_blacklisted_candidates": int(hallmark_blacklisted.sum()),
        "n_final_genes": int(len(selected_symbols)),
        "hallmark_filter_enabled": not args.disable_hallmark_filter,
        "rank_shift_cutoff": rank_cut,
        "pc_loading_cutoff": pc_cut,
        "ood_calibration": {
            "source": "monkey leave-one-out maximum correlation",
            "quantile": args.ood_quantile,
            "threshold": ood_threshold,
            "n_accepted": int(accepted.sum()),
            "coverage": float(accepted.mean()),
        },
        "source_domain_fixed_panel_diagnostic": model_metadata["source_domain_fixed_panel_diagnostic"],
        "model_output": str(args.model_out),
        "top1_distribution": {
            route: dict(Counter(output[f"{route}_top1"]))
            for route in (
                "adapted_raw_tpm",
                "adapted_log1p",
                "adapted_harmonized",
                "adapted_harmonized_calibrated",
            )
        },
        "metrics": metrics,
    }
    (args.outdir / "tumor_adapted_evaluation_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
