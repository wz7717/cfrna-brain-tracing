#!/usr/bin/env python
from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from benchmark.metrics import compute_multiclass_auc  # noqa: E402
from benchmark_runner import (  # noqa: E402
    _compute_confusion_tables,
    _compute_roc_tables,
    _compute_separability_tables,
    _compute_stability_tables,
)
from reporting.benchmark_figure_export import (  # noqa: E402
    export_benchmark_paper_figures,
)


DEFAULT_MATRIX = ROOT / "bo2023 data" / "mfas5_819samples_23605genes_vsd4_rmbatch.xls"
DEFAULT_SAMPLE_INFO = ROOT / "bo2023 data" / "Information of sequenced samples_update_full878_filter819.xlsx"
DEFAULT_OUTDIR = ROOT / "results" / "bo2023_loso_30_vsd_correlation"


def read_vsd_matrix(path: Path) -> pd.DataFrame:
    matrix = pd.read_csv(path, sep="\t", index_col=0)
    matrix.index = matrix.index.astype(str).str.strip()
    matrix.columns = matrix.columns.astype(str).str.strip()
    matrix = matrix[~matrix.index.duplicated(keep="first")]
    matrix = matrix.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    finite = np.isfinite(matrix.to_numpy(dtype=float)).all(axis=1)
    matrix = matrix.loc[finite]
    return matrix


def read_annotations(path: Path, sheet: str, region_col: str) -> pd.DataFrame:
    ann = pd.read_excel(path, sheet_name=sheet)
    required = {"No.", region_col}
    missing = required - set(ann.columns)
    if missing:
        raise ValueError(f"annotation is missing required columns: {sorted(missing)}")
    ann = ann.copy()
    ann["sample_id"] = ann["No."].astype(str).str.strip()
    ann["region_id"] = ann[region_col].astype(str).str.strip()
    ann = ann.drop_duplicates("sample_id")
    return ann[ann["sample_id"].ne("") & ann["region_id"].ne("")]


def softmax(scores: np.ndarray) -> np.ndarray:
    values = np.asarray(scores, dtype=float)
    values = values - np.nanmax(values)
    exp_values = np.exp(values)
    return exp_values / (exp_values.sum() + 1e-12)


def correlation_scores(reference: np.ndarray, sample: np.ndarray, rows: np.ndarray | None = None) -> np.ndarray:
    ref = reference if rows is None else reference[rows, :]
    vec = sample if rows is None else sample[rows]
    ref0 = ref - ref.mean(axis=0, keepdims=True)
    vec0 = vec - vec.mean()
    denom = np.sqrt(np.square(ref0).sum(axis=0) * np.square(vec0).sum() + 1e-12)
    scores = (ref0 * vec0[:, None]).sum(axis=0) / denom
    return np.nan_to_num(scores, nan=0.0, posinf=0.0, neginf=0.0)


def bootstrap_top1_stability(
    reference: np.ndarray,
    sample: np.ndarray,
    predicted_region: str,
    regions: list[str],
    n_bootstrap: int,
    gene_fraction: float,
    seed: int,
) -> float:
    if n_bootstrap <= 0:
        return float("nan")
    rng = np.random.default_rng(seed)
    n_genes = reference.shape[0]
    subset_size = max(2, int(round(n_genes * gene_fraction)))
    hits = 0
    for _ in range(n_bootstrap):
        rows = rng.choice(n_genes, size=subset_size, replace=False)
        top_region = regions[int(np.argmax(correlation_scores(reference, sample, rows)))]
        hits += int(top_region == predicted_region)
    return float(hits / n_bootstrap)


def select_test_samples(
    ann: pd.DataFrame,
    matrix_samples: list[str],
    n_samples: int,
    seed: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    aligned = ann[ann["sample_id"].isin(set(matrix_samples))].copy()
    region_counts = aligned.groupby("region_id")["sample_id"].size()
    singleton_regions = set(region_counts[region_counts < 2].index)
    excluded = aligned[aligned["region_id"].isin(singleton_regions)].copy()
    eligible = aligned[~aligned["region_id"].isin(singleton_regions)].copy()
    if n_samples > len(eligible):
        raise ValueError(f"requested {n_samples} samples but only {len(eligible)} are LOSO-assessable")
    selected = eligible.sample(n=n_samples, random_state=seed).sort_values("sample_id").reset_index(drop=True)
    return selected, excluded


def build_region_reference(
    values: np.ndarray,
    sample_ids: list[str],
    ann: pd.DataFrame,
) -> tuple[np.ndarray, list[str], np.ndarray, dict[str, np.ndarray]]:
    sample_pos = {sample_id: idx for idx, sample_id in enumerate(sample_ids)}
    regions = sorted(ann["region_id"].unique().tolist())
    indices: dict[str, np.ndarray] = {}
    centroids = np.zeros((values.shape[0], len(regions)), dtype=np.float32)
    counts = np.zeros(len(regions), dtype=int)
    for j, region in enumerate(regions):
        ids = ann.loc[ann["region_id"] == region, "sample_id"].tolist()
        idx = np.asarray([sample_pos[x] for x in ids if x in sample_pos], dtype=int)
        if idx.size == 0:
            raise ValueError(f"region {region} has no aligned expression sample")
        indices[region] = idx
        counts[j] = int(idx.size)
        centroids[:, j] = values[:, idx].mean(axis=1, dtype=np.float64).astype(np.float32)
    return centroids, regions, counts, indices


def make_suite(
    detail_df: pd.DataFrame,
    probability_df: pd.DataFrame,
    evaluated_regions: list[str],
    total_classes: int,
) -> dict[str, pd.DataFrame]:
    valid = detail_df[detail_df["abstained"] == 0].copy()
    top1 = float(valid["hit1"].mean()) if len(valid) else float("nan")
    top3 = float(valid["hit3"].mean()) if len(valid) else float("nan")
    balanced = float(valid.groupby("label")["hit1"].mean().mean()) if len(valid) else float("nan")
    auc = compute_multiclass_auc(
        probability_df["label"].astype(str).tolist(),
        probability_df[[c for c in probability_df.columns if c not in {"sample_id", "label"}]],
    )
    summary = {
        "n_total": int(len(detail_df)),
        "n_valid": int(len(valid)),
        "n_classes": int(total_classes),
        "abstain_rate": 0.0,
        "top1_acc": top1,
        "top3_acc": top3,
        "balanced_acc": balanced,
        "auc": auc,
        "mean_top1_confidence": float(valid["top1_confidence"].mean()),
        "mean_decision_margin": float(valid["decision_margin"].mean()),
        "mean_top1_stability": float(valid["top1_stability"].mean()),
    }
    metrics_df = pd.DataFrame(
        [
            {"metric": "N_samples_total", "value": len(detail_df)},
            {"metric": "N_samples_valid", "value": len(valid)},
            {"metric": "N_classes", "value": total_classes},
            {"metric": "Abstain_rate", "value": 0.0},
            {"metric": "Top1_acc_valid", "value": top1},
            {"metric": "Top3_acc_valid", "value": top3},
            {"metric": "Balanced_acc_valid", "value": balanced},
            {"metric": "MacroAUC_ovr_valid", "value": auc},
            {"metric": "Mean_top1_confidence_valid", "value": summary["mean_top1_confidence"]},
            {"metric": "Mean_decision_margin_valid", "value": summary["mean_decision_margin"]},
            {"metric": "Mean_top1_stability_valid", "value": summary["mean_top1_stability"]},
        ]
    )
    confusion_raw, confusion_norm, confusion_long = _compute_confusion_tables(valid)
    roc_curve, roc_summary = _compute_roc_tables(probability_df, evaluated_regions)
    separability, centroids, centroid_dist = _compute_separability_tables(
        probability_df, detail_df, evaluated_regions
    )
    stability_region, stability_bin = _compute_stability_tables(detail_df)
    detail_df.attrs["summary"] = summary
    return {
        "detail_df": detail_df,
        "metrics_df": metrics_df,
        "probability_df": probability_df,
        "confusion_raw_df": confusion_raw,
        "confusion_norm_df": confusion_norm,
        "confusion_long_df": confusion_long,
        "roc_curve_df": roc_curve,
        "roc_summary_df": roc_summary,
        "separability_df": separability,
        "centroid_df": centroids,
        "centroid_distance_df": centroid_dist,
        "stability_region_df": stability_region,
        "stability_bin_df": stability_bin,
    }


def write_chinese_report(
    outdir: Path,
    summary: dict[str, Any],
    selected: pd.DataFrame,
    excluded: pd.DataFrame,
    detail: pd.DataFrame,
    top_confusions: list[dict[str, Any]],
) -> None:
    top1_n = int(detail["hit1"].sum())
    top3_n = int(detail["hit3"].sum())
    med_rank = float(detail["true_rank"].median())
    correct = detail[detail["hit1"] == 1]
    incorrect = detail[detail["hit1"] == 0]
    correct_stability = float(correct["top1_stability"].mean()) if len(correct) else float("nan")
    incorrect_stability = float(incorrect["top1_stability"].mean()) if len(incorrect) else float("nan")
    correct_margin = float(correct["decision_margin"].mean()) if len(correct) else float("nan")
    incorrect_margin = float(incorrect["decision_margin"].mean()) if len(incorrect) else float("nan")
    worst = detail.sort_values("true_rank", ascending=False).head(5)
    confusion_lines = [
        f"- `{x['truth_region']}` -> `{x['pred_region']}`: {int(x['count'])} 次"
        for x in top_confusions
    ] or ["- 本轮无错误混淆。"]
    worst_lines = [
        f"- `{r.sample_id}`: 真实 `{r.label}`，Top1 `{r.pred_top1}`，真实排名 {int(r.true_rank)}"
        for r in worst.itertuples(index=False)
    ]
    content = f"""# Bo2023 30 样本严格 LOSO 脑区回溯验证

## 验证设计

- 数据：Bo2023 VSD + batch-removed 表达矩阵，819 个已知来源样本、{summary['n_genes']} 个基因、{summary['n_regions']} 个 Region。
- 本轮：从可严格留一评估的样本中固定随机种子 `{summary['seed']}` 抽取 30 个样本。
- 不可评估样本：{len(excluded)} 个 singleton Region 样本，留出后训练集中不再存在其真实 Region，因此未进入本轮随机池。
- 每折：仅以另外 818 个样本构建 `gene x region` reference，使用留出样本与每个 region centroid 的 Pearson 相关性排序。
- 表达口径：输入为 VSD normalized-expression pattern；本轮分数代表模式相似度，不代表 TPM 或组织贡献比例。
- Bootstrap：每折对基因抽样 {summary['bootstrap_n']} 次，每次使用 {summary['bootstrap_gene_fraction']:.0%} 基因，评估 Top1 稳定性。

## 核心结果

| 指标 | 结果 |
| --- | ---: |
| Top1 命中 | {top1_n}/30 ({summary['top1_acc']:.1%}) |
| Top3 命中 | {top3_n}/30 ({summary['top3_acc']:.1%}) |
| Macro AUC (OVR) | {summary['macro_auc']:.3f} |
| 真实 Region 中位排名 | {med_rank:.1f} |
| 平均 Top1 bootstrap 稳定性 | {summary['mean_stability']:.3f} |
| 正确 / 错误 Top1 平均稳定性 | {correct_stability:.3f} / {incorrect_stability:.3f} |
| 正确 / 错误 Top1 平均 margin | {correct_margin:.6f} / {incorrect_margin:.6f} |

## 直观解释

本轮 Top1 命中率为 **{summary['top1_acc']:.1%}**，Top3 命中率为 **{summary['top3_acc']:.1%}**。
在每折约 {summary['n_regions']} 个候选 Region 的条件下，随机期望约为 Top1 `{1/summary['n_regions']:.1%}`、Top3 `{3/summary['n_regions']:.1%}`。
因此本结果用于判断 VSD centroid-correlation 是否具备初步脑区筛选能力；30 样本仍是 pilot，不能替代完整 819 样本评估。

需要特别注意：错误 Top1 的平均 bootstrap 稳定性仍达到 **{incorrect_stability:.3f}**，表示部分误判是稳定的相似脑区混淆，而不是随机不稳定造成的。当前 softmax confidence 来自相关系数的相对排序，尚未做概率校准，不应直接当成正确概率解释。

## 主要混淆

{chr(10).join(confusion_lines)}

## 排名最靠后的样本

{chr(10).join(worst_lines)}

## 优化方向

1. 运行全量可评估 LOSO，并将 singleton Region 另列为“无训练参考不可评估”，得到稳定的 region-wise 准确率。
2. 在每一折内部重新筛选 region-specific marker genes，再比较全基因相关性与 marker-weighted 相关性；marker 选择不能使用留出样本。
3. 按相邻/同网络脑区汇总混淆，评估分层分类策略：先定位 lobe/network，再在局部 Region 内细分。
4. 对低 margin 或低 bootstrap stability 的预测设置不确定/人工复核阈值，而不是强行输出 Top1。

## 输出图表

`figures/` 中包含数据库 Benchmark 页面同类型的 Figure1-Figure6：总览指标、混淆矩阵、ROC/AUC、可分离性与稳定性、region 可靠性、failure mode。
"""
    (outdir / "validation_report_cn.md").write_text(content, encoding="utf-8")
    selected.to_csv(outdir / "selected_test_samples.csv", index=False, encoding="utf-8-sig")


def export_loso_dashboard(
    outdir: Path,
    detail: pd.DataFrame,
    n_regions: int,
    title: str = "Bo2023 strict LOSO pilot: 30 held-out samples, VSD centroid correlation",
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    top1 = float(detail["hit1"].mean())
    top3 = float(detail["hit3"].mean())
    chance = [1.0 / n_regions, 3.0 / n_regions]
    fig = plt.figure(figsize=(15, 9))
    grid = fig.add_gridspec(2, 2, width_ratios=[0.9, 1.2], height_ratios=[1.0, 1.0])
    ax1 = fig.add_subplot(grid[0, 0])
    ax2 = fig.add_subplot(grid[0, 1])
    ax3 = fig.add_subplot(grid[1, :])

    x = np.arange(2)
    actual = [top1, top3]
    ax1.bar(x - 0.18, actual, width=0.36, label="LOSO observed", color="#2364aa")
    ax1.bar(x + 0.18, chance, width=0.36, label="Random expectation", color="#bcccdc")
    ax1.set_xticks(x, ["Top1", "Top3"])
    ax1.set_ylim(0, max(0.55, top3 + 0.1))
    ax1.set_ylabel("Hit rate")
    ax1.set_title("A. Accuracy versus chance")
    for i, value in enumerate(actual):
        ax1.text(i - 0.18, value + 0.015, f"{value:.1%}", ha="center")
    ax1.legend(fontsize=9)

    bins = [0.5, 1.5, 3.5, 10.5, 25.5, 50.5, n_regions + 0.5]
    labels = ["1", "2-3", "4-10", "11-25", "26-50", f"51-{n_regions}"]
    rank_group = pd.cut(detail["true_rank"], bins=bins, labels=labels, include_lowest=True)
    rank_counts = rank_group.value_counts(sort=False)
    ax2.bar(labels, rank_counts.values, color="#4c956c")
    ax2.set_title("B. True-region rank distribution")
    ax2.set_ylabel("Samples")
    ax2.tick_params(axis="x", rotation=20)
    for i, value in enumerate(rank_counts.values):
        ax2.text(i, value + 0.15, str(int(value)), ha="center")

    ordered = detail.sort_values("true_rank", ascending=True).reset_index(drop=True)
    colors = np.where(ordered["hit1"].eq(1), "#138a72", np.where(ordered["hit3"].eq(1), "#f4a261", "#d1495b"))
    ax3.scatter(np.arange(len(ordered)), ordered["true_rank"], c=colors, s=55)
    ax3.axhspan(0.5, 3.5, color="#e8f4ea", alpha=0.7)
    ax3.axhline(3.5, color="#4c956c", linestyle="--", linewidth=1)
    ax3.set_xticks(np.arange(len(ordered)), ordered["sample_id"], rotation=65, ha="right", fontsize=8)
    ax3.set_ylim(n_regions + 2, 0)
    ax3.set_ylabel("True-region rank (lower is better)")
    ax3.set_title("C. Per-sample true-region rank; shaded area denotes Top3")
    fig.suptitle(title, fontsize=15)
    fig.tight_layout()
    fig.savefig(outdir / "bo2023_loso30_summary_dashboard.png", dpi=250, bbox_inches="tight")
    fig.savefig(outdir / "bo2023_loso30_summary_dashboard.pdf", bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser(description="Strict Bo2023 leave-one-sample-out brain-region validation.")
    parser.add_argument("--matrix", type=Path, default=DEFAULT_MATRIX)
    parser.add_argument("--sample-info", type=Path, default=DEFAULT_SAMPLE_INFO)
    parser.add_argument("--sample-sheet", default="mfas5_819samples_phenSet4")
    parser.add_argument("--region-col", default="Region")
    parser.add_argument("--n-samples", type=int, default=30)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--bootstrap-n", type=int, default=25)
    parser.add_argument("--bootstrap-gene-fraction", type=float, default=0.7)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument("--save-references", action="store_true")
    args = parser.parse_args()

    outdir = args.outdir
    outdir.mkdir(parents=True, exist_ok=True)
    references_dir = outdir / "fold_references"
    if args.save_references:
        references_dir.mkdir(parents=True, exist_ok=True)

    matrix = read_vsd_matrix(args.matrix)
    ann = read_annotations(args.sample_info, args.sample_sheet, args.region_col)
    matrix_samples = matrix.columns.tolist()
    ann = ann[ann["sample_id"].isin(set(matrix_samples))].copy()
    if len(ann) != matrix.shape[1]:
        raise ValueError(f"expression/annotation alignment failed: matrix={matrix.shape[1]}, annotations={len(ann)}")
    selected, excluded = select_test_samples(ann, matrix_samples, args.n_samples, args.seed)

    values = matrix.to_numpy(dtype=np.float32)
    full_reference, regions, region_counts, region_indices = build_region_reference(values, matrix_samples, ann)
    region_pos = {region: idx for idx, region in enumerate(regions)}
    sample_pos = {sample_id: idx for idx, sample_id in enumerate(matrix_samples)}
    fold_rows: list[dict[str, Any]] = []
    probability_rows: list[dict[str, Any]] = []
    score_rows: list[dict[str, Any]] = []
    fold_manifest: list[dict[str, Any]] = []

    for fold_no, row in enumerate(selected.itertuples(index=False), start=1):
        sample_id = str(row.sample_id)
        truth = str(row.region_id)
        test_idx = sample_pos[sample_id]
        region_idx = region_pos[truth]
        n_before = int(region_counts[region_idx])
        if n_before < 2:
            raise ValueError(f"held-out sample {sample_id} is in singleton region {truth}")
        reference = full_reference.copy()
        source_idx = region_indices[truth]
        train_idx = source_idx[source_idx != test_idx]
        reference[:, region_idx] = values[:, train_idx].mean(axis=1, dtype=np.float64).astype(np.float32)
        sample = values[:, test_idx]
        scores = correlation_scores(reference, sample)
        probs = softmax(scores)
        order = np.argsort(scores)[::-1]
        pred_regions = [regions[int(i)] for i in order]
        true_rank = pred_regions.index(truth) + 1
        stability = bootstrap_top1_stability(
            reference=reference,
            sample=sample,
            predicted_region=pred_regions[0],
            regions=regions,
            n_bootstrap=args.bootstrap_n,
            gene_fraction=args.bootstrap_gene_fraction,
            seed=args.seed + fold_no,
        )
        fold_rows.append(
            {
                "fold": fold_no,
                "sample_id": sample_id,
                "label": truth,
                "pred_top1": pred_regions[0],
                "pred_top2": pred_regions[1],
                "pred_top3": pred_regions[2],
                "true_rank": true_rank,
                "hit1": int(true_rank == 1),
                "hit3": int(true_rank <= 3),
                "top1_score": float(scores[order[0]]),
                "true_region_score": float(scores[region_idx]),
                "top1_confidence": float(probs[order[0]]),
                "top2_confidence": float(probs[order[1]]),
                "decision_margin": float(scores[order[0]] - scores[order[1]]),
                "top1_stability": stability,
                "abstained": 0,
                "traceability": "strict_loso_vsd_correlation",
                "overlap_genes": int(values.shape[0]),
                "truth_region_train_samples": n_before - 1,
            }
        )
        probability = {"sample_id": sample_id, "label": truth}
        probability.update({region: float(probs[j]) for j, region in enumerate(regions)})
        probability_rows.append(probability)
        for rank, j in enumerate(order, start=1):
            score_rows.append(
                {
                    "fold": fold_no,
                    "sample_id": sample_id,
                    "truth_region": truth,
                    "candidate_region": regions[int(j)],
                    "rank": rank,
                    "correlation_score": float(scores[int(j)]),
                    "softmax_confidence": float(probs[int(j)]),
                }
            )
        reference_file = None
        if args.save_references:
            reference_file = references_dir / f"fold_{fold_no:02d}_{sample_id}_reference.npz"
            np.savez_compressed(
                reference_file,
                gene_ids=matrix.index.to_numpy(dtype=str),
                regions=np.asarray(regions, dtype=str),
                reference_matrix=reference,
                training_counts=region_counts - (np.asarray(regions) == truth).astype(int),
                held_out_sample=np.asarray([sample_id], dtype=str),
                truth_region=np.asarray([truth], dtype=str),
            )
        fold_manifest.append(
            {
                "fold": fold_no,
                "held_out_sample": sample_id,
                "truth_region": truth,
                "truth_region_train_samples": n_before - 1,
                "reference_shape": [int(reference.shape[0]), int(reference.shape[1])],
                "reference_file": str(reference_file.relative_to(outdir)) if reference_file else None,
            }
        )

    detail_df = pd.DataFrame(fold_rows)
    probability_df = pd.DataFrame(probability_rows)
    all_scores_df = pd.DataFrame(score_rows)
    evaluated_regions = sorted(detail_df["label"].unique().tolist())
    suite = make_suite(detail_df, probability_df, evaluated_regions, len(regions))
    metrics = {r.metric: float(r.value) for r in suite["metrics_df"].itertuples(index=False)}
    confusion_errors = suite["confusion_long_df"]
    confusion_errors = confusion_errors[
        confusion_errors["truth_region"] != confusion_errors["pred_region"]
    ].sort_values(["count", "row_fraction"], ascending=False)
    top_confusions = confusion_errors.head(5).to_dict(orient="records")

    metadata = {
        "validation_design": "strict leave-one-sample-out",
        "dataset": "Bo2023_WangLab_VSD_region",
        "expression_space": "VSD_batch_removed",
        "method": "correlation",
        "method_description": "Pearson correlation to fold-specific region centroids in raw VSD space",
        "n_sample_total": int(matrix.shape[1]),
        "n_test_samples": int(len(detail_df)),
        "n_regions": int(len(regions)),
        "n_genes": int(matrix.shape[0]),
        "selection_seed": int(args.seed),
        "bootstrap_n": int(args.bootstrap_n),
        "bootstrap_gene_fraction": float(args.bootstrap_gene_fraction),
        "singleton_samples_excluded_from_selection_pool": int(len(excluded)),
    }
    export_benchmark_paper_figures(outdir, suite, metadata=metadata, prefix="bo2023_loso30")
    all_scores_df.to_csv(outdir / "all_candidate_scores.csv", index=False, encoding="utf-8-sig")
    excluded.to_csv(outdir / "singleton_samples_not_loso_assessable.csv", index=False, encoding="utf-8-sig")
    (outdir / "fold_manifest.json").write_text(
        json.dumps(fold_manifest, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    final_summary = {
        "created_at": datetime.now(timezone.utc).isoformat(),
        **metadata,
        "seed": int(args.seed),
        "top1_acc": metrics["Top1_acc_valid"],
        "top3_acc": metrics["Top3_acc_valid"],
        "macro_auc": metrics["MacroAUC_ovr_valid"],
        "mean_stability": metrics["Mean_top1_stability_valid"],
        "top1_hits": int(detail_df["hit1"].sum()),
        "top3_hits": int(detail_df["hit3"].sum()),
    }
    export_loso_dashboard(outdir, detail_df, len(regions))
    (outdir / "validation_summary.json").write_text(
        json.dumps(final_summary, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    write_chinese_report(outdir, final_summary, selected, excluded, detail_df, top_confusions)
    print(json.dumps(final_summary, ensure_ascii=False, indent=2))
    print(f"Outputs written to: {outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
