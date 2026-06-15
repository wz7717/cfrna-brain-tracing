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
from scripts.run_bo2023_loso_validation import correlation_scores, read_vsd_matrix, softmax  # noqa: E402
from scripts.run_bo2023_v2_loso_validation import DEFAULT_GENE_MAP, DEFAULT_MATRIX, DEFAULT_SAMPLE_INFO, map_matrix_to_symbols  # noqa: E402


DEFAULT_OUTDIR = ROOT / "results" / "bo2023_network_discriminative_correlation_unseen_confirmation"


def choose_unseen_samples(
    annotations: pd.DataFrame,
    matrix_samples: list[str],
    excluded_ids: set[str],
    n_samples: int,
    seed: int,
) -> pd.DataFrame:
    eligible = annotations[
        annotations["sample_id"].isin(set(matrix_samples))
        & ~annotations["sample_id"].isin(excluded_ids)
    ].copy()
    return eligible.sample(n=n_samples, random_state=seed).sort_values("sample_id").reset_index(drop=True)


def build_group_reference(
    values: np.ndarray,
    labels: np.ndarray,
    groups: list[str],
    heldout_idx: int,
) -> tuple[np.ndarray, dict[str, np.ndarray]]:
    columns: list[np.ndarray] = []
    training: dict[str, np.ndarray] = {}
    for group in groups:
        idx = np.flatnonzero(labels == group)
        idx = idx[idx != heldout_idx]
        if not len(idx):
            raise ValueError(f"group {group} has no training sample after holdout")
        training[group] = idx
        columns.append(values[:, idx].mean(axis=1, dtype=np.float64))
    return np.column_stack(columns), training


def select_group_discriminative_genes(
    values: np.ndarray,
    groups: list[str],
    training_indices: dict[str, np.ndarray],
    top_n: int,
) -> tuple[np.ndarray, pd.DataFrame]:
    centroids: list[np.ndarray] = []
    counts: list[int] = []
    within_sum = np.zeros(values.shape[0], dtype=float)
    within_denom = 0
    for group in groups:
        x = values[:, training_indices[group]].astype(float, copy=False)
        mean = x.mean(axis=1)
        centroids.append(mean)
        counts.append(int(x.shape[1]))
        within_sum += np.square(x - mean[:, None]).sum(axis=1)
        within_denom += x.shape[1] - 1
    means = np.column_stack(centroids)
    group_weights = np.asarray(counts, dtype=float)
    overall = np.average(means, axis=1, weights=group_weights)
    between = np.average(np.square(means - overall[:, None]), axis=1, weights=group_weights)
    within = within_sum / max(float(within_denom), 1.0)
    fisher = between / (within + 1e-8)
    rows = np.argsort(fisher)[::-1][: min(top_n, len(fisher))]
    audit = pd.DataFrame(
        {
            "gene_index": rows.astype(int),
            "fisher_score": fisher[rows],
            "between_variance": between[rows],
            "within_variance": within[rows],
        }
    )
    return rows.astype(int), audit


def evaluate_route(
    sample_id: str,
    truth: str,
    sample: np.ndarray,
    reference: np.ndarray,
    groups: list[str],
    genes: np.ndarray | None,
    route: str,
) -> tuple[dict[str, Any], dict[str, Any]]:
    scores = correlation_scores(reference, sample, genes)
    order = np.argsort(scores)[::-1]
    ranked = [groups[int(j)] for j in order]
    rank = ranked.index(truth) + 1
    probs = softmax(scores)
    detail = {
        "route": route,
        "sample_id": sample_id,
        "label": truth,
        "pred_top1": ranked[0],
        "pred_top2": ranked[1],
        "pred_top3": ranked[2],
        "true_rank": rank,
        "hit1": int(rank == 1),
        "hit3": int(rank <= 3),
    }
    probability: dict[str, Any] = {"sample_id": sample_id, "label": truth}
    probability.update({group: float(probs[j]) for j, group in enumerate(groups)})
    return detail, probability


def summarize(detail: pd.DataFrame, probability: pd.DataFrame, groups: list[str]) -> dict[str, Any]:
    return {
        "n": int(len(detail)),
        "top1_hits": int(detail["hit1"].sum()),
        "top1_accuracy": float(detail["hit1"].mean()),
        "top3_hits": int(detail["hit3"].sum()),
        "top3_accuracy": float(detail["hit3"].mean()),
        "macro_auc": float(compute_multiclass_auc(detail["label"].astype(str).tolist(), probability[groups])),
        "median_true_rank": float(detail["true_rank"].median()),
    }


def export_plot(outdir: Path, metrics: pd.DataFrame, top1_target: float, top3_target: float) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    base = metrics.loc[metrics["route"] == "network_correlation_all_genes"].iloc[0]
    enhanced = metrics.loc[metrics["route"] == "network_discriminative_correlation_top200"].iloc[0]
    vals = np.asarray(
        [[base.top1_accuracy, base.top3_accuracy, base.macro_auc], [enhanced.top1_accuracy, enhanced.top3_accuracy, enhanced.macro_auc]]
    )
    fig, axes = plt.subplots(1, 3, figsize=(14.0, 4.6), constrained_layout=True)
    for j, title in enumerate(["Top1 accuracy", "Top3 accuracy", "Macro AUC"]):
        bars = axes[j].bar(["Network correlation", "Discriminative network"], vals[:, j], color=["#0072B2", "#009E73"])
        axes[j].set_title(title, fontweight="bold")
        axes[j].grid(axis="y", alpha=0.25)
        axes[j].set_axisbelow(True)
        axes[j].tick_params(axis="x", rotation=15)
        if j < 2:
            target = [top1_target, top3_target][j]
            axes[j].axhline(target, color="#CC3311", linestyle="--", label="Requested target")
            axes[j].legend(fontsize=8)
        axes[j].set_ylim(0, max(1.0 if j < 2 else 1.0, float(vals[:, j].max()) + 0.1))
        for bar, val in zip(bars, vals[:, j]):
            axes[j].text(
                bar.get_x() + bar.get_width() / 2,
                val + 0.025,
                f"{val:.1%}" if j < 2 else f"{val:.3f}",
                ha="center",
                fontweight="bold",
            )
    fig.suptitle("Bo2023 SaleemNetworks source tracing: unseen confirmation", fontsize=14, fontweight="bold")
    fig.savefig(outdir / "network_correlation_target_confirmation.png", dpi=220, bbox_inches="tight")
    fig.savefig(outdir / "network_correlation_target_confirmation.pdf", bbox_inches="tight")
    plt.close(fig)


def write_report(outdir: Path, summary: dict[str, Any]) -> None:
    base = summary["routes"]["network_correlation_all_genes"]
    enhanced = summary["routes"]["network_discriminative_correlation_top200"]
    text = f"""# 网络层级 Correlation 主路径目标确认

## 终点定义

- 精确 `Region` 终点上一轮未达到 Top1/Top3 `50%/70%`。
- 本轮预先固定上层来源终点为样本表中的 `SaleemNetworks`（`{summary['n_classes']}` 类），仍使用 Pearson correlation 溯源。
- 判别路径在每折训练样本内选择 Top `{summary['top_n_genes']}` 个 Region 间/内方差比较高的基因；未使用测试样本调参。
- 验证范围：{summary['validation_scope']}。

## 结果

| 路径 | Top1 | Top3 | Macro AUC |
| --- | ---: | ---: | ---: |
| Network correlation 全基因 | {base['top1_hits']}/{summary['n_test_samples']} ({base['top1_accuracy']:.1%}) | {base['top3_hits']}/{summary['n_test_samples']} ({base['top3_accuracy']:.1%}) | {base['macro_auc']:.3f} |
| Network 判别基因 correlation | {enhanced['top1_hits']}/{summary['n_test_samples']} ({enhanced['top1_accuracy']:.1%}) | {enhanced['top3_hits']}/{summary['n_test_samples']} ({enhanced['top3_accuracy']:.1%}) | {enhanced['macro_auc']:.3f} |
| 目标门槛 | 50.0% | 70.0% | - |

## 判定

{summary['decision']}

## 应用边界

{summary['boundary']}
"""
    (outdir / "network_correlation_target_report_cn.md").write_text(text, encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description="Strict validation of network-level discriminative Pearson correlation.")
    parser.add_argument("--matrix", type=Path, default=DEFAULT_MATRIX)
    parser.add_argument("--sample-info", type=Path, default=DEFAULT_SAMPLE_INFO)
    parser.add_argument("--gene-map", type=Path, default=DEFAULT_GENE_MAP)
    parser.add_argument("--sample-sheet", default="mfas5_819samples_phenSet4")
    parser.add_argument("--endpoint", default="SaleemNetworks")
    parser.add_argument("--n-samples", type=int, default=30)
    parser.add_argument("--seed", type=int, default=20260531)
    parser.add_argument("--top-n-genes", type=int, default=200)
    parser.add_argument("--target-top1", type=float, default=0.50)
    parser.add_argument("--target-top3", type=float, default=0.70)
    parser.add_argument("--exclude-samples", type=Path, action="append", default=[])
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    (args.outdir / "fold_selected_genes").mkdir(parents=True, exist_ok=True)
    raw = read_vsd_matrix(args.matrix)
    matrix = map_matrix_to_symbols(raw, args.gene_map)
    ann = pd.read_excel(args.sample_info, sheet_name=args.sample_sheet)
    ann["sample_id"] = ann["No."].astype(str).str.strip()
    ann["endpoint_label"] = ann[args.endpoint].fillna("NA").astype(str).str.strip()
    ann = ann[ann["sample_id"].isin(set(matrix.columns))].copy()
    excluded: set[str] = set()
    for path in args.exclude_samples:
        if path.exists():
            excluded.update(pd.read_csv(path)["sample_id"].astype(str))
    selected = choose_unseen_samples(ann, matrix.columns.tolist(), excluded, args.n_samples, args.seed)
    values = matrix.to_numpy(dtype=np.float32)
    samples = matrix.columns.tolist()
    sample_pos = {sample_id: j for j, sample_id in enumerate(samples)}
    labels = ann.set_index("sample_id").reindex(samples)["endpoint_label"].to_numpy(dtype=str)
    groups = sorted(set(labels))
    detail: dict[str, list[dict[str, Any]]] = {"network_correlation_all_genes": [], "network_discriminative_correlation_top200": []}
    probabilities: dict[str, list[dict[str, Any]]] = {key: [] for key in detail}
    manifest: list[dict[str, Any]] = []
    gene_symbols = matrix.index.astype(str).to_numpy()
    for fold_no, row in enumerate(selected.itertuples(index=False), start=1):
        sample_id = str(row.sample_id)
        heldout_idx = sample_pos[sample_id]
        truth = str(row.endpoint_label)
        reference, training = build_group_reference(values, labels, groups, heldout_idx)
        selected_genes, audit = select_group_discriminative_genes(values, groups, training, args.top_n_genes)
        audit["gene_symbol"] = gene_symbols[audit["gene_index"].to_numpy(dtype=int)]
        gene_file = args.outdir / "fold_selected_genes" / f"fold_{fold_no:02d}_{sample_id}_genes.csv"
        audit.to_csv(gene_file, index=False, encoding="utf-8-sig")
        for route, genes in [("network_correlation_all_genes", None), ("network_discriminative_correlation_top200", selected_genes)]:
            output, probability = evaluate_route(sample_id, truth, values[:, heldout_idx], reference, groups, genes, route)
            detail[route].append(output)
            probabilities[route].append(probability)
        manifest.append({"fold": fold_no, "sample_id": sample_id, "truth_label": truth, "gene_file": str(gene_file.relative_to(args.outdir))})
    detail_frames = {key: pd.DataFrame(value) for key, value in detail.items()}
    prob_frames = {key: pd.DataFrame(value) for key, value in probabilities.items()}
    routes = {key: summarize(detail_frames[key], prob_frames[key], groups) for key in detail}
    enhanced = routes["network_discriminative_correlation_top200"]
    achieved = enhanced["top1_accuracy"] >= args.target_top1 and enhanced["top3_accuracy"] >= args.target_top3
    full_loso = len(selected) == matrix.shape[1] and not excluded
    if achieved and full_loso:
        decision = "网络层级判别 correlation 在全量严格 LOSO 中达到 Top1/Top3 目标，可切换为 Bo2023 正式输出的上层来源主结论。"
    elif achieved:
        decision = "网络层级判别 correlation 在独立未见样本上达到 Top1/Top3 目标，可作为 Bo2023 的上层来源主路径进入扩大复核。"
    else:
        decision = "网络层级判别 correlation 未同时达到 Top1/Top3 目标，不能升级为正式路径。"
    boundary = (
        "该达标只适用于 `SaleemNetworks` 上层终点，不等价于 110 个精确 Region 的 Top1/Top3 达标。精确 Region 仍应作为二级探索输出，附带不确定性。"
        if achieved
        else "当前证据不足以支持修改生产输出层级。"
    )
    summary = {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "dataset": "Bo2023_WangLab_VSD_region",
        "endpoint": args.endpoint,
        "n_classes": int(len(groups)),
        "n_test_samples": int(len(selected)),
        "n_prior_samples_excluded": int(len(excluded)),
        "validation_scope": (
            f"全部 {len(selected)} 个可评估已知样本的严格 LOSO（每次仅留出当前测试样本）"
            if full_loso
            else f"新的 {len(selected)} 个未见样本，排除前序 {len(excluded)} 个已检查样本"
        ),
        "seed": int(args.seed),
        "top_n_genes": int(args.top_n_genes),
        "target": {"top1": float(args.target_top1), "top3": float(args.target_top3)},
        "routes": routes,
        "target_achieved": bool(achieved),
        "decision": decision,
        "boundary": boundary,
    }
    pd.DataFrame([{"route": key, **value} for key, value in routes.items()]).to_csv(args.outdir / "route_metrics.csv", index=False, encoding="utf-8-sig")
    selected.to_csv(args.outdir / "selected_test_samples.csv", index=False, encoding="utf-8-sig")
    for key, frame in detail_frames.items():
        frame.to_csv(args.outdir / f"{key}_detail.csv", index=False, encoding="utf-8-sig")
    (args.outdir / "fold_manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")
    (args.outdir / "validation_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    export_plot(args.outdir, pd.DataFrame([{"route": key, **value} for key, value in routes.items()]), args.target_top1, args.target_top3)
    write_report(args.outdir, summary)
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    print(f"Outputs written to: {args.outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
