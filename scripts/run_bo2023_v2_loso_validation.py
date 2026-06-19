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

from reporting.benchmark_figure_export import export_benchmark_paper_figures  # noqa: E402
from core.reference_projection import clean_excel_date_gene_symbol  # noqa: E402
from scripts.run_bo2023_loso_validation import (  # noqa: E402
    build_region_reference,
    export_loso_dashboard,
    make_suite,
    read_annotations,
    read_vsd_matrix,
)
from signature_builder import (  # noqa: E402
    DEFAULT_BLOOD_BACKGROUND,
    DEFAULT_HOUSEKEEPING,
    _compute_region_scores,
    _select_signature_genes,
)
from source_tracing_v2 import SourceTracingEngineV2  # noqa: E402


DEFAULT_MATRIX = ROOT / "bo2023 data" / "mfas5_819samples_23605genes_vsd4_rmbatch.xls"
DEFAULT_SAMPLE_INFO = ROOT / "bo2023 data" / "Information of sequenced samples_update_full878_filter819.xlsx"
DEFAULT_GENE_MAP = ROOT / "bo2023_bulk_atlas_buildkit" / "04_expressed_genes_neocortex_plus_subcortical.csv"
DEFAULT_BASELINE = ROOT / "results" / "bo2023_loso_30_vsd_correlation" / "tables" / "benchmark_metrics.csv"
DEFAULT_OUTDIR = ROOT / "results" / "bo2023_loso_30_v2_ensemble_adapted"


def map_matrix_to_symbols(matrix: pd.DataFrame, gene_map_path: Path) -> pd.DataFrame:
    mapping = pd.read_csv(gene_map_path, usecols=["Gene.stable.ID", "Gene.name"])
    mapping["Gene.stable.ID"] = mapping["Gene.stable.ID"].astype(str).str.strip()
    mapping["Gene.name"] = mapping["Gene.name"].map(clean_excel_date_gene_symbol)
    symbol_by_id = mapping.drop_duplicates("Gene.stable.ID").set_index("Gene.stable.ID")["Gene.name"]
    gene_symbols = matrix.index.to_series().map(symbol_by_id).replace("", pd.NA).fillna(matrix.index.to_series())
    out = matrix.groupby(gene_symbols, sort=True).mean()
    out.index.name = "gene_symbol"
    return out


def choose_validation_samples(
    ann: pd.DataFrame,
    matrix_samples: list[str],
    n_samples: int,
    seed: int,
    excluded_sample_ids: set[str],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    aligned = ann[ann["sample_id"].isin(set(matrix_samples))].copy()
    region_counts = aligned.groupby("region_id")["sample_id"].size()
    singleton_regions = set(region_counts[region_counts < 2].index)
    singleton_samples = aligned[aligned["region_id"].isin(singleton_regions)].copy()
    eligible = aligned[
        ~aligned["region_id"].isin(singleton_regions)
        & ~aligned["sample_id"].isin(excluded_sample_ids)
    ].copy()
    if n_samples > len(eligible):
        raise ValueError(f"requested {n_samples} samples but only {len(eligible)} are selectable")
    selected = eligible.sample(n=n_samples, random_state=seed).sort_values("sample_id").reset_index(drop=True)
    return selected, singleton_samples


def select_fold_signature(
    genes: list[str],
    regions: list[str],
    reference: np.ndarray,
    topk_per_region: int,
) -> tuple[np.ndarray, pd.DataFrame]:
    # The Bo2023 matrix is already VSD-normalized. Applying the TPM-oriented
    # log1p transform again would distort fold-local specificity estimates.
    x = reference.astype(float)
    score_matrix = _compute_region_scores(genes, regions, x, "hybrid_specificity")
    signature_rows = _select_signature_genes(
        genes,
        regions,
        score_matrix,
        x,
        DEFAULT_HOUSEKEEPING | DEFAULT_BLOOD_BACKGROUND,
        topk_per_region,
        -np.inf,
        3,
    )
    signature_df = pd.DataFrame(signature_rows, columns=["region_id", "gene_symbol", "raw_specificity_score"])
    signature_df["marker_score"] = signature_df.groupby("region_id")["raw_specificity_score"].rank(pct=True)
    signature_df["specificity_score"] = signature_df["marker_score"]
    signature_df["is_marker"] = 1
    selected = set(signature_df["gene_symbol"].astype(str))
    mask = np.asarray([gene in selected for gene in genes], dtype=bool)
    return mask, signature_df


def build_fold_marker_weight_matrix(
    fold_genes: np.ndarray,
    regions: list[str],
    signature_df: pd.DataFrame,
    specificity_weight: float = 1.25,
    marker_weight: float = 0.75,
    marker_bonus: float = 0.75,
) -> np.ndarray:
    weights = np.ones((len(fold_genes), len(regions)), dtype=float)
    gene_pos = {str(gene): idx for idx, gene in enumerate(fold_genes)}
    region_pos = {region: idx for idx, region in enumerate(regions)}
    for row in signature_df.itertuples(index=False):
        i = gene_pos.get(str(row.gene_symbol))
        j = region_pos.get(str(row.region_id))
        if i is None or j is None:
            continue
        weights[i, j] = (
            1.0
            + specificity_weight * max(float(row.specificity_score), 0.0)
            + marker_weight * max(float(row.marker_score), 0.0)
            + marker_bonus * float(row.is_marker)
        )
    col_means = weights.mean(axis=0, keepdims=True)
    col_means[col_means <= 0] = 1.0
    return weights / col_means


class InMemoryVsdFoldEngine(SourceTracingEngineV2):
    """Run the production V2 scoring path on one leakage-free fold."""

    def __init__(
        self,
        genes: np.ndarray,
        regions: list[str],
        reference: np.ndarray,
        sample_id: str,
        sample_values: np.ndarray,
        marker_weights: np.ndarray | None = None,
        n_weighted_pairs: int = 0,
    ) -> None:
        super().__init__(":memory:")
        self.fold_genes = genes.astype(object)
        self.fold_regions = list(regions)
        self.fold_reference = np.asarray(reference, dtype=float)
        self.fold_sample_id = sample_id
        self.fold_sample_values = np.asarray(sample_values, dtype=float)
        self.fold_marker_weights = (
            np.asarray(marker_weights, dtype=float)
            if marker_weights is not None
            else np.ones_like(self.fold_reference, dtype=float)
        )
        self.fold_weighting_active = marker_weights is not None
        self.fold_n_weighted_pairs = int(n_weighted_pairs)

    def _get_atlas_meta(self, atlas_id: int) -> dict[str, Any]:
        return {
            "atlas_id": atlas_id,
            "atlas_name": "Bo2023_WangLab_VSD_region_LOSO_fold",
            "normalization": "VSD_batch_removed",
            "build_version": "strict_loso_fold",
        }

    def _get_latest_sigset_id(self, atlas_id: int = 1) -> None:
        return None

    def _load_reference_bundle(self, atlas_id: int, sigset_id: int | None, use_value: str, *args: Any, **kwargs: Any):
        a_base = self._apply_value_transform(self.fold_reference, use_value)
        w_ref = self.fold_marker_weights.copy()
        return (
            self.fold_genes.copy(),
            a_base,
            w_ref,
            self.fold_regions.copy(),
            {
                "weighting_active": self.fold_weighting_active,
                "n_weighted_pairs": self.fold_n_weighted_pairs,
                "reference_gene_id_type": "symbol_like",
                "_A_tpm_cache": self.fold_reference.copy(),
            },
        )

    def _load_sample_features(self, sample_id: str, genes: np.ndarray) -> pd.DataFrame:
        if sample_id != self.fold_sample_id:
            return pd.DataFrame()
        values = self.fold_sample_values
        feat = pd.DataFrame(index=pd.Index(genes, name="gene_symbol"))
        feat["tpm_value"] = values
        feat["log_tpm"] = values
        feat["zscore_tpm"] = values
        feat["detected"] = 0.0
        feat["read_count"] = np.nan
        feat["rank_pct"] = feat["tpm_value"].rank(pct=True, method="average").fillna(0.0)
        feat["read_support"] = np.nan
        feat.attrs["sample_gene_id_type"] = "symbol_like"
        return feat


def run_fold(
    fold_no: int,
    sample_id: str,
    truth: str,
    sample: np.ndarray,
    reference: np.ndarray,
    regions: list[str],
    genes: list[str],
    topk_per_region: int,
    bootstrap_n: int,
    bootstrap_gene_fraction: float,
    corr_weight: float | None,
    nnls_weight: float | None,
    rank_weight: float | None,
    ensemble_alpha: float | None,
    method: str,
    fold_marker_weights: bool,
    outdir: Path,
) -> tuple[dict[str, Any], dict[str, Any], pd.DataFrame, dict[str, Any]]:
    mask, signature_df = select_fold_signature(genes, regions, reference, topk_per_region)
    fold_genes = np.asarray(genes, dtype=object)[mask]
    marker_weights = (
        build_fold_marker_weight_matrix(fold_genes, regions, signature_df)
        if fold_marker_weights
        else None
    )
    engine = InMemoryVsdFoldEngine(
        fold_genes,
        regions,
        reference[mask, :],
        sample_id,
        sample[mask],
        marker_weights=marker_weights,
        n_weighted_pairs=(len(signature_df) if fold_marker_weights else 0),
    )
    optional_weights: dict[str, float] = {}
    if corr_weight is not None:
        optional_weights["corr_component_weight"] = corr_weight
    if nnls_weight is not None:
        optional_weights["nnls_component_weight"] = nnls_weight
    if rank_weight is not None:
        optional_weights["rank_component_weight"] = rank_weight
    if ensemble_alpha is not None:
        optional_weights["ensemble_alpha"] = ensemble_alpha
    output = engine.trace(
        sample_id=sample_id,
        method=method,
        atlas_id=4,
        use_value="vsd",
        vsd_compatible=True,
        return_all=True,
        persist=False,
        bootstrap_n=bootstrap_n,
        bootstrap_gene_frac=bootstrap_gene_fraction,
        random_seed=42 + fold_no,
        l2=1e-4,
        min_overlap_genes=5,
        min_overlap_fraction=0.003,
        abstain_on_low_overlap=True,
        **optional_weights,
    )
    ranked = output.get("results", [])
    if not ranked:
        raise RuntimeError(f"fold {fold_no} returned no prediction: {output.get('meta', {})}")
    predictions = [x["region_id"] for x in ranked]
    true_rank = predictions.index(truth) + 1
    result = {
        "fold": fold_no,
        "sample_id": sample_id,
        "label": truth,
        "pred_top1": predictions[0],
        "pred_top2": predictions[1],
        "pred_top3": predictions[2],
        "true_rank": true_rank,
        "hit1": int(true_rank == 1),
        "hit3": int(true_rank <= 3),
        "top1_score": float(ranked[0]["score"]),
        "true_region_score": float(ranked[true_rank - 1]["score"]),
        "top1_confidence": float(ranked[0]["confidence"]),
        "top2_confidence": float(ranked[1]["confidence"]),
        "decision_margin": float(output["meta"]["decision_margin"]),
        "top1_stability": float(ranked[0].get("stability", np.nan)),
        "abstained": 0,
        "traceability": f"strict_loso_v2_{method}_vsd",
        "overlap_genes": int(mask.sum()),
        "n_fold_signature_genes": int(mask.sum()),
        "signal_strategy": output["meta"].get("signal_strategy"),
        "fold_marker_weighting": bool(fold_marker_weights),
    }
    probability = {"sample_id": sample_id, "label": truth}
    probability.update({x["region_id"]: float(x["confidence"]) for x in ranked})
    score_df = pd.DataFrame(
        [
            {
                "fold": fold_no,
                "sample_id": sample_id,
                "truth_region": truth,
                "candidate_region": item["region_id"],
                "rank": item["rank"],
                "score": item["score"],
                "confidence": item["confidence"],
                "fraction": item.get("fraction"),
                "stability": item.get("stability"),
                "corr_component": item.get("corr_component"),
                "nnls_component": item.get("nnls_component"),
                "rank_component": item.get("rank_component"),
                "marker_component": item.get("marker_component"),
                "detect_component": item.get("detect_component"),
                "support_component": item.get("support_component"),
            }
            for item in ranked
        ]
    )
    signature_path = outdir / "fold_signatures" / f"fold_{fold_no:02d}_{sample_id}_signature.csv"
    signature_df.to_csv(signature_path, index=False, encoding="utf-8-sig")
    manifest = {
        "fold": fold_no,
        "sample_id": sample_id,
        "truth_region": truth,
        "n_signature_rows": int(len(signature_df)),
        "n_signature_genes": int(mask.sum()),
        "signature_file": str(signature_path.relative_to(outdir)),
        "signal_strategy": output["meta"].get("signal_strategy"),
        "component_weights": output["meta"].get("signal_component_weights"),
        "use_value": output["meta"].get("use_value"),
        "fold_marker_weighting": bool(fold_marker_weights),
        "reference_weighting": output["meta"].get("reference_weighting"),
    }
    return result, probability, score_df, manifest


def write_report(
    outdir: Path,
    summary: dict[str, Any],
    detail: pd.DataFrame,
    baseline_metrics: pd.DataFrame | None,
) -> None:
    baseline_top1 = baseline_top3 = None
    if baseline_metrics is not None and not baseline_metrics.empty:
        lookup = baseline_metrics.set_index("metric")["value"]
        baseline_top1 = float(lookup.get("Top1_acc_valid", np.nan))
        baseline_top3 = float(lookup.get("Top3_acc_valid", np.nan))
    strategy = detail["signal_strategy"].mode().iloc[0]
    marker_weighting = bool(summary.get("fold_marker_weights", False))
    comparison = ""
    if baseline_top1 is not None and baseline_top3 is not None:
        comparison = (
            f"- 严格 correlation baseline：Top1 `{baseline_top1:.1%}`，Top3 `{baseline_top3:.1%}`。\n"
            f"- 适配后的 V2 ensemble：Top1 `{summary['top1_acc']:.1%}`，Top3 `{summary['top3_acc']:.1%}`。\n"
        )
    text = f"""# Bo2023 图谱与正式 V2 溯源算法适配验证

## 已实施的适配

1. 正式引擎使用 `use_value="vsd"`：对于与 Bo2023 atlas 同尺度的 VSD 样本，直接在 VSD 表达空间计算，避免再次按 TPM/log 路径转换。
2. 本轮折内 marker 权重状态：`{marker_weighting}`。权重仅由每折训练样本构建，不读取留出样本。
3. VSD 输入没有可靠二元 detected 观测，因此 detection 组件保持禁用；marker 与 rank-support 可在有折内权重时进入 ensemble。

## 30 样本严格 LOSO 结果

| 指标 | 适配后的 V2 ensemble |
| --- | ---: |
| Top1 命中 | {int(detail['hit1'].sum())}/30 ({summary['top1_acc']:.1%}) |
| Top3 命中 | {int(detail['hit3'].sum())}/30 ({summary['top3_acc']:.1%}) |
| Macro AUC (OVR) | {summary['macro_auc']:.3f} |
| 真实 Region 中位排名 | {float(detail['true_rank'].median()):.1f} |
| 平均 Top1 稳定性 | {summary['mean_stability']:.3f} |
| 折内 signature gene 数量均值 | {float(detail['n_fold_signature_genes'].mean()):.1f} |

{comparison}
正式 V2 在本次运行实际使用的可辨识信号为 `{strategy}`。

## 结论与下一轮条件

- 若适配后的 V2 Top3 不低于 correlation baseline，且 Top1 有提升或基本持平，则可以开启扩大样本量的确认测试。
- 若 marker-weighted V2 低于同队列 correlation，应保留 correlation 为主路径，并回查 marker 权重构建策略。
- 当前结果只针对 VSD 同尺度的 Bo2023 留出样本；对真实 cfRNA TPM/count 上传样本，仍需要独立的跨平台校准验证。
"""
    (outdir / "v2_adaptation_report_cn.md").write_text(text, encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description="Strict Bo2023 LOSO using the adapted production SourceTracingEngineV2.")
    parser.add_argument("--matrix", type=Path, default=DEFAULT_MATRIX)
    parser.add_argument("--sample-info", type=Path, default=DEFAULT_SAMPLE_INFO)
    parser.add_argument("--gene-map", type=Path, default=DEFAULT_GENE_MAP)
    parser.add_argument("--sample-sheet", default="mfas5_819samples_phenSet4")
    parser.add_argument("--region-col", default="Region")
    parser.add_argument("--n-samples", type=int, default=30)
    parser.add_argument("--method", choices=["ensemble", "correlation"], default="ensemble")
    parser.add_argument("--fold-marker-weights", action="store_true")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--topk-per-region", type=int, default=80)
    parser.add_argument("--bootstrap-n", type=int, default=0)
    parser.add_argument("--bootstrap-gene-fraction", type=float, default=0.7)
    parser.add_argument("--corr-weight", type=float)
    parser.add_argument("--nnls-weight", type=float)
    parser.add_argument("--rank-weight", type=float)
    parser.add_argument("--ensemble-alpha", type=float)
    parser.add_argument("--baseline-metrics", type=Path, default=DEFAULT_BASELINE)
    parser.add_argument("--exclude-samples", type=Path, action="append", default=[])
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    (args.outdir / "fold_signatures").mkdir(parents=True, exist_ok=True)
    raw_matrix = read_vsd_matrix(args.matrix)
    matrix = map_matrix_to_symbols(raw_matrix, args.gene_map)
    ann = read_annotations(args.sample_info, args.sample_sheet, args.region_col)
    ann = ann[ann["sample_id"].isin(set(matrix.columns))].copy()
    excluded_sample_ids: set[str] = set()
    for exclude_path in args.exclude_samples:
        if exclude_path.exists():
            excluded_sample_ids.update(pd.read_csv(exclude_path)["sample_id"].astype(str))
    selected, excluded = choose_validation_samples(
        ann, matrix.columns.tolist(), args.n_samples, args.seed, excluded_sample_ids
    )
    values = matrix.to_numpy(dtype=np.float32)
    full_reference, regions, region_counts, region_indices = build_region_reference(values, matrix.columns.tolist(), ann)
    sample_pos = {sample_id: idx for idx, sample_id in enumerate(matrix.columns)}
    region_pos = {region: idx for idx, region in enumerate(regions)}

    detail_rows: list[dict[str, Any]] = []
    probability_rows: list[dict[str, Any]] = []
    score_frames: list[pd.DataFrame] = []
    fold_manifest: list[dict[str, Any]] = []
    for fold_no, heldout in enumerate(selected.itertuples(index=False), start=1):
        sample_id = str(heldout.sample_id)
        truth = str(heldout.region_id)
        test_idx = sample_pos[sample_id]
        truth_idx = region_pos[truth]
        reference = full_reference.copy()
        train_idx = region_indices[truth][region_indices[truth] != test_idx]
        reference[:, truth_idx] = values[:, train_idx].mean(axis=1, dtype=np.float64).astype(np.float32)
        result, probability, score_df, manifest = run_fold(
            fold_no,
            sample_id,
            truth,
            values[:, test_idx],
            reference,
            regions,
            matrix.index.astype(str).tolist(),
            args.topk_per_region,
            args.bootstrap_n,
            args.bootstrap_gene_fraction,
            args.corr_weight,
            args.nnls_weight,
            args.rank_weight,
            args.ensemble_alpha,
            args.method,
            args.fold_marker_weights,
            args.outdir,
        )
        result["truth_region_train_samples"] = int(region_counts[truth_idx] - 1)
        detail_rows.append(result)
        probability_rows.append(probability)
        score_frames.append(score_df)
        fold_manifest.append(manifest)

    detail_df = pd.DataFrame(detail_rows)
    probability_df = pd.DataFrame(probability_rows)
    suite = make_suite(detail_df, probability_df, sorted(detail_df["label"].unique()), len(regions))
    metric_lookup = suite["metrics_df"].set_index("metric")["value"]
    summary = {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "validation_design": "strict leave-one-sample-out with fold-specific signatures",
        "dataset": "Bo2023_WangLab_VSD_region",
        "engine": "SourceTracingEngineV2",
        "method": args.method,
        "adaptation": "vsd_direct_input_and_fold_specific_marker_weighting",
        "n_sample_total": int(matrix.shape[1]),
        "n_test_samples": int(len(detail_df)),
        "n_regions": int(len(regions)),
        "n_gene_symbols": int(matrix.shape[0]),
        "n_singleton_samples_excluded": int(len(excluded)),
        "n_prior_test_samples_excluded": int(len(excluded_sample_ids)),
        "seed": int(args.seed),
        "topk_per_region": int(args.topk_per_region),
        "bootstrap_n": int(args.bootstrap_n),
        "fold_marker_weights": bool(args.fold_marker_weights),
        "requested_component_weights": {
            "corr": args.corr_weight,
            "nnls": args.nnls_weight,
            "rank": args.rank_weight,
            "ensemble_alpha": args.ensemble_alpha,
        },
        "top1_acc": float(metric_lookup["Top1_acc_valid"]),
        "top3_acc": float(metric_lookup["Top3_acc_valid"]),
        "macro_auc": float(metric_lookup["MacroAUC_ovr_valid"]),
        "mean_stability": float(metric_lookup["Mean_top1_stability_valid"]),
    }
    metadata = {
        "validation_design": summary["validation_design"],
        "dataset": summary["dataset"],
        "method": args.method,
        "method_description": "SourceTracingEngineV2 with VSD-direct input, fold-specific signature genes and optional fold-specific marker weights",
        "use_value": "vsd",
        "n_test_samples": len(detail_df),
        "n_regions": len(regions),
        "bootstrap_n": args.bootstrap_n,
    }
    export_benchmark_paper_figures(args.outdir, suite, metadata=metadata, prefix="bo2023_v2_loso30")
    export_loso_dashboard(
        args.outdir,
        detail_df,
        len(regions),
        title=f"Bo2023 strict LOSO: SourceTracingEngineV2 {args.method}, fold marker weights={args.fold_marker_weights}",
    )
    pd.concat(score_frames, ignore_index=True).to_csv(args.outdir / "all_candidate_scores.csv", index=False, encoding="utf-8-sig")
    selected.to_csv(args.outdir / "selected_test_samples.csv", index=False, encoding="utf-8-sig")
    excluded.to_csv(args.outdir / "singleton_samples_not_loso_assessable.csv", index=False, encoding="utf-8-sig")
    (args.outdir / "fold_manifest.json").write_text(json.dumps(fold_manifest, ensure_ascii=False, indent=2), encoding="utf-8")
    (args.outdir / "validation_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    baseline_metrics = pd.read_csv(args.baseline_metrics) if args.baseline_metrics.exists() else None
    write_report(args.outdir, summary, detail_df, baseline_metrics)
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    print(f"Outputs written to: {args.outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
