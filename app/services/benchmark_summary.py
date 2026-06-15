from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, Optional

import numpy as np
import pandas as pd


@dataclass(frozen=True)
class BenchmarkThresholds:
    top1_high: float = 0.75
    top1_mid: float = 0.50
    topk_gap_good: float = 0.15
    auc_high: float = 0.80
    auc_mid: float = 0.65
    stability_high: float = 0.80
    stability_mid: float = 0.50
    confidence_high: float = 0.70
    confidence_mid: float = 0.45
    margin_high: float = 0.15
    margin_mid: float = 0.05
    failure_high: float = 0.30


THRESHOLDS = BenchmarkThresholds()

METRIC_GUIDE = [
    {"指标": "Top1 accuracy", "说明": "真实来源脑区排在第 1 位的比例。越高越好，最能代表模型的最终精确定位能力。"},
    {"指标": "Top3 accuracy", "说明": "真实来源脑区是否进入前 3 名候选。若明显高于 Top1，说明模型更擅长缩小候选范围。"},
    {"指标": "Rank", "说明": "真实脑区在候选列表中的排序位置。大量样本集中在 rank 1-3，通常说明模型具备实用筛查价值。"},
    {"指标": "Confusion matrix", "说明": "显示真实脑区与预测脑区之间的对应关系。对角线越集中，说明模型越稳定。"},
    {"指标": "ROC / AUC", "说明": "评估各脑区与其他脑区的可区分程度。AUC 越高，排序区分能力越强。"},
    {"指标": "Confidence", "说明": "模型对 Top1 结果的相对把握。适合用于结果分层和人工复核优先级判断。"},
    {"指标": "Margin", "说明": "Top1 与 Top2 分数差距。差距越大，说明第一候选领先越明显。"},
    {"指标": "Stability", "说明": "Bootstrap 后 Top1 是否保持一致。稳定性高说明结果不太依赖少数基因。"},
    {"指标": "Failure mode", "说明": "将错误、弃权、低 overlap、低 confidence 等归类，帮助快速定位主要问题来源。"},
]


def _safe_float(value: Any) -> Optional[float]:
    try:
        if value is None or pd.isna(value):
            return None
        v = float(value)
        return v if np.isfinite(v) else None
    except Exception:
        return None


def _fmt(value: Any, digits: int = 3) -> str:
    v = _safe_float(value)
    return "NA" if v is None else f"{v:.{digits}f}"


def _metric_map(metrics_df: pd.DataFrame) -> Dict[str, Any]:
    if metrics_df is None or metrics_df.empty or not {"metric", "value"}.issubset(metrics_df.columns):
        return {}
    return dict(zip(metrics_df["metric"].astype(str), metrics_df["value"]))


def _level(value: Optional[float], high: float, mid: float) -> str:
    if value is None:
        return "unknown"
    if value >= high:
        return "high"
    if value >= mid:
        return "mid"
    return "low"


def _level_cn(level: str) -> str:
    return {"high": "较高", "mid": "中等", "low": "偏低", "unknown": "暂不可评估"}.get(level, level)


def _is_similar_region_pair(a: Any, b: Any) -> bool:
    x = str(a or "").upper()
    y = str(b or "").upper()
    groups = [
        ("M1", "S1", "PMC", "PM", "SMA"),
        ("V1", "V2", "V3", "VIS"),
        ("PFC", "DLPFC", "OFC", "ACC"),
        ("HIP", "HPC", "CA", "DG", "EC"),
    ]
    return any(any(token in x for token in group) and any(token in y for token in group) for group in groups)


def compute_true_rank_df(suite: Dict[str, pd.DataFrame]) -> pd.DataFrame:
    prob = suite.get("probability_df", pd.DataFrame())
    if prob is None or prob.empty or "label" not in prob.columns:
        return pd.DataFrame()
    class_cols = [c for c in prob.columns if c not in {"sample_id", "label"}]
    rows = []
    for _, row in prob.iterrows():
        label = str(row.get("label"))
        if label not in class_cols:
            rows.append({"sample_id": row.get("sample_id"), "label": label, "true_rank": np.nan, "true_label_score": np.nan})
            continue
        scores = pd.to_numeric(row[class_cols], errors="coerce").fillna(0.0)
        ordered = scores.sort_values(ascending=False)
        rows.append(
            {
                "sample_id": row.get("sample_id"),
                "label": label,
                "true_rank": int(list(ordered.index).index(label) + 1),
                "true_label_score": float(scores[label]),
            }
        )
    return pd.DataFrame(rows)


def build_failure_mode_df(suite: Dict[str, pd.DataFrame]) -> pd.DataFrame:
    detail = suite.get("detail_df", pd.DataFrame())
    if detail is None or detail.empty:
        return pd.DataFrame()
    rows = []
    for _, row in detail.iterrows():
        if int(row.get("abstained", 0) or 0) == 1:
            mode = "abstained / evidence insufficient"
        elif int(row.get("hit1", 0) or 0) == 1:
            mode = "correct Top1"
        elif _safe_float(row.get("overlap_genes")) is not None and float(row.get("overlap_genes")) < 20:
            mode = "wrong with low overlap"
        elif _safe_float(row.get("decision_margin")) is not None and float(row.get("decision_margin")) < THRESHOLDS.margin_mid:
            mode = "wrong with low margin"
        elif _safe_float(row.get("top1_confidence")) is not None and float(row.get("top1_confidence")) < THRESHOLDS.confidence_mid:
            mode = "wrong with low confidence"
        else:
            mode = "wrong despite usable evidence"
        rows.append({"failure_mode": mode})
    counts = pd.DataFrame(rows).value_counts("failure_mode").reset_index(name="n_samples")
    counts["fraction"] = counts["n_samples"] / max(int(counts["n_samples"].sum()), 1)
    return counts


def _top_confusions(conf_long: pd.DataFrame, n: int = 5) -> pd.DataFrame:
    if conf_long is None or conf_long.empty:
        return pd.DataFrame()
    wrong = conf_long[conf_long["truth_region"].astype(str) != conf_long["pred_region"].astype(str)].copy()
    if wrong.empty:
        return wrong
    wrong["similar_region_pair"] = wrong.apply(lambda r: _is_similar_region_pair(r.get("truth_region"), r.get("pred_region")), axis=1)
    return wrong.sort_values(["count", "row_fraction"], ascending=False).head(n)


def build_benchmark_insights(suite: Dict[str, pd.DataFrame], meta: Optional[Dict[str, Any]] = None, k: int = 3) -> Dict[str, Any]:
    meta = meta or {}
    detail = suite.get("detail_df", pd.DataFrame())
    metrics_df = suite.get("metrics_df", pd.DataFrame())
    conf_long = suite.get("confusion_long_df", pd.DataFrame())
    metrics = _metric_map(metrics_df)
    summary = detail.attrs.get("summary", {}) if hasattr(detail, "attrs") else {}

    top1 = _safe_float(summary.get("top1_acc")) or _safe_float(metrics.get("Top1_acc_valid"))
    topk = _safe_float(summary.get(f"top{k}_acc")) or _safe_float(metrics.get(f"Top{k}_acc_valid"))
    auc = _safe_float(summary.get("auc")) or _safe_float(metrics.get("MacroAUC_ovr_valid"))
    stability = _safe_float(summary.get("mean_top1_stability")) or _safe_float(metrics.get("Mean_top1_stability_valid"))
    confidence = _safe_float(summary.get("mean_top1_confidence")) or _safe_float(metrics.get("Mean_top1_confidence_valid"))
    margin = _safe_float(summary.get("mean_decision_margin")) or _safe_float(metrics.get("Mean_decision_margin_valid"))
    abstain = _safe_float(summary.get("abstain_rate")) or _safe_float(metrics.get("Abstain_rate"))

    top_conf = _top_confusions(conf_long)
    similar_confusions = bool(not top_conf.empty and top_conf["similar_region_pair"].mean() >= 0.5)
    if top_conf.empty:
        major_confusion_text = "暂无明显集中混淆"
    else:
        top_pairs = [f"{row.truth_region} -> {row.pred_region}" for row in top_conf.itertuples()]
        major_confusion_text = "；".join(top_pairs[:2])

    top1_level = _level(top1, THRESHOLDS.top1_high, THRESHOLDS.top1_mid)
    auc_level = _level(auc, THRESHOLDS.auc_high, THRESHOLDS.auc_mid)
    stability_level = _level(stability, THRESHOLDS.stability_high, THRESHOLDS.stability_mid)
    margin_level = _level(margin, THRESHOLDS.margin_high, THRESHOLDS.margin_mid)
    topk_gap = (topk - top1) if top1 is not None and topk is not None else None

    if top1_level == "high":
        one_liner = "当前模型已经具备较好的 Top1 脑区识别能力。"
    elif topk_gap is not None and topk_gap >= THRESHOLDS.topk_gap_good:
        one_liner = f"当前模型更擅长缩小候选范围：Top{k} 明显高于 Top1，但精细区分仍有改进空间。"
    elif top1_level == "mid":
        one_liner = "当前模型具备一定识别能力，但相邻或相似脑区之间仍存在稳定混淆。"
    else:
        one_liner = "当前模型的精细脑区定位能力偏弱，建议优先检查 reference、signature 和样本 overlap。"

    if similar_confusions:
        main_problem = f"主要误差集中在相邻或相似脑区之间，例如 {major_confusion_text}。"
        next_step = "建议优先强化这些相近脑区的 marker 权重、signature filtering 或 cortex 内部细分参考。"
    elif abstain is not None and abstain > 0.20:
        main_problem = "当前存在较高的 abstain 比例，说明部分样本证据不足或与参考图谱重叠有限。"
        next_step = "建议检查 gene overlap、输入基因 ID 类型，以及参考图谱对样本的覆盖度。"
    elif top_conf.empty:
        main_problem = "当前未见特别集中的单一错误模式。"
        next_step = "建议保持当前参数，并用独立样本或替代 atlas 继续复核。"
    else:
        main_problem = f"主要误差方向包括 {major_confusion_text}，需要排查是否存在跨系统或跨叶别混淆。"
        next_step = "建议优先检查标签一致性、atlas 层级、以及 signature 与 atlas 的匹配情况。"

    failure_df = build_failure_mode_df(suite)
    failure_high = bool(not failure_df.empty and float(failure_df.loc[failure_df["failure_mode"] != "correct Top1", "fraction"].sum()) >= THRESHOLDS.failure_high)

    return {
        "metrics": {"top1": top1, "topk": topk, "auc": auc, "stability": stability, "confidence": confidence, "margin": margin, "abstain": abstain},
        "levels": {"top1": top1_level, "auc": auc_level, "stability": stability_level, "margin": margin_level},
        "topk_gap": topk_gap,
        "top_confusions": top_conf,
        "similar_confusions": similar_confusions,
        "major_confusion_text": major_confusion_text,
        "failure_df": failure_df,
        "failure_high": failure_high,
        "one_liner": one_liner,
        "main_problem": main_problem,
        "next_step": next_step,
        "method": meta.get("method"),
        "k": k,
    }


def explanation_for(module: str, insights: Dict[str, Any]) -> Dict[str, str]:
    m = insights.get("metrics", {})
    levels = insights.get("levels", {})
    k = int(insights.get("k", 3) or 3)
    topk_gap = insights.get("topk_gap")
    confusion = insights.get("major_confusion_text", "暂无明显集中混淆")
    similar = insights.get("similar_confusions", False)
    failure_high = insights.get("failure_high", False)

    if module == "accuracy":
        interpretation = (
            f"当前 Top{k} 比 Top1 高 {_fmt(topk_gap)}，说明模型通常能把真实来源缩小到较小候选范围。"
            if topk_gap is not None and topk_gap >= THRESHOLDS.topk_gap_good
            else f"当前 Top1 为 {_fmt(m.get('top1'))}，Top{k} 为 {_fmt(m.get('topk'))}，候选排序与最终判定较为一致。"
        )
        return {
            "what": f"Top1 / Top{k} accuracy 衡量真实脑区是否排在第 1 位或前 {k} 位。",
            "how": "先看 Top1 判断最终精细定位能力，再看 TopK 判断候选范围收缩能力。",
            "good": "Top1 >= 0.75 通常较好；如果 Top1 中等但 TopK 较高，说明模型具备较强候选筛查价值。",
            "interpretation": interpretation,
        }
    if module == "confusion":
        interpretation = (
            f"当前主要混淆集中在相邻或转录组相近脑区，例如 {confusion}。"
            if similar
            else f"当前主要混淆方向包括 {confusion}，需要排查是否存在不合理的跨系统误判。"
        )
        return {
            "what": "混淆矩阵展示真实脑区与预测脑区之间的对应关系。",
            "how": "对角线越集中越好；非对角线亮块表示最常见的误判方向。",
            "good": "理想情况是亮度集中在对角线。相邻皮层区域少量混淆通常比跨系统混淆更容易接受。",
            "interpretation": interpretation,
        }
    if module == "rank":
        return {
            "what": "Rank distribution 展示真实脑区在候选列表中的排序位置。",
            "how": "如果大量样本集中在 rank 1-3，说明模型能把真实来源排到较前位置。",
            "good": "大多数样本位于 rank 1 代表精细定位较好；大量位于 rank 1-3 代表候选筛查能力较强。",
            "interpretation": "如果长尾样本较多，建议回查这些样本的 overlap、QC 和标签是否可靠。",
        }
    if module == "roc":
        return {
            "what": "ROC / AUC 衡量每个脑区与其他脑区的可分离程度。",
            "how": "曲线越靠左上、AUC 越接近 1，说明该脑区越容易被正确区分。",
            "good": "AUC >= 0.80 通常较好；0.65-0.80 为中等；低于 0.65 表示区分度有限。",
            "interpretation": f"当前宏平均 AUC 为 {_fmt(m.get('auc'))}，整体区分度属于{_level_cn(levels.get('auc', 'unknown'))}水平。",
        }
    if module == "confidence_margin":
        return {
            "what": "Confidence 反映模型对 Top1 结果的相对把握，margin 反映 Top1 与 Top2 的分数差距。",
            "how": "若正确样本通常有更高 confidence 和更大 margin，说明这些分数适合做结果可信度分层。",
            "good": "平均 confidence 较高且 margin 明显时，通常说明第一候选领先更充分。",
            "interpretation": f"当前平均 confidence 为 {_fmt(m.get('confidence'))}，平均 margin 为 {_fmt(m.get('margin'))}，margin 水平属于{_level_cn(levels.get('margin', 'unknown'))}。",
        }
    if module == "stability":
        return {
            "what": "Bootstrap stability 表示重采样后 Top1 是否保持一致。",
            "how": "稳定性越高，说明结果越不依赖少数基因；低稳定性提示结论需要更谨慎解读。",
            "good": "平均 stability >= 0.80 通常较好；0.50-0.80 为中等；低于 0.50 说明结果偏脆弱。",
            "interpretation": f"当前平均稳定性为 {_fmt(m.get('stability'))}，属于{_level_cn(levels.get('stability', 'unknown'))}水平。",
        }
    if module == "failure":
        return {
            "what": "Failure mode 将错误、弃权、低 overlap、低 margin、低 confidence 等问题拆开来看。",
            "how": "重点看占比最高的问题类型，从而判断瓶颈更接近数据质量、证据不足还是模型区分困难。",
            "good": "理想状态是 correct Top1 占主导，其他问题类型分散且比例较低。",
            "interpretation": "当前存在较集中的错误模式，建议优先处理主导型问题来源。" if failure_high else "当前未见特别强的单一失败模式，可结合混淆矩阵继续定位局部问题。",
        }
    return {"what": "", "how": "", "good": "", "interpretation": ""}


def reviewer_summary(insights: Dict[str, Any]) -> Dict[str, str]:
    m = insights.get("metrics", {})
    levels = insights.get("levels", {})
    k = int(insights.get("k", 3) or 3)
    topk_gap = insights.get("topk_gap")
    confusion = insights.get("major_confusion_text", "暂无明显集中混淆")
    similar = insights.get("similar_confusions", False)

    if topk_gap is not None and topk_gap >= THRESHOLDS.topk_gap_good:
        performance = (
            f"在当前参考图谱与参数设置下，模型的 Top1 accuracy 为 {_fmt(m.get('top1'))}，"
            f"Top{k} accuracy 为 {_fmt(m.get('topk'))}。Top{k} 高于 Top1，提示模型通常能够将真实来源缩小至较小候选范围，"
            "但第一候选的精细定位仍有进一步提升空间。"
        )
    else:
        performance = (
            f"在当前参考图谱与参数设置下，模型的 Top1 accuracy 为 {_fmt(m.get('top1'))}，"
            f"Top{k} accuracy 为 {_fmt(m.get('topk'))}，整体精细定位能力处于{_level_cn(levels.get('top1', 'unknown'))}水平。"
        )

    if similar:
        error = (
            f"误差主要集中于相邻或转录组特征相近的脑区之间，例如 {confusion}。"
            "这类混淆更可能反映局部区域之间的生物学相似性，而非完全随机误判。"
        )
    else:
        error = (
            f"主要误差方向包括 {confusion}。若这些混淆跨越较远脑区或不同系统，"
            "则需要进一步检查标签一致性、参考图谱层级以及 signature 选择。"
        )

    stability = (
        f"Bootstrap 分析显示平均稳定性为 {_fmt(m.get('stability'))}，"
        f"平均 confidence 为 {_fmt(m.get('confidence'))}，平均 decision margin 为 {_fmt(m.get('margin'))}。"
        f"整体稳定性处于{_level_cn(levels.get('stability', 'unknown'))}水平，"
        "支持将稳定性与 margin 作为结果复核和可信度分层的重要依据。"
    )
    return {"performance": performance, "error": error, "stability": stability}
