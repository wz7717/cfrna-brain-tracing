from __future__ import annotations

import io
import json as _json
import traceback
from typing import Dict

import numpy as np
import pandas as pd
import plotly.express as px
import streamlit as st

from app.components.layout import render_kpi_cards, render_panel_header, render_section_band
from app.database_mode import database_label, get_database_mode
from app.i18n import tr
from app.shared import DB_PATH, init_processor, render_page_hero, render_result_hint
from core.methods import METHOD_SPECS, method_choices, method_help_markdown, method_label
from data.dao import get_atlas_options


def _download_df_button(df: pd.DataFrame, label: str, filename: str) -> None:
    if df is None or df.empty:
        return
    buf = io.BytesIO()
    df.to_csv(buf, index=False, encoding="utf-8-sig")
    buf.seek(0)
    st.download_button(label, buf.getvalue(), filename, "text/csv")


def _safe_float(value):
    try:
        if value is None or pd.isna(value):
            return None
        return float(value)
    except Exception:
        return None


def _suite_summary(detail_df: pd.DataFrame, metrics_df: pd.DataFrame, k: int) -> dict:
    summary = detail_df.attrs.get("summary", {}) if hasattr(detail_df, "attrs") else {}
    metric_map = {}
    if metrics_df is not None and not metrics_df.empty and {"metric", "value"}.issubset(metrics_df.columns):
        metric_map = dict(zip(metrics_df["metric"].astype(str), metrics_df["value"]))
    top1 = _safe_float(summary.get("top1_acc")) or _safe_float(metric_map.get("Top1_acc_valid"))
    topk = _safe_float(summary.get(f"top{k}_acc")) or _safe_float(metric_map.get(f"Top{k}_acc_valid"))
    auc = _safe_float(summary.get("auc")) or _safe_float(metric_map.get("MacroAUC_ovr_valid"))
    stability = _safe_float(summary.get("mean_top1_stability")) or _safe_float(metric_map.get("Mean_top1_stability_valid"))
    confidence = _safe_float(summary.get("mean_top1_confidence")) or _safe_float(metric_map.get("Mean_top1_confidence_valid"))
    margin = _safe_float(summary.get("mean_decision_margin")) or _safe_float(metric_map.get("Mean_decision_margin_valid"))
    abstain = _safe_float(summary.get("abstain_rate")) or _safe_float(metric_map.get("Abstain_rate"))
    return {
        "top1": top1,
        "topk": topk,
        "auc": auc,
        "stability": stability,
        "confidence": confidence,
        "margin": margin,
        "abstain": abstain,
    }


def _top_confusion_text(confusion_long_df: pd.DataFrame) -> str:
    if confusion_long_df is None or confusion_long_df.empty:
        return tr("暂未发现集中混淆。", "No concentrated confusion pattern is currently visible.")
    wrong = confusion_long_df[confusion_long_df["truth_region"].astype(str) != confusion_long_df["pred_region"].astype(str)].copy()
    if wrong.empty:
        return tr("主要预测集中在对角线，未见明显错误方向。", "Predictions are largely concentrated on the diagonal with no obvious dominant error direction.")
    sort_cols = [c for c in ["count", "row_fraction"] if c in wrong.columns]
    wrong = wrong.sort_values(sort_cols, ascending=False).head(3)
    pairs = [f"{r.truth_region} -> {r.pred_region}" for r in wrong.itertuples()]
    return " ; ".join(pairs)


def _rule_summary_text(summary: dict, confusion_text: str, k: int) -> dict:
    top1 = summary.get("top1")
    topk = summary.get("topk")
    stability = summary.get("stability")
    auc = summary.get("auc")
    gap = None if top1 is None or topk is None else topk - top1

    if top1 is not None and top1 >= 0.75:
        one_liner = tr(
            "当前模型已具备较好的 Top1 脑区识别能力。",
            "The current model already shows solid Top1 brain-region identification performance.",
        )
    elif gap is not None and gap >= 0.15:
        one_liner = tr(
            f"当前模型更擅长把真实脑区缩小到 Top{k} 候选范围，但精细区分仍有提升空间。",
            f"The current model is better at narrowing the true region into the Top{k} candidate set than making the final fine-grained distinction.",
        )
    else:
        one_liner = tr(
            "当前模型具备一定识别能力，但精细脑区定位仍需进一步优化。",
            "The current model provides usable discrimination, but fine-grained brain-region localization still needs optimization.",
        )

    main_problem = tr(
        f"当前主要混淆方向包括：{confusion_text}",
        f"The main current confusion pattern includes: {confusion_text}",
    )

    if stability is not None and stability < 0.5:
        next_step = tr(
            "建议优先提高稳定性，例如增加更稳健的 signature、检查 overlap 基因数，或降低对弱证据样本的解释强度。",
            "Prioritize improving stability by strengthening the signature set, checking overlap genes, or being more conservative on weak-evidence samples.",
        )
    elif auc is not None and auc < 0.65:
        next_step = tr(
            "建议优先检查 atlas / signature 是否匹配当前样本类型，并重点排查标签一致性。",
            "Prioritize checking whether atlas/signature choices match the current sample type, and review label consistency.",
        )
    else:
        next_step = tr(
            "建议优先针对主要混淆脑区优化 marker 权重或引入更细的区域参考层。",
            "Prioritize refining marker weights for the dominant confusion pairs or introducing a finer regional reference layer.",
        )
    return {"one_liner": one_liner, "main_problem": main_problem, "next_step": next_step}


def _explanation_box(title: str, what_zh: str, what_en: str, how_zh: str, how_en: str, good_zh: str, good_en: str, interp_zh: str, interp_en: str) -> None:
    st.info(
        f"""
**{title}**

**{tr("这是什么", "What")}**  
{tr(what_zh, what_en)}

**{tr("怎么看", "How to read")}**  
{tr(how_zh, how_en)}

**{tr("什么结果算好", "What is good")}**  
{tr(good_zh, good_en)}

**{tr("当前结果说明了什么", "Interpretation")}**  
{tr(interp_zh, interp_en)}
        """
    )


def _render_qc_overview(processor) -> None:
    qc_overview_df = processor.compute_database_cohort_qc()
    st.markdown(f'<div class="result-zone">{tr("结果区：评估前样本质控概览", "Result zone: pre-benchmark cohort QC overview")}</div>', unsafe_allow_html=True)
    render_result_hint(
        tr(
            "建议先看样本总体风险结构，再决定是否仅在 Low risk 样本上运行 Benchmark。",
            "Review the cohort risk structure first, then decide whether Benchmark should be run only on Low-risk samples.",
        )
    )
    if qc_overview_df.empty:
        st.info(tr("当前数据库中没有可用于 cohort QC 校准的样本表达矩阵。", "No cohort-calibratable sample matrices are currently available in the database."))
        return
    total_n = int(len(qc_overview_df))
    low_n = int((qc_overview_df["overall_risk"] == "Low risk").sum())
    moderate_n = int((qc_overview_df["overall_risk"] == "Moderate risk").sum())
    high_n = int((qc_overview_df["overall_risk"] == "High risk").sum())
    uncal_n = int((qc_overview_df["overall_risk"] == "Uncalibrated").sum())
    render_kpi_cards(
        [
            {"icon": "ALL", "label": tr("样本总数", "Total samples"), "value": total_n, "note": tr("参与 cohort QC 评估的样本", "Samples included in cohort QC")},
            {"icon": "LOW", "label": "Low risk", "value": low_n, "note": tr("总体风险较低", "Lower overall risk")},
            {"icon": "MID", "label": "Moderate risk", "value": moderate_n, "note": tr("建议重点复核", "Recommended for closer review")},
            {"icon": "HIGH", "label": "High risk", "value": high_n, "note": tr("解释需更谨慎", "Interpretation should be cautious")},
            {"icon": "UNC", "label": "Uncalibrated", "value": uncal_n, "note": tr("尚缺队列校准", "Still lacks cohort calibration")},
        ]
    )
    focus_only = st.checkbox(tr("只看 Moderate / High risk 样本", "Only show Moderate / High risk samples"), value=False, key="benchmark_qc_focus_only")
    qc_view = qc_overview_df.copy()
    if focus_only:
        qc_view = qc_view[qc_view["overall_risk"].isin(["Moderate risk", "High risk"])].copy()
    st.dataframe(
        qc_view[
            [
                "sample_id",
                "subject_id",
                "overall_risk",
                "rbc_risk",
                "immune_risk",
                "brain_risk",
                "rbc_percentile",
                "immune_percentile",
                "brain_percentile",
                "interpretation",
            ]
        ],
        use_container_width=True,
        hide_index=True,
    )


def display_benchmark_page() -> None:
    db_mode = get_database_mode()
    render_page_hero(
        tr(f"{database_label(db_mode)} - Benchmark 性能中心", f"{database_label(db_mode)} - Benchmark Performance Center"),
        tr(
            "在带标签样本上评估溯源模型的准确性、混淆模式、稳定性和可发表级输出，并让每张图都具备直接可解释性。",
            "Evaluate tracing accuracy, confusion patterns, stability and publication-grade benchmarking outputs on labeled samples with an explanation-first interface.",
        ),
        eyebrow="Benchmark",
        pills=[tr("性能摘要", "Performance summary"), tr("图旁解释", "Explanation-first figures"), tr("论文级输出", "Paper-ready outputs")],
    )
    processor = init_processor()
    atlas_opts = get_atlas_options(DB_PATH, species_mode=db_mode)
    if not atlas_opts:
        st.info(
            tr(
                "当前数据库模式下没有可用 atlas，因此 Benchmark 页面先保持同构布局并显示空态提示。导入对应数据库的 atlas 后即可直接复用同一套评估界面。",
                "No atlas is available in the current database mode. The Benchmark page keeps the same layout but shows an empty-state hint until a matching atlas is imported.",
            )
        )
        return
    _render_qc_overview(processor)

    st.markdown(f'<div class="parameter-zone">{tr("参数区：方法、atlas、signature 与评估设置", "Parameter zone: method, atlas, signature and evaluation settings")}</div>', unsafe_allow_html=True)
    with st.expander(tr("这些 benchmark 指标分别是什么意思？", "What do these benchmark metrics mean?"), expanded=False):
        st.markdown(method_help_markdown(method_choices()))
        guide = pd.DataFrame(
            [
                {tr("指标", "Metric"): "Top1 accuracy", tr("说明", "Description"): tr("真实脑区排在第 1 位的比例。", "Fraction of samples whose true region is ranked first.")},
                {tr("指标", "Metric"): "TopK accuracy", tr("说明", "Description"): tr("真实脑区是否进入前 K 个候选。", "Whether the true region falls within the top-K candidates.")},
                {tr("指标", "Metric"): "Rank", tr("说明", "Description"): tr("真实脑区在候选列表中的排序位置。", "The rank position of the true region in the candidate list.")},
                {tr("指标", "Metric"): "Confusion matrix", tr("说明", "Description"): tr("观察哪些脑区最容易相互混淆。", "Shows which regions are most often confused with each other.")},
                {tr("指标", "Metric"): "ROC / AUC", tr("说明", "Description"): tr("评估脑区区分能力和排序能力。", "Measures discrimination and ranking ability across regions.")},
                {tr("指标", "Metric"): "Confidence / margin", tr("说明", "Description"): tr("反映 Top1 结果的把握程度和领先幅度。", "Reflects confidence in Top1 and its lead over the runner-up.")},
                {tr("指标", "Metric"): "Stability", tr("说明", "Description"): tr("Bootstrap 后 Top1 是否仍保持一致。", "Whether Top1 remains stable after bootstrap resampling.")},
                {tr("指标", "Metric"): "Failure mode", tr("说明", "Description"): tr("拆解主要错误来源，例如低 overlap、低 margin 或错误 Top1。", "Breaks down major failure patterns such as low overlap, low margin or wrong Top1.")},
            ]
        )
        st.dataframe(guide, use_container_width=True, hide_index=True)

    method = st.selectbox(tr("评估方法", "Evaluation method"), method_choices(), format_func=lambda m: f"{m} - {METHOD_SPECS[m].label}", index=0)
    k = st.slider("Top-K", 1, 10, 3)
    limit = st.number_input(tr("最多评估样本数（0 表示不限制）", "Maximum number of samples to evaluate (0 means no limit)"), min_value=0, value=0, step=10)
    label_key = st.text_input(tr("自定义标签字段 key（可选）", "Custom label field key (optional)"), value="")

    st.markdown(f'<div class="form-section">{tr("Atlas / Signature / 输入值", "Atlas / Signature / Input value")}</div>', unsafe_allow_html=True)
    col1, col2, col3 = st.columns(3)
    atlas_id = col1.number_input("atlas_id", min_value=1, value=1, step=1)
    sigset_id_in = col2.text_input("sigset_id", value="")
    use_value = col3.selectbox("use_value", ["log1p", "tpm", "zscore"], index=0)

    st.markdown(f'<div class="form-section">{tr("正则化与 Bootstrap", "Regularization and Bootstrap")}</div>', unsafe_allow_html=True)
    col4, col5, col6 = st.columns(3)
    l2 = col4.number_input("L2", min_value=0.0, value=1e-4, step=1e-4, format="%.6f")
    bootstrap_n = col5.number_input("bootstrap_n", min_value=0, value=50, step=10)
    bootstrap_gene_frac = col6.slider("bootstrap_gene_frac", 0.1, 1.0, 0.7)
    ensemble_alpha = st.slider("ensemble_alpha", 0.0, 1.0, 0.5) if method == "ensemble" else 0.5

    with st.expander(tr("自动调参 / 自动选择权重（ensemble）", "Auto-tuning / automatic weight selection (ensemble)"), expanded=False):
        do_tune = st.checkbox(tr("启用自动选择权重", "Enable automatic weight selection"), value=(method == "ensemble"))
        optimize_metric = st.selectbox(tr("优化目标", "Optimization target"), ["top1_acc", "auc"], index=0)
        alpha_step = st.selectbox(tr("alpha 步长", "alpha step"), [0.05, 0.1, 0.2], index=1)
        tune_l2 = st.checkbox(tr("同时搜索 l2", "Search l2 jointly"), value=False)

    if "benchmark_suite_cache" not in st.session_state:
        st.session_state["benchmark_suite_cache"] = None
    if "benchmark_suite_meta" not in st.session_state:
        st.session_state["benchmark_suite_meta"] = None

    st.markdown(f'<div class="action-zone">{tr("操作区：运行 Benchmark", "Action zone: run Benchmark")}</div>', unsafe_allow_html=True)
    if st.button(tr("运行论文级 Benchmark", "Run paper-grade Benchmark"), type="primary", use_container_width=True):
        with st.spinner(tr("正在运行 Benchmark，请稍候...", "Running Benchmark, please wait...")):
            try:
                from benchmark_runner import auto_tune_ensemble_weights, default_label_extractor, run_paper_grade_benchmark_suite

                if label_key.strip():
                    def extractor(row: Dict):
                        meta = row.get("metadata")
                        if not meta:
                            return None
                        try:
                            obj = _json.loads(meta)
                        except Exception:
                            return None
                        return obj.get(label_key.strip())
                else:
                    extractor = default_label_extractor

                sigset_id = int(sigset_id_in) if sigset_id_in.strip() else None
                eff_l2 = float(l2)
                eff_alpha = float(ensemble_alpha)
                grid_df = None
                best_params = None

                if method == "ensemble" and do_tune:
                    alpha_grid = [round(x, 4) for x in list(np.arange(0.0, 1.0 + 1e-9, float(alpha_step)))]
                    l2_grid = [eff_l2] if not tune_l2 else [0.0, 1e-4, 1e-3, 1e-2]
                    best_params, grid_df = auto_tune_ensemble_weights(
                        db_path=DB_PATH,
                        atlas_id=int(atlas_id),
                        sigset_id=sigset_id,
                        use_value=use_value,
                        l2_grid=l2_grid,
                        alpha_grid=alpha_grid,
                        label_extractor=extractor,
                        limit=(None if int(limit) == 0 else int(limit)),
                        optimize_metric=optimize_metric,
                    )
                    if best_params and "ensemble_alpha" in best_params:
                        eff_alpha = float(best_params["ensemble_alpha"])
                        eff_l2 = float(best_params.get("l2", eff_l2))

                suite = run_paper_grade_benchmark_suite(
                    db_path=DB_PATH,
                    method=method,
                    k=int(k),
                    atlas_id=int(atlas_id),
                    sigset_id=sigset_id,
                    use_value=use_value,
                    l2=eff_l2,
                    ensemble_alpha=eff_alpha,
                    bootstrap_n=int(bootstrap_n),
                    bootstrap_gene_frac=float(bootstrap_gene_frac),
                    label_extractor=extractor,
                    limit=(None if int(limit) == 0 else int(limit)),
                )
                st.session_state["benchmark_suite_cache"] = suite
                st.session_state["benchmark_suite_meta"] = {
                    "method": method,
                    "k": int(k),
                    "atlas_id": int(atlas_id),
                    "sigset_id": sigset_id,
                    "use_value": use_value,
                    "l2": eff_l2,
                    "ensemble_alpha": eff_alpha if method == "ensemble" else None,
                    "bootstrap_n": int(bootstrap_n),
                    "bootstrap_gene_frac": float(bootstrap_gene_frac),
                    "limit": None if int(limit) == 0 else int(limit),
                }
                st.session_state["benchmark_grid_cache"] = grid_df
                st.session_state["benchmark_best_params"] = best_params
            except Exception as exc:
                st.error(tr("Benchmark 运行失败，当前评估未能完成。", "Benchmark execution failed and could not be completed."))
                st.info(f"{tr('原始错误', 'Original error')}: {exc}")
                with st.expander(tr("开发者调试信息", "Developer debug details"), expanded=False):
                    st.code(traceback.format_exc(), language="python")

    suite = st.session_state.get("benchmark_suite_cache")
    meta = st.session_state.get("benchmark_suite_meta") or {}
    grid_df = st.session_state.get("benchmark_grid_cache")
    best_params = st.session_state.get("benchmark_best_params")
    if not isinstance(suite, dict):
        return

    detail_df = suite.get("detail_df", pd.DataFrame())
    metrics_df = suite.get("metrics_df", pd.DataFrame())
    if detail_df is None or detail_df.empty or metrics_df is None or metrics_df.empty:
        st.warning(tr("当前没有可展示的 benchmark 结果。", "There is no benchmark output to display yet."))
        return

    summary = _suite_summary(detail_df, metrics_df, int(meta.get("k", 3) or 3))
    confusion_text = _top_confusion_text(suite.get("confusion_long_df", pd.DataFrame()))
    insight = _rule_summary_text(summary, confusion_text, int(meta.get("k", 3) or 3))

    st.markdown(f'<div class="result-zone">{tr("结果区：Benchmark 总结、图表与自动解释", "Result zone: Benchmark summary, figures and interpretation")}</div>', unsafe_allow_html=True)
    render_result_hint(
        tr(
            "建议先看顶部总结卡片，再依次查看 Top1 / TopK、混淆矩阵、ROC、稳定性和 failure mode，最后再读论文级总结。",
            "Start with the summary cards, then review Top1 / TopK, confusion matrix, ROC, stability and failure mode before reading the paper-style summary.",
        )
    )

    render_kpi_cards(
        [
            {"icon": "T1", "label": "Top1 accuracy", "value": "NA" if summary["top1"] is None else f"{summary['top1']:.3f}", "note": tr("最终脑区定位能力", "Final fine-grained localization ability")},
            {"icon": f"T{meta.get('k', 3)}", "label": f"Top{meta.get('k', 3)} accuracy", "value": "NA" if summary["topk"] is None else f"{summary['topk']:.3f}", "note": tr("候选范围缩小能力", "Candidate-range narrowing ability")},
            {"icon": "AUC", "label": "Macro AUC", "value": "NA" if summary["auc"] is None else f"{summary['auc']:.3f}", "note": tr("整体区分能力", "Overall discrimination ability")},
            {"icon": "STA", "label": tr("平均稳定性", "Mean stability"), "value": "NA" if summary["stability"] is None else f"{summary['stability']:.3f}", "note": tr("Bootstrap 重复性", "Bootstrap repeatability")},
        ]
    )
    st.success(f"{tr('一句话总结', 'One-line summary')}: {insight['one_liner']}")
    left, right = st.columns(2)
    left.warning(f"{tr('当前主要问题', 'Main current issue')}: {insight['main_problem']}")
    right.info(f"{tr('推荐下一步', 'Recommended next step')}: {insight['next_step']}")

    if best_params is not None and isinstance(best_params, dict) and "error" not in best_params:
        st.success(
            tr(
                f"自动调参结果：ensemble_alpha={float(best_params.get('ensemble_alpha', meta.get('ensemble_alpha', 0.5))):.3f}, l2={float(best_params.get('l2', meta.get('l2', 1e-4))):.1e}",
                f"Auto-tuning result: ensemble_alpha={float(best_params.get('ensemble_alpha', meta.get('ensemble_alpha', 0.5))):.3f}, l2={float(best_params.get('l2', meta.get('l2', 1e-4))):.1e}",
            )
        )
    if grid_df is not None and isinstance(grid_df, pd.DataFrame) and not grid_df.empty:
        with st.expander(tr("自动调参网格结果（Top 20）", "Auto-tuning grid results (Top 20)"), expanded=False):
            st.dataframe(grid_df.head(20).replace({np.nan: None}), use_container_width=True, hide_index=True)

    render_section_band(tr("Benchmark 总结", "Benchmark Summary"), tr("核心指标、参数快照和论文级阅读顺序集中展示。", "Core metrics, parameter snapshot and paper-style reading order in one place."))
    left, right = st.columns([0.92, 1.08])
    with left:
        render_panel_header(tr("核心指标表", "Core metrics table"), tr("适合用于汇报或参数组间比较。", "Useful for reporting and side-by-side parameter comparison."))
        st.dataframe(metrics_df.replace({np.nan: None}), use_container_width=True, hide_index=True)
    with right:
        render_panel_header(tr("参数快照", "Parameter snapshot"), tr("保留本次运行的关键设置。", "Captures the key settings used in this run."))
        st.code(_json.dumps(meta, ensure_ascii=False, indent=2), language="json")

    hit_df = pd.DataFrame(
        {
            "metric": ["Top1", f"Top{meta.get('k', 3)}", "Balanced acc", "AUC"],
            "value": [
                detail_df.attrs.get("summary", {}).get("top1_acc", np.nan),
                detail_df.attrs.get("summary", {}).get(f"top{meta.get('k', 3)}_acc", np.nan),
                detail_df.attrs.get("summary", {}).get("balanced_acc", np.nan),
                detail_df.attrs.get("summary", {}).get("auc", np.nan),
            ],
        }
    ).replace({np.nan: 0.0})
    fig = px.bar(hit_df, x="metric", y="value", title=tr("核心性能摘要图", "Publish-grade summary metrics"), color_discrete_sequence=["#2f6df6"])
    fig.update_layout(yaxis=dict(range=[0, 1]), height=360, margin=dict(l=10, r=10, t=60, b=10))
    st.plotly_chart(fig, use_container_width=True)
    _explanation_box(
        "Top1 / TopK",
        "Top1 和 TopK accuracy 分别反映最终定位能力与候选范围缩小能力。",
        "Top1 and TopK accuracy reflect final localization ability and candidate-range narrowing ability.",
        "先看 Top1，再比较 TopK 是否明显更高。",
        "Read Top1 first, then check whether TopK is clearly higher.",
        "Top1 较高，或 TopK 明显高于 Top1，都说明模型具有一定价值。",
        "A high Top1 or a clear TopK > Top1 gap both indicate practical value.",
        insight["one_liner"],
        insight["one_liner"],
    )

    st.markdown("### 1. Confusion Matrix")
    conf_norm = suite.get("confusion_norm_df", pd.DataFrame())
    conf_raw = suite.get("confusion_raw_df", pd.DataFrame())
    t1, t2 = st.tabs([tr("归一化热图", "Normalized heatmap"), tr("原始计数热图", "Raw-count heatmap")])
    with t1:
        if conf_norm is not None and not conf_norm.empty:
            fig = px.imshow(conf_norm, text_auto=".2f", aspect="auto", color_continuous_scale=["#eef4ff", "#9fc0ff", "#2f6df6"])
            fig.update_layout(height=620, margin=dict(l=10, r=10, t=40, b=10))
            st.plotly_chart(fig, use_container_width=True)
    with t2:
        if conf_raw is not None and not conf_raw.empty:
            fig = px.imshow(conf_raw, text_auto=True, aspect="auto", color_continuous_scale=["#eef4ff", "#9fc0ff", "#2f6df6"])
            fig.update_layout(height=620, margin=dict(l=10, r=10, t=40, b=10))
            st.plotly_chart(fig, use_container_width=True)
    _explanation_box(
        "Confusion Matrix",
        "用于观察真实脑区与预测脑区之间的错配方向。",
        "Shows where true and predicted regions disagree.",
        "对角线越集中越好；非对角线亮块提示重点混淆方向。",
        "A stronger diagonal is better; off-diagonal hot spots indicate important confusion directions.",
        "相邻或相似脑区之间少量混淆通常比跨系统混淆更容易接受。",
        "Limited confusion among adjacent or similar regions is usually more acceptable than cross-system confusion.",
        f"当前最值得关注的混淆方向为：{confusion_text}",
        f"The most notable confusion pattern is: {confusion_text}",
    )

    st.markdown("### 2. Rank Distribution")
    prob_df = suite.get("probability_df", pd.DataFrame())
    if prob_df is not None and not prob_df.empty and "label" in prob_df.columns:
        class_cols = [c for c in prob_df.columns if c not in {"sample_id", "label"}]
        rank_rows = []
        for _, row in prob_df.iterrows():
            label = str(row.get("label"))
            if label not in class_cols:
                continue
            ordered = pd.to_numeric(row[class_cols], errors="coerce").fillna(0.0).sort_values(ascending=False)
            rank_rows.append({"true_rank": list(ordered.index).index(label) + 1})
        rank_df = pd.DataFrame(rank_rows)
        if not rank_df.empty:
            counts = rank_df["true_rank"].value_counts().sort_index().reset_index()
            counts.columns = ["true_rank", "n_samples"]
            fig = px.bar(counts, x="true_rank", y="n_samples", text="n_samples", color_discrete_sequence=["#2f6df6"])
            fig.update_layout(height=400, margin=dict(l=10, r=10, t=40, b=10))
            st.plotly_chart(fig, use_container_width=True)
    _explanation_box(
        "Rank Distribution",
        "显示真实脑区在候选排序中的位置分布。",
        "Shows where the true region appears in the ranked candidate list.",
        "大量样本落在 rank 1-3，通常说明模型具备较好的筛选价值。",
        "If many samples fall within rank 1-3, the model usually has useful screening value.",
        "越多样本集中在前几名越好。",
        "More samples concentrated in the top ranks is better.",
        tr("如果长尾样本较多，建议回查这些样本的 overlap、QC 和标签质量。", "A long tail suggests checking overlap, QC and label quality for those samples."),
        tr("如果长尾样本较多，建议回查这些样本的 overlap、QC 和标签质量。", "A long tail suggests checking overlap, QC and label quality for those samples."),
    )

    st.markdown("### 3. ROC / AUC")
    roc_curve_df = suite.get("roc_curve_df", pd.DataFrame())
    if roc_curve_df is not None and not roc_curve_df.empty:
        fig = px.line(roc_curve_df, x="fpr", y="tpr", color="region_id", color_discrete_sequence=px.colors.qualitative.Set2)
        fig.add_shape(type="line", x0=0, y0=0, x1=1, y1=1, line=dict(dash="dash", color="#94a3b8"))
        fig.update_layout(height=560, margin=dict(l=10, r=10, t=40, b=10))
        st.plotly_chart(fig, use_container_width=True)
    roc_summary_df = suite.get("roc_summary_df", pd.DataFrame())
    if roc_summary_df is not None and not roc_summary_df.empty:
        st.dataframe(roc_summary_df.replace({np.nan: None}), use_container_width=True, hide_index=True)
    _explanation_box(
        "ROC / AUC",
        "衡量脑区与其他脑区之间的可区分程度。",
        "Measures how well each region can be separated from all others.",
        "曲线越靠近左上角、AUC 越高，说明排序能力越强。",
        "Curves closer to the top-left and higher AUC indicate stronger ranking ability.",
        "AUC 接近 1 最理想，低于 0.65 通常提示区分度有限。",
        "AUC near 1 is ideal, while values below 0.65 usually indicate limited separation.",
        tr(f"当前 Macro AUC 为 {summary['auc']:.3f}。" if summary["auc"] is not None else "当前 AUC 暂不可评估。", f"Current Macro AUC is {summary['auc']:.3f}." if summary["auc"] is not None else "AUC is not currently evaluable."),
        tr(f"当前 Macro AUC 为 {summary['auc']:.3f}。" if summary["auc"] is not None else "当前 AUC 暂不可评估。", f"Current Macro AUC is {summary['auc']:.3f}." if summary["auc"] is not None else "AUC is not currently evaluable."),
    )

    st.markdown("### 4. Confidence / Margin")
    if detail_df is not None and not detail_df.empty:
        valid = detail_df[detail_df.get("abstained", 0) == 0].copy() if "abstained" in detail_df.columns else detail_df.copy()
        cols = [c for c in ["top1_confidence", "decision_margin"] if c in valid.columns]
        if not valid.empty and cols:
            long_df = valid.melt(id_vars=[c for c in ["sample_id", "label", "hit1"] if c in valid.columns], value_vars=cols, var_name="metric", value_name="value")
            long_df["prediction"] = np.where(long_df.get("hit1", 0).astype(int) == 1, tr("Top1 正确", "Top1 correct"), tr("Top1 错误", "Top1 wrong"))
            fig = px.box(long_df, x="metric", y="value", color="prediction", points="all", color_discrete_sequence=["#1f9d75", "#d43f56"])
            fig.update_layout(height=440, margin=dict(l=10, r=10, t=40, b=10))
            st.plotly_chart(fig, use_container_width=True)
    _explanation_box(
        "Confidence / Margin",
        "分别反映 Top1 结果的把握程度和 Top1 对 Top2 的领先幅度。",
        "These reflect confidence in Top1 and the lead of Top1 over Top2.",
        "正确样本如果整体具有更高 confidence 和更大 margin，说明这些指标有解释价值。",
        "If correct samples show higher confidence and larger margins overall, these metrics are informative.",
        "高 confidence 且 margin 明显时，通常说明第一候选更可靠。",
        "High confidence with a clear margin usually means the first candidate is more reliable.",
        tr("可以把 confidence 和 margin 作为结果分层和人工复核优先级的辅助依据。", "Confidence and margin can help stratify results and prioritize manual review."),
        tr("可以把 confidence 和 margin 作为结果分层和人工复核优先级的辅助依据。", "Confidence and margin can help stratify results and prioritize manual review."),
    )

    st.markdown("### 5. Bootstrap Stability")
    stability_bin_df = suite.get("stability_bin_df", pd.DataFrame())
    if detail_df is not None and not detail_df.empty and "top1_stability" in detail_df.columns:
        valid = detail_df[detail_df.get("abstained", 0) == 0].copy() if "abstained" in detail_df.columns else detail_df.copy()
        if not valid.empty and not valid["top1_stability"].isna().all():
            fig = px.scatter(valid, x="top1_stability", y="top1_confidence", color="label", symbol="hit1")
            fig.update_layout(height=520, margin=dict(l=10, r=10, t=40, b=10))
            st.plotly_chart(fig, use_container_width=True)
        if stability_bin_df is not None and not stability_bin_df.empty:
            fig2 = px.bar(stability_bin_df, x="stability_bin", y="top1_acc", text="n_samples", color_discrete_sequence=["#74a1ff"])
            fig2.update_layout(yaxis=dict(range=[0, 1]), height=360, margin=dict(l=10, r=10, t=40, b=10))
            st.plotly_chart(fig2, use_container_width=True)
    _explanation_box(
        "Bootstrap Stability",
        "衡量重复抽样后 Top1 是否仍保持一致。",
        "Measures whether Top1 remains consistent after resampling.",
        "稳定性越高，说明结果越不依赖少数基因。",
        "Higher stability means the result depends less on a few specific genes.",
        "稳定性较高时，模型结论通常更适合写入正式结果。",
        "Higher stability usually makes the conclusion more suitable for formal reporting.",
        tr(f"当前平均稳定性为 {summary['stability']:.3f}。" if summary["stability"] is not None else "当前稳定性指标暂不可评估。", f"Current mean stability is {summary['stability']:.3f}." if summary["stability"] is not None else "Stability is not currently evaluable."),
        tr(f"当前平均稳定性为 {summary['stability']:.3f}。" if summary["stability"] is not None else "当前稳定性指标暂不可评估。", f"Current mean stability is {summary['stability']:.3f}." if summary["stability"] is not None else "Stability is not currently evaluable."),
    )

    st.markdown("### 6. Failure Mode")
    failure_rows = []
    if detail_df is not None and not detail_df.empty:
        for _, row in detail_df.iterrows():
            if int(row.get("abstained", 0) or 0) == 1:
                mode = tr("弃权 / 证据不足", "abstained / evidence insufficient")
            elif int(row.get("hit1", 0) or 0) == 1:
                mode = tr("Top1 正确", "correct Top1")
            elif _safe_float(row.get("decision_margin")) is not None and float(row.get("decision_margin")) < 0.05:
                mode = tr("错误且 margin 很低", "wrong with low margin")
            elif _safe_float(row.get("top1_confidence")) is not None and float(row.get("top1_confidence")) < 0.45:
                mode = tr("错误且 confidence 很低", "wrong with low confidence")
            else:
                mode = tr("错误但证据仍可用", "wrong despite usable evidence")
            failure_rows.append({"failure_mode": mode})
    failure_df = pd.DataFrame(failure_rows)
    if not failure_df.empty:
        failure_df = failure_df.value_counts("failure_mode").reset_index(name="n_samples")
        failure_df["fraction"] = failure_df["n_samples"] / max(int(failure_df["n_samples"].sum()), 1)
        fig = px.bar(failure_df, x="failure_mode", y="fraction", text="n_samples", color_discrete_sequence=["#7a9ef8"])
        fig.update_layout(height=420, margin=dict(l=10, r=10, t=40, b=10), xaxis_tickangle=-18, yaxis=dict(range=[0, 1]))
        st.plotly_chart(fig, use_container_width=True)
        st.dataframe(failure_df, use_container_width=True, hide_index=True)
    _explanation_box(
        "Failure Mode",
        "把错误类型拆解开，帮助快速定位问题更像样本质量、证据不足，还是模型本身区分困难。",
        "Breaks down error types to show whether the issue is sample quality, insufficient evidence or model discrimination difficulty.",
        "先看占比最高的失败模式，再回到上面的图解释它为什么出现。",
        "Read the dominant failure category first, then interpret it using the plots above.",
        "理想状态是正确 Top1 占主导，其余错误类型分散且比例较低。",
        "Ideally, correct Top1 dominates and other failure classes remain scattered and small.",
        tr("如果某一种失败模式占比过高，通常说明当前瓶颈已经比较集中。", "If one failure category dominates, the current bottleneck is likely concentrated rather than random."),
        tr("如果某一种失败模式占比过高，通常说明当前瓶颈已经比较集中。", "If one failure category dominates, the current bottleneck is likely concentrated rather than random."),
    )

    st.markdown(f"### {tr('样本级明细', 'Sample-level details')}")
    st.dataframe(detail_df.replace({np.nan: None}), use_container_width=True, hide_index=True)

    render_section_band(tr("适合写入论文结果的总结", "Paper-style result summary"), tr("可直接作为论文 Results 或补充说明的起点。", "Can serve as a starting point for manuscript Results or supplementary text."))
    st.markdown(f"**{tr('性能总结', 'Performance summary')}**")
    st.write(insight["one_liner"])
    st.markdown(f"**{tr('错误模式总结', 'Error-pattern summary')}**")
    st.write(insight["main_problem"])
    st.markdown(f"**{tr('稳定性与可信度总结', 'Stability and confidence summary')}**")
    stability_text_zh = (
        f"当前平均稳定性为 {summary['stability']:.3f}。"
        if summary["stability"] is not None
        else "当前稳定性暂不可评估。"
    ) + " 结合 confidence 与 margin，可进一步判断哪些样本更适合被写入正式结论。"
    stability_text_en = (
        f"Current mean stability is {summary['stability']:.3f}."
        if summary["stability"] is not None
        else "Stability is not currently evaluable."
    ) + " Confidence and margin can help decide which samples are most suitable for formal interpretation."
    st.write(
        tr(stability_text_zh, stability_text_en)
    )

    st.markdown(f'<div class="export-zone">{tr("导出区：下载表格、JSON、图包和 PDF 报告", "Export zone: download tables, JSON, figure bundle and PDF report")}</div>', unsafe_allow_html=True)
    c1, c2, c3 = st.columns(3)
    with c1:
        _download_df_button(detail_df, tr("下载样本级明细 CSV", "Download detail CSV"), f"benchmark_detail_{meta.get('method', 'method')}_top{meta.get('k', 3)}.csv")
        _download_df_button(suite.get("confusion_long_df", pd.DataFrame()), tr("下载混淆矩阵长表", "Download confusion long table"), f"benchmark_confusion_long_{meta.get('method', 'method')}.csv")
        _download_df_button(suite.get("roc_curve_df", pd.DataFrame()), tr("下载 ROC 曲线点表", "Download ROC curve points"), f"benchmark_roc_curve_{meta.get('method', 'method')}.csv")
    with c2:
        _download_df_button(metrics_df, tr("下载指标表 CSV", "Download metrics CSV"), f"benchmark_metrics_{meta.get('method', 'method')}.csv")
        _download_df_button(roc_summary_df, tr("下载 ROC 摘要表", "Download ROC summary"), f"benchmark_roc_summary_{meta.get('method', 'method')}.csv")
        _download_df_button(suite.get("stability_region_df", pd.DataFrame()), tr("下载脑区稳定性表", "Download region stability table"), f"benchmark_stability_region_{meta.get('method', 'method')}.csv")
    with c3:
        _download_df_button(stability_bin_df, tr("下载稳定性分层表", "Download stability-bin table"), f"benchmark_stability_bins_{meta.get('method', 'method')}.csv")
        _download_df_button(suite.get("probability_df", pd.DataFrame()), tr("下载概率矩阵", "Download probability matrix"), f"benchmark_probabilities_{meta.get('method', 'method')}.csv")
        st.download_button(
            tr("下载 Benchmark JSON summary", "Download Benchmark JSON summary"),
            _json.dumps({"parameter_snapshot": meta, "summary": summary, "confusion": confusion_text}, ensure_ascii=False, indent=2),
            file_name=f"benchmark_summary_{meta.get('method', 'method')}_top{meta.get('k', 3)}.json",
            mime="application/json",
        )

    try:
        from reporting import build_benchmark_report_bundle_bytes

        bundle = build_benchmark_report_bundle_bytes(
            suite=suite,
            metadata=meta,
            prefix=f"benchmark_{meta.get('method', 'method')}",
        )
        st.download_button(
            tr("下载 Figure1-Figure6 + PDF 图包", "Download Figure1-Figure6 + PDF bundle"),
            bundle,
            file_name=f"benchmark_figure_report_export_{meta.get('method', 'method')}.zip",
            mime="application/zip",
        )
    except Exception as exc:
        st.warning(f"{tr('图包 / PDF 自动导出失败', 'Figure bundle / PDF export failed')}: {exc}")
