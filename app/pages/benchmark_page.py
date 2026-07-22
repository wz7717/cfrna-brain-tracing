from __future__ import annotations

import io
import json as _json
import sqlite3
import traceback
from typing import Dict

import numpy as np
import pandas as pd
import plotly.express as px
import streamlit as st

from app.components.layout import render_kpi_cards, render_panel_header, render_section_band
from app.database_mode import database_label, get_database_mode
from app.i18n import tr
from app.services.benchmark_summary import build_benchmark_insights, explanation_for, reviewer_summary
from app.shared import DB_PATH, init_processor, render_page_hero
from core.methods import METHOD_SPECS, method_choices, method_help_markdown, method_label
from data.dao import get_atlas_options, table_exists


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


def _count_table_rows(table_name: str) -> int:
    if not table_exists(DB_PATH, table_name):
        return 0
    conn = sqlite3.connect(DB_PATH)
    try:
        cur = conn.cursor()
        cur.execute(f"SELECT COUNT(*) FROM {table_name}")
        return int(cur.fetchone()[0])
    finally:
        conn.close()


def _count_default_labeled_samples(limit: int | None = None) -> int:
    if not table_exists(DB_PATH, "braintrace_samples"):
        return 0
    try:
        from benchmark.label_utils import default_label_extractor
    except Exception:
        return 0
    conn = sqlite3.connect(DB_PATH)
    try:
        df = pd.read_sql_query("SELECT metadata FROM braintrace_samples", conn)
    finally:
        conn.close()
    n = 0
    for row in df.to_dict("records"):
        if default_label_extractor(row):
            n += 1
            if limit is not None and n >= int(limit):
                break
    return n


def _render_public_demo_benchmark_boundary(n_atlas: int, n_sigsets: int, n_labeled: int) -> None:
    render_section_band(
        tr("公开 Demo 边界", "Public Demo Boundary"),
        tr(
            "公开云端 demo 只暴露轻量 SaleemNetworks 模型，不下载、不展示完整 Bo2023 reference，也不保存用户上传数据。",
            "The public cloud demo exposes only the lightweight SaleemNetworks model; it does not download or expose the full Bo2023 reference and does not persist uploaded data.",
        ),
    )
    render_kpi_cards(
        [
            {"icon": "ATL", "label": tr("真实 atlas", "Real atlases"), "value": f"{n_atlas:,}", "note": tr("需要完整数据库", "Requires full database")},
            {"icon": "SIG", "label": tr("Signature sets", "Signature sets"), "value": f"{n_sigsets:,}", "note": tr("Benchmark 前置条件", "Benchmark prerequisite")},
            {"icon": "LBL", "label": tr("带标签样本", "Labeled samples"), "value": f"{n_labeled:,}", "note": tr("用于准确率评估", "Used for accuracy evaluation")},
        ]
    )
    st.info(
        tr(
            "完整 Benchmark 需要真实 atlas_versions、signature_sets 和带 ground truth 标签的样本。当前公开 demo 下已禁用运行按钮；可在 Tracing 页面上传 raw counts/logCPM 运行 Network 级轻量溯源。",
            "Full Benchmark requires real atlas_versions, signature_sets and ground-truth labeled samples. The run button is disabled in the current public demo state; use the Tracing page to upload raw counts/logCPM for lightweight Network-level tracing.",
        )
    )
    try:
        from app.pages.tracing_page import _render_network_description_table

        _render_network_description_table()
    except Exception:
        st.caption(
            tr(
                "公开 demo 的输出边界：仅展示 Bo2023 10 个 SaleemNetworks 粗粒度位置/功能候选，不给出完整 reference matrix 或论文级 Benchmark。",
                "Public demo output boundary: only the 10 coarse Bo2023 SaleemNetworks anatomical-functional candidates are shown, without full reference matrices or paper-grade Benchmark.",
            )
        )


def _render_explanation_box(title: str, module: str, insights: dict) -> None:
    explanation = explanation_for(module, insights)
    st.info(
        f"""
**{title}**

**{tr("What", "What")}**
{explanation.get("what", "")}

**{tr("How to read", "How to read")}**
{explanation.get("how", "")}

**{tr("What is good", "What is good")}**
{explanation.get("good", "")}

**{tr("Interpretation", "Interpretation")}**
{explanation.get("interpretation", "")}
        """
    )


def _render_qc_overview(processor) -> None:
    qc_overview_df = processor.compute_database_cohort_qc()
    st.markdown(f'<div class="result-zone">{tr("结果区：评估前样本质控概览", "Result zone: pre-benchmark cohort QC overview")}</div>', unsafe_allow_html=True)
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
        width="stretch",
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
    atlas_opts = get_atlas_options(DB_PATH, species_mode=db_mode, include_legacy=False)
    n_atlas = len(atlas_opts)
    n_sigsets = _count_table_rows("signature_sets")
    n_labeled = _count_default_labeled_samples()
    if n_atlas == 0 or n_sigsets == 0 or n_labeled == 0:
        _render_public_demo_benchmark_boundary(n_atlas, n_sigsets, n_labeled)
        return
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
        st.dataframe(guide, width="stretch", hide_index=True)

    method = st.selectbox(tr("评估方法", "Evaluation method"), method_choices(), format_func=lambda m: f"{m} - {METHOD_SPECS[m].label}", index=0)
    k = st.slider("Top-K", 1, 10, 3)
    limit = st.number_input(tr("最多评估样本数（0 表示不限制）", "Maximum number of samples to evaluate (0 means no limit)"), min_value=0, value=0, step=10)
    label_key = st.text_input(tr("自定义标签字段 key（可选）", "Custom label field key (optional)"), value="")

    st.markdown(f'<div class="form-section">{tr("Atlas / Signature / 输入值", "Atlas / Signature / Input value")}</div>', unsafe_allow_html=True)
    col1, col2, col3 = st.columns(3)
    atlas_id = col1.number_input("atlas_id", min_value=1, value=1, step=1)
    sigset_id_in = col2.text_input("sigset_id", value="")
    use_value_labels = {
        "log1p": tr("logCPM / log expression 主路线", "logCPM / log-expression primary route"),
        "zscore": tr("logCPM 标准化对照", "standardized log-expression control"),
        "tpm": tr("TPM / logTPM fallback", "TPM / logTPM fallback"),
    }
    use_value = col3.selectbox(
        tr("输入表达尺度", "Input expression scale"),
        ["log1p", "zscore", "tpm"],
        index=0,
        format_func=lambda x: use_value_labels.get(x, x),
    )
    st.caption(
        tr(
            "推荐 Benchmark 使用 raw counts/logCPM 派生的 logCPM 表达路线；TPM/logTPM 仅用于兼容旧数据，结果应按 fallback 解释。",
            "Benchmark should prefer the logCPM route derived from raw counts/logCPM inputs; TPM/logTPM is kept only for legacy compatibility and should be interpreted as fallback.",
        )
    )
    if use_value == "tpm":
        st.warning(
            tr(
                "当前选择的是 TPM/logTPM fallback，不等同于当前验证路线中的 Bo2023 logCPM 路线。",
                "The selected TPM/logTPM fallback is not equivalent to the Bo2023 logCPM route used in the current validation path.",
            )
        )

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
    if st.button(tr("运行论文级 Benchmark", "Run paper-grade Benchmark"), type="primary", width="stretch"):
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

    k_value = int(meta.get("k", 3) or 3)
    insights = build_benchmark_insights(suite, meta=meta, k=k_value)
    summary = insights["metrics"]
    confusion_text = insights["major_confusion_text"]
    reviewer = reviewer_summary(insights)
    insight = {
        "one_liner": insights["one_liner"],
        "main_problem": insights["main_problem"],
        "next_step": insights["next_step"],
    }

    st.markdown(f'<div class="result-zone">{tr("结果区：Benchmark 总结、图表与自动解释", "Result zone: Benchmark summary, figures and interpretation")}</div>', unsafe_allow_html=True)

    render_kpi_cards(
        [
            {"icon": "T1", "label": "Top1 accuracy", "value": "NA" if summary["top1"] is None else f"{summary['top1']:.3f}", "note": tr("所选有真值端点的 Top1 命中率", "Top1 hit rate at the selected labeled endpoint")},
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
            st.dataframe(grid_df.head(20).replace({np.nan: None}), width="stretch", hide_index=True)

    render_section_band(tr("Benchmark 总结", "Benchmark Summary"), tr("核心指标、参数快照和论文级阅读顺序集中展示。", "Core metrics, parameter snapshot and paper-style reading order in one place."))
    left, right = st.columns([0.92, 1.08])
    with left:
        render_panel_header(tr("核心指标表", "Core metrics table"), tr("适合用于汇报或参数组间比较。", "Useful for reporting and side-by-side parameter comparison."))
        st.dataframe(metrics_df.replace({np.nan: None}), width="stretch", hide_index=True)
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
    st.plotly_chart(fig, width="stretch")
    _render_explanation_box("Top1 / TopK", "accuracy", insights)

    st.markdown("### 1. Confusion Matrix")
    conf_norm = suite.get("confusion_norm_df", pd.DataFrame())
    conf_raw = suite.get("confusion_raw_df", pd.DataFrame())
    t1, t2 = st.tabs([tr("归一化热图", "Normalized heatmap"), tr("原始计数热图", "Raw-count heatmap")])
    with t1:
        if conf_norm is not None and not conf_norm.empty:
            fig = px.imshow(conf_norm, text_auto=".2f", aspect="auto", color_continuous_scale=["#eef4ff", "#9fc0ff", "#2f6df6"])
            fig.update_layout(height=620, margin=dict(l=10, r=10, t=40, b=10))
            st.plotly_chart(fig, width="stretch")
    with t2:
        if conf_raw is not None and not conf_raw.empty:
            fig = px.imshow(conf_raw, text_auto=True, aspect="auto", color_continuous_scale=["#eef4ff", "#9fc0ff", "#2f6df6"])
            fig.update_layout(height=620, margin=dict(l=10, r=10, t=40, b=10))
            st.plotly_chart(fig, width="stretch")
    _render_explanation_box("Confusion Matrix", "confusion", insights)

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
            st.plotly_chart(fig, width="stretch")
    _render_explanation_box("Rank Distribution", "rank", insights)

    st.markdown("### 3. ROC / AUC")
    roc_curve_df = suite.get("roc_curve_df", pd.DataFrame())
    if roc_curve_df is not None and not roc_curve_df.empty:
        fig = px.line(roc_curve_df, x="fpr", y="tpr", color="region_id", color_discrete_sequence=px.colors.qualitative.Set2)
        fig.add_shape(type="line", x0=0, y0=0, x1=1, y1=1, line=dict(dash="dash", color="#94a3b8"))
        fig.update_layout(height=560, margin=dict(l=10, r=10, t=40, b=10))
        st.plotly_chart(fig, width="stretch")
    roc_summary_df = suite.get("roc_summary_df", pd.DataFrame())
    if roc_summary_df is not None and not roc_summary_df.empty:
        st.dataframe(roc_summary_df.replace({np.nan: None}), width="stretch", hide_index=True)
    _render_explanation_box("ROC / AUC", "roc", insights)

    st.markdown("### 4. Confidence / Margin")
    if detail_df is not None and not detail_df.empty:
        valid = detail_df[detail_df.get("abstained", 0) == 0].copy() if "abstained" in detail_df.columns else detail_df.copy()
        cols = [c for c in ["top1_confidence", "decision_margin"] if c in valid.columns]
        if not valid.empty and cols:
            long_df = valid.melt(id_vars=[c for c in ["sample_id", "label", "hit1"] if c in valid.columns], value_vars=cols, var_name="metric", value_name="value")
            long_df["prediction"] = np.where(long_df.get("hit1", 0).astype(int) == 1, tr("Top1 正确", "Top1 correct"), tr("Top1 错误", "Top1 wrong"))
            fig = px.box(long_df, x="metric", y="value", color="prediction", points="all", color_discrete_sequence=["#1f9d75", "#d43f56"])
            fig.update_layout(height=440, margin=dict(l=10, r=10, t=40, b=10))
            st.plotly_chart(fig, width="stretch")
    _render_explanation_box("Confidence / Margin", "confidence_margin", insights)

    st.markdown("### 5. Bootstrap Stability")
    stability_bin_df = suite.get("stability_bin_df", pd.DataFrame())
    if detail_df is not None and not detail_df.empty and "top1_stability" in detail_df.columns:
        valid = detail_df[detail_df.get("abstained", 0) == 0].copy() if "abstained" in detail_df.columns else detail_df.copy()
        if not valid.empty and not valid["top1_stability"].isna().all():
            fig = px.scatter(valid, x="top1_stability", y="top1_confidence", color="label", symbol="hit1")
            fig.update_layout(height=520, margin=dict(l=10, r=10, t=40, b=10))
            st.plotly_chart(fig, width="stretch")
        if stability_bin_df is not None and not stability_bin_df.empty:
            fig2 = px.bar(stability_bin_df, x="stability_bin", y="top1_acc", text="n_samples", color_discrete_sequence=["#74a1ff"])
            fig2.update_layout(yaxis=dict(range=[0, 1]), height=360, margin=dict(l=10, r=10, t=40, b=10))
            st.plotly_chart(fig2, width="stretch")
    _render_explanation_box("Bootstrap Stability", "stability", insights)

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
        st.plotly_chart(fig, width="stretch")
        st.dataframe(failure_df, width="stretch", hide_index=True)
    _render_explanation_box("Failure Mode", "failure", insights)

    st.markdown(f"### {tr('样本级明细', 'Sample-level details')}")
    st.dataframe(detail_df.replace({np.nan: None}), width="stretch", hide_index=True)

    render_section_band(tr("适合写入论文结果的总结", "Paper-style result summary"), tr("可直接作为论文 Results 或补充说明的起点。", "Can serve as a starting point for manuscript Results or supplementary text."))
    st.markdown(f"**{tr('性能总结', 'Performance summary')}**")
    st.write(reviewer["performance"])
    st.markdown(f"**{tr('错误模式总结', 'Error-pattern summary')}**")
    st.write(reviewer["error"])
    st.markdown(f"**{tr('稳定性与可信度总结', 'Stability and confidence summary')}**")
    st.write(reviewer["stability"])

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
