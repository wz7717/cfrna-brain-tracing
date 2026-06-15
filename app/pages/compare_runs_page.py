from __future__ import annotations

import pandas as pd
import plotly.express as px
import streamlit as st

from app.components.layout import render_kpi_cards, render_panel_header, render_section_band
from app.database_mode import database_label, get_database_mode, matches_species
from app.i18n import tr
from app.shared import DB_PATH, render_page_hero, render_result_hint
from data.dao import get_run_results, list_runs


def display_run_compare():
    db_mode = get_database_mode()
    render_page_hero(
        tr(f"{database_label(db_mode)} - Run 对比", f"{database_label(db_mode)} - Run Comparison"),
        tr(
            "比较同一样本在不同方法或不同参数设置下的溯源结果，判断生物学结论是否足够稳定。",
            "Compare multiple tracing runs from the same sample across methods or parameter settings, and decide whether the biological conclusion is stable enough to report.",
        ),
        eyebrow=tr("对比", "Comparison"),
        pills=[tr("多 Run 复核", "Multi-run review"), tr("参数敏感性", "Parameter sensitivity"), tr("得分稳定性", "Score stability")],
    )
    runs_df = list_runs(DB_PATH, limit=300)
    try:
        from app.shared import init_processor

        samples_df = init_processor().get_all_samples()
        if not samples_df.empty and "species" in samples_df.columns:
            allowed_ids = set(samples_df[samples_df["species"].apply(lambda x: matches_species(x, db_mode))]["sample_id"].astype(str))
            runs_df = runs_df[runs_df["sample_id"].astype(str).isin(allowed_ids)].copy()
    except Exception:
        pass
    if runs_df.empty:
        st.info(tr("当前数据库模式下没有可对比的分析记录。请先运行至少一次 v2 溯源分析。", "There are no comparable analysis runs in the current database mode. Run at least one v2 tracing analysis first."))
        return

    render_kpi_cards(
        [
            {"icon": "RUN", "label": tr("可用 Run", "Available Runs"), "value": f"{len(runs_df):,}", "note": tr("已加载的近期运行记录", "Recent run records loaded")},
            {"icon": "SMP", "label": tr("覆盖样本", "Samples Covered"), "value": f"{runs_df['sample_id'].astype(str).nunique():,}", "note": tr("具有 run 历史的样本数", "Distinct samples with run history")},
            {"icon": "MTH", "label": tr("方法数", "Methods"), "value": f"{runs_df['method'].astype(str).nunique():,}", "note": tr("当前缓存中的方法种类", "Methods represented in the cache")},
        ]
    )

    st.markdown(f'<div class="parameter-zone">{tr("参数区：筛选样本、选择 Run 并设置对比指标", "Parameter zone: filter sample, choose runs and define comparison metric")}</div>', unsafe_allow_html=True)
    sample_ids = [tr("全部样本", "(All samples)")] + sorted(runs_df["sample_id"].astype(str).unique().tolist())
    filt = st.selectbox(tr("按样本筛选", "Filter by sample"), sample_ids, index=0)
    view_df = runs_df if filt == tr("全部样本", "(All samples)") else runs_df[runs_df["sample_id"].astype(str) == filt]
    options = [f"{r.run_id} | sample:{r.sample_id} | {r.method} | {r.created_at}" for r in view_df.itertuples()]
    chosen = st.multiselect(tr("选择要对比的 Runs（建议 2-5 个）", "Choose runs to compare (recommended: 2-5)"), options, default=options[: min(3, len(options))])
    if not chosen:
        st.caption(tr("至少选择一个 Run 后再开始对比。", "Choose at least one run to begin the comparison."))
        return

    run_ids = [item.split(" | ")[0] for item in chosen]
    result_frames = []
    for rid in run_ids:
        df = get_run_results(DB_PATH, rid)
        if not df.empty:
            tmp = df.copy()
            tmp["run_id"] = rid
            result_frames.append(tmp)

    if not result_frames:
        st.warning(tr("所选 Run 当前没有 analysis_results 记录。", "The selected runs do not currently have `analysis_results` entries."))
        return

    res = pd.concat(result_frames, ignore_index=True)
    topn = st.slider(tr("对比 Top-N 脑区", "Top-N regions to compare"), 3, 30, 10)
    metric_choice = st.selectbox(tr("对比指标", "Comparison metric"), ["score", "fraction", "confidence"], index=0)
    if metric_choice not in res.columns:
        st.info(tr(f"所选 Run 不包含 `{metric_choice}` 字段，请更换指标。", f"The selected runs do not contain `{metric_choice}`. Try another metric."))
        return

    res_top = res[res["rank"] <= topn].copy()
    pivot = res_top.pivot_table(index="region_id", columns="run_id", values=metric_choice, aggfunc="first").fillna(0.0)

    render_section_band(
        tr("对比摘要", "Comparison Summary"),
        tr("先看不同 Run 是否持续指向相似脑区，再进一步判断参数敏感性。", "See whether the same brain regions remain dominant across runs."),
    )
    leading_regions = res_top.sort_values(["run_id", "rank"]).groupby("run_id").first().reset_index()
    render_kpi_cards(
        [
            {"icon": "TOP1", "label": tr("Top1 脑区种类数", "Unique Top1 Regions"), "value": f"{leading_regions['region_id'].astype(str).nunique():,}", "note": tr("值越小通常表示跨 Run 一致性越高", "Smaller values often imply stronger cross-run consistency")},
            {"icon": "TOPN", "label": tr("参与对比的脑区", "Compared Regions"), "value": f"{pivot.index.nunique():,}", "note": tr("当前窗口中保留的脑区数", "Regions retained in the comparison window")},
            {"icon": metric_choice.upper(), "label": tr("当前指标", "Metric"), "value": metric_choice, "note": tr("当前用于可视化的比较轴", "Currently visualized comparison axis")},
        ]
    )

    st.markdown(f'<div class="result-zone">{tr("结果区：Run 对比矩阵与得分图", "Result zone: run-level comparison matrix and score plot")}</div>', unsafe_allow_html=True)
    render_result_hint(
        tr(
            "重点看同一脑区在不同 Run 中的 score、fraction 或 confidence 是否同向变化；如果 Top1 频繁变化，通常说明参数较敏感或证据尚不够集中。",
            "Focus on whether the same region changes consistently across runs in score, fraction or confidence. If Top1 flips frequently, the model may be parameter-sensitive or the evidence may still be diffuse.",
        )
    )

    col_table, col_plot = st.columns([1.0, 1.05])
    with col_table:
        render_panel_header(
            tr("Run 对比矩阵", "Run Comparison Matrix"),
            tr("宽表结构便于直接比较每个脑区在不同 Run 中的数值。", "Wide-format view of brain-region metrics across selected runs."),
        )
        st.dataframe(pivot, use_container_width=True)
    with col_plot:
        render_panel_header(
            tr("得分对比图", "Score Comparison Plot"),
            tr("分组条形图可快速观察同一脑区在不同 Run 中的变化。", "Grouped bar plot for direct run-to-run comparison."),
        )
        plot_df = pivot.reset_index().melt(id_vars=["region_id"], var_name="run_id", value_name=metric_choice)
        fig = px.bar(
            plot_df,
            x=metric_choice,
            y="region_id",
            color="run_id",
            orientation="h",
            barmode="group",
            color_discrete_sequence=["#2f6df6", "#74a1ff", "#7cd6c1", "#d9c4ff", "#f2b447"],
        )
        fig.update_layout(height=min(860, 38 * (plot_df["region_id"].nunique() + 5)), margin=dict(l=10, r=10, t=10, b=10))
        st.plotly_chart(fig, use_container_width=True)
