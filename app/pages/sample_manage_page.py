from __future__ import annotations

import io

import streamlit as st

from app.components.layout import render_kpi_cards, render_panel_header, render_section_band
from app.database_mode import database_label, get_database_mode, matches_species
from app.i18n import tr
from app.shared import init_processor, render_page_hero, render_result_hint


def display_sample_list():
    db_mode = get_database_mode()
    render_page_hero(
        tr(f"{database_label(db_mode)} - 样本管理", f"{database_label(db_mode)} - Sample Management"),
        tr(
            "查看已上传的血浆 cfRNA 样本，检查队列级 QC 校准结果，并在同一数据库工作区内完成样本管理。",
            "Review uploaded plasma cfRNA samples, inspect cohort-level QC calibration, and perform careful sample-level management without leaving the database workspace.",
        ),
        eyebrow=tr("样本", "Samples"),
        pills=[tr("样本注册表", "Sample registry"), tr("队列 QC", "Cohort QC"), tr("删除保护", "Deletion safeguards")],
    )
    processor = init_processor()
    samples_df = processor.get_all_samples()
    if not samples_df.empty and "species" in samples_df.columns:
        samples_df = samples_df[samples_df["species"].apply(lambda x: matches_species(x, db_mode))].copy()
    if len(samples_df) == 0:
        st.info(
            tr(
                "当前数据库模式下没有样本。请先在数据提交页面上传对应数据库的表达矩阵。",
                "No samples are available in the current database mode. Upload expression matrices for this workspace first.",
            )
        )
        return

    render_kpi_cards(
        [
            {"icon": "SMP", "label": tr("样本数", "Stored Samples"), "value": f"{len(samples_df):,}", "note": tr("当前 cfrna_samples 中的记录数", "Rows currently indexed in cfrna_samples")},
            {"icon": "SUB", "label": tr("个体数", "Subjects"), "value": f"{samples_df['subject_id'].astype(str).nunique():,}", "note": tr("去重后的个体标识数", "Distinct subject identifiers")},
            {"icon": "QC", "label": tr("QC 状态数", "QC States"), "value": f"{samples_df['qc_status'].astype(str).nunique():,}", "note": tr("当前存在的 QC 标签种类", "Distinct QC labels currently present")},
        ]
    )

    st.markdown(f'<div class="result-zone">{tr("结果区：当前数据库样本列表", "Result zone: current sample registry")}</div>', unsafe_allow_html=True)
    render_result_hint(
        tr(
            "建议先核对 sample ID、subject ID、collection date 和 QC status；删除样本前，请确认它不再参与后续 tracing 或 benchmark。",
            "Check sample ID, subject ID, collection date and QC status first. Before deletion, make sure the sample is no longer needed for tracing or benchmark workflows.",
        )
    )
    st.dataframe(
        samples_df,
        use_container_width=True,
        column_config={
            "sample_id": st.column_config.TextColumn(tr("样本 ID", "Sample ID")),
            "subject_id": st.column_config.TextColumn(tr("个体 ID", "Subject ID")),
            "species": st.column_config.TextColumn(tr("物种", "Species")),
            "diagnosis": st.column_config.TextColumn(tr("诊断", "Diagnosis")),
            "collection_date": st.column_config.TextColumn(tr("采样日期", "Collection date")),
            "qc_status": st.column_config.TextColumn(tr("QC 状态", "QC status")),
        },
    )

    st.markdown(f'<div class="parameter-zone">{tr("参数区：基于当前数据库样本执行 cohort QC 校准", "Parameter zone: run cohort QC calibration on current database samples")}</div>', unsafe_allow_html=True)
    st.caption(
        tr(
            "这个操作会基于当前数据库中的样本重新计算 RBC、免疫背景和脑信号风险，不修改数据库结构。",
            "This recalculates RBC, immune and brain marker risks from the current in-database sample cohort. It does not modify the database schema.",
        )
    )
    if st.button(tr("运行当前数据库样本的 cohort QC 校准", "Run cohort QC calibration for current database samples"), key="run_database_cohort_qc", use_container_width=True):
        qc_df = processor.compute_database_cohort_qc()
        if qc_df.empty:
            st.info(tr("当前没有可用于 cohort QC 校准的表达矩阵。", "No expression matrices are currently available for cohort QC calibration."))
        else:
            render_section_band(
                tr("队列 QC 结果", "Cohort QC Results"),
                tr("这些风险标记是基于队列分布的，分位数极端的样本值得优先复核。", "Risk flags are distribution-aware; percentile extremes deserve closer review."),
            )
            risk_filter = st.checkbox(tr("只看 Moderate / High risk 样本", "Only show Moderate/High risk samples"), value=False, key="cohort_qc_highrisk_only")
            qc_view = qc_df.copy()
            if risk_filter:
                qc_view = qc_view[qc_view["overall_risk"].isin(["Moderate risk", "High risk"])].copy()
                if qc_view.empty:
                    st.info(tr("当前队列 QC 结果中没有 Moderate / High risk 样本。", "No Moderate/High risk samples are present in the current cohort QC result."))

            st.dataframe(
                qc_view,
                use_container_width=True,
                hide_index=True,
                column_config={
                    "sample_id": st.column_config.TextColumn(tr("样本 ID", "Sample ID")),
                    "subject_id": st.column_config.TextColumn(tr("个体 ID", "Subject ID")),
                    "overall_risk": st.column_config.TextColumn(tr("总体风险", "Overall risk")),
                    "rbc_score": st.column_config.NumberColumn(tr("RBC 评分", "RBC score"), format="%.3f"),
                    "rbc_percentile": st.column_config.NumberColumn(tr("RBC 分位数", "RBC percentile"), format="%.1f"),
                    "rbc_risk": st.column_config.TextColumn(tr("RBC 风险", "RBC risk")),
                    "immune_score": st.column_config.NumberColumn(tr("免疫评分", "Immune score"), format="%.3f"),
                    "immune_percentile": st.column_config.NumberColumn(tr("免疫分位数", "Immune percentile"), format="%.1f"),
                    "immune_risk": st.column_config.TextColumn(tr("免疫风险", "Immune risk")),
                    "brain_score": st.column_config.NumberColumn(tr("脑信号评分", "Brain score"), format="%.3f"),
                    "brain_percentile": st.column_config.NumberColumn(tr("脑信号分位数", "Brain percentile"), format="%.1f"),
                    "brain_risk": st.column_config.TextColumn(tr("脑信号风险", "Brain risk")),
                    "hemolysis_mirna_risk": st.column_config.TextColumn(tr("miRNA 溶血风险", "miRNA hemolysis")),
                    "interpretation": st.column_config.TextColumn(tr("解释", "Interpretation")),
                },
            )

            buf = io.BytesIO()
            qc_view.to_csv(buf, index=False, encoding="utf-8-sig")
            buf.seek(0)
            st.download_button(
                tr("下载 cohort QC CSV", "Download cohort QC CSV"),
                buf.getvalue(),
                file_name="cohort_qc_calibration_results.csv",
                mime="text/csv",
            )

    st.markdown(f'<div class="action-zone">{tr("操作区：查看样本详情或执行安全删除", "Action zone: inspect a sample or perform safeguarded deletion")}</div>', unsafe_allow_html=True)
    selected_sample = st.selectbox(tr("选择样本", "Choose sample"), [""] + samples_df["sample_id"].astype(str).tolist())
    if not selected_sample:
        st.caption(tr("选择一个样本后即可查看或管理。", "Choose a sample to inspect or manage it."))
        return

    col_view, col_delete = st.columns(2)
    with col_view:
        if st.button(tr("查看样本元数据", "View sample metadata"), key="view_selected_sample", use_container_width=True):
            info = processor.get_sample_info(selected_sample)
            render_panel_header(
                tr("样本元数据", "Sample Metadata"),
                tr("用于核对、审计和排错的原始存储元数据。", "Raw stored metadata for validation, auditing and troubleshooting."),
            )
            st.json(info)

    confirm_key = f"confirm_delete_{selected_sample}"
    if confirm_key not in st.session_state:
        st.session_state[confirm_key] = False

    with col_delete:
        if st.button(tr("删除样本", "Delete sample"), key="delete_selected_sample", use_container_width=True):
            st.session_state[confirm_key] = True

    if st.session_state[confirm_key]:
        st.markdown(f'<div class="danger-zone">{tr("危险区：样本删除二次确认", "Danger zone: secondary confirmation for sample deletion")}</div>', unsafe_allow_html=True)
        st.warning(
            tr(
                "这个操作会删除样本及其关联表达和分析记录，无法撤销。",
                "This action will remove the sample and its linked expression / analysis records. It cannot be undone.",
            )
        )
        typed_sample_id = st.text_input(
            tr("再次输入样本 ID 以确认删除", "Type the sample ID again to confirm deletion"),
            key=f"delete_verify_{selected_sample}",
            placeholder=selected_sample,
            help=tr(
                "只有当输入内容与选中样本 ID 完全一致时，确认按钮才会启用。",
                "The confirmation button is enabled only when the typed ID exactly matches the selected sample ID.",
            ),
        )
        delete_enabled = typed_sample_id == selected_sample
        if not delete_enabled:
            st.caption(f"{tr('需要输入的确认字符串', 'Required confirmation string')}: {selected_sample}")
        col_yes, col_no = st.columns(2)
        with col_yes:
            if st.button(tr("确认删除", "Confirm deletion"), key="confirm_delete_selected_sample", use_container_width=True, disabled=not delete_enabled):
                processor.delete_sample(selected_sample)
                st.session_state[confirm_key] = False
                st.success(tr("样本已删除。", "The sample has been deleted."))
                st.rerun()
        with col_no:
            if st.button(tr("取消", "Cancel"), key="cancel_delete_selected_sample", use_container_width=True):
                st.session_state[confirm_key] = False
