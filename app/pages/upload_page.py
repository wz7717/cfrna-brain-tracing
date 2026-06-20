from __future__ import annotations

import pandas as pd
import plotly.express as px
import streamlit as st

from app.components.layout import render_kpi_cards, render_panel_header, render_section_band
from app.database_mode import database_label, default_species, get_database_mode
from app.i18n import tr
from app.shared import init_processor, render_page_hero, render_result_hint
from data.qc import compute_sample_qc


SPECIES_OPTIONS = ["Macaca mulatta", "Macaca fascicularis", "Homo sapiens"]
UPLOAD_TEXT_DTYPES = {
    key: "string"
    for key in [
        "sample_id",
        "sample",
        "sampleid",
        "subject_id",
        "subject",
        "subjectid",
        "animal_id",
        "ground_truth_region",
        "source_region",
        "injury_region",
        "label_region",
        "true_source",
        "ground_truth_region_name",
    ]
}


def _species_options_for_mode(db_mode: str) -> list[str]:
    if db_mode == "human":
        return ["Homo sapiens", "Macaca mulatta", "Macaca fascicularis"]
    return SPECIES_OPTIONS


def _upload_mode_defaults(db_mode: str) -> dict:
    if db_mode == "human":
        return {
            "sample_id": "HUM_SAMP0001",
            "subject_id": "HUM_SUB001",
            "age_years": 45.0,
            "diagnosis": "Human brain reference / cfRNA cohort",
            "source_type": "human_brain_transcriptome",
            "surgery_region": "",
            "surgery_side": "Not applicable",
            "post_op_day": 0.0,
            "sample_type": "plasma_cfRNA",
            "library_preparation": "RNA-seq / cfRNA-seq",
            "sequencing_platform": "Illumina NovaSeq",
            "total_reads": 60000000,
            "mapping_rate": 95.0,
        }
    return {
        "sample_id": "SAMP0001",
        "subject_id": "SUB001",
        "age_years": 5.0,
        "diagnosis": "Normal",
        "source_type": "brain_injury",
        "surgery_region": "",
        "surgery_side": "Left",
        "post_op_day": 7.0,
        "sample_type": "plasma",
        "library_preparation": "SMARTer Stranded",
        "sequencing_platform": "Illumina NovaSeq",
        "total_reads": 50000000,
        "mapping_rate": 95.0,
    }


def _safe_select_index(options, value, default=0):
    try:
        return options.index(value)
    except Exception:
        return default


def _read_uploaded_file(uploaded_file):
    if uploaded_file.name.endswith(".csv"):
        return pd.read_csv(uploaded_file, dtype=UPLOAD_TEXT_DTYPES)
    if uploaded_file.name.endswith((".tsv", ".txt")):
        return pd.read_csv(uploaded_file, sep="\t", dtype=UPLOAD_TEXT_DTYPES)
    if uploaded_file.name.endswith(".xlsx"):
        return pd.read_excel(uploaded_file, dtype=UPLOAD_TEXT_DTYPES)
    raise ValueError(tr("暂不支持该文件格式", "Unsupported file format"))


def _render_upload_preview(df: pd.DataFrame, qc: dict) -> None:
    render_section_band(
        tr("提交预览", "Submission Preview"),
        tr(
            "写入数据库前预览表达矩阵结构、元数据抽取结果和上传时 QC。",
            "Preview matrix structure, metadata extraction and upload-time QC before saving.",
        ),
    )
    render_kpi_cards(
        [
            {
                "icon": "GENE",
                "label": tr("检测到的基因", "Detected Genes"),
                "value": f"{int((df['tpm_value'] > 0).sum()):,}",
                "note": tr("表达值大于 0 的基因数", "Genes with non-zero expression"),
            },
            {
                "icon": "EXP",
                "label": tr("表达值中位数", "Median expression"),
                "value": f"{pd.to_numeric(df['tpm_value'], errors='coerce').fillna(0).median():.2f}",
                "note": tr("表达强度中位数", "Median expression intensity"),
            },
            {
                "icon": "RBC",
                "label": tr("RBC 背景评分", "RBC Background Score"),
                "value": "NA" if pd.isna(qc.get("rbc_mrna_score")) else f"{qc['rbc_mrna_score']:.2f}",
                "note": str(qc.get("rbc_mrna_risk", "Unknown")),
            },
            {
                "icon": "IMM",
                "label": tr("免疫背景评分", "Immune Background Score"),
                "value": "NA" if pd.isna(qc.get("immune_mrna_score")) else f"{qc['immune_mrna_score']:.2f}",
                "note": str(qc.get("immune_mrna_risk", "Unknown")),
            },
            {
                "icon": "CNS",
                "label": tr("脑信号评分", "Brain Signal Score"),
                "value": "NA" if pd.isna(qc.get("brain_marker_score")) else f"{qc['brain_marker_score']:.2f}",
                "note": str(qc.get("brain_marker_risk", "Unknown")),
            },
        ]
    )

    col_a, col_b = st.columns([1.05, 0.95])
    with col_a:
        render_panel_header(
            tr("表达分布", "Expression Distribution"),
            tr("查看上传样本的内部表达值分布；推荐输入 raw counts/logCPM，TPM/logTPM 仅作为 fallback。", "Review the internal expression-value distribution; raw counts/logCPM are preferred and TPM/logTPM is fallback only."),
        )
        fig = px.histogram(df, x="tpm_value", nbins=60, color_discrete_sequence=["#2f6df6"])
        fig.update_xaxes(type="log", title=tr("内部表达值（对数坐标）", "Internal expression value (log scale)"))
        fig.update_layout(height=340, margin=dict(l=10, r=10, t=10, b=10), showlegend=False)
        st.plotly_chart(fig, use_container_width=True)
    with col_b:
        render_panel_header(
            tr("高表达基因", "Top Expressed Genes"),
            tr("高表达基因有助于识别血浆背景、污染来源和极端异常值。", "Highest-expression genes are useful for spotting plasma background and outliers."),
        )
        high_expr = df.sort_values("tpm_value", ascending=False).head(20)
        fig2 = px.bar(high_expr, x="tpm_value", y="gene_symbol", orientation="h", color_discrete_sequence=["#74a1ff"])
        fig2.update_layout(height=340, margin=dict(l=10, r=10, t=10, b=10), yaxis={"categoryorder": "total ascending"}, showlegend=False)
        st.plotly_chart(fig2, use_container_width=True)

    render_result_hint(
        tr(
            "建议先看总体风险和 RBC / immune / brain 三类评分，再检查高表达基因是否由血浆背景主导。单样本模式常显示 Uncalibrated，表示仍需要队列参照。",
            "Review the overall risk together with RBC, immune and brain scores first, then inspect whether top expressed genes are dominated by plasma background. Single-sample mode often returns Uncalibrated, which means cohort calibration is still needed.",
        )
    )

    warnings = []
    for key in [
        "interpretation",
        "mir451a_mir23a_ratio_interpretation",
        "rbc_mrna_interpretation",
        "immune_mrna_interpretation",
        "brain_marker_interpretation",
    ]:
        text = str(qc.get(key, "")).strip()
        if text and text not in warnings:
            warnings.append(text)

    for text in warnings[:4]:
        st.info(text)


def display_data_upload():
    db_mode = get_database_mode()
    defaults = _upload_mode_defaults(db_mode)
    species_options = _species_options_for_mode(db_mode)

    if db_mode == "human":
        subtitle = tr(
            "提交 Homo sapiens 血浆 cfRNA 或人脑转录组表达矩阵，核对样本元数据，并在入库前完成上传时 QC。",
            "Submit Homo sapiens plasma cfRNA or human brain transcriptome expression matrices, review metadata, and run upload-time QC before database ingestion.",
        )
        pills = [tr("Homo sapiens", "Homo sapiens"), tr("人脑图谱", "Human brain atlas"), tr("表达矩阵", "Expression matrix"), tr("QC 预览", "QC preview")]
    else:
        subtitle = tr(
            "上传猕猴血浆 cfRNA 表达矩阵，核对元数据完整性，并在写入数据库前完成上传时 QC 预览。",
            "Submit macaque plasma cfRNA expression matrices, review metadata completeness, and perform upload-time QC before saving into the tracing database.",
        )
        pills = [tr("表达矩阵", "Expression matrix"), tr("元数据核对", "Metadata review"), tr("QC 预览", "QC preview"), tr("SQLite 入库", "SQLite-ready ingestion")]

    render_page_hero(
        tr(f"{database_label(db_mode)} - 数据提交门户", f"{database_label(db_mode)} - Data Submission Portal"),
        subtitle,
        eyebrow=tr("提交", "Submission"),
        pills=pills,
    )
    processor = init_processor()

    st.markdown(f'<div class="action-zone">{tr("操作区：上传样本表达矩阵并生成提交预览", "Action zone: upload an expression matrix and generate a submission preview")}</div>', unsafe_allow_html=True)
    uploaded_file = st.file_uploader(
        tr("上传表达矩阵", "Upload expression matrix"),
        type=["csv", "tsv", "txt", "xlsx"],
        help=tr(
            "推荐列：gene_symbol + raw_counts/count/read_count，或 gene_symbol + logCPM。TPM/logTPM 仅作为兼容旧表格的 fallback；不要求用户上传 VSD。",
            "Recommended columns: gene_symbol plus raw_counts/count/read_count, or gene_symbol plus logCPM. TPM/logTPM is accepted only as a legacy fallback; users are not asked to upload VSD.",
        ),
    )

    df = None
    embedded_meta = {}
    is_valid = False
    validation_errors: list[str] = []
    if uploaded_file:
        try:
            df = _read_uploaded_file(uploaded_file)
            embedded_meta = processor.extract_embedded_metadata(df)
            is_valid, validation_errors = processor.validate_expression_data(df)
            if is_valid:
                qc_input = processor.preprocess_expression_data(df, min_tpm=0.0)
                unit = str(qc_input["expression_unit"].iloc[0]) if "expression_unit" in qc_input.columns and len(qc_input) else "unknown"
                if "fallback" in unit.lower():
                    st.warning(
                        tr(
                            f"当前输入被识别为 {unit}。TPM/logTPM 仅用于兼容旧表格，不等同于当前验证路线中的 raw counts/logCPM 输入。",
                            f"Input was detected as {unit}. TPM/logTPM is only a legacy compatibility path and is not equivalent to the raw counts/logCPM route used in current validation.",
                        )
                    )
                _render_upload_preview(qc_input, compute_sample_qc(qc_input))
            else:
                st.error(tr("上传文件未通过校验。", "The uploaded file did not pass validation."))
                for msg in validation_errors:
                    st.warning(msg)
        except Exception as exc:
            st.error(f"{tr('读取上传文件失败', 'Failed to read uploaded file')}: {exc}")

    st.markdown(f'<div class="parameter-zone">{tr("参数区：样本元数据、实验标签与测序信息", "Parameter zone: sample metadata, experimental labels and sequencing information")}</div>', unsafe_allow_html=True)
    with st.form("metadata_form"):
        st.markdown(f'<div class="form-section">{tr("样本身份信息", "Sample Identity")}</div>', unsafe_allow_html=True)
        col1, col2 = st.columns(2)
        with col1:
            sample_id = st.text_input(tr("样本 ID *", "Sample ID *"), value=str(embedded_meta.get("sample_id", defaults["sample_id"])))
            subject_id = st.text_input(tr("个体 ID *", "Subject ID *"), value=str(embedded_meta.get("subject_id", defaults["subject_id"])))
            species_default = embedded_meta.get("species", default_species(db_mode))
            species = st.selectbox(tr("物种 *", "Species *"), species_options, index=_safe_select_index(species_options, species_default))
            age_years = st.number_input(tr("年龄（岁）", "Age (years)"), min_value=0.0, value=float(embedded_meta.get("age_years", defaults["age_years"])))
        with col2:
            sex = st.radio(tr("性别", "Sex"), ["Male", "Female", "Unknown"], index=_safe_select_index(["Male", "Female", "Unknown"], embedded_meta.get("sex", "Unknown"), 2))
            diagnosis = st.text_input(tr("诊断 / 分组", "Diagnosis / Group"), value=str(embedded_meta.get("diagnosis", defaults["diagnosis"])))
            plasma_volume = st.number_input(tr("血浆体积（mL）", "Plasma volume (mL)"), min_value=0.0, value=float(embedded_meta.get("plasma_volume_ml", 2.0)))
            rin_value = st.number_input("RIN", min_value=0.0, max_value=10.0, value=float(embedded_meta.get("rin_value", 8.0)))
        collection_date = st.text_input(tr("采样日期", "Collection date"), value=str(embedded_meta.get("collection_date", "2026-03-12")))

        st.markdown(f'<div class="form-section">{tr("溯源标签", "Tracing Labels")}</div>', unsafe_allow_html=True)
        tag1, tag2, tag3 = st.columns(3)
        with tag1:
            ground_truth_region = st.text_input(tr("真实脑区标签", "Ground-truth region"), value=str(embedded_meta.get("ground_truth_region", "")))
            source_type = st.text_input(tr("来源类型", "Source type"), value=str(embedded_meta.get("source_type", defaults["source_type"])))
        with tag2:
            surgery_region = st.text_input(tr("手术 / 采样脑区", "Surgery / sampling region"), value=str(embedded_meta.get("surgery_region", defaults["surgery_region"])))
            surgery_side = st.text_input(tr("手术侧别", "Surgery side"), value=str(embedded_meta.get("surgery_side", defaults["surgery_side"])))
        with tag3:
            post_op_day = st.number_input(tr("术后天数", "Post-op day"), value=float(embedded_meta.get("post_op_day", defaults["post_op_day"])))
            sample_type = st.text_input(tr("样本类型", "Sample type"), value=str(embedded_meta.get("sample_type", defaults["sample_type"])))

        st.markdown(f'<div class="form-section">{tr("测序元数据", "Sequencing Metadata")}</div>', unsafe_allow_html=True)
        seq1, seq2 = st.columns(2)
        with seq1:
            sequencing_platform = st.text_input(tr("测序平台", "Sequencing platform"), value=str(embedded_meta.get("sequencing_platform", defaults["sequencing_platform"])))
            library_preparation = st.text_input(tr("建库方案", "Library preparation"), value=str(embedded_meta.get("library_preparation", defaults["library_preparation"])))
        with seq2:
            total_reads = st.number_input(tr("总 reads 数", "Total reads"), min_value=0, value=int(embedded_meta.get("total_reads", defaults["total_reads"])))
            mapping_rate = st.number_input(tr("比对率 (%)", "Mapping rate (%)"), min_value=0.0, max_value=100.0, value=float(embedded_meta.get("mapping_rate", defaults["mapping_rate"])))

        st.markdown(f'<div class="action-zone">{tr("操作区：写入样本元数据与表达矩阵", "Action zone: write sample metadata and expression matrix to SQLite")}</div>', unsafe_allow_html=True)
        submit_button = st.form_submit_button(tr("保存样本到 SQLite", "Save sample to SQLite"), type="primary")

    if submit_button:
        if uploaded_file is None or df is None:
            st.error(tr("请先上传表达矩阵。", "Please upload an expression matrix first."))
            return
        if not is_valid:
            st.error(tr("当前矩阵尚未通过校验。", "The current matrix has not passed validation yet."))
            return
        try:
            processed_df = processor.preprocess_expression_data(df)
            metadata = dict(embedded_meta)
            metadata.update(
                {
                    "sample_id": sample_id,
                    "subject_id": subject_id,
                    "species": species,
                    "age_years": age_years,
                    "sex": sex,
                    "diagnosis": diagnosis,
                    "sample_type": sample_type,
                    "plasma_volume_ml": plasma_volume,
                    "collection_date": str(collection_date),
                    "extraction_method": metadata.get("extraction_method", "Unknown"),
                    "rna_concentration_ng_ul": float(metadata.get("rna_concentration_ng_ul", 0.0) or 0.0),
                    "rin_value": rin_value,
                    "library_preparation": library_preparation,
                    "sequencing_platform": sequencing_platform,
                    "total_reads": int(total_reads),
                    "mapped_reads": int(metadata.get("mapped_reads", int(total_reads * mapping_rate / 100))),
                    "mapping_rate": mapping_rate,
                    "ground_truth_region": ground_truth_region,
                    "source_type": source_type,
                    "surgery_region": surgery_region,
                    "surgery_side": surgery_side,
                    "post_op_day": post_op_day,
                    "gene_id_type": embedded_meta.get("gene_id_type", processed_df["gene_id_type"].iloc[0] if len(processed_df) else None),
                    "qc_status": "Pending",
                }
            )
            processor.save_sample_metadata(metadata)
            processor.save_expression_data(sample_id, processed_df)
            report = processor.generate_qc_report(sample_id)

            st.success(tr(f"样本 {sample_id} 已成功写入数据库。", f"Sample {sample_id} has been stored successfully."))
            render_kpi_cards(
                [
                    {"icon": "ID", "label": tr("样本 ID", "Sample ID"), "value": sample_id, "note": tr("提交完成", "Submission completed")},
                    {"icon": "GENE", "label": tr("保存基因数", "Saved Genes"), "value": f"{len(processed_df):,}", "note": tr("过滤后表达记录", "Post-filter expression records")},
                    {"icon": "QC", "label": tr("QC 状态", "QC Status"), "value": str(report.get("status", "Unknown")), "note": tr("当前上传报告状态", "Current upload-time report status")},
                ]
            )
            if report.get("warnings"):
                st.markdown(f'<div class="result-zone">{tr("结果区：QC 报告摘要", "Result zone: QC report summary")}</div>', unsafe_allow_html=True)
                for msg in report["warnings"]:
                    st.info(msg)
        except Exception as exc:
            st.error(tr("样本提交失败，请检查矩阵格式和元数据完整性。", "Sample submission failed. Please check the matrix format and metadata completeness."))
            with st.expander(tr("开发者调试信息", "Developer debug details"), expanded=False):
                st.exception(exc)
