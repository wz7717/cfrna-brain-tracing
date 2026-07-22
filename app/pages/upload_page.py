from __future__ import annotations

import pandas as pd
import streamlit as st

from app.components.layout import render_kpi_cards, render_section_band
from app.database_mode import database_label, default_species, get_database_mode
from app.i18n import tr
from app.shared import init_processor, render_page_hero
from core.network_tracing import load_network_model


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
            "subject_id": "",
            "age_years": None,
            "diagnosis": "",
            "source_type": "",
            "surgery_region": "",
            "surgery_side": "",
            "post_op_day": None,
            "sample_type": "plasma_cfRNA",
            "library_preparation": "",
            "sequencing_platform": "",
            "total_reads": None,
            "mapping_rate": None,
        }
    return {
        "sample_id": "SAMP0001",
        "subject_id": "",
        "age_years": None,
        "diagnosis": "",
        "source_type": "",
        "surgery_region": "",
        "surgery_side": "",
        "post_op_day": None,
        "sample_type": "plasma",
        "library_preparation": "",
        "sequencing_platform": "",
        "total_reads": None,
        "mapping_rate": None,
    }


def _safe_select_index(options, value, default=0):
    try:
        return options.index(value)
    except Exception:
        return default


def _read_uploaded_file(uploaded_file):
    name = str(uploaded_file.name).lower()
    if name.endswith(".csv"):
        return pd.read_csv(uploaded_file, dtype=UPLOAD_TEXT_DTYPES)
    if name.endswith((".tsv", ".txt")):
        return pd.read_csv(uploaded_file, sep="\t", dtype=UPLOAD_TEXT_DTYPES)
    if name.endswith(".xlsx"):
        return pd.read_excel(uploaded_file, dtype=UPLOAD_TEXT_DTYPES)
    raise ValueError(tr("暂不支持该文件格式", "Unsupported file format"))


def _input_audit(df: pd.DataFrame) -> dict:
    unit = str(df["expression_unit"].iloc[0]) if "expression_unit" in df.columns and len(df) else "unknown"
    values = pd.to_numeric(df["tpm_value"], errors="coerce").fillna(0.0)
    model_genes = pd.Index(load_network_model()["genes"].astype(str))
    input_genes = pd.Index(df["gene_symbol"].astype(str).str.strip().unique())
    overlap = int(model_genes.intersection(input_genes).size)
    total_reads = None
    if "read_count" in df.columns and unit == "logCPM_from_raw_counts":
        total_reads = float(pd.to_numeric(df["read_count"], errors="coerce").fillna(0.0).clip(lower=0).sum())
    return {
        "unit": unit,
        "valid_genes": int(len(input_genes)),
        "nonzero_genes": int((values > 0).sum()),
        "total_reads": total_reads,
        "model_overlap": overlap,
        "model_genes": int(len(model_genes)),
    }


def _render_upload_preview(df: pd.DataFrame) -> None:
    audit = _input_audit(df)
    unit_labels = {
        "logCPM_from_raw_counts": "raw counts → logCPM",
        "logCPM": "logCPM",
        "logTPM_fallback": "logTPM fallback",
        "TPM_fallback": "TPM fallback",
    }
    unit_label = unit_labels.get(audit["unit"], audit["unit"])
    reads_value = "—" if audit["total_reads"] is None else f"{audit['total_reads']:,.0f}"
    overlap_fraction = audit["model_overlap"] / audit["model_genes"] if audit["model_genes"] else 0.0
    render_section_band(
        tr("输入完整性与模型覆盖检查", "Input Integrity and Model Coverage"),
        tr(
            "写入数据库前核对输入尺度、有效基因和 Bo2023 Network 模型基因覆盖。",
            "Verify input scale, valid genes and Bo2023 Network model-gene coverage before saving.",
        ),
    )
    render_kpi_cards(
        [
            {
                "icon": "TYPE",
                "label": tr("输入识别类型", "Detected input type"),
                "value": unit_label,
            },
            {
                "icon": "GENE",
                "label": tr("有效基因数", "Valid genes"),
                "value": f"{audit['valid_genes']:,}",
            },
            {
                "icon": "NZ",
                "label": tr("非零基因数", "Non-zero genes"),
                "value": f"{audit['nonzero_genes']:,}",
            },
            {
                "icon": "READ",
                "label": tr("总 reads", "Total reads"),
                "value": reads_value,
            },
            {
                "icon": "COV",
                "label": tr("模型基因覆盖率", "Model-gene coverage"),
                "value": f"{audit['model_overlap']}/{audit['model_genes']} ({overlap_fraction:.1%})",
            },
            {
                "icon": "OK",
                "label": tr("文件校验", "File validation"),
                "value": tr("通过", "Passed"),
            },
        ]
    )


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
    st.warning(
        tr(
            "持久化提示：“保存样本到 SQLite”会将表达矩阵和元数据持久写入应用数据库。不要上传可识别个人身份或敏感的信息。如需移除，请在“样本管理”中删除对应样本及关联记录。直接 Tracing 上传不会通过此提交表单持久化。",
            "Persistence notice: Save sample to SQLite stores the expression matrix and metadata in the application database. Do not upload identifiable or sensitive information. To remove stored data, delete the sample and its linked records in Sample Management. A direct Tracing upload is not persisted through this submission form.",
        )
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
    matrix_total_reads = None
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
                matrix_total_reads = _input_audit(qc_input)["total_reads"]
                unit = str(qc_input["expression_unit"].iloc[0]) if "expression_unit" in qc_input.columns and len(qc_input) else "unknown"
                if "fallback" in unit.lower():
                    st.warning(
                        tr(
                            f"当前输入被识别为 {unit}。TPM/logTPM 仅用于兼容旧表格，不等同于当前验证路线中的 raw counts/logCPM 输入。",
                            f"Input was detected as {unit}. TPM/logTPM is only a legacy compatibility path and is not equivalent to the raw counts/logCPM route used in current validation.",
                        )
                    )
                _render_upload_preview(qc_input)
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
            age_years = st.number_input(tr("年龄（岁）", "Age (years)"), min_value=0.0, value=embedded_meta.get("age_years", defaults["age_years"]), placeholder=tr("未提供", "Not provided"))
        with col2:
            sex = st.radio(tr("性别", "Sex"), ["Male", "Female", "Unknown"], index=_safe_select_index(["Male", "Female", "Unknown"], embedded_meta.get("sex", "Unknown"), 2))
            diagnosis = st.text_input(tr("诊断 / 分组", "Diagnosis / Group"), value=str(embedded_meta.get("diagnosis", defaults["diagnosis"])))
            plasma_volume = st.number_input(tr("血浆体积（mL）", "Plasma volume (mL)"), min_value=0.0, value=embedded_meta.get("plasma_volume_ml"), placeholder=tr("未提供", "Not provided"))
            rin_value = st.number_input("RIN", min_value=0.0, max_value=10.0, value=embedded_meta.get("rin_value"), placeholder=tr("未提供", "Not provided"))
        collection_date = st.text_input(tr("采样日期", "Collection date"), value=str(embedded_meta.get("collection_date", "")), placeholder=tr("未提供", "Not provided"))

        st.markdown(f'<div class="form-section">{tr("溯源标签", "Tracing Labels")}</div>', unsafe_allow_html=True)
        tag1, tag2, tag3 = st.columns(3)
        with tag1:
            ground_truth_region = st.text_input(tr("真实脑区标签", "Ground-truth region"), value=str(embedded_meta.get("ground_truth_region", "")))
            source_type = st.text_input(tr("来源类型", "Source type"), value=str(embedded_meta.get("source_type", defaults["source_type"])))
        with tag2:
            surgery_region = st.text_input(tr("手术 / 采样脑区", "Surgery / sampling region"), value=str(embedded_meta.get("surgery_region", defaults["surgery_region"])))
            surgery_side = st.text_input(tr("手术侧别", "Surgery side"), value=str(embedded_meta.get("surgery_side", defaults["surgery_side"])))
        with tag3:
            post_op_day = st.number_input(tr("术后天数", "Post-op day"), min_value=0.0, value=embedded_meta.get("post_op_day", defaults["post_op_day"]), placeholder=tr("未提供", "Not provided"))
            sample_type = st.text_input(tr("样本类型", "Sample type"), value=str(embedded_meta.get("sample_type", defaults["sample_type"])))

        st.markdown(f'<div class="form-section">{tr("测序元数据", "Sequencing Metadata")}</div>', unsafe_allow_html=True)
        seq1, seq2 = st.columns(2)
        with seq1:
            sequencing_platform = st.text_input(tr("测序平台", "Sequencing platform"), value=str(embedded_meta.get("sequencing_platform", defaults["sequencing_platform"])))
            library_preparation = st.text_input(tr("建库方案", "Library preparation"), value=str(embedded_meta.get("library_preparation", defaults["library_preparation"])))
        with seq2:
            reads_default = embedded_meta.get("total_reads", matrix_total_reads)
            total_reads = st.number_input(tr("总 reads 数", "Total reads"), min_value=0, value=None if reads_default is None else int(reads_default), placeholder=tr("raw counts 自动计算", "Calculated from raw counts"))
            mapping_rate = st.number_input(tr("比对率 (%)", "Mapping rate (%)"), min_value=0.0, max_value=100.0, value=embedded_meta.get("mapping_rate", defaults["mapping_rate"]), placeholder=tr("未提供", "Not provided"))

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
            processed_df = processor.preprocess_expression_data(df, min_tpm=0.0)
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
                    "total_reads": None if total_reads is None else int(total_reads),
                    "mapped_reads": metadata.get("mapped_reads") if metadata.get("mapped_reads") is not None else (None if total_reads is None or mapping_rate is None else int(total_reads * mapping_rate / 100)),
                    "mapping_rate": mapping_rate,
                    "ground_truth_region": ground_truth_region,
                    "source_type": source_type,
                    "surgery_region": surgery_region,
                    "surgery_side": surgery_side,
                    "post_op_day": post_op_day,
                    "gene_id_type": embedded_meta.get("gene_id_type", processed_df["gene_id_type"].iloc[0] if len(processed_df) else None),
                    "qc_status": None,
                }
            )
            # Metadata and expression rows are committed together; failures roll back both.
            # Upload-time legacy RBC/immune/brain marker QC remains intentionally disabled.
            processor.save_sample_with_expression(metadata, processed_df)

            st.success(tr(f"样本 {sample_id} 已成功写入数据库。", f"Sample {sample_id} has been stored successfully."))
            render_kpi_cards(
                [
                    {"icon": "ID", "label": tr("样本 ID", "Sample ID"), "value": sample_id, "note": tr("提交完成", "Submission completed")},
                    {"icon": "GENE", "label": tr("保存基因数", "Saved Genes"), "value": f"{len(processed_df):,}", "note": tr("过滤后表达记录", "Post-filter expression records")},
                ]
            )
        except Exception as exc:
            st.error(tr("样本提交失败，请检查矩阵格式和元数据完整性。", "Sample submission failed. Please check the matrix format and metadata completeness."))
            with st.expander(tr("开发者调试信息", "Developer debug details"), expanded=False):
                st.exception(exc)
