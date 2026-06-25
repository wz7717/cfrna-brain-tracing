from __future__ import annotations

import io
import json
import traceback

import numpy as np
import pandas as pd
import plotly.express as px
import streamlit as st

from app.components.layout import render_kpi_cards, render_panel_header, render_section_band
from app.components.plot_panels import make_fraction_ci_bar, make_score_bar, make_stability_bar
from app.components.result_cards import render_primary_metrics, render_run_meta
from app.database_mode import database_label, get_database_mode, matches_species
from app.i18n import tr
from app.shared import DB_PATH, init_processor, init_tracer, render_page_hero, render_result_hint
from core.bo2023_region_tracing import trace_bo2023_secondary_regions
from core.methods import canonical_method
from core.network_tracing import DEFAULT_BO2023_NETWORK_MODEL, trace_network_expression
from core.region_resolution import annotate_region_candidates
from data.dao import get_atlas_metadata, get_atlas_options, get_sigset_options, table_exists


METHOD_LABELS = {
    "ensemble": ("多信号集成（推荐）", "Ensemble (recommended)"),
    "correlation": ("相关性分析", "Correlation"),
    "nnls_simplex": ("NNLS / simplex 去卷积", "NNLS / simplex deconvolution"),
    "marker": ("标记基因路径", "Marker-gene path"),
}


def _is_vsd_atlas(atlas_meta: dict) -> bool:
    text = " ".join(str(atlas_meta.get(k, "") or "").lower() for k in ["atlas_name", "build_version", "normalization", "notes"])
    return "vsd" in text or "batch_removed" in text or "batch removed" in text


def _is_bo2023_atlas(atlas_meta: dict) -> bool:
    text = " ".join(str(atlas_meta.get(k, "") or "").lower() for k in ["atlas_name", "build_version", "notes"])
    return "bo2023" in text or "wanglab" in text


def _render_vsd_mode_notice(atlas_meta: dict) -> None:
    st.info(
        tr(
            "当前选择的是 VSD + batch-corrected 参考图谱。系统已启用 VSD-compatible tracing mode：结果解释为样本与脑区表达指纹的相似性排序，不解释为 TPM 绝对丰度或真实脑区 RNA 贡献比例。",
            "The selected reference is VSD + batch-corrected. VSD-compatible tracing mode is enabled: results are interpreted as expression-fingerprint similarity rankings, not TPM abundance or biological RNA contribution fractions.",
        )
    )
    st.caption(
        tr(
            f"Atlas normalization: {atlas_meta.get('normalization', 'unknown')}；推荐优先看 correlation、rank、signature 与 bootstrap stability。",
            f"Atlas normalization: {atlas_meta.get('normalization', 'unknown')}; prioritize correlation, rank, signature evidence and bootstrap stability.",
        )
    )


def _render_network_primary(network_out: dict, show_validation_caption: bool = True) -> None:
    rows = network_out.get("results", [])
    meta = network_out.get("meta", {})
    if not rows:
        return
    top = rows[0]
    model_metadata = meta.get("model_metadata", {})
    validation = model_metadata.get("formal_route_validation") or model_metadata.get("full_loso_validation", {})
    render_section_band(
        tr("Network 主结论", "Primary Network Conclusion"),
        tr(
            "经验证的 SaleemNetworks 上层来源预测；精确 Region 排名在下方作为二级候选。",
            "Validated SaleemNetworks-level source prediction; exact Region rankings remain secondary candidates below.",
        ),
    )
    c1, c2, c3 = st.columns(3)
    c1.metric(tr("最可能来源 Network", "Top source Network"), str(top.get("network_id", "NA")))
    c2.metric(tr("Network Top1 置信度", "Network Top1 confidence"), f"{float(top.get('confidence', 0.0)):.3f}")
    c3.metric(tr("Network 模型基因数", "Network model genes"), f"{int(meta.get('n_model_genes', 0))}")
    if validation and show_validation_caption:
        st.caption(
            tr(
                f"固定算法全量 LOSO 验证：Top1 {float(validation.get('top1_accuracy', 0)):.1%}；Top3 {float(validation.get('top3_accuracy', 0)):.1%}。该性能仅代表 SaleemNetworks 上层终点。",
                f"Fixed-algorithm full LOSO validation: Top1 {float(validation.get('top1_accuracy', 0)):.1%}; Top3 {float(validation.get('top3_accuracy', 0)):.1%}. These metrics apply only to the SaleemNetworks endpoint.",
            )
        )
    network_df = pd.DataFrame(rows)
    left, right = st.columns([0.95, 1.05])
    with left:
        st.dataframe(
            network_df.rename(
                columns={
                    "network_id": tr("Network", "Network"),
                    "rank": tr("排名", "Rank"),
                    "score": tr("相关性得分", "Correlation score"),
                    "confidence": tr("置信度", "Confidence"),
                }
            ).head(5),
            use_container_width=True,
            hide_index=True,
        )
    with right:
        figure = px.bar(
            network_df.head(5),
            x="score",
            y="network_id",
            orientation="h",
            color="score",
            color_continuous_scale=["#dbeafe", "#1f7aff"],
        )
        figure.update_layout(yaxis={"categoryorder": "total ascending"}, coloraxis_showscale=False)
        st.plotly_chart(figure, use_container_width=True)


def _render_public_demo_diagnostics(network_out: dict) -> None:
    rows = network_out.get("results", [])
    meta = network_out.get("meta", {})
    scores = [float(row.get("score", 0.0)) for row in rows]
    probs = np.asarray([float(row.get("confidence", 0.0)) for row in rows], dtype=float)
    probs = probs[probs > 0]
    entropy = float(-(probs * np.log(probs)).sum()) if probs.size else 0.0
    margin = float(scores[0] - scores[1]) if len(scores) >= 2 else float("nan")
    n_overlap = int(meta.get("n_overlap_genes", 0))
    n_model = int(meta.get("n_model_genes", 0))
    coverage = float(meta.get("overlap_fraction", 0.0))

    st.info(
        tr(
            "Current manuscript validation context: projected-VSD is used only for Network Top3 beam generation; downstream resolution-group and exploratory exact-region outputs use logCPM-compatible local expression. Network projected-VSD LOSO/LOMO Top3 is 91.58%/91.33%; the three-tier route LOSO/LOMO Network Top3 is 92.19%/91.21%. Network metrics use all 819 samples; region metrics use reference-supported folds only. Exact-region output is exploratory.",
            "Current manuscript validation context: projected-VSD is used only for Network Top3 beam generation; downstream resolution-group and exploratory exact-region outputs use logCPM-compatible local expression. Network projected-VSD LOSO/LOMO Top3 is 91.58%/91.33%; the three-tier route LOSO/LOMO Network Top3 is 92.19%/91.21%. Network metrics use all 819 samples; region metrics use reference-supported folds only. Exact-region output is exploratory.",
        )
    )
    render_section_band(
        tr("Public demo diagnostics", "Public demo diagnostics"),
        tr(
            "Marker coverage, entropy, score margin and scope warnings should be read together before interpreting ranked candidates.",
            "Marker coverage, entropy, score margin and scope warnings should be read together before interpreting ranked candidates.",
        ),
    )
    render_kpi_cards(
        [
            {
                "icon": "COV",
                "label": tr("Marker coverage", "Marker coverage"),
                "value": f"{n_overlap}/{n_model}",
                "note": f"{coverage:.1%} of public Network model genes overlapped",
            },
            {
                "icon": "ENT",
                "label": tr("Entropy", "Entropy"),
                "value": f"{entropy:.3f}",
                "note": "Higher entropy means a flatter, more ambiguous Network ranking.",
            },
            {
                "icon": "MAR",
                "label": tr("Score margin", "Score margin"),
                "value": "NA" if np.isnan(margin) else f"{margin:.4f}",
                "note": "Small Top1-Top2 margin means the leading candidate is weakly separated.",
            },
            {
                "icon": "SCP",
                "label": tr("Scope warning", "Scope warning"),
                "value": tr("Coarse-only", "Coarse-only"),
                "note": "Public demo output is a Network candidate ranking, not deterministic localization.",
            },
        ]
    )
    st.warning(
        tr(
            "Biofluid outputs without independent anatomical truth are transfer stress tests, not localization-validation results.",
            "Biofluid outputs without independent anatomical truth are transfer stress tests, not localization-validation results.",
        )
    )


def _read_demo_expression(uploaded_file) -> tuple[pd.DataFrame, str]:
    name = str(uploaded_file.name).lower()
    if name.endswith((".xlsx", ".xls")):
        df = pd.read_excel(uploaded_file)
    else:
        raw = uploaded_file.getvalue()
        sep = "\t" if name.endswith((".tsv", ".txt")) or raw[:2048].count(b"\t") > raw[:2048].count(b",") else ","
        df = pd.read_csv(io.BytesIO(raw), sep=sep)

    lower_map = {str(col).strip().lower(): col for col in df.columns}
    gene_col = lower_map.get("gene_symbol") or lower_map.get("gene") or lower_map.get("symbol")
    value_col = None
    query_source = ""
    for source, candidates in [
        ("raw_counts", ["raw_counts", "raw_count", "counts", "count", "read_count", "readcount", "reads"]),
        ("logcpm", ["logcpm", "log_cpm", "log2cpm", "log2_cpm"]),
        ("logtpm_fallback", ["logtpm", "log_tpm", "log1p_tpm"]),
        ("tpm_fallback", ["tpm_value", "tpm", "expression", "value"]),
    ]:
        found = next((lower_map.get(candidate) for candidate in candidates if candidate in lower_map), None)
        if found is not None:
            value_col = found
            query_source = source
            break
    if gene_col is None or value_col is None:
        raise ValueError("Input must include gene_symbol/gene and one expression column: raw counts, logCPM, logTPM or TPM.")

    out = df[[gene_col, value_col]].copy()
    out.columns = ["gene_symbol", "query_value"]
    out["gene_symbol"] = out["gene_symbol"].astype(str).str.strip()
    out["query_value"] = pd.to_numeric(out["query_value"], errors="coerce")
    out = out.dropna(subset=["gene_symbol", "query_value"])
    out = out[out["gene_symbol"] != ""]
    if out.empty:
        raise ValueError("No valid expression rows were found.")
    out = out.groupby("gene_symbol", as_index=False)["query_value"].mean()

    if query_source == "raw_counts":
        total = float(out["query_value"].clip(lower=0).sum())
        if total <= 0:
            raise ValueError("Raw counts must sum to a positive value.")
        cpm = out["query_value"].clip(lower=0) / total * 1_000_000.0
        out["tpm_value"] = np.log2(cpm + 1.0)
    elif query_source == "tpm_fallback":
        out["tpm_value"] = np.log1p(out["query_value"].clip(lower=0))
    else:
        out["tpm_value"] = out["query_value"]

    return out[["gene_symbol", "tpm_value"]], query_source


NETWORK_DESCRIPTIONS = [
    {
        "Network": "Cingulate gyrus",
        "Location": "Medial cerebral cortex above the corpus callosum",
        "Function": "Emotion, motivation, pain, attention control and action monitoring",
    },
    {
        "Network": "Frontal (agranular frontal motor areas)",
        "Location": "Posterior frontal motor and premotor cortex",
        "Function": "Motor planning, action execution, eye movement and motor control",
    },
    {
        "Network": "Hippocampal formation",
        "Location": "Medial temporal lobe",
        "Function": "Memory, spatial navigation, episodic learning and contextual encoding",
    },
    {
        "Network": "Lateral Prefrontal Cortex",
        "Location": "Lateral prefrontal cortex",
        "Function": "Executive function, working memory, decision making and cognitive control",
    },
    {
        "Network": "Occipital/Temporal",
        "Location": "Occipital cortex and posterior temporal visual regions",
        "Function": "Visual processing, object recognition and ventral visual-stream functions",
    },
    {
        "Network": "Operculum/Insula",
        "Location": "Insula and frontal/parietal/temporal opercular cortex",
        "Function": "Interoception, taste, pain, somatosensory integration and salience processing",
    },
    {
        "Network": "Orbitomedial Prefrontal Cortex (OMPFC)",
        "Location": "Orbital and medial prefrontal cortex",
        "Function": "Reward, valuation, emotion-guided decisions and social behavior",
    },
    {
        "Network": "Parietal, and Parieto-occipital region",
        "Location": "Parietal cortex and parieto-occipital junction",
        "Function": "Spatial attention, sensory integration, visuospatial processing and action guidance",
    },
    {
        "Network": "Subcortical",
        "Location": "Subcortical structures such as thalamus, striatum and basal ganglia",
        "Function": "Sensory/motor relay, reward, movement regulation and state control",
    },
    {
        "Network": "Temporal",
        "Location": "Lateral and anterior temporal cortex",
        "Function": "Auditory, semantic, memory-related, object and social-information processing",
    },
]


def _render_network_description_table() -> None:
    render_panel_header(
        tr("Bo2023 10 个 Network 说明", "Bo2023 10-Network Guide"),
        tr(
            "公开 demo 只展示这些粗粒度解剖-功能候选来源，不展示完整 Bo2023 表达矩阵。",
            "The public demo shows only these coarse anatomical-functional candidate sources, not the full Bo2023 expression matrix.",
        ),
    )
    st.dataframe(pd.DataFrame(NETWORK_DESCRIPTIONS), use_container_width=True, hide_index=True)


def _render_public_demo_tracing() -> None:
    render_section_band(
        tr("公开 Demo 模式", "Public Demo Mode"),
        tr(
            "云端公开版推荐上传 raw counts 或 logCPM；系统不要求用户上传 VSD，也不下载或公开完整 Bo2023 数据库。",
            "The public cloud demo recommends raw counts or logCPM input; users are not asked to upload VSD, and the full Bo2023 database is not downloaded or exposed.",
        ),
    )
    st.info(
        tr(
            "最佳输入是 RNA-seq raw gene counts，系统会内部计算 logCPM；logCPM 可直接上传。TPM/logTPM 仅作为兼容旧表格的 fallback，精细脑区解释应更谨慎。",
            "The best input is RNA-seq raw gene counts, which are converted internally to logCPM; logCPM can be uploaded directly. TPM/logTPM is accepted only as a legacy fallback and should be interpreted more cautiously for fine regions.",
        )
    )
    st.info(
        tr(
            "诊断信息需联合解读：marker coverage 反映可用模型基因支持度，entropy 和 score margin 反映排序不确定性，scope warnings 提示低置信度或超出参考范围的情况。没有解剖真值的 biofluid 输出只作为 transfer stress test，不作为 localization validation。",
            "Diagnostics should be interpreted together: marker coverage reflects usable model-gene support, entropy and score margin reflect ranking ambiguity, and scope warnings flag low-confidence or out-of-reference cases. Biofluid outputs without anatomical truth are transfer stress tests, not localization-validation results.",
        )
    )
    uploaded = st.file_uploader(
        tr("上传表达矩阵 CSV/TSV/XLSX", "Upload expression matrix CSV/TSV/XLSX"),
        type=["csv", "tsv", "txt", "xlsx", "xls"],
        key="public_demo_expression_upload",
    )
    st.caption(
        tr(
            "至少包含 gene_symbol/gene 和一个表达列。推荐列名：raw_counts/count/read_count 或 logCPM；TPM/logTPM 仅为 fallback。",
            "Requires gene_symbol/gene plus one expression column. Recommended names: raw_counts/count/read_count or logCPM; TPM/logTPM is fallback only.",
        )
    )
    _render_network_description_table()
    if uploaded is None:
        return

    try:
        expr, query_source = _read_demo_expression(uploaded)
    except Exception as exc:
        st.error(f"{tr('无法读取输入文件', 'Could not read input file')}: {exc}")
        return
    if query_source in {"tpm_fallback", "logtpm_fallback"}:
        st.warning(
            tr(
                f"当前输入被识别为 {query_source}。该路径用于兼容旧表格，不等同于当前验证路线中的 Bo2023 logCPM 输入；结果应谨慎解释。",
                f"Input was detected as {query_source}. This is a legacy compatibility path and is not equivalent to the Bo2023 logCPM route used in current validation; interpret results cautiously.",
            )
        )

    render_kpi_cards(
        [
            {"icon": "GENE", "label": tr("有效基因行", "Valid gene rows"), "value": f"{len(expr):,}", "note": tr("用于 Network demo", "Used for Network demo")},
            {"icon": "SRC", "label": tr("Query 来源", "Query source"), "value": query_source, "note": tr("raw counts 会转为 logCPM", "raw counts are converted to logCPM")},
            {"icon": "MODEL", "label": tr("模型", "Model"), "value": "SaleemNetworks", "note": tr("轻量公开模型", "Lightweight public model")},
            {"icon": "DB", "label": tr("完整数据库", "Full database"), "value": tr("未使用", "Not used"), "note": tr("避免公开 Bo2023 数据", "Avoids exposing Bo2023 data")},
        ]
    )
    if st.button(tr("运行公开 Demo 溯源", "Run Public Demo Tracing"), type="primary", use_container_width=True):
        try:
            network_out = trace_network_expression(expr)
            if not network_out.get("results"):
                meta = network_out.get("meta", {})
                st.warning(
                    tr(
                        f"模型基因重叠不足：{meta.get('n_overlap_genes', 0)}/{meta.get('n_model_genes', 0)}。",
                        f"Insufficient model-gene overlap: {meta.get('n_overlap_genes', 0)}/{meta.get('n_model_genes', 0)}.",
                    )
                )
                return
            network_out.setdefault("meta", {})["query_source"] = query_source
            network_out["meta"]["input_recommendation"] = (
                "raw counts/logCPM preferred; TPM/logTPM fallback only; user-uploaded VSD is not required"
            )
            network_out["meta"].pop("model_metadata", None)
            network_out["meta"].pop("pairwise_rescue_validation", None)
            _render_network_primary(network_out, show_validation_caption=False)
            _render_public_demo_diagnostics(network_out)
            result_df = pd.DataFrame(network_out["results"])
            st.download_button(
                tr("下载 Demo JSON", "Download demo JSON"),
                json.dumps(network_out, ensure_ascii=False, indent=2),
                "public_demo_network_tracing.json",
                "application/json",
            )
            st.download_button(
                tr("下载 Demo CSV", "Download demo CSV"),
                result_df.to_csv(index=False).encode("utf-8-sig"),
                "public_demo_network_tracing.csv",
                "text/csv",
            )
        except Exception as exc:
            st.error(f"{tr('公开 Demo 分析失败', 'Public demo analysis failed')}: {exc}")


def _render_v2_results(sample_id: str, out: dict, top_regions: int, network_out: dict | None = None) -> None:
    run_id = out.get("run_id")
    results_rows = out.get("results", [])
    meta = out.get("meta", {})
    render_run_meta(meta)

    for warning in meta.get("warnings", []) or []:
        st.warning(warning)
    if meta.get("traceability") in ("low", "insufficient"):
        st.info(
            tr(
                f"当前样本可溯源性为 {meta.get('traceability')}，建议先查看 overlap 与 QC，再解释脑区排名。",
                f"The current sample traceability is {meta.get('traceability')}. Review overlap and QC before interpreting the region ranking.",
            )
        )

    st.success(tr("分析完成（v2 引擎）。", "Analysis completed (v2 engine)."))
    st.markdown(
        f'<div class="result-zone">{tr("结果区：主要发现、排名、稳定性与导出结果", "Result zone: main findings, rankings, stability and exports")}</div>',
        unsafe_allow_html=True,
    )
    render_result_hint(
        tr(
            "请将 Top 脑区解释为与 cfRNA 样本表达指纹最相近的候选脑区，而不是绝对 RNA 贡献比例。优先结合 confidence、rank/signature 证据和 bootstrap stability。",
            "Interpret Top regions as the brain areas whose reference expression fingerprints best match the cfRNA sample, not as absolute RNA contribution fractions. Combine confidence, rank/signature evidence and bootstrap stability.",
        )
        if meta.get("vsd_compatible_mode")
        else tr(
            "建议先看 Top1 脑区、置信度和稳定性，再结合排名表、fraction CI 和 signature 信号判断结论是否足够稳健。",
            "Start with the Top1 region, confidence and stability, then combine the ranking table, fraction CI and signature evidence to judge whether the conclusion is robust enough.",
        )
    )
    if meta.get("vsd_compatible_mode"):
        st.info(
            tr(
                "本次结果采用 VSD-compatible 解释口径：Top region 代表表达模式最相近的候选脑区；fraction / CI 若存在，仅表示 VSD 表达空间中的拟合权重，不代表真实组织贡献比例。",
                "VSD-compatible interpretation: Top regions are the closest expression-pattern candidates; fraction / CI, when present, are VSD-space fitting weights rather than biological contribution fractions.",
            )
        )

    if network_out:
        _render_network_primary(network_out)
        render_section_band(
            tr("Region 二级候选", "Secondary Region Candidates"),
            tr(
                "精确脑区用于在 Network 主结论下继续探索，不继承 Network 层级的验证准确率。",
                "Exact regions support exploration under the Network conclusion and do not inherit Network-level validation accuracy.",
            ),
        )
        resolution = meta.get("region_resolution_annotation", {})
        if resolution.get("enabled") and resolution.get("manual_review_recommended"):
            st.warning(
                tr(
                    "当前 Region Top1 被标记为低分辨率候选：训练数据无法稳定区分该精确脑区。"
                    f"建议报告候选组 [{resolution.get('top1_group_members', '')}] 并进行人工复核，不将精确 Top1 作为确定结论。",
                    "The current Region Top1 is flagged as low resolution: the training data do not reliably "
                    "separate this exact region. Report the candidate group "
                    f"[{resolution.get('top1_group_members', '')}] and route it for manual review rather than "
                    "treating exact Top1 as definitive.",
                )
            )
    else:
        render_section_band(
            tr("主要读数", "Primary Readout"),
            tr("先看 Top1 结果、置信度和核心支持信息。", "Review Top1, confidence and core support before deeper plots."),
        )
        render_primary_metrics(results_rows, meta)

    df_rank = pd.DataFrame(results_rows)
    if df_rank.empty:
        st.info(tr("当前运行没有返回可展示的脑区结果。", "This run did not return displayable region results."))
        return
    df_rank["region"] = df_rank["region_id"].astype(str)

    col_table, col_plot = st.columns([0.95, 1.05])
    with col_table:
        render_panel_header(
            tr("Region 二级候选表" if network_out else "脑区排名表", "Secondary Region Candidate Table" if network_out else "Region Ranking Table"),
            tr("展示候选脑区的得分、置信度、fraction 和稳定性。", "Top candidate regions with score, confidence, fraction and stability signals."),
        )
        rename = {
            "region_id": tr("脑区", "Region"),
            "rank": tr("排名", "Rank"),
            "score": tr("综合得分", "Integrated score"),
            "confidence": tr("置信度", "Confidence"),
            "fraction": tr("VSD 拟合权重" if meta.get("vsd_compatible_mode") else "贡献比例", "VSD-space fitting weight" if meta.get("vsd_compatible_mode") else "Fraction"),
            "ci_low": tr("CI 下限", "CI low"),
            "ci_high": tr("CI 上限", "CI high"),
            "stability": tr("稳定性", "Stability"),
            "reconstruction_error": tr("重建误差", "Reconstruction error"),
            "resolution_tier": tr("分辨率等级", "Resolution tier"),
            "resolution_group_members": tr("可分辨候选组", "Resolvable candidate group"),
            "manual_review_recommended": tr("人工复核", "Manual review"),
        }
        st.dataframe(df_rank.rename(columns=rename).head(top_regions), use_container_width=True, hide_index=True)
    with col_plot:
        render_panel_header(
            tr("Region 二级候选图" if network_out else "脑区排名图", "Secondary Region Candidate Plot" if network_out else "Source Ranking Plot"),
            tr("用条形图快速比较前列脑区。", "Quick visual comparison of leading source regions."),
        )
        st.plotly_chart(make_score_bar(df_rank.head(max(10, top_regions))), use_container_width=True)

    signature_cols = [c for c in ["marker_component", "support_component", "detect_component"] if c in df_rank.columns]
    if signature_cols:
        render_section_band(
            tr("Signature 证据", "Signature Evidence"),
            tr("观察哪些信号通道在推动特定脑区上升。", "Inspect which signal channels push specific regions upward."),
        )
        sig_viz = df_rank[["region_id"] + signature_cols].melt(
            id_vars="region_id",
            var_name="signature_signal",
            value_name="value",
        )
        fig_sig = px.bar(
            sig_viz,
            x="value",
            y="region_id",
            color="signature_signal",
            orientation="h",
            barmode="group",
            color_discrete_sequence=px.colors.qualitative.Pastel,
        )
        fig_sig.update_layout(height=min(720, 34 * len(df_rank) + 120), yaxis={"categoryorder": "total ascending"})
        st.plotly_chart(fig_sig, use_container_width=True)

    has_ci = {"ci_low", "ci_high", "fraction"}.issubset(df_rank.columns) and df_rank["ci_low"].notna().any()
    has_stab = "stability" in df_rank.columns and df_rank["stability"].notna().any()
    if has_ci or has_stab:
        render_section_band(
            tr("稳定性与拟合权重", "Stability and Fitting Weight") if meta.get("vsd_compatible_mode") else tr("稳定性与 Fraction", "Stability and Fraction"),
            tr(
                "Bootstrap 稳定性优先用于判断 Top 脑区排序是否稳健；VSD 拟合权重只表示标准化表达空间中的拟合支持。",
                "Bootstrap stability should be used to judge whether Top-region rankings are robust; VSD fitting weights indicate support in normalized-expression space only.",
            )
            if meta.get("vsd_compatible_mode")
            else tr("Bootstrap 置信区间与 Top1 稳定性能帮助判断结论是否稳健。", "Bootstrap confidence intervals and Top1 stability help judge robustness."),
        )
        c1, c2 = st.columns(2)
        if has_ci:
            viz_ci = df_rank.dropna(subset=["fraction"]).copy()
            viz_ci["err_plus"] = viz_ci["ci_high"] - viz_ci["fraction"]
            viz_ci["err_minus"] = viz_ci["fraction"] - viz_ci["ci_low"]
            c1.plotly_chart(make_fraction_ci_bar(viz_ci), use_container_width=True)
        else:
            c1.info(tr("当前方法或参数没有生成 fraction CI。", "This method or parameter set did not generate fraction CIs."))
        if has_stab:
            viz_st = df_rank.dropna(subset=["stability"]).copy()
            c2.plotly_chart(make_stability_bar(viz_st), use_container_width=True)
        else:
            c2.info(tr("当前方法或参数没有生成稳定性指标。", "This method or parameter set did not generate stability metrics."))

    st.markdown(f'<div class="export-zone">{tr("导出区：保存当前 Run 结果", "Export zone: save current run outputs")}</div>', unsafe_allow_html=True)
    c1, c2, c3 = st.columns(3)
    with c1:
        buf = io.BytesIO()
        df_rank.to_csv(buf, index=False, encoding="utf-8-sig")
        buf.seek(0)
        st.download_button(tr("下载 CSV", "Download CSV"), buf.getvalue(), f"v2_results_{sample_id}.csv", "text/csv")
    with c2:
        export_out = dict(out)
        if network_out:
            export_out["network_primary"] = network_out
        st.download_button(tr("下载 JSON", "Download JSON"), json.dumps(export_out, ensure_ascii=False, indent=2), f"v2_results_{sample_id}.json", "application/json")
    with c3:
        if run_id and table_exists(DB_PATH, "analysis_runs"):
            try:
                from reporting import build_run_summary, export_run_bundle

                zip_path = export_run_bundle(DB_PATH, run_id)
                with open(zip_path, "rb") as f:
                    st.download_button(tr("下载 Run 图包", "Download run bundle"), f, file_name=f"run_{run_id}.zip", mime="application/zip")
                summary_json = json.dumps(build_run_summary(DB_PATH, run_id), ensure_ascii=False, indent=2)
                st.download_button(tr("下载 Run JSON 摘要", "Download run JSON summary"), summary_json, file_name=f"run_{run_id}_summary.json", mime="application/json")
            except Exception as exc:
                st.warning(f"{tr('导出 Run 报告包失败', 'Failed to export run bundle')}: {exc}")
        else:
            st.info(tr("当前数据库未启用 analysis_runs / analysis_results，无法导出 Run 报告包。", "The current database does not expose analysis_runs / analysis_results, so a run bundle cannot be exported."))


def _render_legacy_results(sample_id: str, results: dict, top_regions: int) -> None:
    st.success(tr("分析完成（legacy 路径）。", "Analysis completed (legacy path)."))
    st.markdown(f'<div class="result-zone">{tr("结果区：legacy 分析结果", "Result zone: legacy analysis results")}</div>', unsafe_allow_html=True)
    render_result_hint(
        tr(
            "legacy 结果适合快速筛查或与 v2 结果对照。建议先看最可能来源和排名表，再决定是否切换到 v2 做更稳健解释。",
            "Legacy results are best used for quick screening or comparison with v2. Read the top source and ranking table first, then decide whether to switch to v2 for a more robust interpretation.",
        )
    )

    top_source = None
    confidence = 0.0
    if "final_ranking" in results and results["final_ranking"]:
        top_source = results["final_ranking"][0][0]
        confidence = float(results["final_ranking"][0][1])
    elif "top_regions" in results and results["top_regions"]:
        top_source = results["top_regions"][0][0]
        confidence = abs(float(results["top_regions"][0][1]["correlation"]))
    elif "components" in results and results["components"]:
        top_source = list(results["components"].items())[0][0]
        confidence = float(list(results["components"].items())[0][1])

    render_kpi_cards(
        [
            {"icon": "TOP", "label": tr("最可能来源", "Top source"), "value": top_source or tr("无法确定", "Undetermined"), "note": tr("legacy Top1 脑区", "Legacy Top1 region")},
            {"icon": "CONF", "label": tr("总体置信度", "Overall confidence"), "value": f"{confidence:.3f}", "note": tr("legacy 相对支持度", "Legacy relative support")},
            {"icon": "GENE", "label": tr("检测基因数", "Input genes"), "value": f"{results.get('n_genes', 'NA')}", "note": tr("本次分析使用的基因数", "Genes used in this analysis")},
        ]
    )

    ranking_df = None
    if "final_ranking" in results and results["final_ranking"]:
        ranking_df = pd.DataFrame(results["final_ranking"][:top_regions], columns=[tr("脑区", "Region"), tr("综合得分", "Integrated score")])
    elif "top_regions" in results and results["top_regions"]:
        ranking_df = pd.DataFrame([(r[0], r[1]["correlation"]) for r in results["top_regions"][:top_regions]], columns=[tr("脑区", "Region"), tr("相关性", "Correlation")])
    elif "components" in results and results["components"]:
        ranking_df = pd.DataFrame(list(results["components"].items())[:top_regions], columns=[tr("脑区", "Region"), tr("贡献度", "Contribution")])

    if ranking_df is not None:
        st.dataframe(ranking_df, use_container_width=True, hide_index=True)
    else:
        st.info(tr("当前 legacy 结果没有可展示的脑区排名。", "No displayable region ranking is available in the legacy result."))

    st.markdown(f'<div class="export-zone">{tr("导出区：保存 legacy 结果", "Export zone: save legacy outputs")}</div>', unsafe_allow_html=True)
    c1, c2 = st.columns(2)
    with c1:
        export_data = ranking_df if ranking_df is not None else pd.DataFrame()
        buf = io.BytesIO()
        export_data.to_csv(buf, index=False, encoding="utf-8-sig")
        buf.seek(0)
        st.download_button(tr("下载 CSV", "Download CSV"), buf.getvalue(), f"tracing_results_{sample_id}.csv", "text/csv")
    with c2:
        st.download_button(tr("下载 JSON", "Download JSON"), json.dumps(results, ensure_ascii=False, indent=2), f"tracing_results_{sample_id}.json", "application/json")


def display_source_tracing() -> None:
    db_mode = get_database_mode()
    render_page_hero(
        tr(f"{database_label(db_mode)} - 溯源分析工作台", f"{database_label(db_mode)} - Tracing Analysis Workspace"),
        tr(
            "在同一科研工作台中对已上传 cfRNA 样本执行脑区溯源、比较方法、检查稳健性并导出可汇报结果。",
            "Run brain-region source tracing on uploaded cfRNA samples, compare methods, inspect robustness and export report-ready outputs from a single scientific analysis workspace.",
        ),
        eyebrow=tr("溯源", "Tracing"),
        pills=[
            tr("样本选择", "Sample selection"),
            tr("图谱与 signature", "Atlas and signature"),
            tr("模型参数", "Model parameters"),
            tr("结果解释", "Result interpretation"),
        ],
    )
    processor = init_processor()
    tracer = init_tracer()
    samples_df = processor.get_all_samples()
    if not samples_df.empty and "species" in samples_df.columns:
        samples_df = samples_df[samples_df["species"].apply(lambda x: matches_species(x, db_mode))].copy()
    if len(samples_df) == 0:
        st.info(tr("当前云端数据库没有样本，已切换到轻量公开 Demo 模式。", "The current cloud database has no samples; switching to lightweight public demo mode."))
        _render_public_demo_tracing()
        return

    st.markdown(f'<div class="action-zone">{tr("操作区：选择待分析样本", "Action zone: choose a sample for tracing")}</div>', unsafe_allow_html=True)
    sample_id = st.selectbox(tr("选择样本", "Choose sample"), samples_df["sample_id"].astype(str).tolist(), index=0)
    cfrna_df = processor.get_sample_expression(sample_id)
    render_kpi_cards(
        [
            {"icon": "SMP", "label": tr("样本 ID", "Sample ID"), "value": sample_id, "note": tr("当前分析样本", "Current analysis sample")},
            {"icon": "GENE", "label": tr("输入基因数", "Input Genes"), "value": f"{len(cfrna_df):,}", "note": tr("样本表达矩阵中的基因数", "Genes in the sample expression matrix")},
            {"icon": "RUN", "label": tr("分析模式", "Analysis Mode"), "value": "v2 preferred", "note": tr("推荐优先使用 v2 集成路径", "v2 integrated path is recommended")},
        ]
    )

    st.markdown(f'<div class="parameter-zone">{tr("参数区：参考图谱、signature 与模型参数", "Parameter zone: atlas, signature and model settings")}</div>', unsafe_allow_html=True)
    atlas_opts = get_atlas_options(DB_PATH, species_mode=db_mode)
    if not atlas_opts:
        st.info(tr("当前数据库模式下没有可用 atlas。请先导入对应物种/数据库的参考图谱。", "No atlas is available in the current database mode. Import the matching reference atlas first."))
        return
    atlas_labels = [x[1] for x in atlas_opts]
    atlas_ids = [x[0] for x in atlas_opts]
    atlas_choice = st.selectbox(tr("Atlas 版本", "Atlas version"), atlas_labels, index=0)
    atlas_id = atlas_ids[atlas_labels.index(atlas_choice)]
    atlas_meta = get_atlas_metadata(DB_PATH, atlas_id)
    is_vsd_atlas = _is_vsd_atlas(atlas_meta)
    if is_vsd_atlas:
        _render_vsd_mode_notice(atlas_meta)
    sig_opts = get_sigset_options(DB_PATH, atlas_id)
    sig_labels = [x[1] for x in sig_opts]
    sig_ids = [x[0] for x in sig_opts]
    sig_choice = st.selectbox(tr("Signature 集合", "Signature set"), sig_labels, index=0)
    sigset_id = sig_ids[sig_labels.index(sig_choice)]
    use_v2 = st.checkbox(tr("使用发布级 v2 引擎", "Use release-grade v2 engine"), value=True)
    if is_vsd_atlas and not use_v2:
        st.warning(
            tr(
                "VSD 图谱建议使用 v2 引擎。legacy 路径不会完整应用 VSD-compatible 权重和解释元数据。",
                "VSD atlases should use the v2 engine. The legacy path does not fully apply VSD-compatible weights or interpretation metadata.",
            )
        )

    with st.expander(tr("高级模型参数", "Advanced model parameters"), expanded=False):
        p1, p2, p3, p4 = st.columns(4)
        norm_options = ["zscore", "vsd", "log1p", "tpm"] if is_vsd_atlas else ["log1p", "tpm", "zscore"]
        use_value = p1.selectbox(
            tr("输入标准化", "Input normalization"),
            norm_options,
            index=0,
            help=tr(
                "VSD 图谱默认使用 z-score，以减少上传样本 TPM/count 与 VSD reference 之间的量纲差异。",
                "VSD references default to z-score for TPM/count-derived uploads. Select vsd only when the uploaded sample is already on the same VSD-normalized scale as the atlas.",
            )
            if is_vsd_atlas
            else None,
        )
        bootstrap_n = p2.slider(tr("Bootstrap 次数", "Bootstrap iterations"), 0, 300, 100, step=25)
        bootstrap_gene_frac = p3.slider(tr("Bootstrap 基因抽样比例", "Bootstrap gene fraction"), 0.3, 1.0, 0.7, step=0.05)
        l2 = p4.number_input("L2", min_value=0.0, value=1e-4, step=1e-4, format="%.5f")
        topk = st.slider(tr("结果返回 TopK", "Return TopK"), 3, 30, 10)

    render_section_band(
        tr("方法选择", "Method Selection"),
        tr("不同路径分别强调集成、相关性、去卷积或 marker 证据。", "Different paths emphasize integration, correlation, deconvolution or marker evidence."),
    )
    method_order = ["correlation", "ensemble", "nnls_simplex", "marker"] if is_vsd_atlas else ["ensemble", "correlation", "nnls_simplex", "marker"]
    if is_vsd_atlas:
        st.info(
            tr(
                "Bo2023 VSD 的正式输出采用两级口径：同尺度 VSD 输入配合相关性分析时，SaleemNetworks 为已验证的主结论，精确 Region 为二级候选。Marker/ensemble 仍仅作探索性对照。",
                "Formal Bo2023 VSD output is two-level: for same-scale VSD input with correlation, SaleemNetworks is the validated primary conclusion and exact Regions are secondary candidates. Marker/ensemble remains exploratory.",
            )
        )
    method_display = st.radio(
        tr("选择溯源算法", "Choose tracing algorithm"),
        [tr(*METHOD_LABELS[code]) for code in method_order],
        horizontal=True,
    )
    method_map = {tr(*v): k for k, v in METHOD_LABELS.items()}
    method_key = method_map[method_display]
    if is_vsd_atlas and method_key == "nnls_simplex":
        st.warning(
            tr(
                "当前 atlas 为 VSD reference。NNLS/simplex 结果只应解释为 VSD 表达空间拟合权重，不应解释为真实脑区贡献比例。建议优先使用多信号集成或相关性分析。",
                "The current atlas is a VSD reference. NNLS/simplex results should be interpreted only as VSD-space fitting weights, not biological contribution fractions. Prefer ensemble or correlation for primary interpretation.",
            )
        )
    c4, c5 = st.columns(2)
    top_regions = c4.slider(tr("结果面板显示 Top N 脑区", "Top N regions shown in results"), 3, 20, 5)
    use_markers = c5.checkbox(
        tr("legacy 路径仅使用 marker genes", "Use marker genes only in legacy mode"),
        value=False,
        help=tr("仅对 legacy 路径生效。v2 推荐保持 signature 模式。", "Only affects the legacy path. v2 is best used with signature mode."),
    )

    st.markdown(f'<div class="action-zone">{tr("操作区：启动溯源分析", "Action zone: start tracing analysis")}</div>', unsafe_allow_html=True)
    if st.button(tr("开始分析", "Start analysis"), type="primary", use_container_width=True):
        with st.spinner(tr("正在进行脑区来源推断，请稍候...", "Running brain-region source tracing, please wait...")):
            try:
                if use_v2 and method_key != "marker":
                    out = tracer.engine_v2.trace(
                        sample_id=sample_id,
                        method=canonical_method(method_key),
                        sigset_id=sigset_id,
                        atlas_id=int(atlas_id),
                        use_value=use_value,
                        bootstrap_n=int(bootstrap_n),
                        bootstrap_gene_frac=float(bootstrap_gene_frac),
                        l2=float(l2),
                        topk=int(max(topk, top_regions)),
                        vsd_compatible=bool(is_vsd_atlas),
                    )
                    network_out = None
                    if is_vsd_atlas and _is_bo2023_atlas(atlas_meta) and method_key == "correlation":
                        if use_value == "vsd" and DEFAULT_BO2023_NETWORK_MODEL.exists():
                            network_out = trace_network_expression(cfrna_df)
                            if not network_out.get("results"):
                                st.warning(
                                    tr(
                                        "当前样本覆盖的 Network 模型基因不足，不能输出经验证的 Network 主结论；下方仅展示 Region 二级候选。",
                                        "The sample does not cover enough Network-model genes for a validated primary Network conclusion; only secondary Region candidates are shown below.",
                                    )
                                )
                                network_out = None
                        elif use_value != "vsd":
                            st.info(
                                tr(
                                    "Network 主结论仅对与 Bo2023 reference 同尺度的 VSD 输入启用；当前输入尺度未经过该路径验证，因此仅展示 Region 候选结果。",
                                    "The primary Network conclusion is enabled only for input on the same VSD scale as the Bo2023 reference. This input scale is not validated for that path, so only Region candidates are shown.",
                                )
                            )
                    if network_out:
                        secondary_out = trace_bo2023_secondary_regions(
                            cfrna_df,
                            network_out,
                            DB_PATH,
                            int(atlas_id),
                            topk=max(int(top_regions), 30),
                        )
                        if secondary_out.get("results"):
                            out = secondary_out
                        else:
                            st.warning(
                                tr(
                                    "Bo2023 专用 Region 二级候选评分未能生成结果，已回退到通用 correlation Region 排名。",
                                    "The Bo2023-specific secondary Region scorer did not return results; falling back to the generic correlation Region ranking.",
                                )
                            )
                        out = annotate_region_candidates(out, network_out)
                    _render_v2_results(sample_id, out, top_regions, network_out=network_out)
                else:
                    if use_markers:
                        markers = tracer.region_signatures[tracer.region_signatures["is_marker"] == 1]["gene_symbol"].unique()
                        cfrna_df_use = cfrna_df[cfrna_df["gene_symbol"].isin(markers)].copy()
                    else:
                        cfrna_df_use = cfrna_df

                    if method_key == "ensemble":
                        results = tracer.integrated_tracing(cfrna_df_use)
                    elif method_key == "correlation":
                        results = tracer.correlation_based_tracing(cfrna_df_use, top_regions)
                    elif method_key == "nnls_simplex":
                        results = tracer.deconvolution_based_tracing(cfrna_df_use, method="NMF")
                    else:
                        results = tracer.marker_gene_based_tracing(cfrna_df_use)

                    tracer.save_tracing_results(sample_id, results)
                    _render_legacy_results(sample_id, results, top_regions)
            except Exception as exc:
                st.error(tr("分析失败：当前样本与参数组合未能完成溯源计算。", "Analysis failed: the current sample and parameter combination could not complete tracing."))
                st.info(f"{tr('原始错误', 'Original error')}: {exc}")
                with st.expander(tr("开发者调试信息", "Developer debug details"), expanded=False):
                    st.code(traceback.format_exc(), language="python")
