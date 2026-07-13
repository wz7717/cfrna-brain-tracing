from __future__ import annotations

import io
import json
import sqlite3
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
from app.shared import DB_PATH, init_processor, is_public_demo_mode, render_page_hero
from core.bo2023_region_tracing import trace_bo2023_secondary_regions
from core.network_tracing import DEFAULT_BO2023_NETWORK_MODEL, trace_network_expression
from data.dao import get_atlas_options, table_exists


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
            width="stretch",
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
        st.plotly_chart(figure, width="stretch")


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

    render_section_band(
        tr("Result confidence checks", "Result confidence checks"),
        tr(
            "Read marker coverage, entropy and score margin together before interpreting the ranked Network candidates.",
            "Read marker coverage, entropy and score margin together before interpreting the ranked Network candidates.",
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
                "label": tr("Interpretation scope", "Interpretation scope"),
                "value": tr("Network-level", "Network-level"),
                "note": "This output is a ranked Network-level candidate list, not deterministic anatomical localization.",
            },
        ]
    )
    st.warning(
        tr(
            "Biofluid outputs without patient-level anatomical truth are projection-feasibility or transfer-stress analyses, not localization-accuracy results.",
            "Biofluid outputs without patient-level anatomical truth are projection-feasibility or transfer-stress analyses, not localization-accuracy results.",
        )
    )


def _read_demo_expression(uploaded_file) -> tuple[pd.DataFrame, str]:
    name = str(uploaded_file.name).lower()
    if name.endswith(".xlsx"):
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
        out["read_count"] = out["query_value"].clip(lower=0)
        if float(out["read_count"].sum()) <= 0:
            raise ValueError("Raw counts must sum to a positive value.")
        return out[["gene_symbol", "read_count"]], query_source
    elif query_source in {"logcpm", "logtpm_fallback"}:
        out["log_tpm"] = out["query_value"]
        return out[["gene_symbol", "log_tpm"]], query_source
    elif query_source == "tpm_fallback":
        out["tpm_value"] = out["query_value"].clip(lower=0)
    else:
        out["log_tpm"] = out["query_value"]
        return out[["gene_symbol", "log_tpm"]], query_source

    return out[["gene_symbol", "tpm_value"]], query_source


def _select_locked_bo2023_atlas(db_mode: str) -> tuple[int, str]:
    atlas_opts = get_atlas_options(DB_PATH, species_mode=db_mode)
    if not atlas_opts:
        raise RuntimeError(
            tr(
                "当前数据库没有可用的 Bo2023 参考图谱，无法运行三层溯源。",
                "No Bo2023 reference atlas is available in the current database; three-tier tracing cannot run.",
            )
        )
    for atlas_id, label in atlas_opts:
        text = str(label).lower()
        if "bo2023" in text or "wanglab" in text or "vsd" in text:
            return int(atlas_id), str(label)
    raise RuntimeError(
        tr(
            "当前数据库的参考图谱列表中未找到 Bo2023 图谱，无法运行三层溯源。",
            "No Bo2023 atlas was found among the current database references; three-tier tracing cannot run.",
        )
    )


def _locked_route_expression(cfrna_df: pd.DataFrame) -> tuple[pd.DataFrame, str]:
    expr = cfrna_df.copy()
    if "gene_symbol" not in expr.columns:
        raise ValueError("Expression data must include gene_symbol.")
    expr["gene_symbol"] = expr["gene_symbol"].astype(str).str.strip()
    expr = expr[expr["gene_symbol"] != ""].copy()
    if expr.empty:
        raise ValueError("Expression data has no valid gene_symbol rows.")

    if "read_count" in expr.columns:
        read_count = pd.to_numeric(expr["read_count"], errors="coerce").fillna(0.0).clip(lower=0.0)
        if float(read_count.sum()) > 0:
            out = pd.DataFrame({"gene_symbol": expr["gene_symbol"], "read_count": read_count})
            return out.groupby("gene_symbol", as_index=False)["read_count"].mean(), "raw_counts"

    if "expression_unit" in expr.columns and "tpm_value" in expr.columns:
        units = expr["expression_unit"].astype(str).str.lower()
        if units.str.contains("logcpm|log_cpm|log2cpm|log2_cpm", regex=True, na=False).any():
            values = pd.to_numeric(expr["tpm_value"], errors="coerce")
            out = pd.DataFrame({"gene_symbol": expr["gene_symbol"], "log_tpm": values})
            return out.dropna(subset=["log_tpm"]).groupby("gene_symbol", as_index=False)["log_tpm"].mean(), "stored_logCPM"

    if "log_tpm" in expr.columns:
        values = pd.to_numeric(expr["log_tpm"], errors="coerce")
        if values.notna().any():
            out = pd.DataFrame({"gene_symbol": expr["gene_symbol"], "log_tpm": values})
            return out.dropna(subset=["log_tpm"]).groupby("gene_symbol", as_index=False)["log_tpm"].mean(), "stored_log"

    if "tpm_value" in expr.columns:
        values = pd.to_numeric(expr["tpm_value"], errors="coerce").fillna(0.0).clip(lower=0.0)
        out = pd.DataFrame({"gene_symbol": expr["gene_symbol"], "tpm_value": values})
        return out.groupby("gene_symbol", as_index=False)["tpm_value"].mean(), "TPM_fallback"

    raise ValueError("Expression data must include read_count, logCPM/log_tpm, or tpm_value.")


def _render_resolution_group_top3(out: dict) -> None:
    group_rows = out.get("meta", {}).get("region_resolution_annotation", {}).get("group_ranking", [])
    if not group_rows:
        return
    render_panel_header(
        tr("Resolution group Top3", "Resolution Group Top3"),
        tr(
            "分辨率组用于报告 Bo2023 训练数据可稳定区分的候选范围。",
            "Resolution groups report the candidate scope that the Bo2023 training data can separate more reliably.",
        ),
    )
    st.dataframe(pd.DataFrame(group_rows).head(3), width="stretch", hide_index=True)


def _run_locked_bo2023_route(expr: pd.DataFrame, atlas_id: int, topk: int = 30) -> tuple[dict, dict]:
    network_out = trace_network_expression(expr, enable_pairwise_rescue=False)
    if not network_out.get("results"):
        meta = network_out.get("meta", {})
        raise ValueError(
            f"Insufficient Network model-gene overlap: {meta.get('n_overlap_genes', 0)}/"
            f"{meta.get('n_model_genes', 0)}."
        )
    out = trace_bo2023_secondary_regions(expr, network_out, DB_PATH, int(atlas_id), topk=max(int(topk), 3))
    if not out.get("results"):
        meta = out.get("meta", {})
        raise ValueError(str(meta.get("error") or "Bo2023 three-tier route returned no region candidates."))
    return network_out, out


def _render_locked_three_tier_results(sample_id: str, out: dict, network_out: dict) -> None:
    network_df = pd.DataFrame(network_out.get("results", [])).head(3)
    exact_df = pd.DataFrame(out.get("results", [])).head(3)
    group_rows = out.get("meta", {}).get("region_resolution_annotation", {}).get("group_ranking", [])
    group_df = pd.DataFrame(group_rows).head(3)

    st.success(tr("三层溯源已完成。", "Three-tier tracing completed."))
    render_section_band(
        tr("Network Top3", "Network Top3"),
        tr("这是论文主路线的上层主结论。", "This is the primary upper-level conclusion in the manuscript route."),
    )
    if network_df.empty:
        st.warning(tr("没有可展示的 Network 候选。", "No Network candidates are available."))
    else:
        st.dataframe(network_df, width="stretch", hide_index=True)

    _render_public_demo_diagnostics(network_out)

    render_section_band(
        tr("Resolution group Top3", "Resolution Group Top3"),
        tr(
            "这是更稳健的可分辨候选组层级，优先用于报告 exact-region 低分辨率时的不确定范围。",
            "This is the more robust resolvable candidate-group tier, useful when exact-region calls are low-resolution.",
        ),
    )
    if group_df.empty:
        st.info(tr("当前结果没有 resolution group ranking。", "No resolution group ranking is available for this result."))
    else:
        st.dataframe(group_df, width="stretch", hide_index=True)

    render_section_band(
        tr("Exact-region exploratory Top3", "Exact-Region Exploratory Top3"),
        tr(
            "精确脑区只作为探索性候选，不继承 Network 层级的验证准确率。",
            "Exact regions are exploratory candidates and do not inherit Network-level validation accuracy.",
        ),
    )
    for warning in out.get("meta", {}).get("warnings", []) or []:
        st.warning(str(warning))
    resolution = out.get("meta", {}).get("region_resolution_annotation", {})
    if resolution.get("manual_review_recommended"):
        st.warning(
            tr(
                f"Exact-region Top1 需要人工复核；建议同时报告候选组 [{resolution.get('top1_group_members', '')}]。",
                f"Exact-region Top1 needs manual review; also report candidate group [{resolution.get('top1_group_members', '')}].",
            )
        )
    if exact_df.empty:
        st.warning(tr("没有可展示的 exact-region 候选。", "No exact-region candidates are available."))
    else:
        st.dataframe(exact_df, width="stretch", hide_index=True)

    export_df = pd.concat(
        [
            network_df.assign(tier="network_top3"),
            group_df.assign(tier="resolution_group_top3"),
            exact_df.assign(tier="exact_region_exploratory_top3"),
        ],
        ignore_index=True,
        sort=False,
    )
    export_out = dict(out)
    export_out["network_primary"] = network_out
    c1, c2 = st.columns(2)
    with c1:
        st.download_button(
            tr("下载 CSV", "Download CSV"),
            export_df.to_csv(index=False).encode("utf-8-sig"),
            f"bo2023_three_tier_top3_{sample_id}.csv",
            "text/csv",
        )
    with c2:
        st.download_button(
            tr("下载 JSON", "Download JSON"),
            json.dumps(export_out, ensure_ascii=False, indent=2),
            f"bo2023_three_tier_{sample_id}.json",
            "application/json",
        )


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
    st.dataframe(pd.DataFrame(NETWORK_DESCRIPTIONS), width="stretch", hide_index=True)


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
        st.dataframe(df_rank.rename(columns=rename).head(top_regions), width="stretch", hide_index=True)
    with col_plot:
        render_panel_header(
            tr("Region 二级候选图" if network_out else "脑区排名图", "Secondary Region Candidate Plot" if network_out else "Source Ranking Plot"),
            tr("用条形图快速比较前列脑区。", "Quick visual comparison of leading source regions."),
        )
        st.plotly_chart(make_score_bar(df_rank.head(max(10, top_regions))), width="stretch")

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
        st.plotly_chart(fig_sig, width="stretch")

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
            c1.plotly_chart(make_fraction_ci_bar(viz_ci), width="stretch")
        else:
            c1.info(tr("当前方法或参数没有生成 fraction CI。", "This method or parameter set did not generate fraction CIs."))
        if has_stab:
            viz_st = df_rank.dropna(subset=["stability"]).copy()
            c2.plotly_chart(make_stability_bar(viz_st), width="stretch")
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
        st.dataframe(ranking_df, width="stretch", hide_index=True)
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


def _get_all_samples_or_empty(processor) -> pd.DataFrame:
    """Treat an unpopulated local SQLite file as an upload-only demo workspace."""
    try:
        return processor.get_all_samples()
    except (sqlite3.Error, pd.errors.DatabaseError):
        return pd.DataFrame()


def _render_public_demo_tracing() -> None:
    render_section_band(
        tr("上传表达矩阵", "Upload Expression Matrix"),
        tr(
            "上传 gene_symbol 加 raw counts 或 logCPM；系统运行三层溯源路线。",
            "Upload gene_symbol plus raw counts or logCPM; the app runs the three-tier tracing route.",
        ),
    )
    st.info(
        tr(
            "主结论看 Network Top3；resolution group Top3 是更稳健的可分辨候选范围；exact-region Top3 只作探索性定位候选。",
            "Read Network Top3 as the primary conclusion; resolution group Top3 is the more robust candidate scope; exact-region Top3 remains exploratory.",
        )
    )
    uploaded = st.file_uploader(
        tr("上传表达矩阵 CSV/TSV/XLSX", "Upload expression matrix CSV/TSV/XLSX"),
        type=["csv", "tsv", "txt", "xlsx"],
        key="public_demo_expression_upload_locked",
    )
    st.caption(
        tr(
            "至少包含 gene_symbol/gene 和一个表达列。推荐列名：raw_counts/count/read_count 或 logCPM。",
            "Requires gene_symbol/gene plus one expression column. Recommended names: raw_counts/count/read_count or logCPM.",
        )
    )
    with st.expander(tr("查看 10 个 Network 候选范围", "View the 10 Network candidate scopes"), expanded=False):
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
                f"当前输入被识别为 {query_source}，仅作为旧表格兼容入口；论文主路线推荐 raw counts 或 logCPM。",
                f"Input was detected as {query_source}; this is only a legacy table compatibility path. The manuscript route prefers raw counts or logCPM.",
            )
        )

    try:
        atlas_id, atlas_label = _select_locked_bo2023_atlas(get_database_mode())
    except RuntimeError as exc:
        st.error(str(exc))
        return
    render_kpi_cards(
        [
            {"icon": "GENE", "label": tr("有效基因行", "Valid gene rows"), "value": f"{len(expr):,}", "note": tr("用于 Bo2023 三层路线", "Used for the Bo2023 three-tier route")},
            {"icon": "SRC", "label": tr("输入口径", "Query source"), "value": query_source, "note": tr("raw counts 会在模型内转 logCPM", "raw counts are converted to logCPM inside the model")},
            {"icon": "REF", "label": tr("参考图谱", "Reference"), "value": "Bo2023", "note": atlas_label},
            {"icon": "SCOPE", "label": tr("输出范围", "Output scope"), "value": tr("三层 Top3", "Three-tier Top3"), "note": "Network / resolution group / exact-region"},
        ]
    )

    if st.button(tr("运行三层溯源", "Run three-tier tracing"), type="primary", width="stretch"):
        try:
            network_out, out = _run_locked_bo2023_route(expr, atlas_id, topk=30)
            network_out.setdefault("meta", {})["query_source"] = query_source
            network_out["meta"]["input_recommendation"] = "raw counts/logCPM preferred; TPM/logTPM fallback only"
            network_out["meta"].pop("model_metadata", None)
            network_out["meta"].pop("pairwise_rescue_validation", None)
            _render_locked_three_tier_results("uploaded_expression", out, network_out)
        except Exception as exc:
            st.error(f"{tr('三层溯源运行失败', 'Three-tier tracing failed')}: {exc}")


def display_source_tracing() -> None:
    db_mode = get_database_mode()
    render_page_hero(
        tr(f"{database_label(db_mode)} - 三层溯源", f"{database_label(db_mode)} - Three-Tier Tracing"),
        tr(
            "上传或选择 cfRNA 表达矩阵后，只运行论文主路线：Network Top3 -> resolution group Top3 -> exact-region exploratory Top3。",
            "Upload or select a cfRNA expression matrix and run only the manuscript route: Network Top3 -> resolution group Top3 -> exact-region exploratory Top3.",
        ),
        eyebrow=tr("溯源分析", "Tracing Analysis"),
        pills=[
            "gene_symbol + raw_counts/logCPM",
            "Network Top3",
            "resolution group Top3",
            "exact-region exploratory Top3",
        ],
    )
    if is_public_demo_mode():
        _render_public_demo_tracing()
        return

    processor = init_processor()
    samples_df = _get_all_samples_or_empty(processor)
    if not samples_df.empty and "species" in samples_df.columns:
        samples_df = samples_df[samples_df["species"].apply(lambda x: matches_species(x, db_mode))].copy()
    if len(samples_df) == 0:
        st.info(tr("当前数据库没有可分析样本，可直接上传表达矩阵运行 Demo。", "No analyzable samples are available; upload an expression matrix to run the demo."))
        _render_public_demo_tracing()
        return

    st.markdown(f'<div class="action-zone">{tr("操作区：选择样本并运行三层溯源", "Action zone: choose a sample and run three-tier tracing")}</div>', unsafe_allow_html=True)
    sample_id = st.selectbox(tr("选择样本", "Choose sample"), samples_df["sample_id"].astype(str).tolist(), index=0)
    cfrna_df = processor.get_sample_expression(sample_id)
    expr, query_source = _locked_route_expression(cfrna_df)
    try:
        atlas_id, atlas_label = _select_locked_bo2023_atlas(db_mode)
    except RuntimeError as exc:
        st.error(str(exc))
        return
    render_kpi_cards(
        [
            {"icon": "SMP", "label": tr("样本 ID", "Sample ID"), "value": sample_id, "note": tr("当前分析样本", "Current analysis sample")},
            {"icon": "GENE", "label": tr("有效基因行", "Valid gene rows"), "value": f"{len(expr):,}", "note": tr("进入锁定路线的基因数", "Genes entering the locked route")},
            {"icon": "SRC", "label": tr("输入口径", "Input scale"), "value": query_source, "note": tr("优先 raw counts/logCPM", "raw counts/logCPM preferred")},
            {"icon": "REF", "label": tr("参考图谱", "Reference"), "value": "Bo2023", "note": atlas_label},
        ]
    )
    if query_source == "TPM_fallback":
        st.warning(tr("该样本未检出 read_count 或 logCPM/log_tpm，当前仅使用 TPM fallback；结果应降级解释。", "This sample has no read_count or logCPM/log_tpm, so TPM fallback is used and interpretation should be downgraded."))

    if st.button(tr("运行三层溯源", "Run three-tier tracing"), type="primary", width="stretch"):
        with st.spinner(tr("正在运行三层候选排名...", "Running three-tier candidate ranking...")):
            try:
                network_out, out = _run_locked_bo2023_route(expr, atlas_id, topk=30)
                network_out.setdefault("meta", {})["query_source"] = query_source
                _render_locked_three_tier_results(sample_id, out, network_out)
            except Exception as exc:
                st.error(tr("分析失败：三层溯源未能完成。", "Analysis failed: three-tier tracing could not complete."))
                st.info(f"{tr('原始错误', 'Original error')}: {exc}")
                with st.expander(tr("开发者调试信息", "Developer debug details"), expanded=False):
                    st.code(traceback.format_exc(), language="python")
