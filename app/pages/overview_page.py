from __future__ import annotations

from html import escape
import sqlite3
from typing import Iterable

import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import streamlit as st

from app.components.layout import (
    render_kpi_cards,
    render_mini_cards,
    render_panel_header,
    render_section_band,
    render_update_list,
)
from app.database_mode import database_label, get_database_mode, matches_species
from app.i18n import tr
from app.shared import DB_PATH, render_page_hero, render_result_hint


def _safe_query(query: str, params: Iterable | None = None) -> pd.DataFrame:
    try:
        with sqlite3.connect(DB_PATH) as conn:
            return pd.read_sql_query(query, conn, params=list(params or []))
    except Exception:
        return pd.DataFrame()


def _atlas_filter_clause(db_mode: str, alias: str = "r") -> tuple[str, list[int]]:
    try:
        with sqlite3.connect(DB_PATH) as conn:
            atlas_ids = _atlas_ids_for_mode(conn, db_mode)
    except Exception:
        atlas_ids = []
    if atlas_ids:
        placeholders = ",".join(["?"] * len(atlas_ids))
        return f" AND {alias}.atlas_id IN ({placeholders})", atlas_ids
    if db_mode == "human":
        return " AND 1=0", []
    return "", []


def _scalar_count(query: str, params: Iterable | None = None) -> int:
    df = _safe_query(query, params)
    if df.empty:
        return 0
    try:
        return int(df.iloc[0, 0] or 0)
    except Exception:
        return 0


def _candidate_values(query: str, params: Iterable | None = None, column: str = "value") -> list[str]:
    df = _safe_query(query, params)
    if df.empty or column not in df.columns:
        return []
    return [str(v) for v in df[column].dropna().astype(str).tolist() if str(v).strip()]


def _build_workbench_search(query: str, db_mode: str) -> dict[str, pd.DataFrame]:
    term = str(query or "").strip()
    if not term:
        return {}
    like = f"%{term}%"
    prefix = f"{term}%"
    atlas_clause, atlas_params = _atlas_filter_clause(db_mode, "r")

    samples = _safe_query(
        """
        SELECT sample_id, subject_id, species, diagnosis, sample_type, qc_status, collection_date
        FROM cfrna_samples
        WHERE sample_id LIKE ?
           OR COALESCE(subject_id, '') LIKE ?
           OR COALESCE(diagnosis, '') LIKE ?
           OR COALESCE(surgery_region, '') LIKE ?
        ORDER BY COALESCE(collection_date, '') DESC, sample_id
        LIMIT 8
        """,
        [like, like, like, like],
    )
    if not samples.empty:
        samples = samples[samples["species"].apply(lambda x: matches_species(x, db_mode))].copy()

    gene_candidates = _candidate_values(
        f"""
        SELECT DISTINCT r.gene_symbol AS value
        FROM reference_expression r
        WHERE (r.gene_symbol = ? COLLATE NOCASE OR r.gene_symbol LIKE ? COLLATE NOCASE)
        {atlas_clause}
        ORDER BY r.gene_symbol
        LIMIT 8
        """,
        [term, prefix, *atlas_params],
    )
    genes = pd.DataFrame()
    if gene_candidates:
        placeholders = ",".join(["?"] * len(gene_candidates))
        genes = _safe_query(
            f"""
            SELECT
                r.gene_symbol,
                MAX(COALESCE(r.gene_name, '')) AS gene_name,
                COUNT(DISTINCT r.region_id) AS region_count,
                COUNT(DISTINCT r.atlas_id) AS atlas_count,
                AVG(r.avg_tpm) AS mean_tpm,
                MAX(r.sample_count) AS sample_count,
                GROUP_CONCAT(DISTINCT COALESCE(a.atlas_name, 'Reference atlas')) AS source_database,
                MAX(COALESCE(a.created_at, '')) AS updated_at
            FROM reference_expression r
            LEFT JOIN atlas_versions a ON a.atlas_id = r.atlas_id
            WHERE r.gene_symbol IN ({placeholders})
            GROUP BY r.gene_symbol
            ORDER BY region_count DESC, mean_tpm DESC
            LIMIT 8
            """,
            gene_candidates,
        )

    region_candidates = _candidate_values(
        f"""
        SELECT DISTINCT r.region_id AS value
        FROM reference_expression r
        WHERE r.region_id LIKE ? COLLATE NOCASE
        {atlas_clause}
        ORDER BY r.region_id
        LIMIT 8
        """,
        [prefix, *atlas_params],
    )
    regions = pd.DataFrame()
    if region_candidates:
        placeholders = ",".join(["?"] * len(region_candidates))
        regions = _safe_query(
            f"""
            SELECT
                r.region_id,
                MAX(COALESCE(r.region_name, r.region_id)) AS region_name,
                COUNT(DISTINCT r.gene_symbol) AS gene_count,
                COUNT(DISTINCT r.atlas_id) AS atlas_count,
                MAX(r.sample_count) AS sample_count,
                GROUP_CONCAT(DISTINCT COALESCE(a.atlas_name, 'Reference atlas')) AS source_database,
                MAX(COALESCE(a.created_at, '')) AS updated_at
            FROM reference_expression r
            LEFT JOIN atlas_versions a ON a.atlas_id = r.atlas_id
            WHERE r.region_id IN ({placeholders})
            GROUP BY r.region_id
            ORDER BY gene_count DESC
            LIMIT 8
            """,
            region_candidates,
        )

    celltype_candidates = _candidate_values(
        f"""
        SELECT DISTINCT r.cell_type_marker AS value
        FROM reference_expression r
        WHERE r.cell_type_marker IS NOT NULL
          AND r.cell_type_marker != ''
          AND r.cell_type_marker LIKE ? COLLATE NOCASE
        {atlas_clause}
        ORDER BY r.cell_type_marker
        LIMIT 8
        """,
        [prefix, *atlas_params],
    )
    celltypes = pd.DataFrame()
    if celltype_candidates:
        placeholders = ",".join(["?"] * len(celltype_candidates))
        celltypes = _safe_query(
            f"""
            SELECT
                r.cell_type_marker AS cell_type,
                COUNT(DISTINCT r.gene_symbol) AS marker_count,
                COUNT(DISTINCT r.region_id) AS region_count,
                COUNT(DISTINCT r.atlas_id) AS atlas_count,
                GROUP_CONCAT(DISTINCT COALESCE(a.atlas_name, 'Reference atlas')) AS source_database,
                MAX(COALESCE(a.created_at, '')) AS updated_at
            FROM reference_expression r
            LEFT JOIN atlas_versions a ON a.atlas_id = r.atlas_id
            WHERE r.cell_type_marker IN ({placeholders})
            GROUP BY r.cell_type_marker
            ORDER BY marker_count DESC
            LIMIT 8
            """,
            celltype_candidates,
        )

    runs = _safe_query(
        """
        SELECT
            ar.run_id,
            ar.sample_id,
            ar.atlas_id,
            ar.method,
            ar.created_at,
            av.atlas_name,
            av.build_version,
            best.region_id AS top_region,
            best.score AS confidence
        FROM analysis_runs ar
        LEFT JOIN atlas_versions av ON av.atlas_id = ar.atlas_id
        LEFT JOIN analysis_results best ON best.run_id = ar.run_id AND best.rank = 1
        WHERE ar.run_id LIKE ?
           OR ar.sample_id LIKE ?
           OR COALESCE(ar.method, '') LIKE ?
           OR COALESCE(best.region_id, '') LIKE ?
        ORDER BY COALESCE(ar.created_at, '') DESC
        LIMIT 8
        """,
        [like, like, like, like],
    )

    return {
        "samples": samples,
        "genes": genes,
        "regions": regions,
        "celltypes": celltypes,
        "runs": runs,
    }


def _result_card(title: str, kind: str, meta: list[str], note: str = "", confidence: str = "") -> str:
    meta_html = "".join(f'<span class="workbench-meta-chip">{escape(str(m))}</span>' for m in meta if str(m).strip())
    confidence_html = f'<div class="workbench-confidence">{escape(confidence)}</div>' if confidence else ""
    note_html = f'<div class="workbench-card-note">{escape(note)}</div>' if note else ""
    return (
        '<div class="workbench-result-card">'
        f'<div class="workbench-card-top"><span class="workbench-kind">{escape(kind)}</span>{confidence_html}</div>'
        f'<div class="workbench-card-title">{escape(title)}</div>'
        f'{note_html}'
        f'<div class="workbench-meta-row">{meta_html}</div>'
        '</div>'
    )


def _render_workbench_results(results: dict[str, pd.DataFrame], query: str) -> None:
    total = sum(len(df) for df in results.values() if isinstance(df, pd.DataFrame) and not df.empty)
    if total == 0:
        st.info(tr("没有找到匹配的样本、基因、脑区、细胞类型或 run。", "No matching samples, genes, regions, cell types or runs were found."))
        return

    st.caption(tr(f"检索到 {total} 条相关结果，已按对象类型分组。", f"{total} related result(s), grouped by entity type."))
    tabs = st.tabs(
        [
            tr(f"样本 {len(results.get('samples', []))}", f"Samples {len(results.get('samples', []))}"),
            tr(f"基因 {len(results.get('genes', []))}", f"Genes {len(results.get('genes', []))}"),
            tr(f"脑区 {len(results.get('regions', []))}", f"Regions {len(results.get('regions', []))}"),
            tr(f"细胞类型 {len(results.get('celltypes', []))}", f"Cell types {len(results.get('celltypes', []))}"),
            tr(f"Runs {len(results.get('runs', []))}", f"Runs {len(results.get('runs', []))}"),
        ]
    )

    with tabs[0]:
        blocks = []
        for _, row in results.get("samples", pd.DataFrame()).iterrows():
            blocks.append(
                _result_card(
                    str(row.get("sample_id", "NA")),
                    tr("样本", "Sample"),
                    [
                        f"{tr('来源', 'Source')}: cfrna_samples",
                        f"{tr('物种', 'Species')}: {row.get('species', 'NA')}",
                        f"QC: {row.get('qc_status', 'Unknown') or 'Unknown'}",
                        f"{tr('更新时间', 'Updated')}: {str(row.get('collection_date', ''))[:10] or 'NA'}",
                    ],
                    note=f"{row.get('subject_id', '') or ''} {row.get('diagnosis', '') or ''} {row.get('sample_type', '') or ''}".strip(),
                )
            )
        st.markdown(f'<div class="workbench-results-grid">{"".join(blocks)}</div>', unsafe_allow_html=True) if blocks else st.info(tr("没有样本匹配。", "No sample matches."))

    with tabs[1]:
        blocks = []
        for _, row in results.get("genes", pd.DataFrame()).iterrows():
            blocks.append(
                _result_card(
                    str(row.get("gene_symbol", "NA")),
                    tr("基因", "Gene"),
                    [
                        f"{tr('来源数据库', 'Source database')}: {row.get('source_database', 'Reference atlas')}",
                        f"{tr('样本数', 'Samples')}: {int(row.get('sample_count') or 0)}",
                        f"{tr('参考版本', 'Reference versions')}: {int(row.get('atlas_count') or 0)}",
                        f"{tr('更新时间', 'Updated')}: {str(row.get('updated_at', ''))[:10] or 'NA'}",
                    ],
                    note=f"{row.get('gene_name', '') or ''} | {tr('覆盖脑区', 'Regions')}: {int(row.get('region_count') or 0)}",
                )
            )
        st.markdown(f'<div class="workbench-results-grid">{"".join(blocks)}</div>', unsafe_allow_html=True) if blocks else st.info(tr("没有基因匹配。", "No gene matches."))

    with tabs[2]:
        blocks = []
        for _, row in results.get("regions", pd.DataFrame()).iterrows():
            blocks.append(
                _result_card(
                    str(row.get("region_id", "NA")),
                    tr("脑区", "Brain region"),
                    [
                        f"{tr('来源数据库', 'Source database')}: {row.get('source_database', 'Reference atlas')}",
                        f"{tr('样本数', 'Samples')}: {int(row.get('sample_count') or 0)}",
                        f"{tr('参考版本', 'Reference versions')}: {int(row.get('atlas_count') or 0)}",
                        f"{tr('更新时间', 'Updated')}: {str(row.get('updated_at', ''))[:10] or 'NA'}",
                    ],
                    note=f"{row.get('region_name', '') or ''} | {tr('基因记录', 'Gene records')}: {int(row.get('gene_count') or 0)}",
                )
            )
        st.markdown(f'<div class="workbench-results-grid">{"".join(blocks)}</div>', unsafe_allow_html=True) if blocks else st.info(tr("没有脑区匹配。", "No region matches."))

    with tabs[3]:
        blocks = []
        for _, row in results.get("celltypes", pd.DataFrame()).iterrows():
            blocks.append(
                _result_card(
                    str(row.get("cell_type", "NA")),
                    tr("细胞类型", "Cell type"),
                    [
                        f"{tr('来源数据库', 'Source database')}: {row.get('source_database', 'Reference atlas')}",
                        f"{tr('marker 数', 'Markers')}: {int(row.get('marker_count') or 0)}",
                        f"{tr('覆盖脑区', 'Regions')}: {int(row.get('region_count') or 0)}",
                        f"{tr('更新时间', 'Updated')}: {str(row.get('updated_at', ''))[:10] or 'NA'}",
                    ],
                )
            )
        st.markdown(f'<div class="workbench-results-grid">{"".join(blocks)}</div>', unsafe_allow_html=True) if blocks else st.info(tr("当前 reference_expression 中没有匹配的细胞类型 marker。", "No matching cell type marker is available in reference_expression."))

    with tabs[4]:
        blocks = []
        for _, row in results.get("runs", pd.DataFrame()).iterrows():
            confidence = ""
            if pd.notna(row.get("confidence")):
                confidence = f"{tr('置信度', 'Confidence')}: {float(row.get('confidence')):.3f}"
            blocks.append(
                _result_card(
                    str(row.get("run_id", "NA")),
                    "Run",
                    [
                        f"{tr('来源数据库', 'Source database')}: analysis_runs",
                        f"{tr('样本数', 'Samples')}: 1",
                        f"{tr('参考版本', 'Reference version')}: {row.get('atlas_name', 'NA')} {row.get('build_version', '') or ''}",
                        f"{tr('更新时间', 'Updated')}: {str(row.get('created_at', ''))[:10] or 'NA'}",
                    ],
                    note=f"{tr('样本', 'Sample')}: {row.get('sample_id', 'NA')} | Top: {row.get('top_region', 'NA')}",
                    confidence=confidence,
                )
            )
        st.markdown(f'<div class="workbench-results-grid">{"".join(blocks)}</div>', unsafe_allow_html=True) if blocks else st.info(tr("没有 run 匹配。", "No run matches."))


def _render_global_workbench_search(compact: bool = False) -> None:
    title = tr("研究工作台检索", "Research Workspace Search")
    subtitle = tr(
        "输入样本 ID、基因名、脑区、cell type 或 run ID，直接定位到统一结果视图。",
        "Search sample IDs, genes, brain regions, cell types or run IDs from one entry point.",
    )
    render_panel_header(title, subtitle)
    query = st.text_input(
        tr("全局搜索", "Global search"),
        placeholder=tr("例如 RBFOX3 / GFAP / CTX_M1 / sample ID / run ID", "e.g. RBFOX3 / GFAP / CTX_M1 / sample ID / run ID"),
        key="workbench_global_search_compact" if compact else "workbench_global_search",
        label_visibility="collapsed",
    )
    if query.strip():
        _render_workbench_results(_build_workbench_search(query, get_database_mode()), query)


def _atlas_ids_for_mode(conn: sqlite3.Connection, db_mode: str) -> list[int]:
    try:
        df = pd.read_sql_query("SELECT atlas_id, species FROM atlas_versions", conn)
    except Exception:
        return []
    if df.empty:
        return [1] if db_mode == "rhesus" else []
    df = df[df["species"].apply(lambda x: matches_species(x, db_mode))].copy()
    if df.empty:
        return []
    return df["atlas_id"].astype(int).tolist()


def _build_dashboard_stats(db_mode: str) -> dict:
    stats = {
        "reference_samples": 0,
        "detected_genes": 0,
        "brain_regions": 0,
        "sequencing_runs": 0,
        "species": 0,
        "uploaded_samples": 0,
    }
    try:
        with sqlite3.connect(DB_PATH) as conn:
            cur = conn.cursor()
            stats["reference_samples"] = int(
                pd.read_sql_query("SELECT species FROM cfrna_samples", conn)["species"].apply(lambda x: matches_species(x, db_mode)).sum()
            ) if True else 0
            atlas_ids = _atlas_ids_for_mode(conn, db_mode)
            if atlas_ids:
                placeholders = ",".join(["?"] * len(atlas_ids))
                stats["detected_genes"] = int(
                    pd.read_sql_query(
                        f"SELECT COUNT(DISTINCT gene_symbol) AS n FROM reference_expression WHERE atlas_id IN ({placeholders})",
                        conn,
                        params=atlas_ids,
                    ).iloc[0]["n"]
                )
                stats["brain_regions"] = int(
                    pd.read_sql_query(
                        f"SELECT COUNT(DISTINCT region_id) AS n FROM reference_expression WHERE atlas_id IN ({placeholders})",
                        conn,
                        params=atlas_ids,
                    ).iloc[0]["n"]
                )
                v2_runs = int(
                    pd.read_sql_query(
                        f"SELECT COUNT(*) AS n FROM analysis_runs WHERE atlas_id IN ({placeholders})",
                        conn,
                        params=atlas_ids,
                    ).iloc[0]["n"]
                )
            else:
                stats["detected_genes"] = 0
                stats["brain_regions"] = 0
                v2_runs = 0
            legacy_runs = int(cur.execute("SELECT COUNT(*) FROM source_tracing_results").fetchone()[0]) if db_mode == "rhesus" else 0
            stats["sequencing_runs"] = legacy_runs + v2_runs
            stats["uploaded_samples"] = stats["reference_samples"]
            species_df = pd.read_sql_query("SELECT DISTINCT species FROM cfrna_samples WHERE species IS NOT NULL AND species != ''", conn)
            species_df = species_df[species_df["species"].apply(lambda x: matches_species(x, db_mode))]
            stats["species"] = int(len(species_df))
    except Exception:
        pass
    return stats


def _tissue_composition_df(db_mode: str) -> pd.DataFrame:
    atlas_filter = ""
    params = []
    try:
        with sqlite3.connect(DB_PATH) as conn:
            atlas_ids = _atlas_ids_for_mode(conn, db_mode)
        if atlas_ids:
            placeholders = ",".join(["?"] * len(atlas_ids))
            atlas_filter = f"WHERE atlas_id IN ({placeholders})"
            params = atlas_ids
        elif db_mode == "human":
            return pd.DataFrame()
    except Exception:
        if db_mode == "human":
            return pd.DataFrame()
    df = _safe_query(
        f"""
        SELECT region_id, COUNT(*) AS gene_records
        FROM reference_expression
        {atlas_filter}
        GROUP BY region_id
        ORDER BY gene_records DESC
        LIMIT 8
        """,
        params,
    )
    if df.empty:
        return df
    total = df["gene_records"].sum()
    df["fraction"] = df["gene_records"] / total
    return df


def _expression_overview_df(db_mode: str) -> pd.DataFrame:
    atlas_filter = ""
    params = []
    try:
        with sqlite3.connect(DB_PATH) as conn:
            atlas_ids = _atlas_ids_for_mode(conn, db_mode)
        if atlas_ids:
            placeholders = ",".join(["?"] * len(atlas_ids))
            atlas_filter = f"WHERE atlas_id IN ({placeholders})"
            params = atlas_ids
        elif db_mode == "human":
            return pd.DataFrame()
    except Exception:
        if db_mode == "human":
            return pd.DataFrame()
    return _safe_query(
        f"""
        SELECT gene_symbol, AVG(avg_tpm) AS avg_tpm, AVG(sample_count) AS sample_count
        FROM reference_expression
        {atlas_filter}
        GROUP BY gene_symbol
        HAVING avg_tpm IS NOT NULL
        LIMIT 4000
        """,
        params,
    )


def _qc_distribution_df() -> pd.DataFrame:
    df = _safe_query(
        """
        SELECT COALESCE(qc_status, 'Unknown') AS qc_status, COUNT(*) AS n
        FROM cfrna_samples
        GROUP BY COALESCE(qc_status, 'Unknown')
        """
    )
    if df.empty:
        return pd.DataFrame({"qc_status": ["Unknown"], "n": [0], "fraction": [0.0]})
    total = max(int(df["n"].sum()), 1)
    df["fraction"] = df["n"] / total
    return df


def _latest_updates() -> list[dict]:
    updates = []
    recent_runs = _safe_query(
        """
        SELECT run_id, sample_id, method, created_at
        FROM analysis_runs
        ORDER BY created_at DESC
        LIMIT 3
        """
    )
    for _, row in recent_runs.iterrows():
        updates.append(
            {
                "title": tr(f"样本 {str(row.get('sample_id', 'NA'))} 的分析已完成", f"Run for sample {str(row.get('sample_id', 'NA'))} completed"),
                "note": tr(
                    f"方法: {row.get('method', 'unknown')} | Run ID: {str(row.get('run_id', ''))[:8]}...",
                    f"Method: {row.get('method', 'unknown')} | Run ID: {str(row.get('run_id', ''))[:8]}...",
                ),
                "date": str(row.get("created_at", ""))[:10],
            }
        )

    if not updates:
        updates = [
            {
                "title": tr("已启用面向评估的 Dashboard", "Benchmark-ready dashboard enabled"),
                "note": tr("已统一 KPI 卡片、图表面板与快速检索入口。", "Unified KPI cards, chart panels and quick search entry points."),
                "date": "2026-04-27",
            },
            {
                "title": tr("支持队列化 QC 校准", "Cohort QC calibration support"),
                "note": tr("RBC、免疫背景和脑信号风险支持基于队列分布的校准。", "RBC, immune and brain-signal risks now support cohort-distribution calibration."),
                "date": "2026-04-27",
            },
            {
                "title": tr("Bo2023 浏览器缺表保护已增强", "Bo2023 browser fallback refined"),
                "note": tr("缺失 atlas 表时会显示友好提示，而不是原始报错。", "Missing atlas tables now return graceful hints instead of raw traceback."),
                "date": "2026-04-26",
            },
        ]
    return updates


def _overview_copy(db_mode: str) -> dict:
    if db_mode == "human":
        return {
            "subtitle": tr(
                "整合 HPA 与 Allen Human Brain Atlas 的人脑转录组参考图谱，用于基因表达检索、脑区对照、样本 QC 与后续溯源评估。",
                "A human brain transcriptome workspace integrating HPA and Allen Human Brain Atlas references for gene search, region comparison, sample QC and downstream tracing evaluation.",
            ),
            "glass_note": tr(
                "当前工作区聚焦 Homo sapiens 脑区转录组图谱。建议先查看 atlas 覆盖、基因表达分布和可用样本，再进入数据提交、图谱浏览和 Benchmark。",
                "This workspace focuses on Homo sapiens brain transcriptome references. Review atlas coverage, expression landscape and available samples before submission, atlas browsing and Benchmark workflows.",
            ),
            "species_icon": "HUM",
            "sample_note": tr("当前人脑模式下的入库样本", "Uploaded samples in human mode"),
        }
    return {
        "subtitle": tr(
            "用于血浆 cfRNA 溯源、脑区去卷积、Benchmark 评估与猕猴多组学图谱整合的综合平台。",
            "A comprehensive platform for plasma cfRNA tracing, brain-region deconvolution, benchmark evaluation, and macaque multi-omics atlas integration.",
        ),
        "glass_note": tr(
            "这个首页是整个平台的科研入口：先看整体状态，再进入数据提交、图谱浏览、溯源分析、Run 对比和论文级 Benchmark 解释。",
            "This dashboard acts as the scientific front door of the platform: overview first, then data submission, atlas browsing, tracing analysis, run comparison and publication-grade Benchmark interpretation.",
        ),
        "species_icon": "NHP",
        "sample_note": tr("当前猕猴血浆cfRNA损伤溯源数据库中的入库样本", "Uploaded samples in macaque plasma cfRNA tracing mode"),
    }


def _render_tissue_panel(df: pd.DataFrame) -> None:
    render_panel_header(
        tr("组织组成", "Tissue Composition"),
        tr("基于脑区参考表达记录的 atlas 组成概览。", "Reference atlas composition derived from region-level expression records."),
    )
    if df.empty:
        st.info(tr("当前没有可用的参考组成数据。", "No reference composition data is currently available."))
        return
    fig = px.pie(
        df,
        values="gene_records",
        names="region_id",
        hole=0.62,
        color_discrete_sequence=["#2f6df6", "#74a1ff", "#9bc0ff", "#7cd6c1", "#b4d9a7", "#d9c4ff", "#f5d597", "#d5dce6"],
    )
    fig.update_traces(textinfo="none", hovertemplate="%{label}<br>%{value} records<br>%{percent}")
    fig.update_layout(height=360, margin=dict(l=10, r=10, t=10, b=10), showlegend=True)
    st.plotly_chart(fig, use_container_width=True)


def _render_expression_panel(df: pd.DataFrame) -> None:
    render_panel_header(
        tr("基因表达概览", "Gene Expression Overview"),
        tr("查看参考表达强度与脑区覆盖范围之间的关系。", "Average reference expression versus regional coverage across genes."),
    )
    if df.empty:
        st.info(tr("当前参考表无法计算表达概览。", "No expression overview could be computed from the current reference tables."))
        return
    plot_df = df.copy()
    plot_df["avg_tpm"] = pd.to_numeric(plot_df["avg_tpm"], errors="coerce").fillna(0.0) + 1e-6
    plot_df["sample_count"] = pd.to_numeric(plot_df["sample_count"], errors="coerce").fillna(0.0)
    fig = px.scatter(
        plot_df,
        x="avg_tpm",
        y="sample_count",
        opacity=0.45,
        render_mode="webgl",
        color="sample_count",
        color_continuous_scale=["#dbe7ff", "#75a0ff", "#2f6df6"],
        labels={"avg_tpm": tr("平均 TPM", "Average TPM"), "sample_count": tr("脑区样本覆盖", "Regional sample count")},
        hover_data=["gene_symbol"],
    )
    fig.update_xaxes(type="log")
    fig.update_layout(height=360, margin=dict(l=10, r=10, t=10, b=10), coloraxis_showscale=False)
    st.plotly_chart(fig, use_container_width=True)


def _render_qc_panel(df: pd.DataFrame) -> None:
    render_panel_header(
        tr("样本质控", "Sample Quality Control"),
        tr("数据库内已上传 cfRNA 样本的当前 QC 风险分布。", "Current database-level QC distribution for uploaded plasma cfRNA samples."),
    )
    if df.empty:
        st.info(tr("当前还没有可用的 QC 分布数据。", "No QC distribution data is available yet."))
        return
    color_map = {
        "Low risk": "#1f9d75",
        "Moderate risk": "#f2b447",
        "High risk": "#d43f56",
        "Unknown": "#94a3b8",
        "Uncalibrated": "#7f90a5",
        "Pass": "#1f9d75",
        "Warning": "#f2b447",
        "Fail": "#d43f56",
    }
    fig = go.Figure()
    for _, row in df.iterrows():
        label = str(row["qc_status"])
        fig.add_trace(
            go.Bar(
                x=[row["fraction"]],
                y=[tr("QC 概况", "QC profile")],
                orientation="h",
                name=label,
                marker=dict(color=color_map.get(label, "#94a3b8")),
                hovertemplate=f"{label}: {int(row['n'])} samples (%{{x:.1%}})<extra></extra>",
            )
        )
    fig.update_layout(
        barmode="stack",
        height=240,
        margin=dict(l=10, r=10, t=10, b=10),
        xaxis=dict(range=[0, 1], tickformat=".0%"),
        yaxis=dict(showticklabels=False),
        legend=dict(orientation="h", yanchor="bottom", y=1.02, x=0),
    )
    st.plotly_chart(fig, use_container_width=True)
    mini = [{"title": str(row["qc_status"]), "note": f"{int(row['n'])} samples | {row['fraction']:.1%}"} for _, row in df.iterrows()]
    render_mini_cards(mini)


def _render_quick_search() -> None:
    _render_global_workbench_search(compact=True)


def display_database_overview():
    db_mode = get_database_mode()
    stats = _build_dashboard_stats(db_mode)
    copy = _overview_copy(db_mode)
    render_page_hero(
        tr(f"{database_label(db_mode)}", f"{database_label(db_mode)}"),
        copy["subtitle"],
        eyebrow=tr("仪表盘", "Dashboard"),
        pills=[
            tr("参考图谱", "Reference atlas"),
            tr("溯源分析", "Tracing analysis"),
            tr("Benchmark", "Benchmark"),
            tr("带 QC 的结果解释", "QC-aware interpretation"),
        ],
    )
    render_result_hint(copy["glass_note"])
    _render_global_workbench_search()

    render_kpi_cards(
        [
            {"icon": "REF", "label": tr("参考样本", "Reference Samples"), "value": f"{stats['reference_samples']:,}", "note": tr("当前入库样本总数", "Total uploaded samples")},
            {"icon": "DNA", "label": tr("检测基因", "Detected Genes"), "value": f"{stats['detected_genes']:,}", "note": tr("参考表达表中的去重基因", "Distinct genes in reference tables")},
            {"icon": "REG", "label": tr("脑区数", "Brain Regions"), "value": f"{stats['brain_regions']:,}", "note": tr("Atlas 覆盖脑区", "Atlas region coverage")},
            {"icon": "RUN", "label": tr("分析记录", "Sequencing Runs"), "value": f"{stats['sequencing_runs']:,}", "note": tr("可回溯的分析记录", "Available tracing records")},
            {"icon": copy["species_icon"], "label": tr("物种", "Species"), "value": f"{stats['species']:,}", "note": tr("当前数据库中的物种标签", "Species labels currently present")},
        ]
    )

    render_section_band(
        tr("数据库仪表盘", "Database Dashboard"),
        tr("统一查看参考组成、表达图景、样本 QC 状态和快捷入口。", "Reference composition, expression landscape, QC status and navigation shortcuts."),
    )
    col1, col2, col3 = st.columns([1.05, 1.25, 1.0])
    with col1:
        _render_tissue_panel(_tissue_composition_df(db_mode))
    with col2:
        _render_expression_panel(_expression_overview_df(db_mode))
    with col3:
        _render_qc_panel(_qc_distribution_df())

    render_result_hint(
        tr(
            "优先看 KPI、QC 分布和参考图谱覆盖范围；如果样本数或 atlas 覆盖不足，建议先完成数据上传或补充参考资源。",
            "Start with the KPI cards, QC distribution and atlas coverage. If sample count or atlas coverage is limited, expand the uploaded data or reference resources first.",
        )
    )

    render_section_band(
        tr("更新与快捷操作", "Updates & Quick Actions"),
        tr("在一个视图里查看最新变化、常用流程和首轮检索入口。", "Latest changes, common workflows and first-step searches in one view."),
    )
    col4, col5, col6 = st.columns([1.0, 1.05, 1.05])
    with col4:
        render_panel_header(
            tr("最近更新", "Latest Updates"),
            tr("在解读结果前，先确认最近运行和系统更新。", "Recent runs and system-level changes worth checking before interpretation."),
        )
        render_update_list(_latest_updates())
    with col5:
        render_panel_header(
            tr("常用分析入口", "Popular Analyses"),
            tr("快速进入最常见的科研使用流程。", "Shortcut entry points for the most common review workflows."),
        )
        quick_cols = st.columns(2)
        buttons = [
            (tr("溯源分析", "Tracing Analysis"), "tracing"),
            (tr("Benchmark", "Benchmark"), "benchmark"),
            (tr("参考图谱", "Reference Atlas"), "atlas"),
            (tr("Run 对比", "Run Compare"), "compare"),
        ]
        for idx, (label, target) in enumerate(buttons):
            with quick_cols[idx % 2]:
                if st.button(label, key=f"overview_jump_{target}", use_container_width=True):
                    st.session_state.page = target
                    st.rerun()
    with col6:
        _render_quick_search()
