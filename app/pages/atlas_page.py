from __future__ import annotations

import sqlite3
from typing import Iterable

import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import streamlit as st

from app.components.layout import render_kpi_cards, render_panel_header, render_section_band
from app.database_mode import DATABASE_MODES, database_label, get_database_mode, matches_species
from app.i18n import tr
from app.shared import DB_PATH, render_page_hero, render_result_hint
from data.dao import get_atlas_options


@st.cache_data(show_spinner=False)
def _load_table(query: str, params: tuple = ()) -> pd.DataFrame:
    conn = sqlite3.connect(DB_PATH, timeout=15)
    try:
        return pd.read_sql_query(query, conn, params=params)
    finally:
        conn.close()


@st.cache_data(show_spinner=False)
def _sqlite_object_exists(name: str) -> bool:
    conn = sqlite3.connect(DB_PATH, timeout=15)
    try:
        cur = conn.cursor()
        cur.execute("SELECT name FROM sqlite_master WHERE name=? AND type IN ('table', 'view')", (name,))
        return cur.fetchone() is not None
    finally:
        conn.close()


def _missing_bo2023_objects(names: Iterable[str]) -> list[str]:
    return [name for name in names if not _sqlite_object_exists(name)]


def _show_missing_objects(names: Iterable[str], label: str) -> bool:
    missing = _missing_bo2023_objects(names)
    if missing:
        st.info(
            tr(
                f"{label} 暂不可用：当前数据库缺少 {', '.join(missing)}。请先导入对应 Bo2023 数据层。",
                f"{label} is currently unavailable because the database is missing {', '.join(missing)}. Please import the required Bo2023 layer first.",
            )
        )
        return True
    return False


@st.cache_data(show_spinner=False)
def _get_bo2023_table_names() -> list[str]:
    df = _load_table("SELECT name FROM sqlite_master WHERE type='table' AND name LIKE 'bo2023_%' ORDER BY name")
    return df["name"].tolist() if not df.empty else []


@st.cache_data(show_spinner=False)
def _get_bo2023_views() -> list[str]:
    df = _load_table("SELECT name FROM sqlite_master WHERE type='view' AND name LIKE 'bo2023_%' ORDER BY name")
    return df["name"].tolist() if not df.empty else []


@st.cache_data(show_spinner=False)
def _load_bo2023_regions() -> pd.DataFrame:
    if _missing_bo2023_objects(["bo2023_regions", "bo2023_region_sample_coverage", "bo2023_region_qc_overview"]):
        return pd.DataFrame()
    return _load_table(
        """
        SELECT
            r.region,
            r.full_name,
            COALESCE(r.dictionary_lobe, c.lobe, r.lobe) AS lobe,
            r.regional_map,
            COALESCE(c.neocortexregion, r.neocortexregion) AS neocortexregion,
            c.total_sample_count,
            q.rin,
            q.uniquely_mapped_percent,
            q.pct_mrna_bases,
            q.pct_usable_bases,
            q.pct_correct_strand_reads
        FROM bo2023_regions r
        LEFT JOIN bo2023_region_sample_coverage c ON r.region = c.region
        LEFT JOIN bo2023_region_qc_overview q ON r.region = q.region
        ORDER BY r.region
        """
    )


@st.cache_data(show_spinner=False)
def _load_bo2023_qc() -> pd.DataFrame:
    if _missing_bo2023_objects(["bo2023_region_qc_overview"]):
        return pd.DataFrame()
    return _load_table(
        """
        SELECT region, lobe, neocortexregion, rin, uniquely_mapped_percent,
               pct_pf_reads_aligned, pct_mrna_bases, pct_usable_bases,
               pct_correct_strand_reads, at_dropout, gc_dropout
        FROM bo2023_region_qc_overview
        ORDER BY region
        """
    )


@st.cache_data(show_spinner=False)
def _load_bo2023_deg_summary() -> pd.DataFrame:
    if _missing_bo2023_objects(["bo2023_region_deg_summary"]):
        return pd.DataFrame()
    return _load_table("SELECT type, region, number, lobe, percentage FROM bo2023_region_deg_summary ORDER BY number DESC")


@st.cache_data(show_spinner=False)
def _load_bo2023_ct_catalog() -> pd.DataFrame:
    if _missing_bo2023_objects(["bo2023_ct_gene_catalog"]):
        return pd.DataFrame()
    return _load_table("SELECT gene_name, aggregate, age2, age3, age4, age5, age6, age7, age8 FROM bo2023_ct_gene_catalog ORDER BY aggregate DESC")


@st.cache_data(show_spinner=False)
def _search_ct_gene(keyword: str) -> pd.DataFrame:
    if _missing_bo2023_objects(["bo2023_ct_gene_catalog"]):
        return pd.DataFrame()
    kw = f"%{keyword.strip()}%"
    return _load_table(
        """
        SELECT gene_name, aggregate, age2, age3, age4, age5, age6, age7, age8
        FROM bo2023_ct_gene_catalog
        WHERE gene_name LIKE ?
        ORDER BY aggregate DESC
        LIMIT 100
        """,
        (kw,),
    )


@st.cache_data(show_spinner=False)
def _search_module_gene(keyword: str) -> pd.DataFrame:
    if _missing_bo2023_objects(["bo2023_wgcna_cortical_expressed_genes_and_modules"]):
        return pd.DataFrame()
    kw = f"%{keyword.strip()}%"
    return _load_table(
        """
        SELECT gene_name, geneid, gene_type, m20, m4
        FROM bo2023_wgcna_cortical_expressed_genes_and_modules
        WHERE gene_name LIKE ? OR geneid LIKE ?
        ORDER BY gene_name
        LIMIT 200
        """,
        (kw, kw),
    )


@st.cache_data(show_spinner=False)
def _load_custom_table_preview(table_name: str, limit: int = 200) -> pd.DataFrame:
    safe_name = "".join(ch for ch in table_name if ch.isalnum() or ch == "_")
    if not safe_name.startswith("bo2023_"):
        raise ValueError(tr("仅允许预览 bo2023_ 表或视图。", "Only bo2023_ tables or views can be previewed."))
    if not _sqlite_object_exists(safe_name):
        return pd.DataFrame()
    return _load_table(f"SELECT * FROM {safe_name} LIMIT {int(limit)}")


@st.cache_data(show_spinner=False)
def _load_bo2023_expression_atlases() -> pd.DataFrame:
    if _missing_bo2023_objects(["atlas_versions", "reference_expression"]):
        return pd.DataFrame()
    return _load_table(
        """
        SELECT a.atlas_id, a.atlas_name, a.build_version, a.normalization, COUNT(r.id) AS n_expression_rows
        FROM atlas_versions a
        JOIN reference_expression r ON r.atlas_id = a.atlas_id
        WHERE lower(COALESCE(a.atlas_name, '') || ' ' || COALESCE(a.build_version, '')) LIKE '%bo2023%'
           OR lower(COALESCE(a.atlas_name, '') || ' ' || COALESCE(a.build_version, '')) LIKE '%wang%'
        GROUP BY a.atlas_id, a.atlas_name, a.build_version, a.normalization
        ORDER BY a.atlas_id DESC
        """
    )


@st.cache_data(show_spinner=False)
def _load_imported_bo2023_region_overview(atlas_id: int) -> pd.DataFrame:
    return _load_table(
        """
        SELECT
            m.region_id,
            m.region_name,
            m.region_acronym,
            m.parent_region_id,
            COUNT(DISTINCT r.gene_symbol) AS n_genes,
            MAX(r.sample_count) AS sample_count
        FROM macaque_brain_atlas m
        LEFT JOIN reference_expression r
          ON r.atlas_id = m.atlas_id
         AND r.region_id = m.region_id
        WHERE m.atlas_id = ?
        GROUP BY m.region_id, m.region_name, m.region_acronym, m.parent_region_id
        ORDER BY m.region_id
        """,
        (int(atlas_id),),
    )


@st.cache_data(show_spinner=False)
def _load_bo2023_gene_expression(atlas_id: int, gene_symbol: str) -> pd.DataFrame:
    return _load_table(
        """
        SELECT region_id, region_name, avg_tpm, median_tpm, expression_class, sample_count
        FROM reference_expression
        WHERE atlas_id = ? AND upper(gene_symbol) = upper(?)
        ORDER BY avg_tpm DESC
        """,
        (int(atlas_id), gene_symbol.strip()),
    )


@st.cache_data(show_spinner=False)
def _load_bo2023_m1_s1_expression(atlas_id: int, gene_symbols: tuple[str, ...]) -> pd.DataFrame:
    if not gene_symbols:
        return pd.DataFrame()
    placeholders = ",".join(["?"] * len(gene_symbols))
    params = [int(atlas_id)] + [g.upper() for g in gene_symbols]
    return _load_table(
        f"""
        SELECT gene_symbol, region_id, region_name, avg_tpm
        FROM reference_expression
        WHERE atlas_id = ?
          AND upper(gene_symbol) IN ({placeholders})
          AND (
            upper(region_id) LIKE 'M1%' OR upper(region_name) LIKE '%M1%'
            OR upper(region_id) LIKE 'S1%' OR upper(region_name) LIKE '%S1%'
          )
        ORDER BY gene_symbol, region_id
        """,
        tuple(params),
    )


@st.cache_data(show_spinner=False)
def _load_atlas_catalog() -> pd.DataFrame:
    if _missing_bo2023_objects(["atlas_versions"]):
        return pd.DataFrame()
    return _load_table(
        """
        SELECT atlas_id, atlas_name, species, level, build_version, gene_id_type,
               normalization, created_at, notes
        FROM atlas_versions
        ORDER BY atlas_id DESC
        """
    )


@st.cache_data(show_spinner=False)
def _load_atlas_region_options(atlas_id: int) -> pd.DataFrame:
    return _load_table(
        """
        SELECT
            region_id,
            MAX(region_name) AS region_name,
            COUNT(DISTINCT gene_symbol) AS n_genes,
            MAX(sample_count) AS sample_count
        FROM reference_expression
        WHERE atlas_id = ?
        GROUP BY region_id
        ORDER BY region_id
        """,
        (int(atlas_id),),
    )


@st.cache_data(show_spinner=False)
def _load_atlas_celltype_options(atlas_id: int) -> pd.DataFrame:
    return _load_table(
        """
        SELECT cell_type_marker, COUNT(DISTINCT gene_symbol) AS n_genes
        FROM reference_expression
        WHERE atlas_id = ?
          AND cell_type_marker IS NOT NULL
          AND cell_type_marker != ''
        GROUP BY cell_type_marker
        ORDER BY n_genes DESC, cell_type_marker
        LIMIT 200
        """,
        (int(atlas_id),),
    )


def _expression_filter_sql(region_ids: list[str], celltype: str | None) -> tuple[str, list]:
    clauses = []
    params: list = []
    if region_ids:
        placeholders = ",".join(["?"] * len(region_ids))
        clauses.append(f"region_id IN ({placeholders})")
        params.extend(region_ids)
    if celltype:
        clauses.append("cell_type_marker = ?")
        params.append(celltype)
    return (" AND " + " AND ".join(clauses) if clauses else ""), params


@st.cache_data(show_spinner=False)
def _load_atlas_region_ranking(atlas_id: int, region_ids: tuple[str, ...], celltype: str | None) -> pd.DataFrame:
    extra_sql, extra_params = _expression_filter_sql(list(region_ids), celltype)
    return _load_table(
        f"""
        SELECT
            region_id,
            MAX(region_name) AS region_name,
            COUNT(DISTINCT gene_symbol) AS gene_count,
            AVG(avg_tpm) AS mean_tpm,
            MAX(sample_count) AS sample_count,
            SUM(CASE WHEN expression_class IS NOT NULL AND expression_class != '' THEN 1 ELSE 0 END) AS classified_genes
        FROM reference_expression
        WHERE atlas_id = ?
        {extra_sql}
        GROUP BY region_id
        ORDER BY mean_tpm DESC, gene_count DESC
        LIMIT 80
        """,
        tuple([int(atlas_id)] + extra_params),
    )


@st.cache_data(show_spinner=False)
def _load_atlas_gene_candidates(atlas_id: int, region_ids: tuple[str, ...], celltype: str | None, limit: int = 40) -> list[str]:
    extra_sql, extra_params = _expression_filter_sql(list(region_ids), celltype)
    df = _load_table(
        f"""
        SELECT gene_symbol, AVG(avg_tpm) AS mean_tpm
        FROM reference_expression
        WHERE atlas_id = ?
        {extra_sql}
        GROUP BY gene_symbol
        HAVING mean_tpm IS NOT NULL
        ORDER BY mean_tpm DESC
        LIMIT ?
        """,
        tuple([int(atlas_id)] + extra_params + [int(limit)]),
    )
    if df.empty:
        return []
    return df["gene_symbol"].astype(str).tolist()


@st.cache_data(show_spinner=False)
def _load_atlas_expression_matrix(
    atlas_id: int,
    region_ids: tuple[str, ...],
    gene_symbols: tuple[str, ...],
    celltype: str | None,
) -> pd.DataFrame:
    genes = [g.strip().upper() for g in gene_symbols if g.strip()]
    extra_sql, extra_params = _expression_filter_sql(list(region_ids), celltype)
    gene_sql = ""
    gene_params: list = []
    if genes:
        placeholders = ",".join(["?"] * len(genes))
        gene_sql = f" AND upper(gene_symbol) IN ({placeholders})"
        gene_params = genes
    return _load_table(
        f"""
        SELECT gene_symbol, region_id, region_name, avg_tpm, median_tpm,
               sample_count, expression_class, cell_type_marker
        FROM reference_expression
        WHERE atlas_id = ?
        {extra_sql}
        {gene_sql}
        ORDER BY region_id, gene_symbol
        LIMIT 5000
        """,
        tuple([int(atlas_id)] + extra_params + gene_params),
    )


@st.cache_data(show_spinner=False)
def _load_marker_evidence(atlas_id: int, region_ids: tuple[str, ...], celltype: str | None, limit: int = 200) -> pd.DataFrame:
    sigset = _load_table(
        """
        SELECT sigset_id, method, topk_per_region, created_at
        FROM signature_sets
        WHERE atlas_id = ?
        ORDER BY sigset_id DESC
        LIMIT 1
        """,
        (int(atlas_id),),
    )
    if not sigset.empty:
        sigset_id = int(sigset.iloc[0]["sigset_id"])
        region_clause = ""
        params: list = [sigset_id, int(atlas_id)]
        if region_ids:
            placeholders = ",".join(["?"] * len(region_ids))
            region_clause = f" AND sg.region_id IN ({placeholders})"
            params.extend(region_ids)
        celltype_clause = ""
        if celltype:
            celltype_clause = " AND r.cell_type_marker = ?"
            params.append(celltype)
        return _load_table(
            f"""
            SELECT
                sg.region_id,
                MAX(r.region_name) AS region_name,
                sg.gene_symbol,
                sg.weight,
                MAX(r.avg_tpm) AS avg_tpm,
                MAX(r.expression_class) AS expression_class,
                MAX(r.cell_type_marker) AS cell_type_marker
            FROM signature_genes sg
            LEFT JOIN reference_expression r
              ON r.atlas_id = ?
             AND r.region_id = sg.region_id
             AND r.gene_symbol = sg.gene_symbol
            WHERE sg.sigset_id = ?
            {region_clause}
            {celltype_clause}
            GROUP BY sg.region_id, sg.gene_symbol, sg.weight
            ORDER BY sg.region_id, sg.weight DESC, avg_tpm DESC
            LIMIT ?
            """,
            tuple([int(atlas_id), sigset_id] + params[2:] + [int(limit)]),
        )

    extra_sql, extra_params = _expression_filter_sql(list(region_ids), celltype)
    return _load_table(
        f"""
        SELECT region_id, region_name, gene_symbol, 1.0 AS weight, avg_tpm,
               expression_class, cell_type_marker
        FROM reference_expression
        WHERE atlas_id = ?
        {extra_sql}
        ORDER BY avg_tpm DESC
        LIMIT ?
        """,
        tuple([int(atlas_id)] + extra_params + [int(limit)]),
    )


def _parse_gene_input(text: str) -> tuple[str, ...]:
    cleaned = str(text or "").replace("\n", ",").replace(";", ",")
    genes = []
    for item in cleaned.split(","):
        value = item.strip()
        if value:
            genes.append(value)
    return tuple(dict.fromkeys(genes))


def _render_marker_overlap(marker_df: pd.DataFrame, query_genes: tuple[str, ...]) -> None:
    if marker_df.empty:
        st.info(tr("当前筛选下没有可用 marker 证据。", "No marker evidence is available for the current filters."))
        return
    overlap = marker_df.copy()
    query_set = {g.upper() for g in query_genes}
    if query_set:
        overlap["overlap"] = overlap["gene_symbol"].astype(str).str.upper().isin(query_set)
        summary = overlap.groupby("region_id", as_index=False).agg(
            marker_count=("gene_symbol", "nunique"),
            overlap_count=("overlap", "sum"),
        )
    else:
        summary = overlap.groupby("region_id", as_index=False).agg(marker_count=("gene_symbol", "nunique"))
        summary["overlap_count"] = 0
    fig = px.bar(
        summary.sort_values(["overlap_count", "marker_count"], ascending=False).head(40),
        x="region_id",
        y=["marker_count", "overlap_count"],
        barmode="group",
        title=tr("Marker overlap / evidence coverage", "Marker overlap / evidence coverage"),
        color_discrete_sequence=["#7aa2ff", "#1f9d75"],
    )
    fig.update_layout(height=360, margin=dict(l=10, r=10, t=60, b=10))
    st.plotly_chart(fig, use_container_width=True)


def _render_atlas_explorer(default_mode: str) -> None:
    catalog = _load_atlas_catalog()
    if catalog.empty:
        st.info(tr("当前数据库没有可浏览的 atlas_versions。", "No browsable atlas_versions are available in the current database."))
        return

    species_labels = {
        "rhesus": tr("猴 / macaque", "Macaque"),
        "human": tr("人 / human", "Human"),
    }
    mode_options = ["rhesus", "human"]
    default_index = mode_options.index(default_mode) if default_mode in mode_options else 0
    left, right = st.columns([0.92, 2.35], gap="large")

    with left:
        render_panel_header(
            tr("筛选器", "Filters"),
            tr("按物种、脑区、细胞类型和数据来源收窄 atlas 视图。", "Narrow the atlas view by species, region, cell type and source."),
        )
        selected_mode = st.radio(
            tr("物种", "Species"),
            options=mode_options,
            index=default_index,
            format_func=lambda m: species_labels[m],
            key="atlas_explorer_species",
        )
        mode_catalog = catalog[catalog["species"].apply(lambda x: matches_species(x, selected_mode))].copy()
        if mode_catalog.empty:
            st.warning(tr("该物种下没有可用 atlas。", "No atlas is available for this species."))
            return

        source_options = [tr("全部来源", "All sources")] + sorted(mode_catalog["atlas_name"].dropna().astype(str).unique().tolist())
        source = st.selectbox(tr("数据来源", "Data source"), source_options, key="atlas_explorer_source")
        if source != tr("全部来源", "All sources"):
            mode_catalog = mode_catalog[mode_catalog["atlas_name"].astype(str) == source].copy()

        atlas_labels = [
            f"{int(r.atlas_id)} | {r.atlas_name} | {r.build_version} | {r.level}"
            for r in mode_catalog.itertuples()
        ]
        atlas_label = st.selectbox(tr("参考 atlas", "Reference atlas"), atlas_labels, key="atlas_explorer_atlas")
        atlas_id = int(atlas_label.split(" | ")[0])
        atlas_meta = mode_catalog[mode_catalog["atlas_id"].astype(int) == atlas_id].iloc[0].to_dict()

        regions = _load_atlas_region_options(atlas_id)
        region_options = regions["region_id"].astype(str).tolist() if not regions.empty else []
        selected_regions = st.multiselect(
            tr("脑区", "Brain regions"),
            region_options,
            default=region_options[: min(8, len(region_options))],
            key="atlas_explorer_regions",
        )

        celltypes = _load_atlas_celltype_options(atlas_id)
        celltype_options = [tr("全部细胞类型", "All cell types")]
        if not celltypes.empty:
            celltype_options += celltypes["cell_type_marker"].astype(str).tolist()
        celltype_label = st.selectbox(tr("细胞类型", "Cell type"), celltype_options, key="atlas_explorer_celltype")
        celltype = None if celltype_label == tr("全部细胞类型", "All cell types") else celltype_label

        gene_text = st.text_area(
            tr("基因列表", "Gene list"),
            placeholder=tr("逗号分隔，例如 RBFOX3, GFAP, MBP", "Comma-separated, e.g. RBFOX3, GFAP, MBP"),
            height=86,
            key="atlas_explorer_genes",
        )
        top_gene_count = st.slider(tr("自动补充 Top genes", "Auto-fill top genes"), 10, 60, 30, step=10, key="atlas_explorer_top_genes")

    query_genes = _parse_gene_input(gene_text)
    region_tuple = tuple(selected_regions)
    if not query_genes:
        query_genes = tuple(_load_atlas_gene_candidates(atlas_id, region_tuple, celltype, top_gene_count))
    matrix_df = _load_atlas_expression_matrix(atlas_id, region_tuple, query_genes, celltype)
    ranking_df = _load_atlas_region_ranking(atlas_id, region_tuple, celltype)
    marker_df = _load_marker_evidence(atlas_id, region_tuple, celltype)

    with right:
        render_panel_header(
            tr("矩阵 / 图谱视图", "Matrix / Atlas View"),
            tr("主视图包含表达热图、脑区排名和 marker overlap。", "The main view combines expression heatmap, region ranking and marker overlap."),
        )
        render_kpi_cards(
            [
                {"icon": "ATL", "label": tr("Atlas", "Atlas"), "value": str(atlas_id), "note": str(atlas_meta.get("atlas_name", ""))},
                {"icon": "REG", "label": tr("脑区", "Regions"), "value": f"{len(region_tuple) or (regions['region_id'].nunique() if not regions.empty else 0):,}", "note": tr("当前筛选覆盖", "Current selection")},
                {"icon": "GEN", "label": tr("基因", "Genes"), "value": f"{len(query_genes):,}", "note": tr("热图候选基因", "Heatmap candidate genes")},
            ]
        )

        chart_tabs = st.tabs([tr("热图", "Heatmap"), tr("脑区排名", "Region ranking"), tr("Marker overlap", "Marker overlap")])
        with chart_tabs[0]:
            if matrix_df.empty:
                st.info(tr("当前筛选下没有表达矩阵记录。", "No expression matrix records match the current filters."))
            else:
                heat = matrix_df.pivot_table(index="region_id", columns="gene_symbol", values="avg_tpm", aggfunc="mean").fillna(0)
                fig = go.Figure(
                    data=go.Heatmap(
                        z=heat.values,
                        x=heat.columns.tolist(),
                        y=heat.index.tolist(),
                        colorscale="Blues",
                        colorbar=dict(title="avg TPM"),
                    )
                )
                fig.update_layout(height=max(420, min(900, 24 * len(heat.index) + 120)), margin=dict(l=10, r=10, t=20, b=10))
                st.plotly_chart(fig, use_container_width=True)
        with chart_tabs[1]:
            if ranking_df.empty:
                st.info(tr("当前筛选下没有脑区排名。", "No region ranking is available for the current filters."))
            else:
                fig = px.bar(
                    ranking_df.sort_values("mean_tpm", ascending=True).tail(35),
                    x="mean_tpm",
                    y="region_id",
                    color="gene_count",
                    orientation="h",
                    hover_data=["region_name", "sample_count", "classified_genes"],
                    color_continuous_scale=["#dbe7ff", "#7aa2ff", "#2f6df6"],
                    title=tr("脑区平均表达排名", "Region ranking by mean expression"),
                )
                fig.update_layout(height=520, margin=dict(l=10, r=10, t=60, b=10), coloraxis_colorbar=dict(title="genes"))
                st.plotly_chart(fig, use_container_width=True)
        with chart_tabs[2]:
            _render_marker_overlap(marker_df, query_genes)

    render_section_band(
        tr("详情、证据与导出", "Details, Evidence & Export"),
        tr("查看基因列表、marker 证据、引用信息，并下载当前筛选结果。", "Inspect gene lists, marker evidence, citations and download the current selection."),
    )
    detail_tabs = st.tabs([tr("基因列表", "Gene list"), tr("Marker 证据", "Marker evidence"), tr("引用 / 版本", "References / Version"), tr("下载", "Download")])
    with detail_tabs[0]:
        if matrix_df.empty:
            st.info(tr("没有可显示的基因表达明细。", "No gene-expression details are available."))
        else:
            st.dataframe(matrix_df, use_container_width=True, hide_index=True)
    with detail_tabs[1]:
        if marker_df.empty:
            st.info(tr("没有可显示的 marker 证据。", "No marker evidence is available."))
        else:
            st.dataframe(marker_df, use_container_width=True, hide_index=True)
    with detail_tabs[2]:
        ref_df = pd.DataFrame([atlas_meta])
        st.dataframe(ref_df, use_container_width=True, hide_index=True)
        notes = str(atlas_meta.get("notes") or "").strip()
        if notes:
            st.info(notes)
    with detail_tabs[3]:
        st.download_button(
            tr("下载表达矩阵 CSV", "Download expression matrix CSV"),
            matrix_df.to_csv(index=False).encode("utf-8-sig"),
            file_name=f"atlas_{atlas_id}_expression_selection.csv",
            mime="text/csv",
            disabled=matrix_df.empty,
        )
        st.download_button(
            tr("下载 marker 证据 CSV", "Download marker evidence CSV"),
            marker_df.to_csv(index=False).encode("utf-8-sig"),
            file_name=f"atlas_{atlas_id}_marker_evidence.csv",
            mime="text/csv",
            disabled=marker_df.empty,
        )


def _display_imported_bo2023_region_browser() -> None:
    atlases = _load_bo2023_expression_atlases()
    if atlases.empty:
        st.info(tr("尚未检测到已入库的 Bo2023 / Wang 表达 atlas。", "No imported Bo2023 / Wang expression atlas was detected."))
        return
    labels = [f"{int(r.atlas_id)} | {r.atlas_name} | {r.build_version} | rows={int(r.n_expression_rows):,}" for r in atlases.itertuples()]
    selected = st.selectbox(tr("选择已入库 atlas", "Choose imported atlas"), labels, key="bo2023_region_atlas")
    atlas_id = int(selected.split(" | ")[0])
    overview = _load_imported_bo2023_region_overview(atlas_id)
    if overview.empty:
        st.info(tr("当前 atlas 没有可浏览的脑区注释。", "The selected atlas has no browsable region annotation."))
        return

    render_kpi_cards(
        [
            {"icon": "REG", "label": tr("可浏览脑区", "Browsable regions"), "value": f"{overview['region_id'].nunique():,}", "note": tr("来自已入库 Bo2023 atlas", "From imported Bo2023 atlas")},
            {"icon": "GENE", "label": tr("覆盖基因", "Covered genes"), "value": f"{int(overview['n_genes'].max()):,}", "note": tr("每个脑区的 gene-symbol 覆盖", "Gene-symbol coverage per region")},
            {"icon": "SMP", "label": tr("最大样本数", "Max samples"), "value": f"{int(overview['sample_count'].max()):,}", "note": tr("单个脑区聚合样本数", "Samples aggregated in one region")},
        ]
    )
    render_result_hint(
        tr(
            "这里展示的是已经由 819 个样本聚合得到的 Bo2023_WangLab_VSD_region atlas，可直接浏览 110 个猕猴脑区。",
            "This table shows the imported Bo2023_WangLab_VSD_region atlas aggregated from 819 samples; 110 macaque brain regions are directly browsable.",
        )
    )
    keyword = st.text_input(tr("筛选脑区", "Filter regions"), placeholder=tr("例如 F1, V1, thalamus", "e.g. F1, V1, thalamus"), key="bo2023_region_filter")
    view = overview.copy()
    if keyword.strip():
        q = keyword.strip().lower()
        view = view[
            view["region_id"].astype(str).str.lower().str.contains(q)
            | view["region_name"].astype(str).str.lower().str.contains(q)
            | view["parent_region_id"].astype(str).str.lower().str.contains(q)
        ]
    st.dataframe(view, use_container_width=True, hide_index=True)
    if not view.empty:
        fig = px.bar(
            view.sort_values("sample_count", ascending=False).head(40),
            x="region_id",
            y="sample_count",
            color="parent_region_id",
            hover_data=["region_name", "n_genes"],
            title=tr("Bo2023 脑区样本覆盖", "Bo2023 regional sample coverage"),
            color_discrete_sequence=px.colors.qualitative.Set3,
        )
        fig.update_layout(height=460, margin=dict(l=10, r=10, t=60, b=10))
        st.plotly_chart(fig, use_container_width=True)


def _display_legacy_atlas() -> None:
    render_section_band(tr("参考图谱", "Reference Atlas"), tr("浏览 legacy atlas 的脑区注释和参考表达。", "Browse legacy atlas annotations and reference expression."))
    render_result_hint(
        tr(
            "建议先确认脑区 ID、全称和层级关系，再输入基因查看该基因在不同脑区中的参考表达差异。",
            "Review region IDs, full names and hierarchy first, then inspect how a gene varies across regions.",
        )
    )
    regions_df = _load_table(
        """
        SELECT region_id, region_name, region_acronym, parent_region_id
        FROM macaque_brain_atlas
        ORDER BY region_id
        """
    )
    ref_count_df = _load_table("SELECT COUNT(*) AS n FROM reference_expression")
    ref_rows = int(ref_count_df.iloc[0]["n"]) if not ref_count_df.empty else 0
    render_kpi_cards(
        [
            {"icon": "REG", "label": tr("脑区数", "Brain Regions"), "value": f"{regions_df['region_id'].nunique() if not regions_df.empty else 0}", "note": tr("legacy atlas 脑区", "Legacy atlas regions")},
            {"icon": "REF", "label": tr("参考记录", "Reference Rows"), "value": f"{ref_rows:,}", "note": tr("参考表达记录数", "Reference expression records")},
        ]
    )
    st.dataframe(regions_df, use_container_width=True, hide_index=True)

    render_panel_header(tr("Legacy 基因检索", "Legacy Gene Search"), tr("查看某个基因在 legacy atlas 各脑区中的分布。", "Inspect a gene across legacy atlas regions."))
    gene_input = st.text_input(tr("输入基因符号", "Enter gene symbol"), placeholder=tr("例如 GAD1, SLC17A7", "e.g. GAD1, SLC17A7"), key="legacy_gene_input")
    if gene_input:
        gene_df = _load_table(
            """
            SELECT gene_symbol, gene_name, region_id, region_name, avg_tpm
            FROM reference_expression
            WHERE upper(gene_symbol) = upper(?)
            ORDER BY avg_tpm DESC
            """,
            (gene_input.strip(),),
        )
        if gene_df.empty:
            st.warning(tr(f"未找到基因 {gene_input} 的参考表达数据。", f"No reference-expression data was found for {gene_input}."))
        else:
            st.dataframe(gene_df, use_container_width=True, hide_index=True)
            fig = px.bar(gene_df, x="region_id", y="avg_tpm", color="avg_tpm", title=tr(f"{gene_input} 在 legacy atlas 中的平均表达", f"Average expression of {gene_input} across the legacy atlas"), color_continuous_scale=["#eef4ff", "#8fb5ff", "#2f6df6"])
            fig.update_layout(height=420, margin=dict(l=10, r=10, t=60, b=10))
            st.plotly_chart(fig, use_container_width=True)


def _display_bo2023_expression_browser() -> None:
    atlases = _load_bo2023_expression_atlases()
    render_panel_header(tr("Bo2023 表达浏览器", "Bo2023 Expression Browser"), tr("查询已导入的 Bo2023 / Wang gene-by-region 表达矩阵。", "Query imported Bo2023 / Wang gene-by-region matrices."))
    if atlases.empty:
        st.info(tr("尚未检测到可用的 Bo2023 gene × region 表达矩阵。请先导入重建矩阵；注释层浏览仍可正常使用。", "No Bo2023 gene × region expression matrix was detected. Please import a reconstructed matrix first; annotation-layer browsing still works."))
        return
    atlas_labels = [f"{int(r.atlas_id)} | {r.atlas_name} | {r.build_version} | rows={int(r.n_expression_rows):,}" for r in atlases.itertuples()]
    selected = st.selectbox(tr("选择表达 atlas", "Choose expression atlas"), atlas_labels, key="bo2023_expr_atlas")
    atlas_id = int(selected.split(" | ")[0])
    gene = st.text_input(tr("按 gene 查询脑区表达", "Query gene expression by region"), value="GAD1", key="bo2023_expr_gene").strip()
    if not gene:
        return
    expr = _load_bo2023_gene_expression(atlas_id, gene)
    if expr.empty:
        st.warning(tr(f"当前表达 atlas 中未找到 gene: {gene}", f"The current expression atlas does not contain gene: {gene}"))
        return
    render_kpi_cards(
        [
            {"icon": "REG", "label": tr("覆盖脑区", "Covered regions"), "value": f"{expr['region_id'].nunique()}", "note": tr("该基因有表达记录的脑区数", "Number of regions with expression records")},
            {"icon": "MAX", "label": tr("最高 avg TPM", "Maximum avg TPM"), "value": f"{float(expr['avg_tpm'].max()):.2f}", "note": tr("最强表达脑区", "Strongest-expression region")},
            {"icon": "MED", "label": tr("中位 avg TPM", "Median avg TPM"), "value": f"{float(expr['avg_tpm'].median()):.2f}", "note": tr("整体表达水平", "Overall expression level")},
        ]
    )
    render_result_hint(tr("先看 Top region 排名，再看热图确认这个基因是否具有清晰脑区偏好。", "Review the top region ranking first, then use the heatmap to confirm whether the gene has a clear regional preference."))
    topn = st.slider(tr("显示 Top N region", "Show Top N regions"), 5, min(80, len(expr)), min(20, len(expr)), key="bo2023_expr_topn")
    top_expr = expr.head(topn)
    col1, col2 = st.columns(2)
    with col1:
        fig = px.bar(top_expr, x="avg_tpm", y="region_id", color="expression_class", orientation="h", hover_data=["region_name"], title=f"{gene} region ranking", color_discrete_sequence=px.colors.qualitative.Set2)
        fig.update_layout(height=max(420, 24 * len(top_expr)), yaxis={"categoryorder": "total ascending"})
        st.plotly_chart(fig, use_container_width=True)
    with col2:
        heat = expr.copy()
        heat["gene_symbol"] = gene.upper()
        fig2 = px.density_heatmap(heat, x="gene_symbol", y="region_id", z="avg_tpm", title=f"{gene} expression heatmap")
        fig2.update_layout(height=max(420, 16 * len(heat)), margin=dict(l=10, r=10, t=60, b=10))
        st.plotly_chart(fig2, use_container_width=True)
    st.dataframe(expr, use_container_width=True, hide_index=True)

    st.markdown(f"#### {tr('M1 vs S1 精细比较', 'M1 vs S1 fine comparison')}")
    gene_list = st.text_input(tr("M1 / S1 比较基因（逗号分隔）", "M1 / S1 comparison genes (comma-separated)"), value=gene, key="bo2023_m1s1_genes")
    genes = tuple(g.strip() for g in gene_list.split(",") if g.strip())
    m1s1 = _load_bo2023_m1_s1_expression(atlas_id, genes)
    if m1s1.empty:
        st.info(tr("未找到匹配的 M1 / S1 表达记录，请确认矩阵中的 region 注释包含 M1 / S1 标识。", "No matching M1 / S1 records were found. Please confirm that the matrix annotations include M1 / S1 labels."))
    else:
        fig3 = px.bar(m1s1, x="region_id", y="avg_tpm", color="gene_symbol", barmode="group", title=tr("M1 vs S1 表达比较", "M1 vs S1 expression comparison"))
        fig3.update_layout(height=420, margin=dict(l=10, r=10, t=60, b=10))
        st.plotly_chart(fig3, use_container_width=True)
        st.dataframe(m1s1, use_container_width=True, hide_index=True)


def _display_human_atlas_browser() -> None:
    atlases = _load_bo2023_expression_atlases()
    human_opts = get_atlas_options(DB_PATH, species_mode="human")
    if not human_opts:
        st.info(
            tr(
                "当前数据库尚未导入 human atlas。请先执行人脑参考导入脚本。",
                "No human atlas has been imported into the current database yet. Import a human reference atlas first.",
            )
        )
        return
    render_section_band(
        tr("人脑图谱浏览器", "Human Atlas Browser"),
        tr(
            "使用与猕猴血浆cfRNA损伤溯源数据库一致的交互方式浏览 human brain reference atlas。",
            "Browse the human brain reference atlas using the same interaction pattern as the macaque plasma cfRNA tracing database.",
        ),
    )
    labels = [label for _, label in human_opts]
    ids = [atlas_id for atlas_id, _ in human_opts]
    atlas_choice = st.selectbox(tr("选择人脑 atlas", "Choose human atlas"), labels, index=0, key="human_atlas_choice")
    atlas_id = ids[labels.index(atlas_choice)]
    gene = st.text_input(tr("查询基因", "Query gene"), value="GFAP", key="human_gene_query").strip()
    if not gene:
        return
    expr = _load_bo2023_gene_expression(atlas_id, gene)
    if expr.empty:
        st.info(tr(f"当前 human atlas 中未找到 {gene}。", f"{gene} was not found in the current human atlas."))
        return
    render_kpi_cards(
        [
            {"icon": "REG", "label": tr("覆盖脑区", "Covered regions"), "value": expr["region_id"].nunique(), "note": tr("该基因覆盖的脑区数", "Number of covered regions")},
            {"icon": "MAX", "label": tr("最高 avg TPM", "Maximum avg TPM"), "value": f"{float(expr['avg_tpm'].max()):.2f}", "note": tr("最强表达脑区", "Strongest region")},
            {"icon": "MED", "label": tr("中位 avg TPM", "Median avg TPM"), "value": f"{float(expr['avg_tpm'].median()):.2f}", "note": tr("整体表达水平", "Overall level")},
        ]
    )
    col1, col2 = st.columns(2)
    with col1:
        fig = px.bar(expr.head(20), x="avg_tpm", y="region_id", color="expression_class", orientation="h", title=tr(f"{gene} 脑区排名", f"{gene} region ranking"))
        fig.update_layout(height=max(420, 24 * min(20, len(expr))), yaxis={"categoryorder": "total ascending"})
        st.plotly_chart(fig, use_container_width=True)
    with col2:
        heat = expr.copy()
        heat["gene_symbol"] = gene.upper()
        fig2 = px.density_heatmap(heat, x="gene_symbol", y="region_id", z="avg_tpm", title=tr(f"{gene} 表达热图", f"{gene} expression heatmap"))
        fig2.update_layout(height=max(420, 16 * len(heat)), margin=dict(l=10, r=10, t=60, b=10))
        st.plotly_chart(fig2, use_container_width=True)
    st.dataframe(expr, use_container_width=True, hide_index=True)


def _display_bo2023_browser() -> None:
    if not _get_bo2023_table_names():
        st.warning(tr("当前数据库中未检测到 bo2023_* 表。", "No bo2023_* tables were detected in the current database."))
        return

    regions_df = _load_bo2023_regions()
    qc_df = _load_bo2023_qc()
    deg_df = _load_bo2023_deg_summary()
    ct_df = _load_bo2023_ct_catalog()

    render_section_band(tr("Bo2023 图谱浏览器", "Bo2023 Atlas Browser"), tr("浏览脑区注释、QC、DEG、细胞类型相关基因和表达矩阵。", "Browse region annotation, QC, DEG summaries, cell-type-linked genes and expression matrices."))
    render_kpi_cards(
        [
            {"icon": "REG", "label": tr("Bo2023 脑区", "Bo2023 Regions"), "value": f"{regions_df['region'].nunique() if not regions_df.empty else 0}", "note": tr("已注释脑区", "Annotated regions")},
            {"icon": "LOB", "label": tr("脑叶类别", "Lobes"), "value": f"{regions_df['lobe'].nunique() if not regions_df.empty else 0}", "note": tr("脑叶分类数", "Lobe categories")},
            {"icon": "GEN", "label": tr("CT 基因", "CT Genes"), "value": f"{ct_df['gene_name'].nunique() if not ct_df.empty else 0}", "note": tr("细胞类型相关基因", "Cell-type related genes")},
            {"icon": "OBJ", "label": tr("对象数", "Objects"), "value": f"{len(_get_bo2023_table_names()) + len(_get_bo2023_views())}", "note": tr("Bo2023 表和视图", "Bo2023 tables and views")},
        ]
    )
    render_result_hint(tr("根据问题选择标签页：脑区总览看注释与覆盖，QC 看质量，表达浏览器看表达偏好，全表预览用于核对导入内容。", "Pick tabs by task: use Region Overview for annotation and coverage, QC for quality, Expression Browser for regional preference, and Table Preview for import auditing."))

    tab1, tab2, tab3, tab4, tab5 = st.tabs([tr("脑区总览", "Region Overview"), tr("QC 与覆盖", "QC & Coverage"), tr("基因 / 模块检索", "Gene / Module Search"), tr("表达浏览器", "Expression Browser"), tr("全表预览", "Table Preview")])

    with tab1:
        render_panel_header(
            tr("已入库 Bo2023 表达 atlas 脑区", "Imported Bo2023 expression-atlas regions"),
            tr("直接浏览由 Wang lab 819 个样本聚合得到的 gene x region atlas。", "Browse the gene-by-region atlas aggregated from 819 Wang-lab samples."),
        )
        _display_imported_bo2023_region_browser()
        if not _show_missing_objects(["bo2023_regions", "bo2023_region_sample_coverage", "bo2023_region_qc_overview"], "Region Overview"):
            render_panel_header(tr("脑区卡片", "Region Cards"), tr("浏览 Bo2023 脑区及其样本覆盖和基础 QC。", "Browse Bo2023 regions with coverage and basic QC fields."))
            lobe_options = [tr("全部", "All")] + sorted([x for x in regions_df["lobe"].dropna().unique().tolist() if str(x).strip()])
            selected_lobe = st.selectbox(tr("按 lobe 筛选", "Filter by lobe"), lobe_options, key="bo2023_lobe")
            filtered = regions_df.copy()
            if selected_lobe != tr("全部", "All"):
                filtered = filtered[filtered["lobe"] == selected_lobe]
            st.dataframe(filtered, use_container_width=True, hide_index=True)
            if not filtered.empty:
                fig = px.bar(filtered.sort_values("total_sample_count", ascending=False), x="region", y="total_sample_count", color="lobe", hover_data=["full_name", "regional_map", "rin"], title=tr("脑区样本覆盖", "Regional sample coverage"), color_discrete_sequence=px.colors.qualitative.Set3)
                fig.update_layout(height=460, margin=dict(l=10, r=10, t=60, b=10))
                st.plotly_chart(fig, use_container_width=True)

    with tab2:
        if not _show_missing_objects(["bo2023_region_qc_overview"], "QC & Coverage"):
            render_panel_header(tr("QC 概览", "QC Overview"), tr("查看不同脑叶的 region-level 质量指标。", "Inspect region-level quality metrics across lobes."))
            qc_lobes = [tr("全部", "All")] + sorted([x for x in qc_df["lobe"].dropna().unique().tolist() if str(x).strip()])
            qc_lobe = st.selectbox(tr("QC lobe 筛选", "QC lobe filter"), qc_lobes, key="bo2023_qc_lobe")
            qc_view = qc_df.copy()
            if qc_lobe != tr("全部", "All"):
                qc_view = qc_view[qc_view["lobe"] == qc_lobe]
            st.dataframe(qc_view, use_container_width=True, hide_index=True)
            metric_name = st.selectbox(tr("选择 QC 指标", "Choose QC metric"), ["rin", "uniquely_mapped_percent", "pct_pf_reads_aligned", "pct_mrna_bases", "pct_usable_bases", "pct_correct_strand_reads", "at_dropout", "gc_dropout"], key="bo2023_qc_metric")
            fig = px.box(qc_view, x="lobe", y=metric_name, color="neocortexregion", points="all", title=f"{metric_name} by lobe")
            fig.update_layout(height=460, margin=dict(l=10, r=10, t=60, b=10))
            st.plotly_chart(fig, use_container_width=True)

    with tab3:
        render_panel_header(tr("基因 / 模块检索", "Gene / Module Search"), tr("搜索 CT 基因、模块基因和 DEG 汇总。", "Search CT-linked genes, module genes and DEG summaries."))
        gene_keyword = st.text_input(tr("搜索基因（CT gene / module gene）", "Search gene (CT gene / module gene)"), placeholder=tr("例如 FEZF2, GAD1, CRYM", "e.g. FEZF2, GAD1, CRYM"), key="bo2023_gene_kw")
        c1, c2 = st.columns(2)
        if gene_keyword:
            with c1:
                ct_hit = _search_ct_gene(gene_keyword)
                st.markdown(f"#### {tr('CT 相关基因', 'CT-related genes')}")
                st.dataframe(ct_hit, use_container_width=True, hide_index=True)
                if not ct_hit.empty:
                    long_ct = ct_hit.melt(id_vars=["gene_name"], value_vars=[c for c in ct_hit.columns if c.startswith("age")], var_name="age", value_name="score")
                    fig = px.line(long_ct, x="age", y="score", color="gene_name", markers=True, title=tr("年龄相关 CT 基因变化", "Age-dependent CT gene profile"))
                    fig.update_layout(height=380, margin=dict(l=10, r=10, t=60, b=10))
                    st.plotly_chart(fig, use_container_width=True)
            with c2:
                module_hit = _search_module_gene(gene_keyword)
                st.markdown(f"#### {tr('模块 / WGCNA 命中', 'Module / WGCNA hits')}")
                st.dataframe(module_hit, use_container_width=True, hide_index=True)
        else:
            st.info(tr("输入基因名后，可同时检索 CT 基因目录和 WGCNA 模块基因。", "Enter a gene name to search both the CT-gene catalog and WGCNA module genes."))

        st.markdown(f"#### {tr('脑区 DEG 概览', 'Regional DEG overview')}")
        if deg_df.empty:
            st.info(tr("当前数据库中没有 bo2023_region_deg_summary。", "The current database does not contain bo2023_region_deg_summary."))
        else:
            deg_type = st.selectbox(tr("DEG 类型", "DEG type"), sorted(deg_df["type"].dropna().unique().tolist()), key="bo2023_deg_type")
            top_n = st.slider(tr("显示 Top N 脑区", "Show Top N regions"), 10, 60, 20, step=5, key="bo2023_deg_n")
            deg_view = deg_df[deg_df["type"] == deg_type].sort_values("number", ascending=False).head(top_n)
            st.dataframe(deg_view, use_container_width=True, hide_index=True)
            fig = px.bar(deg_view, x="region", y="number", color="lobe", title=f"{deg_type} DEG count by region", color_discrete_sequence=px.colors.qualitative.Pastel)
            fig.update_layout(height=420, margin=dict(l=10, r=10, t=60, b=10))
            st.plotly_chart(fig, use_container_width=True)

    with tab4:
        _display_bo2023_expression_browser()

    with tab5:
        render_panel_header(tr("Bo2023 全表预览", "Bo2023 Table Preview"), tr("预览导入的 bo2023_* 表并下载当前切片。", "Preview imported bo2023_* tables and download the current slice."))
        all_sources = _get_bo2023_views() + _get_bo2023_table_names()
        if not all_sources:
            st.info(tr("当前数据库中没有可预览的 bo2023_* 表或视图。", "There are no previewable bo2023_* tables or views in the current database."))
        else:
            selected_table = st.selectbox(tr("选择 bo2023_* 表或视图", "Choose bo2023_* table or view"), all_sources, key="bo2023_custom_table")
            preview_limit = st.slider(tr("预览行数", "Preview rows"), 20, 500, 100, 20, key="bo2023_preview_limit")
            preview_df = _load_custom_table_preview(selected_table, preview_limit)
            st.dataframe(preview_df, use_container_width=True, hide_index=True)
            csv_bytes = preview_df.to_csv(index=False).encode("utf-8-sig")
            st.download_button(tr("下载当前预览 CSV", "Download current preview CSV"), csv_bytes, file_name=f"{selected_table}_preview.csv", mime="text/csv")


def _display_human_atlas_browser_v2() -> None:
    human_opts = get_atlas_options(DB_PATH, species_mode="human")
    if not human_opts:
        st.info(
            tr(
                "当前数据库尚未导入 human atlas。请先执行人脑参考图谱导入脚本。",
                "No human atlas has been imported into the current database yet. Import a human reference atlas first.",
            )
        )
        return

    render_section_band(
        tr("人脑图谱浏览器", "Human Atlas Browser"),
        tr(
            "浏览 HPA 与 Allen Human Brain Atlas 的 gene-by-region 表达矩阵，支持基因查询、脑区排名和表达热图。",
            "Browse HPA and Allen Human Brain Atlas gene-by-region matrices with gene search, region ranking and expression heatmaps.",
        ),
    )
    labels = [label for _, label in human_opts]
    ids = [atlas_id for atlas_id, _ in human_opts]
    atlas_choice = st.selectbox(tr("选择人脑 atlas", "Choose human atlas"), labels, index=0, key="human_atlas_choice_v2")
    atlas_id = ids[labels.index(atlas_choice)]
    gene = st.text_input(tr("查询基因", "Query gene"), value="GFAP", key="human_gene_query_v2").strip()
    if not gene:
        return

    expr = _load_bo2023_gene_expression(atlas_id, gene)
    if expr.empty:
        st.info(tr(f"当前 human atlas 中未找到 {gene}。", f"{gene} was not found in the current human atlas."))
        return

    render_kpi_cards(
        [
            {"icon": "REG", "label": tr("覆盖脑区", "Covered regions"), "value": expr["region_id"].nunique(), "note": tr("该基因覆盖的脑区数", "Number of covered regions")},
            {"icon": "MAX", "label": tr("最高 avg TPM", "Maximum avg TPM"), "value": f"{float(expr['avg_tpm'].max()):.2f}", "note": tr("最强表达脑区", "Strongest region")},
            {"icon": "MED", "label": tr("中位 avg TPM", "Median avg TPM"), "value": f"{float(expr['avg_tpm'].median()):.2f}", "note": tr("整体表达水平", "Overall level")},
        ]
    )
    render_result_hint(
        tr(
            "优先查看 region ranking 判断最高表达脑区，再用热图确认该基因是否具有人脑区域偏好。",
            "Start with the region ranking to identify the strongest regions, then use the heatmap to confirm whether the gene has a human brain regional preference.",
        )
    )

    col1, col2 = st.columns(2)
    with col1:
        fig = px.bar(
            expr.head(20),
            x="avg_tpm",
            y="region_id",
            color="expression_class",
            orientation="h",
            title=tr(f"{gene} 脑区排名", f"{gene} region ranking"),
        )
        fig.update_layout(height=max(420, 24 * min(20, len(expr))), yaxis={"categoryorder": "total ascending"})
        st.plotly_chart(fig, use_container_width=True)
    with col2:
        heat = expr.copy()
        heat["gene_symbol"] = gene.upper()
        fig2 = px.density_heatmap(
            heat,
            x="gene_symbol",
            y="region_id",
            z="avg_tpm",
            title=tr(f"{gene} 表达热图", f"{gene} expression heatmap"),
        )
        fig2.update_layout(height=max(420, 16 * len(heat)), margin=dict(l=10, r=10, t=60, b=10))
        st.plotly_chart(fig2, use_container_width=True)
    st.dataframe(expr, use_container_width=True, hide_index=True)


def display_atlas_browser() -> None:
    db_mode = get_database_mode()
    current_view = st.session_state.get("atlas_view", "legacy")
    if db_mode == "human":
        render_page_hero(
            tr(f"{database_label(db_mode)} - 人脑参考图谱", f"{database_label(db_mode)} - Human Reference Atlas"),
            tr(
                "浏览 Homo sapiens 脑区转录组参考图谱，支持 HPA 与 Allen Human Brain Atlas 的基因表达检索和脑区对照。",
                "Browse Homo sapiens brain transcriptome references with HPA and Allen Human Brain Atlas gene-expression search and region comparison.",
            ),
            eyebrow="Atlas",
            pills=[
                tr("HPA Brain", "HPA Brain"),
                tr("Allen HBA", "Allen HBA"),
                tr("基因检索", "Gene search"),
                tr("脑区表达", "Regional expression"),
            ],
        )
        _render_atlas_explorer("human")
        return

    render_page_hero(
        tr(f"{database_label(db_mode)} - Atlas Explorer", f"{database_label(db_mode)} - Atlas Explorer"),
        tr(
            "以工作台方式浏览猴 / 人脑参考图谱：左侧筛选物种、脑区、细胞类型和参考来源，右侧查看表达热图、脑区排名和 marker overlap，下方核对证据并导出。",
            "Explore macaque and human brain references as a workspace: filter species, regions, cell types and source on the left; inspect heatmaps, region ranking and marker overlap on the right; review evidence and export below.",
        ),
        eyebrow="Atlas",
        pills=[
            tr("左侧筛选", "Left filters"),
            tr("表达热图", "Expression heatmap"),
            tr("脑区排名", "Region ranking"),
            tr("Marker overlap", "Marker overlap"),
            tr("证据导出", "Evidence export"),
        ],
    )
    _render_atlas_explorer("rhesus")

    if current_view != "legacy":
        with st.expander(tr("Bo2023 深度浏览 / 导入审计", "Bo2023 deep browser / import audit"), expanded=False):
            _display_bo2023_browser()
