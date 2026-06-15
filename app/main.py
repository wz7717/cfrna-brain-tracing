from __future__ import annotations

import streamlit as st

from app.components.layout import render_top_toolbar
from app.database_mode import DATABASE_MODES, database_label, mode_ready_message, set_database_mode
from app.i18n import get_language_mode, set_language_mode, tr
from app.pages.atlas_page import display_atlas_browser
from app.pages.benchmark_page import display_benchmark_page
from app.pages.compare_runs_page import display_run_compare
from app.pages.overview_page import display_database_overview
from app.pages.sample_manage_page import display_sample_list
from app.pages.tracing_page import display_source_tracing
from app.pages.upload_page import display_data_upload
from app.shared import DB_PATH, init_database, inject_global_style, render_startup_check_summary
from data.dao import get_system_metrics


PAGES = {
    "overview": {"zh": "数据总览", "en": "Data Overview", "icon": "DB", "section": "OVERVIEW", "func": display_database_overview},
    "upload": {"zh": "数据提交", "en": "Data Submission", "icon": "UP", "section": "DATA BROWSER", "func": display_data_upload},
    "samples": {"zh": "样本管理", "en": "Sample Management", "icon": "SM", "section": "DATA BROWSER", "func": display_sample_list},
    "tracing": {"zh": "溯源分析", "en": "Tracing Analysis", "icon": "TR", "section": "ANALYSIS", "func": display_source_tracing},
    "compare": {"zh": "Run 对比", "en": "Run Comparison", "icon": "RC", "section": "ANALYSIS", "func": display_run_compare},
    "benchmark": {"zh": "性能评估", "en": "Performance", "icon": "BM", "section": "BENCHMARK", "func": display_benchmark_page},
    "atlas": {"zh": "参考图谱", "en": "Reference Atlas", "icon": "AT", "section": "ATLAS", "func": display_atlas_browser},
    "bo2023": {"zh": "Bo2023 图谱浏览器", "en": "Bo2023 Atlas Browser", "icon": "BO", "section": "ATLAS", "func": display_atlas_browser},
}

NAV_ORDER = ["OVERVIEW", "DATA BROWSER", "ANALYSIS", "BENCHMARK", "ATLAS"]

SECTION_LABELS = {
    "OVERVIEW": {"zh": "数据总览", "en": "Overview"},
    "DATA BROWSER": {"zh": "数据浏览", "en": "Data Browser"},
    "ANALYSIS": {"zh": "分析", "en": "Analysis"},
    "BENCHMARK": {"zh": "性能评估", "en": "Benchmark"},
    "ATLAS": {"zh": "图谱", "en": "Atlas"},
}


def _render_sidebar() -> None:
    with st.sidebar:
        st.markdown(f"### {tr('显示语言', 'Language')}")
        current_language = get_language_mode()
        lang_col_zh, lang_col_en = st.columns(2)
        with lang_col_zh:
            if st.button("中文", key="language_zh", use_container_width=True, type="primary" if current_language == "zh" else "secondary"):
                set_language_mode("zh")
                st.rerun()
        with lang_col_en:
            if st.button("English", key="language_en", use_container_width=True, type="primary" if current_language == "en" else "secondary"):
                set_language_mode("en")
                st.rerun()

        current_db_mode = st.session_state.get("database_mode", "rhesus")
        st.markdown(f"### {tr('数据库选择', 'Database Workspace')}")
        for mode_key in ["rhesus", "human"]:
            meta = DATABASE_MODES[mode_key]
            if st.button(
                tr(meta["zh"], meta["en"]),
                key=f"db_mode_{mode_key}",
                use_container_width=True,
                type="primary" if current_db_mode == mode_key else "secondary",
            ):
                set_database_mode(mode_key)
                st.session_state.atlas_view = "legacy"
                st.session_state.page = "overview"
                st.rerun()
        st.caption(mode_ready_message(current_db_mode))

        current_page = st.session_state.get("page", "overview")
        for section in NAV_ORDER:
            section_pages = [(page_key, meta) for page_key, meta in PAGES.items() if meta["section"] == section]
            is_current_section = any(page_key == current_page for page_key, _ in section_pages)
            section_label = tr(SECTION_LABELS[section]["zh"], SECTION_LABELS[section]["en"])
            if len(section_pages) == 1:
                page_key, meta = section_pages[0]
                if st.button(
                    section_label,
                    key=f"nav_section_{section}",
                    use_container_width=True,
                    type="primary" if is_current_section else "secondary",
                ):
                    st.session_state.page = page_key
                    st.rerun()
                continue

            with st.expander(section_label, expanded=is_current_section):
                for page_key, meta in section_pages:
                    is_active = page_key == current_page
                    if st.button(
                        f'{meta["icon"]}  {tr(meta["zh"], meta["en"])}',
                        key=f"nav_{page_key}",
                        use_container_width=True,
                        type="primary" if is_active else "secondary",
                    ):
                        st.session_state.page = page_key
                        if page_key == "atlas":
                            st.session_state.atlas_view = "legacy"
                        elif page_key == "bo2023":
                            st.session_state.atlas_view = "bo2023"
                        st.rerun()

        st.markdown("---")
        st.markdown(f"### {tr('系统概览', 'System Snapshot')}")
        st.caption(f"{tr('当前数据库', 'Current database')}: {database_label()}")
        try:
            metrics = get_system_metrics(DB_PATH)
        except Exception as exc:
            metrics = {}
            st.warning(f"{tr('系统指标读取失败', 'System metrics unavailable')}: {exc}")
        st.metric(tr("样本数", "Samples"), metrics.get("n_samples", 0))
        st.metric(tr("分析数", "Analyses"), metrics.get("n_analyses", 0))
        render_startup_check_summary(expanded=False)
        st.markdown(
            f"""
            <div class="sidebar-footnote">
                <strong>{tr("推荐流程", "Recommended workflow")}</strong><br>
                {tr(
                    "建议从 Dashboard 开始，先查看样本与数据库状态，再上传矩阵、运行溯源分析，最后用 Benchmark 和 atlas 解释结果。",
                    "Start from Dashboard, review samples and database status, then upload matrices, run tracing analysis, and finish with Benchmark plus atlas-based interpretation.",
                )}
            </div>
            """,
            unsafe_allow_html=True,
        )


def main() -> None:
    st.set_page_config(
        page_title="cfRNA Brain Injury Source Tracing Database",
        page_icon=":brain:",
        layout="wide",
        initial_sidebar_state="expanded",
    )
    inject_global_style()

    try:
        init_database()
    except Exception as exc:
        st.error(f"{tr('数据库初始化失败', 'Database initialization failed')}: {exc}")
        st.info(
            tr(
                "请检查 SQLite 路径、文件权限以及项目依赖是否完整。",
                "Please check the SQLite path, file permissions, and project dependencies.",
            )
        )
        return

    if "page" not in st.session_state or st.session_state.page not in PAGES:
        st.session_state.page = "overview"

    _render_sidebar()
    current = PAGES[st.session_state.page]
    render_top_toolbar(tr(current["zh"], current["en"]))
    current["func"]()


if __name__ == "__main__":
    main()
