from __future__ import annotations

from importlib import import_module

import streamlit as st

from app.i18n import get_language_mode, set_language_mode, tr
from app.shared import init_database, inject_global_style


PAGES = {
    "overview": {
        "zh": "溯源分析",
        "en": "Tracing Analysis",
        "func": "app.pages.tracing_page:display_source_tracing",
    },
    "upload": {
        "zh": "样本上传",
        "en": "Sample Upload",
        "func": "app.pages.upload_page:display_data_upload",
    },
    "samples": {
        "zh": "样本管理",
        "en": "Sample Management",
        "func": "app.pages.sample_manage_page:display_sample_list",
    },
}
NAV_ORDER = ["upload", "samples", "overview"]


def resolve_page_func(func_path: str):
    module_name, func_name = func_path.split(":", 1)
    return getattr(import_module(module_name), func_name)


def _set_page(page_key: str) -> None:
    if page_key in PAGES:
        st.session_state["page"] = page_key


def _render_sidebar() -> None:
    with st.sidebar:
        st.markdown("### BrainTrace")
        current_language = get_language_mode()
        lang_col_zh, lang_col_en = st.columns(2)
        with lang_col_zh:
            st.button(
                "中文",
                key="language_zh",
                width="stretch",
                type="primary" if current_language == "zh" else "secondary",
                on_click=set_language_mode,
                args=("zh",),
            )
        with lang_col_en:
            st.button(
                "English",
                key="language_en",
                width="stretch",
                type="primary" if current_language == "en" else "secondary",
                on_click=set_language_mode,
                args=("en",),
            )

        st.divider()
        st.caption(tr("工作流程", "Workflow"))
        current_page = st.session_state.get("page", "overview")
        for page_key in NAV_ORDER:
            meta = PAGES[page_key]
            st.button(
                tr(meta["zh"], meta["en"]),
                key=f"nav_{page_key}",
                width="stretch",
                type="primary" if current_page == page_key else "secondary",
                on_click=_set_page,
                args=(page_key,),
            )

        st.divider()
        st.caption("Network Top3 → resolution group Top3 → exact-region exploratory Top3")


def main() -> None:
    st.set_page_config(
        page_title="BrainTrace",
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
    resolve_page_func(current["func"])()


if __name__ == "__main__":
    main()
