from __future__ import annotations

import streamlit as st


LANGUAGE_OPTIONS = {
    "zh": "中文",
    "en": "English",
}


def get_language_mode() -> str:
    mode = st.session_state.get("ui_language_mode", "zh")
    return mode if mode in LANGUAGE_OPTIONS else "zh"


def set_language_mode(mode: str) -> None:
    st.session_state["ui_language_mode"] = mode if mode in LANGUAGE_OPTIONS else "zh"


def tr(zh: str, en: str) -> str:
    mode = get_language_mode()
    if mode == "zh":
        return zh
    if mode == "en":
        return en
    return zh
