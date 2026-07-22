from __future__ import annotations

import streamlit as st

from app.i18n import tr


DATABASE_MODES = {
    "rhesus": {
        "zh": "BrainTrace",
        "en": "BrainTrace",
        "species_tokens": ["rhesus", "macaca", "macaque", "macaca mulatta", "macaca fascicularis"],
        "default_species": "Macaca mulatta",
    },
    "human": {
        "zh": "人脑转录组图谱",
        "en": "Human Brain Transcriptome Atlas",
        "species_tokens": ["human", "homo sapiens"],
        "default_species": "Homo sapiens",
    },
}


def get_database_mode() -> str:
    return st.session_state.get("database_mode", "rhesus")


def set_database_mode(mode: str) -> None:
    st.session_state["database_mode"] = mode if mode in DATABASE_MODES else "rhesus"


def database_label(mode: str | None = None) -> str:
    mode = mode or get_database_mode()
    meta = DATABASE_MODES.get(mode, DATABASE_MODES["rhesus"])
    return tr(meta["zh"], meta["en"])


def species_tokens(mode: str | None = None) -> list[str]:
    mode = mode or get_database_mode()
    return list(DATABASE_MODES.get(mode, DATABASE_MODES["rhesus"]).get("species_tokens", []))


def matches_species(value: str | None, mode: str | None = None) -> bool:
    text = str(value or "").strip().lower()
    if not text:
        return False
    return any(token in text for token in species_tokens(mode))


def default_species(mode: str | None = None) -> str:
    mode = mode or get_database_mode()
    return str(DATABASE_MODES.get(mode, DATABASE_MODES["rhesus"]).get("default_species", "Macaca mulatta"))


def mode_ready_message(mode: str | None = None) -> str:
    mode = mode or get_database_mode()
    if mode == "human":
        return tr(
            "当前已切换到人脑转录组图谱模式。页面会优先展示 Homo sapiens atlas、样本和分析结果。",
            "The workspace is now in Human Brain Transcriptome Atlas mode and prioritizes Homo sapiens atlases, samples and analysis outputs.",
        )
    return tr(
        "当前为 BrainTrace 工作区。",
        "Current workspace: BrainTrace.",
    )
