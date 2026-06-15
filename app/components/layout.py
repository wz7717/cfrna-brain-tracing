from __future__ import annotations

from html import escape
from typing import Iterable

import streamlit as st

from app.i18n import tr


def render_top_toolbar(current_page: str) -> None:
    st.markdown(
        f'<div class="top-toolbar">'
        f'<div class="top-toolbar-meta">'
        f'<span>{escape(tr("cfRNA 脑损伤溯源数据库", "cfRNA Brain Injury Source Tracing Database"))}</span>'
        f'<span>|</span>'
        f'<span>{escape(current_page)}</span>'
        f'</div>'
        f'<div class="top-toolbar-links">'
        f'<span class="top-toolbar-chip">{escape(tr("首页", "Home"))}</span>'
        f'<span class="top-toolbar-chip">{escape(tr("帮助", "Help"))}</span>'
        f'<span class="top-toolbar-chip">{escape(tr("下载", "Download"))}</span>'
        f'<span class="top-toolbar-chip">{escape(tr("关于", "About"))}</span>'
        f'</div>'
        f'</div>',
        unsafe_allow_html=True,
    )


def render_kpi_cards(cards: Iterable[dict]) -> None:
    blocks = []
    for card in cards:
        icon = escape(str(card.get("icon", "•")))
        label = escape(str(card.get("label", "")))
        value = escape(str(card.get("value", "")))
        blocks.append(
            f'<div class="kpi-card">'
            f'<div class="kpi-card-top">'
            f'<div class="kpi-card-icon">{icon}</div>'
            f'<div class="kpi-card-body">'
            f'<div class="kpi-card-label">{label}</div>'
            f'<div class="kpi-card-value">{value}</div>'
            f'</div>'
            f'</div>'
            f'</div>'
        )
    if blocks:
        st.markdown(f'<section class="kpi-grid">{"".join(blocks)}</section>', unsafe_allow_html=True)


def render_section_band(title: str, note: str = "") -> None:
    st.markdown(
        f'<div class="section-band"><div class="section-band-title">{escape(title)}</div></div>',
        unsafe_allow_html=True,
    )


def render_panel_header(title: str, subtitle: str = "") -> None:
    st.markdown(f'<div class="panel-title">{escape(title)}</div>', unsafe_allow_html=True)


def render_mini_cards(cards: Iterable[dict]) -> None:
    blocks = []
    for card in cards:
        title = escape(str(card.get("title", "")))
        blocks.append(f'<div class="mini-card"><div class="mini-card-title">{title}</div></div>')
    if blocks:
        st.markdown(f'<div class="mini-grid">{"".join(blocks)}</div>', unsafe_allow_html=True)


def render_update_list(items: Iterable[dict]) -> None:
    blocks = []
    for item in items:
        title = escape(str(item.get("title", "")))
        date = escape(str(item.get("date", "")))
        blocks.append(
            f'<div class="update-item">'
            f'<div class="update-dot"></div>'
            f'<div><div class="update-item-title">{title}</div></div>'
            f'<div class="update-item-date">{date}</div>'
            f'</div>'
        )
    if blocks:
        st.markdown(f'<div class="update-list">{"".join(blocks)}</div>', unsafe_allow_html=True)
