from __future__ import annotations

from pathlib import Path
from typing import Dict
from html import escape
import os
import sqlite3
import types
import pandas as pd
import streamlit as st
from app.i18n import tr
from database_init import CSFRNASourceDatabase
from source_tracing import CSFRNASourceTracer
from data_processor import DataProcessor
from data.migrations import run_migrations
from data.bo2023_buildkit import import_buildkit_dir
from data.system_check import run_system_check

PROJECT_ROOT = Path(__file__).resolve().parents[1]
DB_PATH = str(Path(os.environ.get("CFRNA_DB_PATH", PROJECT_ROOT / "cfrna_source_tracing.db")).resolve())
CSS_PATH = PROJECT_ROOT / "app" / "styles" / "main.css"
GLOBAL_STYLE = """
<style>
:root {
    --cf-bg: #f5f8fc;
    --cf-surface: rgba(255, 255, 255, 0.9);
    --cf-surface-soft: #eef4fb;
    --cf-border: #dce6f2;
    --cf-text: #18395f;
    --cf-muted: #6b7a90;
    --cf-teal: #1f7aff;
    --cf-teal-dark: #1668e3;
    --cf-blue: #1f7aff;
    --cf-navy: #12385c;
    --cf-navy-dark: #0f2d4a;
    --cf-amber: #f2a93b;
    --cf-red: #b42318;
    --cf-shadow: 0 18px 40px rgba(24, 57, 95, 0.08);
}

.stApp {
    background:
        radial-gradient(circle at 52% 48%, rgba(255,255,255,.92), rgba(245,248,252,.82) 34%, rgba(239,245,253,.9) 100%),
        var(--cf-bg);
    color: var(--cf-text);
    font-family: "Inter", "PingFang SC", "Microsoft YaHei", Arial, sans-serif;
}

.block-container {
    padding-top: 1.3rem;
    padding-bottom: 3rem;
    max-width: 1560px;
}

section[data-testid="stSidebar"] {
    background: linear-gradient(180deg, var(--cf-navy) 0%, #102f4e 46%, #0b2946 100%);
    border-right: 1px solid rgba(255,255,255,.08);
}

section[data-testid="stSidebar"] * {
    color: #dcecff !important;
}

section[data-testid="stSidebar"] [data-testid="stMetric"] {
    background: rgba(255,255,255,.08);
    border: 1px solid rgba(209,229,250,.16);
    border-radius: 8px;
    padding: 0.75rem;
    box-shadow: none;
}

section[data-testid="stSidebar"] hr {
    border-color: rgba(255,255,255,.08);
}

section[data-testid="stSidebar"] label,
section[data-testid="stSidebar"] p,
section[data-testid="stSidebar"] span {
    color: #dcecff !important;
}

section[data-testid="stSidebar"] [role="radiogroup"] label {
    background: transparent;
    border: 1px solid transparent;
    border-radius: 8px;
    padding: 0.45rem 0.65rem;
    margin-bottom: 0.35rem;
    box-shadow: none;
}

section[data-testid="stSidebar"] [role="radiogroup"] label:hover {
    background: rgba(255,255,255,.08);
    border-color: rgba(209,229,250,.24);
    transform: translateX(1px);
}

section[data-testid="stSidebar"] [role="radiogroup"] label[data-checked="true"],
section[data-testid="stSidebar"] [role="radiogroup"] label:has(input:checked) {
    background: linear-gradient(135deg, #2384ff, #1d6ff2);
    border-color: transparent;
    box-shadow: 0 10px 22px rgba(20,102,221,.28);
    color: #ffffff !important;
}

section[data-testid="stSidebar"] h1,
section[data-testid="stSidebar"] h2,
section[data-testid="stSidebar"] h3 {
    color: #f3f8ff !important;
}

.main-header {
    font-size: 2.15rem;
    color: #12365f;
    text-align: left;
    margin-bottom: 1rem;
    font-weight: 900;
    letter-spacing: 0.2px;
}

.sub-header {
    font-size: 1.45rem;
    color: #12365f;
    margin-top: 1.7rem;
    margin-bottom: 1rem;
    padding: 0.9rem 1.05rem;
    border-left: 6px solid var(--cf-blue);
    border-bottom: 1px solid var(--cf-border);
    background: linear-gradient(90deg, rgba(31,122,255,.1) 0%, rgba(255,255,255,0) 100%);
    border-radius: 10px;
    font-weight: 820;
}

.page-hero {
    position: relative;
    overflow: hidden;
    margin: 0.15rem 0 1.1rem 0;
    padding: 1.2rem 1.25rem 1.1rem 1.25rem;
    border-radius: 12px;
    border: 1px solid rgba(205, 220, 240, 0.95);
    background:
        linear-gradient(135deg, rgba(255,255,255,.94) 0%, rgba(244,248,253,.88) 58%, rgba(236,243,252,.92) 100%);
    box-shadow: 0 18px 36px rgba(24,57,95,.08);
    backdrop-filter: blur(12px);
}

.page-hero::before {
    content: "";
    position: absolute;
    inset: 0 auto 0 0;
    width: 6px;
    background: linear-gradient(180deg, #3f97ff 0%, #1f7aff 55%, #1358cb 100%);
}

.page-hero::after {
    content: "";
    position: absolute;
    right: -84px;
    top: -84px;
    width: 230px;
    height: 230px;
    border-radius: 999px;
    background: radial-gradient(circle, rgba(31,122,255,.16) 0%, rgba(31,122,255,.05) 42%, rgba(31,122,255,0) 76%);
    pointer-events: none;
}

.page-hero-eyebrow {
    display: inline-flex;
    align-items: center;
    gap: 0.4rem;
    margin-bottom: 0.35rem;
    color: #52759d;
    font-size: 0.82rem;
    font-weight: 800;
    text-transform: uppercase;
    letter-spacing: 0.08em;
}

.page-hero-title {
    margin: 0;
    color: #12365f;
    font-size: 2rem;
    line-height: 1.12;
    font-weight: 900;
    letter-spacing: 0;
}

.page-hero-subtitle {
    margin: 0.55rem 0 0 0;
    max-width: 980px;
    color: #5f7591;
    font-size: 1rem;
    line-height: 1.7;
}

.page-hero-pills {
    display: flex;
    flex-wrap: wrap;
    gap: 0.55rem;
    margin-top: 0.9rem;
}

.page-hero-pill {
    display: inline-flex;
    align-items: center;
    padding: 0.36rem 0.72rem;
    border-radius: 999px;
    border: 1px solid rgba(31,122,255,.12);
    background: rgba(255,255,255,.78);
    color: #32557e;
    font-size: 0.84rem;
    font-weight: 700;
    box-shadow: inset 0 1px 0 rgba(255,255,255,.78);
}

.status-strip {
    display: grid;
    grid-template-columns: repeat(auto-fit, minmax(180px, 1fr));
    gap: 0.8rem;
    margin: -0.1rem 0 1.15rem 0;
}

.status-card {
    border-radius: 12px;
    border: 1px solid rgba(215, 226, 241, 0.96);
    background: rgba(255,255,255,.82);
    box-shadow: 0 12px 26px rgba(24,57,95,.06);
    backdrop-filter: blur(10px);
    padding: 0.85rem 0.95rem 0.8rem 0.95rem;
}

.status-card-label {
    color: #6a7d97;
    font-size: 0.8rem;
    font-weight: 700;
    margin-bottom: 0.28rem;
}

.status-card-value {
    color: #173a60;
    font-size: 1.3rem;
    font-weight: 900;
    line-height: 1.15;
}

.status-card-note {
    color: #6b7a90;
    font-size: 0.82rem;
    line-height: 1.45;
    margin-top: 0.28rem;
}

.metric-card {
    background: var(--cf-teal);
    border-radius: 8px;
    padding: 1.4rem;
    color: white;
    text-align: center;
    box-shadow: var(--cf-shadow);
}

.info-box,
.success-box,
.warning-box {
    border-radius: 8px;
    padding: 1rem 1.1rem;
    margin: 0.75rem 0;
    box-shadow: 0 2px 8px rgba(22,43,38,0.03);
}

.info-box {
    background-color: rgba(31,122,255,.08);
    border-left: 5px solid var(--cf-blue);
}

.success-box {
    background-color: #e7f6ed;
    border-left: 5px solid #238a50;
}

.warning-box {
    background-color: #fff8ea;
    border-left: 5px solid var(--cf-amber);
}

.action-zone {
    margin: 1.6rem 0 1rem 0;
    padding: 1rem 1.1rem 1rem 1.25rem;
    border: 1px solid var(--cf-border);
    border-left: 8px solid var(--cf-amber);
    border-radius: 8px;
    background: linear-gradient(135deg, #fff7ea, #fffdf7);
    color: var(--cf-text);
    font-weight: 800;
    position: relative;
    box-shadow: 0 12px 24px rgba(24,57,95,0.06);
    outline: 1px solid rgba(183,121,31,0.12);
}

.parameter-zone,
.result-zone,
.export-zone,
.danger-zone {
    margin: 1.6rem 0 1rem 0;
    padding: 1rem 1.1rem 1rem 1.25rem;
    border-radius: 8px;
    font-weight: 800;
    border: 1px solid var(--cf-border);
    position: relative;
    box-shadow: 0 12px 24px rgba(24,57,95,0.06);
    outline: 1px solid rgba(24,57,95,0.04);
}

.action-zone::before,
.parameter-zone::before,
.result-zone::before,
.export-zone::before,
.danger-zone::before,
.form-section::before {
    display: inline-flex;
    align-items: center;
    justify-content: center;
    width: 1.35rem;
    height: 1.35rem;
    margin-right: 0.45rem;
    border-radius: 999px;
    color: #ffffff;
    font-size: 0.82rem;
    font-weight: 850;
    vertical-align: -0.08rem;
}

.action-zone::before {
    content: "▶";
    background: var(--cf-amber);
}

.parameter-zone::before {
    content: "⚙";
    background: var(--cf-teal);
}

.result-zone::before {
    content: "✓";
    background: var(--cf-blue);
}

.export-zone::before {
    content: "↓";
    background: var(--cf-amber);
}

.danger-zone::before {
    content: "!";
    background: var(--cf-red);
}

.parameter-zone {
    border-left: 8px solid var(--cf-teal);
    background: linear-gradient(135deg, rgba(31,122,255,.08), rgba(255,255,255,.96));
}

.result-zone {
    border-left: 8px solid var(--cf-blue);
    background: linear-gradient(135deg, rgba(35,120,247,.12), rgba(255,255,255,.97));
}

.export-zone {
    border-left: 8px solid var(--cf-amber);
    background: linear-gradient(135deg, #fff1cf 0%, #fffaf0 100%);
}

.danger-zone {
    border-left: 8px solid var(--cf-red);
    background: linear-gradient(90deg, #ffe4e2 0%, #fff3f2 100%);
    color: #7a1b14;
    box-shadow: 0 6px 18px rgba(180,35,24,0.08);
}

div[data-testid="stForm"],
div[data-testid="stExpander"],
div[data-testid="stDataFrame"],
div[data-testid="stTable"],
div[data-testid="stMetric"] {
    border-radius: 8px;
}

div[data-testid="stForm"] {
    background: rgba(255,255,255,.84);
    border: 1px solid var(--cf-border);
    box-shadow: 0 12px 30px rgba(24,57,95,.06);
    padding: 1.15rem 1.2rem;
    backdrop-filter: blur(10px);
}

div[data-testid="stForm"] > div {
    gap: 0.85rem;
}

div[data-testid="stExpander"] {
    border: 1px solid var(--cf-border);
    background: rgba(255,255,255,.84);
    box-shadow: 0 8px 22px rgba(24,57,95,.05);
    backdrop-filter: blur(10px);
}

div[data-testid="stDataFrame"],
div[data-testid="stTable"] {
    border: 1px solid var(--cf-border);
    box-shadow: 0 8px 22px rgba(24,57,95,0.04);
    background: rgba(255,255,255,.88);
}

div[data-testid="stPlotlyChart"],
div[data-testid="stImage"],
div[data-testid="stJson"] {
    border-radius: 12px;
    border: 1px solid rgba(220,230,242,.96);
    background: rgba(255,255,255,.86);
    box-shadow: 0 12px 28px rgba(24,57,95,.06);
    padding: 0.8rem 0.9rem;
    backdrop-filter: blur(10px);
}

div[data-testid="stMetric"] {
    background: rgba(255,255,255,.92);
    border: 1px solid var(--cf-border);
    padding: 0.95rem 1rem;
    box-shadow: 0 10px 24px rgba(24,57,95,0.06);
    position: relative;
    overflow: hidden;
    border-radius: 12px;
}

div[data-testid="stMetric"]::before {
    content: "";
    position: absolute;
    inset: 0 0 auto 0;
    height: 4px;
    background: linear-gradient(90deg, #2a86ff, #1d6ff2);
}

div[data-testid="stMetric"] label,
div[data-testid="stMetric"] [data-testid="stMetricLabel"] {
    color: var(--cf-muted) !important;
    font-weight: 700;
}

div[data-testid="stMetric"] [data-testid="stMetricValue"] {
    color: var(--cf-text) !important;
    font-weight: 800;
    letter-spacing: 0;
}

div[data-testid="stMetric"] [data-testid="stMetricDelta"] {
    font-weight: 700;
}

div[data-testid="stMetric"]:hover {
    transform: translateY(-1px);
    box-shadow: 0 7px 18px rgba(22,43,38,0.07);
    transition: all 140ms ease;
}

.stTabs [data-baseweb="tab-list"] {
    gap: 0.6rem;
    border-bottom: 1px solid var(--cf-border);
}

.stTabs [data-baseweb="tab"] {
    background: rgba(255,255,255,.68);
    border: 1px solid var(--cf-border);
    border-bottom: none;
    border-radius: 10px 10px 0 0;
    color: #273a56;
    padding: 0.62rem 1rem;
    font-weight: 800;
}

.stTabs [aria-selected="true"] {
    background: #ffffff;
    color: var(--cf-blue);
    border-top: 3px solid var(--cf-blue);
    font-weight: 700;
    box-shadow: 0 6px 16px rgba(24,57,95,.05);
}

div.stButton > button,
div[data-testid="stFormSubmitButton"] > button {
    border-radius: 8px;
    border: 1px solid #d9e3f0;
    background: #ffffff;
    color: #344966;
    font-weight: 700;
    min-height: 2.45rem;
    box-shadow: 0 8px 18px rgba(24,57,95,.08);
    transition: all 150ms ease;
}

div.stButton > button:hover,
div[data-testid="stFormSubmitButton"] > button:hover {
    border-color: rgba(31,122,255,.36);
    background: #ffffff;
    color: var(--cf-blue);
    transform: translateY(-1px);
    box-shadow: 0 12px 24px rgba(24,57,95,.1);
}

div.stButton > button:active,
div[data-testid="stFormSubmitButton"] > button:active {
    transform: translateY(0);
    box-shadow: 0 2px 6px rgba(15,118,110,0.12);
}

button[kind="primary"],
div[data-testid="stFormSubmitButton"] button[kind="primary"] {
    background: linear-gradient(135deg, #2484ff, #116de9) !important;
    color: white !important;
    border-color: var(--cf-blue) !important;
    box-shadow: 0 10px 20px rgba(31,122,255,.2) !important;
}

button[kind="primary"]:hover,
div[data-testid="stFormSubmitButton"] button[kind="primary"]:hover {
    background: linear-gradient(135deg, #1f7aff, #0f5dd1) !important;
    color: white !important;
}

div[data-testid="stDownloadButton"] > button {
    border-radius: 8px;
    border: 1px solid #d9e3f0;
    background: rgba(255,255,255,.92);
    color: #28405f;
    font-weight: 700;
    min-height: 2.45rem;
    box-shadow: 0 8px 18px rgba(24,57,95,.08);
}

div[data-testid="stDownloadButton"] > button:hover {
    background: #ffffff;
    color: var(--cf-blue);
    border-color: rgba(31,122,255,.36);
    transform: translateY(-1px);
}

div.st-key-delete_selected_sample button,
div.st-key-confirm_delete_selected_sample button {
    background: #fff0ef !important;
    border: 1.5px solid var(--cf-red) !important;
    color: #8f1d16 !important;
    box-shadow: 0 2px 8px rgba(180,35,24,0.12) !important;
}

div.st-key-delete_selected_sample button:hover,
div.st-key-confirm_delete_selected_sample button:hover {
    background: var(--cf-red) !important;
    border-color: #8f1d16 !important;
    color: white !important;
    box-shadow: 0 5px 16px rgba(180,35,24,0.24) !important;
}

input,
textarea,
div[data-baseweb="select"] > div {
    border-radius: 8px !important;
}

div[data-testid="stTextInput"],
div[data-testid="stNumberInput"],
div[data-testid="stSelectbox"],
div[data-testid="stRadio"] {
    background: rgba(255,255,255,.78);
    border: 1px solid #d5e1ef;
    border-radius: 8px;
    padding: 0.65rem 0.75rem;
    margin-bottom: 0.35rem;
    box-shadow: inset 0 1px 0 rgba(255,255,255,.85);
}

div[data-testid="stTextInput"] label,
div[data-testid="stNumberInput"] label,
div[data-testid="stSelectbox"] label,
div[data-testid="stRadio"] label {
    color: var(--cf-text) !important;
    font-weight: 700;
}

div[data-testid="stTextInput"]:focus-within,
div[data-testid="stNumberInput"]:focus-within,
div[data-testid="stSelectbox"]:focus-within,
div[data-testid="stRadio"]:focus-within {
    border-color: rgba(31,122,255,.6);
    background: #ffffff;
    box-shadow: 0 0 0 4px rgba(31,122,255,.08);
}

.form-section {
    margin: 1.05rem 0 0.75rem 0;
    padding: 0.72rem 0.9rem 0.72rem 1rem;
    border-left: 5px solid var(--cf-blue);
    border-radius: 8px;
    background: rgba(255,255,255,.75);
    color: #17375d;
    font-weight: 750;
    box-shadow: inset 0 0 0 1px rgba(31,122,255,0.06);
    border-top: 1px solid rgba(31,122,255,0.16);
    border-bottom: 1px solid rgba(31,122,255,0.10);
    position: relative;
}

.form-section::before {
    content: "•";
    width: 1.1rem;
    height: 1.1rem;
    background: var(--cf-blue);
    font-size: 0.7rem;
}

.result-hint {
    margin: -0.2rem 0 1rem 0;
    padding: 0.75rem 0.9rem;
    border: 1px solid #c9dcff;
    border-left: 5px solid var(--cf-blue);
    border-radius: 8px;
    background: rgba(255,255,255,.78);
    color: #244264;
    font-size: 0.95rem;
    line-height: 1.55;
    box-shadow: 0 8px 18px rgba(24,57,95,.05);
}

div[data-testid="stAlert"] {
    border-radius: 8px;
    border: 1px solid rgba(24,57,95,0.08);
}
</style>
"""

def inject_global_style() -> None:
    css = GLOBAL_STYLE
    try:
        if CSS_PATH.exists():
            css = f"<style>\n{CSS_PATH.read_text(encoding='utf-8')}\n</style>"
    except Exception:
        css = GLOBAL_STYLE
    st.markdown(css, unsafe_allow_html=True)


def render_page_hero(title: str, subtitle: str, eyebrow: str = "Workflow", pills: list[str] | None = None) -> None:
    pill_html = ""
    if pills:
        pill_html = '<div class="page-hero-pills">' + "".join(
            f'<span class="page-hero-pill">{escape(str(p))}</span>' for p in pills if str(p).strip()
        ) + "</div>"
    st.markdown(
        f"""
        <section class="page-hero">
            <div class="page-hero-layout">
                <div class="page-hero-content">
                    <div class="page-hero-eyebrow">{escape(str(eyebrow))}</div>
                    <h1 class="page-hero-title">{escape(str(title))}</h1>
                    {pill_html}
                </div>
                <div class="page-hero-visual" aria-hidden="true">
                    <div class="page-hero-dots"></div>
                    <div class="page-hero-wave"></div>
                    <div class="page-hero-macaque">macaque atlas</div>
                    <div class="page-hero-brain">🧠</div>
                    <div class="page-hero-dna">⟲</div>
                    <div class="page-hero-tube">🧪</div>
                </div>
            </div>
        </section>
        """,
        unsafe_allow_html=True,
    )


def render_status_strip(cards: list[dict]) -> None:
    blocks = []
    for card in cards:
        label = str(card.get("label", "")).strip()
        value = str(card.get("value", "")).strip()
        if not (label or value):
            continue
        blocks.append(
            f'<div class="status-card">'
            f'<div class="status-card-label">{escape(label)}</div>'
            f'<div class="status-card-value">{escape(value)}</div>'
            f"</div>"
        )
    if blocks:
        st.markdown(f'<section class="status-strip">{"".join(blocks)}</section>', unsafe_allow_html=True)


def render_result_hint(text: str) -> None:
    return None


@st.cache_resource
def init_database() -> str:
    db_path = Path(DB_PATH)
    db = CSFRNASourceDatabase(DB_PATH)
    if not db_path.exists(): db.initialize_database(); db.close()
    else: db.connect(); db.create_database_schema(); db.close()
    run_migrations(DB_PATH)
    buildkit_dir = PROJECT_ROOT / "bo2023_bulk_atlas_buildkit"
    try:
        if buildkit_dir.exists():
            conn = sqlite3.connect(DB_PATH)
            try:
                cur = conn.cursor()
                cur.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='bo2023_buildkit_catalog'")
                has_catalog = cur.fetchone() is not None
            finally:
                conn.close()
            if not has_catalog:
                import_buildkit_dir(DB_PATH, buildkit_dir)
    except Exception:
        pass
    return DB_PATH


@st.cache_data(show_spinner=False, ttl=60)
def get_startup_check() -> Dict:
    return run_system_check(DB_PATH, PROJECT_ROOT)


def render_startup_check_summary(expanded: bool = False) -> None:
    check = get_startup_check()
    items = check.get("items", [])
    n_errors = sum(1 for i in items if getattr(i, "status", "") == "error")
    n_warnings = sum(1 for i in items if getattr(i, "status", "") == "warning")
    if n_errors:
        st.error(f"Startup self-check found {n_errors} blocking issue(s).")
    elif n_warnings:
        st.warning(f"Startup self-check passed with {n_warnings} warning(s).")
    else:
        st.success("Startup self-check passed.")

    with st.expander("Startup self-check details", expanded=expanded or bool(n_errors)):
        for item in items:
            status = getattr(item, "status", "warning")
            text = f"**{getattr(item, 'name', 'Check')}**: {getattr(item, 'detail', '')}"
            if status == "ok":
                st.success(text)
            elif status == "error":
                st.error(text)
            else:
                st.warning(text)


@st.cache_resource
def init_tracer() -> CSFRNASourceTracer:
    init_database(); tracer = CSFRNASourceTracer(DB_PATH); tracer.load_reference_data(); return tracer


def _compute_database_cohort_qc_fallback(processor: DataProcessor):
    samples_df = processor.get_all_samples()
    if samples_df.empty:
        return samples_df

    sample_map = {}
    for sample_id in samples_df["sample_id"].astype(str).tolist():
        expr_df = processor.get_sample_expression(sample_id)
        if expr_df is not None and not expr_df.empty:
            sample_map[sample_id] = expr_df

    if not sample_map:
        return samples_df.iloc[0:0].copy()

    from data.qc import compute_cohort_qc

    cohort_qc = compute_cohort_qc(sample_map)
    rows = []
    for sample_id, qc in cohort_qc.items():
        rows.append(
            {
                "sample_id": sample_id,
                "overall_risk": qc.get("overall_risk"),
                "gene_id_type": qc.get("gene_id_type"),
                "rbc_score": qc.get("rbc_mrna_score"),
                "rbc_percentile": qc.get("rbc_mrna_percentile"),
                "rbc_risk": qc.get("rbc_mrna_risk"),
                "immune_score": qc.get("immune_mrna_score"),
                "immune_percentile": qc.get("immune_mrna_percentile"),
                "immune_risk": qc.get("immune_mrna_risk"),
                "brain_score": qc.get("brain_marker_score"),
                "brain_percentile": qc.get("brain_marker_percentile"),
                "brain_risk": qc.get("brain_marker_risk"),
                "hemolysis_mirna_risk": qc.get("hemolysis_mirna_risk"),
                "interpretation": qc.get("interpretation"),
            }
        )

    qc_df = pd.DataFrame(rows)
    if qc_df.empty:
        return qc_df
    return samples_df.merge(qc_df, on="sample_id", how="left")


@st.cache_resource
def init_processor() -> DataProcessor:
    init_database()
    processor = DataProcessor(DB_PATH)
    if not hasattr(processor, "compute_database_cohort_qc"):
        processor.compute_database_cohort_qc = types.MethodType(_compute_database_cohort_qc_fallback, processor)
    return processor

def get_database_stats(conn: sqlite3.Connection) -> Dict:
    stats = {}; cursor = conn.cursor();
    for table in ['macaque_brain_atlas', 'reference_expression', 'cfrna_samples', 'cfrna_expression', 'source_tracing_results', 'region_gene_signature']:
        cursor.execute(f"SELECT COUNT(*) FROM {table}"); stats[table] = cursor.fetchone()[0]
    return stats
