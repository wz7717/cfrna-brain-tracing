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
DB_PATH = str(Path(os.environ.get("BRAINTRACE_DB_PATH", PROJECT_ROOT / "braintrace_source_tracing.db")).resolve())
CSS_PATH = PROJECT_ROOT / "app" / "styles" / "main.css"
GLOBAL_STYLE = """
<style>
:root {
    --app-bg: #f8fafc;
    --app-text: #0f172a;
    --app-muted: #64748b;
    --app-border: #dbe3ef;
    --app-primary: #1f7aff;
}

html, body, [data-testid="stAppViewContainer"] {
    background: var(--app-bg);
    color: var(--app-text);
}

.top-toolbar,
.kpi-grid,
.status-strip,
.mini-grid {
    display: flex;
    gap: 0.75rem;
    flex-wrap: wrap;
}

.top-toolbar,
.kpi-card,
.status-card,
.mini-card,
.section-band {
    border: 1px solid var(--app-border);
    border-radius: 8px;
    background: #ffffff;
}

.top-toolbar {
    align-items: center;
    justify-content: space-between;
    padding: 0.7rem 1rem;
    margin-bottom: 1rem;
}

.top-toolbar-meta,
.top-toolbar-links {
    display: flex;
    align-items: center;
    gap: 0.5rem;
    flex-wrap: wrap;
}

.top-toolbar-chip {
    color: var(--app-muted);
    font-size: 0.85rem;
}

.kpi-card,
.status-card,
.mini-card {
    padding: 0.8rem;
}

.section-band {
    padding: 0.75rem 1rem;
    margin: 1rem 0 0.75rem;
}

.panel-title {
    font-weight: 700;
    margin: 0.5rem 0;
}

.action-zone,
.parameter-zone,
.result-zone,
.export-zone,
.form-section {
    color: var(--app-muted);
    font-weight: 700;
    margin: 1rem 0 0.5rem;
}
</style>
"""


def is_public_demo_mode() -> bool:
    return os.environ.get("BRAINTRACE_PUBLIC_DEMO", "").strip().lower() in {"1", "true", "yes", "on"}


def inject_global_style() -> None:
    try:
        if CSS_PATH.exists():
            css = f"<style>\n{CSS_PATH.read_text(encoding='utf-8')}\n</style>"
        else:
            css = GLOBAL_STYLE
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
    for table in ['macaque_brain_atlas', 'reference_expression', 'braintrace_samples', 'braintrace_expression', 'source_tracing_results', 'region_gene_signature']:
        cursor.execute(f"SELECT COUNT(*) FROM {table}"); stats[table] = cursor.fetchone()[0]
    return stats
