from __future__ import annotations

from pathlib import Path

import pandas as pd
import streamlit as st

from src.reference_tracing.scoring import run_reference_tracing


st.set_page_config(page_title="Reference Tracing v1", layout="wide")
st.title("cfRNA Brain Tracing Reference v1")
st.caption("Research exploration output: evidence, signal, risk, and candidate source. Not a clinical decision model.")

uploaded = st.file_uploader("Upload expression matrix (TSV or CSV)", type=["tsv", "txt", "csv"])
expr_type = st.selectbox("Expression type", ["raw counts", "TPM", "CPM", "FPKM"])

if uploaded is not None:
    sep = "," if uploaded.name.lower().endswith(".csv") else "\t"
    expr = pd.read_csv(uploaded, sep=sep)
    gene_col = st.selectbox("Gene column", list(expr.columns))
    st.dataframe(expr.head(20), use_container_width=True)

    if st.button("Run cfRNA brain tracing v1", type="primary"):
        outdir = Path("results") / "cfrna_tracing_v1" / "streamlit_last_run"
        with st.spinner("Scoring samples against reference v1..."):
            result = run_reference_tracing(expr, gene_col, expr_type, outdir, make_plots=True)

        overall = result["sample_overall_tracing_summary"]
        contam = result["sample_contamination_scores"]
        brain = result["sample_brain_signal_scores"]
        cell = result["sample_celltype_scores"]
        region = result["sample_region_top_hits"]
        injury = result["sample_injury_state_scores"]

        st.subheader("Overall summary")
        st.dataframe(overall, use_container_width=True, hide_index=True)

        if (overall["rbc_risk"].eq("High") | overall["immune_risk"].eq("High")).any():
            st.warning("High peripheral contamination risk detected; blood/immune background may confound brain signal interpretation.")

        c1, c2 = st.columns(2)
        with c1:
            st.subheader("Brain signal")
            st.dataframe(brain, use_container_width=True, hide_index=True)
            st.subheader("Top celltype evidence")
            st.dataframe(cell, use_container_width=True, hide_index=True)
        with c2:
            st.subheader("Contamination")
            st.dataframe(contam, use_container_width=True, hide_index=True)
            st.subheader("Injury-state score")
            st.dataframe(injury, use_container_width=True, hide_index=True)

        st.subheader("Top region evidence")
        st.info("Region score reflects similarity to reference brain-region expression signatures and should not be interpreted as definitive anatomical localization.")
        st.dataframe(region, use_container_width=True, hide_index=True)

        downloads = {
            "sample_overall_tracing_summary.tsv": overall,
            "sample_celltype_scores.tsv": cell,
            "sample_region_top_hits.tsv": region,
            "sample_injury_state_scores.tsv": injury,
        }
        st.subheader("Downloads")
        for fname, df in downloads.items():
            st.download_button(fname, df.to_csv(sep="\t", index=False).encode("utf-8"), file_name=fname)
        report_path = outdir / "cfrna_tracing_report.md"
        if report_path.exists():
            st.download_button(
                "cfrna_tracing_report.md",
                report_path.read_text(encoding="utf-8").encode("utf-8"),
                file_name="cfrna_tracing_report.md",
            )
else:
    st.info("Upload a gene x sample expression matrix to run reference tracing v1.")
