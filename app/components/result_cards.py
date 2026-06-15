from __future__ import annotations

import streamlit as st

from app.i18n import tr


def render_run_meta(meta: dict) -> None:
    sigset_id = meta.get("sigset_id")
    st.caption(
        f"{tr('运行参数', 'Run parameters')}: "
        f"atlas_id={meta.get('atlas_id')} | "
        f"sigset_id={sigset_id} | "
        f"use_value={meta.get('use_value')} | "
        f"bootstrap_n={meta.get('bootstrap_n')}"
    )
    if meta.get("atlas_normalization"):
        mode_label = (
            tr("VSD 兼容表达模式相似性", "VSD-compatible pattern similarity")
            if meta.get("vsd_compatible_mode")
            else tr("TPM 兼容溯源", "TPM-compatible tracing")
        )
        st.caption(
            f"{tr('参考图谱标准化', 'Reference normalization')}: "
            f"{meta.get('atlas_normalization')} | "
            f"{tr('解释模式', 'Interpretation mode')}: "
            f"{mode_label}"
        )
    if sigset_id is None:
        st.warning(
            tr(
                "当前分析未启用 signature 集合，结果可能更依赖全基因或 marker 模式，区分度会相对保守。",
                "No signature set was used in this run. The result may rely more on full-gene or marker-only evidence and can be more conservative.",
            )
        )


def render_primary_metrics(results_rows, meta) -> None:
    c1, c2, c3 = st.columns(3)
    if results_rows:
        top_row = results_rows[0]
        c1.metric(
            tr("最可能来源脑区", "Top source region"),
            str(top_row.get("region_id", "NA")),
            delta=f"score: {float(top_row.get('score', 0.0)):.3f}",
        )
        c2.metric(
            tr("Top1 置信度", "Top1 confidence"),
            f"{float(top_row.get('confidence', 0.0)):.3f}",
        )
    else:
        c1.metric(tr("最可能来源脑区", "Top source region"), tr("无法确定", "Undetermined"))
        c2.metric(tr("Top1 置信度", "Top1 confidence"), "0.000")
    c3.metric(tr("使用基因数", "Genes used"), f"{int(meta.get('n_genes', 0))}")
