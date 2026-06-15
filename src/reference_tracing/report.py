from __future__ import annotations

from pathlib import Path

import pandas as pd


def _top_lines(df: pd.DataFrame, cols: list[str], n: int = 5) -> str:
    if df.empty:
        return "_No rows._"
    return df[cols].head(n).to_markdown(index=False) if hasattr(df, "to_markdown") else df[cols].head(n).to_csv(sep="\t", index=False)


def _table_md(df: pd.DataFrame, cols: list[str], n: int = 5) -> str:
    view = df[cols].head(n).fillna("").astype(str)
    lines = ["| " + " | ".join(cols) + " |", "| " + " | ".join(["---"] * len(cols)) + " |"]
    for _, row in view.iterrows():
        lines.append("| " + " | ".join(str(row[c]).replace("|", "\\|") for c in cols) + " |")
    return "\n".join(lines)


def write_report(path: str | Path, input_meta: dict, outputs: dict) -> None:
    overall = outputs["sample_overall_tracing_summary"]
    contam = outputs["sample_contamination_scores"]
    brain = outputs["sample_brain_signal_scores"]
    cell = outputs["sample_celltype_scores"]
    region = outputs["sample_region_top_hits"]
    injury = outputs["sample_injury_state_scores"]
    lines = [
        "# cfRNA Brain Injury Source-Tracing Report",
        "",
        "## 1. Input summary",
        "",
        f"- sample number: {input_meta.get('sample_number')}",
        f"- gene number: {input_meta.get('gene_number')}",
        f"- detected genes: {input_meta.get('detected_genes')}",
        f"- matched reference genes: {input_meta.get('matched_reference_genes')}",
        f"- transform: {input_meta.get('expression_transform')}",
        f"- small cohort warning: {input_meta.get('small_cohort_warning')}",
        "",
        "## 2. Contamination assessment",
        "",
        "RBC and immune scores are marker mean log1p expression summaries. High contamination risk means peripheral contamination may confound brain signal interpretation.",
        "",
        _table_md(contam, ["sample", "rbc_score", "rbc_risk", "immune_score", "immune_risk"]),
        "",
        "## 3. Brain-enriched signal",
        "",
        "Brain signal scores are relative evidence scores based on GTEx brain-enriched markers; they are not absolute brain-source proportions.",
        "",
        _table_md(brain, ["sample", "brain_signal_score", "brain_signal_level", "number_of_detected_brain_markers"]),
        "",
        "## 4. Cell-type source evidence",
        "",
        _table_md(cell, ["sample", "top_celltype_1", "top_celltype_2", "top_celltype_3"]),
        "",
        "## 5. Brain-region similarity",
        "",
        "Region score reflects similarity to reference brain-region expression signatures and should not be interpreted as definitive anatomical localization.",
        "",
        _table_md(region, [c for c in ["sample", "top_region_1", "top_region_2", "top_region_3"] if c in region.columns]),
        "",
        "## 6. Injury-state evidence",
        "",
        "TBI/injury response score is candidate injury-state evidence. Garza TBI grouping is inferred from sample names and has metadata limitations.",
        "",
        _table_md(injury, ["sample", "injury_state_score", "injury_state_level"]),
        "",
        "## 7. Overall interpretation",
        "",
        _table_md(overall, ["sample", "brain_signal_level", "rbc_risk", "immune_risk", "injury_state_level", "overall_interpretation", "warning_flags"], n=20),
        "",
        "## 8. Limitations",
        "",
        "- cfRNA expression is affected by RNA stability, release mechanisms, and peripheral background.",
        "- The v1 reference is a research exploration tool, not a clinical decision model.",
        "- Garza TBI metadata has limitations; TBI/control grouping is inferred from sample names.",
        "- HPA aggregate lacks detection fraction, so celltype markers are candidate markers rather than a final gold standard.",
        "- Brain-region source resolution should be interpreted conservatively.",
        "- Interpretation should be combined with longitudinal personal baselines plus experimental and clinical context.",
    ]
    Path(path).write_text("\n".join(lines), encoding="utf-8")
