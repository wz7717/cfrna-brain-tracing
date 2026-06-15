from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def _save(fig, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def _bar(df: pd.DataFrame, x: str, y: str, title: str, path: Path, color: str = "#36688D") -> None:
    fig, ax = plt.subplots(figsize=(max(6, len(df) * 0.7), 4))
    ax.bar(df[x].astype(str), pd.to_numeric(df[y], errors="coerce").fillna(0), color=color)
    ax.set_title(title)
    ax.set_ylabel(y.replace("_", " "))
    ax.tick_params(axis="x", rotation=45, labelsize=8)
    _save(fig, path)


def _heatmap(matrix: pd.DataFrame, title: str, path: Path) -> None:
    if matrix.empty:
        fig, ax = plt.subplots(figsize=(6, 2))
        ax.text(0.5, 0.5, "No data", ha="center", va="center")
        ax.axis("off")
        _save(fig, path)
        return
    values = matrix.apply(pd.to_numeric, errors="coerce").fillna(0).to_numpy()
    fig, ax = plt.subplots(figsize=(max(6, matrix.shape[1] * 0.55), max(3, matrix.shape[0] * 0.35)))
    im = ax.imshow(values, aspect="auto", cmap="viridis")
    ax.set_title(title)
    ax.set_yticks(range(matrix.shape[0]))
    ax.set_yticklabels(matrix.index.astype(str), fontsize=8)
    ax.set_xticks(range(matrix.shape[1]))
    ax.set_xticklabels(matrix.columns.astype(str), rotation=45, ha="right", fontsize=8)
    fig.colorbar(im, ax=ax, fraction=0.025)
    _save(fig, path)


def make_all_plots(outputs: dict, figdir: str | Path) -> None:
    figdir = Path(figdir)
    brain = outputs["sample_brain_signal_scores"]
    contam = outputs["sample_contamination_scores"]
    cell = outputs["sample_celltype_scores"]
    region = outputs["sample_region_scores"]
    injury = outputs["sample_injury_state_scores"]
    overall = outputs["sample_overall_tracing_summary"]

    _bar(brain, "sample", "brain_signal_score", "Brain-enriched cfRNA signal", figdir / "sample_brain_signal_barplot.png", "#2F6F4E")

    cmat = contam.set_index("sample")[["rbc_score", "immune_score"]]
    _heatmap(cmat, "RBC and immune contamination scores", figdir / "sample_contamination_scores.png")

    ccols = [c for c in cell.columns if c not in {"sample", "top_celltype_1", "top_celltype_2", "top_celltype_3", "top_celltype_scores"}]
    _heatmap(cell.set_index("sample")[ccols], "Brain cell-type evidence scores", figdir / "sample_celltype_heatmap.png")

    if not region.empty:
        top_regions = (
            region.sort_values("combined_region_score", ascending=False)
            .groupby("sample")
            .head(10)
            .pivot_table(index="sample", columns="region", values="combined_region_score", aggfunc="max")
            .fillna(0)
        )
    else:
        top_regions = pd.DataFrame()
    _heatmap(top_regions, "Top brain-region similarity scores", figdir / "sample_region_top_hits_heatmap.png")

    _bar(injury, "sample", "injury_state_score", "TBI/injury-state candidate evidence", figdir / "sample_injury_state_scores.png", "#9A4C36")

    risk_map = {"Low": 0, "Moderate": 1, "High": 2}
    risk = overall.set_index("sample")[["brain_signal_level", "rbc_risk", "immune_risk", "injury_state_level"]].apply(lambda col: col.map(risk_map).fillna(0))
    _heatmap(risk, "Overall evidence/risk levels", figdir / "sample_overall_risk_heatmap.png")
