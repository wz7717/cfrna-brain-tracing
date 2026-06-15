from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pandas as pd

from .db import DEFAULT_DB, read_marker_candidates, read_table
from .io import ensure_outdir
from .preprocessing import level_from_percentile, percentile_series, prepare_expression_matrix, zscore_rows


RBC_GENES = ["HBA1", "HBA2", "HBB", "ALAS2", "CA1", "CA2", "SLC4A1", "GYPA"]
IMMUNE_GENES = ["PTPRC", "LST1", "TYROBP", "CD74", "HLA-DRA", "S100A8", "S100A9"]
CORE_GROUPS = {
    "neuron": ["RBFOX3", "SNAP25", "SYT1", "TUBB3", "MAP2"],
    "astrocyte": ["GFAP", "AQP4", "ALDH1L1", "SLC1A3"],
    "oligodendrocyte": ["MBP", "PLP1", "MOBP", "MAG", "MOG"],
    "OPC": ["PDGFRA", "VCAN", "CSPG4"],
    "microglia": ["P2RY12", "CX3CR1", "TMEM119", "AIF1"],
    "endothelial": ["CLDN5", "PECAM1", "VWF", "KDR"],
}
META = {"gene_symbol", "gene_symbol_raw", "ensembl_id", "source_atlas", "doi", "recommended_use", "reference_subtype"}


def _present(expr: pd.DataFrame, genes: list[str]) -> list[str]:
    return [g for g in genes if g in expr.index]


def _mean_score(expr: pd.DataFrame, genes: list[str]) -> pd.Series:
    genes = _present(expr, genes)
    if not genes:
        return pd.Series(0.0, index=expr.columns)
    return expr.loc[genes].mean(axis=0)


def _risk_table(expr: pd.DataFrame, genes: list[str], label: str) -> pd.DataFrame:
    score = _mean_score(expr, genes)
    percentile = percentile_series(score)
    z = (score - score.mean()) / (score.std() if score.std() else np.nan)
    z = z.fillna(0)
    return pd.DataFrame(
        {
            "sample": score.index,
            f"{label}_score": score.values,
            f"{label}_percentile": percentile.values,
            f"{label}_zscore": z.values,
            f"{label}_risk": [level_from_percentile(p, s) for p, s in zip(percentile, score)],
            f"{label}_detected_markers": len(_present(expr, genes)),
        }
    )


def score_contamination(expr: pd.DataFrame) -> pd.DataFrame:
    rbc = _risk_table(expr, RBC_GENES, "rbc")
    immune = _risk_table(expr, IMMUNE_GENES, "immune")
    return rbc.merge(immune, on="sample", how="outer")


def score_brain_signal(expr: pd.DataFrame, markers: pd.DataFrame) -> pd.DataFrame:
    brain = markers[markers["marker_type_normalized"].eq("brain_enriched")]
    genes = sorted(set(brain["gene_symbol"].dropna().astype(str).str.upper()))
    genes = _present(expr, genes)
    if not genes:
        score = pd.Series(0.0, index=expr.columns)
        frac = 0.0
        detected = 0
    else:
        z = zscore_rows(expr.loc[genes])
        median_z = z.median(axis=0)
        detected_by_sample = (expr.loc[genes] > 0).sum(axis=0)
        frac_by_sample = detected_by_sample / max(len(genes), 1)
        score = median_z + 0.5 * frac_by_sample
        frac = None
        detected = None
    percentile = percentile_series(score)
    background = expr.drop(index=genes, errors="ignore").mean(axis=0) if len(expr.index) > len(genes) else pd.Series(0.0, index=expr.columns)
    return pd.DataFrame(
        {
            "sample": expr.columns,
            "brain_marker_mean": _mean_score(expr, genes).reindex(expr.columns).values,
            "background_gene_mean": background.reindex(expr.columns).values,
            "brain_signal_score": score.reindex(expr.columns).values,
            "brain_marker_zscore": score.reindex(expr.columns).values,
            "brain_marker_percentile": percentile.reindex(expr.columns).values,
            "number_of_detected_brain_markers": (expr.loc[genes] > 0).sum(axis=0).reindex(expr.columns).values if genes else 0,
            "detected_marker_fraction": ((expr.loc[genes] > 0).sum(axis=0) / max(len(genes), 1)).reindex(expr.columns).values if genes else 0,
            "brain_signal_level": [level_from_percentile(p, s) for p, s in zip(percentile.reindex(expr.columns), score.reindex(expr.columns))],
        }
    )


def normalize_celltype(name: str) -> str | None:
    text = str(name).upper()
    if any(k in text for k in ["NEURON", "EXCITATORY", "INHIBITORY", "GLUTAMATERGIC", "GABAERGIC", "INTERNEURON"]):
        return "neuron"
    if "ASTRO" in text:
        return "astrocyte"
    if "OPC" in text or "OLIGODENDROCYTE_PRECURSOR" in text or "OLIGODENDROCYTE PRECURSOR" in text:
        return "OPC"
    if "OLIGO" in text or "OLIGODENDROCYTE" in text:
        return "oligodendrocyte"
    if "MICROGLIA" in text or "MACROPHAGE" in text:
        return "microglia"
    if any(k in text for k in ["ENDOTHELIAL", "VASCULAR", "PERICYTE"]):
        return "endothelial"
    return None


def score_celltypes(expr: pd.DataFrame, markers: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    cell = markers[markers["marker_type_normalized"].eq("brain_celltype_marker")].copy()
    rows = []
    for raw_celltype, sub in cell.groupby("celltype", dropna=True):
        broad = normalize_celltype(raw_celltype)
        if not broad:
            continue
        genes = _present(expr, sorted(set(sub["gene_symbol"].astype(str).str.upper())))
        if len(genes) < 2:
            continue
        score = expr.loc[genes].mean(axis=0)
        for sample, value in score.items():
            rows.append({"sample": sample, "celltype": broad, "raw_celltype": raw_celltype, "score": value, "marker_count": len(genes)})
    if not rows:
        for broad, genes0 in CORE_GROUPS.items():
            genes = _present(expr, genes0)
            score = _mean_score(expr, genes)
            for sample, value in score.items():
                rows.append({"sample": sample, "celltype": broad, "raw_celltype": "curated_core_marker_fallback", "score": value, "marker_count": len(genes)})
    detail = pd.DataFrame(rows)
    broad = detail.groupby(["sample", "celltype"], as_index=False).agg(score=("score", "mean"), marker_count=("marker_count", "max"))
    wide = broad.pivot(index="sample", columns="celltype", values="score").fillna(0).reset_index()
    top_rows = []
    for sample, sub in broad.sort_values("score", ascending=False).groupby("sample"):
        names = sub.head(3)["celltype"].tolist()
        scores = sub.head(3)["score"].tolist()
        top_rows.append(
            {
                "sample": sample,
                "top_celltype_1": names[0] if len(names) > 0 else "",
                "top_celltype_2": names[1] if len(names) > 1 else "",
                "top_celltype_3": names[2] if len(names) > 2 else "",
                "top_celltype_scores": ";".join(f"{n}:{s:.3g}" for n, s in zip(names, scores)),
            }
        )
    return wide.merge(pd.DataFrame(top_rows), on="sample", how="left"), detail


def _cosine(a: np.ndarray, b: np.ndarray) -> float:
    denom = np.linalg.norm(a) * np.linalg.norm(b)
    return float(np.dot(a, b) / denom) if denom else float("nan")


def _spearman(a: np.ndarray, b: np.ndarray) -> float:
    if len(a) < 30:
        return float("nan")
    ar = pd.Series(a).rank().to_numpy()
    br = pd.Series(b).rank().to_numpy()
    return _cosine(ar - ar.mean(), br - br.mean())


def score_regions(expr: pd.DataFrame, markers: pd.DataFrame, db_path: str | Path | None = None) -> tuple[pd.DataFrame, pd.DataFrame]:
    region_markers = markers[markers["marker_type_normalized"].eq("brain_region_marker")]
    refs = read_table("brain_region_reference", db_path)
    region_cols = [c for c in refs.columns if c not in META and c != "reference_subtype"] if not refs.empty else []
    rows = []
    for region, sub in region_markers.groupby("region", dropna=True):
        genes = _present(expr, sorted(set(sub["gene_symbol"].astype(str).str.upper())))
        marker_mean = _mean_score(expr, genes)
        ref_col = region if region in region_cols else str(region).replace(" ", "_").replace("-", "_")
        for sample in expr.columns:
            cosine = float("nan")
            spearman = float("nan")
            if not refs.empty and ref_col in refs.columns:
                ref_sub = refs[["gene_symbol", ref_col]].dropna()
                ref_sub["gene_symbol"] = ref_sub["gene_symbol"].astype(str).str.upper()
                common = [g for g in ref_sub["gene_symbol"] if g in expr.index]
                if len(common) >= 10:
                    rvec = ref_sub.set_index("gene_symbol").loc[common, ref_col].astype(float).to_numpy()
                    svec = expr.loc[common, sample].astype(float).to_numpy()
                    cosine = _cosine(svec, rvec)
                    spearman = _spearman(svec, rvec)
            rows.append({"sample": sample, "region": region, "marker_mean_score": marker_mean.get(sample, 0.0), "cosine_similarity": cosine, "spearman_correlation": spearman, "marker_count": len(genes)})
    scores = pd.DataFrame(rows)
    if scores.empty:
        return scores, scores
    scores["combined_region_score"] = scores[["marker_mean_score", "cosine_similarity"]].apply(pd.to_numeric, errors="coerce").mean(axis=1)
    top_rows = []
    for sample, sub in scores.sort_values("combined_region_score", ascending=False).groupby("sample"):
        top = sub.head(3)
        row = {"sample": sample}
        for i, (_, r) in enumerate(top.iterrows(), start=1):
            row[f"top_region_{i}"] = r["region"]
            row[f"top_region_{i}_score"] = r["combined_region_score"]
        top_rows.append(row)
    return scores, pd.DataFrame(top_rows)


def score_injury(expr: pd.DataFrame, markers: pd.DataFrame) -> pd.DataFrame:
    injury = markers[markers["marker_type_normalized"].eq("tbi_injury_response")]
    genes = _present(expr, sorted(set(injury["gene_symbol"].astype(str).str.upper())))
    score = _mean_score(expr, genes)
    percentile = percentile_series(score)
    immune_score = _mean_score(expr, IMMUNE_GENES)
    oligo_genes = set(CORE_GROUPS["oligodendrocyte"])
    oligo_injury_genes = [g for g in genes if g in oligo_genes or "OLIG" in g]
    oligo_score = _mean_score(expr, oligo_injury_genes)
    return pd.DataFrame(
        {
            "sample": expr.columns,
            "injury_state_score": score.reindex(expr.columns).values,
            "injury_state_percentile": percentile.reindex(expr.columns).values,
            "injury_state_level": [level_from_percentile(p, s) for p, s in zip(percentile.reindex(expr.columns), score.reindex(expr.columns))],
            "oligodendrocyte_injury_related_score": oligo_score.reindex(expr.columns).values,
            "immune_inflammation_response_score": immune_score.reindex(expr.columns).values,
            "injury_marker_count": len(genes),
            "metadata_limit_note": "Garza TBI grouping is inferred from sample names; score is candidate injury-state evidence, not clinical diagnosis.",
        }
    )


def build_overall(contam: pd.DataFrame, brain: pd.DataFrame, cell: pd.DataFrame, region_top: pd.DataFrame, injury: pd.DataFrame, small_warning: bool) -> pd.DataFrame:
    out = brain.merge(contam, on="sample", how="outer").merge(cell, on="sample", how="left").merge(region_top, on="sample", how="left").merge(injury, on="sample", how="left")
    for col in ["top_celltype_1", "top_celltype_2", "top_celltype_3", "top_region_1", "top_region_2", "top_region_3"]:
        if col not in out.columns:
            out[col] = ""
    interpretations = []
    warnings = []
    for _, row in out.iterrows():
        flags = []
        if small_warning:
            flags.append("small cohort warning: percentile ranks are not reliable")
        if row.get("rbc_risk") == "High" or row.get("immune_risk") == "High":
            flags.append("peripheral contamination may confound brain signal interpretation")
        brain_level = row.get("brain_signal_level", "Low")
        injury_level = row.get("injury_state_level", "Low")
        if brain_level == "High" and injury_level == "High":
            interp = "strong molecular evidence of brain/injury-associated signal"
        elif brain_level == "Moderate" and (row.get("rbc_risk") == "High" or row.get("immune_risk") == "High"):
            interp = "possible brain signal, but interpretation limited by blood/immune background"
        elif brain_level == "Low":
            interp = "no robust brain-enriched cfRNA signal detected"
        else:
            interp = f"{brain_level.lower()} brain-associated signal evidence with {injury_level.lower()} injury-state evidence"
        interpretations.append(interp)
        warnings.append("; ".join(flags))
    out["overall_interpretation"] = interpretations
    out["warning_flags"] = warnings
    cols = [
        "sample", "brain_signal_score", "brain_signal_level", "rbc_score", "rbc_risk", "immune_score", "immune_risk",
        "top_celltype_1", "top_celltype_2", "top_celltype_3", "top_region_1", "top_region_2", "top_region_3",
        "injury_state_score", "injury_state_level", "overall_interpretation", "warning_flags",
    ]
    return out[cols]


def run_reference_tracing(expr_df: pd.DataFrame, gene_col: str, expression_type: str, outdir: str | Path, db_path: str | Path | None = None, make_plots: bool = True) -> dict[str, pd.DataFrame | dict | Path]:
    outdir = ensure_outdir(outdir)
    db_path = db_path or DEFAULT_DB
    expr, input_meta = prepare_expression_matrix(expr_df, gene_col, expression_type)
    markers = read_marker_candidates(db_path)
    input_meta["matched_reference_genes"] = int(len(set(expr.index) & set(read_table("gene_id_map", db_path).get("gene_symbol", pd.Series(dtype=str)).astype(str))))
    input_meta["small_cohort_warning"] = expr.shape[1] < 5
    contam = score_contamination(expr)
    brain = score_brain_signal(expr, markers)
    cell_scores, cell_detail = score_celltypes(expr, markers)
    region_scores, region_top = score_regions(expr, markers, db_path)
    injury = score_injury(expr, markers)
    overall = build_overall(contam, brain, cell_scores, region_top, injury, expr.shape[1] < 5)
    outputs = {
        "sample_contamination_scores": contam,
        "sample_brain_signal_scores": brain,
        "sample_celltype_scores": cell_scores,
        "sample_celltype_marker_detail": cell_detail,
        "sample_region_scores": region_scores,
        "sample_region_top_hits": region_top,
        "sample_injury_state_scores": injury,
        "sample_overall_tracing_summary": overall,
    }
    for name, df in outputs.items():
        df.to_csv(outdir / f"{name}.tsv", sep="\t", index=False)
    if make_plots:
        from .plotting import make_all_plots

        make_all_plots(outputs, outdir / "figures")
    from .report import write_report

    write_report(outdir / "cfrna_tracing_report.md", input_meta, outputs)
    return {"input_meta": input_meta, "outdir": outdir, **outputs}
