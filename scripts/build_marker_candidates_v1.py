#!/usr/bin/env python
"""Build first-version marker candidate tables from processed references."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
V1 = ROOT / "data" / "processed_reference" / "v1"
OUT = V1 / "marker_candidates"
META = {"gene_symbol", "gene_symbol_raw", "ensembl_id", "source_atlas", "doi", "recommended_use"}


RBC_GENES = ["HBA1", "HBA2", "HBB", "ALAS2", "CA1", "CA2", "SLC4A1", "GYPA"]
IMMUNE_GENES = ["PTPRC", "LST1", "TYROBP", "CD74", "HLA-DRA", "S100A8", "S100A9"]


def read_matrix(name: str) -> pd.DataFrame:
    return pd.read_csv(V1 / name, sep="\t")


def expr_cols(df: pd.DataFrame) -> list[str]:
    return [c for c in df.columns if c not in META]


def add_expression_scale(df: pd.DataFrame, cols: list[str]) -> pd.DataFrame:
    out = df.copy()
    out[cols] = np.expm1(out[cols].apply(pd.to_numeric, errors="coerce"))
    return out


def top_markers(df: pd.DataFrame, group_cols: list[str], group_label: str, min_expr: float = 0.5) -> pd.DataFrame:
    raw = add_expression_scale(df, group_cols)
    rows = []
    all_values = raw[group_cols]
    for group in group_cols:
        target = pd.to_numeric(raw[group], errors="coerce")
        others = all_values[[c for c in group_cols if c != group]].mean(axis=1)
        rows.append(
            pd.DataFrame(
                {
                    "gene_symbol": raw["gene_symbol"],
                    "gene_symbol_raw": raw["gene_symbol_raw"],
                    "ensembl_id": raw.get("ensembl_id", ""),
                    group_label: group,
                    "target_mean_expression": target,
                    "background_mean_expression": others,
                    "log2FC": np.log2((target + 0.1) / (others + 0.1)),
                    "detection_fraction": np.nan,
                    "source_atlas": raw.get("source_atlas", ""),
                    "doi": raw.get("doi", ""),
                    "recommended_use": raw.get("recommended_use", ""),
                }
            )
        )
    out = pd.concat(rows, ignore_index=True)
    out = out[(out["target_mean_expression"] >= min_expr) & (out["log2FC"] > 1)]
    return out.sort_values([group_label, "log2FC"], ascending=[True, False])


def build_brain_enriched(report: list[str]) -> pd.DataFrame:
    gtex = read_matrix("peripheral_background_gtex_median_tpm.tsv")
    cols = expr_cols(gtex)
    brain_cols = [c for c in cols if c.lower().startswith("brain -")]
    peripheral_cols = [c for c in cols if c not in brain_cols]
    raw = add_expression_scale(gtex, cols)
    raw["brain_median_tpm"] = raw[brain_cols].median(axis=1)
    raw["max_peripheral_median_tpm"] = raw[peripheral_cols].max(axis=1)
    raw["max_peripheral_tissue"] = raw[peripheral_cols].idxmax(axis=1)
    raw["brain_vs_peripheral_ratio"] = raw["brain_median_tpm"] / raw["max_peripheral_median_tpm"].replace(0, np.nan)
    out = raw[
        (raw["brain_median_tpm"] > 1)
        & (raw["brain_vs_peripheral_ratio"] > 5)
        & (raw["max_peripheral_tissue"].astype(str) != "Whole Blood")
    ][
        [
            "gene_symbol",
            "gene_symbol_raw",
            "ensembl_id",
            "brain_median_tpm",
            "max_peripheral_median_tpm",
            "max_peripheral_tissue",
            "brain_vs_peripheral_ratio",
            "source_atlas",
            "doi",
            "recommended_use",
        ]
    ].sort_values("brain_vs_peripheral_ratio", ascending=False)
    out.to_csv(OUT / "brain_enriched_genes_vs_gtex.tsv", sep="\t", index=False)
    report.append(f"- Brain-enriched genes vs GTEx: {len(out)} candidates.")
    return out


def build_contamination(report: list[str]) -> pd.DataFrame:
    rows = [{"gene_symbol": g, "contamination_class": "rbc", "source_atlas": "curated", "doi": "", "recommended_use": "contamination_marker"} for g in RBC_GENES]
    rows += [{"gene_symbol": g, "contamination_class": "immune", "source_atlas": "curated", "doi": "", "recommended_use": "contamination_marker"} for g in IMMUNE_GENES]
    out = pd.DataFrame(rows)
    out.to_csv(OUT / "blood_immune_rbc_contamination_genes.tsv", sep="\t", index=False)
    report.append(f"- Contamination markers: {len(out)} curated RBC/immune genes.")
    return out


def build_region(report: list[str]) -> pd.DataFrame:
    frames = []
    for fname in ["brain_region_reference_hpa.tsv", "brain_region_reference_allen_hba.tsv"]:
        path = V1 / fname
        if not path.exists():
            continue
        df = read_matrix(fname)
        frames.append(top_markers(df, expr_cols(df), "region", min_expr=0.5))
    out = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    out.to_csv(OUT / "brain_region_marker_candidates.tsv", sep="\t", index=False)
    report.append(f"- Brain-region marker candidates: {len(out)} rows.")
    return out


def build_celltype(report: list[str]) -> pd.DataFrame:
    df = read_matrix("celltype_reference_hpa_single_nuclei_brain.tsv")
    out = top_markers(df, expr_cols(df), "celltype", min_expr=0.5)
    out.to_csv(OUT / "brain_celltype_marker_candidates.tsv", sep="\t", index=False)
    report.append(f"- Brain cell-type marker candidates: {len(out)} rows from HPA single-nuclei brain.")
    return out


def build_tbi(report: list[str]) -> pd.DataFrame:
    df = read_matrix("injury_state_reference_garza2023_tbi.tsv")
    cols = expr_cols(df)
    raw = add_expression_scale(df, cols)
    tbi_cols = [c for c in cols if c.startswith("TBI|") and "ctrl" not in c.lower()]
    control_cols = [c for c in cols if "ctrl" in c.lower() or c.startswith("hGPC|")]
    if tbi_cols and control_cols:
        target = raw[tbi_cols].mean(axis=1)
        background = raw[control_cols].mean(axis=1)
        out = pd.DataFrame(
            {
                "gene_symbol": raw["gene_symbol"],
                "gene_symbol_raw": raw["gene_symbol_raw"],
                "ensembl_id": raw.get("ensembl_id", ""),
                "target_mean_expression": target,
                "background_mean_expression": background,
                "log2FC": np.log2((target + 0.1) / (background + 0.1)),
                "source_atlas": raw.get("source_atlas", ""),
                "doi": raw.get("doi", ""),
                "recommended_use": raw.get("recommended_use", ""),
                "comparison_note": "TBI non-control samples vs TBI control and hGPC columns inferred from sample names",
            }
        )
        out = out[(out["target_mean_expression"] >= 0.5) & (out["log2FC"] > 1)].sort_values("log2FC", ascending=False)
        report.append(f"- TBI injury-response candidates: {len(out)} rows using sample-name inferred TBI/control groups.")
    else:
        raw["mean_expression"] = raw[cols].mean(axis=1)
        out = raw.sort_values("mean_expression", ascending=False).head(5000)
        out["comparison_note"] = "No usable TBI/control grouping detected; high-expression candidates only."
        report.append("- TBI injury-response candidates: no metadata grouping available; high-expression fallback used.")
    out.to_csv(OUT / "tbi_injury_response_candidates.tsv", sep="\t", index=False)
    return out


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    report = ["# Marker Build Report", ""]
    build_brain_enriched(report)
    build_region(report)
    build_celltype(report)
    build_tbi(report)
    build_contamination(report)
    report.extend(
        [
            "",
            "## Notes",
            "",
            "- GTEx median TPM is used as peripheral background.",
            "- Region and cell-type marker tables use group-vs-other log2FC from processed log1p matrices after converting back to expression scale.",
            "- Detection fraction is left blank unless the upstream atlas provides detection statistics.",
        ]
    )
    (OUT / "marker_build_report.md").write_text("\n".join(report), encoding="utf-8")
    print(f"Wrote marker candidates to {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
