from __future__ import annotations

from typing import Any, Dict, Iterable, Mapping, Sequence

import numpy as np
import pandas as pd

from core.gene_utils import guess_gene_id_type

# ---------------------------------------------------------------------------
# QC risk framework notes
#
# 1) For plasma/serum cell-free miRNA, miR-451a / miR-23a-3p ratio and
#    DeltaCq(miR-23a-3p - miR-451a) are commonly used hemolysis indicators.
#    The literature often uses 5 and 7-8 as practical risk boundaries.
#
# 2) For mRNA / cfRNA-seq RBC, immune, and brain marker scores, there is no
#    universal clinical TPM cutoff. These signals should be calibrated using
#    project-specific negative controls, hemolysis-positive controls, and the
#    empirical cohort distribution.
#
# 3) The outputs here are exploratory QC risk flags. They are not clinical
#    diagnoses and should not be the sole basis for sample exclusion.
# ---------------------------------------------------------------------------

RBC_PANEL = ["HBA1", "HBA2", "HBB", "ALAS2", "CA1", "CA2", "SLC4A1", "GYPA"]
IMMUNE_PANEL = ["PTPRC", "LST1", "TYROBP", "CD74", "HLA-DRA", "S100A8", "S100A9"]
BRAIN_PANEL = ["RBFOX3", "SNAP25", "SLC17A7", "GAD1", "GFAP", "AQP4", "MBP", "MOG", "P2RY12", "TMEM119", "CLDN5"]

LOW_RISK = "Low risk"
MODERATE_RISK = "Moderate risk"
HIGH_RISK = "High risk"
UNKNOWN_RISK = "Unknown"
UNCALIBRATED_RISK = "Uncalibrated"


def _normalize_expression_df(df: pd.DataFrame) -> pd.DataFrame:
    tmp = df.copy()
    if "gene_symbol" in tmp.columns:
        tmp["gene_symbol"] = tmp["gene_symbol"].astype(str).str.strip().str.upper()
    if "tpm_value" in tmp.columns:
        tmp["tpm_value"] = pd.to_numeric(tmp["tpm_value"], errors="coerce").fillna(0.0)
    return tmp


def _aggregate_expression_df(df: pd.DataFrame) -> pd.DataFrame:
    tmp = _normalize_expression_df(df)
    if "gene_symbol" not in tmp.columns or "tpm_value" not in tmp.columns:
        return tmp
    return tmp.groupby("gene_symbol", as_index=False)["tpm_value"].mean()


def _percentile_against_reference(raw_score: float, reference: Sequence[float]) -> float:
    vals = pd.to_numeric(pd.Series(reference), errors="coerce").dropna().astype(float)
    if vals.empty or not np.isfinite(raw_score):
        return float("nan")
    return float((vals <= raw_score).mean() * 100.0)


def _zscore_against_reference(raw_score: float, reference: Sequence[float]) -> float:
    vals = pd.to_numeric(pd.Series(reference), errors="coerce").dropna().astype(float)
    if vals.empty or not np.isfinite(raw_score):
        return float("nan")
    std = float(vals.std(ddof=0))
    if std == 0:
        return 0.0
    return float((raw_score - float(vals.mean())) / std)


def _cohort_reference_from_qc(reference: Any, score_key: str) -> list[float]:
    if reference is None:
        return []
    if isinstance(reference, pd.DataFrame):
        if score_key in reference.columns:
            return pd.to_numeric(reference[score_key], errors="coerce").dropna().astype(float).tolist()
        return []
    if isinstance(reference, Mapping):
        if score_key in reference:
            return pd.to_numeric(pd.Series(reference[score_key]), errors="coerce").dropna().astype(float).tolist()
        vals = []
        for item in reference.values():
            if isinstance(item, Mapping):
                value = item.get(score_key)
                if value is not None and pd.notna(value):
                    vals.append(float(value))
        return vals
    if isinstance(reference, Sequence) and not isinstance(reference, (str, bytes)):
        vals = []
        for item in reference:
            if isinstance(item, Mapping):
                value = item.get(score_key)
                if value is not None and pd.notna(value):
                    vals.append(float(value))
        return vals
    return []


def _risk_from_percentile(percentile: float, *, high_bad: bool) -> str:
    if not np.isfinite(percentile):
        return UNCALIBRATED_RISK
    if high_bad:
        if percentile >= 90:
            return HIGH_RISK
        if percentile >= 75:
            return MODERATE_RISK
        return LOW_RISK
    if percentile < 10:
        return HIGH_RISK
    if percentile < 25:
        return MODERATE_RISK
    return LOW_RISK


def _panel_interpretation(panel_label: str, risk: str, *, high_bad: bool) -> str:
    if risk == UNKNOWN_RISK:
        return f"{panel_label} risk could not be assessed."
    if risk == UNCALIBRATED_RISK:
        return f"{panel_label} raw score is available, but cohort calibration is still needed."
    if high_bad:
        if risk == HIGH_RISK:
            return f"{panel_label} is in the extreme high tail of the cohort distribution."
        if risk == MODERATE_RISK:
            return f"{panel_label} is elevated relative to the cohort and deserves attention."
        return f"{panel_label} is within the lower-risk range of the current cohort."
    if risk == HIGH_RISK:
        return f"{panel_label} is very low relative to the cohort, so source interpretation should be cautious."
    if risk == MODERATE_RISK:
        return f"{panel_label} is somewhat weak relative to the cohort and should be interpreted cautiously."
    return f"{panel_label} is within an acceptable range for the current cohort."


def compute_marker_panel_score(
    df: pd.DataFrame,
    panel: Sequence[str],
    *,
    score_name: str,
    cohort_reference: Any = None,
    higher_percentile_means_higher_risk: bool = True,
) -> Dict[str, Any]:
    tmp = _aggregate_expression_df(df)
    out: Dict[str, Any] = {
        "panel_name": score_name,
        "panel_size": int(len(panel)),
        "matched_markers": [],
        "raw_values": {},
        "score_raw_mean": float("nan"),
        "score_raw_median": float("nan"),
        "score_raw_sum": float("nan"),
        "score_log1p_mean": float("nan"),
        "score_log1p_median": float("nan"),
        "score_log1p_sum": float("nan"),
        "detected_fraction": float("nan"),
        "percentile": float("nan"),
        "z_score": float("nan"),
        "risk": UNKNOWN_RISK,
        "interpretation": "",
    }

    if "gene_symbol" not in tmp.columns or "tpm_value" not in tmp.columns:
        out["interpretation"] = f"{score_name} could not be computed because required columns are missing."
        return out

    panel_upper = [str(g).strip().upper() for g in panel]
    s = tmp.set_index("gene_symbol")["tpm_value"]
    values = pd.Series({g: float(s.get(g, 0.0)) for g in panel_upper}, dtype=float)
    log_vals = np.log1p(values)

    out["matched_markers"] = [g for g in panel_upper if float(values.get(g, 0.0)) > 0]
    out["raw_values"] = values.to_dict()
    out["score_raw_mean"] = float(values.mean())
    out["score_raw_median"] = float(values.median())
    out["score_raw_sum"] = float(values.sum())
    out["score_log1p_mean"] = float(log_vals.mean())
    out["score_log1p_median"] = float(log_vals.median())
    out["score_log1p_sum"] = float(log_vals.sum())
    out["detected_fraction"] = float((values > 0).mean())

    cohort_scores = _cohort_reference_from_qc(cohort_reference, score_name)
    raw_score = out["score_log1p_mean"]
    if cohort_scores:
        out["percentile"] = _percentile_against_reference(raw_score, cohort_scores)
        out["z_score"] = _zscore_against_reference(raw_score, cohort_scores)
        out["risk"] = _risk_from_percentile(out["percentile"], high_bad=higher_percentile_means_higher_risk)
    else:
        out["risk"] = UNCALIBRATED_RISK if np.isfinite(raw_score) else UNKNOWN_RISK

    out["interpretation"] = _panel_interpretation(score_name, out["risk"], high_bad=higher_percentile_means_higher_risk)
    return out


def _extract_mirna_anchor_metrics(df: pd.DataFrame) -> Dict[str, Any]:
    tmp = _aggregate_expression_df(df)
    out: Dict[str, Any] = {
        "mir451a_mir23a_ratio": float("nan"),
        "mir451a_mir23a_ratio_risk": UNKNOWN_RISK,
        "mir451a_mir23a_ratio_interpretation": "",
        "delta_cq_mir23a_minus_mir451a": float("nan"),
        "delta_cq_mir23a_minus_mir451a_risk": UNKNOWN_RISK,
        "delta_cq_mir23a_minus_mir451a_interpretation": "",
    }

    if "gene_symbol" not in tmp.columns or "tpm_value" not in tmp.columns:
        return out

    s = tmp.set_index("gene_symbol")["tpm_value"]
    mir451a = float(s.get("MIR451A", s.get("MIR-451A", 0.0)))
    mir23a = float(s.get("MIR23A", s.get("MIR-23A-3P", s.get("MIR23A-3P", 0.0))))

    if mir451a > 0 or mir23a > 0:
        ratio = float(mir451a / (mir23a + 1e-8))
        out["mir451a_mir23a_ratio"] = ratio
        if ratio > 7:
            risk = HIGH_RISK
        elif ratio >= 5:
            risk = MODERATE_RISK
        else:
            risk = LOW_RISK
        out["mir451a_mir23a_ratio_risk"] = risk
        out["mir451a_mir23a_ratio_interpretation"] = (
            "Literature-anchored miR-451a / miR-23a-3p ratio suggests low hemolysis risk."
            if risk == LOW_RISK
            else "Literature-anchored miR-451a / miR-23a-3p ratio suggests moderate hemolysis risk."
            if risk == MODERATE_RISK
            else "Literature-anchored miR-451a / miR-23a-3p ratio suggests high hemolysis risk."
        )

    if "delta_cq" in df.columns and "gene_symbol" in df.columns:
        cq = df[["gene_symbol", "delta_cq"]].copy()
        cq["gene_symbol"] = cq["gene_symbol"].astype(str).str.strip().str.upper()
        cq["delta_cq"] = pd.to_numeric(cq["delta_cq"], errors="coerce")
        cq_s = cq.groupby("gene_symbol", as_index=True)["delta_cq"].mean()
        val = float(cq_s.get("MIR23A-3P", np.nan) - cq_s.get("MIR451A", np.nan)) if ("MIR23A-3P" in cq_s.index and "MIR451A" in cq_s.index) else float("nan")
        if np.isfinite(val):
            out["delta_cq_mir23a_minus_mir451a"] = val
            if val > 7:
                risk = HIGH_RISK
            elif val > 5:
                risk = MODERATE_RISK
            else:
                risk = LOW_RISK
            out["delta_cq_mir23a_minus_mir451a_risk"] = risk
            out["delta_cq_mir23a_minus_mir451a_interpretation"] = (
                "Literature-anchored DeltaCq suggests low hemolysis risk."
                if risk == LOW_RISK
                else "Literature-anchored DeltaCq suggests possible hemolysis."
                if risk == MODERATE_RISK
                else "Literature-anchored DeltaCq suggests high hemolysis risk."
            )
    return out


def grade_sample_qc_risk(qc: Mapping[str, Any]) -> str:
    if not int(qc.get("qc_applicable", 0)):
        return UNKNOWN_RISK

    has_cohort = bool(qc.get("has_cohort_reference", False))
    if not has_cohort:
        direct_mirna_risk = qc.get("hemolysis_mirna_risk")
        if direct_mirna_risk in {LOW_RISK, MODERATE_RISK, HIGH_RISK}:
            return direct_mirna_risk
        return UNCALIBRATED_RISK

    rbc_risk = qc.get("rbc_mrna_risk", UNKNOWN_RISK)
    immune_risk = qc.get("immune_mrna_risk", UNKNOWN_RISK)
    brain_risk = qc.get("brain_marker_risk", UNKNOWN_RISK)

    if HIGH_RISK in {rbc_risk, immune_risk, brain_risk}:
        return HIGH_RISK
    if MODERATE_RISK in {rbc_risk, immune_risk, brain_risk}:
        return MODERATE_RISK
    if LOW_RISK in {rbc_risk, immune_risk, brain_risk}:
        return LOW_RISK
    return UNKNOWN_RISK


def interpret_qc_risk(qc: Mapping[str, Any]) -> str:
    overall = qc.get("overall_risk", grade_sample_qc_risk(qc))
    parts = []

    hem = qc.get("hemolysis_mirna_risk")
    if hem in {LOW_RISK, MODERATE_RISK, HIGH_RISK}:
        parts.append(f"miRNA-based hemolysis anchor: {hem}.")
    else:
        parts.append(f"RBC mRNA score: {qc.get('rbc_mrna_risk', UNKNOWN_RISK)}.")

    parts.append(f"Immune background: {qc.get('immune_mrna_risk', UNKNOWN_RISK)}.")
    parts.append(f"Brain signal: {qc.get('brain_marker_risk', UNKNOWN_RISK)}.")

    if overall == HIGH_RISK:
        parts.append("Overall QC risk is high; sample interpretation should be conservative.")
    elif overall == MODERATE_RISK:
        parts.append("Overall QC risk is moderate; key results should be interpreted with added caution.")
    elif overall == LOW_RISK:
        parts.append("Overall QC risk is low under the current cohort calibration.")
    elif overall == UNCALIBRATED_RISK:
        parts.append("Only raw scores are available; cohort calibration is still needed before strong QC calls.")
    else:
        parts.append("QC risk could not be determined from the available identifiers or fields.")

    return " ".join(parts)


def compute_sample_qc(df: pd.DataFrame, cohort_reference: Any = None) -> Dict[str, Any]:
    out: Dict[str, Any] = {
        "hemolysis_hbb_hba_ratio": float("nan"),
        "immune_ptprc": float("nan"),
        "albumin_alb": float("nan"),
        "brain_signal_score": float("nan"),
        "n_total_genes": float(len(df)) if df is not None else float("nan"),
        "n_detected_genes": float("nan"),
        "gene_id_type": "unknown",
        "qc_applicable": 0.0,
        "has_cohort_reference": False,
        "rbc_mrna_score": float("nan"),
        "rbc_mrna_detected_fraction": float("nan"),
        "rbc_mrna_percentile": float("nan"),
        "rbc_mrna_z_score": float("nan"),
        "rbc_mrna_risk": UNKNOWN_RISK,
        "rbc_mrna_interpretation": "",
        "immune_mrna_score": float("nan"),
        "immune_mrna_detected_fraction": float("nan"),
        "immune_mrna_percentile": float("nan"),
        "immune_mrna_z_score": float("nan"),
        "immune_mrna_risk": UNKNOWN_RISK,
        "immune_mrna_interpretation": "",
        "brain_marker_score": float("nan"),
        "brain_marker_detected_fraction": float("nan"),
        "brain_marker_percentile": float("nan"),
        "brain_marker_z_score": float("nan"),
        "brain_marker_risk": UNKNOWN_RISK,
        "brain_marker_interpretation": "",
        "hemolysis_mirna_risk": UNKNOWN_RISK,
        "overall_risk": UNKNOWN_RISK,
        "interpretation": "",
    }

    if df is None or "gene_symbol" not in df.columns or "tpm_value" not in df.columns:
        return out

    tmp = _aggregate_expression_df(df)
    out["n_total_genes"] = float(len(tmp))
    out["n_detected_genes"] = float(int((pd.to_numeric(tmp["tpm_value"], errors="coerce").fillna(0) > 0).sum()))
    out["gene_id_type"] = guess_gene_id_type(tmp["gene_symbol"].astype(str).tolist())
    if out["gene_id_type"] != "symbol_like":
        out["interpretation"] = interpret_qc_risk(out)
        return out

    out["qc_applicable"] = 1.0
    out["has_cohort_reference"] = bool(_cohort_reference_from_qc(cohort_reference, "rbc_mrna_score"))

    s = tmp.set_index("gene_symbol")["tpm_value"]
    out["hemolysis_hbb_hba_ratio"] = float(float(s.get("HBB", 0.0)) / (float(s.get("HBA1", 0.0)) + float(s.get("HBA2", 0.0)) + 1e-8))
    out["immune_ptprc"] = float(s.get("PTPRC", 0.0))
    out["albumin_alb"] = float(s.get("ALB", 0.0))

    rbc_panel = compute_marker_panel_score(
        tmp,
        RBC_PANEL,
        score_name="rbc_mrna_score",
        cohort_reference=cohort_reference,
        higher_percentile_means_higher_risk=True,
    )
    immune_panel = compute_marker_panel_score(
        tmp,
        IMMUNE_PANEL,
        score_name="immune_mrna_score",
        cohort_reference=cohort_reference,
        higher_percentile_means_higher_risk=True,
    )
    brain_panel = compute_marker_panel_score(
        tmp,
        BRAIN_PANEL,
        score_name="brain_marker_score",
        cohort_reference=cohort_reference,
        higher_percentile_means_higher_risk=False,
    )

    out["rbc_mrna_score"] = rbc_panel["score_log1p_mean"]
    out["rbc_mrna_detected_fraction"] = rbc_panel["detected_fraction"]
    out["rbc_mrna_percentile"] = rbc_panel["percentile"]
    out["rbc_mrna_z_score"] = rbc_panel["z_score"]
    out["rbc_mrna_risk"] = rbc_panel["risk"]
    out["rbc_mrna_interpretation"] = rbc_panel["interpretation"]

    out["immune_mrna_score"] = immune_panel["score_log1p_mean"]
    out["immune_mrna_detected_fraction"] = immune_panel["detected_fraction"]
    out["immune_mrna_percentile"] = immune_panel["percentile"]
    out["immune_mrna_z_score"] = immune_panel["z_score"]
    out["immune_mrna_risk"] = immune_panel["risk"]
    out["immune_mrna_interpretation"] = immune_panel["interpretation"]

    out["brain_marker_score"] = brain_panel["score_log1p_mean"]
    out["brain_marker_detected_fraction"] = brain_panel["detected_fraction"]
    out["brain_marker_percentile"] = brain_panel["percentile"]
    out["brain_marker_z_score"] = brain_panel["z_score"]
    out["brain_marker_risk"] = brain_panel["risk"]
    out["brain_marker_interpretation"] = brain_panel["interpretation"]
    out["brain_signal_score"] = out["brain_marker_score"]

    mirna_metrics = _extract_mirna_anchor_metrics(df)
    out.update(mirna_metrics)
    if mirna_metrics.get("mir451a_mir23a_ratio_risk") in {LOW_RISK, MODERATE_RISK, HIGH_RISK}:
        out["hemolysis_mirna_risk"] = mirna_metrics["mir451a_mir23a_ratio_risk"]
    elif mirna_metrics.get("delta_cq_mir23a_minus_mir451a_risk") in {LOW_RISK, MODERATE_RISK, HIGH_RISK}:
        out["hemolysis_mirna_risk"] = mirna_metrics["delta_cq_mir23a_minus_mir451a_risk"]

    out["overall_risk"] = grade_sample_qc_risk(out)
    out["interpretation"] = interpret_qc_risk(out)
    return out


def compute_cohort_qc(samples: Any) -> Any:
    if isinstance(samples, Mapping):
        raw = {sample_id: compute_sample_qc(df) for sample_id, df in samples.items()}
        ref = {
            "rbc_mrna_score": [v.get("rbc_mrna_score") for v in raw.values()],
            "immune_mrna_score": [v.get("immune_mrna_score") for v in raw.values()],
            "brain_marker_score": [v.get("brain_marker_score") for v in raw.values()],
        }
        return {sample_id: compute_sample_qc(df, cohort_reference=ref) for sample_id, df in samples.items()}

    if isinstance(samples, Sequence) and not isinstance(samples, (str, bytes)):
        raw_list = [compute_sample_qc(df) for df in samples]
        ref = {
            "rbc_mrna_score": [v.get("rbc_mrna_score") for v in raw_list],
            "immune_mrna_score": [v.get("immune_mrna_score") for v in raw_list],
            "brain_marker_score": [v.get("brain_marker_score") for v in raw_list],
        }
        return [compute_sample_qc(df, cohort_reference=ref) for df in samples]

    raise TypeError("samples must be a mapping of sample_id -> DataFrame or a sequence of DataFrame objects.")


def grade_sample_qc(qc: Dict[str, Any]) -> str:
    return grade_sample_qc_risk(qc)
