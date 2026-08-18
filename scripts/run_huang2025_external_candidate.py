#!/usr/bin/env python
"""Run the provenance-remediated Huang 2025 cfRNA external-domain audit.

The public Huang expression matrix is used only as a 159-profile computational
stress-test resource.  Its labels encode disease group and biofluid, not a
validated patient identifier or CSF--plasma correspondence.  This runner
therefore analyses all 77 CSF and 82 plasma profiles as separate,
fluid-specific profile cohorts and contains no patient-paired or synthetic-mixture
calculation.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, spearmanr


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from core.bo2023_region_tracing import ROUTE_NAME, trace_bo2023_secondary_regions  # noqa: E402
from core.network_tracing import (  # noqa: E402
    DEFAULT_BO2023_NETWORK_MODEL,
    DEFAULT_BO2023_REFERENCE_PROJECTOR,
    trace_network_expression,
)


DEFAULT_INPUT = ROOT / "external_inputs" / "huang2025_pmc12041490" / "41698_2025_909_MOESM2_ESM.csv"
DEFAULT_OUTDIR = ROOT / "reproducibility" / "huang_2025"
DEFAULT_DB_PATH = ROOT / "braintrace_source_tracing.db"
SAMPLE_RE = re.compile(r"^(?P<disease>GLI|MEN|NOR)_(?P<fluid>CSF|plasma)\d+$")
DISEASE_LABEL = {"GLI": "glioma", "MEN": "meningioma", "NOR": "control"}
EXPECTED_COUNTS = {
    ("GLI", "CSF"): 16,
    ("GLI", "plasma"): 18,
    ("MEN", "CSF"): 43,
    ("MEN", "plasma"): 46,
    ("NOR", "CSF"): 18,
    ("NOR", "plasma"): 18,
}
EXPECTED_TOTAL = 159
EXPECTED_CSF = 77
EXPECTED_PLASMA = 82
N_BOOTSTRAP = 10_000
RANDOM_SEED = 20250719

SOURCE_DOI = "10.1038/s41698-025-00909-6"
SOURCE_ARTICLE_URL = "https://doi.org/10.1038/s41698-025-00909-6"
SOURCE_QC_NOTE = (
    "The source study reports aggregate sequencing-QC exclusions for its clinical analyses "
    "(five CSF and one plasma profile); no public per-profile QC-status mapping accompanies "
    "the expression matrix."
)
PATIENT_ID_STATUS = "unknown_not_supplied_in_public_expression_matrix"
OMPFC_NETWORK = "Orbitomedial Prefrontal Cortex (OMPFC)"
PLATELET_MARKERS = ["PF4", "PPBP", "RGS18", "GP9", "ITGA2B", "TUBB1", "SELP", "NRGN"]
EV_MARKERS = ["CD9", "CD63", "CD81", "TSG101", "PDCD6IP", "SDCBP", "FLOT1", "FLOT2"]


def json_default(value: Any) -> Any:
    """Convert numpy/pandas values without silently changing their numerical value."""

    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, np.floating):
        return None if math.isnan(float(value)) else float(value)
    if isinstance(value, np.bool_):
        return bool(value)
    if value is pd.NA:
        return None
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, (pd.Timestamp, datetime)):
        return value.isoformat()
    raise TypeError(f"Cannot serialise {type(value).__name__}")


def write_json(path: Path, payload: Any) -> None:
    path.write_text(
        json.dumps(payload, indent=2, ensure_ascii=False, default=json_default) + "\n",
        encoding="utf-8",
        newline="\n",
    )


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def asset_record(path: Path | None, label: str) -> dict[str, Any]:
    """Describe a supplied provenance asset without inventing a missing input."""

    if path is None:
        return {"label": label, "provided": False}
    resolved = path.resolve()
    try:
        locator = resolved.relative_to(ROOT.resolve()).as_posix()
        path_kind = "repository_relative"
    except ValueError:
        locator = resolved.name
        path_kind = "external_basename"
    result: dict[str, Any] = {
        "label": label,
        "provided": True,
        "path": locator,
        "path_kind": path_kind,
        "exists": resolved.exists(),
    }
    if resolved.exists() and resolved.is_file():
        result.update({"bytes": resolved.stat().st_size, "sha256": sha256_file(resolved)})
    return result


def parse_sample_id(sample_id: str) -> dict[str, Any]:
    """Parse only the public disease and fluid labels; never infer patient identity."""

    match = SAMPLE_RE.fullmatch(str(sample_id))
    if not match:
        raise ValueError(f"Unexpected Huang 2025 sample identifier: {sample_id!r}")
    disease_code = match.group("disease")
    fluid = match.group("fluid")
    return {
        "sample_id": str(sample_id),
        "disease_code": disease_code,
        "disease_group": DISEASE_LABEL[disease_code],
        "fluid": fluid,
        "tumor_status": "control" if disease_code == "NOR" else "tumour",
        "patient_id": pd.NA,
        "patient_id_status": PATIENT_ID_STATUS,
    }


def source_log2_rpm_to_log1p_rpm(values: np.ndarray) -> np.ndarray:
    numeric = np.asarray(values, dtype=np.float64)
    if not np.isfinite(numeric).all() or (numeric < 0).any():
        raise ValueError("Huang matrix contains negative or non-finite values")
    return (numeric * math.log(2.0)).astype(np.float32)


def normalized_entropy(probabilities: list[float]) -> float:
    values = np.asarray(probabilities, dtype=float)
    values = np.clip(values, 0.0, None)
    total = float(values.sum())
    if total <= 0 or len(values) <= 1:
        return float("nan")
    values = values / total
    nonzero = values[values > 0]
    return float(-(nonzero * np.log(nonzero)).sum() / np.log(len(values)))


def cliffs_delta(left: np.ndarray, right: np.ndarray) -> float:
    x = np.asarray(left, dtype=float)
    y = np.asarray(right, dtype=float)
    return float((np.greater.outer(x, y).sum() - np.less.outer(x, y).sum()) / (len(x) * len(y)))


def bootstrap_median_difference(
    left: np.ndarray, right: np.ndarray, *, n_resamples: int, seed: int
) -> tuple[float, float, float]:
    x = np.asarray(left, dtype=float)
    y = np.asarray(right, dtype=float)
    rng = np.random.default_rng(seed)
    differences = np.empty(n_resamples, dtype=float)
    for index in range(n_resamples):
        differences[index] = float(
            np.median(rng.choice(x, size=len(x), replace=True))
            - np.median(rng.choice(y, size=len(y), replace=True))
        )
    return (
        float(np.median(x) - np.median(y)),
        float(np.quantile(differences, 0.025)),
        float(np.quantile(differences, 0.975)),
    )


def benjamini_hochberg(p_values: Iterable[float]) -> list[float]:
    values = np.asarray(list(p_values), dtype=float)
    order = np.argsort(values)
    adjusted = np.empty(len(values), dtype=float)
    running = 1.0
    for reverse_rank, index in enumerate(order[::-1], start=1):
        rank = len(values) - reverse_rank + 1
        running = min(running, float(values[index]) * len(values) / rank)
        adjusted[index] = min(running, 1.0)
    return adjusted.tolist()


def top_ids(rows: list[dict[str, Any]], key: str, k: int = 3) -> list[str]:
    values = [str(row.get(key, "")) for row in rows[:k] if str(row.get(key, "")).strip()]
    if len(values) < k:
        raise ValueError(f"Expected {k} nonempty {key} values, found {len(values)}")
    return values


def read_expression(path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Read the public matrix and retain the locked model and marker genes."""

    header = pd.read_csv(path, nrows=0).columns.astype(str).tolist()
    if len(header) != 160:
        raise ValueError(f"Expected 160 columns, found {len(header)}")
    sample_ids = header[1:]
    metadata = pd.DataFrame([parse_sample_id(sample_id) for sample_id in sample_ids])
    observed = metadata.groupby(["disease_code", "fluid"]).size().to_dict()
    if observed != EXPECTED_COUNTS:
        raise ValueError(f"Unexpected disease-by-fluid sample counts: {observed!r}")
    if len(metadata) != EXPECTED_TOTAL:
        raise ValueError(f"Expected {EXPECTED_TOTAL} profiles, found {len(metadata)}")
    if int(metadata["fluid"].eq("CSF").sum()) != EXPECTED_CSF or int(metadata["fluid"].eq("plasma").sum()) != EXPECTED_PLASMA:
        raise ValueError("Unexpected aggregate CSF/plasma count")

    with np.load(DEFAULT_BO2023_REFERENCE_PROJECTOR, allow_pickle=False) as payload:
        projector_genes = set(payload["genes"].astype(str))
    with np.load(DEFAULT_BO2023_NETWORK_MODEL, allow_pickle=False) as payload:
        network_genes = set(payload["genes"].astype(str))
    required_genes = projector_genes | network_genes | set(PLATELET_MARKERS) | set(EV_MARKERS)

    raw = pd.read_csv(path, dtype={header[0]: str}, low_memory=False)
    raw = raw.rename(columns={header[0]: "gene_symbol"})
    if raw.shape != (83_929, 160):
        raise ValueError(f"Expected 83,929 x 160 matrix, found {raw.shape}")
    raw["gene_symbol"] = raw["gene_symbol"].astype(str).str.strip()
    if raw["gene_symbol"].eq("").any() or raw["gene_symbol"].duplicated().any():
        raise ValueError("Gene symbols must be nonblank and unique")
    selected = raw.loc[raw["gene_symbol"].isin(required_genes)].copy()
    values = selected[sample_ids].to_numpy(dtype=np.float64)
    selected.loc[:, sample_ids] = source_log2_rpm_to_log1p_rpm(values)
    return selected.set_index("gene_symbol").astype(np.float32), metadata


def compare_tumour_control(sample_df: pd.DataFrame) -> pd.DataFrame:
    """Run six two-sided Mann--Whitney U tests at the profile level when pairing metadata are unavailable."""

    rows: list[dict[str, Any]] = []
    metrics = ["atlas_fit_score", "network_margin", "network_entropy"]
    counter = 0
    for fluid in ["CSF", "plasma"]:
        subset = sample_df.loc[sample_df["fluid"].eq(fluid)]
        tumour = subset.loc[subset["tumor_status"].eq("tumour")]
        control = subset.loc[subset["tumor_status"].eq("control")]
        for metric in metrics:
            counter += 1
            x = tumour[metric].to_numpy(dtype=float)
            y = control[metric].to_numpy(dtype=float)
            statistic, raw_p = mannwhitneyu(x, y, alternative="two-sided")
            difference, ci_low, ci_high = bootstrap_median_difference(
                x, y, n_resamples=N_BOOTSTRAP, seed=RANDOM_SEED + counter
            )
            rows.append(
                {
                    "fluid": fluid,
                    "comparison": "tumour_vs_control",
                    "test": "two-sided Mann-Whitney U (profile-level; pairing unavailable)",
                    "metric": metric,
                    "n_tumour": len(x),
                    "n_control": len(y),
                    "median_tumour": float(np.median(x)),
                    "median_control": float(np.median(y)),
                    "median_difference_tumour_minus_control": difference,
                    "bootstrap_ci95_low": ci_low,
                    "bootstrap_ci95_high": ci_high,
                    "mann_whitney_u": float(statistic),
                    "raw_p": float(raw_p),
                    "cliffs_delta": cliffs_delta(x, y),
                }
            )
    frame = pd.DataFrame(rows)
    frame["bh_fdr"] = benjamini_hochberg(frame["raw_p"].tolist())
    return frame


def marker_correlations(sample_df: pd.DataFrame) -> pd.DataFrame:
    """Calculate descriptive, fluid-specific marker correlations with OMPFC score."""

    rows: list[dict[str, Any]] = []
    for fluid in ["CSF", "plasma"]:
        subset = sample_df.loc[sample_df["fluid"].eq(fluid)]
        for marker_class, markers in [
            ("platelet-associated", PLATELET_MARKERS),
            ("extracellular-vesicle-associated", EV_MARKERS),
        ]:
            for marker in markers:
                marker_column = f"expr_{marker}"
                pair = subset[["OMPFC_network_score", marker_column]].dropna()
                if len(pair) >= 3 and pair.iloc[:, 0].nunique() > 1 and pair.iloc[:, 1].nunique() > 1:
                    rho, raw_p = spearmanr(pair.iloc[:, 0], pair.iloc[:, 1])
                    status = "estimated"
                else:
                    rho, raw_p, status = np.nan, np.nan, "not_estimable"
                rows.append(
                    {
                        "fluid": fluid,
                        "marker_class": marker_class,
                        "marker": marker,
                        "target": "OMPFC network score",
                        "n": int(len(pair)),
                        "spearman_rho": float(rho) if pd.notna(rho) else np.nan,
                        "raw_p": float(raw_p) if pd.notna(raw_p) else np.nan,
                        "status": status,
                    }
                )
    result = pd.DataFrame(rows)
    estimable = result["raw_p"].notna()
    result.loc[estimable, "bh_fdr"] = benjamini_hochberg(result.loc[estimable, "raw_p"].tolist())
    return result


def distribution_table(sample_df: pd.DataFrame, rank_depth: int) -> pd.DataFrame:
    """Report network identities within predeclared independent cohorts."""

    cohorts = {
        "all_159_profiles": pd.Series(True, index=sample_df.index),
        "all_CSF_profiles": sample_df["fluid"].eq("CSF"),
        "all_plasma_profiles": sample_df["fluid"].eq("plasma"),
        "CSF_tumour_profiles": sample_df["fluid"].eq("CSF") & sample_df["tumor_status"].eq("tumour"),
        "CSF_control_profiles": sample_df["fluid"].eq("CSF") & sample_df["tumor_status"].eq("control"),
        "plasma_tumour_profiles": sample_df["fluid"].eq("plasma") & sample_df["tumor_status"].eq("tumour"),
        "plasma_control_profiles": sample_df["fluid"].eq("plasma") & sample_df["tumor_status"].eq("control"),
    }
    rank_columns = [f"network_top{rank}" for rank in range(1, rank_depth + 1)]
    rows: list[dict[str, Any]] = []
    for cohort, mask in cohorts.items():
        subset = sample_df.loc[mask]
        n_profiles = len(subset)
        networks = pd.concat([subset[column] for column in rank_columns], ignore_index=True)
        for network, count in networks.value_counts(dropna=True).items():
            rows.append(
                {
                    "cohort": cohort,
                    "rank_depth": rank_depth,
                    "network": network,
                    "count": int(count),
                    "denominator_profiles": n_profiles,
                    "percent_profiles": float(100 * count / n_profiles),
                }
            )
    return pd.DataFrame(rows).sort_values(["cohort", "count", "network"], ascending=[True, False, True]).reset_index(drop=True)


def fluid_summary(sample_df: pd.DataFrame) -> pd.DataFrame:
    """Summarise all/full-fluid and tumour/control distributions without pairing."""

    cohorts = {
        "all_159_profiles": pd.Series(True, index=sample_df.index),
        "all_CSF_profiles": sample_df["fluid"].eq("CSF"),
        "all_plasma_profiles": sample_df["fluid"].eq("plasma"),
        "CSF_tumour_profiles": sample_df["fluid"].eq("CSF") & sample_df["tumor_status"].eq("tumour"),
        "CSF_control_profiles": sample_df["fluid"].eq("CSF") & sample_df["tumor_status"].eq("control"),
        "plasma_tumour_profiles": sample_df["fluid"].eq("plasma") & sample_df["tumor_status"].eq("tumour"),
        "plasma_control_profiles": sample_df["fluid"].eq("plasma") & sample_df["tumor_status"].eq("control"),
    }
    rows: list[dict[str, Any]] = []
    for cohort, mask in cohorts.items():
        subset = sample_df.loc[mask]
        n_profiles = int(len(subset))
        top1 = int(subset["network_top1"].eq(OMPFC_NETWORK).sum())
        top3 = int(subset[["network_top1", "network_top2", "network_top3"]].eq(OMPFC_NETWORK).any(axis=1).sum())
        rows.append(
            {
                "cohort": cohort,
                "n_profiles": n_profiles,
                "OMPFC_top1_numerator": top1,
                "OMPFC_top1_denominator": n_profiles,
                "OMPFC_top1_percent": float(100 * top1 / n_profiles),
                "OMPFC_top3_numerator": top3,
                "OMPFC_top3_denominator": n_profiles,
                "OMPFC_top3_percent": float(100 * top3 / n_profiles),
                "median_atlas_fit_score": float(subset["atlas_fit_score"].median()),
                "median_network_margin": float(subset["network_margin"].median()),
                "median_network_entropy": float(subset["network_entropy"].median()),
                "median_OMPFC_network_score": float(subset["OMPFC_network_score"].median()),
            }
        )
    return pd.DataFrame(rows)


def build_ledger(metadata: pd.DataFrame, sample_df: pd.DataFrame) -> pd.DataFrame:
    """Produce the prescribed per-profile inclusion ledger."""

    available = sample_df.set_index("sample_id")["BrainTrace_output_available"]
    ledger = metadata.set_index("sample_id").copy()
    ledger["expression_available"] = True
    ledger["BrainTrace_output_available"] = available
    ledger["included_in_full159_audit"] = True
    ledger["included_in_CSF_analysis"] = ledger["fluid"].eq("CSF")
    ledger["included_in_plasma_analysis"] = ledger["fluid"].eq("plasma")
    ledger["included_in_tumour_control_analysis"] = True
    ledger["source_QC_status_if_known"] = "not_publicly_mapped_per_profile"
    ledger["source_QC_note"] = SOURCE_QC_NOTE
    return ledger.reset_index().loc[
        :,
        [
            "sample_id",
            "fluid",
            "disease_group",
            "tumor_status",
            "patient_id",
            "patient_id_status",
            "expression_available",
            "BrainTrace_output_available",
            "included_in_full159_audit",
            "included_in_CSF_analysis",
            "included_in_plasma_analysis",
            "included_in_tumour_control_analysis",
            "source_QC_status_if_known",
            "source_QC_note",
        ],
    ]


def validate_canonical_invariants(
    metadata: pd.DataFrame, ledger: pd.DataFrame, sample_df: pd.DataFrame, comparisons: pd.DataFrame
) -> None:
    """Fail closed on the provenance and denominator commitments of this audit."""

    if len(metadata) != EXPECTED_TOTAL or len(ledger) != EXPECTED_TOTAL or len(sample_df) != EXPECTED_TOTAL:
        raise ValueError("The canonical audit requires all 159 published profiles.")
    if int(metadata["fluid"].eq("CSF").sum()) != EXPECTED_CSF or int(metadata["fluid"].eq("plasma").sum()) != EXPECTED_PLASMA:
        raise ValueError("The canonical audit requires 77 CSF and 82 plasma profiles.")
    if not ledger["patient_id"].isna().all() or not ledger["patient_id_status"].eq(PATIENT_ID_STATUS).all():
        raise ValueError("Patient identity must remain unavailable in the public-matrix audit.")
    if not sample_df["BrainTrace_output_available"].all():
        raise ValueError("Every published profile must yield a traceable canonical output.")
    if len(comparisons) != 6 or not comparisons["test"].eq("two-sided Mann-Whitney U (profile-level; pairing unavailable)").all():
        raise ValueError("Exactly six independent tumour-control tests are required.")
    if set(comparisons["n_tumour"].astype(int)) != {59, 64} or set(comparisons["n_control"].astype(int)) != {18}:
        raise ValueError("Tumour-control denominators are inconsistent with the complete fluid cohorts.")


def make_summary(
    metadata: pd.DataFrame,
    sample_df: pd.DataFrame,
    fluid: pd.DataFrame,
    comparisons: pd.DataFrame,
    correlations: pd.DataFrame,
) -> dict[str, Any]:
    return {
        "protocol_status": "huang_2025_provenance_remediated",
        "source_doi": SOURCE_DOI,
        "source_article_url": SOURCE_ARTICLE_URL,
        "n_profiles": int(len(metadata)),
        "n_csf": int(metadata["fluid"].eq("CSF").sum()),
        "n_plasma": int(metadata["fluid"].eq("plasma").sum()),
        "n_traceable_outputs": int(sample_df["BrainTrace_output_available"].sum()),
        "traceable_output_percent": float(100 * sample_df["BrainTrace_output_available"].mean()),
        "patient_id_metadata": "not_available_in_public_expression_matrix",
        "patient_paired_analysis": "NOT_SUPPORTED",
        "synthetic_matched_admixture": "REMOVED_FROM_CANONICAL_ANALYSIS",
        "analysis_unit": "profile-level observation within fluid-specific cohort; patient-level dependence cannot be assessed from the public matrix",
        "source_clinical_qc_note": SOURCE_QC_NOTE,
        "interpretation": {
            "supported_claim": "technical portability and domain-shift audit",
            "not_supported_claims": [
                "patient-level CSF-plasma correspondence",
                "synthetic matched-mixture behavior",
                "anatomical localization accuracy",
                "clinical validity",
            ],
        },
        "fluid_summary": fluid.to_dict(orient="records"),
        "independent_tumour_control_tests": comparisons.to_dict(orient="records"),
        "minimum_bh_fdr": float(comparisons["bh_fdr"].min()),
        "marker_correlations": correlations.to_dict(orient="records"),
    }


def write_results_markdown(
    path: Path,
    summary: dict[str, Any],
    fluid: pd.DataFrame,
    comparisons: pd.DataFrame,
    correlations: pd.DataFrame,
) -> None:
    """Write a prose report bounded to supported technical-transfer claims."""

    csf = fluid.loc[fluid["cohort"].eq("all_CSF_profiles")].iloc[0]
    plasma = fluid.loc[fluid["cohort"].eq("all_plasma_profiles")].iloc[0]
    lines = [
        "# Huang 2025 cfRNA: provenance-remediated external domain audit",
        "",
        f"Source: Huang et al. (2025), DOI [{SOURCE_DOI}]({SOURCE_ARTICLE_URL}).",
        "",
        "## Scope and provenance",
        "",
        "The public expression matrix was used as a computational stress-test resource. Although the source study applied its own sequencing-QC exclusions for its clinical analyses, all 159 profiles available in the published matrix were considered here for technical transfer auditing; this analysis does not attempt to reproduce the source study’s QC-filtered clinical analysis.",
        "",
        "No patient-level CSF-plasma correspondence was assumed. CSF and plasma were analysed as independent fluid-specific cohorts. Sample-label suffixes were not interpreted as patient identifiers.",
        "",
        "This audit supports a narrow technical-portability/domain-shift statement. It does not establish patient correspondence, synthetic mixture behavior, anatomical localization accuracy, tumour-source discrimination, or clinical validity.",
        "",
        "## Cohort accounting",
        "",
        f"- Published-matrix audit universe: {summary['n_profiles']} profiles ({summary['n_csf']} CSF; {summary['n_plasma']} plasma).",
        f"- Traceable BrainTrace outputs: {summary['n_traceable_outputs']}/{summary['n_profiles']} ({summary['traceable_output_percent']:.1f}%).",
        f"- OMPFC Top1: CSF {int(csf['OMPFC_top1_numerator'])}/{int(csf['OMPFC_top1_denominator'])} ({csf['OMPFC_top1_percent']:.1f}%); plasma {int(plasma['OMPFC_top1_numerator'])}/{int(plasma['OMPFC_top1_denominator'])} ({plasma['OMPFC_top1_percent']:.1f}%).",
        "",
        "## Independent tumour-control diagnostics",
        "",
        f"Six two-sided Mann-Whitney U tests were run across fluid and metric; the smallest Benjamini-Hochberg FDR was {summary['minimum_bh_fdr']:.6f}. All tests are profile-level analyses with pairing unavailable; patient-level dependence cannot be verified from the public matrix.",
        "",
        "| Fluid | Metric | tumour n | control n | raw P | BH-FDR |",
        "|---|---|---:|---:|---:|---:|",
    ]
    for _, row in comparisons.iterrows():
        lines.append(f"| {row['fluid']} | {row['metric']} | {int(row['n_tumour'])} | {int(row['n_control'])} | {row['raw_p']:.6f} | {row['bh_fdr']:.6f} |")
    lines.extend(
        [
            "",
            "## Exploratory marker correlations",
            "",
            "Marker associations are descriptive, fluid-specific Spearman correlations with the OMPFC network score; they are not matched-biofluid comparisons; patient-level dependence cannot be assessed from the public matrix.",
            "",
            "| Fluid | Marker class | Marker | n | rho | raw P | BH-FDR |",
            "|---|---|---|---:|---:|---:|---:|",
        ]
    )
    for _, row in correlations.iterrows():
        rho = "NA" if pd.isna(row["spearman_rho"]) else f"{row['spearman_rho']:.4f}"
        raw_p = "NA" if pd.isna(row["raw_p"]) else f"{row['raw_p']:.6f}"
        fdr = "NA" if pd.isna(row.get("bh_fdr")) else f"{row['bh_fdr']:.6f}"
        lines.append(f"| {row['fluid']} | {row['marker_class']} | {row['marker']} | {int(row['n'])} | {rho} | {raw_p} | {fdr} |")
    lines.extend(
        [
            "",
            "## Canonical outputs",
            "",
            "The sample ledger, per-profile outputs, rankings, independent-cohort distributions, statistics, machine-readable summaries, and audit manifest in this directory are the canonical Huang 2025 remediation outputs. Pseudo-paired CSF-plasma and synthetic-mixture outputs are intentionally absent.",
            "",
        ]
    )
    path.write_text("\n".join(lines), encoding="utf-8", newline="\n")


def run(args: argparse.Namespace) -> int:
    if not args.input_csv.exists():
        raise FileNotFoundError(f"Published Huang expression matrix not found: {args.input_csv}")
    if not args.db_path.exists():
        raise FileNotFoundError(f"BrainTrace source-tracing database not found: {args.db_path}")
    outdir = args.outdir.resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    expression, metadata = read_expression(args.input_csv)
    sample_ids = metadata["sample_id"].astype(str).tolist()
    sample_rows: list[dict[str, Any]] = []
    network_rows: list[dict[str, Any]] = []
    region_rows: list[dict[str, Any]] = []

    for index, meta in enumerate(metadata.to_dict("records"), start=1):
        sample_id = str(meta["sample_id"])
        frame = pd.DataFrame(
            {"gene_symbol": expression.index.astype(str), "log_tpm": expression[sample_id].to_numpy(dtype=np.float32)}
        )
        network_out = trace_network_expression(
            frame, min_overlap_fraction=0.50, project_to_vsd=True, enable_pairwise_rescue=False
        )
        region_out = trace_bo2023_secondary_regions(
            frame, network_out, str(args.db_path), int(args.atlas_id), topk=30
        )
        network_ranked = network_out.get("results", [])
        region_ranked = region_out.get("results", [])
        if len(network_ranked) != 10 or not region_ranked:
            raise RuntimeError(f"Untraceable locked-route output for {sample_id}")

        network_top3 = top_ids(network_ranked, "network_id", 3)
        region_top3 = top_ids(region_ranked, "region_id", 3)
        scores_by_id = {str(row["network_id"]): float(row["score"]) for row in network_ranked}
        confidences = [float(row["confidence"]) for row in network_ranked]
        marker_expression = {
            f"expr_{marker}": float(expression.loc[marker, sample_id]) if marker in expression.index else np.nan
            for marker in PLATELET_MARKERS + EV_MARKERS
        }
        sample_row = {
            **meta,
            "route": ROUTE_NAME,
            "input_scale": "source_log2_rpm_plus1_times_ln2",
            "network_traceability": network_out["meta"].get("traceability"),
            "region_traceability": region_out["meta"].get("traceability"),
            "BrainTrace_output_available": bool(
                network_out["meta"].get("traceability") == "high" and region_out["meta"].get("traceability") == "high"
            ),
            "n_input_genes": int(expression.shape[0]),
            "n_projector_genes": network_out["meta"]["reference_projection"].get("n_projector_genes"),
            "n_projector_overlap_genes": network_out["meta"]["reference_projection"].get("n_input_projector_overlap_genes"),
            "n_network_model_genes": network_out["meta"].get("n_model_genes"),
            "n_network_overlap_genes": network_out["meta"].get("n_overlap_genes"),
            "network_overlap_fraction": network_out["meta"].get("overlap_fraction"),
            "atlas_fit_score": float(network_ranked[0]["score"]),
            "network_margin": float(network_ranked[0]["score"] - network_ranked[1]["score"]),
            "network_entropy": normalized_entropy(confidences),
            "network_top1": network_top3[0],
            "network_top2": network_top3[1],
            "network_top3": network_top3[2],
            "OMPFC_network_score": scores_by_id.get(OMPFC_NETWORK, np.nan),
            "region_top1": region_top3[0],
            "region_top2": region_top3[1],
            "region_top3": region_top3[2],
            "n_candidate_regions": region_out["meta"].get("n_candidate_regions"),
            "n_region_overlap_genes": region_out["meta"].get("n_overlap_genes"),
            **marker_expression,
        }
        sample_rows.append(sample_row)
        for rank, row in enumerate(network_ranked, start=1):
            network_rows.append(
                {"sample_id": sample_id, "disease_group": meta["disease_group"], "fluid": meta["fluid"], "rank": rank, **dict(row)}
            )
        for rank, row in enumerate(region_ranked, start=1):
            region_rows.append(
                {"sample_id": sample_id, "disease_group": meta["disease_group"], "fluid": meta["fluid"], "rank": rank, **dict(row)}
            )
        print(f"[{index:03d}/{len(sample_ids):03d}] {sample_id}: {network_top3[0]} -> {region_top3[0]}")

    sample_df = pd.DataFrame(sample_rows)
    network_df = pd.DataFrame(network_rows)
    region_df = pd.DataFrame(region_rows)
    ledger = build_ledger(metadata, sample_df)
    comparisons = compare_tumour_control(sample_df)
    correlations = marker_correlations(sample_df)
    top1_distribution = distribution_table(sample_df, 1)
    top3_distribution = distribution_table(sample_df, 3)
    fluid = fluid_summary(sample_df)
    validate_canonical_invariants(metadata, ledger, sample_df, comparisons)
    summary = make_summary(metadata, sample_df, fluid, comparisons, correlations)

    ledger.to_csv(outdir / "huang_2025_sample_ledger.csv", index=False, lineterminator="\n")
    sample_df.to_csv(outdir / "huang_2025_sample_outputs.csv", index=False, lineterminator="\n")
    network_df.to_csv(outdir / "huang_2025_network_rankings.csv", index=False, lineterminator="\n")
    region_df.to_csv(outdir / "huang_2025_exact_region_rankings.csv", index=False, lineterminator="\n")
    top1_distribution.to_csv(outdir / "huang_2025_network_top1_distribution.csv", index=False, lineterminator="\n")
    top3_distribution.to_csv(outdir / "huang_2025_network_top3_distribution.csv", index=False, lineterminator="\n")
    fluid.to_csv(outdir / "huang_2025_fluid_summary.csv", index=False, lineterminator="\n")
    correlations.to_csv(outdir / "huang_2025_marker_correlations.csv", index=False, lineterminator="\n")
    comparisons.to_csv(outdir / "huang_2025_tumour_control_comparisons.csv", index=False, lineterminator="\n")
    pd.DataFrame(
        [
            {
                "protocol_status": summary["protocol_status"],
                "source_doi": SOURCE_DOI,
                "n_profiles": summary["n_profiles"],
                "n_csf": summary["n_csf"],
                "n_plasma": summary["n_plasma"],
                "n_traceable_outputs": summary["n_traceable_outputs"],
                "traceable_output_percent": summary["traceable_output_percent"],
                "patient_id_metadata": summary["patient_id_metadata"],
                "patient_paired_analysis": summary["patient_paired_analysis"],
                "synthetic_matched_admixture": summary["synthetic_matched_admixture"],
                "minimum_bh_fdr": summary["minimum_bh_fdr"],
            }
        ]
    ).to_csv(
        outdir / "huang_2025_canonical_summary.csv",
        index=False,
        lineterminator="\n",
    )

    manifest = {
        "protocol_status": summary["protocol_status"],
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "source": {
            "citation": "Huang et al., 2025",
            "doi": SOURCE_DOI,
            "article_url": SOURCE_ARTICLE_URL,
            "published_matrix_cohort": {
                "profiles": EXPECTED_TOTAL,
                "csf": EXPECTED_CSF,
                "plasma": EXPECTED_PLASMA,
                "patients_reported_by_source": 85,
                "disease_groups_reported_by_source": {"glioma": 18, "meningioma": 46, "control": 21},
            },
            "source_clinical_qc_note": SOURCE_QC_NOTE,
            "matrix_scale_interpretation": "audited as log2(RPM+1); converted to ln(RPM+1) by multiplying values by ln(2) before BrainTrace inference",
        },
        "provenance_guardrails": {
            "patient_id_metadata": "not_available_in_public_expression_matrix",
            "patient_paired_analysis": "NOT_SUPPORTED",
            "synthetic_matched_admixture": "REMOVED_FROM_CANONICAL_ANALYSIS",
            "sample_suffix_interpreted_as_patient_id": False,
            "CSF_to_plasma_sample_name_substitution": "PROHIBITED",
        },
        "analysis_scope": {
            "audit_universe": "all 159 published-matrix profiles",
            "analysis_unit": "profile-level observation within fluid-specific cohort; patient-level dependence cannot be assessed from the public matrix",
            "supported_claim": "technical portability and domain-shift audit",
            "not_supported_claims": summary["interpretation"]["not_supported_claims"],
        },
        "route": ROUTE_NAME,
        "pairwise_rescue": False,
        "model_assets": [asset_record(DEFAULT_BO2023_NETWORK_MODEL, "locked_network_model"), asset_record(DEFAULT_BO2023_REFERENCE_PROJECTOR, "locked_reference_projector")],
        "input_assets": [
            asset_record(args.input_csv, "published_expression_matrix"),
            asset_record(args.source_xlsb, "source_supplementary_xlsb"),
            asset_record(args.db_path, "braintrace_source_tracing_database"),
        ],
        "model_overlap": {
            "n_source_matrix_features": 83929,
            "n_selected_input_genes": int(sample_df["n_input_genes"].iloc[0]),
            "n_projector_genes": int(sample_df["n_projector_genes"].iloc[0]),
            "n_projector_overlap_genes": int(sample_df["n_projector_overlap_genes"].iloc[0]),
            "n_network_model_genes": int(sample_df["n_network_model_genes"].iloc[0]),
            "n_network_overlap_genes": int(sample_df["n_network_overlap_genes"].iloc[0]),
        },
        "outputs": [
            "huang_2025_sample_ledger.csv",
            "huang_2025_sample_outputs.csv",
            "huang_2025_network_rankings.csv",
            "huang_2025_exact_region_rankings.csv",
            "huang_2025_network_top1_distribution.csv",
            "huang_2025_network_top3_distribution.csv",
            "huang_2025_fluid_summary.csv",
            "huang_2025_marker_correlations.csv",
            "huang_2025_tumour_control_comparisons.csv",
            "huang_2025_canonical_summary.csv",
            "huang_2025_canonical_summary.json",
            "HUANG_2025_RESULTS.md",
        ],
    }
    write_json(outdir / "huang_2025_audit_manifest.json", manifest)
    write_json(outdir / "huang_2025_canonical_summary.json", summary)
    write_results_markdown(outdir / "HUANG_2025_RESULTS.md", summary, fluid, comparisons, correlations)
    print(json.dumps({key: summary[key] for key in ("n_profiles", "n_csf", "n_plasma", "n_traceable_outputs", "minimum_bh_fdr")}, indent=2))
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-csv", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--source-xlsb", type=Path, default=None)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument("--db-path", type=Path, default=DEFAULT_DB_PATH)
    parser.add_argument("--atlas-id", type=int, default=1)
    return run(parser.parse_args())


if __name__ == "__main__":
    raise SystemExit(main())
