from __future__ import annotations

import json
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd


@dataclass(frozen=True)
class ProjectorFit:
    genes: np.ndarray
    slope: np.ndarray
    intercept: np.ndarray
    clip_low: np.ndarray
    clip_high: np.ndarray
    logcpm_mean: np.ndarray
    logcpm_sd: np.ndarray
    vsd_mean: np.ndarray
    vsd_sd: np.ndarray
    fallback_reason: np.ndarray


def _json_default(value: Any) -> Any:
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, np.floating):
        return float(value)
    if isinstance(value, np.ndarray):
        return value.tolist()
    return str(value)


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(payload, indent=2, ensure_ascii=False, default=_json_default) + "\n",
        encoding="utf-8",
    )


def read_bo2023_gene_matrix(path: Path, dtype: str | type = "float32") -> pd.DataFrame:
    """Read Bo2023 gene x sample matrix whose header omits the gene-id column."""
    with path.open("rt", encoding="utf-8") as handle:
        samples = handle.readline().rstrip("\n\r").split("\t")
    names = ["gene_id", *samples]
    matrix = pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=names,
        skiprows=1,
        index_col=0,
        low_memory=False,
    )
    matrix.index = matrix.index.astype(str).str.strip()
    matrix.columns = matrix.columns.astype(str).str.strip()
    matrix = matrix[~matrix.index.duplicated(keep="first")]
    return matrix.apply(pd.to_numeric, errors="coerce").fillna(0.0).astype(dtype)


def read_gene_map(path: Path) -> pd.DataFrame:
    mapping = pd.read_csv(path, usecols=["Gene.stable.ID", "Gene.name"])
    mapping = mapping.rename(columns={"Gene.stable.ID": "gene_id", "Gene.name": "gene_symbol"})
    mapping["gene_id"] = mapping["gene_id"].astype(str).str.strip()
    mapping["gene_symbol"] = mapping["gene_symbol"].map(clean_excel_date_gene_symbol)
    mapping.loc[mapping["gene_symbol"].isin(["", "nan", "None"]), "gene_symbol"] = pd.NA
    return mapping.dropna(subset=["gene_id"]).drop_duplicates(["gene_id", "gene_symbol"])


def map_index_to_symbols(matrix: pd.DataFrame, gene_map: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    symbol_by_id = (
        gene_map.dropna(subset=["gene_symbol"])
        .drop_duplicates("gene_id")
        .set_index("gene_id")["gene_symbol"]
    )
    raw = matrix.index.to_series()
    symbols = raw.map(symbol_by_id).replace("", pd.NA).fillna(raw)
    mapped = matrix.groupby(symbols, sort=True).mean()
    mapped.index.name = "gene_symbol"
    audit = pd.DataFrame({"gene_id": raw.to_numpy(), "gene_symbol": symbols.to_numpy()})
    audit["mapped_to_symbol"] = audit["gene_id"].ne(audit["gene_symbol"])
    return mapped, audit


def compute_logcpm(counts: pd.DataFrame) -> pd.DataFrame:
    numeric = counts.apply(pd.to_numeric, errors="coerce").fillna(0.0).astype("float64")
    lib_size = numeric.sum(axis=0)
    safe_lib_size = lib_size.mask(lib_size <= 0, np.nan)
    cpm = numeric.divide(safe_lib_size, axis=1) * 1_000_000.0
    return np.log1p(cpm).fillna(0.0).astype("float32")


def align_matrices(left: pd.DataFrame, right: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, list[str], list[str]]:
    common_genes = sorted(set(left.index.astype(str)) & set(right.index.astype(str)))
    common_samples = [sample for sample in left.columns.astype(str) if sample in set(right.columns.astype(str))]
    return (
        left.loc[common_genes, common_samples],
        right.loc[common_genes, common_samples],
        common_genes,
        common_samples,
    )


def fit_linear_projector(
    logcpm: pd.DataFrame,
    vsd: pd.DataFrame,
    min_nonzero_samples: int = 10,
    min_logcpm_sd: float = 1e-8,
    clip_quantiles: tuple[float, float] = (0.005, 0.995),
) -> tuple[ProjectorFit, pd.DataFrame]:
    if not logcpm.index.equals(vsd.index) or not logcpm.columns.equals(vsd.columns):
        raise ValueError("logCPM and VSD matrices must be identically aligned before fitting.")

    x = logcpm.to_numpy(dtype=np.float64)
    y = vsd.to_numpy(dtype=np.float64)
    genes = logcpm.index.astype(str).to_numpy()
    n = x.shape[1]

    x_mean = x.mean(axis=1)
    y_mean = y.mean(axis=1)
    x_centered = x - x_mean[:, None]
    y_centered = y - y_mean[:, None]
    x_var_sum = np.square(x_centered).sum(axis=1)
    y_var_sum = np.square(y_centered).sum(axis=1)
    xy_sum = (x_centered * y_centered).sum(axis=1)

    x_sd = x.std(axis=1, ddof=1)
    y_sd = y.std(axis=1, ddof=1)
    nonzero = (x > 0).sum(axis=1)

    slope = np.divide(xy_sum, x_var_sum, out=np.zeros_like(xy_sum), where=x_var_sum > 0)
    intercept = y_mean - slope * x_mean
    pred = slope[:, None] * x + intercept[:, None]

    residual = y - pred
    residual_sd = np.sqrt(np.square(residual).sum(axis=1) / np.maximum(n - 2, 1))
    r = np.divide(
        xy_sum,
        np.sqrt(x_var_sum * y_var_sum),
        out=np.zeros_like(xy_sum),
        where=(x_var_sum > 0) & (y_var_sum > 0),
    )
    r2 = np.square(r)
    spearman = [
        float(pd.Series(x[i, :]).corr(pd.Series(y[i, :]), method="spearman"))
        if x_sd[i] > 0 and y_sd[i] > 0
        else float("nan")
        for i in range(x.shape[0])
    ]

    clip_low = np.quantile(y, clip_quantiles[0], axis=1)
    clip_high = np.quantile(y, clip_quantiles[1], axis=1)
    fallback = np.full(x.shape[0], "", dtype=object)
    low_nonzero = nonzero < int(min_nonzero_samples)
    low_sd = x_sd <= float(min_logcpm_sd)
    fallback[low_nonzero] = "low_nonzero_count_samples"
    fallback[low_sd] = np.where(fallback[low_sd] == "", "low_logcpm_sd", fallback[low_sd] + ";low_logcpm_sd")
    if np.any(fallback != ""):
        slope[fallback != ""] = 0.0
        intercept[fallback != ""] = y_mean[fallback != ""]
        pred[fallback != "", :] = y_mean[fallback != "", None]
        residual = y - pred
        residual_sd = np.sqrt(np.square(residual).sum(axis=1) / np.maximum(n - 1, 1))

    fit = ProjectorFit(
        genes=genes,
        slope=slope.astype("float32"),
        intercept=intercept.astype("float32"),
        clip_low=clip_low.astype("float32"),
        clip_high=clip_high.astype("float32"),
        logcpm_mean=x_mean.astype("float32"),
        logcpm_sd=x_sd.astype("float32"),
        vsd_mean=y_mean.astype("float32"),
        vsd_sd=y_sd.astype("float32"),
        fallback_reason=fallback.astype(str),
    )
    params = pd.DataFrame(
        {
            "gene_symbol": genes,
            "n_train_samples": n,
            "n_nonzero_count_samples": nonzero,
            "logcpm_mean": x_mean,
            "logcpm_sd": x_sd,
            "vsd_mean": y_mean,
            "vsd_sd": y_sd,
            "slope": slope,
            "intercept": intercept,
            "r2": r2,
            "spearman_r": spearman,
            "residual_sd": residual_sd,
            "fallback_reason": fallback,
            "clip_low": clip_low,
            "clip_high": clip_high,
        }
    )
    return fit, params


def apply_projector(fit: ProjectorFit, logcpm: pd.DataFrame) -> pd.DataFrame:
    aligned = logcpm.reindex(fit.genes).fillna(0.0).astype("float64")
    x = aligned.to_numpy(dtype=np.float64)
    projected = fit.slope[:, None] * x + fit.intercept[:, None]
    projected = np.clip(projected, fit.clip_low[:, None], fit.clip_high[:, None])
    return pd.DataFrame(projected.astype("float32"), index=fit.genes, columns=aligned.columns)


def save_projector_npz(path: Path, fit: ProjectorFit, metadata: dict[str, Any] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        path,
        genes=fit.genes.astype(str),
        slope=fit.slope,
        intercept=fit.intercept,
        clip_low=fit.clip_low,
        clip_high=fit.clip_high,
        logcpm_mean=fit.logcpm_mean,
        logcpm_sd=fit.logcpm_sd,
        vsd_mean=fit.vsd_mean,
        vsd_sd=fit.vsd_sd,
        fallback_reason=fit.fallback_reason.astype(str),
        metadata=json.dumps(metadata or {}, ensure_ascii=False),
    )


def load_projector_npz(path: Path) -> ProjectorFit:
    with np.load(path, allow_pickle=False) as data:
        return ProjectorFit(
            genes=data["genes"].astype(str),
            slope=data["slope"].astype("float32"),
            intercept=data["intercept"].astype("float32"),
            clip_low=data["clip_low"].astype("float32"),
            clip_high=data["clip_high"].astype("float32"),
            logcpm_mean=data["logcpm_mean"].astype("float32"),
            logcpm_sd=data["logcpm_sd"].astype("float32"),
            vsd_mean=data["vsd_mean"].astype("float32"),
            vsd_sd=data["vsd_sd"].astype("float32"),
            fallback_reason=data["fallback_reason"].astype(str),
        )


def summarize_fit(projected: pd.DataFrame, native: pd.DataFrame, params: pd.DataFrame) -> dict[str, Any]:
    x = projected.to_numpy(dtype=np.float64)
    y = native.to_numpy(dtype=np.float64)
    sample_corr = [
        float(pd.Series(x[:, i]).corr(pd.Series(y[:, i]), method="pearson"))
        for i in range(x.shape[1])
    ]
    gene_corr = [
        float(pd.Series(x[i, :]).corr(pd.Series(y[i, :]), method="pearson"))
        for i in range(x.shape[0])
    ]
    delta = x - y
    return {
        "n_genes": int(projected.shape[0]),
        "n_samples": int(projected.shape[1]),
        "global_pearson": float(pd.Series(x.ravel()).corr(pd.Series(y.ravel()), method="pearson")),
        "mae": float(np.mean(np.abs(delta))),
        "rmse": float(np.sqrt(np.mean(np.square(delta)))),
        "median_sample_pearson": float(np.nanmedian(sample_corr)),
        "p10_sample_pearson": float(np.nanquantile(sample_corr, 0.10)),
        "median_gene_pearson": float(np.nanmedian(gene_corr)),
        "p10_gene_pearson": float(np.nanquantile(gene_corr, 0.10)),
        "median_gene_r2": float(params["r2"].median()),
        "p10_gene_r2": float(params["r2"].quantile(0.10)),
        "n_fallback_genes": int(params["fallback_reason"].astype(str).ne("").sum()),
    }


def missing_items(expected: Iterable[str], observed: Iterable[str]) -> list[str]:
    seen = set(map(str, observed))
    return sorted(str(x) for x in expected if str(x) not in seen)
EXCEL_DATE_SYMBOL_RE = re.compile(r"^\d{4}-(\d{2})-(\d{2})(?: 00:00:00)?$")


def clean_excel_date_gene_symbol(value: object) -> str:
    text = "" if pd.isna(value) else str(value).strip()
    match = EXCEL_DATE_SYMBOL_RE.match(text)
    if not match:
        return text
    month = int(match.group(1))
    day = int(match.group(2))
    prefix = {3: "MARCH", 9: "SEPT", 12: "DEC"}.get(month)
    if not prefix:
        return text
    return f"{prefix}{day}"
