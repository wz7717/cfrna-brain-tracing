"""Public API for reporting and export helpers."""

from __future__ import annotations

import json
import os
import sqlite3
import tempfile
import zipfile
from typing import Optional

import pandas as pd

from core.methods import method_label


def _json_safe_records(df: pd.DataFrame):
    return json.loads(df.to_json(orient="records", force_ascii=False))


def build_run_summary(db_path: str, run_id: str) -> dict:
    db_path = os.path.abspath(db_path)
    conn = sqlite3.connect(db_path)
    try:
        run = pd.read_sql_query("SELECT * FROM analysis_runs WHERE run_id=?", conn, params=[run_id])
        res = pd.read_sql_query(
            "SELECT * FROM analysis_results WHERE run_id=? ORDER BY rank ASC",
            conn,
            params=[run_id],
        )
    finally:
        conn.close()
    if run.empty:
        raise ValueError(f"run_id not found: {run_id}")
    meta = run.iloc[0].to_dict()
    try:
        params = json.loads(meta.get("params_json") or "{}")
    except Exception:
        params = {}
    top_hit = res.iloc[0].to_dict() if not res.empty else {}
    return {
        "run_id": run_id,
        "sample_id": meta.get("sample_id"),
        "method": meta.get("method"),
        "method_label": method_label(meta.get("method")),
        "atlas_id": meta.get("atlas_id"),
        "sigset_id": meta.get("sigset_id"),
        "created_at": meta.get("created_at"),
        "parameter_snapshot": params,
        "atlas_name": params.get("atlas_name"),
        "atlas_normalization": params.get("atlas_normalization"),
        "vsd_compatible_mode": params.get("vsd_compatible_mode", False),
        "result_interpretation": params.get("result_interpretation"),
        "recommended_interpretation": params.get("recommended_interpretation"),
        "top_region": top_hit.get("region_id"),
        "top_score": top_hit.get("score"),
        "top_confidence": top_hit.get("confidence"),
        "n_results": int(len(res)),
        "results": _json_safe_records(res),
    }


def export_run_bundle(db_path: str, run_id: str, out_zip_path: Optional[str] = None) -> str:
    db_path = os.path.abspath(db_path)
    conn = sqlite3.connect(db_path)
    try:
        run = pd.read_sql_query(
            "SELECT * FROM analysis_runs WHERE run_id=?",
            conn,
            params=[run_id],
        )
        res = pd.read_sql_query(
            "SELECT * FROM analysis_results WHERE run_id=? ORDER BY rank ASC",
            conn,
            params=[run_id],
        )
        if run.empty or res.empty:
            raise ValueError(f"run_id not found or has no results: {run_id}")
        meta = run.iloc[0].to_dict()
        try:
            meta["params"] = json.loads(meta.get("params_json") or "{}")
        except Exception:
            meta["params"] = meta.get("params_json")
        meta["method_label"] = method_label(meta.get("method"))
        summary = build_run_summary(db_path, run_id)
        tmpdir = tempfile.mkdtemp(prefix="cfrna_tracing_")
        results_tsv = os.path.join(tmpdir, "results.tsv")
        meta_json = os.path.join(tmpdir, "meta.json")
        summary_json = os.path.join(tmpdir, "summary.json")
        raw_json = os.path.join(tmpdir, "raw_results.json")
        res.to_csv(results_tsv, sep="\t", index=False)
        with open(meta_json, "w", encoding="utf-8") as f:
            json.dump(meta, f, ensure_ascii=False, indent=2)
        with open(summary_json, "w", encoding="utf-8") as f:
            json.dump(summary, f, ensure_ascii=False, indent=2)
        with open(raw_json, "w", encoding="utf-8") as f:
            json.dump(res.to_dict(orient="records"), f, ensure_ascii=False, indent=2)
        if out_zip_path is None:
            out_zip_path = os.path.join(tmpdir, f"run_{run_id}.zip")
        out_zip_path = os.path.abspath(out_zip_path)
        os.makedirs(os.path.dirname(out_zip_path), exist_ok=True)
        with zipfile.ZipFile(out_zip_path, "w", compression=zipfile.ZIP_DEFLATED) as z:
            z.write(results_tsv, arcname="results.tsv")
            z.write(meta_json, arcname="meta.json")
            z.write(summary_json, arcname="summary.json")
            z.write(raw_json, arcname="raw_results.json")
        return out_zip_path
    finally:
        conn.close()


build_report_bundle = export_run_bundle


def export_benchmark_paper_figures(*args, **kwargs):
    from reporting.benchmark_figure_export import export_benchmark_paper_figures as _fn
    return _fn(*args, **kwargs)


def build_benchmark_figure_bundle_bytes(*args, **kwargs):
    from reporting.benchmark_figure_export import build_benchmark_figure_bundle_bytes as _fn
    return _fn(*args, **kwargs)


def export_benchmark_report_pdf(*args, **kwargs):
    from reporting.benchmark_figure_export import export_benchmark_report_pdf as _fn
    return _fn(*args, **kwargs)


def build_benchmark_report_bundle_bytes(*args, **kwargs):
    from reporting.benchmark_figure_export import build_benchmark_report_bundle_bytes as _fn
    return _fn(*args, **kwargs)


__all__ = [
    "export_run_bundle",
    "build_run_summary",
    "build_report_bundle",
    "export_benchmark_paper_figures",
    "build_benchmark_figure_bundle_bytes",
    "export_benchmark_report_pdf",
    "build_benchmark_report_bundle_bytes",
]
