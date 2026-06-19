from __future__ import annotations

import importlib.util
import sqlite3
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List


REQUIRED_TABLES = [
    "cfrna_samples",
    "cfrna_expression",
    "macaque_brain_atlas",
    "reference_expression",
    "source_tracing_results",
    "atlas_versions",
    "signature_sets",
    "signature_genes",
    "analysis_runs",
    "analysis_results",
    "sample_qc",
]

BO2023_CORE_TABLES = [
    "bo2023_buildkit_catalog",
    "bo2023_regions",
    "bo2023_region_qc_overview",
    "bo2023_region_sample_coverage",
]

REQUIRED_MODULES = [
    "streamlit",
    "pandas",
    "numpy",
    "plotly",
    "sklearn",
    "scipy",
    "openpyxl",
    "matplotlib",
]


@dataclass
class CheckItem:
    name: str
    status: str
    detail: str


def _object_exists(conn: sqlite3.Connection, name: str) -> bool:
    cur = conn.cursor()
    cur.execute(
        "SELECT name FROM sqlite_master WHERE name=? AND type IN ('table', 'view')",
        (name,),
    )
    return cur.fetchone() is not None


def _count_rows(conn: sqlite3.Connection, name: str) -> int | None:
    try:
        cur = conn.cursor()
        cur.execute(f"SELECT COUNT(*) FROM {name}")
        return int(cur.fetchone()[0])
    except Exception:
        return None


def _check_modules(modules: Iterable[str]) -> CheckItem:
    missing = [m for m in modules if importlib.util.find_spec(m) is None]
    if missing:
        return CheckItem("Python dependencies", "error", "Missing modules: " + ", ".join(missing))
    return CheckItem("Python dependencies", "ok", "Required runtime modules are importable.")


def run_system_check(db_path: str, project_root: str | Path | None = None) -> Dict[str, object]:
    root = Path(project_root) if project_root is not None else Path.cwd()
    db = Path(db_path)
    if not db.is_absolute():
        db = root / db

    items: List[CheckItem] = [_check_modules(REQUIRED_MODULES)]
    if not db.exists():
        items.append(CheckItem("SQLite database", "warning", f"Database file not found yet: {db}"))
        return {
            "ok": False,
            "has_database": False,
            "items": items,
            "missing_required_tables": REQUIRED_TABLES,
            "missing_bo2023_objects": BO2023_CORE_TABLES,
            "bo2023_expression_ready": False,
        }

    items.append(CheckItem("SQLite database", "ok", f"Database file found: {db}"))
    missing_required: List[str] = []
    missing_bo2023: List[str] = []
    counts: Dict[str, int | None] = {}
    bo2023_expression_rows = 0
    bo2023_expression_atlases = 0

    try:
        conn = sqlite3.connect(str(db))
        try:
            for table in REQUIRED_TABLES:
                if not _object_exists(conn, table):
                    missing_required.append(table)
                else:
                    counts[table] = _count_rows(conn, table)
            for obj in BO2023_CORE_TABLES:
                if not _object_exists(conn, obj):
                    missing_bo2023.append(obj)
                else:
                    counts[obj] = _count_rows(conn, obj)
            if _object_exists(conn, "atlas_versions") and _object_exists(conn, "reference_expression"):
                cur = conn.cursor()
                cur.execute(
                    """
                    SELECT COUNT(*), COUNT(DISTINCT r.atlas_id)
                    FROM reference_expression r
                    JOIN atlas_versions a ON a.atlas_id = r.atlas_id
                    WHERE lower(COALESCE(a.atlas_name, '') || ' ' || COALESCE(a.build_version, '')) LIKE '%bo2023%'
                       OR lower(COALESCE(a.atlas_name, '') || ' ' || COALESCE(a.build_version, '')) LIKE '%wang%'
                    """
                )
                row = cur.fetchone()
                bo2023_expression_rows = int(row[0] or 0)
                bo2023_expression_atlases = int(row[1] or 0)
        finally:
            conn.close()
    except Exception as exc:
        items.append(CheckItem("SQLite connection", "error", f"Could not inspect database: {exc}"))
        return {
            "ok": False,
            "has_database": True,
            "items": items,
            "missing_required_tables": REQUIRED_TABLES,
            "missing_bo2023_objects": BO2023_CORE_TABLES,
            "bo2023_expression_ready": False,
        }

    if missing_required:
        items.append(CheckItem("Required SQLite tables", "error", "Missing: " + ", ".join(missing_required)))
    else:
        items.append(CheckItem("Required SQLite tables", "ok", "All required application tables are present."))

    if missing_bo2023:
        items.append(CheckItem("Bo2023 annotation layer", "warning", "Missing optional objects: " + ", ".join(missing_bo2023)))
    else:
        items.append(CheckItem("Bo2023 annotation layer", "ok", "Bo2023 annotation/QC browser objects are present."))

    if bo2023_expression_rows > 0:
        items.append(
            CheckItem(
                "Bo2023 expression matrix",
                "ok",
                f"Expression-ready atlas rows: {bo2023_expression_rows:,} across {bo2023_expression_atlases} atlas version(s).",
            )
        )
    else:
        items.append(
            CheckItem(
                "Bo2023 expression matrix",
                "warning",
                "No imported Bo2023 gene x region expression matrix was detected in reference_expression.",
            )
        )

    ok = not any(i.status == "error" for i in items)
    return {
        "ok": ok,
        "has_database": True,
        "items": items,
        "counts": counts,
        "missing_required_tables": missing_required,
        "missing_bo2023_objects": missing_bo2023,
        "bo2023_expression_ready": bo2023_expression_rows > 0,
        "bo2023_expression_rows": bo2023_expression_rows,
    }
