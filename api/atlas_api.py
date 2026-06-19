from __future__ import annotations

import os
import sqlite3
from pathlib import Path
from typing import Any

from fastapi import Depends, FastAPI, Header, HTTPException, Query


DB_PATH = Path(os.environ.get("CFRNA_DB_PATH", "cfrna_source_tracing.db")).resolve()
API_KEY = os.environ.get("CFRNA_API_KEY", "").strip()

app = FastAPI(title="cfRNA-BrainTrace Atlas API", version="1.0.0")


def require_api_key(x_api_key: str | None = Header(default=None)) -> None:
    if API_KEY and x_api_key != API_KEY:
        raise HTTPException(status_code=401, detail="Invalid API key")


def query_rows(sql: str, params: tuple[Any, ...] = ()) -> list[dict[str, Any]]:
    if not DB_PATH.exists():
        raise HTTPException(status_code=503, detail=f"Database not found: {DB_PATH}")
    conn = sqlite3.connect(DB_PATH, timeout=15)
    conn.row_factory = sqlite3.Row
    try:
        rows = conn.execute(sql, params).fetchall()
        return [dict(row) for row in rows]
    finally:
        conn.close()


def expression_filter(region_id: list[str], celltype: str | None) -> tuple[str, list[Any]]:
    clauses: list[str] = []
    params: list[Any] = []
    clean_regions = [str(item).strip() for item in region_id if str(item).strip()]
    if clean_regions:
        placeholders = ",".join(["?"] * len(clean_regions))
        clauses.append(f"region_id IN ({placeholders})")
        params.extend(clean_regions)
    if celltype:
        clauses.append("cell_type_marker = ?")
        params.append(celltype)
    return (" AND " + " AND ".join(clauses) if clauses else ""), params


@app.get("/health", dependencies=[Depends(require_api_key)])
def health() -> dict[str, Any]:
    exists = DB_PATH.exists()
    return {"ok": exists, "db_path": str(DB_PATH), "db_exists": exists}


@app.get("/atlas_versions", dependencies=[Depends(require_api_key)])
def atlas_versions() -> dict[str, Any]:
    items = query_rows(
        """
        SELECT atlas_id, atlas_name, species, level, build_version, gene_id_type,
               normalization, created_at, notes
        FROM atlas_versions
        ORDER BY atlas_id DESC
        """
    )
    return {"items": items}


@app.get("/atlases/{atlas_id}/regions", dependencies=[Depends(require_api_key)])
def atlas_regions(atlas_id: int) -> dict[str, Any]:
    items = query_rows(
        """
        SELECT
            region_id,
            MAX(region_name) AS region_name,
            COUNT(DISTINCT gene_symbol) AS n_genes,
            MAX(sample_count) AS sample_count
        FROM reference_expression
        WHERE atlas_id = ?
        GROUP BY region_id
        ORDER BY region_id
        """,
        (atlas_id,),
    )
    return {"items": items}


@app.get("/atlases/{atlas_id}/celltypes", dependencies=[Depends(require_api_key)])
def atlas_celltypes(atlas_id: int) -> dict[str, Any]:
    items = query_rows(
        """
        SELECT cell_type_marker, COUNT(DISTINCT gene_symbol) AS n_genes
        FROM reference_expression
        WHERE atlas_id = ?
          AND cell_type_marker IS NOT NULL
          AND cell_type_marker != ''
        GROUP BY cell_type_marker
        ORDER BY n_genes DESC, cell_type_marker
        LIMIT 200
        """,
        (atlas_id,),
    )
    return {"items": items}


@app.get("/atlases/{atlas_id}/region-ranking", dependencies=[Depends(require_api_key)])
def atlas_region_ranking(
    atlas_id: int,
    region_id: list[str] = Query(default=[]),
    celltype: str = "",
) -> dict[str, Any]:
    extra_sql, extra_params = expression_filter(region_id, celltype or None)
    items = query_rows(
        f"""
        SELECT
            region_id,
            MAX(region_name) AS region_name,
            COUNT(DISTINCT gene_symbol) AS gene_count,
            AVG(avg_tpm) AS mean_tpm,
            MAX(sample_count) AS sample_count,
            SUM(CASE WHEN expression_class IS NOT NULL AND expression_class != '' THEN 1 ELSE 0 END) AS classified_genes
        FROM reference_expression
        WHERE atlas_id = ?
        {extra_sql}
        GROUP BY region_id
        ORDER BY mean_tpm DESC, gene_count DESC
        LIMIT 80
        """,
        tuple([atlas_id] + extra_params),
    )
    return {"items": items}


@app.get("/atlases/{atlas_id}/gene-candidates", dependencies=[Depends(require_api_key)])
def atlas_gene_candidates(
    atlas_id: int,
    region_id: list[str] = Query(default=[]),
    celltype: str = "",
    limit: int = Query(default=40, ge=1, le=100),
) -> dict[str, Any]:
    extra_sql, extra_params = expression_filter(region_id, celltype or None)
    rows = query_rows(
        f"""
        SELECT gene_symbol, AVG(avg_tpm) AS mean_tpm
        FROM reference_expression
        WHERE atlas_id = ?
        {extra_sql}
        GROUP BY gene_symbol
        HAVING mean_tpm IS NOT NULL
        ORDER BY mean_tpm DESC
        LIMIT ?
        """,
        tuple([atlas_id] + extra_params + [limit]),
    )
    return {"genes": [row["gene_symbol"] for row in rows]}


@app.get("/atlases/{atlas_id}/expression", dependencies=[Depends(require_api_key)])
def atlas_expression(
    atlas_id: int,
    region_id: list[str] = Query(default=[]),
    gene: list[str] = Query(default=[]),
    celltype: str = "",
) -> dict[str, Any]:
    extra_sql, extra_params = expression_filter(region_id, celltype or None)
    genes = [item.strip().upper() for item in gene if item.strip()]
    gene_sql = ""
    gene_params: list[Any] = []
    if genes:
        placeholders = ",".join(["?"] * len(genes))
        gene_sql = f" AND upper(gene_symbol) IN ({placeholders})"
        gene_params = genes
    items = query_rows(
        f"""
        SELECT gene_symbol, region_id, region_name, avg_tpm, median_tpm,
               sample_count, expression_class, cell_type_marker
        FROM reference_expression
        WHERE atlas_id = ?
        {extra_sql}
        {gene_sql}
        ORDER BY region_id, gene_symbol
        LIMIT 5000
        """,
        tuple([atlas_id] + extra_params + gene_params),
    )
    return {"items": items}


@app.get("/atlases/{atlas_id}/marker-evidence", dependencies=[Depends(require_api_key)])
def marker_evidence(
    atlas_id: int,
    region_id: list[str] = Query(default=[]),
    celltype: str = "",
    limit: int = Query(default=200, ge=1, le=1000),
) -> dict[str, Any]:
    sigsets = query_rows(
        """
        SELECT sigset_id, method, topk_per_region, created_at
        FROM signature_sets
        WHERE atlas_id = ?
        ORDER BY sigset_id DESC
        LIMIT 1
        """,
        (atlas_id,),
    )
    if sigsets:
        sigset_id = int(sigsets[0]["sigset_id"])
        region_clause = ""
        params: list[Any] = [atlas_id, sigset_id]
        clean_regions = [str(item).strip() for item in region_id if str(item).strip()]
        if clean_regions:
            placeholders = ",".join(["?"] * len(clean_regions))
            region_clause = f" AND sg.region_id IN ({placeholders})"
            params.extend(clean_regions)
        celltype_clause = ""
        if celltype:
            celltype_clause = " AND r.cell_type_marker = ?"
            params.append(celltype)
        params.append(limit)
        items = query_rows(
            f"""
            SELECT
                sg.region_id,
                MAX(r.region_name) AS region_name,
                sg.gene_symbol,
                sg.weight,
                MAX(r.avg_tpm) AS avg_tpm,
                MAX(r.expression_class) AS expression_class,
                MAX(r.cell_type_marker) AS cell_type_marker
            FROM signature_genes sg
            LEFT JOIN reference_expression r
              ON r.atlas_id = ?
             AND r.region_id = sg.region_id
             AND r.gene_symbol = sg.gene_symbol
            WHERE sg.sigset_id = ?
            {region_clause}
            {celltype_clause}
            GROUP BY sg.region_id, sg.gene_symbol, sg.weight
            ORDER BY sg.region_id, sg.weight DESC, avg_tpm DESC
            LIMIT ?
            """,
            tuple(params),
        )
        return {"items": items}

    extra_sql, extra_params = expression_filter(region_id, celltype or None)
    items = query_rows(
        f"""
        SELECT region_id, region_name, gene_symbol, 1.0 AS weight, avg_tpm,
               expression_class, cell_type_marker
        FROM reference_expression
        WHERE atlas_id = ?
        {extra_sql}
        ORDER BY avg_tpm DESC
        LIMIT ?
        """,
        tuple([atlas_id] + extra_params + [limit]),
    )
    return {"items": items}
