from __future__ import annotations

import json
import os
import urllib.parse
import urllib.request
from typing import Iterable

import pandas as pd
import streamlit as st


def _config(name: str) -> str | None:
    value = os.environ.get(name)
    if value:
        return value.strip()
    try:
        secret_value = st.secrets.get(name)
    except Exception:
        secret_value = None
    return str(secret_value).strip() if secret_value else None


def is_configured() -> bool:
    return bool(_config("CFRNA_API_URL"))


def _url(path: str, params: dict | None = None) -> str:
    base = (_config("CFRNA_API_URL") or "").rstrip("/")
    query = urllib.parse.urlencode(params or {}, doseq=True)
    return f"{base}{path}" + (f"?{query}" if query else "")


def _get_json(path: str, params: dict | None = None) -> dict:
    headers = {"Accept": "application/json"}
    api_key = _config("CFRNA_API_KEY")
    if api_key:
        headers["X-API-Key"] = api_key
    req = urllib.request.Request(_url(path, params), headers=headers, method="GET")
    with urllib.request.urlopen(req, timeout=60) as response:
        return json.loads(response.read().decode("utf-8"))


def _records_to_frame(payload: dict, key: str = "items") -> pd.DataFrame:
    return pd.DataFrame(payload.get(key, []))


@st.cache_data(show_spinner=False, ttl=300)
def atlas_catalog() -> pd.DataFrame:
    if not is_configured():
        return pd.DataFrame()
    return _records_to_frame(_get_json("/atlas_versions"))


@st.cache_data(show_spinner=False, ttl=300)
def atlas_regions(atlas_id: int) -> pd.DataFrame:
    if not is_configured():
        return pd.DataFrame()
    return _records_to_frame(_get_json(f"/atlases/{int(atlas_id)}/regions"))


@st.cache_data(show_spinner=False, ttl=300)
def atlas_celltypes(atlas_id: int) -> pd.DataFrame:
    if not is_configured():
        return pd.DataFrame()
    return _records_to_frame(_get_json(f"/atlases/{int(atlas_id)}/celltypes"))


@st.cache_data(show_spinner=False, ttl=300)
def atlas_region_ranking(atlas_id: int, region_ids: tuple[str, ...], celltype: str | None) -> pd.DataFrame:
    if not is_configured():
        return pd.DataFrame()
    return _records_to_frame(
        _get_json(
            f"/atlases/{int(atlas_id)}/region-ranking",
            {"region_id": list(region_ids), "celltype": celltype or ""},
        )
    )


@st.cache_data(show_spinner=False, ttl=300)
def atlas_gene_candidates(atlas_id: int, region_ids: tuple[str, ...], celltype: str | None, limit: int) -> list[str]:
    if not is_configured():
        return []
    payload = _get_json(
        f"/atlases/{int(atlas_id)}/gene-candidates",
        {"region_id": list(region_ids), "celltype": celltype or "", "limit": int(limit)},
    )
    return [str(item) for item in payload.get("genes", [])]


@st.cache_data(show_spinner=False, ttl=300)
def atlas_expression_matrix(
    atlas_id: int,
    region_ids: tuple[str, ...],
    gene_symbols: tuple[str, ...],
    celltype: str | None,
) -> pd.DataFrame:
    if not is_configured():
        return pd.DataFrame()
    return _records_to_frame(
        _get_json(
            f"/atlases/{int(atlas_id)}/expression",
            {
                "region_id": list(region_ids),
                "gene": [g for g in gene_symbols if g],
                "celltype": celltype or "",
            },
        )
    )


@st.cache_data(show_spinner=False, ttl=300)
def marker_evidence(atlas_id: int, region_ids: tuple[str, ...], celltype: str | None, limit: int) -> pd.DataFrame:
    if not is_configured():
        return pd.DataFrame()
    return _records_to_frame(
        _get_json(
            f"/atlases/{int(atlas_id)}/marker-evidence",
            {"region_id": list(region_ids), "celltype": celltype or "", "limit": int(limit)},
        )
    )
