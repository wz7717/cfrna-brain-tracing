from __future__ import annotations

from typing import Dict, Optional, Sequence
import re
import pandas as pd

ENSEMBL_RE = re.compile(r"^(ENS[A-Z]*G\d+)$", re.I)
SYMBOL_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]{1,30}$")
REGION_ALIASES: Dict[str, str] = {
    "CTX_V4": "CTX_ACC",
    "V4": "ACC",
    "CTX_ACC": "CTX_ACC",
    "ACC": "ACC",
}

def normalize_region_label(label: Optional[str]) -> Optional[str]:
    if label is None:
        return None
    s = str(label).strip()
    if not s:
        return None
    return REGION_ALIASES.get(s, s)

def is_mixed_region_label(label: Optional[str]) -> bool:
    if label is None:
        return False
    s = str(label)
    return any(x in s for x in ['+', ',', ';', '|', '/'])

def guess_gene_id_type(genes: Sequence[str]) -> str:
    vals = [str(g).strip() for g in genes if str(g).strip()]
    if not vals:
        return 'unknown'
    n = min(200, len(vals))
    vals = vals[:n]
    ens = sum(bool(ENSEMBL_RE.match(v)) for v in vals)
    sym = sum(bool(SYMBOL_RE.match(v)) and not bool(ENSEMBL_RE.match(v)) for v in vals)
    if ens / max(n, 1) >= 0.6:
        return 'ensembl_like'
    if sym / max(n, 1) >= 0.6:
        return 'symbol_like'
    return 'mixed'

def resolve_gene_alias(gene: str, alias_map: Optional[Dict[str, str]] = None) -> str:
    g = str(gene).strip()
    if alias_map and g in alias_map:
        return alias_map[g]
    return g
