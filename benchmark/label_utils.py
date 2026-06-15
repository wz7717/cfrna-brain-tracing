from __future__ import annotations

import json
from typing import Dict, Optional
from core.gene_utils import normalize_region_label, is_mixed_region_label

def safe_soft_label(label: Optional[str]) -> Optional[str]:
    if label is None: return None
    return normalize_region_label(str(label).strip())

def default_label_extractor(sample_row: Dict) -> Optional[str]:
    meta = sample_row.get('metadata')
    if not meta: return None
    try: obj = json.loads(meta)
    except Exception: return None
    for k in ('ground_truth_region', 'source_region', 'injury_region', 'label_region', 'true_source', 'surgery_region'):
        if k in obj and obj[k] is not None: return safe_soft_label(obj[k])
    return None

def is_mixed_label(label: Optional[str]) -> bool:
    return is_mixed_region_label(label)
