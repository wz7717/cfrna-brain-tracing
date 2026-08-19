from __future__ import annotations

import json
import numpy as np
from pathlib import Path

from scripts.verify_regenerated_npz import compare


def test_npz_comparison_checks_values_not_compressed_bytes(tmp_path: Path) -> None:
    canonical = tmp_path / "canonical.npz"
    regenerated = tmp_path / "regenerated.npz"
    np.savez_compressed(canonical, numbers=np.array([1.0, 2.0]), labels=np.array(["a", "b"]))
    np.savez_compressed(regenerated, numbers=np.array([1.0, 2.0]), labels=np.array(["a", "b"]))
    assert compare(canonical, regenerated, "canonical", "regenerated")["status"] == "PASS"
    np.savez_compressed(regenerated, numbers=np.array([1.0, 3.0]), labels=np.array(["a", "b"]))
    assert compare(canonical, regenerated, "canonical", "regenerated")["status"] == "BLOCKED: SCIENTIFIC_OUTPUT_DRIFT"


def test_npz_comparison_ignores_only_generation_timestamp_metadata(tmp_path: Path) -> None:
    canonical = tmp_path / "canonical.npz"
    regenerated = tmp_path / "regenerated.npz"
    np.savez_compressed(canonical, metadata=json.dumps({"created_at_utc": "old", "panel": "fixed"}))
    np.savez_compressed(regenerated, metadata=json.dumps({"created_at_utc": "new", "panel": "fixed"}))
    assert compare(canonical, regenerated, "canonical", "regenerated")["status"] == "PASS"
    np.savez_compressed(regenerated, metadata=json.dumps({"created_at_utc": "new", "panel": "changed"}))
    assert compare(canonical, regenerated, "canonical", "regenerated")["status"] == "BLOCKED: SCIENTIFIC_OUTPUT_DRIFT"
