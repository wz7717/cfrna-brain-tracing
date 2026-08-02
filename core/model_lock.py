from __future__ import annotations

import hashlib
import json
from pathlib import Path
from types import MappingProxyType
from typing import Any, Mapping


ROOT = Path(__file__).resolve().parents[1]
MODEL_LOCK_PATH = ROOT / "data" / "models" / "canonical110_model_lock.json"
EXPECTED_MODEL_LOCK_ID = "canonical110-v0.1.10-20260802"

LOCKED_PRODUCTION_PARAMETERS: Mapping[str, Any] = MappingProxyType(
    {
        "route_name": "projected_vsd_network_top3_logcpm_resolution_local_exact",
        "region_count": 110,
        "network_count": 10,
        "beam_count": 120,
        "network_top_k": 3,
        "network_gene_count": 200,
        "network_min_overlap_fraction": 0.50,
        "project_to_vsd": True,
        "enable_pairwise_rescue": False,
        "region_local_top_n_genes": 200,
        "region_min_overlap_genes": 20,
        "exact_top50_gene_count": 50,
        "exact_top100_gene_count": 100,
        "exact_top50_weight": 0.25,
        "resolution_group_mean_weight": 0.10,
        "allow_development_fallback": False,
    }
)

LOCKED_ARTIFACT_PATHS = frozenset(
    {
        "data/models/bo2023_saleem_network_top200_model.npz",
        "data/models/bo2023_saleem_network_top200_model.json",
        "data/models/bo2023_saleem_network_top200_model_genes.csv",
        "data/models/bo2023_saleem_network_pairwise_rescue_model.json",
        "data/models/bo2023_reference_projector_linear_full.npz",
        "data/models/bo2023_formal_region_logcpm_reference_matrix.npz",
        "data/models/bo2023_formal_region_beam_gene_panels.json",
        "data/models/bo2023_region_resolution_groups.json",
    }
)


class ModelLockError(RuntimeError):
    """Raised when the frozen production model differs from its lock."""


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_model_lock(manifest_path: Path = MODEL_LOCK_PATH) -> dict[str, Any]:
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ModelLockError(f"model lock manifest is unavailable or invalid: {exc}") from exc
    if manifest.get("lock_id") != EXPECTED_MODEL_LOCK_ID:
        raise ModelLockError("model lock ID differs from the frozen canonical 110 release")
    if manifest.get("status") != "frozen":
        raise ModelLockError("model lock status is not frozen")
    if manifest.get("production_parameters") != dict(LOCKED_PRODUCTION_PARAMETERS):
        raise ModelLockError("production parameters differ from the frozen canonical 110 route")
    artifacts = manifest.get("artifacts")
    if not isinstance(artifacts, dict) or set(artifacts) != LOCKED_ARTIFACT_PATHS:
        raise ModelLockError("locked artifact inventory is incomplete or contains unexpected paths")
    return manifest


def verify_locked_artifact(path: Path, expected: Mapping[str, Any]) -> None:
    if not path.is_file():
        raise ModelLockError(f"locked model artifact is missing: {path}")
    expected_size = int(expected.get("size", -1))
    if path.stat().st_size != expected_size:
        raise ModelLockError(f"locked model artifact size changed: {path}")
    expected_hash = str(expected.get("sha256", "")).lower()
    if len(expected_hash) != 64 or _sha256(path) != expected_hash:
        raise ModelLockError(f"locked model artifact SHA-256 changed: {path}")


def verify_locked_model_bundle(
    root: Path = ROOT,
    manifest_path: Path = MODEL_LOCK_PATH,
) -> dict[str, Any]:
    manifest = load_model_lock(manifest_path)
    resolved_root = root.resolve()
    for relative_path, expected in manifest["artifacts"].items():
        path = (resolved_root / relative_path).resolve()
        if resolved_root not in path.parents:
            raise ModelLockError(f"locked artifact escapes repository root: {relative_path}")
        verify_locked_artifact(path, expected)
    return manifest


def assert_locked_production_parameters(actual: Mapping[str, Any]) -> None:
    actual_parameters = dict(actual)
    expected_parameters = dict(LOCKED_PRODUCTION_PARAMETERS)
    if actual_parameters != expected_parameters:
        changed = sorted(
            key
            for key in set(actual_parameters) | set(expected_parameters)
            if actual_parameters.get(key) != expected_parameters.get(key)
        )
        raise ModelLockError(
            "production model parameters differ from the frozen canonical 110 route: "
            + ", ".join(changed)
        )
