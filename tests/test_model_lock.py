from pathlib import Path
import hashlib

import pytest

from core.model_lock import (
    LOCKED_PRODUCTION_PARAMETERS,
    ModelLockError,
    assert_locked_production_parameters,
    load_model_lock,
    verify_locked_artifact,
    verify_locked_model_bundle,
)


def test_committed_canonical110_model_bundle_matches_lock() -> None:
    manifest = verify_locked_model_bundle()

    assert manifest["status"] == "frozen"
    assert manifest["production_parameters"] == dict(LOCKED_PRODUCTION_PARAMETERS)


def test_locked_artifact_rejects_byte_change(tmp_path: Path) -> None:
    artifact = tmp_path / "model.bin"
    artifact.write_bytes(b"locked")
    expected = {
        "size": artifact.stat().st_size,
        "sha256": hashlib.sha256(artifact.read_bytes()).hexdigest(),
    }
    artifact.write_bytes(b"edited")

    with pytest.raises(ModelLockError, match="SHA-256 changed"):
        verify_locked_artifact(artifact, expected)


def test_lock_manifest_contains_exact_frozen_inventory() -> None:
    manifest = load_model_lock()

    assert len(manifest["artifacts"]) == 8
    assert all(len(item["sha256"]) == 64 for item in manifest["artifacts"].values())


def test_production_parameter_drift_is_rejected() -> None:
    changed = dict(LOCKED_PRODUCTION_PARAMETERS)
    changed["exact_top50_weight"] = 0.30

    with pytest.raises(ModelLockError, match="exact_top50_weight"):
        assert_locked_production_parameters(changed)
