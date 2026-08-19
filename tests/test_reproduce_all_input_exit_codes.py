from __future__ import annotations

import hashlib
import json
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def _sha256(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def _write_manifest(tmp_path: Path, item: dict[str, object]) -> Path:
    manifest = {
        "schema": "braintrace.external_input_manifest.v1",
        "canonical_external_root": "external_data",
        "profiles": {"full": {}, "portable": {}},
        "inputs": [item],
    }
    path = tmp_path / "manifest.json"
    path.write_text(json.dumps(manifest), encoding="utf-8")
    return path


def _file_item(expected: str) -> dict[str, object]:
    return {
        "alias": "fixture",
        "storage": "external",
        "kind": "file",
        "canonical_path": "Fixture/input.txt",
        "sha256": expected,
        "required": True,
        "profiles": ["full"],
        "classification": "downloadable public",
        "role": "test fixture",
        "accession_or_source": "test",
    }


def _run(tmp_path: Path, manifest: Path, root: Path) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [
            sys.executable,
            str(ROOT / "reproduce_all.py"),
            "--profile",
            "full",
            "--verify-only",
            "--input-manifest",
            str(manifest),
            "--external-data-root",
            str(root),
            "--output-dir",
            str(tmp_path / "audit"),
        ],
        cwd=ROOT,
        text=True,
        capture_output=True,
        check=False,
    )


def test_runner_missing_input_returns_nonzero(tmp_path: Path) -> None:
    result = _run(tmp_path, _write_manifest(tmp_path, _file_item("0" * 64)), tmp_path / "external_data")
    assert result.returncode != 0


def test_runner_bad_hash_returns_nonzero(tmp_path: Path) -> None:
    root = tmp_path / "external_data"
    path = root / "Fixture" / "input.txt"
    path.parent.mkdir(parents=True)
    path.write_bytes(b"actual")
    result = _run(tmp_path, _write_manifest(tmp_path, _file_item(_sha256(b"different"))), root)
    assert result.returncode != 0


def test_runner_incomplete_directory_returns_nonzero(tmp_path: Path) -> None:
    root = tmp_path / "external_data"
    (root / "Fixture").mkdir(parents=True)
    item = {
        "alias": "directory_fixture",
        "storage": "external",
        "kind": "directory",
        "canonical_path": "Fixture",
        "tree_sha256": "0" * 64,
        "allow_extra_files": False,
        "required_files": [{"path": "required.txt", "size": 3, "sha256": _sha256(b"yes")}],
        "required": True,
        "profiles": ["full"],
        "classification": "downloadable public",
        "role": "test fixture",
        "accession_or_source": "test",
    }
    result = _run(tmp_path, _write_manifest(tmp_path, item), root)
    assert result.returncode != 0


def test_runner_valid_input_returns_zero(tmp_path: Path) -> None:
    root = tmp_path / "external_data"
    path = root / "Fixture" / "input.txt"
    payload = b"valid input\n"
    path.parent.mkdir(parents=True)
    path.write_bytes(payload)
    result = _run(tmp_path, _write_manifest(tmp_path, _file_item(_sha256(payload))), root)
    assert result.returncode == 0, result.stdout + result.stderr
