from __future__ import annotations

import hashlib
import json
from pathlib import Path

from core.external_inputs import (
    deterministic_tree_hash,
    select_external_root,
    verify_inputs,
)


def _sha256(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def _write_manifest(tmp_path: Path, item: dict) -> Path:
    manifest = {
        "schema": "braintrace.external_input_manifest.v1",
        "canonical_external_root": "external_data",
        "profiles": {"full": {}, "portable": {}},
        "inputs": [item],
    }
    path = tmp_path / "manifest.json"
    path.write_text(json.dumps(manifest), encoding="utf-8")
    return path


def _file_item(expected: str) -> dict:
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


def test_external_root_priority_is_cli_then_environment_then_repository(tmp_path: Path) -> None:
    cli_root = tmp_path / "cli"
    env_root = tmp_path / "env"
    repository = tmp_path / "repository"
    assert select_external_root(cli_root, repo_root=repository, environ={"BRAINTRACE_EXTERNAL_DATA_ROOT": str(env_root)}).source == "cli"
    assert select_external_root(None, repo_root=repository, environ={"BRAINTRACE_EXTERNAL_DATA_ROOT": str(env_root)}).source == "environment"
    default = select_external_root(None, repo_root=repository, environ={})
    assert default.source == "repository_default"
    assert default.path == (repository / "external_data").resolve()


def test_missing_required_file_is_fail_closed(tmp_path: Path) -> None:
    manifest = _write_manifest(tmp_path, _file_item("0" * 64))
    report = verify_inputs(
        profile="full",
        external_data_root=tmp_path / "external_data",
        manifest_path=manifest,
        repo_root=tmp_path,
    )
    assert report["status"] == "FAIL"
    assert report["inputs"][0]["status"] == "MISSING"


def test_bad_hash_is_fail_closed(tmp_path: Path) -> None:
    root = tmp_path / "external_data"
    path = root / "Fixture" / "input.txt"
    path.parent.mkdir(parents=True)
    path.write_text("actual", encoding="utf-8")
    manifest = _write_manifest(tmp_path, _file_item(_sha256(b"different")))
    report = verify_inputs(profile="full", external_data_root=root, manifest_path=manifest, repo_root=tmp_path)
    assert report["status"] == "FAIL"
    assert report["inputs"][0]["status"] == "HASH_MISMATCH"


def test_invalid_directory_is_fail_closed(tmp_path: Path) -> None:
    root = tmp_path / "external_data"
    directory = root / "Fixture"
    directory.mkdir(parents=True)
    required = {"path": "required.txt", "size": 3, "sha256": _sha256(b"yes")}
    item = {
        "alias": "directory_fixture",
        "storage": "external",
        "kind": "directory",
        "canonical_path": "Fixture",
        "tree_sha256": "0" * 64,
        "allow_extra_files": False,
        "required_files": [required],
        "required": True,
        "profiles": ["full"],
        "classification": "downloadable public",
        "role": "test fixture",
        "accession_or_source": "test",
    }
    manifest = _write_manifest(tmp_path, item)
    report = verify_inputs(profile="full", external_data_root=root, manifest_path=manifest, repo_root=tmp_path)
    assert report["status"] == "FAIL"
    assert report["inputs"][0]["status"] == "INCOMPLETE_DIRECTORY"


def test_valid_file_is_pass_and_report_uses_portable_locator(tmp_path: Path) -> None:
    root = tmp_path / "external_data"
    path = root / "Fixture" / "input.txt"
    path.parent.mkdir(parents=True)
    payload = b"valid input\n"
    path.write_bytes(payload)
    manifest = _write_manifest(tmp_path, _file_item(_sha256(payload)))
    report = verify_inputs(profile="full", external_data_root=root, manifest_path=manifest, repo_root=tmp_path)
    assert report["status"] == "PASS"
    assert report["inputs"][0]["status"] == "PASS"
    assert report["inputs"][0]["locator"] == "external_source::fixture/input.txt"
    assert str(root) not in json.dumps(report)


def test_repository_text_hash_uses_canonical_lf_across_checkout_line_endings(tmp_path: Path) -> None:
    tracked = tmp_path / "tracked" / "input.csv"
    tracked.parent.mkdir(parents=True)
    tracked.write_bytes(b"header\r\nvalue\r\n")
    item = {
        "alias": "tracked_text",
        "storage": "repository",
        "kind": "file",
        "canonical_path": "tracked/input.csv",
        "sha256": _sha256(b"header\nvalue\n"),
        "hash_mode": "utf8_lf",
        "required": True,
        "profiles": ["portable"],
        "classification": "repository tracked",
        "role": "test tracked text",
        "accession_or_source": "test",
    }
    manifest = _write_manifest(tmp_path, item)
    report = verify_inputs(profile="portable", manifest_path=manifest, repo_root=tmp_path)
    assert report["status"] == "PASS"
    assert report["inputs"][0]["hash_mode"] == "utf8_lf"


def test_directory_tree_hash_is_deterministic_and_rejects_extra_files(tmp_path: Path) -> None:
    root = tmp_path / "external_data"
    directory = root / "Fixture"
    directory.mkdir(parents=True)
    required_path = directory / "required.txt"
    required_path.write_bytes(b"yes")
    tree_hash, _ = deterministic_tree_hash(directory)
    item = {
        "alias": "directory_fixture",
        "storage": "external",
        "kind": "directory",
        "canonical_path": "Fixture",
        "tree_sha256": tree_hash,
        "allow_extra_files": False,
        "required_files": [{"path": "required.txt", "size": 3, "sha256": _sha256(b"yes")}],
        "required": True,
        "profiles": ["full"],
        "classification": "downloadable public",
        "role": "test fixture",
        "accession_or_source": "test",
    }
    manifest = _write_manifest(tmp_path, item)
    assert verify_inputs(profile="full", external_data_root=root, manifest_path=manifest, repo_root=tmp_path)["status"] == "PASS"
    (directory / "unexpected.txt").write_text("no", encoding="utf-8")
    report = verify_inputs(profile="full", external_data_root=root, manifest_path=manifest, repo_root=tmp_path)
    assert report["status"] == "FAIL"
    assert report["inputs"][0]["status"] == "UNEXPECTED_REQUIRED_FILE"


def test_release_gate_rejects_legacy_fallback(tmp_path: Path) -> None:
    legacy = tmp_path / "legacy" / "input.txt"
    legacy.parent.mkdir(parents=True)
    payload = b"legacy input"
    legacy.write_bytes(payload)
    item = _file_item(_sha256(payload))
    item["legacy_path"] = "legacy/input.txt"
    manifest = _write_manifest(tmp_path, item)
    warning_report = verify_inputs(profile="full", manifest_path=manifest, repo_root=tmp_path)
    assert warning_report["inputs"][0]["path_mode"] == "legacy_fallback"
    assert warning_report["status"] == "PASS"
    release_report = verify_inputs(profile="full", manifest_path=manifest, repo_root=tmp_path, release_gate=True)
    assert release_report["status"] == "FAIL"
    assert release_report["inputs"][0]["status"] == "SOURCE_CONTRACT_ERROR"
