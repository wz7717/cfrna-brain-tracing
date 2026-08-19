"""Canonical external-input contract for BrainTrace manuscript reproduction.

The repository deliberately does not distribute controlled or large raw inputs.
This module is the single resolver and verifier for those inputs.  It records
portable logical locators in audit output; physical paths never enter public
provenance reports.
"""

from __future__ import annotations

import hashlib
import json
import os
import warnings
import zipfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

from core.provenance_hashes import sha256_utf8_lf_text


REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MANIFEST = REPO_ROOT / "reproducibility" / "external_input_manifest.json"
ENV_EXTERNAL_DATA_ROOT = "BRAINTRACE_EXTERNAL_DATA_ROOT"

FAILURE_STATUSES = frozenset(
    {
        "MISSING",
        "HASH_MISMATCH",
        "NOT_IN_MANIFEST",
        "WRONG_TYPE",
        "INVALID_ARCHIVE",
        "INCOMPLETE_DIRECTORY",
        "SOURCE_CONTRACT_ERROR",
        "UNEXPECTED_REQUIRED_FILE",
        "INVALID_SCHEMA",
    }
)


class InputContractError(ValueError):
    """Raised when the tracked input contract itself is malformed."""


@dataclass(frozen=True)
class RootSelection:
    path: Path
    source: str


def sha256_file(path: Path) -> str:
    """Return the lower-case SHA-256 of a regular file."""

    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def deterministic_tree_hash(path: Path) -> tuple[str, list[dict[str, Any]]]:
    """Hash a directory as sorted ``relative-path<TAB>size<TAB>sha256`` rows."""

    if not path.is_dir():
        raise ValueError(f"not a directory: {path}")
    records: list[dict[str, Any]] = []
    for child in sorted((item for item in path.rglob("*") if item.is_file()), key=lambda item: item.as_posix()):
        relative = child.relative_to(path).as_posix()
        records.append(
            {
                "path": relative,
                "size": child.stat().st_size,
                "sha256": sha256_file(child),
            }
        )
    payload = "".join(
        f"{record['path']}\t{record['size']}\t{record['sha256']}\n" for record in records
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest(), records


def load_manifest(path: Path = DEFAULT_MANIFEST) -> dict[str, Any]:
    """Load and validate the tracked, machine-readable external-input manifest."""

    try:
        manifest = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise InputContractError(f"cannot read external-input manifest: {exc}") from exc
    if not isinstance(manifest, dict) or manifest.get("schema") != "braintrace.external_input_manifest.v1":
        raise InputContractError("external-input manifest schema must be braintrace.external_input_manifest.v1")
    if not isinstance(manifest.get("inputs"), list) or not manifest["inputs"]:
        raise InputContractError("external-input manifest must contain a non-empty inputs list")
    aliases: set[str] = set()
    for item in manifest["inputs"]:
        if not isinstance(item, dict):
            raise InputContractError("each external-input manifest item must be an object")
        alias = item.get("alias")
        if not isinstance(alias, str) or not alias or alias in aliases:
            raise InputContractError("external-input aliases must be unique non-empty strings")
        aliases.add(alias)
        if item.get("storage") not in {"external", "repository"}:
            raise InputContractError(f"{alias}: storage must be external or repository")
        if item.get("kind") not in {"file", "directory"}:
            raise InputContractError(f"{alias}: kind must be file or directory")
        if not isinstance(item.get("canonical_path"), str) or not item["canonical_path"]:
            raise InputContractError(f"{alias}: canonical_path is required")
        if not isinstance(item.get("profiles"), list) or not item["profiles"]:
            raise InputContractError(f"{alias}: profiles is required")
        if item["kind"] == "file" and not _valid_digest(item.get("sha256")):
            raise InputContractError(f"{alias}: file SHA-256 is required")
        hash_mode = item.get("hash_mode", "raw")
        if hash_mode not in {"raw", "utf8_lf"}:
            raise InputContractError(f"{alias}: hash_mode must be raw or utf8_lf")
        if item["kind"] == "directory":
            if not _valid_digest(item.get("tree_sha256")):
                raise InputContractError(f"{alias}: directory tree_sha256 is required")
            if not isinstance(item.get("required_files"), list) or not item["required_files"]:
                raise InputContractError(f"{alias}: directory required_files is required")
            for required in item["required_files"]:
                if not isinstance(required, dict) or not isinstance(required.get("path"), str):
                    raise InputContractError(f"{alias}: malformed directory required file")
                if not isinstance(required.get("size"), int) or not _valid_digest(required.get("sha256")):
                    raise InputContractError(f"{alias}: directory required file needs size and SHA-256")
    return manifest


def _valid_digest(value: object) -> bool:
    return isinstance(value, str) and len(value) == 64 and all(char in "0123456789abcdefABCDEF" for char in value)


def select_external_root(
    external_data_root: Path | str | None = None,
    *,
    repo_root: Path = REPO_ROOT,
    environ: dict[str, str] | None = None,
) -> RootSelection:
    """Resolve the canonical-root priority: CLI, environment, repository default."""

    environment = os.environ if environ is None else environ
    if external_data_root is not None:
        return RootSelection(Path(external_data_root).expanduser().resolve(), "cli")
    environment_root = environment.get(ENV_EXTERNAL_DATA_ROOT, "").strip()
    if environment_root:
        return RootSelection(Path(environment_root).expanduser().resolve(), "environment")
    return RootSelection((repo_root / "external_data").resolve(), "repository_default")


def input_locator(item: dict[str, Any]) -> str:
    """Return a portable, stable source locator suitable for public artifacts."""

    alias = str(item["alias"])
    filename = Path(str(item["canonical_path"])).name
    if item["storage"] == "repository":
        return f"repository::{item['canonical_path']}"
    return f"external_source::{alias}/{filename}"


def resolve_input_path(
    item: dict[str, Any],
    selection: RootSelection,
    *,
    repo_root: Path = REPO_ROOT,
    allow_legacy_fallback: bool = True,
) -> tuple[Path, str]:
    """Resolve an input without persisting the physical path in any report.

    Legacy locations are only considered after a canonical external-root lookup
    misses.  They remain intentionally marked so release gates can reject them.
    """

    if item["storage"] == "repository":
        return repo_root / str(item["canonical_path"]), "repository"
    canonical = selection.path / str(item["canonical_path"])
    if canonical.exists() or not allow_legacy_fallback:
        return canonical, "canonical"
    legacy_path = item.get("legacy_path")
    if isinstance(legacy_path, str) and legacy_path:
        legacy = repo_root / legacy_path
        if legacy.exists():
            return legacy, "legacy_fallback"
    return canonical, "canonical"


def _verify_archive(path: Path, archive_type: str | None) -> str | None:
    if archive_type != "zip":
        return None
    if not zipfile.is_zipfile(path):
        return "INVALID_ARCHIVE"
    try:
        with zipfile.ZipFile(path) as archive:
            return "INVALID_ARCHIVE" if archive.testzip() is not None else None
    except (OSError, zipfile.BadZipFile):
        return "INVALID_ARCHIVE"


def _verify_file(item: dict[str, Any], path: Path) -> tuple[str, dict[str, Any]]:
    if not path.exists():
        return "MISSING", {}
    if not path.is_file():
        return "WRONG_TYPE", {}
    expected = item.get("sha256")
    if not _valid_digest(expected):
        return "NOT_IN_MANIFEST", {}
    hash_mode = str(item.get("hash_mode", "raw"))
    observed = (
        sha256_utf8_lf_text(path).lower()
        if hash_mode == "utf8_lf"
        else sha256_file(path)
    )
    if observed.lower() != str(expected).lower():
        return "HASH_MISMATCH", {"observed_sha256": observed, "expected_sha256": str(expected).lower()}
    archive_status = _verify_archive(path, item.get("archive"))
    if archive_status:
        return archive_status, {"sha256": observed}
    return "PASS", {"sha256": observed, "size": path.stat().st_size, "hash_mode": hash_mode}


def _verify_directory(item: dict[str, Any], path: Path) -> tuple[str, dict[str, Any]]:
    if not path.exists():
        return "MISSING", {}
    if not path.is_dir():
        return "WRONG_TYPE", {}
    required_files = item.get("required_files")
    tree_sha256 = item.get("tree_sha256")
    if not isinstance(required_files, list) or not _valid_digest(tree_sha256):
        return "NOT_IN_MANIFEST", {}
    observed_required: set[str] = set()
    for required in required_files:
        if not isinstance(required, dict) or not isinstance(required.get("path"), str):
            return "INVALID_SCHEMA", {}
        relative = required["path"].replace("\\", "/")
        candidate = path / relative
        observed_required.add(relative)
        if not candidate.exists() or not candidate.is_file():
            return "INCOMPLETE_DIRECTORY", {"missing_required_file": relative}
        if candidate.stat().st_size != required.get("size"):
            return "INCOMPLETE_DIRECTORY", {"size_mismatch_file": relative}
        observed = sha256_file(candidate)
        if observed.lower() != str(required.get("sha256", "")).lower():
            return "HASH_MISMATCH", {"hash_mismatch_file": relative, "observed_sha256": observed}
        archive_status = _verify_archive(candidate, required.get("archive"))
        if archive_status:
            return archive_status, {"invalid_archive_file": relative}
    observed_tree, records = deterministic_tree_hash(path)
    actual_files = {str(record["path"]) for record in records}
    if item.get("allow_extra_files") is False and actual_files != observed_required:
        return "UNEXPECTED_REQUIRED_FILE", {
            "unexpected_files": sorted(actual_files - observed_required),
            "missing_files": sorted(observed_required - actual_files),
        }
    if observed_tree.lower() != str(tree_sha256).lower():
        return "HASH_MISMATCH", {"observed_tree_sha256": observed_tree, "expected_tree_sha256": str(tree_sha256).lower()}
    return "PASS", {"tree_sha256": observed_tree, "n_files": len(records)}


def verify_inputs(
    *,
    profile: str,
    external_data_root: Path | str | None = None,
    manifest_path: Path = DEFAULT_MANIFEST,
    release_gate: bool = False,
    repo_root: Path = REPO_ROOT,
    environ: dict[str, str] | None = None,
) -> dict[str, Any]:
    """Verify every required input for a profile and return a portable audit.

    The result is deliberately data-only: physical root paths, host names and
    user names are omitted so it can be safely committed or archived.
    """

    try:
        manifest = load_manifest(manifest_path)
    except InputContractError as exc:
        return {
            "schema": "braintrace.external_input_verification.v1",
            "profile": profile,
            "status": "FAIL",
            "failure_count": 1,
            "inputs": [{"alias": "manifest", "status": "INVALID_SCHEMA", "detail": str(exc)}],
        }
    declared_profiles = manifest.get("profiles", {})
    if profile not in declared_profiles:
        return {
            "schema": "braintrace.external_input_verification.v1",
            "profile": profile,
            "status": "FAIL",
            "failure_count": 1,
            "inputs": [{"alias": "profile", "status": "SOURCE_CONTRACT_ERROR", "detail": "unknown profile"}],
        }
    selection = select_external_root(external_data_root, repo_root=repo_root, environ=environ)
    records: list[dict[str, Any]] = []
    for item in manifest["inputs"]:
        if profile not in item["profiles"]:
            continue
        path, path_mode = resolve_input_path(item, selection, repo_root=repo_root, allow_legacy_fallback=True)
        if path_mode == "legacy_fallback" and release_gate:
            status, details = "SOURCE_CONTRACT_ERROR", {"detail": "legacy fallback is prohibited by the release gate"}
        elif item["kind"] == "file":
            status, details = _verify_file(item, path)
        else:
            status, details = _verify_directory(item, path)
        records.append(
            {
                "alias": item["alias"],
                "locator": input_locator(item),
                "storage": item["storage"],
                "kind": item["kind"],
                "required": bool(item.get("required", False)),
                "path_mode": path_mode,
                "status": status,
                **details,
            }
        )
    failures = [record for record in records if record["status"] in FAILURE_STATUSES]
    return {
        "schema": "braintrace.external_input_verification.v1",
        "profile": profile,
        "status": "PASS" if not failures else "FAIL",
        "failure_count": len(failures),
        "root_selection": selection.source,
        "canonical_root_locator": "external_data/",
        "legacy_fallback_used": [record["alias"] for record in records if record["path_mode"] == "legacy_fallback"],
        "warnings": [
            f"WARNING: {record['alias']} used the deprecated legacy fallback; canonical external_data is required for release."
            for record in records
            if record["path_mode"] == "legacy_fallback"
        ],
        "inputs": records,
    }


def input_by_alias(alias: str, manifest_path: Path = DEFAULT_MANIFEST) -> dict[str, Any]:
    """Return a manifest item by its stable alias."""

    for item in load_manifest(manifest_path)["inputs"]:
        if item["alias"] == alias:
            return item
    raise KeyError(alias)


def resolve_alias(
    alias: str,
    *,
    external_data_root: Path | str | None = None,
    release_gate: bool = False,
    repo_root: Path = REPO_ROOT,
    manifest_path: Path = DEFAULT_MANIFEST,
) -> Path:
    """Resolve one source for scripts that take explicit input-path arguments."""

    item = input_by_alias(alias, manifest_path=manifest_path)
    path, mode = resolve_input_path(
        item,
        select_external_root(external_data_root, repo_root=repo_root),
        repo_root=repo_root,
        allow_legacy_fallback=not release_gate,
    )
    if mode == "legacy_fallback" and release_gate:
        raise InputContractError(f"{alias}: legacy fallback is prohibited by the release gate")
    if mode == "legacy_fallback":
        warnings.warn(
            f"WARNING: {alias} resolved via deprecated legacy fallback; use canonical external_data instead.",
            RuntimeWarning,
            stacklevel=2,
        )
    return path


def portable_origin_locator(path: Path, *, alias: str) -> str:
    """Convert an external historical source into a non-machine-specific locator."""

    try:
        return f"repository::{path.resolve().relative_to(REPO_ROOT.resolve()).as_posix()}"
    except ValueError:
        return f"external_source::{alias}/{path.name}"
