#!/usr/bin/env python
"""Manifest-driven, one-shot BrainTrace release metadata finalization.

The script is intentionally unable to finalize an engineering pre-release.
``--finalize`` requires a separately recorded all-green release gate, then
updates every version-bearing release surface from one manifest in one run.
"""

from __future__ import annotations

import argparse
import json
import re
import tomllib
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MANIFEST = ROOT / "release" / "v0.1.18" / "release_manifest.json"
SOFTWARE_DOI = "10.5281/zenodo.22022058"
FULL_REPRO_DOI = "10.5281/zenodo.22022059"
SOFTWARE_CONCEPT_DOI = "10.5281/zenodo.20773674"
FULL_REPRO_CONCEPT_DOI = "10.5281/zenodo.21920696"
PREVIOUS_RELEASE_SHA = "2cceaea5c42979689e96f97f75b66eaaeafa4629"
PROVENANCE_BODY_MARKERS = {
    "DATA_PROVENANCE.md": "## Source Datasets",
    "reproducibility/DATA_PROVENANCE.md": "## 1. Primary Atlas Data",
}
PUBLIC_EXAMPLE_OUTPUTS = (
    "examples/expected_output_counts.json",
    "examples/expected_output_logcpm.json",
)


def load_manifest(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def validate_manifest(manifest: dict[str, Any]) -> list[str]:
    errors: list[str] = []
    version = str(manifest.get("version", ""))
    if not re.fullmatch(r"0\.1\.\d+", version):
        errors.append("version must use the 0.1.x form")
    if manifest.get("tag") != f"v{version}":
        errors.append("tag must equal v<version>")
    if manifest.get("software_version_doi") != SOFTWARE_DOI:
        errors.append(f"software DOI is not the reserved v{version} draft DOI")
    if manifest.get("full_repro_version_doi") != FULL_REPRO_DOI:
        errors.append(f"full-repro DOI is not the reserved v{version} draft DOI")
    if manifest.get("software_concept_doi") != SOFTWARE_CONCEPT_DOI:
        errors.append("software concept DOI does not match the existing BrainTrace series")
    if manifest.get("full_repro_concept_doi") != FULL_REPRO_CONCEPT_DOI:
        errors.append("full-repro concept DOI does not match the existing BrainTrace series")
    previous = manifest.get("previous_release")
    if not isinstance(previous, dict) or previous.get("version") != "0.1.17":
        errors.append("previous release must identify immutable v0.1.17")
    elif (
        previous.get("tag") != "v0.1.17"
        or previous.get("git_sha") != PREVIOUS_RELEASE_SHA
        or previous.get("software_version_doi") != "10.5281/zenodo.22006038"
        or previous.get("full_repro_version_doi") != "10.5281/zenodo.22005947"
        or previous.get("immutable") is not True
    ):
        errors.append("immutable v0.1.17 identity changed")
    if manifest.get("scientific_frozen") is not True:
        errors.append("scientific_frozen must be true")
    if not isinstance(manifest.get("creators"), list) or not manifest["creators"]:
        errors.append("creators must be a non-empty authoritative list")
    if not str(manifest.get("title", "")).strip() or not str(manifest.get("description", "")).strip():
        errors.append("title and description are required")
    return errors


def _replace_once(text: str, pattern: str, replacement: str, label: str) -> str:
    updated, count = re.subn(pattern, replacement, text, count=1, flags=re.MULTILINE)
    if count != 1:
        raise ValueError(f"could not uniquely update {label}")
    return updated


def _replace_release_state_text(text: str, previous: str, finalized: str, label: str) -> str:
    """Apply one idempotent, generator-owned release-state rewrite."""
    if previous in text:
        return text.replace(previous, finalized)
    if finalized in text:
        return text
    raise ValueError(f"could not locate {label} for release-state finalization")


def _write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8", newline="\n")


def _creators(manifest: dict[str, Any]) -> list[dict[str, str]]:
    return [{"name": str(record["name"])} for record in manifest["creators"]]


def _zenodo_payload(manifest: dict[str, Any]) -> dict[str, Any]:
    return {
        "title": manifest.get("software_record_title", manifest["title"]),
        "description": manifest["description"],
        "creators": _creators(manifest),
        "keywords": list(manifest.get("keywords", [])),
        "license": manifest["license"],
        "version": f"v{manifest['version']}",
        "related_identifiers": [
            {
                "identifier": manifest["repository"],
                "relation": "isSupplementTo",
                "scheme": "url",
            },
            {
                "identifier": f"https://doi.org/{manifest['software_version_doi']}",
                "relation": "isIdenticalTo",
                "scheme": "url",
            },
            {
                "identifier": f"https://doi.org/{manifest['full_repro_version_doi']}",
                "relation": "isSupplementTo",
                "scheme": "url",
            },
        ],
    }


def _release_block(manifest: dict[str, Any]) -> str:
    return (
        "<!-- BRAINTRACE:CURRENT_RELEASE:START -->\n"
        "## Current release\n\n"
        f"BrainTrace v{manifest['version']} is the current software release. "
        f"Software DOI: [https://doi.org/{manifest['software_version_doi']}](https://doi.org/{manifest['software_version_doi']}). "
        f"Full reproducibility DOI: [https://doi.org/{manifest['full_repro_version_doi']}](https://doi.org/{manifest['full_repro_version_doi']}).\n"
        "<!-- BRAINTRACE:CURRENT_RELEASE:END -->\n"
    )


def _data_provenance_header(manifest: dict[str, Any]) -> str:
    return (
        f"# Data Provenance - BrainTrace v{manifest['version']}\n\n"
        "This manifest documents every external dataset used to build, validate, and benchmark the "
        "BrainTrace hierarchical brain-origin candidate ranking tool. The current executable metadata and immutable "
        f"software release are v{manifest['version']}, archived at version DOI "
        f"[{manifest['software_version_doi']}](https://doi.org/{manifest['software_version_doi']}). "
        "Its full reproducibility archive, including the materialized required Git LFS payload, is archived separately at "
        f"[{manifest['full_repro_version_doi']}](https://doi.org/{manifest['full_repro_version_doi']}); "
        f"the persistent software concept DOI is [{manifest['software_concept_doi']}](https://doi.org/{manifest['software_concept_doi']}). "
        "Public source data are available under the cited accessions except where explicitly noted below."
    )


def _release_notes(manifest: dict[str, Any]) -> str:
    return (
        f"# BrainTrace v{manifest['version']}\n\n"
        f"- Software DOI: [https://doi.org/{manifest['software_version_doi']}](https://doi.org/{manifest['software_version_doi']})\n"
        f"- Full reproducibility DOI: [https://doi.org/{manifest['full_repro_version_doi']}](https://doi.org/{manifest['full_repro_version_doi']})\n"
        f"- Model lock: `{manifest['model_lock_id']}`\n"
        "- Release integrity, portable reproduction, and public-package compliance gates were strengthened.\n"
        "- Frozen model artifacts and formal scientific outputs are unchanged from v0.1.17.\n"
        "- External input verification remains fail-closed under the canonical `external_data/` layout.\n"
    )


def _check_release_gate(gate_path: Path) -> list[str]:
    """Validate facts that exist before the one-shot metadata finalization.

    The exact release-archive checksum is intentionally excluded here: the
    archive can only be built from a clean tree *after* finalization.  It is
    verified as a separate, mandatory pre-merge release-artifact gate.
    """
    gate = json.loads(gate_path.read_text(encoding="utf-8"))
    expected = {
        "status": "PASS",
        "scientific_drift": 0,
        "app_docker": "PASS",
        "repro_docker": "PASS",
        "github_actions": "GREEN",
    }
    return [f"release gate {key!r} is not {value!r}" for key, value in expected.items() if gate.get(key) != value]


def _finalize_readme_release_state(text: str, manifest: dict[str, Any]) -> str:
    version = str(manifest["version"])
    text = _replace_release_state_text(
        text,
        "### Huang 2025 cfRNA provenance remediation (unreleased working package)",
        "### Huang 2025 cfRNA provenance remediation",
        "Huang section heading",
    )
    text = _replace_release_state_text(
        text,
        "`reproducibility/huang_2025/` is a provenance-remediated external-audit package\n"
        "included in the v0.1.17 release.",
        f"`reproducibility/huang_2025/` is a provenance-remediated external-audit package\n"
        f"included in the v{version} release.",
        "Huang release state",
    )
    text = _replace_release_state_text(
        text,
        "The package supports technical portability/domain-shift statements only; it\n"
        "does not validate anatomical localization, tumour-source discrimination,\n"
        "patient-level stability, or clinical performance. Its formal release requires\n"
        "explicit scientific approval.",
        "The package supports technical portability/domain-shift statements only; it\n"
        "does not validate anatomical localization, tumour-source discrimination,\n"
        "patient-level stability, or clinical performance.",
        "Huang approval wording",
    )
    text = _replace_release_state_text(
        text,
        "BrainTrace v0.1.17 is the current software release. It does not alter the\n"
        "frozen model or the formal scientific results retained from v0.1.16.",
        f"BrainTrace v{version} is the current software release. It does not alter the\n"
        "frozen model or the formal scientific results retained from v0.1.17.",
        "current release status",
    )
    text = _replace_release_state_text(
        text,
        "The current executable software metadata is BrainTrace v0.1.17, the\n"
        "release-engineering and full-reproducibility finalization for the\n"
        "Bioinformatics Application Note. It does not change the frozen model,\n"
        "ontology, formal prediction set, Network Top1/Top3, resolution-group or\n"
        "exact-region endpoints. The production model remains locked under\n"
        "`canonical110-v0.1.12-20260813`.",
        f"The current executable software metadata is BrainTrace v{version}, the\n"
        "public-package integrity and full-reproducibility patch for the\n"
        "Bioinformatics Application Note. It does not change the frozen model,\n"
        "ontology, formal prediction set, Network Top1/Top3, resolution-group or\n"
        "exact-region endpoints. The production model remains locked under\n"
        "`canonical110-v0.1.12-20260813`.",
        "README current metadata status",
    )
    text = _replace_release_state_text(
        text,
        "primary endpoint or benchmark result changed in v0.1.17.",
        f"primary endpoint or benchmark result changed in v{version}.",
        "README frozen-result version",
    )
    text = _replace_release_state_text(
        text,
        "BrainTrace v0.1.16 is the current immutable GitHub/Zenodo software release under\n"
        "version DOI `https://doi.org/10.5281/zenodo.21974954`. BrainTrace v0.1.15 remains\n"
        "the previous immutable release under version DOI\n"
        "`https://doi.org/10.5281/zenodo.21970252`. The frozen v0.1.12 scientific\n"
        "release remains archived under `https://doi.org/10.5281/zenodo.21911532`; the\n"
        "v0.1.14 historical software record remains at\n"
        "`https://doi.org/10.5281/zenodo.21920261`; the persistent Zenodo concept DOI is\n"
        "`https://doi.org/10.5281/zenodo.20773674`.\n\n"
        "Because Zenodo's GitHub-generated source archive does not materialize Git LFS\n"
        "objects, the v0.1.16 full reproducibility archive, including the unchanged\n"
        "164,161,292-byte payload, is deposited separately under\n"
        "`https://doi.org/10.5281/zenodo.21974991`. The same materialized archive is\n"
        "available as a checksum-matched asset on the\n"
        "[GitHub v0.1.16 release](https://github.com/wz7717/cfrna-brain-tracing/releases/tag/v0.1.16).\n"
        "Its size is 198,284,868 bytes and its SHA-256 is\n"
        "`742c0390aa413deb18a12d6ae6164df52c5ed06e45e5c0e4218cd996231ca24e`.",
        f"BrainTrace v{version} is the current immutable GitHub/Zenodo software release\n"
        f"under software DOI `https://doi.org/{manifest['software_version_doi']}`; its\n"
        "separate materialized full reproducibility archive is under\n"
        f"`https://doi.org/{manifest['full_repro_version_doi']}`. Exact filenames,\n"
        f"inventories and SHA-256 values are released alongside the GitHub v{version}\n"
        "release and verified against the corresponding Zenodo records. BrainTrace\n"
        "v0.1.17 remains a historical immutable software/full-reproducibility release\n"
        "under `https://doi.org/10.5281/zenodo.22006038` /\n"
        "`https://doi.org/10.5281/zenodo.22005947`. BrainTrace v0.1.16 remains a\n"
        "historical immutable software/full-reproducibility release\n"
        "under `https://doi.org/10.5281/zenodo.21974954` /\n"
        "`https://doi.org/10.5281/zenodo.21974991`. BrainTrace v0.1.15 remains the\n"
        "previous immutable release under version DOI\n"
        "`https://doi.org/10.5281/zenodo.21970252`. The frozen v0.1.12 scientific\n"
        "release remains archived under `https://doi.org/10.5281/zenodo.21911532`; the\n"
        "v0.1.14 historical software record remains at\n"
        "`https://doi.org/10.5281/zenodo.21920261`; the persistent Zenodo concept DOI is\n"
        "`https://doi.org/10.5281/zenodo.20773674`.",
        "README archive status",
    )
    return text


def _finalize_provenance_release_state(relative: str, text: str, manifest: dict[str, Any]) -> str:
    version = str(manifest["version"])
    if relative == "DATA_PROVENANCE.md":
        return _replace_release_state_text(
            text,
            "The v0.1.17 release records AHBA\n",
            f"The v{version} release records AHBA\n",
            "DATA_PROVENANCE current release state",
        )
    text = _replace_release_state_text(
        text,
        "## 4. Manuscript and Candidate Provenance Artifacts (Dynamically Generated)",
        "## 4. Manuscript and Provenance Artifacts (Dynamically Generated)",
        "reproducibility provenance heading",
    )
    text = _replace_release_state_text(
        text,
        "The candidate AHBA, TCGA/BraTS,\n",
        "The current AHBA, TCGA/BraTS,\n",
        "reproducibility provenance description",
    )
    text = _replace_release_state_text(
        text,
        "the candidate provenance set, run\n",
        "the current provenance set, run\n",
        "reproducibility regeneration wording",
    )
    return _replace_release_state_text(
        text,
        "| BrainTrace archive | Zenodo | Current immutable v0.1.17 software/full archive: 10.5281/zenodo.22006038 / 10.5281/zenodo.22005947; v0.1.16, v0.1.15 and v0.1.14 remain immutable historical records; v0.1.12 scientific release 10.5281/zenodo.21911532; software concept 10.5281/zenodo.20773674 | Public |",
        f"| BrainTrace archive | Zenodo | Current immutable v{version} software/full archive: {manifest['software_version_doi']} / {manifest['full_repro_version_doi']}; v0.1.17, v0.1.16, v0.1.15 and v0.1.14 remain immutable historical records; v0.1.12 scientific release 10.5281/zenodo.21911532; software concept {manifest['software_concept_doi']} | Public |",
        "reproducibility archive table",
    )


def _sync_public_example_versions(root: Path, version: str) -> list[Path]:
    """Update only the version-bearing metadata in locked public fixtures."""
    paths: list[Path] = []
    for relative in PUBLIC_EXAMPLE_OUTPUTS:
        path = root / relative
        payload = json.loads(path.read_text(encoding="utf-8"))
        metadata = payload.get("meta")
        if not isinstance(metadata, dict) or "braintrace_version" not in metadata:
            raise ValueError(f"public example fixture lacks braintrace_version: {relative}")
        metadata["braintrace_version"] = version
        _write(path, json.dumps(payload, ensure_ascii=False, indent=2) + "\n")
        paths.append(path)
    return paths


def finalize(root: Path, manifest_path: Path, gate_path: Path) -> list[Path]:
    manifest = load_manifest(manifest_path)
    errors = [*validate_manifest(manifest), *_check_release_gate(gate_path)]
    if errors:
        raise ValueError("; ".join(errors))
    if manifest.get("release_state") not in {"engineering_pre_finalization", "final"}:
        raise ValueError("release manifest is not in a finalization-compatible state")

    version = str(manifest["version"])
    pyproject = root / "pyproject.toml"
    _write(pyproject, _replace_once(pyproject.read_text(encoding="utf-8"), r'^version = "[^"]+"$', f'version = "{version}"', "pyproject version"))
    cli = root / "core" / "cli.py"
    _write(cli, _replace_once(cli.read_text(encoding="utf-8"), r'^CURRENT_SOFTWARE_VERSION = "[^"]+"$', f'CURRENT_SOFTWARE_VERSION = "{version}"', "CLI version"))

    docker = root / "Dockerfile"
    docker_text = docker.read_text(encoding="utf-8")
    docker_text = re.sub(r"braintrace:v0\.1\.\d+", f"braintrace:v{version}", docker_text)
    docker_text = _replace_once(docker_text, r'LABEL org\.opencontainers\.image\.version="v[^"]+"', f'LABEL org.opencontainers.image.version="v{version}"', "application Docker version")
    _write(docker, docker_text)
    repro = root / "Dockerfile.repro"
    repro_text = repro.read_text(encoding="utf-8")
    repro_text = re.sub(r"braintrace-repro:v0\.1\.\d+", f"braintrace-repro:v{version}", repro_text)
    repro_text = _replace_once(repro_text, r"ARG BRAINTRACE_VERSION=[^\n]+", f"ARG BRAINTRACE_VERSION={version}", "repro Docker version")
    _write(repro, repro_text)

    readme = root / "README.md"
    readme_text = readme.read_text(encoding="utf-8")
    block = _release_block(manifest)
    if "<!-- BRAINTRACE:CURRENT_RELEASE:START -->" in readme_text:
        readme_text = re.sub(r"<!-- BRAINTRACE:CURRENT_RELEASE:START -->.*?<!-- BRAINTRACE:CURRENT_RELEASE:END -->\n", block, readme_text, flags=re.DOTALL)
    else:
        readme_text = _replace_once(readme_text, r"(?=## Interfaces)", block + "\n", "README current-release block")
    readme_text = re.sub(r"braintrace:v0\.1\.\d+", f"braintrace:v{version}", readme_text)
    readme_text = _finalize_readme_release_state(readme_text, manifest)
    _write(readme, readme_text)

    for relative, marker in PROVENANCE_BODY_MARKERS.items():
        path = root / relative
        text = _finalize_provenance_release_state(relative, path.read_text(encoding="utf-8"), manifest)
        index = text.find(marker)
        if index < 0:
            raise ValueError(f"could not find source-dataset section in {relative}")
        _write(path, _data_provenance_header(manifest) + "\n" + text[index:])

    _write(root / ".zenodo.json", json.dumps(_zenodo_payload(manifest), ensure_ascii=False, indent=2) + "\n")
    release_dir = root / "release" / f"v{version}"
    _write(release_dir / "RELEASE_NOTES.md", _release_notes(manifest))
    _write(release_dir / "RELEASE_CHECKLIST.md", "# Release checklist\n\n- [x] Version metadata synchronized from `release_manifest.json`.\n- [ ] Upload exact checked archives to the two reserved Zenodo drafts.\n- [ ] Verify remote file checksums and metadata before publication.\n")
    manifest["release_state"] = "final"
    _write(manifest_path, json.dumps(manifest, ensure_ascii=False, indent=2) + "\n")
    example_outputs = _sync_public_example_versions(root, version)
    return [pyproject, cli, docker, repro, readme, root / "DATA_PROVENANCE.md", root / "reproducibility" / "DATA_PROVENANCE.md", root / ".zenodo.json", release_dir / "RELEASE_NOTES.md", release_dir / "RELEASE_CHECKLIST.md", manifest_path, *example_outputs]


def check(root: Path, manifest_path: Path) -> dict[str, Any]:
    manifest = load_manifest(manifest_path)
    errors = validate_manifest(manifest)
    pyproject_authors = tomllib.loads((root / "pyproject.toml").read_text(encoding="utf-8")).get("project", {}).get("authors", [])
    author_names = [str(author.get("name", "")) for author in pyproject_authors if isinstance(author, dict)]
    if author_names != [record["name"] for record in manifest.get("creators", []) if isinstance(record, dict)]:
        errors.append("release-manifest creator order does not match pyproject authors")
    return {
        "schema": "braintrace.release_metadata_sync_audit.v1",
        "version": manifest.get("version"),
        "release_state": manifest.get("release_state"),
        "status": "PASS" if not errors else "FAIL",
        "errors": errors,
        "finalization_ready": manifest.get("release_state") == "engineering_pre_finalization" and not errors,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--finalize", action="store_true")
    parser.add_argument("--release-gate-record", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.check and args.finalize:
        parser.error("--check and --finalize are mutually exclusive")
    if args.finalize:
        if args.release_gate_record is None:
            parser.error("--finalize requires --release-gate-record")
        changed = finalize(ROOT, args.manifest.resolve(), args.release_gate_record.resolve())
        payload: dict[str, Any] = {"status": "PASS", "changed_files": [path.relative_to(ROOT).as_posix() for path in changed]}
    else:
        payload = check(ROOT, args.manifest.resolve())
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8", newline="\n")
    print(json.dumps(payload, ensure_ascii=False))
    return 0 if payload["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
