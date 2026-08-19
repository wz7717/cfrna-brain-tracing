#!/usr/bin/env python
"""Strict, manifest-driven BrainTrace manuscript reproduction gate.

``python reproduce_all.py`` is the full, fail-closed manuscript gate.
``python reproduce_all.py --verify-only`` strictly verifies every required
external input. ``--profile portable`` is the repository-only CI gate.
"""

from __future__ import annotations

import argparse
import importlib.metadata
import json
import os
import platform
import re
import subprocess
import sys
import time
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Callable, Iterable


ROOT = Path(__file__).resolve().parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from core.external_inputs import (  # noqa: E402
    DEFAULT_MANIFEST as DEFAULT_EXTERNAL_MANIFEST,
    input_locator,
    load_manifest,
    resolve_alias,
    select_external_root,
    sha256_file,
    verify_inputs,
)


AUDIT_SCHEMA = "braintrace.full_reproduction_audit.v1"
DEFAULT_AUDIT_DIR = ROOT / "reproducibility_audit"
PRIVATE_PATH_RE = re.compile(
    r"(?:[A-Za-z]:\\[^\s\"']+|/Users/[^\s\"']+|/home/(?!braintrace(?:/|$))[^\s\"']+)",
    flags=re.IGNORECASE,
)


@dataclass(frozen=True)
class ReproductionStep:
    """One executable verification or regeneration step."""

    name: str
    description: str
    profiles: frozenset[str]
    command: Callable[["RunContext"], list[str]]
    display_command: Callable[["RunContext"], str]
    required_outputs: Callable[["RunContext"], list[Path]] = lambda _context: []
    timeout_seconds: int = 3600


@dataclass(frozen=True)
class RunContext:
    profile: str
    output_dir: Path
    external_root: Path | None
    release_gate: bool
    manifest_path: Path = DEFAULT_EXTERNAL_MANIFEST

    def external(self, alias: str) -> Path:
        if self.external_root is None:
            raise RuntimeError(f"{alias} requires the full external-data profile")
        return resolve_alias(
            alias,
            external_data_root=self.external_root,
            release_gate=self.release_gate,
            manifest_path=self.manifest_path,
        )


def _python(*args: str | Path) -> list[str]:
    return [sys.executable, *(str(arg) for arg in args)]


def _safe_tail(text: str, context: RunContext, *, limit: int = 4000) -> str:
    """Keep diagnostics useful without allowing private paths into an audit."""

    redacted = text
    for known_path in (ROOT, context.output_dir, context.external_root):
        if known_path is not None:
            redacted = redacted.replace(str(known_path), "<workspace-path>")
    return PRIVATE_PATH_RE.sub("<private-local-path>", redacted)[-limit:]


def _git_value(*args: str) -> str:
    try:
        result = subprocess.run(["git", *args], cwd=ROOT, text=True, capture_output=True, check=False)
    except OSError:
        return "unavailable"
    return result.stdout.strip() if result.returncode == 0 else "unavailable"


def environment_record() -> dict[str, object]:
    """Capture reproducibility state without host or user identifiers."""

    packages: dict[str, str] = {}
    for package in ("numpy", "pandas", "scipy", "scikit-learn", "matplotlib", "openpyxl", "nibabel", "pytest", "coverage"):
        try:
            packages[package] = importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:
            packages[package] = "not-installed"
    build_git_sha = os.environ.get("BRAINTRACE_GIT_SHA", "").strip()
    build_git_clean = os.environ.get("BRAINTRACE_GIT_CLEAN", "").strip().lower()
    git_status = _git_value("status", "--porcelain")
    return {
        "git_sha": build_git_sha or _git_value("rev-parse", "HEAD"),
        "git_clean": (
            True
            if build_git_clean == "true"
            else False
            if build_git_clean == "false"
            else None
            if git_status == "unavailable"
            else git_status == ""
        ),
        "python_version": platform.python_version(),
        "python_implementation": platform.python_implementation(),
        "platform_system": platform.system(),
        "dependency_versions": packages,
    }


def _output_hashes(paths: Iterable[Path], context: RunContext) -> tuple[dict[str, str], list[str]]:
    hashes: dict[str, str] = {}
    missing: list[str] = []
    for path in paths:
        try:
            logical = f"generated::{path.relative_to(context.output_dir).as_posix()}"
        except ValueError:
            logical = f"generated::{path.name}"
        if not path.exists() or not path.is_file():
            missing.append(logical)
        else:
            hashes[logical] = sha256_file(path)
    return hashes, missing


def execute_step(step: ReproductionStep, context: RunContext) -> dict[str, object]:
    """Run one command; nonzero and missing outputs both fail closed."""

    started = time.perf_counter()
    try:
        result = subprocess.run(
            step.command(context),
            cwd=ROOT,
            text=True,
            encoding="utf-8",
            errors="replace",
            capture_output=True,
            timeout=step.timeout_seconds,
            env={**os.environ, "PYTHONUTF8": "1"},
            check=False,
        )
        exit_code = result.returncode
        output = (result.stdout or "") + "\n" + (result.stderr or "")
    except subprocess.TimeoutExpired as exc:
        exit_code = 124
        stdout = exc.stdout.decode("utf-8", errors="replace") if isinstance(exc.stdout, bytes) else (exc.stdout or "")
        stderr = exc.stderr.decode("utf-8", errors="replace") if isinstance(exc.stderr, bytes) else (exc.stderr or "")
        output = stdout + "\n" + stderr + "\nSTEP_TIMEOUT"
    duration = time.perf_counter() - started
    hashes, missing = _output_hashes(step.required_outputs(context), context)
    return {
        "name": step.name,
        "description": step.description,
        "command": step.display_command(context),
        "exit_code": exit_code,
        "duration_seconds": round(duration, 6),
        "status": "PASS" if exit_code == 0 and not missing else "FAIL",
        "missing_outputs": missing,
        "generated_artifact_sha256": hashes,
        "output_tail": _safe_tail(output, context),
    }


def _portable_steps() -> list[ReproductionStep]:
    return [
        ReproductionStep(
            name="archived_benchmark_provenance",
            description="Verify the immutable archived GSE189919 benchmark manifest and its recorded origin/staged SHA-256 chain.",
            profiles=frozenset({"portable", "full"}),
            command=lambda ctx: _python(
                "scripts/verify_archived_benchmark_provenance.py",
                "--manifest", ROOT / "reproducibility" / "formal_real_input_performance_manifest.json",
                "--provenance", ROOT / "reproducibility" / "formal_real_input_performance_provenance.json",
                "--output", ctx.output_dir / "ARCHIVED_BENCHMARK_PROVENANCE_VERIFICATION.json",
            ),
            display_command=lambda _ctx: "python scripts/verify_archived_benchmark_provenance.py --manifest reproducibility/formal_real_input_performance_manifest.json --provenance reproducibility/formal_real_input_performance_provenance.json",
            required_outputs=lambda ctx: [ctx.output_dir / "ARCHIVED_BENCHMARK_PROVENANCE_VERIFICATION.json"], timeout_seconds=300,
        ),
        ReproductionStep(
            name="pytest",
            description="Run the complete test and scientific-contract suite.",
            profiles=frozenset({"portable", "full"}),
            command=lambda _ctx: _python("-m", "pytest", "tests", "-q"),
            display_command=lambda _ctx: "python -m pytest tests -q",
            timeout_seconds=1800,
        ),
        ReproductionStep(
            name="model_lock",
            description="Verify frozen model artifact hashes and lock identifier.",
            profiles=frozenset({"portable", "full"}),
            command=lambda _ctx: _python("scripts/verify_locked_artifacts.py"),
            display_command=lambda _ctx: "python scripts/verify_locked_artifacts.py",
        ),
        ReproductionStep(
            name="current_scientific_provenance",
            description="Verify current AHBA, TCGA, LOMO, tier, orthology and sign-flip provenance.",
            profiles=frozenset({"portable", "full"}),
            command=lambda _ctx: _python("scripts/verify_nonhuang_scientific_provenance.py"),
            display_command=lambda _ctx: "python scripts/verify_nonhuang_scientific_provenance.py",
        ),
        ReproductionStep(
            name="scientific_freeze",
            description="Compare all frozen canonical artifacts and current endpoint contracts to the accepted scientific baseline.",
            profiles=frozenset({"portable", "full"}),
            command=lambda ctx: _python("scripts/audit_scientific_freeze.py", "--output", ctx.output_dir / "SCIENTIFIC_FREEZE_AUDIT.json"),
            display_command=lambda _ctx: "python scripts/audit_scientific_freeze.py --output <audit>/SCIENTIFIC_FREEZE_AUDIT.json",
            required_outputs=lambda ctx: [ctx.output_dir / "SCIENTIFIC_FREEZE_AUDIT.json"],
            timeout_seconds=600,
        ),
        ReproductionStep(
            name="lomo_exact_evidence",
            description="Verify the staged LOMO Exact input and derived F1 evidence.",
            profiles=frozenset({"portable", "full"}),
            command=lambda _ctx: _python("scripts/generate_lomo_exact_f1_evidence.py", "--verify-only"),
            display_command=lambda _ctx: "python scripts/generate_lomo_exact_f1_evidence.py --verify-only",
        ),
        ReproductionStep(
            name="csv_assets",
            description="Validate tracked manuscript CSV schemas and integrity.",
            profiles=frozenset({"portable", "full"}),
            command=lambda _ctx: _python("scripts/validate_csv_assets.py"),
            display_command=lambda _ctx: "python scripts/validate_csv_assets.py",
        ),
        ReproductionStep(
            name="public_provenance_paths",
            description="Reject private machine paths from current public provenance artifacts.",
            profiles=frozenset({"portable", "full"}),
            command=lambda _ctx: _python("scripts/audit_public_provenance_paths.py"),
            display_command=lambda _ctx: "python scripts/audit_public_provenance_paths.py",
        ),
        ReproductionStep(
            name="release_state",
            description="Semantically audit current/pre-final release wording and stale current-version references.",
            profiles=frozenset({"portable", "full"}),
            command=lambda _ctx: _python("scripts/audit_release_state.py"),
            display_command=lambda _ctx: "python scripts/audit_release_state.py",
        ),
        ReproductionStep(
            name="regeneration_comparison",
            description="Regenerate deterministic evidence in a temporary location and compare it to the frozen archive.",
            profiles=frozenset({"portable", "full"}),
            command=lambda ctx: _python("scripts/verify_regeneration.py", "--output", ctx.output_dir / "REGENERATION_COMPARISON.json"),
            display_command=lambda _ctx: "python scripts/verify_regeneration.py --output <audit>/REGENERATION_COMPARISON.json",
            required_outputs=lambda ctx: [ctx.output_dir / "REGENERATION_COMPARISON.json"],
        ),
    ]


def _full_steps() -> list[ReproductionStep]:
    """Current full workflow: all raw paths are resolved via canonical aliases."""

    def build_db(ctx: RunContext) -> list[str]:
        return _python(
            "scripts/build_bo2023_atlas_from_wang_matrix.py",
            "--db", ctx.output_dir / "build" / "braintrace_source_tracing.db",
            "--matrix", ctx.external("bo2023_vsd"),
            "--counts", ctx.external("bo2023_counts"),
            "--sample-info", ctx.external("bo2023_sample_metadata"),
            "--gene-map", ctx.external("bo2023_gene_map"),
            "--outdir", ctx.output_dir / "build" / "atlas",
        )

    def initialize_db(ctx: RunContext) -> list[str]:
        return _python(
            "scripts/initialize_reproduction_database.py",
            "--db", ctx.output_dir / "build" / "braintrace_source_tracing.db",
        )

    def build_model(ctx: RunContext) -> list[str]:
        return _python(
            "scripts/build_bo2023_network_model.py",
            "--matrix", ctx.external("bo2023_vsd"),
            "--sample-info", ctx.external("bo2023_sample_metadata"),
            "--gene-map", ctx.external("bo2023_gene_map"),
            "--out", ctx.output_dir / "build" / "bo2023_saleem_network_top200_model.npz",
        )

    def verify_model(ctx: RunContext) -> list[str]:
        return _python(
            "scripts/verify_regenerated_npz.py",
            "--canonical", ROOT / "data" / "models" / "bo2023_saleem_network_top200_model.npz",
            "--regenerated", ctx.output_dir / "build" / "bo2023_saleem_network_top200_model.npz",
            "--canonical-label", "repository::data/models/bo2023_saleem_network_top200_model.npz",
            "--regenerated-label", "generated::build/bo2023_saleem_network_top200_model.npz",
            "--output", ctx.output_dir / "NETWORK_MODEL_REGENERATION_COMPARISON.json",
        )

    def build_projector(ctx: RunContext) -> list[str]:
        return _python(
            "scripts/build_bo2023_reference_projector.py",
            "--counts", ctx.external("bo2023_counts"),
            "--vsd", ctx.external("bo2023_vsd"),
            "--sample-info", ctx.external("bo2023_sample_metadata"),
            "--gene-map", ctx.external("bo2023_gene_map"),
            "--outdir", ctx.output_dir / "build" / "reference_projector",
        )

    def verify_projector(ctx: RunContext) -> list[str]:
        return _python(
            "scripts/verify_regenerated_npz.py",
            "--canonical", ROOT / "data" / "models" / "bo2023_reference_projector_linear_full.npz",
            "--regenerated", ctx.output_dir / "build" / "reference_projector" / "bo2023_reference_projector_linear_full.npz",
            "--canonical-label", "repository::data/models/bo2023_reference_projector_linear_full.npz",
            "--regenerated-label", "generated::build/reference_projector/bo2023_reference_projector_linear_full.npz",
            "--output", ctx.output_dir / "REFERENCE_PROJECTOR_REGENERATION_COMPARISON.json",
        )

    def sparse(ctx: RunContext) -> list[str]:
        return _python(
            "scripts/run_p0_4_sparse_domain_shift_sensitivity.py",
            "--counts", ctx.external("bo2023_counts"),
            "--vsd", ctx.external("bo2023_vsd"),
            "--sample-info", ctx.external("bo2023_sample_metadata"),
            "--gene-map", ctx.external("bo2023_gene_map"),
            "--replicates", "30", "--seed", "20260711", "--bootstrap-seed", "20260716",
            "--outdir", ctx.output_dir / "sparse_30_repeat",
        )

    def verify_sparse(ctx: RunContext) -> list[str]:
        return _python(
            "scripts/verify_sparse_sensitivity_run.py",
            "--methods", ctx.output_dir / "sparse_30_repeat" / "P0_4_METHODS.json",
            "--per-repeat", ctx.output_dir / "sparse_30_repeat" / "p0_8_sparse_sensitivity_per_repeat.csv",
            "--detail", ctx.output_dir / "sparse_30_repeat" / "p0_4_sparse_sensitivity_sample_detail.csv",
            "--baseline-gate", ctx.output_dir / "sparse_30_repeat" / "baseline_gate.json",
            "--output", ctx.output_dir / "SPARSE_SENSITIVITY_VERIFICATION.json",
        )

    def stage_brats(ctx: RunContext) -> list[str]:
        return _python(
            "scripts/prepare_brats_tcga_lgg_bundle.py",
            "--archive",
            ctx.external("brats_training_bundle") / "PKG - BraTS-TCGA-LGG" / "BraTS-TCGA-LGG" / "Pre-operative_TCGA_LGG_NIfTI_and_Segmentations.zip",
            "--outdir",
            ctx.output_dir / "brats_staging",
        )

    def tcga_tracing(ctx: RunContext) -> list[str]:
        return _python(
            "scripts/score_tcga_gbm_lgg_sample_tracing_with_mri_labels.py",
            "--matrix", ctx.external("tcga_expression_sample_mean"),
            "--manifest", ctx.external("tcga_gdc_manifest"),
            "--db", ctx.output_dir / "build" / "braintrace_source_tracing.db",
            "--atlas-id", "1",
            "--no-label-evaluation",
            "--outdir", ctx.output_dir / "tcga_tracing",
        )

    def tcga_mri_truth(ctx: RunContext) -> list[str]:
        return _python(
            "scripts/evaluate_brats_tcga_lgg_65_mri_truth.py",
            "--audit", ctx.output_dir / "brats_staging" / "brats_tcga_lgg_training_65_patient_audit.csv",
            "--network-predictions", ctx.output_dir / "tcga_tracing" / "tcga_gbm_lgg_sample_network_tracing.csv",
            "--atlas", ctx.external("sri24_tzo116_atlas"),
            "--atlas-lut", ctx.external("sri24_tzo116_lut"),
            "--outdir", ctx.output_dir / "tcga_brats_truth",
        )

    def verify_tcga_mri_truth(ctx: RunContext) -> list[str]:
        return _python(
            "scripts/verify_tcga_brats_truth_basis_run.py",
            "--input", ctx.output_dir / "tcga_brats_truth" / "brats_tcga_lgg_65_mri_truth_and_predictions.csv",
            "--output", ctx.output_dir / "TCGA_BRATS_TRUTH_BASIS_VERIFICATION.json",
        )

    def ahba(ctx: RunContext) -> list[str]:
        return _python(
            "scripts/run_ahba_projected_vsd_formal_three_tier_external.py",
            "--zip-dir", ctx.external("ahba_rnaseq_donor_2001").parent,
            "--bo-counts", ctx.external("bo2023_counts"),
            "--bo-vsd", ctx.external("bo2023_vsd"),
            "--sample-info", ctx.external("bo2023_sample_metadata"),
            "--gene-map", ctx.external("bo2023_gene_map"),
            "--outdir", ctx.output_dir / "ahba",
        )

    def verify_ahba(ctx: RunContext) -> list[str]:
        return _python(
            "scripts/verify_ahba_endpoint_run.py",
            "--input", ctx.output_dir / "ahba" / "ahba_formal_three_tier_sample_detail.csv",
            "--output", ctx.output_dir / "AHBA_ENDPOINT_VERIFICATION.json",
        )

    def gse_benchmark(ctx: RunContext) -> list[str]:
        return _python(
            "scripts/benchmark_real_input_inference.py",
            "--data-dir", ctx.external("gse189919_counts").parent,
            "--outdir", ctx.output_dir / "gse189919_benchmark",
            "--formal-workload", "51", "--authorize-formal",
        )

    def gse_route(ctx: RunContext) -> list[str]:
        return _python(
            "scripts/run_gse189919_latest_main_route.py",
            "--data-dir", ctx.external("gse189919_counts").parent,
            "--db-path", ctx.output_dir / "build" / "braintrace_source_tracing.db",
            "--outdir", ctx.output_dir / "gse189919_route",
        )

    def verify_gse_benchmark(ctx: RunContext) -> list[str]:
        return _python(
            "scripts/verify_gse189919_benchmark_run.py",
            "--manifest", ctx.output_dir / "gse189919_benchmark" / "manifest.json",
            "--timing", ctx.output_dir / "gse189919_benchmark" / "timing.csv",
            "--bo2023-counts", ctx.external("bo2023_counts"),
            "--output", ctx.output_dir / "GSE189919_BENCHMARK_VERIFICATION.json",
        )

    def huang(ctx: RunContext) -> list[str]:
        return _python(
            "scripts/run_huang2025_external_candidate.py",
            "--input-csv", ctx.external("huang2025_matrix_csv"),
            "--source-xlsb", ctx.external("huang2025_matrix_xlsb"),
            "--db-path", ctx.output_dir / "build" / "braintrace_source_tracing.db",
            "--outdir", ctx.output_dir / "huang2025",
        )

    def verify_huang(ctx: RunContext) -> list[str]:
        return _python(
            "scripts/verify_huang2025_run.py",
            "--input", ctx.output_dir / "huang2025" / "huang_2025_canonical_summary.json",
            "--output", ctx.output_dir / "HUANG2025_RUN_VERIFICATION.json",
        )

    return [
        ReproductionStep(
            name="reproduction_database_schema",
            description="Create a clean temporary SQLite schema without demo rows before loading verified Bo2023 data.",
            profiles=frozenset({"full"}), command=initialize_db,
            display_command=lambda _ctx: "python scripts/initialize_reproduction_database.py --db <audit>/build/braintrace_source_tracing.db",
            required_outputs=lambda ctx: [ctx.output_dir / "build" / "braintrace_source_tracing.db"], timeout_seconds=300,
        ),
        ReproductionStep(
            name="atlas_database",
            description="Rebuild the Bo2023 SQLite atlas in the audit directory.",
            profiles=frozenset({"full"}), command=build_db,
            display_command=lambda _ctx: "python scripts/build_bo2023_atlas_from_wang_matrix.py --db <audit>/build/braintrace_source_tracing.db --matrix external_source::bo2023_vsd/...",
            required_outputs=lambda ctx: [ctx.output_dir / "build" / "braintrace_source_tracing.db"], timeout_seconds=7200,
        ),
        ReproductionStep(
            name="network_model_regeneration",
            description="Regenerate the network model in the audit directory for content comparison.",
            profiles=frozenset({"full"}), command=build_model,
            display_command=lambda _ctx: "python scripts/build_bo2023_network_model.py --matrix external_source::bo2023_vsd/... --out <audit>/build/bo2023_saleem_network_top200_model.npz",
            required_outputs=lambda ctx: [ctx.output_dir / "build" / "bo2023_saleem_network_top200_model.npz"], timeout_seconds=3600,
        ),
        ReproductionStep(
            name="network_model_regeneration_comparison",
            description="Compare regenerated network-model arrays to the locked artifact without relying on ZIP serialization bytes.",
            profiles=frozenset({"full"}), command=verify_model,
            display_command=lambda _ctx: "python scripts/verify_regenerated_npz.py --canonical repository::data/models/bo2023_saleem_network_top200_model.npz --regenerated <audit>/build/bo2023_saleem_network_top200_model.npz",
            required_outputs=lambda ctx: [ctx.output_dir / "NETWORK_MODEL_REGENERATION_COMPARISON.json"], timeout_seconds=300,
        ),
        ReproductionStep(
            name="reference_projector_regeneration",
            description="Regenerate the count-to-VSD reference projector in the audit directory.",
            profiles=frozenset({"full"}), command=build_projector,
            display_command=lambda _ctx: "python scripts/build_bo2023_reference_projector.py --counts external_source::bo2023_counts/... --outdir <audit>/build/reference_projector",
            required_outputs=lambda ctx: [ctx.output_dir / "build" / "reference_projector" / "bo2023_reference_projector_linear_full.npz"], timeout_seconds=7200,
        ),
        ReproductionStep(
            name="reference_projector_regeneration_comparison",
            description="Compare regenerated reference-projector arrays to the locked artifact without relying on ZIP serialization bytes.",
            profiles=frozenset({"full"}), command=verify_projector,
            display_command=lambda _ctx: "python scripts/verify_regenerated_npz.py --canonical repository::data/models/bo2023_reference_projector_linear_full.npz --regenerated <audit>/build/reference_projector/bo2023_reference_projector_linear_full.npz",
            required_outputs=lambda ctx: [ctx.output_dir / "REFERENCE_PROJECTOR_REGENERATION_COMPARISON.json"], timeout_seconds=300,
        ),
        ReproductionStep(
            name="brats_bundle_staging",
            description="Safely extract the verified read-only BraTS bundle into the audit directory and derive its portable patient audit.",
            profiles=frozenset({"full"}), command=stage_brats,
            display_command=lambda _ctx: "python scripts/prepare_brats_tcga_lgg_bundle.py --archive external_source::brats_training_bundle/... --outdir <audit>/brats_staging",
            required_outputs=lambda ctx: [
                ctx.output_dir / "brats_staging" / "brats_tcga_lgg_training_65_patient_audit.csv",
                ctx.output_dir / "brats_staging" / "BRATS_BUNDLE_STAGING_AUDIT.json",
            ], timeout_seconds=3600,
        ),
        ReproductionStep(
            name="tcga_frozen_tracing",
            description="Regenerate the current frozen TCGA Network predictions without a legacy label-file dependency.",
            profiles=frozenset({"full"}), command=tcga_tracing,
            display_command=lambda _ctx: "python scripts/score_tcga_gbm_lgg_sample_tracing_with_mri_labels.py --no-label-evaluation --matrix external_source::tcga_expression_sample_mean/... --db <audit>/build/braintrace_source_tracing.db",
            required_outputs=lambda ctx: [ctx.output_dir / "tcga_tracing" / "tcga_gbm_lgg_sample_network_tracing.csv"], timeout_seconds=7200,
        ),
        ReproductionStep(
            name="tcga_brats_truth_basis",
            description="Derive the current MRI truth basis from staged BraTS NIfTI and the verified SRI24/TZO116+ atlas.",
            profiles=frozenset({"full"}), command=tcga_mri_truth,
            display_command=lambda _ctx: "python scripts/evaluate_brats_tcga_lgg_65_mri_truth.py --audit <audit>/brats_staging/... --network-predictions <audit>/tcga_tracing/... --atlas external_source::sri24_tzo116_atlas/...",
            required_outputs=lambda ctx: [ctx.output_dir / "tcga_brats_truth" / "brats_tcga_lgg_65_mri_truth_and_predictions.csv"], timeout_seconds=7200,
        ),
        ReproductionStep(
            name="tcga_brats_truth_basis_comparison",
            description="Compare regenerated TCGA/BraTS strict Top3 truth-basis values to the frozen scientific contract.",
            profiles=frozenset({"full"}), command=verify_tcga_mri_truth,
            display_command=lambda _ctx: "python scripts/verify_tcga_brats_truth_basis_run.py --input <audit>/tcga_brats_truth/brats_tcga_lgg_65_mri_truth_and_predictions.csv",
            required_outputs=lambda ctx: [ctx.output_dir / "TCGA_BRATS_TRUTH_BASIS_VERIFICATION.json"], timeout_seconds=300,
        ),
        ReproductionStep(
            name="sparse_sensitivity_30_repeats",
            description="Run the current frozen-seed 30-repeat sparse/domain-shift protocol.",
            profiles=frozenset({"full"}), command=sparse,
            display_command=lambda _ctx: "python scripts/run_p0_4_sparse_domain_shift_sensitivity.py --replicates 30 --seed 20260711 --bootstrap-seed 20260716 --outdir <audit>/sparse_30_repeat",
            required_outputs=lambda ctx: [
                ctx.output_dir / "sparse_30_repeat" / "p0_8_sparse_sensitivity_per_repeat.csv",
                ctx.output_dir / "sparse_30_repeat" / "P0_4_METHODS.json",
            ], timeout_seconds=21600,
        ),
        ReproductionStep(
            name="sparse_sensitivity_30_repeats_comparison",
            description="Verify all four nonbaseline sparse scenarios retain 30 fixed-seed independent repeats and the locked baseline endpoint.",
            profiles=frozenset({"full"}), command=verify_sparse,
            display_command=lambda _ctx: "python scripts/verify_sparse_sensitivity_run.py --methods <audit>/sparse_30_repeat/P0_4_METHODS.json --per-repeat <audit>/sparse_30_repeat/p0_8_sparse_sensitivity_per_repeat.csv",
            required_outputs=lambda ctx: [ctx.output_dir / "SPARSE_SENSITIVITY_VERIFICATION.json"], timeout_seconds=300,
        ),
        ReproductionStep(
            name="ahba_endpoint_evaluability",
            description="Run the current AHBA projected-VSD three-tier endpoint workflow.",
            profiles=frozenset({"full"}), command=ahba,
            display_command=lambda _ctx: "python scripts/run_ahba_projected_vsd_formal_three_tier_external.py --zip-dir external_source::ahba_rnaseq_donor_2001/... --outdir <audit>/ahba",
            required_outputs=lambda ctx: [
                ctx.output_dir / "ahba" / "ahba_formal_three_tier_sample_detail.csv",
                ctx.output_dir / "ahba" / "ahba_formal_three_tier_summary.json",
            ], timeout_seconds=21600,
        ),
        ReproductionStep(
            name="ahba_endpoint_comparison",
            description="Compare the regenerated AHBA endpoint counts to the frozen 231/223/88 scientific contract.",
            profiles=frozenset({"full"}), command=verify_ahba,
            display_command=lambda _ctx: "python scripts/verify_ahba_endpoint_run.py --input <audit>/ahba/ahba_formal_three_tier_sample_detail.csv",
            required_outputs=lambda ctx: [ctx.output_dir / "AHBA_ENDPOINT_VERIFICATION.json"], timeout_seconds=300,
        ),
        ReproductionStep(
            name="gse189919_benchmark",
            description="Run the frozen 51-profile, 28,415-gene benchmark with three warm repeats.",
            profiles=frozenset({"full"}), command=gse_benchmark,
            display_command=lambda _ctx: "python scripts/benchmark_real_input_inference.py --data-dir external_source::gse189919_counts/... --formal-workload 51 --authorize-formal",
            required_outputs=lambda ctx: [ctx.output_dir / "gse189919_benchmark" / "manifest.json"], timeout_seconds=3600,
        ),
        ReproductionStep(
            name="gse189919_frozen_route",
            description="Run the formal GSE189919 frozen no-pairwise route against the regenerated atlas.",
            profiles=frozenset({"full"}), command=gse_route,
            display_command=lambda _ctx: "python scripts/run_gse189919_latest_main_route.py --data-dir external_source::gse189919_counts/... --db-path <audit>/build/braintrace_source_tracing.db",
            required_outputs=lambda ctx: [ctx.output_dir / "gse189919_route" / "manifest.json"], timeout_seconds=7200,
        ),
        ReproductionStep(
            name="gse189919_benchmark_contract",
            description="Verify the new GSE189919 51-profile/three-warm-repeat benchmark schema without equating hardware-dependent timings.",
            profiles=frozenset({"full"}), command=verify_gse_benchmark,
            display_command=lambda _ctx: "python scripts/verify_gse189919_benchmark_run.py --manifest <audit>/gse189919_benchmark/manifest.json --timing <audit>/gse189919_benchmark/timing.csv --bo2023-counts external_source::bo2023_counts/...",
            required_outputs=lambda ctx: [ctx.output_dir / "GSE189919_BENCHMARK_VERIFICATION.json"], timeout_seconds=300,
        ),
        ReproductionStep(
            name="huang2025_external_audit",
            description="Run the current Huang 2025 provenance-remediated external-domain audit.",
            profiles=frozenset({"full"}), command=huang,
            display_command=lambda _ctx: "python scripts/run_huang2025_external_candidate.py --input-csv external_source::huang2025_matrix_csv/... --source-xlsb external_source::huang2025_matrix_xlsb/...",
            required_outputs=lambda ctx: [ctx.output_dir / "huang2025" / "huang_2025_canonical_summary.json"], timeout_seconds=7200,
        ),
        ReproductionStep(
            name="huang2025_external_comparison",
            description="Verify the frozen Huang2025 profile counts and non-pairing provenance guardrails.",
            profiles=frozenset({"full"}), command=verify_huang,
            display_command=lambda _ctx: "python scripts/verify_huang2025_run.py --input <audit>/huang2025/huang_2025_canonical_summary.json",
            required_outputs=lambda ctx: [ctx.output_dir / "HUANG2025_RUN_VERIFICATION.json"], timeout_seconds=300,
        ),
    ]


def _write_audit(path: Path, payload: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8", newline="\n")


def verify_raw_data(
    external_data_root: Path | str | None = None,
    *,
    profile: str = "full",
    release_gate: bool | None = None,
    manifest_path: Path = DEFAULT_EXTERNAL_MANIFEST,
) -> dict[str, object]:
    """Backward-compatible strict raw-input verification API."""

    strict_release = profile == "full" if release_gate is None else release_gate
    return verify_inputs(
        profile=profile,
        external_data_root=external_data_root,
        release_gate=strict_release,
        manifest_path=manifest_path,
    )


def _report_input_hashes(input_report: dict[str, object]) -> dict[str, str]:
    hashes: dict[str, str] = {}
    for item in input_report.get("inputs", []):
        if isinstance(item, dict):
            digest = item.get("sha256") or item.get("tree_sha256") or item.get("observed_sha256")
            if isinstance(digest, str) and len(digest) == 64:
                hashes[str(item.get("alias", "unknown"))] = digest
    return hashes


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--profile", choices=("full", "portable"), default="full")
    parser.add_argument("--verify-only", action="store_true", help="Strictly verify the selected input profile and write an audit.")
    parser.add_argument("--external-data-root", type=Path, help="Canonical external_data root; takes highest priority.")
    parser.add_argument(
        "--input-manifest",
        type=Path,
        default=DEFAULT_EXTERNAL_MANIFEST,
        help="External-input manifest to verify (defaults to the tracked release contract).",
    )
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_AUDIT_DIR, help="Directory for FULL_REPRODUCTION_AUDIT.json and regenerated outputs.")
    parser.add_argument("--step", action="append", help="Run only named step(s), after strict input verification.")
    parser.add_argument("--allow-legacy-fallback", action="store_true", help="Development-only compatibility switch; prohibited for full reproduction.")
    parser.add_argument("--release-gate", action="store_true", help="Reject all legacy input fallbacks; implied by the full profile.")
    args = parser.parse_args(argv)
    if args.profile == "full" and args.allow_legacy_fallback:
        parser.error("the full manuscript gate requires canonical external_data; legacy fallback is not permitted")

    output_dir = args.output_dir.resolve()
    manifest_path = args.input_manifest.resolve()
    release_gate = bool(args.release_gate or args.profile == "full")
    selection = select_external_root(args.external_data_root)
    context = RunContext(
        args.profile,
        output_dir,
        selection.path if args.profile == "full" else None,
        release_gate,
        manifest_path,
    )
    input_report = verify_raw_data(
        args.external_data_root,
        profile=args.profile,
        release_gate=release_gate,
        manifest_path=manifest_path,
    )
    audit: dict[str, object] = {
        "schema": AUDIT_SCHEMA,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "profile": args.profile,
        "release_gate": release_gate,
        "environment": environment_record(),
        "canonical_external_root": "external_data/",
        "external_root_selection": selection.source,
        "input_contract": input_report,
        "canonical_source_aliases": [
            {"alias": item["alias"], "locator": input_locator(item)}
            for item in load_manifest(manifest_path)["inputs"] if args.profile in item["profiles"]
        ],
        "input_sha256": _report_input_hashes(input_report),
        "steps": [],
        "scientific_freeze": {"status": "NOT_RUN"},
        "status": "FAIL",
    }
    output_dir.mkdir(parents=True, exist_ok=True)
    audit_path = output_dir / "FULL_REPRODUCTION_AUDIT.json"

    for warning in input_report.get("warnings", []):
        if isinstance(warning, str):
            print(warning, file=sys.stderr)

    if input_report.get("status") != "PASS":
        audit["scientific_freeze"] = {"status": "NOT_RUN", "reason": "input contract failed closed"}
        _write_audit(audit_path, audit)
        print(json.dumps({"status": "FAIL", "reason": "input_contract", "audit": "FULL_REPRODUCTION_AUDIT.json"}))
        return 2
    if args.verify_only:
        audit["scientific_freeze"] = {"status": "NOT_RUN", "reason": "verify-only"}
        audit["status"] = "PASS"
        _write_audit(audit_path, audit)
        print(json.dumps({"status": "PASS", "audit": "FULL_REPRODUCTION_AUDIT.json"}))
        return 0

    steps = [step for step in [*_portable_steps(), *_full_steps()] if args.profile in step.profiles]
    if args.step:
        requested = set(args.step)
        unknown = requested - {step.name for step in steps}
        if unknown:
            parser.error(f"unknown step(s): {', '.join(sorted(unknown))}")
        steps = [step for step in steps if step.name in requested]
    records = [execute_step(step, context) for step in steps]
    audit["steps"] = records
    scientific_names = {
        "archived_benchmark_provenance",
        "model_lock",
        "network_model_regeneration_comparison",
        "reference_projector_regeneration_comparison",
        "ahba_endpoint_comparison",
        "huang2025_external_comparison",
        "current_scientific_provenance",
        "scientific_freeze",
        "lomo_exact_evidence",
        "regeneration_comparison",
        "tcga_brats_truth_basis_comparison",
        "sparse_sensitivity_30_repeats_comparison",
    }
    scientific_records = [record for record in records if record["name"] in scientific_names]
    audit["scientific_freeze"] = {
        "status": "PASS" if scientific_records and all(record["status"] == "PASS" for record in scientific_records) else "FAIL",
        "checks": [{"name": record["name"], "status": record["status"]} for record in scientific_records],
    }
    audit["status"] = "PASS" if records and all(record["status"] == "PASS" for record in records) else "FAIL"
    _write_audit(audit_path, audit)
    print(json.dumps({"status": audit["status"], "audit": "FULL_REPRODUCTION_AUDIT.json"}))
    return 0 if audit["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
