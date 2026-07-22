#!/usr/bin/env python
"""Smoke-test real-input inference performance without emitting clinical metadata.

This harness is deliberately restricted to the approved three-sample smoke gate.
The registered 1/8/51 workloads are metadata only until a formal benchmark run is
separately authorized.
"""

from __future__ import annotations

import argparse
import csv
import ctypes
import hashlib
import importlib.metadata
import json
import os
import platform
import subprocess
import sys
import time
import tracemalloc
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

import pandas as pd


PROCESS_ENTRY_NS = time.perf_counter_ns()
ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from core import bo2023_region_tracing  # noqa: E402
from core.bo2023_region_tracing import ROUTE_NAME, trace_bo2023_secondary_regions  # noqa: E402
from core.network_tracing import trace_network_expression  # noqa: E402


COUNT_FILE_NAME = "GSE189919_count.csv.gz"
SMOKE_SAMPLE_COUNT = 3
FORMAL_WORKLOADS = [1, 8, 51]
FORMAL_WARM_REPEATS = 3
CONCURRENCY = 1
PACKAGED_REFERENCE_SOURCE = "packaged_region_logcpm_reference"


@dataclass(frozen=True)
class ProcessMemory:
    working_set_bytes: int | None
    peak_working_set_bytes: int | None


class _PROCESS_MEMORY_COUNTERS(ctypes.Structure):
    _fields_ = [
        ("cb", ctypes.c_ulong),
        ("PageFaultCount", ctypes.c_ulong),
        ("PeakWorkingSetSize", ctypes.c_size_t),
        ("WorkingSetSize", ctypes.c_size_t),
        ("QuotaPeakPagedPoolUsage", ctypes.c_size_t),
        ("QuotaPagedPoolUsage", ctypes.c_size_t),
        ("QuotaPeakNonPagedPoolUsage", ctypes.c_size_t),
        ("QuotaNonPagedPoolUsage", ctypes.c_size_t),
        ("PagefileUsage", ctypes.c_size_t),
        ("PeakPagefileUsage", ctypes.c_size_t),
    ]


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def anonymize_sample_id(sample_id: str) -> str:
    return hashlib.sha256(str(sample_id).encode("utf-8")).hexdigest()[:16]


def process_memory() -> ProcessMemory:
    if os.name != "nt":
        return ProcessMemory(None, None)
    kernel32 = ctypes.WinDLL("kernel32", use_last_error=True)
    psapi = ctypes.WinDLL("psapi", use_last_error=True)
    kernel32.GetCurrentProcess.argtypes = []
    kernel32.GetCurrentProcess.restype = ctypes.c_void_p
    psapi.GetProcessMemoryInfo.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_ulong]
    psapi.GetProcessMemoryInfo.restype = ctypes.c_int
    counters = _PROCESS_MEMORY_COUNTERS()
    counters.cb = ctypes.sizeof(counters)
    ok = psapi.GetProcessMemoryInfo(kernel32.GetCurrentProcess(), ctypes.byref(counters), counters.cb)
    if not ok:
        raise ctypes.WinError()
    return ProcessMemory(int(counters.WorkingSetSize), int(counters.PeakWorkingSetSize))


def total_physical_memory_bytes() -> int | None:
    if os.name != "nt":
        return None

    class MEMORYSTATUSEX(ctypes.Structure):
        _fields_ = [
            ("dwLength", ctypes.c_ulong),
            ("dwMemoryLoad", ctypes.c_ulong),
            ("ullTotalPhys", ctypes.c_ulonglong),
            ("ullAvailPhys", ctypes.c_ulonglong),
            ("ullTotalPageFile", ctypes.c_ulonglong),
            ("ullAvailPageFile", ctypes.c_ulonglong),
            ("ullTotalVirtual", ctypes.c_ulonglong),
            ("ullAvailVirtual", ctypes.c_ulonglong),
            ("ullAvailExtendedVirtual", ctypes.c_ulonglong),
        ]

    status = MEMORYSTATUSEX()
    status.dwLength = ctypes.sizeof(status)
    kernel32 = ctypes.WinDLL("kernel32", use_last_error=True)
    kernel32.GlobalMemoryStatusEx.argtypes = [ctypes.c_void_p]
    kernel32.GlobalMemoryStatusEx.restype = ctypes.c_int
    if not kernel32.GlobalMemoryStatusEx(ctypes.byref(status)):
        raise ctypes.WinError()
    return int(status.ullTotalPhys)


def git_commit() -> str:
    result = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=ROOT,
        check=True,
        capture_output=True,
        text=True,
    )
    return result.stdout.strip()


def software_versions() -> dict[str, str]:
    packages = ["numpy", "pandas", "scipy", "scikit-learn", "streamlit"]
    versions: dict[str, str] = {}
    for package in packages:
        try:
            versions[package] = importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:
            versions[package] = "not-installed"
    return versions


def environment_manifest() -> dict[str, Any]:
    return {
        "os": platform.platform(),
        "machine": platform.machine(),
        "processor": platform.processor(),
        "logical_cpu_count": os.cpu_count(),
        "total_physical_memory_bytes": total_physical_memory_bytes(),
        "python": sys.version.replace("\n", " "),
        "packages": software_versions(),
        "git_commit": git_commit(),
    }


def read_counts_subset(count_path: Path, sample_count: int) -> tuple[pd.DataFrame, list[str], int]:
    header = pd.read_csv(count_path, nrows=0)
    if "Geneid" not in header.columns:
        raise ValueError("GSE189919 count matrix is missing Geneid")
    all_sample_ids = [str(column) for column in header.columns if str(column) != "Geneid"]
    if len(all_sample_ids) < sample_count:
        raise ValueError(f"requested {sample_count} samples, found {len(all_sample_ids)}")
    selected = all_sample_ids[:sample_count]
    frame = pd.read_csv(count_path, usecols=["Geneid", *selected]).set_index("Geneid")
    if frame.index.duplicated().any():
        raise ValueError("GSE189919 count matrix contains duplicate gene identifiers")
    frame = frame.apply(pd.to_numeric, errors="raise")
    if frame.isna().any().any() or (frame < 0).any().any():
        raise ValueError("GSE189919 count matrix contains missing or negative values")
    return frame.astype(float), selected, len(all_sample_ids)


def expression_for_sample(counts: pd.DataFrame, sample_id: str) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "gene_symbol": counts.index.astype(str),
            "read_count": counts[sample_id].to_numpy(dtype=float),
        }
    )


def assert_canonical_contract() -> dict[str, Any]:
    matrix, region_network, source = bo2023_region_tracing._load_packaged_region_reference_matrix()  # noqa: SLF001
    beams = bo2023_region_tracing._load_formal_beam_gene_panels()  # noqa: SLF001
    bo2023_region_tracing._validate_formal_beam_gene_panels(beams, matrix, region_network)  # noqa: SLF001
    networks = sorted(set(region_network.values()))
    contract = {
        "canonical_regions": int(matrix.shape[1]),
        "networks": int(len(networks)),
        "network_top3_beams": int(len(beams)),
        "reference_genes": int(matrix.shape[0]),
        "reference_source": source,
    }
    expected = {
        "canonical_regions": 110,
        "networks": 10,
        "network_top3_beams": 120,
        "reference_genes": 21668,
        "reference_source": PACKAGED_REFERENCE_SOURCE,
    }
    for key, value in expected.items():
        if contract[key] != value:
            raise AssertionError(f"canonical contract mismatch for {key}: {contract[key]} != {value}")
    return contract


def run_frozen_route(expression: pd.DataFrame) -> tuple[dict[str, Any], dict[str, Any]]:
    network_output = trace_network_expression(
        expression,
        min_overlap_fraction=0.50,
        project_to_vsd=True,
        enable_pairwise_rescue=False,
    )
    rescue = network_output.get("meta", {}).get("pairwise_rescue", {})
    if rescue.get("enabled") or rescue.get("switched"):
        raise AssertionError(f"pairwise rescue must remain disabled: {rescue}")
    if len(network_output.get("results", [])) != 10:
        raise AssertionError("frozen Network route did not return exactly 10 candidates")
    region_output = trace_bo2023_secondary_regions(
        expression,
        network_output,
        "__packaged_reference_only__",
        None,
        topk=30,
    )
    if not region_output.get("results"):
        raise AssertionError(f"frozen region route returned no candidates: {region_output.get('meta', {})}")
    if len(region_output.get("meta", {}).get("network_beam", [])) != 3:
        raise AssertionError("frozen region route did not retain a three-Network beam")
    return network_output, region_output


def _timed_workload(
    counts: pd.DataFrame,
    sample_ids: list[str],
    *,
    phase: str,
    repeat: int,
) -> tuple[list[dict[str, Any]], float]:
    rows: list[dict[str, Any]] = []
    wall_started = time.perf_counter_ns()
    rescue_switches = 0
    for index, sample_id in enumerate(sample_ids, start=1):
        expression = expression_for_sample(counts, sample_id)
        sample_started = time.perf_counter_ns()
        network_output, region_output = run_frozen_route(expression)
        elapsed_seconds = (time.perf_counter_ns() - sample_started) / 1_000_000_000.0
        rescue = network_output.get("meta", {}).get("pairwise_rescue", {})
        rescue_switches += int(bool(rescue.get("switched")))
        rows.append(
            {
                "phase": phase,
                "repeat": repeat,
                "sample_index": index,
                "sample_hash": anonymize_sample_id(sample_id),
                "elapsed_seconds": elapsed_seconds,
                "network_candidates": len(network_output.get("results", [])),
                "region_candidates": len(region_output.get("results", [])),
                "pairwise_rescue_switched": bool(rescue.get("switched")),
            }
        )
    wall_seconds = (time.perf_counter_ns() - wall_started) / 1_000_000_000.0
    if rescue_switches:
        raise AssertionError("pairwise rescue switched during frozen-route benchmark")
    return rows, wall_seconds


def _timing_summary(rows: list[dict[str, Any]], wall_seconds: float) -> dict[str, Any]:
    elapsed = pd.Series([float(row["elapsed_seconds"]) for row in rows], dtype=float)
    return {
        "samples": len(rows),
        "wall_total_seconds": wall_seconds,
        "wall_seconds_per_sample": wall_seconds / len(rows),
        "sample_time_sum_seconds": float(elapsed.sum()),
        "sample_time_p50_seconds": float(elapsed.quantile(0.50)),
        "sample_time_p95_seconds": float(elapsed.quantile(0.95)),
        "pairwise_rescue_switches": int(
            sum(bool(row["pairwise_rescue_switched"]) for row in rows)
        ),
    }


def _memory_snapshot() -> dict[str, Any]:
    current_py, peak_py = tracemalloc.get_traced_memory()
    return {
        "process": asdict(process_memory()),
        "python_current_bytes": current_py,
        "python_peak_bytes": peak_py,
    }


def run_formal_workload(
    data_dir: Path,
    outdir: Path,
    workload: int,
    *,
    authorized: bool,
) -> dict[str, Any]:
    """Run one explicitly authorized formal workload in its own process.

    Each command invocation measures one workload. The cold measure starts at
    module entry (after the Python interpreter has started), includes remaining
    imports, input hashing/loading, contract validation, and one full workload.
    It is followed by one untimed single-sample warmup and three timed repeats.
    """
    if not authorized:
        raise PermissionError("formal benchmark requires explicit --authorize-formal")
    if workload not in FORMAL_WORKLOADS:
        raise ValueError(f"workload must be one of {FORMAL_WORKLOADS}")
    if not data_dir.is_absolute():
        raise ValueError("--data-dir must be an absolute read-only input path")
    count_path = data_dir / COUNT_FILE_NAME
    if not count_path.is_file():
        raise FileNotFoundError(count_path)
    outdir.mkdir(parents=True, exist_ok=False)

    memory_start = process_memory()
    tracemalloc.start()
    input_sha256 = sha256_file(count_path)
    counts, sample_ids, total_input_samples = read_counts_subset(count_path, workload)
    contract = assert_canonical_contract()

    cold_rows, cold_inference_wall = _timed_workload(
        counts,
        sample_ids,
        phase="cold",
        repeat=0,
    )
    cold_process_entry_seconds = (time.perf_counter_ns() - PROCESS_ENTRY_NS) / 1_000_000_000.0
    cold_memory = _memory_snapshot()

    # Untimed warmup is deliberately excluded from all warm summaries.
    run_frozen_route(expression_for_sample(counts, sample_ids[0]))

    timing_rows = list(cold_rows)
    warm_repeats: list[dict[str, Any]] = []
    for repeat in range(1, FORMAL_WARM_REPEATS + 1):
        tracemalloc.reset_peak()
        rows, wall_seconds = _timed_workload(
            counts,
            sample_ids,
            phase="warm",
            repeat=repeat,
        )
        timing_rows.extend(rows)
        warm_repeats.append(
            {
                "repeat": repeat,
                "timing": _timing_summary(rows, wall_seconds),
                "memory": _memory_snapshot(),
            }
        )
    tracemalloc.stop()

    timing_path = outdir / "timing.csv"
    with timing_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(timing_rows[0]))
        writer.writeheader()
        writer.writerows(timing_rows)

    warm_rows = [row for row in timing_rows if row["phase"] == "warm"]
    warm_wall_total = sum(
        float(item["timing"]["wall_total_seconds"]) for item in warm_repeats
    )
    manifest = {
        "schema_version": 2,
        "status": "formal_real_input_performance_gate",
        "interpretation_boundary": (
            "Engineering performance on public GSE189919 count inputs only; "
            "this benchmark makes no accuracy, localization, or clinical claim."
        ),
        "route": ROUTE_NAME,
        "input": {
            "file_name": COUNT_FILE_NAME,
            "sha256": input_sha256,
            "bytes": count_path.stat().st_size,
            "total_samples_in_header": total_input_samples,
            "workload_samples": workload,
            "fixed_header_order": True,
            "clinical_metadata_read": False,
            "clinical_metadata_copied": False,
            "persistent_database_read": False,
        },
        "preregistration": {
            "formal_workload_sizes": FORMAL_WORKLOADS,
            "formal_warm_repeats": FORMAL_WARM_REPEATS,
            "formal_execution_authorized": True,
            "concurrency": CONCURRENCY,
            "independent_process_per_workload": True,
            "untimed_warmup_samples": 1,
        },
        "measurement_definitions": {
            "cold_process_entry_to_completion_seconds": (
                "Elapsed from this Python module's timer initialization through remaining "
                "imports, input SHA/read, canonical contract validation, and one complete "
                "workload; excludes interpreter launch before module execution."
            ),
            "cold_inference_wall_seconds": (
                "Wall time for the first complete workload after input and contract setup."
            ),
            "warm_wall_total_seconds": (
                "Wall time around each complete workload after one untimed warmup."
            ),
            "sample_time_p50_p95_seconds": (
                "Nearest-linear pandas quantiles of per-sample frozen-route inference "
                "durations; aggregate warm quantiles pool all three fixed-order repeats."
            ),
            "memory": (
                "Windows PeakWorkingSetSize is the process-lifetime peak and includes native "
                "allocations; tracemalloc records Python allocations and can exclude NumPy native memory."
            ),
        },
        "cold": {
            "process_entry_to_completion_seconds": cold_process_entry_seconds,
            "timing": _timing_summary(cold_rows, cold_inference_wall),
            "memory": cold_memory,
        },
        "warm": {
            "repeats": warm_repeats,
            "aggregate": _timing_summary(warm_rows, warm_wall_total),
        },
        "memory_at_process_start": asdict(memory_start),
        "canonical_contract": contract,
        "environment": environment_manifest(),
        "outputs": {"timing_csv": timing_path.name},
    }
    manifest_path = outdir / "manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )

    # Output artifacts must not expose raw sample identifiers.
    output_text = timing_path.read_text(encoding="utf-8") + manifest_path.read_text(encoding="utf-8")
    leaked = [sample_id for sample_id in sample_ids if sample_id in output_text]
    if leaked:
        raise AssertionError("raw sample identifiers leaked into benchmark artifacts")
    return manifest


def run_smoke(data_dir: Path, outdir: Path) -> dict[str, Any]:
    if not data_dir.is_absolute():
        raise ValueError("--data-dir must be an absolute read-only input path")
    count_path = data_dir / COUNT_FILE_NAME
    if not count_path.is_file():
        raise FileNotFoundError(count_path)
    outdir.mkdir(parents=True, exist_ok=True)

    memory_start = process_memory()
    tracemalloc.start()
    cold_started = time.perf_counter_ns()
    input_sha256 = sha256_file(count_path)
    counts, sample_ids, total_input_samples = read_counts_subset(count_path, SMOKE_SAMPLE_COUNT)
    contract = assert_canonical_contract()
    warmup_expression = expression_for_sample(counts, sample_ids[0])
    warmup_network, _ = run_frozen_route(warmup_expression)
    warmup_rescue = warmup_network.get("meta", {}).get("pairwise_rescue", {})
    cold_seconds = (time.perf_counter_ns() - cold_started) / 1_000_000_000.0
    cold_current_py, cold_peak_py = tracemalloc.get_traced_memory()
    cold_memory = process_memory()

    tracemalloc.reset_peak()
    timing_rows: list[dict[str, Any]] = []
    warm_started = time.perf_counter_ns()
    rescue_switches = 0
    for index, sample_id in enumerate(sample_ids, start=1):
        expression = expression_for_sample(counts, sample_id)
        sample_started = time.perf_counter_ns()
        network_output, region_output = run_frozen_route(expression)
        elapsed_seconds = (time.perf_counter_ns() - sample_started) / 1_000_000_000.0
        rescue = network_output.get("meta", {}).get("pairwise_rescue", {})
        rescue_switches += int(bool(rescue.get("switched")))
        timing_rows.append(
            {
                "sample_index": index,
                "sample_hash": anonymize_sample_id(sample_id),
                "elapsed_seconds": elapsed_seconds,
                "network_candidates": len(network_output.get("results", [])),
                "region_candidates": len(region_output.get("results", [])),
                "pairwise_rescue_switched": bool(rescue.get("switched")),
            }
        )
    warm_seconds = (time.perf_counter_ns() - warm_started) / 1_000_000_000.0
    warm_current_py, warm_peak_py = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    warm_memory = process_memory()
    if rescue_switches != 0 or warmup_rescue.get("switched"):
        raise AssertionError("pairwise rescue switched during frozen-route smoke")

    timing_path = outdir / "timing.csv"
    with timing_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(timing_rows[0]))
        writer.writeheader()
        writer.writerows(timing_rows)

    manifest = {
        "schema_version": 1,
        "status": "three-sample-smoke-only_not_formal_benchmark",
        "route": ROUTE_NAME,
        "input": {
            "file_name": COUNT_FILE_NAME,
            "sha256": input_sha256,
            "bytes": count_path.stat().st_size,
            "total_samples_in_header": total_input_samples,
            "smoke_samples": SMOKE_SAMPLE_COUNT,
            "clinical_metadata_copied": False,
        },
        "preregistration": {
            "formal_workload_sizes": FORMAL_WORKLOADS,
            "formal_plan": "1 cold run plus 3 warm repeats at each workload",
            "formal_execution_authorized": False,
            "concurrency": CONCURRENCY,
            "fixed_input_order": True,
            "warmup_samples_excluded_from_timing": 1,
        },
        "smoke": {
            "cold_startup_seconds_including_input_contract_and_untimed_warmup": cold_seconds,
            "warm_timed_samples": len(timing_rows),
            "warm_total_seconds": warm_seconds,
            "warm_seconds_per_sample": warm_seconds / len(timing_rows),
            "pairwise_rescue_switches": rescue_switches,
            "memory": {
                "process_start": asdict(memory_start),
                "after_cold_startup": asdict(cold_memory),
                "after_warm_workload": asdict(warm_memory),
                "cold_python_current_bytes": cold_current_py,
                "cold_python_peak_bytes": cold_peak_py,
                "warm_python_current_bytes": warm_current_py,
                "warm_python_peak_bytes": warm_peak_py,
                "note": "PeakWorkingSetSize is the Windows process lifetime peak; tracemalloc covers Python allocations and may exclude native NumPy allocations.",
            },
        },
        "canonical_contract": contract,
        "environment": environment_manifest(),
        "outputs": {"timing_csv": timing_path.name},
    }
    manifest_path = outdir / "manifest.json"
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return manifest


def main() -> int:
    parser = argparse.ArgumentParser(description="Run real-input performance gates.")
    parser.add_argument("--data-dir", type=Path, required=True, help="Absolute read-only GSE189919 input directory.")
    parser.add_argument("--outdir", type=Path, required=True, help="Directory for anonymized smoke evidence.")
    parser.add_argument(
        "--formal-workload",
        type=int,
        choices=FORMAL_WORKLOADS,
        help="Run one registered formal workload in this independent process.",
    )
    parser.add_argument(
        "--authorize-formal",
        action="store_true",
        help="Explicitly acknowledge coordinator authorization for a formal workload.",
    )
    args = parser.parse_args()
    if args.formal_workload is None:
        if args.authorize_formal:
            parser.error("--authorize-formal requires --formal-workload")
        manifest = run_smoke(args.data_dir, args.outdir)
    else:
        manifest = run_formal_workload(
            args.data_dir,
            args.outdir,
            args.formal_workload,
            authorized=args.authorize_formal,
        )
    print(json.dumps({"status": manifest["status"], "outdir": str(args.outdir.resolve())}, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
