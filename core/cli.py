from __future__ import annotations

import argparse
import json
import sys
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from typing import Any, Sequence

import pandas as pd

import core.inference as locked_inference
import core.model_lock as model_lock
import core.production_route as production_route
import core.query_input as query_input


def _package_version() -> str:
    try:
        return version("braintrace")
    except PackageNotFoundError:
        return "0.1.14"


def _json_default(value: Any) -> Any:
    if hasattr(value, "item"):
        return value.item()
    if isinstance(value, Path):
        return str(value)
    raise TypeError(f"Object of type {type(value).__name__} is not JSON serializable")


def _query_payload(network_out: dict, region_out: dict, query_source: str) -> dict[str, Any]:
    parameters = dict(production_route.production_implementation_parameters())
    network_meta = network_out.get("meta", {})
    region_meta = region_out.get("meta", {})
    group_rows = region_meta.get("region_resolution_annotation", {}).get("group_ranking", [])
    lock_meta = network_meta.get("model_lock") or region_meta.get("model_lock", {})
    return {
        "network_top3": network_out.get("results", [])[:3],
        "resolution_group_top3": group_rows[:3],
        "exact_region_exploratory_top3": region_out.get("results", [])[:3],
        "meta": {
            "braintrace_version": _package_version(),
            "model_lock": lock_meta,
            "route_name": parameters["route_name"],
            "query_source": query_source,
            "reference_expression_source": region_meta.get("reference_expression_source"),
            "network_overlap": {
                "n_overlap_genes": network_meta.get("n_overlap_genes"),
                "n_model_genes": network_meta.get("n_model_genes"),
                "overlap_fraction": network_meta.get("overlap_fraction"),
            },
            "network_top3": [row.get("network_id") for row in network_out.get("results", [])[:3]],
            "resolution_group_ranking": group_rows,
            "exact_region_exploratory_ranking": region_out.get("results", []),
        },
    }


def _write_query_output(payload: dict[str, Any], output: Path, output_format: str) -> None:
    if output_format == "json":
        output.write_text(
            json.dumps(payload, indent=2, ensure_ascii=False, default=_json_default) + "\n",
            encoding="utf-8",
        )
        return

    frames = []
    for key in ("network_top3", "resolution_group_top3", "exact_region_exploratory_top3"):
        frame = pd.DataFrame(payload[key])
        if not frame.empty:
            frame.insert(0, "tier", key)
            frames.append(frame)
    pd.concat(frames, ignore_index=True, sort=False).to_csv(output, index=False)


def command_query(args: argparse.Namespace) -> int:
    expression, query_source = query_input.read_expression_file(args.input)
    network_out, region_out = locked_inference.run_locked_three_tier_route(expression)
    payload = _query_payload(network_out, region_out, query_source)
    _write_query_output(payload, args.output, args.format)
    print(f"BrainTrace query complete: {args.output}")
    return 0


def command_models(args: argparse.Namespace) -> int:
    parameters, manifest = production_route.verify_production_route()
    print(f"BrainTrace version: {_package_version()}")
    print(f"Lock ID: {manifest['lock_id']}")
    print(f"Lock status: {manifest['status']}")
    for key in (
        "route_name",
        "network_count",
        "region_count",
        "beam_count",
        "network_gene_count",
        "region_local_top_n_genes",
        "exact_top50_gene_count",
        "exact_top100_gene_count",
        "project_to_vsd",
        "allow_development_fallback",
    ):
        print(f"{key} = {str(parameters[key]).lower() if isinstance(parameters[key], bool) else parameters[key]}")
    print("Locked artifacts:")
    for relative_path in sorted(manifest["artifacts"]):
        print(f"- {relative_path}")
    return 0


def command_validate(args: argparse.Namespace) -> int:
    _, manifest = production_route.verify_production_route()
    total = len(manifest["artifacts"])
    print("BrainTrace frozen model bundle: PASS")
    print(f"Lock ID: {manifest['lock_id']}")
    print(f"Status: {manifest['status']}")
    print(f"Artifacts verified: {total}/{total}")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="braintrace",
        description="BrainTrace locked three-tier brain-origin candidate ranking",
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {_package_version()}")
    subparsers = parser.add_subparsers(dest="command", required=True)

    query = subparsers.add_parser("query", help="run one expression table through the locked production route")
    query.add_argument("--input", type=Path, required=True, help="CSV, TSV, TXT or XLSX expression table")
    query.add_argument("--output", type=Path, required=True, help="output JSON or CSV path")
    query.add_argument("--format", choices=("json", "csv"), default="json", help="output format (default: json)")
    query.set_defaults(handler=command_query)

    models = subparsers.add_parser("models", help="show the frozen production model and artifact inventory")
    models.set_defaults(handler=command_models)

    validate = subparsers.add_parser("validate", help="verify the installed frozen production model bundle")
    validate.set_defaults(handler=command_validate)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return int(args.handler(args))
    except (OSError, ValueError, model_lock.ModelLockError) as exc:
        print(f"braintrace: error: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
