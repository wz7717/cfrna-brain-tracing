from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

from benchmark_runner import run_paper_grade_benchmark_suite, default_label_extractor
from core.methods import method_choices
from reporting import (
    build_benchmark_report_bundle_bytes,
    export_benchmark_paper_figures,
    export_benchmark_report_pdf,
    export_run_bundle,
)
from data.bo2023_buildkit import import_buildkit_dir
from data.bo2023_region_matrix import import_region_matrix
from signature_builder import build_signature_set


def _make_label_extractor(label_key: str | None):
    if not label_key:
        return default_label_extractor

    def extractor(row):
        meta = row.get('metadata')
        if not meta:
            return None
        try:
            obj = json.loads(meta)
        except Exception:
            return None
        return obj.get(label_key)

    return extractor


def cmd_benchmark(args: argparse.Namespace) -> int:
    suite = run_paper_grade_benchmark_suite(
        db_path=args.db,
        method=args.method,
        k=args.topk,
        atlas_id=args.atlas_id,
        sigset_id=args.sigset_id,
        use_value=args.use_value,
        l2=args.l2,
        ensemble_alpha=args.ensemble_alpha,
        bootstrap_n=args.bootstrap_n,
        bootstrap_gene_frac=args.bootstrap_gene_frac,
        label_extractor=_make_label_extractor(args.label_key),
        limit=args.limit,
    )
    outdir = Path(args.output_dir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    meta = {
        'db': str(Path(args.db).resolve()),
        'method': args.method,
        'topk': args.topk,
        'atlas_id': args.atlas_id,
        'sigset_id': args.sigset_id,
        'use_value': args.use_value,
        'l2': args.l2,
        'ensemble_alpha': args.ensemble_alpha,
        'bootstrap_n': args.bootstrap_n,
        'bootstrap_gene_frac': args.bootstrap_gene_frac,
        'label_key': args.label_key,
    }
    export_benchmark_paper_figures(outdir, suite=suite, metadata=meta, prefix=args.prefix)
    report_pdf = outdir / f'{args.prefix}_benchmark_report.pdf'
    export_benchmark_report_pdf(report_pdf, suite=suite, metadata=meta)

    metrics_path = outdir / 'tables' / 'benchmark_metrics.csv'
    detail_path = outdir / 'tables' / 'benchmark_detail.csv'
    print(f'Benchmark export completed: {outdir}')
    if metrics_path.exists():
        print(f'Metrics: {metrics_path}')
    if detail_path.exists():
        print(f'Detail: {detail_path}')
    print(f'Report PDF: {report_pdf}')
    return 0


def cmd_benchmark_bundle(args: argparse.Namespace) -> int:
    suite = run_paper_grade_benchmark_suite(
        db_path=args.db,
        method=args.method,
        k=args.topk,
        atlas_id=args.atlas_id,
        sigset_id=args.sigset_id,
        use_value=args.use_value,
        l2=args.l2,
        ensemble_alpha=args.ensemble_alpha,
        bootstrap_n=args.bootstrap_n,
        bootstrap_gene_frac=args.bootstrap_gene_frac,
        label_extractor=_make_label_extractor(args.label_key),
        limit=args.limit,
    )
    payload = build_benchmark_report_bundle_bytes(
        suite=suite,
        metadata={
            'db': str(Path(args.db).resolve()),
            'method': args.method,
            'topk': args.topk,
            'atlas_id': args.atlas_id,
            'sigset_id': args.sigset_id,
            'use_value': args.use_value,
        },
        prefix=args.prefix,
    )
    out = Path(args.output_zip).resolve()
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_bytes(payload)
    print(f'Benchmark bundle written to: {out}')
    return 0




def cmd_import_bo2023_buildkit(args: argparse.Namespace) -> int:
    result = import_buildkit_dir(args.db, args.buildkit_dir)
    print(f"Imported Bo2023 buildkit into SQLite: {result['tables_imported']} tables")
    return 0





def cmd_build_signature(args: argparse.Namespace) -> int:
    sigset_id = build_signature_set(
        db_path=args.db,
        atlas_id=args.atlas_id,
        method=args.method,
        topk_per_region=args.topk_per_region,
        remove_housekeeping=not args.keep_housekeeping,
        remove_blood_background=not args.keep_blood_background,
    )
    print(json.dumps({'sigset_id': sigset_id, 'atlas_id': args.atlas_id}, ensure_ascii=False, indent=2))
    return 0

def cmd_import_bo2023_region_matrix(args: argparse.Namespace) -> int:
    result = import_region_matrix(
        db_path=args.db,
        matrix_path=args.matrix,
        annotation_path=args.annotation,
        atlas_name=args.atlas_name,
        build_version=args.build_version,
        gene_id_type=args.gene_id_type,
        normalization=args.normalization,
        notes=args.notes,
    )
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return 0

def cmd_export_run(args: argparse.Namespace) -> int:
    out = export_run_bundle(args.db, args.run_id, args.output_zip)
    print(out)
    return 0


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog='cfrna-tracing', description='cfRNA source tracing CLI')
    sub = p.add_subparsers(dest='command', required=True)

    common = argparse.ArgumentParser(add_help=False)
    common.add_argument('--db', required=True, help='Path to SQLite database')
    common.add_argument('--method', default='ensemble', choices=method_choices())
    common.add_argument('--topk', type=int, default=3)
    common.add_argument('--atlas-id', type=int, default=1)
    common.add_argument('--sigset-id', type=int, default=None)
    common.add_argument('--use-value', default='log1p', choices=['log1p', 'tpm', 'zscore'])
    common.add_argument('--l2', type=float, default=1e-4)
    common.add_argument('--ensemble-alpha', type=float, default=0.5)
    common.add_argument('--bootstrap-n', type=int, default=50)
    common.add_argument('--bootstrap-gene-frac', type=float, default=0.7)
    common.add_argument('--label-key', default='')
    common.add_argument('--limit', type=int, default=None)
    common.add_argument('--prefix', default='benchmark_papergrade')

    p_bench = sub.add_parser('benchmark', parents=[common], help='Run benchmark and export Figure1-Figure6 plus PDF report')
    p_bench.add_argument('--output-dir', required=True)
    p_bench.set_defaults(func=cmd_benchmark)

    p_bundle = sub.add_parser('benchmark-bundle', parents=[common], help='Run benchmark and write a zipped figure/report bundle')
    p_bundle.add_argument('--output-zip', required=True)
    p_bundle.set_defaults(func=cmd_benchmark_bundle)

    p_imp = sub.add_parser('import-bo2023-buildkit', help='Import Bo2023 supplementary buildkit CSVs into SQLite')
    p_imp.add_argument('--db', required=True)
    p_imp.add_argument('--buildkit-dir', required=True)
    p_imp.set_defaults(func=cmd_import_bo2023_buildkit)

    p_mat = sub.add_parser('import-bo2023-region-matrix', help='Import reconstructed Bo2023 gene×region TPM matrix into SQLite as a selectable atlas')
    p_mat.add_argument('--db', required=True)
    p_mat.add_argument('--matrix', required=True, help='Path to bo2023_gene_by_region_mean_tpm.tsv.gz')
    p_mat.add_argument('--annotation', required=True, help='Path to sample_annotation_master_auto_brain_region.tsv')
    p_mat.add_argument('--atlas-name', default='WangLab Bo2023 reconstructed bulk atlas')
    p_mat.add_argument('--build-version', default='reconstructed_from_PRJNA905082')
    p_mat.add_argument('--gene-id-type', default='gene_symbol')
    p_mat.add_argument('--normalization', default='TPM')
    p_mat.add_argument('--notes', default='')
    p_mat.set_defaults(func=cmd_import_bo2023_region_matrix)

    p_sig = sub.add_parser('build-signature-set', help='Build signature genes for a chosen atlas')
    p_sig.add_argument('--db', required=True)
    p_sig.add_argument('--atlas-id', type=int, required=True)
    p_sig.add_argument('--method', default='hybrid_specificity')
    p_sig.add_argument('--topk-per-region', type=int, default=120)
    p_sig.add_argument('--keep-housekeeping', action='store_true')
    p_sig.add_argument('--keep-blood-background', action='store_true')
    p_sig.set_defaults(func=cmd_build_signature)

    p_run = sub.add_parser('export-run', help='Export one analysis run from the database')
    p_run.add_argument('--db', required=True)
    p_run.add_argument('--run-id', required=True)
    p_run.add_argument('--output-zip', default=None)
    p_run.set_defaults(func=cmd_export_run)
    return p


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return int(args.func(args))


if __name__ == '__main__':
    raise SystemExit(main())
