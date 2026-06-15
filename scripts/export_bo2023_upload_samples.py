#!/usr/bin/env python
from __future__ import annotations

import argparse
import json
import re
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.run_bo2023_loso_validation import read_annotations, read_vsd_matrix  # noqa: E402
from scripts.run_bo2023_v2_loso_validation import (  # noqa: E402
    DEFAULT_GENE_MAP,
    DEFAULT_MATRIX,
    DEFAULT_SAMPLE_INFO,
    map_matrix_to_symbols,
)


DEFAULT_OUTDIR = ROOT / "exports" / "bo2023_vsd_upload_samples_819_20260525"
UPLOAD_MIN_VALUE = 0.100001


def safe_filename_token(value: object) -> str:
    token = str(value).strip()
    token = re.sub(r'[<>:"/\\|?*]+', "-", token)
    token = re.sub(r"\s+", "_", token).strip(" ._-")
    return token or "NA"


def cell_text(row: pd.Series, column: str) -> str:
    value = row.get(column, "")
    return "" if pd.isna(value) else str(value).strip()


def main() -> None:
    parser = argparse.ArgumentParser(description="Export Bo2023 reference samples as directly uploadable CSV files.")
    parser.add_argument("--matrix", type=Path, default=DEFAULT_MATRIX)
    parser.add_argument("--sample-info", type=Path, default=DEFAULT_SAMPLE_INFO)
    parser.add_argument("--gene-map", type=Path, default=DEFAULT_GENE_MAP)
    parser.add_argument("--sample-sheet", default="mfas5_819samples_phenSet4")
    parser.add_argument("--region-col", default="Region")
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    args = parser.parse_args()

    outdir = args.outdir
    samples_dir = outdir / "samples"
    if samples_dir.exists() and any(samples_dir.iterdir()):
        raise FileExistsError(f"Refusing to mix exports into non-empty directory: {samples_dir}")
    samples_dir.mkdir(parents=True, exist_ok=True)

    raw_matrix = read_vsd_matrix(args.matrix)
    matrix = map_matrix_to_symbols(raw_matrix, args.gene_map)
    ann = read_annotations(args.sample_info, args.sample_sheet, args.region_col)
    ann = ann.set_index("sample_id", drop=False)

    matrix_samples = matrix.columns.astype(str).tolist()
    missing_annotations = sorted(set(matrix_samples) - set(ann.index))
    extra_annotations = sorted(set(ann.index) - set(matrix_samples))
    if missing_annotations or extra_annotations or len(matrix_samples) != 819:
        raise ValueError(
            "Sample alignment failed: "
            f"matrix={len(matrix_samples)}, missing_annotations={len(missing_annotations)}, "
            f"extra_annotations={len(extra_annotations)}"
        )

    expression_values = matrix.to_numpy(dtype=float)
    if not np.isfinite(expression_values).all():
        raise ValueError("The upload schema requires finite expression values.")

    manifest_rows: list[dict[str, object]] = []
    gene_symbols = matrix.index.astype(str).to_numpy()
    for sample_id in matrix_samples:
        info = ann.loc[sample_id]
        region_id = cell_text(info, "region_id")
        filename_region = safe_filename_token(region_id)
        filename = f"{safe_filename_token(sample_id)}__{filename_region}.csv"
        raw_values = matrix[sample_id].to_numpy(dtype=float)
        raw_min = float(raw_values.min())
        upload_shift = max(0.0, UPLOAD_MIN_VALUE - raw_min)
        upload_values = raw_values + upload_shift

        upload_table = pd.DataFrame(
            {
                "gene_symbol": gene_symbols,
                "tpm_value": upload_values,
                "sample_id": "",
                "ground_truth_region": "",
            }
        )
        upload_table.loc[0, "sample_id"] = sample_id
        upload_table.loc[0, "ground_truth_region"] = region_id
        upload_table.to_csv(samples_dir / filename, index=False, encoding="utf-8-sig")

        manifest_rows.append(
            {
                "sample_id": sample_id,
                "source_region": region_id,
                "filename_region": filename_region,
                "filename": filename,
                "relative_path": f"samples/{filename}",
                "lobe": cell_text(info, "Lobe"),
                "saleem_network": cell_text(info, "SaleemNetworks"),
                "n_genes": int(matrix.shape[0]),
                "expression_value_column": "tpm_value",
                "expression_scale": "Bo2023 VSD batch-removed with per-sample upload shift",
                "raw_vsd_min": raw_min,
                "upload_shift": upload_shift,
                "uploaded_min_value": float(upload_values.min()),
                "tracing_input_normalization": "vsd",
            }
        )

    manifest = pd.DataFrame(manifest_rows)
    if manifest["filename"].duplicated().any():
        duplicates = manifest.loc[manifest["filename"].duplicated(keep=False), "filename"].tolist()
        raise ValueError(f"Generated duplicate filenames: {duplicates[:5]}")
    manifest.to_csv(outdir / "upload_manifest.csv", index=False, encoding="utf-8-sig")

    summary = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "source_matrix": str(args.matrix),
        "source_annotations": str(args.sample_info),
        "gene_map": str(args.gene_map),
        "n_samples": int(matrix.shape[1]),
        "n_gene_symbols_per_file": int(matrix.shape[0]),
        "raw_min_expression_value": float(expression_values.min()),
        "raw_max_expression_value": float(expression_values.max()),
        "upload_min_value_floor": UPLOAD_MIN_VALUE,
        "n_samples_shifted_for_upload": int((manifest["upload_shift"] > 0).sum()),
        "max_upload_shift": float(manifest["upload_shift"].max()),
        "output_directory": str(outdir),
        "sample_files_directory": str(samples_dir),
        "filename_pattern": "<sample_id>__<filename-safe source_region>.csv",
        "filename_note": "The source Region remains unchanged in ground_truth_region and the manifest.",
        "upload_columns": ["gene_symbol", "tpm_value", "sample_id", "ground_truth_region"],
        "tracing_input_normalization": "vsd",
        "upload_shift_note": "Each sample is shifted by one constant only when needed; Pearson correlation rankings are unchanged.",
    }
    with (outdir / "export_summary.json").open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, ensure_ascii=False, indent=2)
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
