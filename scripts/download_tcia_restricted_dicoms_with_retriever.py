from __future__ import annotations

import argparse
import getpass
import os
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_RETRIEVER = ROOT / "tools" / "tcia_data_retriever" / "extracted" / "TCIA_Data_Retriever.exe"
DEFAULT_INPUT = (
    ROOT
    / "data"
    / "tcia_tcga_glioma_mri"
    / "tcia_data_retriever_manifests"
    / "tcga_73_complete_brats4_retriever_series_uid_only.csv"
)
DEFAULT_OUTPUT = ROOT / "data" / "tcia_tcga_glioma_mri" / "dicom_raw_73_nbia"


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Download restricted TCIA/NBIA DICOM series with TCIA Data Retriever."
    )
    parser.add_argument("--retriever", type=Path, default=DEFAULT_RETRIEVER)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--user", default=os.environ.get("TCIA_USER", ""))
    parser.add_argument("--processes", default="2")
    parser.add_argument("--max-connections", default="4")
    args = parser.parse_args()

    if not args.user:
        raise SystemExit(
            "Set TCIA_USER or pass --user before running. Do not paste credentials into chat."
        )
    password = os.environ.get("TCIA_PASS", "") or getpass.getpass("TCIA/NBIA password: ")
    if not password:
        raise SystemExit("Password is required.")

    cmd = [
        str(args.retriever),
        "--cli",
        "--input",
        str(args.input),
        "--output",
        str(args.output),
        "--image-url",
        "https://services.cancerimagingarchive.net/nbia-api/services/v2/getImage",
        "--meta-url",
        "https://services.cancerimagingarchive.net/nbia-api/services/v2/getSeriesMetaData",
        "--token-url",
        "https://services.cancerimagingarchive.net/nbia-api/oauth/token",
        "--user",
        args.user,
        "--passwd",
        password,
        "--processes",
        str(args.processes),
        "--max-connections",
        str(args.max_connections),
        "--server-friendly",
        "--skip-existing",
        "--verbose",
        "--save-log",
    ]

    args.output.mkdir(parents=True, exist_ok=True)
    proc = subprocess.run(cmd, cwd=ROOT, text=True)
    raise SystemExit(proc.returncode)


if __name__ == "__main__":
    main()
