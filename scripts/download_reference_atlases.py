#!/usr/bin/env python
"""Download and organize reference atlases for human brain injury cfRNA tracing.

The script favors processed matrices, aggregate expression tables, and metadata.
FASTQ/BAM/CRAM/SRA files are excluded unless --raw is provided, and very large
single-cell matrices are excluded unless --large is provided.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import logging
import os
import re
import shutil
import subprocess
import sys
import time
import zipfile
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable
from urllib.parse import unquote, urljoin, urlparse

import pandas as pd
import requests
try:
    from tqdm import tqdm
except ModuleNotFoundError:
    class tqdm:  # type: ignore[no-redef]
        def __init__(self, *args, **kwargs):
            self.total = kwargs.get("total")

        def __enter__(self):
            return self

        def __exit__(self, exc_type, exc, tb):
            return False

        def update(self, n):
            return None

try:
    import yaml
except ModuleNotFoundError:  # keep the downloader usable before optional deps are installed
    yaml = None

try:
    from bs4 import BeautifulSoup
except ModuleNotFoundError:
    BeautifulSoup = None


ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "data" / "reference_atlases"
MANIFEST_DIR = BASE / "00_manifest"
LOG_DIR = BASE / "99_logs"

RAW_EXTENSIONS = (".fastq", ".fq", ".bam", ".cram", ".sra", ".bigwig", ".bw")
LARGE_HINTS = (
    "cellxgene",
    "full",
    "raw_counts",
    "gene_tpm.gct.gz",
    "gene_reads.gct.gz",
    "transcript_",
    ".h5ad",
    ".loom",
)

MANIFEST_FIELDS = [
    "atlas_name",
    "publication",
    "doi",
    "source_url",
    "accession",
    "file_name",
    "file_type",
    "expected_format",
    "download_status",
    "md5_or_size_if_available",
    "notes",
    "recommended_use",
]


@dataclass
class AtlasFile:
    atlas_name: str
    publication: str
    doi: str
    source_url: str
    accession: str
    file_name: str
    file_type: str
    expected_format: str
    recommended_use: str
    destination: Path
    notes: str = ""
    requires_large: bool = False
    requires_raw: bool = False
    manual: bool = False
    download_status: str = "planned"
    md5_or_size_if_available: str = ""

    def manifest_row(self) -> dict[str, str]:
        row = asdict(self)
        row.pop("destination", None)
        row.pop("requires_large", None)
        row.pop("requires_raw", None)
        row.pop("manual", None)
        return {k: str(row.get(k, "")) for k in MANIFEST_FIELDS}


def ensure_dirs() -> None:
    for path in [
        MANIFEST_DIR,
        LOG_DIR,
        BASE / "01_normal_human_brain_celltype",
        BASE / "02_normal_human_brain_region",
        BASE / "03_brain_injury_state",
        BASE / "04_peripheral_background",
    ]:
        path.mkdir(parents=True, exist_ok=True)


def setup_logging() -> None:
    ensure_dirs()
    log_path = LOG_DIR / f"download_reference_atlases_{time.strftime('%Y%m%d_%H%M%S')}.log"
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        handlers=[logging.FileHandler(log_path, encoding="utf-8"), logging.StreamHandler(sys.stdout)],
    )
    logging.info("Log file: %s", log_path)


def make_session() -> requests.Session:
    session = requests.Session()
    session.headers.update(
        {
            "User-Agent": "cfRNA-reference-atlas-downloader/1.0 (+https://github.com/; research use)",
            "Accept": "text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8",
        }
    )
    return session


def safe_filename_from_url(url: str, fallback: str = "downloaded_file") -> str:
    parsed = urlparse(url)
    name = unquote(Path(parsed.path).name)
    return name or fallback


def looks_raw(url_or_name: str) -> bool:
    value = url_or_name.lower()
    return any(ext in value for ext in RAW_EXTENSIONS)


def looks_large(url_or_name: str) -> bool:
    value = url_or_name.lower()
    return any(hint in value for hint in LARGE_HINTS)


def infer_file_type(name: str) -> tuple[str, str]:
    lower = name.lower()
    if lower.endswith(".h5ad"):
        return "single_cell_matrix", "h5ad"
    if lower.endswith(".loom"):
        return "single_cell_matrix", "loom"
    if lower.endswith(".mtx") or lower.endswith(".mtx.gz"):
        return "sparse_matrix", "mtx"
    if lower.endswith(".gct") or lower.endswith(".gct.gz"):
        return "expression_matrix", "gct"
    if lower.endswith(".csv") or lower.endswith(".csv.gz"):
        return "table", "csv"
    if lower.endswith(".tsv") or lower.endswith(".tsv.gz") or lower.endswith(".txt") or lower.endswith(".txt.gz"):
        return "table", "tsv"
    if lower.endswith(".zip"):
        return "compressed_archive", "zip"
    if lower.endswith(".gz"):
        return "compressed_file", "gz"
    return "file", Path(lower).suffix.lstrip(".") or "unknown"


def append_failed(row: AtlasFile, reason: str) -> None:
    path = LOG_DIR / "failed_downloads.tsv"
    write_header = not path.exists()
    with path.open("a", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["time", "atlas_name", "source_url", "file_name", "reason"], delimiter="\t")
        if write_header:
            writer.writeheader()
        writer.writerow(
            {
                "time": time.strftime("%Y-%m-%d %H:%M:%S"),
                "atlas_name": row.atlas_name,
                "source_url": row.source_url,
                "file_name": row.file_name,
                "reason": reason,
            }
        )


def append_failed_url(atlas_name: str, url: str, file_name: str, reason: str) -> None:
    path = LOG_DIR / "failed_downloads.tsv"
    write_header = not path.exists()
    with path.open("a", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["time", "atlas_name", "source_url", "file_name", "reason"], delimiter="\t")
        if write_header:
            writer.writeheader()
        writer.writerow(
            {
                "time": time.strftime("%Y-%m-%d %H:%M:%S"),
                "atlas_name": atlas_name,
                "source_url": url,
                "file_name": file_name,
                "reason": reason,
            }
        )


def fetch_html(session: requests.Session, url: str, timeout: int = 45) -> str | None:
    try:
        response = session.get(url, timeout=timeout)
        response.raise_for_status()
    except requests.RequestException as exc:
        logging.warning("Could not fetch %s: %s", url, exc)
        append_failed_url("page_parse", url, "__page_parse__", str(exc))
        return None
    return response.text


def link_text(anchor) -> str:
    return " ".join(anchor.get_text(" ", strip=True).split())


def parse_links(
    session: requests.Session,
    page_url: str,
    keywords: Iterable[str],
    extensions: Iterable[str] = (".zip", ".gz", ".tsv", ".csv", ".txt", ".h5ad", ".loom", ".mtx"),
) -> list[tuple[str, str, str]]:
    html = fetch_html(session, page_url)
    if html is None:
        return []
    keyword_re = re.compile("|".join(re.escape(k) for k in keywords), re.I) if keywords else None
    found: list[tuple[str, str, str]] = []
    if BeautifulSoup is not None:
        soup = BeautifulSoup(html, "html.parser")
        anchors = [(anchor["href"].strip(), link_text(anchor)) for anchor in soup.find_all("a", href=True)]
    else:
        anchors = []
        for match in re.finditer(r"<a\b[^>]*href=[\"']([^\"']+)[\"'][^>]*>(.*?)</a>", html, flags=re.I | re.S):
            href = match.group(1).strip()
            text = re.sub(r"<[^>]+>", " ", match.group(2))
            text = " ".join(text.split())
            anchors.append((href, text))
    for href, text in anchors:
        context = f"{text} {href}"
        full_url = urljoin(page_url, href)
        lower_url = full_url.lower()
        has_ext = any(ext in lower_url for ext in extensions)
        has_keyword = bool(keyword_re.search(context)) if keyword_re else True
        if has_ext and has_keyword:
            found.append((full_url, text or safe_filename_from_url(full_url), context[:500]))
    dedup: dict[str, tuple[str, str, str]] = {}
    for item in found:
        dedup[item[0]] = item
    return list(dedup.values())


def github_raw_url(url: str) -> str:
    return url.replace("https://github.com/", "https://raw.githubusercontent.com/").replace("/blob/", "/")


def create_manual_note(path: Path, title: str, lines: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    body = "\n".join([f"# {title}", "", *lines, ""])
    if path.exists() and path.read_text(encoding="utf-8") == body:
        return
    path.write_text(body, encoding="utf-8")
    logging.info("Wrote manual download note: %s", path)


def aws_cli_available() -> bool:
    return shutil.which("aws") is not None


def list_public_s3(prefix: str) -> list[str]:
    if not aws_cli_available():
        return []
    try:
        result = subprocess.run(
            ["aws", "s3", "ls", "--no-sign-request", prefix],
            capture_output=True,
            text=True,
            timeout=120,
            check=False,
        )
    except Exception as exc:  # noqa: BLE001
        logging.warning("Could not run aws s3 ls for %s: %s", prefix, exc)
        return []
    if result.returncode != 0:
        logging.warning("aws s3 ls failed for %s: %s", prefix, result.stderr.strip())
        return []
    urls = []
    for line in result.stdout.splitlines():
        parts = line.split()
        if not parts:
            continue
        name = parts[-1]
        if name.endswith("/"):
            continue
        urls.append(prefix.rstrip("/") + "/" + name)
    return urls


def write_outputs(rows: list[AtlasFile], source_name: str | None = None) -> None:
    ensure_dirs()
    table = pd.DataFrame([row.manifest_row() for row in rows], columns=MANIFEST_FIELDS)
    if source_name:
        table.to_csv(MANIFEST_DIR / f"{source_name}_manifest.tsv", sep="\t", index=False)
    table.to_csv(MANIFEST_DIR / "reference_atlas_manifest.tsv", sep="\t", index=False)
    table.to_json(MANIFEST_DIR / "reference_atlas_manifest.json", orient="records", indent=2, force_ascii=False)
    with (MANIFEST_DIR / "reference_atlas_manifest.yaml").open("w", encoding="utf-8") as handle:
        records = table.to_dict(orient="records")
        if yaml is not None:
            yaml.safe_dump(records, handle, allow_unicode=True, sort_keys=False)
        else:
            handle.write(json.dumps(records, indent=2, ensure_ascii=False))


def md5sum(path: Path, chunk_size: int = 1024 * 1024) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def download_file(session: requests.Session, row: AtlasFile, dry_run: bool = False) -> AtlasFile:
    row.destination.mkdir(parents=True, exist_ok=True)
    target = row.destination / row.file_name
    if row.manual:
        row.download_status = "manual_required"
        return row
    if dry_run:
        row.download_status = "dry_run_planned"
        return row
    if target.exists() and target.stat().st_size > 0:
        row.download_status = "skipped_exists"
        row.md5_or_size_if_available = f"size={target.stat().st_size}"
        logging.info("Skip existing file: %s", target)
        return row

    temp = target.with_suffix(target.suffix + ".part")
    headers = {}
    existing = temp.stat().st_size if temp.exists() else 0
    if existing:
        headers["Range"] = f"bytes={existing}-"

    logging.info("Downloading %s -> %s", row.source_url, target)
    try:
        with session.get(row.source_url, stream=True, timeout=120, headers=headers) as response:
            response.raise_for_status()
            if existing and response.status_code != 206:
                existing = 0
                temp.unlink(missing_ok=True)
            total = response.headers.get("Content-Length")
            total_size = int(total) + existing if total and total.isdigit() else None
            mode = "ab" if existing else "wb"
            with temp.open(mode) as handle, tqdm(
                total=total_size,
                initial=existing,
                unit="B",
                unit_scale=True,
                desc=row.file_name[:40],
            ) as bar:
                for chunk in response.iter_content(chunk_size=1024 * 1024):
                    if chunk:
                        handle.write(chunk)
                        bar.update(len(chunk))
        temp.replace(target)
        row.download_status = "downloaded"
        row.md5_or_size_if_available = f"size={target.stat().st_size};md5={md5sum(target)}"
        extract_if_needed(target)
        quick_check(target)
    except Exception as exc:  # noqa: BLE001 - continue all remaining downloads
        row.download_status = "failed"
        row.notes = f"{row.notes}; download failed: {exc}".strip("; ")
        append_failed(row, str(exc))
        logging.exception("Failed download: %s", row.source_url)
    return row


def extract_if_needed(path: Path) -> None:
    lower = path.name.lower()
    try:
        if lower.endswith(".zip"):
            out_dir = path.with_suffix("")
            out_dir.mkdir(exist_ok=True)
            with zipfile.ZipFile(path) as archive:
                archive.extractall(out_dir)
            logging.info("Extracted zip to %s", out_dir)
        elif lower.endswith(".gz") and not lower.endswith((".h5ad.gz", ".mtx.gz")):
            out_path = path.with_suffix("")
            if not out_path.exists():
                with gzip.open(path, "rb") as src, out_path.open("wb") as dst:
                    shutil.copyfileobj(src, dst)
                logging.info("Decompressed gz to %s", out_path)
    except Exception as exc:  # noqa: BLE001
        logging.warning("Could not extract %s: %s", path, exc)


def quick_check(path: Path) -> None:
    lower = path.name.lower()
    candidate = path
    if lower.endswith(".gz") and not lower.endswith((".h5ad.gz", ".mtx.gz")):
        candidate = path.with_suffix("")
        lower = candidate.name.lower()
    try:
        if lower.endswith(".gct"):
            with candidate.open("r", encoding="utf-8", errors="replace") as handle:
                _ = handle.readline()
                _ = handle.readline()
                header = handle.readline().rstrip("\n").split("\t")
            logging.info("GCT header columns=%d first=%s", len(header), header[:8])
        elif lower.endswith((".tsv", ".csv", ".txt")):
            sep = "," if lower.endswith(".csv") else "\t"
            preview = pd.read_csv(candidate, sep=sep, nrows=5)
            columns = [str(c).lower() for c in preview.columns]
            kind = "expression_matrix" if len(columns) > 20 else "metadata_or_annotation"
            logging.info("Preview %s: shape=%s inferred=%s", candidate.name, preview.shape, kind)
        elif lower.endswith(".h5ad"):
            logging.info("H5AD exists, size=%d; not loading during download", candidate.stat().st_size)
    except Exception as exc:  # noqa: BLE001
        logging.warning("Preview check failed for %s: %s", candidate, exc)


def atlas_file_from_url(
    url: str,
    destination: Path,
    atlas_name: str,
    publication: str,
    doi: str,
    accession: str,
    recommended_use: str,
    notes: str = "",
    file_name: str | None = None,
) -> AtlasFile:
    name = file_name or safe_filename_from_url(url)
    file_type, expected_format = infer_file_type(name)
    return AtlasFile(
        atlas_name=atlas_name,
        publication=publication,
        doi=doi,
        source_url=url,
        accession=accession,
        file_name=name,
        file_type=file_type,
        expected_format=expected_format,
        recommended_use=recommended_use,
        destination=destination,
        notes=notes,
        requires_large=looks_large(url) or looks_large(name),
        requires_raw=looks_raw(url) or looks_raw(name),
    )


def source_siletti(session: requests.Session) -> list[AtlasFile]:
    dest = BASE / "01_normal_human_brain_celltype" / "siletti2023_adult_human_brain"
    rows: list[AtlasFile] = []
    create_manual_note(
        dest / "manual_download_siletti_abc_s3.md",
        "Siletti 2023 / Allen ABC Atlas S3 Manual Download",
        [
            "- GitHub: https://github.com/linnarsson-lab/adult-human-brain",
            "- HCA atlas: https://data.humancellatlas.org/hca-bio-networks/nervous-system/atlases/brain-v1-0",
            "- CELLxGENE collection: https://cellxgene.cziscience.com/collections/283d65eb-dd53-496d-adb7-7570c7caa443",
            "- Allen WHB dataset: https://alleninstitute.github.io/abc_atlas_access/descriptions/WHB_dataset.html",
            "- Allen WHB 10Xv3: https://alleninstitute.github.io/abc_atlas_access/descriptions/WHB-10Xv3.html",
            "- Default first-version downloads should target metadata, taxonomy, cluster annotation, and gene metadata.",
            "- List metadata with: `aws s3 ls --no-sign-request s3://allen-brain-cell-atlas/metadata/WHB-10Xv3/20241115/`",
            "- Large expression matrices are about 70GB and should only be copied with `--large`: `aws s3 ls --no-sign-request s3://allen-brain-cell-atlas/expression_matrices/WHB-10Xv3/20240330/`",
            "- MapMyCells output is about 8GB and should only be copied with `--large`: `aws s3 ls --no-sign-request s3://allen-brain-cell-atlas/mapmycells/WHB-10Xv3/20240831/`",
            "- Do not copy raw FASTQ/BAM by default.",
        ],
    )
    s3_metadata = "s3://allen-brain-cell-atlas/metadata/WHB-10Xv3/20241115/"
    s3_expression = "s3://allen-brain-cell-atlas/expression_matrices/WHB-10Xv3/20240330/"
    s3_mapmycells = "s3://allen-brain-cell-atlas/mapmycells/WHB-10Xv3/20240831/"
    s3_urls = list_public_s3(s3_metadata)
    if not s3_urls:
        rows.append(
            AtlasFile(
                "Siletti2023_adult_human_whole_brain",
                "Siletti et al., Science 2023",
                "10.1126/science.add7046",
                s3_metadata,
                "s3://allen-brain-cell-atlas/metadata/WHB-10Xv3/20241115/",
                "manual_download_siletti_abc_s3.md",
                "manual_download_instruction",
                "markdown",
                "celltype_reference",
                dest,
                notes="Install/configure AWS CLI or use manual S3 browser. Metadata/taxonomy/annotations should be downloaded before full matrices.",
                manual=True,
            )
        )
    for s3_url in s3_urls:
        name = safe_filename_from_url(s3_url)
        file_type, expected_format = infer_file_type(name)
        rows.append(
            AtlasFile(
                "Siletti2023_adult_human_whole_brain",
                "Siletti et al., Science 2023",
                "10.1126/science.add7046",
                s3_url,
                "s3://allen-brain-cell-atlas/metadata/WHB-10Xv3/20241115/",
                name,
                file_type,
                expected_format,
                "celltype_reference",
                dest,
                notes="Public Allen ABC Atlas WHB-10Xv3 metadata; download with `aws s3 cp --no-sign-request` or manual S3 access.",
                manual=True,
            )
        )
    for large_prefix, label in [(s3_expression, "70GB expression matrices"), (s3_mapmycells, "8GB MapMyCells output")]:
        rows.append(
            AtlasFile(
                "Siletti2023_adult_human_whole_brain",
                "Siletti et al., Science 2023",
                "10.1126/science.add7046",
                large_prefix,
                large_prefix,
                f"siletti_{label.replace(' ', '_')}_s3_prefix",
                "s3_prefix",
                "s3",
                "celltype_reference",
                dest,
                notes=f"{label}; list with `aws s3 ls --no-sign-request {large_prefix}` and download only with --large.",
                requires_large=True,
                manual=True,
            )
        )
    for url, text, context in parse_links(
        session,
        "https://github.com/linnarsson-lab/adult-human-brain",
        ["h5ad", "loom", "cluster", "taxonomy", "metadata", "annotation", "expression"],
    ):
        if "github.com" in url and "/blob/" in url:
            url = github_raw_url(url)
        rows.append(
            atlas_file_from_url(
                url,
                dest,
                "Siletti2023_adult_human_whole_brain",
                "Siletti et al., Science 2023",
                "10.1126/science.add7046",
                "CELLxGENE:283d65eb-dd53-496d-adb7-7570c7caa443; s3://allen-brain-cell-atlas",
                "celltype_reference",
                notes=f"Parsed from GitHub: {text}; {context}",
            )
        )
    for page in [
        "https://data.humancellatlas.org/hca-bio-networks/nervous-system/atlases/brain-v1-0",
        "https://alleninstitute.github.io/abc_atlas_access/descriptions/WHB_dataset.html",
        "https://alleninstitute.github.io/abc_atlas_access/descriptions/WHB-10Xv3.html",
    ]:
        for url, text, context in parse_links(
            session,
            page,
            ["metadata", "taxonomy", "cluster", "annotation", "gene", "WHB", "10Xv3"],
        ):
            row = atlas_file_from_url(
                url,
                dest,
                "Siletti2023_adult_human_whole_brain",
                "Siletti et al., Science 2023",
                "10.1126/science.add7046",
                "HCA/CELLxGENE/Allen ABC Atlas WHB-10Xv3",
                "celltype_reference",
                notes=f"Parsed from {page}: {text}; {context}",
            )
            if "expression_matrices" in url or "mapmycells" in url:
                row.requires_large = True
            rows.append(row)
    rows.append(
        AtlasFile(
            "Siletti2023_adult_human_whole_brain",
            "Siletti et al., Science 2023",
            "10.1126/science.add7046",
            "https://cellxgene.cziscience.com/collections/283d65eb-dd53-496d-adb7-7570c7caa443",
            "CELLxGENE:283d65eb-dd53-496d-adb7-7570c7caa443",
            "cellxgene_collection_manual_or_large_h5ad",
            "single_cell_matrix",
            "h5ad",
            "celltype_reference",
            dest,
            notes="Full CELLxGENE h5ad is large; use --large only if a direct asset URL is discovered manually.",
            requires_large=True,
            manual=True,
        )
    )
    return rows


def source_brain_cell_atlas(session: requests.Session) -> list[AtlasFile]:
    dest = BASE / "01_normal_human_brain_celltype" / "brain_cell_atlas_chen2024"
    pages = ["https://www.braincellatlas.org/dataSet", "https://www.braincellatlas.org", "https://www.nature.com/articles/s41591-024-03150-z"]
    rows: list[AtlasFile] = []
    for page in pages:
        for url, text, context in parse_links(session, page, ["human", "brain", "adult", "healthy", "h5ad", "metadata", "count", "matrix", "sparse"]):
            rows.append(
                atlas_file_from_url(
                    url,
                    dest,
                    "Chen2024_Brain_Cell_Atlas",
                    "Chen et al., Nature Medicine 2024",
                    "10.1038/s41591-024-03150-z",
                    "Brain Cell Atlas portal",
                    "celltype_reference",
                    notes=f"Parsed from {page}: {text}; {context}",
                )
            )
    create_manual_note(
        dest.parent / "brain_cell_atlas_manual_download.md",
        "Brain Cell Atlas Manual Download",
        [
            "- Portal: https://www.braincellatlas.org",
            "- Data page: https://www.braincellatlas.org/dataSet",
            "- Article: https://www.nature.com/articles/s41591-024-03150-z",
            "- Priority files: human adult healthy h5ad, raw counts or normalized counts, cell metadata, cell type annotation, and brain region annotation.",
            "- First-version recommendation: find aggregated or cell-type average matrices before downloading the full cell-level matrix.",
            "- Do not download the full 11M-level cell atlas by default.",
            "- If only h5ad is available, later processing should use `anndata.read_h5ad(..., backed=\"r\")` and aggregate by cell type and brain region.",
            "- Save files under `data/reference_atlases/01_normal_human_brain_celltype/brain_cell_atlas_chen2024/`.",
        ],
    )
    if not rows:
        rows.append(
            AtlasFile(
                "Chen2024_Brain_Cell_Atlas",
                "Chen et al., Nature Medicine 2024",
                "10.1038/s41591-024-03150-z",
                "https://www.braincellatlas.org/dataSet",
                "Brain Cell Atlas portal",
                "brain_cell_atlas_manual_download.md",
                "manual_download_instruction",
                "markdown",
                "celltype_reference",
                dest.parent,
                notes="Dataset page may require JavaScript or portal interaction.",
                manual=True,
            )
        )
    return rows


def source_allen_human_brain_atlas(session: requests.Session) -> list[AtlasFile]:
    dest = BASE / "02_normal_human_brain_region" / "allen_human_brain_atlas"
    page = "https://human.brain-map.org/static/download"
    rows: list[AtlasFile] = []
    donor_ids = ["H0351.2001", "H0351.2002", "H0351.1009", "H0351.1012", "H0351.1015", "H0351.1016"]
    seed_downloads = [
        ("H0351.2001_normalized_microarray.zip", "https://human.brain-map.org/api/v2/well_known_file_download/178238387"),
        ("H0351.2002_normalized_microarray.zip", "https://human.brain-map.org/api/v2/well_known_file_download/178238373"),
        ("H0351.1009_normalized_microarray.zip", "https://human.brain-map.org/api/v2/well_known_file_download/178238359"),
        ("H0351.1012_normalized_microarray.zip", "https://human.brain-map.org/api/v2/well_known_file_download/178238316"),
        ("H0351.1015_normalized_microarray.zip", "https://human.brain-map.org/api/v2/well_known_file_download/178238266"),
        ("H0351.1016_normalized_microarray.zip", "https://human.brain-map.org/api/v2/well_known_file_download/178236545"),
        ("H0351.2001_rnaseq.zip", "https://human.brain-map.org/api/v2/well_known_file_download/278447594"),
        ("H0351.2002_rnaseq.zip", "https://human.brain-map.org/api/v2/well_known_file_download/278448166"),
    ]
    for file_name, url in seed_downloads:
        rows.append(
            atlas_file_from_url(
                url,
                dest,
                "Allen_Human_Brain_Atlas_Hawrylycz2012",
                "Hawrylycz et al., Nature 2012",
                "10.1038/nature11405",
                ",".join(donor_ids),
                "brain_region_reference",
                notes="Seeded from Allen AHBA static download page; archives include expression matrix plus probe/sample/anatomical metadata or RNA-seq counts/TPM.",
                file_name=file_name,
            )
        )
    for url, text, context in parse_links(session, page, donor_ids + ["RNA-Sequencing", "microarray", "normalized", "TPM", "counts", "metadata"]):
        if any(donor in f"{url} {text} {context}" for donor in donor_ids):
            rows.append(
                atlas_file_from_url(
                    url,
                    dest,
                    "Allen_Human_Brain_Atlas_Hawrylycz2012",
                    "Hawrylycz et al., Nature 2012",
                    "10.1038/nature11405",
                    ",".join(donor_ids),
                    "brain_region_reference",
                    notes=f"Parsed from Allen AHBA download page: {text}; {context}",
                )
            )
    if not rows:
        rows.append(
            AtlasFile(
                "Allen_Human_Brain_Atlas_Hawrylycz2012",
                "Hawrylycz et al., Nature 2012",
                "10.1038/nature11405",
                page,
                ",".join(donor_ids),
                "allen_human_brain_atlas_download_page_manual_check",
                "download_page",
                "html",
                "brain_region_reference",
                dest,
                notes="No direct links parsed; use Allen download page for six microarray donors and two RNA-seq donors.",
                manual=True,
            )
        )
    return rows


def source_hpa_brain(session: requests.Session) -> list[AtlasFile]:
    dest = BASE / "02_normal_human_brain_region" / "human_protein_atlas_brain"
    rows: list[AtlasFile] = []
    seed_files = [
        "rna_single_nuclei_brain_datasets.tsv.zip",
        "rna_single_nuclei_cluster_type.tsv.zip",
        "rna_single_nuclei_cluster_type_cluster_types.tsv.zip",
        "rna_single_nuclei_brain_cluster.tsv.zip",
        "rna_single_nuclei_brain_clusters.tsv.zip",
        "rna_brain_region_hpa.tsv.zip",
        "rna_brain_region_hpa_brain_regions.tsv.zip",
        "rna_brain_hpa.tsv.zip",
        "rna_pfc_brain_hpa.tsv.zip",
        "rna_pfc_brain_hpa_subregions.tsv.zip",
        "transcript_rna_brain.tsv.zip",
        "transcript_rna_pfcbrain.tsv.zip",
    ]
    for name in seed_files:
        rows.append(
            atlas_file_from_url(
                f"https://www.proteinatlas.org/download/{name}",
                dest,
                "Human_Protein_Atlas_Brain_Sjostedt2020",
                "Sjöstedt et al., Science 2020",
                "10.1126/science.aay5947",
                "Human Protein Atlas v25 brain downloads",
                "brain_region_reference;celltype_reference",
                notes="Seeded HPA download filename; script also parses HPA pages for related files.",
            )
        )
    for page in ["https://www.proteinatlas.org/about/download", "https://www.proteinatlas.org/humanproteome/brain/data", "https://www.proteinatlas.org/humanproteome/single%2Bcell/single%2Bnuclei%2Bbrain/data"]:
        for url, text, context in parse_links(session, page, ["brain", "single nuclei brain", "rna expression", "cluster type", "cluster"]):
            rows.append(
                atlas_file_from_url(
                    url,
                    dest,
                    "Human_Protein_Atlas_Brain_Sjostedt2020",
                    "Sjöstedt et al., Science 2020",
                    "10.1126/science.aay5947",
                    "Human Protein Atlas brain downloads",
                    "brain_region_reference;celltype_reference",
                    notes=f"Parsed from {page}: {text}; {context}",
                )
            )
    return rows


def source_hodge2019(session: requests.Session) -> list[AtlasFile]:
    dest = BASE / "01_normal_human_brain_celltype" / "hodge2019_human_cortex"
    rows: list[AtlasFile] = []
    pages = [
        "https://celltypes.brain-map.org/rnaseq",
        "https://celltypes.brain-map.org/rnaseq/human",
        "https://brain-map.org/our-research/cell-types-taxonomies/cell-types-database-rna-seq-data",
        "https://celltypes.brain-map.org/",
        "https://portal.brain-map.org/atlases-and-data/rnaseq#Human_Cortex",
        "https://biccn.org",
    ]
    for page in pages:
        for url, text, context in parse_links(
            session,
            page,
            [
                "human",
                "MTG",
                "middle temporal gyrus",
                "SMART-seq",
                "single nucleus",
                "RNA-seq",
                "exonic",
                "intronic",
                "raw counts",
                "metadata",
                "cluster",
                "taxonomy",
                "matrix",
            ],
        ):
            rows.append(
                atlas_file_from_url(
                    url,
                    dest,
                    "Hodge2019_human_cortex_MTG",
                    "Hodge et al., Nature 2019",
                    "10.1038/s41586-019-1506-7",
                    "Allen Cell Types / BICCN",
                    "celltype_reference",
                    notes=f"Parsed from {page}: {text}; {context}",
                )
            )
    create_manual_note(
        dest / "hodge2019_manual_download.md",
        "Hodge 2019 Human Cortex MTG Manual Download",
        [
            "- Primary Allen Cell Types RNA-seq entry: https://celltypes.brain-map.org/rnaseq",
            "- Human RNA-seq entry: https://celltypes.brain-map.org/rnaseq/human",
            "- Allen RNA-seq portal: https://portal.brain-map.org/atlases-and-data/rnaseq#Human_Cortex",
            "- Recommended download: gene expression matrix for human MTG / middle temporal gyrus.",
            "- Recommended download: sample metadata.",
            "- Recommended download: gene metadata.",
            "- Recommended download: cluster annotation / taxonomy.",
            "- Search terms: human, MTG, middle temporal gyrus, SMART-seq, single nucleus, RNA-seq, exonic, intronic, raw counts, metadata, cluster, taxonomy.",
            "- This data is used as a human cortex cell type reference.",
            "- Save processed matrices and annotations under `data/reference_atlases/01_normal_human_brain_celltype/hodge2019_human_cortex/`.",
        ],
    )
    if not rows:
        rows.append(
            AtlasFile(
                "Hodge2019_human_cortex_MTG",
                "Hodge et al., Nature 2019",
                "10.1038/s41586-019-1506-7",
                pages[0],
                "Allen Cell Types / BICCN",
                "hodge2019_manual_download.md",
                "manual_download_instruction",
                "markdown",
                "celltype_reference",
                dest,
                notes="No direct downloadable files parsed from public pages.",
                manual=True,
            )
        )
    return rows


def source_garza2023_tbi(session: requests.Session) -> list[AtlasFile]:
    dest = BASE / "03_brain_injury_state" / "garza2023_tbi_gse209552"
    page = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE209552"
    rows: list[AtlasFile] = []
    supp_base = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE209nnn/GSE209552/suppl/"
    seed_files = [
        "GSE209552_TBI_gene_count_matrix_2.csv.gz",
        "GSE209552_hGPC_gene_count_matrix_2.csv.gz",
    ]
    for name in seed_files:
        rows.append(
            atlas_file_from_url(
                supp_base + name,
                dest,
                "Garza2023_human_TBI_snRNAseq",
                "Garza et al., Cell Reports 2023",
                "10.1016/j.celrep.2023.113395",
                "GSE209552",
                "injury_state_reference",
                notes="Seeded GEO supplementary processed count matrix; no SRA raw reads by default.",
            )
        )
    for url, text, context in parse_links(session, page, ["supplementary", "TBI", "hGPC", "metadata", "matrix", "count", "cluster", "csv", "txt"]):
        if "ftp.ncbi.nlm.nih.gov" in url or "download" in url.lower():
            rows.append(
                atlas_file_from_url(
                    url,
                    dest,
                    "Garza2023_human_TBI_snRNAseq",
                    "Garza et al., Cell Reports 2023",
                    "10.1016/j.celrep.2023.113395",
                    "GSE209552",
                    "injury_state_reference",
                    notes=f"Parsed from GEO: {text}; {context}",
                )
            )
    return rows


def source_allen_aging_tbi(session: requests.Session) -> list[AtlasFile]:
    dest = BASE / "03_brain_injury_state" / "allen_aging_dementia_tbi_gse104687"
    rows: list[AtlasFile] = []
    pages = ["https://aging.brain-map.org/download/index", "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE104687"]
    for page in pages:
        for url, text, context in parse_links(session, page, ["FPKM", "TPM", "normalized", "unnormalized", "expression", "metadata", "pathology", "RNA-seq", "csv", "zip", "txt"]):
            rows.append(
                atlas_file_from_url(
                    url,
                    dest,
                    "Allen_Aging_Dementia_TBI_GSE104687",
                    "Allen Aging, Dementia and TBI Study",
                    "GSE104687",
                    "GSE104687",
                    "injury_state_reference",
                    notes=f"Parsed from {page}: {text}; {context}",
                )
            )
    if not rows:
        rows.append(
            AtlasFile(
                "Allen_Aging_Dementia_TBI_GSE104687",
                "Allen Aging, Dementia and TBI Study",
                "GSE104687",
                pages[0],
                "GSE104687",
                "allen_aging_tbi_manual_download_page",
                "download_page",
                "html",
                "injury_state_reference",
                dest,
                notes="No direct links parsed; download FPKM/TPM/metadata from Allen Aging Brain or GEO manually.",
                manual=True,
            )
        )
    return rows


def source_gtex(session: requests.Session) -> list[AtlasFile]:
    del session
    dest = BASE / "04_peripheral_background" / "gtex_v8"
    base_url = "https://storage.googleapis.com/adult-gtex/bulk-gex/v8/rna-seq/"
    files = [
        ("GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz", False),
        ("GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz", True),
        ("GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz", True),
    ]
    rows = []
    for name, large in files:
        row = atlas_file_from_url(
            base_url + name,
            dest,
            "GTEx_v8_bulk_RNAseq",
            "GTEx Consortium, Science 2020",
            "10.1126/science.aaz1776",
            "GTEx v8",
            "peripheral_background",
            notes="Use median TPM by tissue first for peripheral/background filtering; full TPM/read count require --large.",
        )
        row.requires_large = large
        rows.append(row)
    return rows


SOURCE_BUILDERS = {
    "siletti": source_siletti,
    "brain_cell_atlas": source_brain_cell_atlas,
    "allen_human_brain_atlas": source_allen_human_brain_atlas,
    "hpa_brain": source_hpa_brain,
    "hodge2019": source_hodge2019,
    "garza2023_tbi": source_garza2023_tbi,
    "allen_aging_tbi": source_allen_aging_tbi,
    "gtex": source_gtex,
}


def select_sources(args: argparse.Namespace) -> list[str]:
    if args.all:
        return list(SOURCE_BUILDERS)
    if args.source:
        return [args.source]
    raise SystemExit("Use --all or --source SOURCE. See --help.")


def deduplicate(rows: list[AtlasFile]) -> list[AtlasFile]:
    seen: set[tuple[str, str]] = set()
    out = []
    for row in rows:
        key = (row.source_url, row.file_name)
        if key not in seen:
            seen.add(key)
            out.append(row)
    return out


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--all", action="store_true", help="Process all configured atlas sources.")
    parser.add_argument("--source", choices=sorted(SOURCE_BUILDERS), help="Process one source.")
    parser.add_argument("--large", action="store_true", help="Allow large full matrices such as h5ad/loom/full GTEx TPM.")
    parser.add_argument("--raw", action="store_true", help="Allow raw sequencing files. Off by default.")
    parser.add_argument("--dry-run", action="store_true", help="Build manifests and manual notes without downloading.")
    args = parser.parse_args(argv)

    setup_logging()
    session = make_session()
    selected = select_sources(args)
    all_rows: list[AtlasFile] = []
    for name in selected:
        logging.info("Collecting source candidates: %s", name)
        rows = deduplicate(SOURCE_BUILDERS[name](session))
        processed: list[AtlasFile] = []
        for row in rows:
            if row.requires_raw and not args.raw:
                row.download_status = "skipped_raw_requires_--raw"
                processed.append(row)
                continue
            if row.requires_large and not args.large:
                row.download_status = "skipped_large_requires_--large"
                processed.append(row)
                continue
            processed.append(download_file(session, row, dry_run=args.dry_run))
        write_outputs(processed, name)
        all_rows.extend(processed)
    write_outputs(all_rows)
    logging.info("Finished. Manifest: %s", MANIFEST_DIR / "reference_atlas_manifest.tsv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
