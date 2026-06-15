# Human Brain Injury cfRNA Reference Atlases

This project area builds a reproducible local cache under `data/reference_atlases/` for normal human brain cell-type references, brain-region references, injury-state references, and peripheral background correction.

## Directory Layout

```text
data/reference_atlases/
  00_manifest/
  01_normal_human_brain_celltype/
  02_normal_human_brain_region/
  03_brain_injury_state/
  04_peripheral_background/
  99_logs/
```

## Install

```bash
pip install -r requirements_reference_atlases.txt
```

## Download Commands

```bash
python scripts/download_reference_atlases.py --all
python scripts/download_reference_atlases.py --source siletti
python scripts/download_reference_atlases.py --source brain_cell_atlas
python scripts/download_reference_atlases.py --source allen_human_brain_atlas
python scripts/download_reference_atlases.py --source hpa_brain
python scripts/download_reference_atlases.py --source hodge2019
python scripts/download_reference_atlases.py --source garza2023_tbi
python scripts/download_reference_atlases.py --source allen_aging_tbi
python scripts/download_reference_atlases.py --source gtex
python scripts/download_reference_atlases.py --all --large
python scripts/download_reference_atlases.py --all --raw
```

Use `--dry-run` to build manifests and manual download notes without downloading files.

## Inspect and Rebuild Manifest

```bash
python scripts/inspect_reference_atlases.py
python scripts/build_reference_manifest.py
```

`inspect_reference_atlases.py` writes `data/reference_atlases/00_manifest/reference_file_inspection.tsv`. `build_reference_manifest.py` merges per-source manifests and local files into:

- `data/reference_atlases/00_manifest/reference_atlas_manifest.tsv`
- `data/reference_atlases/00_manifest/reference_atlas_manifest.json`
- `data/reference_atlases/00_manifest/reference_atlas_manifest.yaml`

## Data Sources and Intended Use

| Source | DOI / accession | Use | Default behavior | Large/raw behavior |
|---|---|---|---|---|
| Siletti et al. adult human whole brain snRNA-seq | DOI: `10.1126/science.add7046`; CELLxGENE collection `283d65eb-dd53-496d-adb7-7570c7caa443`; Allen ABC Atlas S3 `s3://allen-brain-cell-atlas` | Normal human brain cell-type reference | Downloads Siletti-derived HPA single-nuclei brain aggregate tables when available; parses GitHub for processed annotations/matrices | Full h5ad/loom matrices are skipped unless `--large`; raw NeMO reads are not downloaded |
| Chen et al. Brain Cell Atlas | DOI: `10.1038/s41591-024-03150-z` | Integrated human brain cell-type reference | Parses Brain Cell Atlas pages for public h5ad/metadata/annotation links and writes manual instructions if the portal requires JavaScript | Full 11.3M-cell data is skipped unless `--large` |
| Hodge et al. human cortex / MTG atlas | DOI: `10.1038/s41586-019-1506-7` | Human cortex cell-type reference | Parses Allen Cell Types, legacy Cell Types, BICCN, and Allen RNA-seq pages; writes manual instructions when direct files are not exposed | Large full matrices require `--large` |
| Allen Human Brain Atlas / Hawrylycz et al. | DOI: `10.1038/nature11405` | Normal brain-region reference | Parses `https://human.brain-map.org/static/download` for six normalized adult microarray donors and two adult RNA-seq donors | Raw sequencing files are not used |
| Human Protein Atlas Brain Atlas / Sjöstedt et al. | DOI: `10.1126/science.aay5947` | Brain-region reference and cell-type aggregate reference | Downloads/parses brain region RNA, prefrontal cortex RNA, and single-nuclei brain cluster/type ZIP/TSV files; ZIPs are extracted while originals are retained | Very large archives are controlled by `--large` if detected |
| Garza et al. human TBI snRNA-seq | DOI: `10.1016/j.celrep.2023.113395`; GEO `GSE209552` | Acute TBI injury-state reference | Uses GEO/NCBI supplementary files and seed processed count matrix names; skips SRA raw reads | SRA/raw files require `--raw` |
| Allen Aging, Dementia and TBI Study | GEO `GSE104687` | Chronic TBI / aging brain background reference | Parses Allen Aging Brain and GEO for RNA-seq expression, sample metadata, and pathology metadata | BAM/bigWig/raw files require `--raw` |
| GTEx v8 | DOI: `10.1126/science.aaz1776` | Peripheral tissue background reference | Downloads median TPM by tissue from GTEx v8 Google Storage | Full gene TPM and raw reads GCT files require `--large` |

## Recommended cfRNA Tracing Workflow

Do not use a single atlas as the final source-of-origin answer. A robust brain injury cfRNA workflow should combine references in stages:

1. Use GTEx v8 median TPM to filter non-brain peripheral tissue signals and common blood/immune/liver/lung/kidney/muscle background genes.
2. Use Siletti, Chen Brain Cell Atlas, and Hodge/Allen cortex references to score brain cell-type signatures.
3. Use Allen Human Brain Atlas and HPA brain region expression to estimate brain-region origin.
4. Overlay Garza TBI and Allen Aging/Dementia/TBI signatures to detect injury-state, chronic injury, aging, or neurodegeneration-related programs.
5. Report agreement and disagreement across atlases rather than forcing one atlas-specific label.

## Manual or Controlled Access Notes

Some portals expose metadata through JavaScript interfaces, controlled access systems, or interactive download widgets. The downloader records these pages in the manifest as `manual_required` instead of forcing access.

- Brain Cell Atlas may require browser interaction to list the latest dataset files.
- CELLxGENE collections may require manual export/API selection for full h5ad assets.
- NeMO/dbGaP/Synapse or other controlled-access raw reads must be obtained under the relevant data-use agreement and are not downloaded by default.
- Allen/BICCN legacy pages may require manual selection if direct file URLs are not present in page HTML.

## Failure Handling

Each file is downloaded independently. Existing files are skipped. Failed URLs are appended to `data/reference_atlases/99_logs/failed_downloads.tsv` with timestamp, atlas name, URL, file name, and error. The script continues to the next file and still writes manifest outputs, so a stale URL does not block the whole reference-library build.

## Next Run Order

Recommended processing:

```bash
python scripts/build_processed_reference_v1.py
python scripts/build_marker_candidates_v1.py
python scripts/build_sqlite_reference_db_v1.py
```

QC:

```bash
python scripts/qc_reference_downloads.py
```

Optional download repair:

```bash
python scripts/download_reference_atlases.py --source siletti
python scripts/download_reference_atlases.py --source hodge2019
python scripts/download_reference_atlases.py --source brain_cell_atlas
```

Large-file downloads:

```bash
python scripts/download_reference_atlases.py --source siletti --large
python scripts/download_reference_atlases.py --source gtex --large
python scripts/download_reference_atlases.py --source hpa_brain --large
```
