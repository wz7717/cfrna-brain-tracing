# BrainTrace

BrainTrace is a Python and Streamlit application for hierarchical
brain-origin candidate inference from RNA expression profiles using a primate
transcriptomic reference.

The current submission workflow is:

```text
Query logCPM/logTPM-compatible expression
-> reference-fitted projection to Bo2023-like VSD
-> 10-class Network Top3 beam generation
-> logCPM-compatible resolution-group and exploratory exact-region reranking
```

Projected VSD is used only for broad Network beam generation. Downstream
resolution-group and exploratory exact-region candidates are reranked within
the retained Network beam using logCPM-compatible local expression. Exact-region
output is a candidate ranking, not a deterministic localization endpoint.

The software reports candidates at three main levels:

1. 10-class anatomical Network;
2. resolution group;
3. exploratory exact region.

Current evidence supports coarse candidate ranking more strongly than exact
localization. Biofluid cohorts without patient-level anatomical truth are
treated as external transfer stress tests, not localization-accuracy
validation.

## Interfaces

- Streamlit application: `streamlit_app.py`
- Locked Network inference: `core/network_tracing.py`
- Versioned models: `data/models/`
- Validation and export scripts: `scripts/`
- Tests: `tests/`

## Installation and entry points

Choose one of the following workflows. The application environment is the
smallest supported installation for interactive inference; it is not the full
manuscript-reproduction environment.

### Run the application

```bash
python -m venv .venv
# Windows PowerShell: .venv\Scripts\Activate.ps1
# macOS/Linux: source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
python -m streamlit run streamlit_app.py
```

### Fully reproduce the paper

Obtain the external source datasets and place them at the paths documented in
[`DATA_PROVENANCE.md`](DATA_PROVENANCE.md) and `SHA256SUMS.txt`. Then run the
checksum gate followed by the complete build, validation and CSV-generation
pipeline:

```bash
python -m venv .venv-repro
# Windows PowerShell: .venv-repro\Scripts\Activate.ps1
# macOS/Linux: source .venv-repro/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements_reproducible.txt
python reproduce_all.py --verify-only
python reproduce_all.py --output-dir reproducibility_audit
```

`requirements-repro.txt` adds testing extras to the complete reproducibility
environment; it is intended for CI/development checks rather than as the main
installation entry point.

## Input format

Recommended query input is a gene-level RNA-seq expression table with:

```text
gene_symbol
raw_counts
```

or a pre-normalized logCPM table:

```text
gene_symbol
logCPM
```

Raw counts are converted internally to logCPM. The validated route uses
logCPM-derived query expression with projected VSD for Network-level candidate
generation, followed by logCPM-based resolution-group and local reranking where
the full private reference is available.

TPM or logTPM tables remain accepted for backward compatibility, but they are
treated as fallback inputs and fine-region interpretations should be reported
more cautiously. Users are not expected to upload VSD; VSD-like query
expression is generated internally by the projector.

Optional sample, subject, diagnosis and anatomical metadata can be included and
are retained in exported reports.

## Reproducibility and data policy

This deployment repository contains application code, lightweight model
artifacts and tests. It intentionally excludes:

- patient-level MRI and clinical data;
- TCGA, GEO, AHBA and other downloaded expression matrices;
- the local SQLite database;
- large generated result directories;
- third-party executables and vendored environments.

Users must obtain source datasets from their original repositories under the
applicable access and licensing conditions. External TPM-like inputs are used
in a cross-scale correlation stress test; `log1p(TPM)` is not described as a
full conversion to Bo2023 VSD.

Public derived reproducibility assets include:

- `reproducibility/crosswalks/P2_Bo2023_Saleem_crosswalk.csv`, the complete
  110-region locked-Network/Saleem-style crosswalk;
- `reproducibility/independent_enrichment/`, containing the frozen GO:BP,
  KEGG and broad cell-type enrichment outputs, protocol and hash manifest;
- `reproducibility/s18_s19/`, containing the portable S18-S19 runner, frozen
  output and portable input manifest.

Run the S18-S19 analysis from the repository root with:

```bash
python reproducibility/s18_s19/run_s18_s19_sensitivity.py
```

Publisher-hosted source workbooks and raw expression matrices remain excluded;
the public enrichment directory contains only derived marker sets and results.

## Validation summary

Canonical 110 locked route (rechecked 2026-07-13):

- The Bo2023 reference retains the paper's 110 post-QC anatomical region IDs. Each region has one canonical parent Network; the two discordant assay-level labels are normalized as `10m -> Orbitomedial Prefrontal Cortex (OMPFC)` and `V2 -> Occipital/Temporal` without modifying the source workbook.
- Internal LOSO Network Top1/Top3: 58.97% (`483/819`) / 91.94% (`753/819`).
- Internal LOMO Network Top1/Top3 on the frozen no-pairwise route: 55.56% (`455/819`) / 91.58% (`750/819`). The former 59.10% (`484/819`) Top1 value used pairwise rescue and is deprecated for the locked production route.
- Internal LOSO resolution-group Top1/Top3: 45.21% (`368/814`) / 72.48% (`590/814`).
- Internal LOMO resolution-group Top1/Top3: 42.36% (`344/812`) / 70.07% (`569/812`).
- Internal LOSO exact-region Top1/Top3: 22.36% (`182/814`) / 45.21% (`368/814`).
- Internal LOMO exact-region Top1/Top3: 21.80% (`177/812`) / 42.61% (`346/812`).
- AHBA technical-replicate-collapsed, network-qualified mapped-label Network Top1/Top3: 73.99% (`165/223`) / 94.62% (`211/223`).
- AHBA technical-replicate-collapsed, network-qualified mapped-label resolution-group Top1/Top3: 42.05% (`37/88`) / 68.18% (`60/88`).
- AHBA technical-replicate-collapsed, network-qualified mapped-label exact-region Top1/Top3: 27.27% (`24/88`) / 45.45% (`40/88`).
- TCGA/BraTS Network Top3 (strict): 31.25% (`20/64` evaluable cases, binomial p=0.4602 vs 30% uniform null).
- TCGA/BraTS broad-anatomy Top3: 79.69% (`51/64` evaluable cases, binomial p<0.001 vs 20% uniform null).
- GSE189919 projector gene overlap: 15,622 / 21,668 (72.10%); projection feasibility / transfer stress test only, not localization accuracy.

The repository exposes only the locked three-tier submission route. Historical
baseline and tumour-adapted routes are not distributed in this release.

The aligned endpoint definitions use projected VSD only for the Network Top3
beam, an independent logCPM/Fisher-Top200 ranking for resolution groups, and an
independent `0.25 x Top50 + 0.75 x Top100` logCPM fusion for exact regions.
Resolution-group ordering does not rewrite the exact-region ranking. The static
deployment hierarchy does not define validation groups: LOSO/LOMO groups are
fold-local, and AHBA groups are constructed from the locked macaque training
reference.

Network metrics include all 819 samples. Region-level metrics are restricted
to folds in which the held-out truth region remains represented in the
training reference; five LOSO samples and seven LOMO samples are therefore
excluded only from resolution-group and exact-region evaluation.

The 91.94% LOSO Network Top3 value uses all 819 Network-evaluable samples as
the denominator. Region-level LOSO metrics use 814 reference-supported samples
because five samples lacked a truth-region reference after fold construction.
Region-level LOMO metrics use 812 reference-supported samples because seven
samples lacked a truth-region reference after fold construction.

## Status

This repository contains the v0.1.10 revised public submission release for the
Bioinformatics Application Note describing BrainTrace. The software is
intended for research use in hierarchical brain-origin candidate ranking and
resolution-limit auditing. It is not a clinical diagnostic device and does not
provide stand-alone clinical localization from unlabeled biofluid RNA. The
v0.1.10 release uses the audited canonical 110-region assets, enforces the exact
110-region / 10-Network / 120-beam production contract, cleans the AHBA
resolution-group detail export, and reports the formally rerun validation
metrics. v0.1.10 supersedes v0.1.9 and includes the fourth-round reviewer
analyses, synchronized reproducibility metadata, and automated CI checks. It is
archived at Zenodo with version DOI
`https://doi.org/10.5281/zenodo.21756113`. The persistent Zenodo concept DOI is
`https://doi.org/10.5281/zenodo.20773674`.
