# cfRNA-BrainTrace

cfRNA-BrainTrace is a Python and Streamlit application for hierarchical
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
- Command-line interface: `cli.py`
- Locked Network inference: `core/network_tracing.py`
- Versioned models: `data/models/`
- Validation and export scripts: `scripts/`
- Tests: `tests/`

## Installation

```bash
python -m venv .venv
.venv\Scripts\activate
pip install -r requirements.txt
pip install -e .
```

## Run the application

```bash
streamlit run streamlit_app.py
```

## Command-line use

```bash
cfrna-tracing --help
```

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

This repository contains code, lightweight model artifacts, tests,
documentation and manuscript assets. It intentionally excludes:

- patient-level MRI and clinical data;
- TCGA, GEO, AHBA and other downloaded expression matrices;
- the local SQLite database;
- large generated result directories;
- third-party executables and vendored environments.

Users must obtain source datasets from their original repositories under the
applicable access and licensing conditions. External TPM-like inputs are used
in a cross-scale correlation stress test; `log1p(TPM)` is not described as a
full conversion to Bo2023 VSD.

## Validation summary

Current submission route:

- Network projected-VSD LOSO Top1/Top3: 58.00% / 91.58%.
- Network projected-VSD LOMO Top1/Top3: 53.72% / 91.33%.
- Locked three-tier route LOSO Network Top1/Top3: 58.24% / 92.19% (`n=819`).
- Locked three-tier route LOMO Network Top3: 91.21%.
- Resolution-group Top3 LOSO/LOMO: 72.36% (`n=814`) / 69.09% (`n=812`).
- Exact-region Top3 LOSO/LOMO: 45.33% (`n=814`) / 42.36% (`n=812`).
- AHBA mapped-label Network Top1/Top3: 74.68% / 94.42%.
- AHBA resolution-group Top1/Top3: 36.26% / 67.03%.
- AHBA exact-region Top1/Top3: 24.18% / 42.86%.
- TCGA/BraTS Network Top3: 40.00%.
- TCGA/BraTS broad-anatomy Top3: 64.62%.

Earlier baseline routes were used during development but are not part of the
current submission route or reported validation results.

Network metrics include all 819 samples. Region-level metrics are restricted
to folds in which the held-out truth region remains represented in the
training reference; five LOSO samples and seven LOMO samples are therefore
excluded only from resolution-group and exact-region evaluation.

See `manuscript/` for the Bioinformatics Application Note draft and
supplementary material.

## Status

This repository contains the v0.1.6 public submission release for the
Bioinformatics Application Note describing cfRNA-BrainTrace. The software is
intended for research use in hierarchical brain-origin candidate ranking and
resolution-limit auditing. It is not a clinical diagnostic device and does not
provide stand-alone clinical localization from unlabeled biofluid RNA. The
manuscript-associated v0.1.6 release is archived at Zenodo with version DOI
`https://doi.org/10.5281/zenodo.20780280` and concept DOI
`https://doi.org/10.5281/zenodo.20773674`.
