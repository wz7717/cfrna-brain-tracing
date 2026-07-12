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

- Locked three-tier route LOSO Network Top1/Top3: 58.24% / 92.19% (`n=819`).
- Locked three-tier route LOMO Network Top3: 91.21%.
- Resolution-group Top3 LOSO/LOMO: 72.36% (`n=814`) / 69.09% (`n=812`).
- Exact-region Top3 LOSO/LOMO: 45.33% (`n=814`) / 42.36% (`n=812`).
- AHBA technical-replicate-collapsed, network-qualified mapped-label Network Top1/Top3: 73.99% / 94.62% (`n=223`).
- AHBA technical-replicate-collapsed, network-qualified mapped-label resolution-group Top1/Top3: 40.91% / 62.50% (`n=88`).
- AHBA technical-replicate-collapsed, network-qualified mapped-label exact-region Top1/Top3: 23.86% / 46.59% (`n=88`).
- TCGA/BraTS Network Top3: 40.63% (`26/64` evaluable cases).
- TCGA/BraTS broad-anatomy Top3: 65.63% (`42/64` evaluable cases).
- GSE189919 projector gene overlap: 15,622 / 21,668 (72.10%); projection feasibility / transfer stress test only, not localization accuracy.

The repository exposes only the locked three-tier submission route. Historical
baseline and tumour-adapted routes are not distributed in this release.

Network metrics include all 819 samples. Region-level metrics are restricted
to folds in which the held-out truth region remains represented in the
training reference; five LOSO samples and seven LOMO samples are therefore
excluded only from resolution-group and exact-region evaluation.

The 92.19% LOSO Network Top3 value uses all 819 Network-evaluable samples as
the denominator. Region-level LOSO metrics use 814 reference-supported samples
because five samples lacked a truth-region reference after fold construction.
Region-level LOMO metrics use 812 reference-supported samples because seven
samples lacked a truth-region reference after fold construction.

See `manuscript/` for the Bioinformatics Application Note draft and
supplementary material.

## Status

This repository contains the v0.1.7 revised public submission release for the
Bioinformatics Application Note describing cfRNA-BrainTrace. The software is
intended for research use in hierarchical brain-origin candidate ranking and
resolution-limit auditing. It is not a clinical diagnostic device and does not
provide stand-alone clinical localization from unlabeled biofluid RNA. The
v0.1.7 supersedes v0.1.6; its Zenodo version DOI is assigned when the release
archive is published. The persistent Zenodo concept DOI is
`https://doi.org/10.5281/zenodo.20773674`.
