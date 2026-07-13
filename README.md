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
- Locked Network inference: `core/network_tracing.py`
- Versioned models: `data/models/`
- Validation and export scripts: `scripts/`
- Tests: `tests/`

## Installation

```bash
python -m venv .venv
.venv\Scripts\activate
pip install -r requirements.txt
```

## Run the application

```bash
streamlit run streamlit_app.py
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

## Validation summary

Latest aligned locked route (rechecked 2026-07-13):

- Internal LOSO Network Top1/Top3: 58.24% (`477/819`) / 92.19% (`755/819`).
- Internal LOMO Network Top1/Top3: 57.75% (`473/819`) / 91.21% (`747/819`).
- Internal LOSO resolution-group Top1/Top3: 44.47% (`362/814`) / 72.36% (`589/814`).
- Internal LOMO resolution-group Top1/Top3: 41.38% (`336/812`) / 69.09% (`561/812`).
- Internal LOSO exact-region Top1/Top3: 22.48% (`183/814`) / 45.33% (`369/814`).
- Internal LOMO exact-region Top1/Top3: 22.17% (`180/812`) / 42.36% (`344/812`).
- AHBA technical-replicate-collapsed, network-qualified mapped-label Network Top1/Top3: 73.99% (`165/223`) / 94.62% (`211/223`).
- AHBA technical-replicate-collapsed, network-qualified mapped-label resolution-group Top1/Top3: 40.91% (`36/88`) / 65.91% (`58/88`).
- AHBA technical-replicate-collapsed, network-qualified mapped-label exact-region Top1/Top3: 23.86% (`21/88`) / 46.59% (`41/88`).
- TCGA/BraTS Network Top3: 40.63% (`26/64` evaluable cases).
- TCGA/BraTS broad-anatomy Top3: 65.63% (`42/64` evaluable cases).
- GSE189919 projector gene overlap: 15,622 / 21,668 (72.10%); projection feasibility / transfer stress test only, not localization accuracy.

The repository exposes only the locked three-tier submission route. Historical
baseline and tumour-adapted routes are not distributed in this release.

The aligned endpoint definitions use projected VSD only for the Network Top3
beam, an independent logCPM/Fisher-Top200 ranking for resolution groups, and an
independent `0.25 x Top50 + 0.75 x Top100` logCPM fusion for exact regions.
Resolution-group ordering does not rewrite the exact-region ranking. The
previous AHBA resolution-group Top3 value was produced by deriving group order
from the exact-region fusion and is superseded by the aligned independent-group
result reported above.

Network metrics include all 819 samples. Region-level metrics are restricted
to folds in which the held-out truth region remains represented in the
training reference; five LOSO samples and seven LOMO samples are therefore
excluded only from resolution-group and exact-region evaluation.

The 92.19% LOSO Network Top3 value uses all 819 Network-evaluable samples as
the denominator. Region-level LOSO metrics use 814 reference-supported samples
because five samples lacked a truth-region reference after fold construction.
Region-level LOMO metrics use 812 reference-supported samples because seven
samples lacked a truth-region reference after fold construction.

## Status

This repository contains the v0.1.7 revised public submission release for the
Bioinformatics Application Note describing cfRNA-BrainTrace. The software is
intended for research use in hierarchical brain-origin candidate ranking and
resolution-limit auditing. It is not a clinical diagnostic device and does not
provide stand-alone clinical localization from unlabeled biofluid RNA. The
v0.1.7 supersedes v0.1.6; its Zenodo version DOI is assigned when the release
archive is published. The persistent Zenodo concept DOI is
`https://doi.org/10.5281/zenodo.20773674`.
