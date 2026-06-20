# cfRNA-BrainTrace

cfRNA-BrainTrace is a Python and Streamlit application for hierarchical
brain-origin candidate inference from RNA expression profiles using a primate
transcriptomic reference.

The locked workflow is:

```text
Bo2023 VSD reference
-> fold-selected 200 genes
-> Pearson correlation
-> Top-3 pairwise rescue
```

The software reports candidates at four levels:

1. lobe;
2. broad anatomy;
3. 10-class anatomical Network;
4. exploratory exact region.

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

- Bo2023 strict LOSO Network Top1/Top3: 55.8% / 88.0%.
- Bo2023 leave-one-monkey-out Top1/Top3: 53.2% / 86.7%.
- AHBA normal-human Network Top1/Top3: 32.6% / 55.4%.
- Paired TCGA-LGG/BraTS broad-anatomy Top3 strict/tolerant:
  75.4% / 83.1%.

See `manuscript/` for the Bioinformatics Application Note draft and
supplementary material.

## Status

This is a research software release candidate. It is not a clinical diagnostic
device. The public repository URL, archived release DOI and OSI-approved
license must be finalized before journal submission.
