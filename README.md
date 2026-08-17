# BrainTrace

BrainTrace is a Python and Streamlit application for hierarchical
brain-origin candidate inference from RNA expression profiles using a primate
transcriptomic reference.

The current submission workflow is:

```text
Query logCPM/logTPM-compatible expression
-> reference-fitted projection to Bo2023-like VSD
-> fixed 10-class Network Top3 candidate filtering
-> logCPM-compatible resolution-group and exploratory exact-region reranking
```

Projected VSD is used only for broad Network candidate filtering. Downstream
resolution-group and exploratory exact-region candidates are reranked within
the retained Network candidate set using logCPM-compatible local expression. Exact-region
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
- Command-line interface: `core/cli.py`
- Shared locked production inference: `core/inference.py`
- Locked Network inference: `core/network_tracing.py`
- Versioned models: `data/models/`
- Validation and export scripts: `scripts/`
- Tests: `tests/`

## Quick tutorial

### Option 1 — Live demo

1. Open <https://brain-cfrna-tracing.streamlit.app/>.
2. On **Tracing Analysis**, click **Load example data**.
3. Click **Run example**.
4. Inspect Network Top3, Resolution Group Top3 and Exact-Region Exploratory
   Top3 together with coverage, entropy and score margin.

The loaded file is a synthetic software smoke test, not a biological sample
or localization-validation dataset. It is not saved to Sample Management or
SQLite.

### Option 2 — CLI

```bash
pip install -e .
braintrace validate
braintrace query \
  --input examples/braintrace_example_counts.tsv \
  --output example_result.json
python examples/verify_examples.py
```

### Interpreting output

- **Network Top3** is the primary validated candidate tier.
- **Resolution Group Top3** is the recommended resolvable candidate tier.
- **Exact-Region Top3** is exploratory and is not a deterministic localization
  call.

## Installation and entry points

Choose one of the following workflows. The application environment is the
smallest supported installation for interactive inference; it is not the full
manuscript-reproduction environment.

### Docker (recommended for quick start)

```bash
docker build -t braintrace:v0.1.16 .
docker run -p 8501:8501 braintrace:v0.1.16
# Open http://localhost:8501 in your browser

# CLI mode
docker run --rm --entrypoint braintrace braintrace:v0.1.16 --help

# CLI query with the current directory mounted at /work
docker run --rm --entrypoint braintrace \
  -v "$PWD:/work" \
  braintrace:v0.1.16 \
  query --input /work/sample_counts.tsv --output /work/result.json
```

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

#### Requirements files at a glance

| File | Purpose | Packages | Install for |
|------|---------|----------|-------------|
| `requirements.txt` | Core application (Streamlit + inference) | 8 (pinned) | Interactive use |
| `requirements_reproducible.txt` | Full manuscript reproduction | 24 (pinned) | Regenerating all results |
| `requirements-repro.txt` | Testing extras | includes reproducible + pytest | CI and development |
| `environment.yml` | Conda mirror of reproducible | 24 | Conda users |

### Command-line interface

The package installs the `braintrace` console entry point:

```bash
pip install -e .

# Show help
braintrace --help

# Run a single-sample query from the command line
braintrace query --input sample_counts.tsv --output result.json

# List available models
braintrace models

# Verify the frozen production model bundle
braintrace validate
```

Both `braintrace query` and the Streamlit web application call the same locked
production function in `core/inference.py`. The web application
(`streamlit run app/main.py`) provides the complete released three-tier route:
Network Top3, Resolution Group Top3 and Exact-Region Exploratory Top3. Exact
regions remain exploratory candidates, not validated localization calls.

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
generation, followed by logCPM-based resolution-group and exploratory
exact-region reranking using the frozen region-reference artifacts distributed
with the release.

The non-redistributed Bo2023 author-package matrices are required only for
rebuilding the frozen reference artifacts and reproducing
reference-construction analyses from source; they are not required for
inference with the released model bundle.

TPM or logTPM tables remain accepted for backward compatibility, but they are
treated as fallback inputs and fine-region interpretations should be reported
more cautiously. Users are not expected to upload VSD; VSD-like query
expression is generated internally by the projector.

**Gene identifier format**: The tool accepts HGNC gene symbols as the primary
identifier. Macaque Ensembl gene IDs (prefix `ENSMFAG`) are also recognized
when present in the frozen reference lookup. Entrez Gene IDs and Ensembl
human gene IDs are not directly supported; users should convert these to
HGNC symbols before upload. Unrecognized gene symbols are set to zero in the
projection step and do not contribute to ranking.

Optional sample, subject, diagnosis and anatomical metadata can be included and
are retained in exported reports.

## Reproducibility and data policy

Released-model inference reproducibility requires the public code and frozen
model artifacts distributed with this release. Full reference-construction
reproduction additionally requires the external or non-redistributed source
datasets identified in [`DATA_PROVENANCE.md`](DATA_PROVENANCE.md). These are
distinct reproducibility levels: absence of the author-package matrices limits
rebuilding from source, not inference with the released frozen bundle.

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
- `reproducibility/round5_analysis/`, containing the Network-layer direct
  projection ablation, projector OLS-fit audit, LOSO/LOMO MRR and NDCG@3
  summaries, and the P0-3 donor-Network pseudobulk DESeq2 audit of the frozen
  Network Top200 panel, plus the P0-4 same-version three-background GO:BP/KEGG
  stability audit, with protocols and hash manifests.

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
- TCGA/BraTS Network Top3 (strict, current tracer): 23.81% (`15/63` primary edema-comparator cases; exploratory one-sided exact-binomial p=0.8888 vs the 30% uniform-network null).
- TCGA/BraTS broad-anatomy Top3 (strict, current tracer): 82.54% (`52/63` primary edema-comparator cases; descriptive only, with no valid prespecified Top3 null).
- TCGA/BraTS Network candidate-set any-hit Top3: 36.51% (`23/63`; descriptive sensitivity; candidate set is Networks with >=20% edema overlap, not anatomical adjacency).
- The primary edema comparator excludes TCGA-HT-7686 (no label-2 edema voxels) and TCGA-HT-7680 (cerebellar/out-of-scope edema); other MRI truth bases retain n=65.
- The former `20/64` and `51/64` values belong to the archived 2026-06-09 endpoint and are not current release metrics.
- The current reproducibility endpoint is generated by `score_tcga_gbm_lgg_sample_tracing_with_mri_labels.py` and `evaluate_brats_tcga_lgg_65_mri_truth.py` under `results/tcga_brats_current/`; the archived 2026-06-09 files are historical inputs only.
- GSE189919 projector gene overlap: 15,622 / 21,668 (72.10%); projection feasibility / transfer stress test only, not localization accuracy.

The repository exposes only the locked three-tier submission route. Historical
baseline and tumour-adapted routes are not distributed in this release.

The aligned endpoint definitions use projected VSD only for the fixed Network
Top3 candidate filter, an independent logCPM/Fisher-Top200 ranking for resolution groups, and an
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

## Bo2023 source-data placement

Bo2023 raw sequencing is available under NCBI SRA BioProject `PRJNA905082`.
The exact processed count, VSD and metadata files used for the frozen projector
were obtained from the paper's released supplementary/author package; they are
not redistributed here and are not implied to be reproducible from SRA without
the authors' processing workflow. Place them at the paths below before running
the full reproduction pipeline. The checksum gate rejects non-matching files.

| Expected relative path | SHA-256 |
|---|---|
| `external_data/Bo2023/mfas5_819samples_28415genes_featurecounts_counts.txt` | `1fb3a512da11ab0c327c07c114da3b9c38cab0a504682f2c7c036eedb3c7561a` |
| `external_data/Bo2023/mfas5_819samples_23605genes_vsd4_rmbatch.xls` | `286aeab66b21b7fa012fac8ceaa24497894327e0736f9f6b200334c57089a1b3` |
| `external_data/Bo2023/Information of sequenced samples_update_full878_filter819.xlsx` | `9a2fe2bec1475f6ad613883d0ff5925b1e6ba36e800caa922c35d4f8ae7d3645` |

### Public Bo2023 file-level links

The publisher's direct file links are recorded here so that the public source
package is not confused with the exact author-package inputs above:

- Supplementary Data 1–15 ZIP: <https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41467-023-37246-w/MediaObjects/41467_2023_37246_MOESM3_ESM.zip>
- Source Data workbook: <https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41467-023-37246-w/MediaObjects/41467_2023_37246_MOESM5_ESM.xlsx>
- Article landing page: <https://www.nature.com/articles/s41467-023-37246-w>
- SRA BioProject: <https://www.ncbi.nlm.nih.gov/bioproject/PRJNA905082>
- Analysis-code archive cited by Bo et al.: <https://doi.org/10.5281/zenodo.7641873>

The publisher ZIP contains related supplementary workbooks (including a
Supplementary Data 1 metadata workbook), but it does not expose the two exact
SHA-256-matched count/VSD filenames listed above. No stable public file URL for
those processed author-package matrices was identified in the publisher or SRA
records checked on 2026-08-10; obtain them from the corresponding authors and
verify the hashes rather than substituting a similarly named workbook.

## Availability and maintenance

The authors commit to maintaining the public BrainTrace web service for at
least three years after publication. Versioned source code and archived
software releases will remain available through GitHub and Zenodo.

The software source is open under the MIT license, and archived versions remain
citable. Web-service availability is supplementary to the archived executable
code: a temporary service interruption does not remove the versioned CLI,
source, examples or immutable archives.

## Continuous-integration coverage

CI runs the full test suite with branch coverage for the `core/` package. The
v0.1.16 submission-readiness run produced **120 passed, 2 skipped**. Coverage.py 7.14.3
reported **71.55% line coverage**, **53.93% branch coverage**, and **68% combined
terminal coverage** (1,406 statements; 382 branches). The exact text and XML
reports are `reproducibility/coverage_report.txt` and
`reproducibility/coverage.xml`; the workflow also uploads the XML report as a
CI artifact. This percentage is a software-test coverage metric, not a
statistical confidence or biological validation measure.

## Status

The current executable software metadata is BrainTrace v0.1.16, a
documentation, public-example and web-usability patch prepared for the
Bioinformatics Application Note. It does not change the frozen model,
ontology, formal prediction set, Network Top1/Top3, resolution-group or
exact-region endpoints. The production model remains locked under
`canonical110-v0.1.12-20260813`.

No model artifact, learned parameter, anatomical ontology, formal prediction,
primary endpoint or benchmark result changed in v0.1.16. The
software is intended for research use in hierarchical brain-origin candidate
ranking and resolution-limit auditing. It is not a clinical diagnostic device
and does not provide stand-alone clinical localization from unlabeled biofluid
RNA.

BrainTrace v0.1.15 remains the previous immutable GitHub/Zenodo release under version
DOI `https://doi.org/10.5281/zenodo.21970252`. The frozen v0.1.12 scientific
release remains archived under `https://doi.org/10.5281/zenodo.21911532`; the
v0.1.14 historical software record remains at
`https://doi.org/10.5281/zenodo.21920261`; the persistent Zenodo concept DOI is
`https://doi.org/10.5281/zenodo.20773674`.

v0.1.16 release metadata will be synchronized with its immutable archive after
archival. No v0.1.16 version DOI is claimed before that archive exists.

The official Zenodo v0.1.15 software record is
`https://doi.org/10.5281/zenodo.21970252`. Because Zenodo's GitHub-generated
source archive does not materialize Git LFS objects, the v0.1.15 full
reproducibility archive, including the unchanged 164 MB payload, is deposited
separately under `https://doi.org/10.5281/zenodo.21970278`. The same payload is
also
available as a checksum-matched asset on the
[GitHub v0.1.15 release](https://github.com/wz7717/cfrna-brain-tracing/releases/tag/v0.1.15).

The v0.1.14 LOMO Network F1 discrepancy arose from stale rounded class-level
recall values and a non-formal historical route being used to estimate summary
metrics. The v0.1.15 values are regenerated from the same locked formal
prediction-level source: Top1 remains 455/819, and micro-F1 is exactly 455/819.

The 164 MB `permutation_fg_max.npz` reproducibility intermediate is stored with
Git LFS. Clone with Git LFS enabled (or run `git lfs pull`) to retrieve its full
contents; its SHA-256 is included in the release checksum manifests.
