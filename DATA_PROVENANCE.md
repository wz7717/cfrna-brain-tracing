# Data Provenance — BrainTrace v0.1.10

This manifest documents every external dataset used to build, validate, and
benchmark the BrainTrace hierarchical brain-origin candidate ranking tool
(v0.1.10; persistent Zenodo concept DOI:
[10.5281/zenodo.20773674](https://doi.org/10.5281/zenodo.20773674)).
All datasets are publicly available under the cited accessions and licensing
terms. No patient-level clinical identifiers are redistributed.

## Source Datasets

| # | Dataset | Accession / DOI | Samples | Role |
|---|---------|----------------|---------|------|
| 1 | Bo2023 Macaque Brain Atlas | [10.1038/s41593-023-01379-0](https://doi.org/10.1038/s41593-023-01379-0) | 819 (9 donors, 110 regions) | Training reference |
| 2 | Allen Human Brain Atlas (AHBA) | [https://human.brain-map.org](https://human.brain-map.org) | 2 donors, 231 mapped samples | Cross-species external validation |
| 3 | TCGA-GBM / TCGA-LGG | [https://portal.gdc.cancer.gov](https://portal.gdc.cancer.gov) | 20 evaluable cases | Glioma domain-shift test |
| 4 | BraTS-TCGA-LGG | [10.5281/zenodo.3718921](https://doi.org/10.5281/zenodo.3718921) | 65 cases (MRI truth) | Anatomical validation |
| 5 | GSE189919 (GEO) | [https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE189919](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE189919) | 72,108 gene panel | Engineering benchmark |
| 6 | Huang2025 cfRNA | [10.1186/s41698-025-00909-0](https://doi.org/10.1186/s41698-025-00909-0) | 159 CSF/plasma profiles | cfRNA domain-shift audit |

## Data Access & Licensing

- **Bo2023**: Published with the original paper; raw expression matrices and
  sample metadata are available from the authors or supplementary material.
  The BrainTrace repository contains derived model artifacts (SHA-256 locked)
  but not the raw expression matrices.
- **AHBA**: Available via the Allen Institute API; requires acceptance of
  Allen Institute terms of use. BrainTrace distributes only aggregated
  validation results, not raw AHBA expression data.
- **TCGA**: Available through the NCI Genomic Data Commons; requires dbGaP
  authorization for controlled-access tiers. BrainTrace uses only open-access
  gene expression data.
- **BraTS-TCGA-LGG**: Zenodo-hosted imaging dataset (CC BY 4.0). The nested
  NIfTI zip is automatically extracted by `reproduce_all.py` if present.
- **GSE189919**: GEO public dataset; no access restrictions.
- **Huang2025**: Supplementary material distributed with the published article
  (PMC12041490); no additional access restrictions.

## Immutable Snapshot Policy

The validation results archived in `reproducibility/v4_p*.csv` and summarized
in the manuscript correspond to frozen snapshots of these datasets. Users
re-running `reproduce_all.py` must obtain the same data versions from their
original repositories to reproduce the exact numerical results. The software
version v0.1.10 is archived as a distinct version under the Zenodo concept DOI above.

## SHA-256 Integrity

The `data/models/canonical110_model_lock.json` manifest records SHA-256
checksums for all locked model artifacts. Production inference verifies these
hashes at load time and fails closed on any mismatch. See `MODEL_LOCK.md` for
the complete lock specification.
