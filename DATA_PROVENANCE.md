# Data Provenance — BrainTrace v0.1.10

This manifest documents every external dataset used to build, validate, and
benchmark the BrainTrace hierarchical brain-origin candidate ranking tool
(v0.1.10, version DOI: [10.5281/zenodo.21756113](https://doi.org/10.5281/zenodo.21756113); concept DOI: [10.5281/zenodo.20773674](https://doi.org/10.5281/zenodo.20773674)).
All datasets are publicly available under the cited accessions and licensing
terms. No patient-level clinical identifiers are redistributed.

## Source Datasets

| # | Dataset | Accession / DOI | Samples | Role |
|---|---------|----------------|---------|------|
| 1 | Bo2023 Macaque Brain Atlas | [10.1038/s41593-023-01379-0](https://doi.org/10.1038/s41593-023-01379-0) | 819 (9 donors, 110 regions) | Training reference |
| 2 | Allen Human Brain Atlas (AHBA) | [https://human.brain-map.org](https://human.brain-map.org) | 2 donors, 231 mapped samples | Cross-species external validation |
| 3 | TCGA-GBM / TCGA-LGG | [NCI Genomic Data Commons](https://portal.gdc.cancer.gov) | 65 expression cases linked to the imaging cohort; 64 edema-evaluable after one cerebellar case was excluded as out of scope | Glioma domain-shift test |
| 4 | BraTS-TCGA-LGG | [10.5281/zenodo.3718921](https://doi.org/10.5281/zenodo.3718921) | 65 MRI cases; 64 edema-evaluable | Imaging-derived anatomical truth for the linked TCGA cases |
| 5 | GSE189919 (GEO) | [GSE189919](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE189919) | 51 samples; 72,108 expression rows; 15,622/21,668 frozen-projector genes overlap | Engineering benchmark and input-domain audit |
| 6 | Huang2025 cfRNA | [10.1186/s41698-025-00909-0](https://doi.org/10.1186/s41698-025-00909-0) | 159 CSF/plasma profiles | cfRNA domain-shift audit |
| 7 | Gene Ontology, KEGG, and g:Profiler | [Gene Ontology](https://geneontology.org), [KEGG](https://www.kegg.jp), [g:Profiler](https://biit.cs.ut.ee/gprofiler) | Frozen 200-gene panel; 179 mapped by g:Profiler; 21,668-gene model-space background | Independent GO:BP/KEGG annotation; not model training or predictive validation |
| 8 | Chiou2023 rhesus macaque single-cell atlas | [10.1126/sciadv.adh1914](https://doi.org/10.1126/sciadv.adh1914) | Published Tables S3/S7; seven prespecified broad marker families | Primary independent cell-type annotation-bias analysis |
| 9 | Siletti2023 adult human brain cell atlas | [10.1126/science.add7046](https://doi.org/10.1126/science.add7046) | Published cluster annotations; the same seven broad marker families | Human-reference sensitivity analysis |

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
- **GO/KEGG/g:Profiler**: Database content is versioned and may change. The
  submission calculation package records g:Profiler version
  `e114_eg62_p19_27110d83`, query date 2026-07-31, requested sources, mapped
  identifiers, background, and response metadata.
- **Chiou2023/Siletti2023**: BrainTrace redistributes derived marker-family
  intersections and enrichment summaries only. The source spreadsheets remain
  subject to their publishers' and authors' terms.

## Frozen Labels and Analysis Roles

The formal Bo2023 scoring truth is frozen before evaluation. Two audited label
normalizations are applied before scoring: sample `19R470` is normalized from
the raw Parietal label to `Occipital/Temporal`, and sample
`xz-L-10m-10o_BRRL200004691-1A` is normalized from raw Lateral Prefrontal Cortex
to `Orbitomedial Prefrontal Cortex (OMPFC)`. These are anatomical crosswalk
corrections, not post-hoc changes based on predictions.

The fixed comparator uses exactly the locked 200 Network genes and leave-one-
macaque-out evaluation. Its row-level result file contains 819 samples and
sums to 389 Top1 and 680 Top3 hits. Historical k=500 RF results are retained
only as a separately labelled legacy comparator.

The 3′ analysis is a transcript-coordinate/detection proxy sensitivity test;
it is not molecule-level read generation and must not be described as a true
3′-end sequencing simulation. GO/KEGG and cell-type analyses annotate the
frozen panel and do not establish cell of origin, mechanism, causality, or
clinical predictive performance.

## Immutable Snapshot Policy

The validation results archived in `reproducibility/v4_p*.csv` and summarized
in the manuscript correspond to frozen snapshots of these datasets. Users
re-running `reproduce_all.py` must obtain the same data versions from their
original repositories to reproduce the exact numerical results. The software
version v0.1.10 is permanently archived at Zenodo under version DOI
10.5281/zenodo.21756113. The concept DOI resolves to the latest Zenodo version.
Post-archive metadata corrections and newly public reproducibility assets are
part of a later repository state unless a new immutable release archive
explicitly includes them.

## SHA-256 Integrity

The `data/models/canonical110_model_lock.json` manifest records SHA-256
checksums for all locked model artifacts. Production inference verifies these
hashes at load time and fails closed on any mismatch. See `MODEL_LOCK.md` for
the complete lock specification.

The frozen 200-gene RF row-level result has SHA-256
`d27b67c5ca7fa79b9186a79b24d94cf2624c4007c95bb53016043f554ceba3ac`.
The independent enrichment manifest records the input SHA-256 values for the
locked panel, 21,668-gene projector background, and the Chiou2023 and
Siletti2023 source spreadsheets. Derived enrichment outputs, marker-family
intersections, protocol and manifest are public under
`reproducibility/independent_enrichment/`; publisher-hosted source workbooks
are not redistributed. The Saleem crosswalk is public under
`reproducibility/crosswalks/`, and the portable S18-S19 runner and input hashes
are public under `reproducibility/s18_s19/`.
