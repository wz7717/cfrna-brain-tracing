# Data Provenance - BrainTrace v0.1.16

This manifest documents every external dataset used to build, validate, and
benchmark the BrainTrace hierarchical brain-origin candidate ranking tool.
The current executable metadata and immutable software release are v0.1.16,
archived at version DOI
[10.5281/zenodo.21974954](https://doi.org/10.5281/zenodo.21974954). Its full
reproducibility archive, including the materialized 164 MB Git LFS payload, is
archived separately at
[10.5281/zenodo.21974991](https://doi.org/10.5281/zenodo.21974991); the persistent
software concept DOI is
[10.5281/zenodo.20773674](https://doi.org/10.5281/zenodo.20773674). The frozen
v0.1.12 scientific release remains archived at version DOI
[10.5281/zenodo.21911532](https://doi.org/10.5281/zenodo.21911532).
Public source data are available under the cited accessions except where
explicitly noted below. The exact Bo2023 processed author-package matrices used
by the frozen projector have no identified stable public URL; they must be
obtained from the original authors and verified against the reported SHA-256
hashes. No patient-level clinical identifiers are redistributed.

## Source Datasets

| # | Dataset | Accession / DOI | Samples | Role |
|---|---------|----------------|---------|------|
| 1 | Bo2023 Macaque Brain Atlas | [10.1038/s41467-023-37246-w](https://doi.org/10.1038/s41467-023-37246-w); SRA [PRJNA905082](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA905082) | 819 (9 donors, 110 regions) | Training reference |
| 2 | Allen Human Brain Atlas (AHBA) | [https://human.brain-map.org](https://human.brain-map.org) | 2 donors, 231 mapped samples | Cross-species external validation |
| 3 | TCGA-GBM / TCGA-LGG | [NCI Genomic Data Commons](https://portal.gdc.cancer.gov) | 65 expression cases linked to the imaging cohort; 63 primary edema-comparator cases after excluding TCGA-HT-7686 (no label-2 edema voxels) and TCGA-HT-7680 (cerebellar/out of scope) | Glioma domain-shift test |
| 4 | BraTS-TCGA-LGG | [10.5281/zenodo.3718921](https://doi.org/10.5281/zenodo.3718921) | 65 MRI cases; 63 primary edema-comparator cases under the same exclusions | Imaging-derived anatomical truth for the linked TCGA cases |
| 5 | GSE189919 (GEO) | [GSE189919](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE189919) | 51 samples; 72,108 expression rows; 15,622/21,668 frozen-projector genes overlap | Engineering benchmark and input-domain audit |
| 6 | Huang2025 cfRNA | [10.1038/s41698-025-00909-6](https://doi.org/10.1038/s41698-025-00909-6) | 159 published-matrix profiles (77 CSF, 82 plasma) | Independent-fluid technical portability/domain-shift audit; not localization validation |
| 7 | Gene Ontology, KEGG, and g:Profiler | [Gene Ontology](https://geneontology.org), [KEGG](https://www.kegg.jp), [g:Profiler](https://biit.cs.ut.ee/gprofiler) | Frozen 200-gene panel; 179 mapped by g:Profiler; 21,668-gene model-space background | Independent GO:BP/KEGG annotation; not model training or predictive validation |
| 8 | Chiou2023 rhesus macaque single-cell atlas | [10.1126/sciadv.adh1914](https://doi.org/10.1126/sciadv.adh1914) | Published Tables S3/S7; seven prespecified broad marker families | Primary independent cell-type annotation-bias analysis |
| 9 | Siletti2023 adult human brain cell atlas | [10.1126/science.add7046](https://doi.org/10.1126/science.add7046) | Published cluster annotations; the same seven broad marker families | Human-reference sensitivity analysis |

## Data Access & Licensing

- **Bo2023**: Raw sequencing is deposited under SRA BioProject PRJNA905082.
  The exact processed count, VSD and sample-metadata files used here came from
  the paper's released supplementary/author package. The repository contains
  SHA-256-locked derived model artifacts but does not redistribute those source
  files or claim that the processed matrices can be reconstructed from SRA
  without the authors' processing workflow. For full reproduction, place the
  files at `external_data/Bo2023/` using the exact names and hashes below:

  | File | SHA-256 |
  |---|---|
  | `mfas5_819samples_28415genes_featurecounts_counts.txt` | `1fb3a512da11ab0c327c07c114da3b9c38cab0a504682f2c7c036eedb3c7561a` |
  | `mfas5_819samples_23605genes_vsd4_rmbatch.xls` | `286aeab66b21b7fa012fac8ceaa24497894327e0736f9f6b200334c57089a1b3` |
  | `Information of sequenced samples_update_full878_filter819.xlsx` | `9a2fe2bec1475f6ad613883d0ff5925b1e6ba36e800caa922c35d4f8ae7d3645` |
  Public file-level links are the [Supplementary Data 1–15 ZIP](https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41467-023-37246-w/MediaObjects/41467_2023_37246_MOESM3_ESM.zip), the [Source Data workbook](https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41467-023-37246-w/MediaObjects/41467_2023_37246_MOESM5_ESM.xlsx), the [article landing page](https://www.nature.com/articles/s41467-023-37246-w), the [SRA BioProject](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA905082), and the authors' [Zenodo code archive](https://doi.org/10.5281/zenodo.7641873). The publisher ZIP contains related supplementary workbooks but does not expose the two exact SHA-matched count/VSD author-package filenames; no stable public file URL for those matrices was identified in the publisher or SRA records checked on 2026-08-10. The author package remains required for those inputs and must be hash-verified.
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
  (PMC12041490); no additional access restrictions. The published expression
  matrix contains 159 profiles (77 CSF and 82 plasma). The source article
  reports five CSF and one plasma sequencing-QC exclusions for its own clinical
  analyses, but does not supply a public per-profile QC-status map or a
  patient-level CSF-plasma correspondence. The remediation audit therefore
  uses all 159 published-matrix profiles only as unpaired, fluid-specific
  profile-level technical stress-test observations; patient-level dependence
  cannot be assessed from the public matrix. It does not infer patient identifiers,
  substitute a plasma sample from a CSF sample name, construct synthetic
  mixtures, reproduce the source clinical analysis, or validate localization.
  For the documented full-pipeline path, place the source files at
  `external_inputs/huang2025_pmc12041490/41698_2025_909_MOESM2_ESM.csv`
  (SHA-256 `ef0c72c17d65a0293ec4089880716ca3db1ad74764f43fe3bbe828b3e62ea6a3`)
  and `external_inputs/huang2025_pmc12041490/41698_2025_909_MOESM2_ESM.xlsb`
  (SHA-256 `4ca9519933b75c52b4b6650405971cd167a4ecfaacd8d48d1247f604a800a7df`).
- **GO/KEGG/g:Profiler**: Database content is versioned and may change. The
  public `reproducibility/independent_enrichment/` package records g:Profiler version
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

The 110 Bo2023 regions are used as a dense, repeated-donor primate anatomical
reference, not as direct human anatomical identity. Broad transcriptional
gradients and network-linked cortical expression reported across primates
support a cautious cross-species organizational bridge (Krienen et al., 2016;
Burt et al., 2018), while AHBA mapping, orthology loss and species-specific
association-cortex differences are audited explicitly. Saleem nomenclature is
therefore a naming crosswalk rather than proof of one-to-one homology.

The 3′ analysis is a transcript-coordinate/detection proxy sensitivity test;
it is not molecule-level read generation and must not be described as a true
3′-end sequencing simulation. GO/KEGG and cell-type analyses annotate the
frozen panel and do not establish cell of origin, mechanism, causality, or
clinical predictive performance.

The Round-5 P0-3 formal differential-expression audit uses the immutable
Bo2023 feature-count matrix. Raw counts are summed within each observed
donor-by-Network cell (74 pseudobulk units from 9 macaques) before fitting a
DESeq2 likelihood-ratio test with full design `~ donor_id + network` and
reduced design `~ donor_id`. This avoids treating multiple tissue dissections
from one animal as independent biological replication. The audit tests the
frozen Network Top200; it does not retrain or replace the production panel.

The Round-5 P0-4 enrichment sensitivity audit compares the same locked Top200
query under three universes in one g:Profiler database version: the 21,668-gene
model space, the g:Profiler annotated domain, and 21,375 unique nonblank gene
symbols among all genes tested in the primary donor-by-Network pseudobulk
DESeq2 analysis (with no significance filter). Full returned GO:BP/KEGG term
tables, stability metrics, service metadata and SHA-256 values are retained in
`reproducibility/round5_analysis/p0_4_background_sensitivity/`.

## Immutable Snapshot Policy

The validation results archived in `reproducibility/v4_p*.csv` and summarized
in the manuscript correspond to frozen snapshots of these datasets. Users
re-running `reproduce_all.py` must obtain the same data versions from their
original repositories to reproduce the exact numerical results. Software
version v0.1.16 is the current immutable Zenodo release under version DOI
10.5281/zenodo.21974954; the separate full reproducibility archive is
10.5281/zenodo.21974991 and contains the materialized 164 MB Git LFS payload
because Zenodo's GitHub-generated source archive does not materialize Git LFS
objects. The previous v0.1.15 software/full-archive records remain immutable at
10.5281/zenodo.21970252 and 10.5281/zenodo.21970278. The immutable v0.1.14
historical software/full-archive records remain available at 10.5281/zenodo.21920261
and 10.5281/zenodo.21920697. The concept DOI resolves to the latest Zenodo
version, and the v0.1.12 scientific release remains available at
10.5281/zenodo.21911532.

## SHA-256 Integrity

The `data/models/canonical110_model_lock.json` manifest records SHA-256
checksums for all locked model artifacts. Production inference verifies these
hashes at load time and fails closed on any mismatch. See `MODEL_LOCK.md` for
the complete lock specification.

The frozen 200-gene RF row-level result has SHA-256
`d27b67c5ca7fa79b9186a79b24d94cf2624c4007c95bb53016043f554ceba3ac`.
The formal LOMO Network prediction-level source is
`reproducibility/p2_publication_completeness/formal_lomo_network_detail.csv`;
its route-family filter, SHA-256, integer TP/FP/FN accounting and derived
macro/weighted/micro F1 values are recorded in
`reproducibility/lomo_network_f1_provenance.json` and the generated
`reproducibility/LOMO_NETWORK_F1_PROVENANCE.md` report.
The independent enrichment manifest records the input SHA-256 values for the
locked panel, 21,668-gene projector background, and the Chiou2023 and
Siletti2023 source spreadsheets. Derived enrichment outputs, marker-family
intersections, protocol and manifest are public under
`reproducibility/independent_enrichment/`; publisher-hosted source workbooks
are not redistributed. The Saleem crosswalk is public under
`reproducibility/crosswalks/`, and the portable S18-S19 runner and input hashes
are public under `reproducibility/s18_s19/`.

The Round-5 audit outputs are under `reproducibility/round5_analysis/`. They
separate the Network-layer projection ablation, paired-reference OLS fit
quality, LOSO/LOMO ranking metrics, and the P0-3 donor-Network pseudobulk
DESeq2 marker audit, plus the P0-4 three-universe GO:BP/KEGG background
sensitivity audit. Source counts and the publisher workbook are not
redistributed; filenames and hashes are recorded in the P0-3 input manifest.
The OLS statistics are in-sample engineering-fit diagnostics and are not
evidence of DESeq2-VST equivalence or external calibration.
The P1 E4 mapping-impact output is
`reproducibility/p1_e4_ahba_network_mapping_impact.csv`, with its derivation
manifest in `p1_e4_ahba_network_mapping_impact_manifest.json`. CI branch
coverage is archived in `reproducibility/coverage_report.txt` and
`reproducibility/coverage.xml`.
