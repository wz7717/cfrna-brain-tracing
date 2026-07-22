# MB single-cell / developmental-origin signature sources

Generated: 2026-07-09

This folder contains downloaded supplementary tables and extracted gene-set candidates for medulloblastoma cell-state and developmental-origin analyses. The immediate target is downstream scoring of GSE189919 CSF RNA-seq samples and correlation against BrainTrace Network scores.

## Extracted output

- `mb_literature_signature_gene_sets.csv`
  Unified long-format gene-set table with columns:
  `source`, `table`, `signature`, `gene_symbol`, `rank`, `score`, `note`.
- `mb_literature_signature_gene_sets_summary.csv`
  Per-source, per-signature gene counts.

## Source 1: Hovestadt et al., Nature 2019

Article: "Resolving medulloblastoma cellular architecture by single-cell genomics"

Local files:

- `Hovestadt2019_Nature_supp_table2_transcriptional_programs.xlsx`
- `Hovestadt2019_Nature_supp_table3_neuronal_like_comparison.xlsx`

Useful extracted signatures:

- `WNT-A`, `WNT-B`, `WNT-C`, `WNT-D`
- `SHH-A`, `SHH-B`, `SHH-C`
- `Group3/4-A`, `Group3/4-B`, `Group3/4-C`
- `Neuron-like_WNT`, `Neuron-like_SHH`, `Neuron-like_Group 3/4`, `Neuron-like_shared`

Rationale:

This is the most directly useful source for MB cell-state programs. The paper states that WNT, SHH and Group 3 tumors contain undifferentiated and differentiated neuronal-like malignant populations; Group 4 tumors are differentiated neuronal-like; Group 3/4 spans primitive progenitor-like to mature neuronal-like states.

Recommended use:

- For broad MB cell-state scoring, use Hovestadt meta-programs first.
- For neuronal-like scoring, use Supplementary Table 3 gene sets.
- For compact scoring, consider top 30 genes per Hovestadt meta-program, matching the original paper's cell-scoring description.

## Source 2: Hendrikse et al., Nature 2022

Article: "Failure of human rhombic lip differentiation underlies medulloblastoma formation"

Local files:

- `Hendrikse2022_Nature_supp_table4_G3_G4_subtype_genes.xlsx`
- `Hendrikse2022_Nature_supp_table7_glutamatergic_cell_markers.xlsx`

Useful extracted signatures:

- Group 3/4 subtype DE gene sets:
  - `Group3_alpha`, `Group3_beta`, `Group3_gamma`
  - `Group4_alpha`, `Group4_beta`, `Group4_gamma`
- Human glutamatergic developmental markers:
  - `RL-VZ`, `RL-SVZ`
  - `Early_UBC`, `Late_UBC`
  - `GCP`, `Early_GN`, `Late_GN`

Rationale:

This is the strongest source for Group 3/4 subtype genes and human cerebellar glutamatergic lineage markers. It can support testing whether GSE189919 Network scores correlate with RL/UBC/GCP/GN-like programs.

Recommended use:

- Use Table 7 for compact developmental lineage signatures.
- Use Table 4 subtype DEGs carefully: the gene sets are large, so for sample scoring prefer top-ranked genes by log2 fold change or use rank/GSEA-style methods.

## Source 3: Smith et al., Nature 2022

Article: "Unified rhombic lip origins of group 3 and group 4 medulloblastoma"

Local file:

- `Smith2022_Nature_supp_table5_human_cerebellar_signatures.xlsx`

Useful extracted signatures:

- `RL_svz`, `RL_vz`, `EGL`, `PCL`
- `Photoreceptor_gene_set:Gene set: Descartes Photoreceptor cells`
- `Photoreceptor_gene_set:Gene set: RLsvz-Photoreceptor cells`
- `UBC_gene_set:Gene set: Descartes Unipolar Brush Cells`
- `UBC_gene_set:Gene set: RLsvz-Unipolar brush cells`

Rationale:

This is the strongest source for rhombic-lip-derived UBC and photoreceptor-like programs in Group 3/4 MB. It is useful for testing whether MB CSF Network scores align with rhombic lip, UBC, or photoreceptor-like signatures.

Recommended use:

- Use UBC and photoreceptor gene sets as targeted hypothesis tests.
- Use RL_svz/RL_vz/EGL/PCL as developmental compartment scores.

## Suggested analysis plan for GSE189919

1. Load GSE189919 raw counts and compute logCPM.
2. Filter signature genes to those present in GSE189919.
3. Score each sample using:
   - mean z-score or median z-score for compact gene sets;
   - rank-based AUCell/ssGSEA-like scoring for large DEG lists.
4. Join sample-level signature scores with:
   - `results/gse189919_latest_main_route_20260708/gse189919_latest_main_route_network_rankings.csv`
   - `results/gse189919_latest_main_route_20260708/gse189919_latest_main_route_sample_summary.csv`
5. Compute Spearman correlations:
   - signature score vs each Network score;
   - signature score vs selected Network ranks or Top3 membership.
6. Report FDR-corrected correlation tables and a heatmap.

## Interpretation boundary

These signatures come from tumor tissue or developmental atlases, whereas GSE189919 is CSF RNA-seq. A positive association should be interpreted as transcriptomic-program similarity in biofluid RNA, not direct tumor-cell fraction, anatomical origin, or localization accuracy.
