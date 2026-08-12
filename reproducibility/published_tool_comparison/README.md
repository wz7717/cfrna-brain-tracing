# Published-tool comparison (P0-2)

This directory separates an actually executed TissueEnrich analysis from a
historically prepared but explicitly discontinued CIBERSORTx branch. The two tools do
not estimate BrainTrace's 110-region hierarchy, so their outputs are treated as
complementary tissue/cell-composition evidence rather than Top1/Top3 localization
competitors.

## TissueEnrich: completed

TissueEnrich 1.32.0 was actually run under R 4.6.1/Bioconductor 3.23. The query
was the frozen 200-gene Network panel and the custom background was the frozen
21,668-gene projector space. HPA and GTEx were evaluated independently using
all TissueEnrich tissue-specific gene categories and BH correction within each
dataset.

The strongest result in HPA was Cerebral Cortex (76 query genes, 4.94-fold,
BH-adjusted P = 4.38e-33). The strongest result in GTEx was Brain (33 query
genes, 4.14-fold, BH-adjusted P = 7.89e-11). These findings establish that the
locked panel is enriched for brain-associated genes in an independent published
tool; they do not establish brain-region localization accuracy or cell of origin.

Run from the workspace root after installing the recorded R environment:

```text
external_tools/R-4.6.1/bin/Rscript.exe \
  github_main_sync/reproducibility/published_tool_comparison/run_tissueenrich.R \
  github_main_sync \
  github_main_sync/reproducibility/published_tool_comparison/outputs
```

## CIBERSORTx: input prepared, execution discontinued

A local-only, upload-ready TPM mixture was prepared from the same 65 frozen
TCGA/BraTS cases used by BrainTrace. The mixture, its manifest and the licensed
LM22 input are deliberately not redistributed in the public repository.

CIBERSORTx was not run, and the user has explicitly terminated this optimization
branch. The official web service requires user login and acceptance of Stanford's
non-commercial terms, and the official fractions container requires CIBERSORTx
credentials. The former analysis plan is retained here only as provenance; it is
not an active pending task. No CIBERSORTx result, immune-composition estimate or
domain-shift comparison is available. LM22 must not be described as a brain-cell
or anatomical localization estimate.

No binary marker matrix has been fabricated from Chiou/Siletti marker lists:
CIBERSORTx requires quantitative cell-type expression profiles, and ranked
marker names alone are not a valid signature expression matrix.

## Status rule

P0-2 CIBERSORTx branch: discontinued by user decision. TissueEnrich is complete
and reproducible; CIBERSORTx was not executed. Do not resume this branch unless a
future user request explicitly reopens it.
