# P0-3 DESeq2 marker audit results

## Primary donor-Network pseudobulk analysis

Raw counts from the 819 Bo2023 tissue dissections were summed within each observed donor-by-Network cell, yielding 74 pseudobulk units from 9 macaques and 10 canonical Networks. The donor-blocked design was full rank (18/18 columns): full `~ donor_id + network` versus reduced `~ donor_id`. After the prespecified filter of at least 10 pseudobulk counts in at least 3 units, 23,331 of 28,415 genes were tested.

The global DESeq2 likelihood-ratio test identified 16,024/23,331 genes (68.68%) with an any-Network effect at BH-FDR < 0.05. All 200 locked Network-panel genes mapped to the raw-count/gene-annotation space, were tested, and passed BH-FDR < 0.05. The least significant locked gene still had adjusted P = 6.23e-7. The median descriptive highest-Network-versus-other-Network mean effect among the locked genes was 0.628 log2 units (range 0.175-7.992).

The historical Fisher-like predictive ranking had weak concordance with formal DESeq2 evidence across the locked panel (Spearman rho = 0.098 with both the LRT statistic and -log10 adjusted P). This is not a contradiction: the Fisher-like score was optimized to rank a fixed-size predictive panel on VSD-scale data, whereas the DESeq2 LRT tests whether any Network term improves a negative-binomial count model after donor blocking. The DESeq2 result supports differential expression of the panel but does not retroactively redefine its predictive ordering.

## Tissue-level sensitivity analysis

A donor-fixed-effect LRT retaining all 819 tissue dissections identified 18,698/19,735 tested genes (94.74%) and supported all 200 locked genes. Because residual tissue dissections from the same animal are not independent donor replicates, this analysis can overstate precision and is retained only as a sensitivity analysis. The manuscript-facing inference uses the more conservative donor-Network pseudobulk result.

## Limitations

- The donor-by-Network table is incomplete. The design is estimable, but not every donor contributes every Network.
- Hippocampal formation has only three pseudobulk donors (8 original tissue samples), so its Network-specific descriptive effects have limited donor-level precision.
- Pseudobulk units contain 1-63 tissue dissections. DESeq2 size factors address library-scale differences, but the breadth of anatomical sampling within a Network is not identical across donors.
- The global LRT provides formal evidence for an any-Network expression effect. The reported highest-Network effect is descriptive and does not constitute a second family of Network-specific hypothesis tests.
- This is an inferential audit of the frozen Network Top200, not a replacement panel, cross-species validation, cfRNA validation, or new localization-accuracy experiment.
