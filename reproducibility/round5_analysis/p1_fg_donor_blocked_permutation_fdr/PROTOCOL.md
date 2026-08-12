# P1 E3 donor-blocked F_g permutation FDR

This run is the formal donor-aware sensitivity for the descriptive F_g
Network-ranking score. It does not replace the donor-by-Network DESeq2 LRT.

- Input: 819 Bo2023 tissue samples, 28415 raw genes.
- Blocking: repeated dissections were summed within each observed donor x
  region block, producing 459 blocks from
  9 donors and 110 regions.
- Transformation: logCPM on the block-level raw counts, followed by gene-wise
  centering within donor.
- Statistic: F_g = between-Network variance / within-Network variance, with
  df_between=9 and df_within=449.
- Null: Network labels were permuted within donor, preserving each donor's
  block count and observed Network composition.
- Multiplicity: 5000 permutations, seed 20260809, empirical p-values
  with the add-one correction, then BH across the 21,668-gene
  frozen model universe.
- Production panel: the frozen Top200 was not reselected or reordered.

The earlier `p1_fg_permutation_fdr` directory is retained as a historical
110-region, non-donor-aware exploratory sensitivity and must not be described
as donor-aware.
