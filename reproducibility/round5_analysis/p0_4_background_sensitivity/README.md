# Round-5 P0-4 background-sensitivity audit

This package tests whether GO:BP and KEGG conclusions for the immutable
Network Top200 panel depend materially on the enrichment universe.

All three queries were executed together against g:Profiler database version
`e114_eg62_p19_27110d83` with the same 200-gene query, FDR correction and
GO:BP/KEGG sources:

1. the prespecified 21,668-gene frozen model space;
2. the g:Profiler annotated domain; and
3. 21,375 unique nonblank gene symbols among all 23,331 genes tested in the
   primary donor-by-Network pseudobulk DESeq2 audit, without a significance
   filter.

The expression background produced 438 significant GO:BP terms and 10 KEGG
pathways, compared with 446/11 for the model background and 410/8 for the
annotated domain. Relative to the model background, it retained 438/446
GO:BP terms (98.2%) and 10/11 KEGG terms (90.9%). Spearman correlations of
`-log10(FDR)` over all common tested terms were 0.99995 and 0.99988,
respectively. Across all three universes, 406 GO:BP and 8 KEGG terms were
significant. Thus the broad annotation pattern is stable, while individual
near-threshold terms remain background-sensitive.

`outputs/` contains all returned terms, pairwise and three-way stability
metrics, and the prespecified representative-theme audit. `p0_4_manifest.json`
records inputs, hashes, background construction and service metadata;
`SHA256SUMS.txt` protects all derived files. Re-running the public script may
change results if the upstream annotation database changes.

These are annotation-bias analyses, not evidence of mechanism, cell of origin,
causality or predictive validity.
