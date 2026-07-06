# Reference audit report v6

Final status: Ready

Scope audited:
- Main manuscript DOCX and Markdown v6 references-checked outputs.
- Current submission package manuscript Markdown/DOCX copies.
- Supplementary DOCX, PDF, Markdown source and embedded supplementary text.
- Package text assets with extensions md, txt, csv, tsv, json, yml and yaml.
- Figure 1 caption and alt text in the main manuscript.

## Main manuscript references retained

1. Bakas,S. et al. (2017) Advancing The Cancer Genome Atlas glioma MRI collections with expert segmentation labels and radiomic features. Sci. Data, 4, 170117.
Citation location: main Markdown line 48.
Citation sentence: In TCGA/BraTS glioma tissue RNA-seq with MRI-derived labels from the TCGA glioma MRI collections (Bakas et al., 2017), results support only coarse anatomical consistency because MRI truth labels are human atlas labels rather than Bo2023 macaque exact-region identifiers.

2. Bo,T. et al. (2023) Brain-wide and cell-specific transcriptomic insights into MRI-derived cortical morphology in macaque monkeys. Nat. Commun., 14, 1499.
Citation location: main Markdown line 28.
Citation sentence: Network denotes the 10-class macaque functional-anatomical source space used by the Bo2023 macaque transcriptomic reference (Bo et al., 2023): Cingulate gyrus; Frontal agranular motor areas; Hippocampal formation; Lateral Prefrontal Cortex; Occipital/Temporal; Operculum/Insula; Orbitomedial Prefrontal Cortex; Parietal and Parieto-occipital region; Subcortical; and Temporal.

3. Hawrylycz,M.J. et al. (2012) An anatomically comprehensive atlas of the adult human brain transcriptome. Nature, 489, 391-399.
Citation location: main Markdown line 48.
Citation sentence: In Allen Human Brain Atlas (AHBA) mapped-label external validation (Hawrylycz et al., 2012), the locked production route achieved Network Top1/Top3 accuracy of 74.68% (95% CI 68.73-79.83%)/94.42% (90.69-96.71%), resolution-group Top1/Top3 accuracy of 36.26% (27.13-46.51%)/67.03% (56.86-75.83%) and exact Top1/Top3 accuracy of 24.18% (16.54-33.90%)/42.86% (33.18-53.11%).

4. Vorperian,S.K. et al. (2022) Cell types of origin of the cell-free transcriptome. Nat. Biotechnol., 40, 855-861.
Citation location: main Markdown line 24.
Citation sentence: RNA expression profiles, including cell-free RNA profiles, can retain tissue- and cell-of-origin information (Vorperian et al., 2022), but within-brain tracing is limited by regional similarity, atlas granularity and domain shift.

## Main manuscript references removed

Removed from the main Reference list because they were not cited in the main text and were not necessary for the short Application Note core argument:
- Avila Cobos,F. et al. incorrect 2021 Brief. Bioinform. bbab265 mixed entry.
- GTEx Consortium. (2020).
- Jain,A. and Tuteja,G. (2019).
- Kang,H.J. et al. (2011).
- Miller,J.A. et al. (2014).
- Newman,A.M. et al. (2019).
- Papatheodorou,I. et al. (2020).
- Tabula Sapiens Consortium. (2022).
- Uhlen,M. et al. (2015).
- Wilson,E.B. (1927), moved out of the main manuscript because the main text reports confidence intervals but does not specify the Wilson method.

## Erroneous or mismatched references

Avila Cobos mismatch found: Yes.
Action taken: Deleted from the main manuscript. The manuscript and supplement do not substantively discuss transcriptomic deconvolution benchmarking, so the corrected Avila Cobos 2020 Nat. Commun. reference was not added.

## Supplementary references retained

1. Wilson,E.B. (1927) Probable inference, the law of succession, and statistical inference. J. Am. Stat. Assoc., 22, 209-212.
Citation location: Supplementary Markdown line 71.
Citation sentence: Wilson confidence intervals (Wilson, 1927) and statistical tests were generated for the formal P0 evidence package.
Reference list location: Supplementary Markdown lines 276-278 under Supplementary references.

No other Supplementary reference list entries were retained. Dataset names such as Bo2023, AHBA and TCGA/BraTS remain dataset/source-space names in the supplementary narrative unless they are formal author-year citations.

## Citation-reference consistency checks

Main manuscript author-year citations detected:
- Bakas et al., 2017: matched to retained main Reference list.
- Bo et al., 2023: matched to retained main Reference list.
- Hawrylycz et al., 2012: matched to retained main Reference list.
- Vorperian et al., 2022: matched to retained main Reference list.

Main Reference entries without main-text citation: None detected.
Main-text author-year citations without Reference entry: None detected.
Supplementary author-year citations without Supplementary reference entry: None detected.
Supplementary reference entries without Supplementary citation: None detected.

## Markdown and LaTeX residue check

Checked patterns included visible Markdown/LaTeX residue markers: starred et al., backticks, escaped parentheses, and escaped lambda notation.
Main v6 Markdown: no residue detected.
Main v6 DOCX extracted text: no residue detected.
Supplementary v6 DOCX extracted text: no residue detected.
Supplementary package Markdown: no residue detected.
Package text assets scanned: no Avila Cobos, bbab265, starred et al., backticks, escaped parentheses or escaped lambda patterns detected.

## Validation and content boundary

No validation numbers, software route descriptions, Figure 1 artwork, Figure 1 caption values, claim boundary, data availability links, AI-assisted editing disclosure, author order or corresponding-author order were changed during this reference cleanup.

## Render and package QA

Main v6 DOCX was exported to PDF for QA and rendered to 10 page PNGs.
Supplementary v6 DOCX was exported to the requested PDF and rendered to 35 page PNGs.
Visual contact-sheet inspection found complete page sequences and no obvious text overflow or broken tables caused by the reference edits.

Final state: Ready.
