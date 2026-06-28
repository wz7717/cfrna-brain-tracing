# cfRNA-BrainTrace submission-ready final v5 changelog

Date: 2026-06-28

## Figure 1 replacement

- Figure 1 was replaced with the new two-panel workflow + validation version.
- Panel A shows the cfRNA-BrainTrace workflow.
- Panel B shows the validation summary.
- Final figure files were exported/copied to:
  - `/mnt/data/Figure1_cfRNA_BrainTrace_final.svg`
  - `/mnt/data/Figure1_cfRNA_BrainTrace_final.png`
  - `/mnt/data/Figure1_cfRNA_BrainTrace_final.pdf`
- The DOCX files embed the 600 dpi PNG version and retain the SVG/PDF as vector source files.

## Manuscript text

- `Fig. 1A` was added in Section 2 System and methods after the operational workflow description.
- `Fig. 1B` was added in Section 4 Validation after the end-to-end validation summary.
- The old validation-only Figure 1 caption was removed.
- The Figure 1 caption was replaced with the workflow + validation evidence caption requested for v5.
- Figure 1 alt text was updated to the requested two-panel description.

## Validation-number consistency

- The Figure 1 Panel B values were checked against the manuscript text and alt text:
  - Internal LOSO: Network Top3 92.19%, resolution-group Top3 72.36%, exact-region Top3 45.33%.
  - Internal LOMO: Network Top3 91.21%, resolution-group Top3 69.09%, exact-region Top3 42.36%.
  - AHBA mapped-label: Network Top3 94.42%, resolution-group Top3 67.03%, exact-region Top3 42.86%.
  - TCGA/BraTS coarse consistency: Network Top3 40.00%, broad-anatomy Top3 64.62%.
- TCGA/BraTS broad anatomy was kept as a separate category, labelled in the figure legend as `Broad anatomy (TCGA/BraTS only)`.
- No validation numbers were changed in the manuscript body.

## Claim boundary

- Biofluid boundary language was retained in the figure note, caption and alt text.
- The manuscript still treats biofluid datasets without patient-level anatomical truth as projection-feasibility or transfer-stress analyses only.
- TCGA/BraTS remains described only as coarse anatomical consistency, not macaque exact-region validation.
- Network remains the primary endpoint, resolution group remains secondary, and exact region remains exploratory.
- No new claims of deterministic brain-region prediction, clinical cfRNA localization, exact-region localization endpoint, biofluid localization accuracy or GSE189919 localization accuracy were introduced.

## QA

- No old Figure 1 or old validation-only caption remains in the v5 manuscript outputs.
- No duplicate Figure 1 is present in either DOCX.
- No standalone `Figure legends` heading remains.
- Word core metadata titles were cleaned in both v5 DOCX files to remove inherited `v4 review format` / `v4 clean` labels.
- The review-format DOCX retains continuous line numbering and rendered to 8 PNG pages with LibreOffice/soffice + Poppler; all pages were visually inspected.
- The clean DOCX rendered to 5 PNG pages with LibreOffice/soffice + Poppler; all pages were visually inspected.
- Panel A, Panel B, A/B labels, bar labels, legend and bottom biofluid boundary note are visible in the rendered DOCX outputs.
- Poppler emitted non-blocking `nameToUnicode` warnings caused by the local Chinese user-profile path; PNG rendering completed successfully.

## Output files

- Review-format Word: `/mnt/data/Bioinformatics_Application_Note_cfRNA_BrainTrace_submission_ready_final_v5_review_format.docx`
- Clean Word: `/mnt/data/Bioinformatics_Application_Note_cfRNA_BrainTrace_submission_ready_final_v5_clean.docx`
- Markdown: `/mnt/data/Bioinformatics_Application_Note_cfRNA_BrainTrace_submission_ready_final_v5.md`
- Changelog: `/mnt/data/Bioinformatics_Application_Note_cfRNA_BrainTrace_submission_ready_final_v5_changelog.md`

## Final status

Ready after author confirmation.
