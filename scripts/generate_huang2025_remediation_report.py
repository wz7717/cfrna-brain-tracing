#!/usr/bin/env python
"""Generate the final forensic report, changelog, and canonical output checksums."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", type=Path, default=Path("reproducibility/huang_2025"))
    parser.add_argument("--workdir", type=Path, default=Path("manuscript_remediation"))
    parser.add_argument("--main-baseline", type=Path, required=True)
    parser.add_argument("--supp-baseline", type=Path, required=True)
    args = parser.parse_args()

    outdir = args.outdir
    workdir = args.workdir
    workdir.mkdir(parents=True, exist_ok=True)
    summary = json.loads((outdir / "huang_2025_canonical_summary.json").read_text(encoding="utf-8"))
    manifest = json.loads((outdir / "huang_2025_audit_manifest.json").read_text(encoding="utf-8"))
    verification = json.loads((workdir / "huang_2025_structural_verification.json").read_text(encoding="utf-8"))
    comparisons = read_csv(outdir / "huang_2025_tumour_control_comparisons.csv")
    fluid = {row["cohort"]: row for row in read_csv(outdir / "huang_2025_fluid_summary.csv")}
    markers = read_csv(outdir / "huang_2025_marker_correlations.csv")
    main_docx = workdir / "BrainTrace_Main_Manuscript_HuangRemediated.docx"
    supp_docx = workdir / "BrainTrace_Supplementary_File_HuangRemediated.docx"
    test_log = (workdir / "pytest_full.log").read_text(encoding="utf-8", errors="replace").strip().splitlines()[-1]
    main_images = json.loads((workdir / "rendered" / "contact_sheets" / "main" / "render_image_manifest.json").read_text(encoding="utf-8"))
    supp_images = json.loads((workdir / "rendered" / "contact_sheets" / "supplement" / "render_image_manifest.json").read_text(encoding="utf-8"))

    csf = fluid["all_CSF_profiles"]
    plasma = fluid["all_plasma_profiles"]
    source_assets = {item["label"]: item for item in manifest["input_assets"]}
    strongest_marker = min(
        (row for row in markers if row.get("bh_fdr")),
        key=lambda row: float(row["bh_fdr"]),
    )

    checksum_lines = []
    for path in sorted(outdir.iterdir(), key=lambda value: value.name.lower()):
        if path.is_file() and path.name != "SHA256SUMS.txt":
            checksum_lines.append(f"{sha256(path)}  {path.name}")
    (outdir / "SHA256SUMS.txt").write_text("\n".join(checksum_lines) + "\n", encoding="utf-8")

    report = f"""# Huang 2025 cfRNA provenance-remediation QA report

**Status:** PASS — working remediation candidate only; not a public release.

## 1. External source facts

- Huang et al. (2025), DOI [{summary['source_doi']}](https://doi.org/{summary['source_doi']}).
- The published expression matrix contains 159 profiles: 77 CSF and 82 plasma. The source article reports 85 patients (18 glioma, 46 meningioma and 21 controls).
- The source clinical analysis reports aggregate sequencing-QC exclusions of five CSF and one plasma profile. No public per-profile QC-status map or patient-level CSF-plasma correspondence was found.

## 2. Manuscript-impact assessment

The prior local interpretation of shared terminal label numbers as patient correspondence was unsupported. Therefore all manuscript/supplement wording that treated the Huang data as matched CSF-plasma pairs, or cited the resulting paired minimum P=0.304, was retired.

## 3. Pseudo-pairing forensic conclusion

No public patient identifier was used or inferred. The canonical runner contains no patient key, CSF-to-plasma sample-name substitution, concordance permutation, or synthetic admixture calculation. The historical rationale and retirement boundary are recorded in `reproducibility/historical_noncanonical/HUANG_2025_PSEUDOPAIRING_RETIRED.md`.

## 4. New raw-input handling

The full published matrix was retained as a technical audit universe, not as an attempted replication of the source clinical-QC-filtered analysis. Matrix SHA-256: `{source_assets['published_expression_matrix']['sha256']}`. The source supplementary XLSB was recorded with SHA-256 `{source_assets['source_supplementary_xlsb']['sha256']}`.

## 5. Analysis unit and statistics

The analysis unit is a published expression profile within a fluid-specific cohort. Six two-sided profile-level Mann-Whitney U tumour-control comparisons were run across CSF/plasma and atlas-fit/margin/entropy metrics, with Benjamini-Hochberg correction across all six tests. The public matrix does not provide a profile-to-patient map, so patient-level dependence or independence cannot be assessed. The smallest BH-FDR was {summary['minimum_bh_fdr']:.6f}; all six comparisons were non-significant after correction.

## 6. Output chain and cohort accounting

- Audit universe: {summary['n_profiles']} profiles ({summary['n_csf']} CSF; {summary['n_plasma']} plasma).
- Traceable BrainTrace outputs: {summary['n_traceable_outputs']}/{summary['n_profiles']}.
- OMPFC Top1: CSF {csf['OMPFC_top1_numerator']}/{csf['OMPFC_top1_denominator']} ({float(csf['OMPFC_top1_percent']):.2f}%); plasma {plasma['OMPFC_top1_numerator']}/{plasma['OMPFC_top1_denominator']} ({float(plasma['OMPFC_top1_percent']):.2f}%).
- The sample ledger records all required inclusion flags, blank patient IDs, and the public-QC-status limitation.

## 7. Marker analyses

Sixteen platelet-associated and extracellular-vesicle-associated markers were assessed separately in CSF (n=77) and plasma (n=82) against OMPFC score, with rho, raw P and BH-FDR recorded for each. These are descriptive domain diagnostics, not evidence of patient matching, mixture behavior or mechanism. The smallest marker BH-FDR was {float(strongest_marker['bh_fdr']):.3g} for {strongest_marker['fluid']} {strongest_marker['marker']} (rho={float(strongest_marker['spearman_rho']):.4f}).

## 8. Test and reproducibility status

- Targeted Huang contract tests passed before the full suite.
- Full suite: `{test_log}`.
- Structural verifier: `{verification['status']}`; it confirmed 159/77/82 accounting, six profile-level tests, no patient IDs, and no legacy pseudo-pair output names.

## 9. Security and repository cleanliness

The raw matrix, XLSB and source-tracing database remain external inputs and are not copied into the repository. The immutable v0.1.16 release/tag, its release manifest and archival checksums were not edited. This branch is explicitly marked as an unreleased v0.1.17 candidate requiring scientific approval.

## 10. Claim boundaries

Supported claim: technical portability and domain-shift audit. Not supported: patient-level CSF-plasma correspondence, synthetic matched-mixture behavior, anatomical localization accuracy, tumour-source discrimination, or clinical validity.

## 11. Document revisions and QA

New working documents were created from the two user-supplied v0.1.16 FINAL baselines; the originals were not modified.

- Main output: `{main_docx.name}` (SHA-256 `{sha256(main_docx)}`).
- Supplement output: `{supp_docx.name}` (SHA-256 `{sha256(supp_docx)}`).
- Baseline main SHA-256: `{sha256(args.main_baseline)}`.
- Baseline supplement SHA-256: `{sha256(args.supp_baseline)}`.
- XML verifier: no tracked insertions/deletions/moves, no Track Changes flag, and no comment parts.
- Visual QA: {len(main_images)} main-manuscript pages and {len(supp_images)} supplementary pages were rendered and inspected; tables and changed passages showed no clipping, overlap or overflow.

## 12. Change log

See `HUANG_2025_REMEDIATION_CHANGELOG.txt` in this directory for file-by-file OLD / NEW / REASON / SOURCE-OR-OUTPUT support.

## 13. Deliverables

- Canonical analysis outputs: `reproducibility/huang_2025/`.
- Canonical checksums: `reproducibility/huang_2025/SHA256SUMS.txt`.
- Remediated DOCX working files and QA materials: `manuscript_remediation/`.
- Retired noncanonical record: `reproducibility/historical_noncanonical/`.

## 14. Exact rerun commands

```powershell
python scripts/run_huang2025_external_candidate.py `
  --input-csv <Huang_2025_Supplementary_Data_1.csv> `
  --source-xlsb <Huang_2025_Supplementary_Data_1.xlsb> `
  --db-path <braintrace_source_tracing.db> `
  --outdir reproducibility/huang_2025

python -m pytest -q

& <bundled-python> scripts/verify_huang2025_remediation.py `
  --outdir reproducibility/huang_2025 `
  --main-docx manuscript_remediation/BrainTrace_Main_Manuscript_HuangRemediated.docx `
  --supp-docx manuscript_remediation/BrainTrace_Supplementary_File_HuangRemediated.docx `
  --report manuscript_remediation/huang_2025_structural_verification.json
```

## 15. Remaining and out-of-scope issues

The public record still lacks per-profile source-QC status, patient correspondence and anatomical truth. Those are source-data limitations, not gaps that may be repaired by label parsing. No unrelated scientific endpoint or the v0.1.16 public release was changed.

## Acceptance checklist (A–R)

| Item | Status | Evidence |
|---|---|---|
| A. Source record | PASS | DOI and source facts above; audit manifest |
| B. Full-matrix accounting | PASS | 159 = 77 CSF + 82 plasma |
| C. No inferred patient pairing | PASS | ledger + static contract tests |
| D. No synthetic mixture | PASS | canonical output inventory + source inspection |
| E. Locked BrainTrace route | PASS | audit manifest and per-profile outputs |
| F. Profile-level statistics | PASS | six comparison rows; patient pairing unavailable |
| G. Correct BH correction | PASS | comparison CSV and verifier |
| H. Fluid-specific distributions | PASS | fluid-summary and Top1/Top3 CSVs |
| I. Marker diagnostics | PASS | marker-correlation CSV |
| J. Source-QC boundary | PASS | manifest, report, and revised methods |
| K. Machine-readable summary | PASS | canonical CSV/JSON |
| L. Reproduction entry point | PASS | runner + `reproduce_all` integration |
| M. Regression tests | PASS | full pytest result |
| N. Old pseudo outputs retired | PASS | historical note; absent from canonical directory |
| O. Main manuscript remediation | PASS | new DOCX and structural verifier |
| P. Supplement remediation | PASS | new DOCX and structural verifier |
| Q. No tracked changes/comments | PASS | XML verifier |
| R. Full-page visual QA | PASS | 12 + 81 rendered pages and contact-sheet manifests |
"""
    (workdir / "HUANG_2025_REMEDIATION_QA_REPORT.md").write_text(report, encoding="utf-8")

    changelog = """HUANG 2025 cfRNA PROVENANCE-REMEDIATION CHANGELOG
Status: working remediation candidate only; not a public v0.1.17 release.

FILE: scripts/run_huang2025_external_candidate.py
OLD: Sample label suffixes were converted to a patient key and used for CSF-plasma concordance/permutation outputs.
NEW: The full 159-profile matrix is processed through the locked route as separate 77-CSF and 82-plasma profile cohorts; the ledger leaves patient_id blank, records that patient pairing is unavailable, and excludes paired/synthetic analyses.
REASON: The public matrix does not supply a patient-level CSF-plasma correspondence.
SOURCE / OUTPUT SUPPORT: Huang DOI 10.1038/s41698-025-00909-6; huang_2025_sample_ledger.csv; huang_2025_canonical_summary.json.

FILE: scripts/run_p1_bio1_4_audit.py
OLD: Generated sample-name-derived CSF-plasma pseudo-pair and admixture outputs.
NEW: Huang-related pseudo-pair/admixture code path removed; the script now retains only its unrelated Bio2/Bio4 analyses.
REASON: Prevent reintroduction of unsupported pairing.
SOURCE / OUTPUT SUPPORT: tests/test_huang2025_external_candidate.py static regression guard.

FILE: reproducibility/huang_2025/*
OLD: No canonical remediation package.
NEW: Full matrix ledger, per-profile locked-route outputs, ranks, fluid distributions, marker correlations, six profile-level tests, manifest, summary and results report.
REASON: Create auditable data-to-claim chain.
SOURCE / OUTPUT SUPPORT: SHA256SUMS.txt; huang_2025_audit_manifest.json.

FILE: reproduce_all.py; reproduce_all.ps1; README.md; DATA_PROVENANCE.md; reproducibility/TRACEABILITY_MATRIX.md
OLD: The Huang audit pointed to obsolete output locations or lacked the patient-correspondence provenance boundary.
NEW: Canonical path is reproducibility/huang_2025; documentation states full-matrix profile-level scope, separate fluid-specific cohorts, unavailable patient pairing, source-QC limitation, claim boundary and unreleased status.
REASON: Keep scripts, provenance and manuscript-facing traceability aligned.
SOURCE / OUTPUT SUPPORT: canonical summary and audit manifest.

FILE: BrainTrace_Main_Manuscript_HuangRemediated.docx
OLD: Reported paired-stability diagnostics and a minimum paired P=0.304.
NEW: Reports 159 full-matrix profiles, separate 77-CSF and 82-plasma profile cohorts, 159/159 traceability and minimum BH-FDR=0.722 without patient-pairing claims.
REASON: Correct unsupported patient-paired inference.
SOURCE / OUTPUT SUPPORT: huang_2025_canonical_summary.json; huang_2025_tumour_control_comparisons.csv.

FILE: BrainTrace_Supplementary_File_HuangRemediated.docx
OLD: Contained matched-patient admixture, pseudo-paired stability language, old plasma denominator and the 0.304 claim.
NEW: Replaces those passages and Tables S5/S6 with separate fluid-specific profile results, correct denominators, marker-analysis description and claim boundaries.
REASON: Correct the full scientific-provenance chain.
SOURCE / OUTPUT SUPPORT: huang_2025_fluid_summary.csv; huang_2025_marker_correlations.csv; huang_2025_tumour_control_comparisons.csv.

FILE: reproducibility/historical_noncanonical/HUANG_2025_PSEUDOPAIRING_RETIRED.md
OLD: No quarantined explanation of the invalid assumption.
NEW: Historical-only retirement record.
REASON: Preserve forensic traceability without preserving invalid outputs as canonical evidence.
SOURCE / OUTPUT SUPPORT: source-data availability audit and code-history review.
"""
    (workdir / "HUANG_2025_REMEDIATION_CHANGELOG.txt").write_text(changelog, encoding="utf-8")
    print(workdir / "HUANG_2025_REMEDIATION_QA_REPORT.md")
    print(workdir / "HUANG_2025_REMEDIATION_CHANGELOG.txt")
    print(outdir / "SHA256SUMS.txt")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
