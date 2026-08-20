# Public release content policy

BrainTrace public releases are limited to scientific software, frozen model
artifacts, provenance records, reproducibility workflows, validation code and
read-only scientific QA material.

The following material is outside every public release surface, including the
Git repository, source archives, full-reproducibility archives and container
build contexts:

- Main manuscript, supplementary-material and cover-letter documents;
- utilities that write, revise, format or lay out submission documents;
- utilities that generate proposed wording, change lists, response text or
  submission-compression checklists;
- rendered manuscript-figure revision bundles derived from submission-document
  media;
- local submission assembly, rendering and office-review artifacts.

Scientific computation remains in scope when it generates data, statistics,
figures, provenance, machine-readable evidence or reproducibility reports.
Read-only validators may inspect scientific artifacts without modifying a
submission document.

Before a public change is accepted, run:

```bash
python scripts/audit_public_release_content.py
python scripts/generate_package_integrity.py --check
```

The first command fails closed on prohibited public-tree content. The second
checks that `PACKAGE_MANIFEST.csv` and `SHA256SUMS.txt` exactly describe the
current intended public tree. Release construction repeats the content audit,
and CI runs both checks.
