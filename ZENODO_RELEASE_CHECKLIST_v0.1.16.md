# Zenodo release checklist — BrainTrace v0.1.16

Do not insert or infer a v0.1.16 DOI before Zenodo creates the immutable
version record.

After the GitHub v0.1.16 release is created:

- [ ] Create a new Zenodo version from the v0.1.16 GitHub release.
- [ ] Confirm title, creators, MIT license, repository URL, publication date
      and version `v0.1.16`.
- [ ] Confirm the archived source contains `examples/`, Help & Tutorial,
      `VERSION_STRING_AUDIT_v0.1.16.txt` and
      `SUBMISSION_ONLINE_RESOURCE_QA_v0.1.16.md`.
- [ ] Upload/materialize the full reproducibility archive, including the
      approximately 164 MB Git LFS payload, if the GitHub-generated Zenodo
      source archive does not materialize LFS objects.
- [ ] Verify the full archive against its release checksum manifest.
- [ ] Record the real Zenodo v0.1.16 version DOI; never reuse the v0.1.15 DOI
      as though it identified v0.1.16.
- [ ] Add `RELEASE_INTEGRITY_v0.1.16.json` with the real Git tag, GitHub
      release URL, version DOI, full-archive DOI/asset and SHA-256 values.
- [ ] Synchronize the real version DOI and archive wording in `README.md` and
      `DATA_PROVENANCE.md`.
- [ ] Synchronize the DOI in the manuscript and Supplementary Material where
      the current software release is cited.
- [ ] Confirm the Zenodo concept DOI still resolves to the latest version.
- [ ] Re-run `braintrace validate`, `python examples/verify_examples.py`, the
      full test suite and package checksum verification after metadata edits.

The immutable v0.1.15 records remain historical and must not be edited.
