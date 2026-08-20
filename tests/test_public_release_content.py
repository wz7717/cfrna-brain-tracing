from __future__ import annotations

from pathlib import Path

import pytest

from scripts.audit_public_release_content import ROOT, audit, classify
from scripts.generate_package_integrity import verify


def categories(relative: str, source: str = "") -> set[str]:
    return {finding["category"] for finding in classify(relative, source.encode("utf-8"))}


def test_current_public_tree_passes_content_audit() -> None:
    report = audit(ROOT)
    assert report["status"] == "PASS", report["findings"]


def test_submission_authoring_filename_is_rejected() -> None:
    relative = "scripts/" + "update_main_" + "manuscript.py"
    assert "PROHIBITED_MANUSCRIPT_AUTHORING_SCRIPT" in categories(relative, "print('draft')\n")


def test_main_supplement_docx_writer_is_rejected() -> None:
    source = "\n".join(
        [
            "from " + "docx import Document",
            "MAIN_" + "SOURCE = 'draft.docx'",
            "document = Document()",
            "document." + "save('draft.docx')",
        ]
    )
    assert "PROHIBITED_MAIN_SUPPLEMENT_DOCX_WRITER" in categories(
        "scripts/scientific_transform.py", source
    )


def test_manuscript_wording_generator_is_rejected() -> None:
    phrase = "Suggested " + "manuscript wording"
    source = "Path('proposal.md')." + f"write_text({phrase!r}, encoding='utf-8')"
    assert "PROHIBITED_MANUSCRIPT_WORDING_GENERATOR" in categories(
        "scripts/proposal_export.py", source
    )


def test_publication_table_and_legend_authoring_is_rejected() -> None:
    phrase = "Publication-ready " + "tables for Bioinformatics Draft"
    source = "Path('tables_and_legends.md')." + f"write_text({phrase!r}, encoding='utf-8')"
    assert "PROHIBITED_MANUSCRIPT_WORDING_GENERATOR" in categories(
        "scripts/table_export.py", source
    )


def test_manuscript_media_figure_compositor_is_rejected() -> None:
    media = "manuscript/" + "_extracted_imgs/word/media/image1.png"
    source = "\n".join(
        [
            f"original = Image.open({media!r})",
            "figure." + "savefig('figure1.png')",
        ]
    )
    assert "PROHIBITED_MANUSCRIPT_FIGURE_LAYOUT_SCRIPT" in categories(
        "reproducibility/build_figure.py", source
    )


def test_read_only_document_validator_is_allowed() -> None:
    target = "Main_" + "Manuscript.docx"
    source = "\n".join(
        [
            "from " + "docx import Document",
            f"document = Document({target!r})",
            "print(len(document.paragraphs))",
        ]
    )
    assert categories("scripts/validate_submission_snapshot.py", source) == set()


def test_reproducibility_report_writer_is_allowed() -> None:
    source = "\n".join(
        [
            "from " + "docx import Document",
            "document = Document()",
            "document." + "save('Reproducibility_Audit_Report.docx')",
        ]
    )
    assert categories("reproducibility/generate_audit_report.py", source) == set()


def test_submission_document_artifact_is_rejected() -> None:
    relative = "deliverables/" + "Main_" + "Manuscript.docx"
    assert "PROHIBITED_MANUSCRIPT_ARTIFACT" in categories(relative)


def test_rendered_figure_revision_artifact_is_rejected() -> None:
    relative = "reproducibility/" + "figure1b_" + "revision/figure1_revised.png"
    assert "PROHIBITED_MANUSCRIPT_REVISION_RENDER" in categories(relative)


def test_private_manuscript_source_metadata_is_rejected() -> None:
    media = "manuscript/" + "_extracted_imgs/word/media/image1.png"
    source = '{"original_figure": ' + repr(media) + "}"
    assert "PROHIBITED_MANUSCRIPT_SOURCE_METADATA" in categories(
        "reproducibility/figure_revision_metadata.json", source
    )


def test_package_integrity_manifests_match_public_tree() -> None:
    if not (ROOT / ".git").exists():
        pytest.skip("package manifests describe the complete public source checkout")
    passed, stale = verify(Path(ROOT))
    assert passed, stale


def test_docker_context_excludes_local_manuscript_production_material() -> None:
    if not (ROOT / ".dockerignore").is_file():
        pytest.skip("Docker does not copy its context-control file into the runtime image")
    patterns = set((ROOT / ".dockerignore").read_text(encoding="utf-8").splitlines())
    required = {
        "manuscript/",
        "manuscript_remediation/",
        "reproducibility/round5_analysis/figure1b_revision/",
        "scripts/*manuscript*.py",
        "scripts/*docx*.py",
        "scripts/*change_list*.py",
        "scripts/*checklist*.py",
        "scripts/*remediation_report*.py",
        "scripts/create_render_contact_sheets.py",
        "*.docx",
        "*.docm",
    }
    assert required <= patterns
