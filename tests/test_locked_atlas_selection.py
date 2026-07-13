from __future__ import annotations

import pytest

from app.pages import tracing_page


def test_locked_atlas_selection_rejects_missing_atlas(monkeypatch) -> None:
    monkeypatch.setattr(tracing_page, "get_atlas_options", lambda *args, **kwargs: [])
    monkeypatch.setattr(tracing_page, "packaged_formal_region_assets_available", lambda: False)

    with pytest.raises(RuntimeError, match="Bo2023"):
        tracing_page._select_locked_bo2023_atlas("rhesus")


def test_locked_atlas_selection_rejects_non_bo2023_atlas(monkeypatch) -> None:
    monkeypatch.setattr(
        tracing_page,
        "get_atlas_options",
        lambda *args, **kwargs: [(7, "Unrelated reference")],
    )
    monkeypatch.setattr(tracing_page, "packaged_formal_region_assets_available", lambda: False)

    with pytest.raises(RuntimeError, match="Bo2023"):
        tracing_page._select_locked_bo2023_atlas("rhesus")


def test_locked_atlas_selection_uses_matching_reference(monkeypatch) -> None:
    monkeypatch.setattr(
        tracing_page,
        "get_atlas_options",
        lambda *args, **kwargs: [(2, "Other atlas"), (4, "Bo2023_WangLab_VSD_region")],
    )

    assert tracing_page._select_locked_bo2023_atlas("rhesus") == (4, "Bo2023_WangLab_VSD_region")


def test_locked_atlas_selection_uses_packaged_formal_reference_without_db_atlas(monkeypatch) -> None:
    monkeypatch.setattr(tracing_page, "get_atlas_options", lambda *args, **kwargs: [])
    monkeypatch.setattr(tracing_page, "packaged_formal_region_assets_available", lambda: True)

    atlas_id, label = tracing_page._select_locked_bo2023_atlas("rhesus")

    assert atlas_id is None
    assert "Bo2023" in label
