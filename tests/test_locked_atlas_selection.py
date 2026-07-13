from __future__ import annotations

import pytest

from app.pages import tracing_page


def test_locked_atlas_selection_rejects_missing_atlas(monkeypatch) -> None:
    monkeypatch.setattr(tracing_page, "get_atlas_options", lambda *args, **kwargs: [])

    with pytest.raises(RuntimeError, match="Bo2023"):
        tracing_page._select_locked_bo2023_atlas("rhesus")


def test_locked_atlas_selection_rejects_non_bo2023_atlas(monkeypatch) -> None:
    monkeypatch.setattr(
        tracing_page,
        "get_atlas_options",
        lambda *args, **kwargs: [(7, "Unrelated reference")],
    )

    with pytest.raises(RuntimeError, match="Bo2023"):
        tracing_page._select_locked_bo2023_atlas("rhesus")


def test_locked_atlas_selection_uses_matching_reference(monkeypatch) -> None:
    monkeypatch.setattr(
        tracing_page,
        "get_atlas_options",
        lambda *args, **kwargs: [(2, "Other atlas"), (4, "Bo2023_WangLab_VSD_region")],
    )

    assert tracing_page._select_locked_bo2023_atlas("rhesus") == (4, "Bo2023_WangLab_VSD_region")
