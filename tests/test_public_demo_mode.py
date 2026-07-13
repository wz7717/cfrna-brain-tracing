from __future__ import annotations

import importlib
import inspect
import sys

import pandas as pd


def test_main_exposes_tracing_and_sample_upload_with_lazy_imports(monkeypatch) -> None:
    monkeypatch.setenv("CFRNA_PUBLIC_DEMO", "1")
    sys.modules.pop("app.main", None)
    sys.modules.pop("app.pages.tracing_page", None)

    main = importlib.import_module("app.main")

    assert sorted(main.PAGES) == ["overview", "samples", "upload"]
    assert main.NAV_ORDER == ["upload", "samples", "overview"]
    assert "app.pages.tracing_page" not in sys.modules

    func = main.resolve_page_func(main.PAGES["overview"]["func"])
    assert func.__name__ == "display_source_tracing"
    upload_func = main.resolve_page_func(main.PAGES["upload"]["func"])
    assert upload_func.__name__ == "display_data_upload"
    samples_func = main.resolve_page_func(main.PAGES["samples"]["func"])
    assert samples_func.__name__ == "display_sample_list"


def test_sidebar_navigation_uses_a_single_widget_rerun(monkeypatch) -> None:
    monkeypatch.setenv("CFRNA_PUBLIC_DEMO", "1")
    main = importlib.import_module("app.main")

    source = inspect.getsource(main._render_sidebar)
    assert "st.rerun" not in source
    assert "on_click=_set_page" in source

    main.st.session_state.clear()
    main._set_page("samples")
    assert main.st.session_state["page"] == "samples"

    main._set_page("unknown")
    assert main.st.session_state["page"] == "samples"


def test_public_demo_tracing_enters_demo_before_backend_init(monkeypatch) -> None:
    monkeypatch.setenv("CFRNA_PUBLIC_DEMO", "1")
    tracing_page = importlib.import_module("app.pages.tracing_page")
    called = {"demo": False}

    monkeypatch.setattr(tracing_page, "render_page_hero", lambda *args, **kwargs: None)
    monkeypatch.setattr(tracing_page, "get_database_mode", lambda: "rhesus")
    monkeypatch.setattr(tracing_page, "database_label", lambda *args, **kwargs: "Rhesus")
    monkeypatch.setattr(tracing_page, "_render_public_demo_tracing", lambda: called.__setitem__("demo", True))
    monkeypatch.setattr(
        tracing_page,
        "init_processor",
        lambda: (_ for _ in ()).throw(AssertionError("public demo should not initialize processor")),
    )
    tracing_page.display_source_tracing()

    assert called["demo"]


def test_missing_sample_table_falls_back_to_empty_sample_list() -> None:
    tracing_page = importlib.import_module("app.pages.tracing_page")

    class EmptyDatabase:
        @staticmethod
        def get_all_samples():
            raise pd.errors.DatabaseError("no such table: cfrna_samples")

    assert tracing_page._get_all_samples_or_empty(EmptyDatabase()).empty
