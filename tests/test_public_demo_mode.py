from __future__ import annotations

import importlib
import sys


def test_public_demo_main_keeps_only_demo_pages_and_lazy_imports(monkeypatch) -> None:
    monkeypatch.setenv("CFRNA_PUBLIC_DEMO", "1")
    sys.modules.pop("app.main", None)
    sys.modules.pop("app.pages.tracing_page", None)

    main = importlib.import_module("app.main")

    assert sorted(main.PAGES) == ["tracing"]
    assert main.NAV_ORDER == ["ANALYSIS"]
    assert "app.pages.tracing_page" not in sys.modules

    func = main.resolve_page_func(main.PAGES["tracing"]["func"])
    assert func.__name__ == "display_source_tracing"


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
