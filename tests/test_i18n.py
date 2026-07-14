from types import SimpleNamespace

from app import i18n


def test_english_is_the_default_language(monkeypatch):
    monkeypatch.setattr(i18n, "st", SimpleNamespace(session_state={}))

    assert i18n.get_language_mode() == "en"
    assert i18n.tr("中文", "English") == "English"


def test_invalid_language_falls_back_to_english(monkeypatch):
    state = {}
    monkeypatch.setattr(i18n, "st", SimpleNamespace(session_state=state))

    i18n.set_language_mode("unsupported")

    assert state["ui_language_mode"] == "en"
    assert i18n.get_language_mode() == "en"
