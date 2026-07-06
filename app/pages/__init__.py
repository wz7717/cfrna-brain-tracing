from __future__ import annotations

from importlib import import_module

_EXPORTS = {
    "display_database_overview": "app.pages.overview_page",
    "display_data_upload": "app.pages.upload_page",
    "display_sample_list": "app.pages.sample_manage_page",
    "display_source_tracing": "app.pages.tracing_page",
    "display_atlas_browser": "app.pages.atlas_page",
    "display_run_compare": "app.pages.compare_runs_page",
    "display_benchmark_page": "app.pages.benchmark_page",
}

__all__ = sorted(_EXPORTS)


def __getattr__(name: str):
    if name not in _EXPORTS:
        raise AttributeError(name)
    return getattr(import_module(_EXPORTS[name]), name)
