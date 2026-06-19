from __future__ import annotations

from .dao import get_atlas_options, get_sigset_options
from signature_builder import get_latest_sigset_id as _get_latest_sigset_id


def list_atlases(db_path: str):
    return get_atlas_options(db_path)


def list_signature_sets(db_path: str, atlas_id: int):
    return get_sigset_options(db_path, atlas_id)


def get_latest_sigset_id(db_path: str, atlas_id: int = 1):
    return _get_latest_sigset_id(db_path, atlas_id)
