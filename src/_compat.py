"""Shared deprecation notice for legacy ``src`` imports."""

from __future__ import annotations

import warnings


def warn_legacy_module(old_name: str, new_name: str) -> None:
    """Warn that a legacy module will be removed in version 3."""
    warnings.warn(
        f"{old_name} is deprecated; import {new_name} instead. "
        "Legacy src modules will be removed in Marker Gene Finder 3.0.",
        DeprecationWarning,
        stacklevel=2,
    )
