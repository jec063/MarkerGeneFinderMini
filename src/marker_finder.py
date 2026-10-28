"""Compatibility entry point for the marker_gene_finder v2 implementation."""

from __future__ import annotations

import sys

if __name__ == "__main__":
    from marker_gene_finder._core import main

    main()
else:
    from marker_gene_finder import _core

    # Preserve historical imports and monkeypatch behavior while the old module
    # path remains supported during the v2 deprecation window.
    sys.modules[__name__] = _core
