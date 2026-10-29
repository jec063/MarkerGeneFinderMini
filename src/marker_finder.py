"""Compatibility entry point for the marker_gene_finder v2 implementation."""

from __future__ import annotations

import sys
import warnings

if __name__ == "__main__":
    warnings.warn(
        "python src/marker_finder.py is deprecated; use the installed "
        "marker-gene-finder command. The legacy command will be removed in 3.0.",
        FutureWarning,
        stacklevel=1,
    )
    from marker_gene_finder._core import main

    main()
else:
    from ._compat import warn_legacy_module

    warn_legacy_module("src.marker_finder", "marker_gene_finder.analysis")
    from marker_gene_finder import _core

    # Preserve historical imports and monkeypatch behavior while the old module
    # path remains supported during the v2 deprecation window.
    sys.modules[__name__] = _core
