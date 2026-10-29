"""Compatibility import for the v2 package version."""

from ._compat import warn_legacy_module

warn_legacy_module("src.version", "marker_gene_finder.version")

from marker_gene_finder.version import APP_VERSION

__all__ = ["APP_VERSION"]
