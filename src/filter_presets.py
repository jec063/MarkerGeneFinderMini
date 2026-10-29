"""Compatibility imports for the v2 filter-preset API."""

from ._compat import warn_legacy_module

warn_legacy_module("src.filter_presets", "marker_gene_finder.filter_presets")

from marker_gene_finder.filter_presets import *  # noqa: F403
