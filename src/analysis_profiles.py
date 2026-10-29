"""Compatibility imports for the v2 analysis-profile API."""

from ._compat import warn_legacy_module

warn_legacy_module("src.analysis_profiles", "marker_gene_finder.profiles")

from marker_gene_finder.profiles import *  # noqa: F403
