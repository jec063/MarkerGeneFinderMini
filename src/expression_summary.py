"""Compatibility imports for the v2 expression-summary API."""

from ._compat import warn_legacy_module

warn_legacy_module("src.expression_summary", "marker_gene_finder.expression_summary")

from marker_gene_finder.expression_summary import *  # noqa: F403
