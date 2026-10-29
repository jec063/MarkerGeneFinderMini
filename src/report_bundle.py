"""Compatibility imports for the v2 report-bundle API."""

from ._compat import warn_legacy_module

warn_legacy_module("src.report_bundle", "marker_gene_finder.report_bundle")

from marker_gene_finder.report_bundle import *  # noqa: F403
