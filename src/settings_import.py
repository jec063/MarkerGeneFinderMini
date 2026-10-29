"""Compatibility imports for the v2 settings-import API."""

from ._compat import warn_legacy_module

warn_legacy_module("src.settings_import", "marker_gene_finder.settings")

from marker_gene_finder.settings import *  # noqa: F403
