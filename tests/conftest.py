"""Exercise the primary suite as if the v3 legacy removals had happened."""
import importlib.abc
import os
import sys


class NoLegacyImports(importlib.abc.MetaPathFinder):
    def find_spec(self, fullname, path=None, target=None):
        if fullname == "src" or fullname.startswith("src.") or fullname in {
            "marker_gene_finder._core", "marker_gene_finder._legacy_cli",
        }:
            raise ImportError(f"Legacy module blocked by migration check: {fullname}")
        return None


if os.environ.get("MGF_DISABLE_LEGACY") == "1":
    sys.meta_path.insert(0, NoLegacyImports())
