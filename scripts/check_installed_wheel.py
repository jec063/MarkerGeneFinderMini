"""Run from outside the checkout with a wheel installed and dependencies present.

This script and its inputs must first be copied outside the source checkout.
"""
import argparse
import importlib.abc
from importlib import metadata
import os
from pathlib import Path
import runpy
import subprocess
import sys


class NoLegacyImports(importlib.abc.MetaPathFinder):
    def find_spec(self, fullname, path=None, target=None):
        if fullname == "src" or fullname.startswith("src.") or fullname in {
            "marker_gene_finder._core", "marker_gene_finder._legacy_cli",
        }:
            raise ImportError(f"Installed-wheel check forbids legacy import: {fullname}")
        return None


def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--checkout", type=Path, required=True)
    args = parser.parse_args()
    checkout = args.checkout.resolve()
    assert not Path.cwd().resolve().is_relative_to(checkout)
    sys.meta_path.insert(0, NoLegacyImports())
    import marker_gene_finder
    package = Path(marker_gene_finder.__file__).resolve()
    assert not package.is_relative_to(checkout), package
    assert metadata.version("marker-gene-finder-mini") == marker_gene_finder.APP_VERSION
    root = Path(__file__).resolve().parent
    runpy.run_path(str(root / "migration.py"))["run_demo"](root / "results")
    from streamlit.testing.v1 import AppTest
    app = AppTest.from_file(str(root / "app.py"), default_timeout=30).run()
    assert not app.exception, app.exception
    app.checkbox(key="use_example_data").check().run()
    next(button for button in app.button if button.label == "Find marker genes").click().run()
    assert not app.exception and not app.error
    assert not app.session_state["marker_results"].empty
    assert not any(name == "src" or name.startswith("src.") for name in sys.modules)
    # Apply the same guard to fresh Python processes, including the console script.
    (root / "sitecustomize.py").write_text(
        "import sys\nfrom check_installed_wheel import NoLegacyImports\n"
        "sys.meta_path.insert(0, NoLegacyImports())\n"
    )
    environment = dict(os.environ)
    environment["PYTHONPATH"] = str(root) + os.pathsep + environment.get("PYTHONPATH", "")
    for command in ([sys.executable, "-m", "marker_gene_finder", "--version"],
                    ["marker-gene-finder", "list-engines"]):
        subprocess.run(command, check=True, env=environment)
    print(f"Installed wheel, migration examples, and Streamlit passed: {package}")


if __name__ == "__main__":
    main()
