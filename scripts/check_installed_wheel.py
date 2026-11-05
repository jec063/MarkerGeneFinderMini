"""Run from outside the checkout with a wheel installed and dependencies present.

This script and its inputs must first be copied outside the source checkout.
"""
import argparse
import importlib.util
from importlib import metadata
from pathlib import Path
import runpy
import subprocess
import sys


def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--checkout", type=Path, required=True)
    args = parser.parse_args()
    checkout = args.checkout.resolve()
    assert not Path.cwd().resolve().is_relative_to(checkout)
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
    for module in ("src", "marker_gene_finder._core", "marker_gene_finder._legacy_cli"):
        assert importlib.util.find_spec(module) is None, module
    files = metadata.files("marker-gene-finder-mini")
    assert files is not None
    for path in files:
        assert path.parts[0] != "src", path
        assert str(path) not in {"marker_gene_finder/_core.py", "marker_gene_finder/_legacy_cli.py"}, path
    for command in ([sys.executable, "-m", "marker_gene_finder", "--version"],
                    ["marker-gene-finder", "list-engines"]):
        subprocess.run(command, check=True)
    rejected = subprocess.run(
        ["marker-gene-finder", "--input", "missing.csv", "--output", "unwanted.csv"],
        capture_output=True, text=True,
    )
    assert rejected.returncode == 2 and not Path("unwanted.csv").exists()
    print(f"Installed wheel, migration examples, and Streamlit passed: {package}")


if __name__ == "__main__":
    main()
