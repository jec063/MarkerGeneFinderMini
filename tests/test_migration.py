"""Executable public migration examples and the frozen v3 surface."""
import ast
import json
from pathlib import Path
import runpy

ROOT = Path(__file__).resolve().parents[1]


def test_migration_examples(tmp_path):
    runpy.run_path(str(ROOT / "examples/migration.py"))["run_demo"](tmp_path)


def test_planned_v3_contract_matches_reviewed_snapshot():
    actual = runpy.run_path(str(ROOT / "scripts/public_contract.py"))["contract"]()
    # JSON normalizes tuple defaults to the portable array representation.
    actual = json.loads(json.dumps(actual))
    assert actual == json.loads((ROOT / "docs/public-api-v3.json").read_text())


def test_primary_code_and_examples_do_not_import_legacy_modules():
    paths = [ROOT / "app.py", *(ROOT / "examples").glob("*.py"),
             *(ROOT / "tests").glob("test_*.py")]
    for path in paths:
        for node in ast.walk(ast.parse(path.read_text())):
            if isinstance(node, ast.Import):
                names = [alias.name for alias in node.names]
            elif isinstance(node, ast.ImportFrom):
                names = [node.module or ""]
                if node.module == "marker_gene_finder":
                    names += ["marker_gene_finder." + alias.name for alias in node.names]
            else:
                continue
            assert not any(name == "src" or name.startswith("src.") or name in {
                "marker_gene_finder._core", "marker_gene_finder._legacy_cli",
            } for name in names), path
