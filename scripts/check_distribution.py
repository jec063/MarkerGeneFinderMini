"""Reject legacy modules in built wheel and source archives."""
from pathlib import Path, PurePosixPath
import sys
import tarfile
from zipfile import ZipFile


def check(path):
    if path.suffix == ".whl":
        with ZipFile(path) as archive:
            names = archive.namelist()
    else:
        with tarfile.open(path) as archive:
            names = archive.getnames()
    for name in names:
        parts = PurePosixPath(name).parts
        assert "src" not in parts, (path, name)
        assert not name.endswith(("/marker_gene_finder/_core.py", "/marker_gene_finder/_legacy_cli.py")), (path, name)
        assert name not in {"marker_gene_finder/_core.py", "marker_gene_finder/_legacy_cli.py"}, (path, name)
    assert any(name.endswith("marker_gene_finder/__main__.py") for name in names)
    print(f"Verified v3 distribution has no legacy modules: {path.name}")


if __name__ == "__main__":
    for argument in sys.argv[1:]:
        check(Path(argument))
