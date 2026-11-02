"""Publish complete files using a unique sibling and an atomic replacement.

Atomicity is per file on filesystems supporting os.replace; this is not a
multi-file transaction or a guarantee against power loss. A hard process kill
may leave a hidden temporary file, but never publishes it as the destination.
"""
from contextlib import contextmanager
import os
from pathlib import Path
import tempfile


@contextmanager
def atomic_path(path):
    """Yield a staging path; replace the destination only after success."""
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    descriptor, name = tempfile.mkstemp(
        prefix=f".{destination.name}.", suffix=".tmp", dir=destination.parent
    )
    os.close(descriptor)
    temporary = Path(name)
    try:
        yield temporary
        os.replace(temporary, destination)
    finally:
        temporary.unlink(missing_ok=True)


def write_text(path, contents, encoding="utf-8"):
    with atomic_path(path) as temporary:
        temporary.write_text(contents, encoding=encoding)


def write_bytes(path, contents):
    with atomic_path(path) as temporary:
        temporary.write_bytes(contents)


def write_csv(path, markers):
    with atomic_path(path) as temporary:
        markers.to_csv(temporary, index=False)
