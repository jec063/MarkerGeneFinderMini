"""Shared cluster progress contract for Native, Scanpy, and legacy CLI.

Callbacks run on the caller thread, once per completed cluster in selection
order, with one-based completion counts. Callback exceptions propagate.
"""
from collections.abc import Callable
import sys

ProgressCallback = Callable[[int, int, object], None]


def report_progress(callback, completed, total, cluster):
    if callback is not None:
        callback(completed, total, cluster)


def console_progress(completed: int, total: int, cluster: object) -> None:
    print(f"[{completed}/{total}] Completed cluster {cluster}",
          file=sys.stderr, flush=True)
