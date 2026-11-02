"""Reproducible v2.8 release benchmark; run from the repository root.

Each scenario runs in a fresh process, warms once, then measures five runs.
Peak RSS includes interpreter, dependencies, input generation, and warmup.
"""
import argparse
import json
import os
from pathlib import Path
import platform
import resource
import statistics
import subprocess
import sys
import tempfile
import time

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

from marker_gene_finder.analysis import run
from marker_gene_finder.batch import run_batch
from marker_gene_finder.provenance import dependency_versions


def dataset(cells=600, genes=300):
    rng = np.random.default_rng(280)
    labels = np.array(["A", "B", "C"])[np.arange(cells) % 3]
    values = rng.poisson(0.2, (cells, genes)).astype(float)
    for index, label in enumerate(["A", "B", "C"]):
        values[labels == label, index * 4:index * 4 + 4] += 4
    return ad.AnnData(
        sparse.csr_matrix(values),
        obs=pd.DataFrame({"cluster": pd.Categorical(labels)},
                         index=[f"cell{i}" for i in range(cells)]),
        var=pd.DataFrame(index=[f"gene{i:03}" for i in range(genes)]),
    )


def measure(case, cells=600, genes=300):
    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory)
        inputs = root / "inputs"
        inputs.mkdir()
        data = dataset(cells, genes)
        h5ad = inputs / "expression.h5ad"
        data.write_h5ad(h5ad)
        frame = pd.DataFrame(data.X.toarray(), columns=data.var_names)
        frame.insert(0, "cluster", data.obs.cluster.to_numpy())
        csv = inputs / "expression.csv"
        frame.to_csv(csv, index=False)
        def execute():
            if case == "batch":
                records = run_batch(inputs, root / "batch", top_n=5)
                assert all(row["status"] == "complete" for row in records)
            else:
                run(csv if case == "csv" else h5ad, root / "markers.csv",
                    engine="scanpy" if case == "scanpy" else "native", top_n=5)
        execute()
        durations = []
        for _ in range(5):
            start = time.perf_counter()
            execute()
            durations.append(time.perf_counter() - start)
        return {"median_seconds": statistics.median(durations),
                "seconds": durations,
                "peak_rss_kib": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--case", choices=["csv", "h5ad", "scanpy", "batch"])
    parser.add_argument("--output", type=Path)
    parser.add_argument("--cells", type=int, default=600)
    parser.add_argument("--genes", type=int, default=300)
    args = parser.parse_args()
    if args.case:
        print(json.dumps(measure(args.case, args.cells, args.genes)))
    else:
        report = {"python": platform.python_version(), "platform": platform.platform(),
                  "cells": args.cells, "genes": args.genes, "seed": 280, "repeats": 5,
                  "dependencies": dependency_versions(), "cases": {}}
        for case in ["csv", "h5ad", "scanpy", "batch"]:
            result = subprocess.run([sys.executable, __file__, "--case", case,
                                     "--cells", str(args.cells), "--genes", str(args.genes)],
                                    env={**os.environ, "NUMBA_NUM_THREADS": "1",
                                         "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1"},
                                    capture_output=True, text=True, check=True)
            report["cases"][case] = json.loads(result.stdout)
        args.output.write_text(json.dumps(report, indent=2) + "\n")
