"""Real PBMC measurements match each engine's independent v2.9 reference."""
from dataclasses import replace
import hashlib
import json
from pathlib import Path
import zipfile

import anndata as ad
import pandas as pd
import pytest
import scanpy as sc

from marker_gene_finder import AnalysisConfig, analyze_file, InputValidationError
from marker_gene_finder.batch import run_batch, verify_batch_manifest
from marker_gene_finder.cli import main

FIXTURES = Path(__file__).parent / 'fixtures/pbmc68k'


@pytest.fixture(scope='module')
def biological_inputs(tmp_path_factory):
    root = tmp_path_factory.mktemp('pbmc')
    source = Path(sc.__file__).parent / 'datasets/10x_pbmc68k_reduced.h5ad'
    provenance = json.loads((FIXTURES / 'provenance.json').read_text())
    assert hashlib.sha256(source.read_bytes()).hexdigest() == provenance['sha256'], (
        'Scanpy PBMC fixture changed; review provenance and scientific references.'
    )
    original = sc.datasets.pbmc68k_reduced()
    data = ad.AnnData(original.raw.X.astype('float64'),
                      obs=original.obs[['bulk_labels']].copy(), var=original.raw.var.copy())
    assert data.shape == (700, 765)
    data.layers['expression'] = data.X.copy()
    data.raw = data.copy()
    data.write_h5ad(root / 'pbmc.h5ad')
    frame = pd.DataFrame(data.X.toarray(), columns=data.var_names)
    frame.insert(0, 'bulk_labels', data.obs.bulk_labels.to_numpy())
    frame.to_csv(root / 'pbmc.csv', index=False)
    return root


@pytest.mark.parametrize('engine', ['native', 'scanpy'])
def test_pbmc_workflows(engine, biological_inputs, tmp_path):
    reference = pd.read_csv(FIXTURES / f'{engine}.csv')
    config = AnalysisConfig(engine=engine, cluster_column='bulk_labels', top_n=5)
    for label, options in [('X', {}), ('layer', {'layer': 'expression'}), ('raw', {'use_raw': True})]:
        result = analyze_file(biological_inputs / 'pbmc.h5ad', tmp_path / f'{label}.csv',
                              replace(config, **options))
        pd.testing.assert_frame_equal(result.markers, reference, check_dtype=False,
                                      rtol=1e-6, atol=1e-10)
    if engine == 'native':
        result = analyze_file(biological_inputs / 'pbmc.csv', tmp_path / 'csv.csv', config)
        pd.testing.assert_frame_equal(result.markers, reference, check_dtype=False,
                                      rtol=1e-6, atol=1e-10)
    else:
        with pytest.raises(InputValidationError, match='does not support csv'):
            analyze_file(biological_inputs / 'pbmc.csv', tmp_path / 'unsupported.csv', config)
        assert not (tmp_path / 'unsupported.csv').exists()
    config.to_json(tmp_path / 'config.json')
    assert AnalysisConfig.from_json(tmp_path / 'config.json') == config
    assert main(['analyze', str(biological_inputs / 'pbmc.h5ad'), str(tmp_path / 'cli.csv'),
                 '--config', str(tmp_path / 'config.json'), '--bundle', str(tmp_path / 'bundle.zip')]) == 0
    pd.testing.assert_frame_equal(pd.read_csv(tmp_path / 'cli.csv'), reference,
                                  check_dtype=False, rtol=1e-6, atol=1e-10)
    with zipfile.ZipFile(tmp_path / 'bundle.zip') as bundle:
        assert bundle.testzip() is None
    options = dict(engine=engine, cluster_column='bulk_labels', top_n=5,
                   include_patterns=['*.h5ad'] if engine == 'scanpy' else [])
    records = run_batch(biological_inputs, tmp_path / 'batch', **options)
    assert len(records) == (1 if engine == 'scanpy' else 2)
    assert all(row['status'] == 'complete' for row in records)
    resumed = run_batch(biological_inputs, tmp_path / 'batch', resume=True, **options)
    assert all(row['status'] == 'skipped' for row in resumed)
    assert verify_batch_manifest(tmp_path / 'batch/batch_manifest.json')['valid']
