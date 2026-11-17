"""Export-only regression checks; no expression data or network required."""
import json
from pathlib import Path
import shutil

from compare import compare

ROOT=Path(__file__).resolve().parent


def test_original_threshold_and_identical():
    result=compare(ROOT/'sample/A',ROOT/'sample/B')
    assert compare(ROOT/'extended/A.zip',ROOT/'sample/B.zip')==result
    assert result['counts']==dict(added=7,removed=7,reranked=12,unchanged=1)
    assert result['explainable_changes']==dict(added=0,removed=7,reranked=12)
    unchanged=compare(ROOT/'sample/A',ROOT/'sample/A')
    assert unchanged['counts']==dict(added=0,removed=0,reranked=0,unchanged=20)


def test_top_n_and_combined():
    top=compare(ROOT/'extended/A.zip',ROOT/'extended/top_n/B.zip')
    assert top['counts']==dict(added=0,removed=10,reranked=0,unchanged=10)
    assert top['explainable_changes']['removed']==10
    both=compare(ROOT/'extended/A.zip',ROOT/'extended/combined/B.zip')
    assert both['counts']==dict(added=0,removed=10,reranked=9,unchanged=1)
    assert both['explainable_changes']==dict(added=0,removed=10,reranked=9)
    removed=[r for r in both['markers'] if r['change']=='removed']
    assert sum(r['evidence'].startswith('threshold crossing') for r in removed)==7
    assert sum(r['evidence'].startswith('top_n truncation') for r in removed)==3


def test_mismatched_context():
    result=compare(ROOT/'sample/A',ROOT/'extended/context_mismatch/B')
    assert result['context_status']=='different'
    assert 'input_sha256' in result['context_changes']
    assert not any(r['explainable'] for r in result['markers'])


def test_grouping_and_missing_metadata_block_causes(tmp_path):
    copied=tmp_path/'B';shutil.copytree(ROOT/'sample/B',copied)
    path=copied/'analysis_config.json'
    config=json.loads(path.read_text());config['cluster_column']='different_grouping'
    path.write_text(json.dumps(config))
    result=compare(ROOT/'sample/A',copied)
    assert 'cluster_column' in result['context_changes']
    assert not any(r['explainable'] for r in result['markers'])
    shutil.rmtree(copied);shutil.copytree(ROOT/'sample/B',copied)
    p=copied/'provenance.json';provenance=json.loads(p.read_text());del provenance['engine_version'];p.write_text(json.dumps(provenance))
    assert not any(r['explainable'] for r in compare(ROOT/'sample/A',copied)['markers'])


def test_inconsistent_exported_settings_block_causes(tmp_path):
    copied=tmp_path/'B';shutil.copytree(ROOT/'extended/top_n/B',copied)
    p=copied/'provenance.json';provenance=json.loads(p.read_text());provenance['effective_config']['top_n']=15;p.write_text(json.dumps(provenance))
    assert not any(r['explainable'] for r in compare(ROOT/'sample/A',copied)['markers'])
