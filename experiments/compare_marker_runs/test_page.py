"""Read-only upload/page tests using committed bundles, without expression data."""
import io
import json
from pathlib import Path
import zipfile

import pytest
import streamlit as st
from streamlit.testing.v1 import AppTest

from upload import BundleError,compare_uploads,dataset_notice
from compare import compare

ROOT=Path(__file__).resolve().parent


def bundle(directory,changes=None,omit=()):
    output=io.BytesIO()
    with zipfile.ZipFile(output,'w') as archive:
        for p in directory.iterdir():
            if p.name not in omit:
                archive.writestr(p.name,(changes or {}).get(p.name,p.read_bytes()))
    return output.getvalue()


@pytest.mark.parametrize('later,expected',[
    ('sample/B',dict(added=7,removed=7,reranked=12,unchanged=1)),
    ('extended/top_n/B',dict(added=0,removed=10,reranked=0,unchanged=10)),
    ('extended/combined/B',dict(added=0,removed=10,reranked=9,unchanged=1)),
])
def test_upload_parity(later,expected):
    result,shared=compare_uploads(bundle(ROOT/'sample/A'),bundle(ROOT/later))
    assert result==compare(ROOT/'sample/A',ROOT/later)
    assert result['counts']==expected and shared==['B cells']


def test_metadata_and_validation():
    earlier=bundle(ROOT/'sample/A')
    for invalid in [b'not zip',bundle(ROOT/'sample/B',omit=['marker_gene_results.csv']),
                    bundle(ROOT/'sample/B',changes={'analysis_config.json':b'[]'})]:
        with pytest.raises(BundleError):compare_uploads(earlier,invalid)
    malicious=io.BytesIO()
    with zipfile.ZipFile(malicious,'w') as archive:archive.writestr('../bad','bad')
    with pytest.raises(BundleError,match='flat'):compare_uploads(earlier,malicious.getvalue())
    nohash=json.loads((ROOT/'sample/B/provenance.json').read_text());nohash.pop('input_sha256')
    report,_=compare_uploads(earlier,bundle(ROOT/'sample/B',{'provenance.json':json.dumps(nohash)}))
    assert 'cannot be verified' in dataset_notice(report)
    assert not any(row['explainable'] for row in report['markers'])


def test_page_counts_and_stale_clear(monkeypatch):
    uploads={'Earlier run':bundle(ROOT/'sample/A'),'Later run':bundle(ROOT/'sample/B')}
    def uploaded(label,**kwargs):
        x=io.BytesIO(uploads[label]);x.name=label+'.zip';return x
    monkeypatch.setattr(st,'file_uploader',uploaded)
    page=AppTest.from_file(ROOT/'streamlit_app.py').run(timeout=30)
    assert not page.exception and not page.error
    assert page.table[-1].value.iloc[0].to_dict()==dict(Added=7,Removed=7,Reranked=12,Unchanged=1)
    added=page.dataframe[0].value
    assert len(added)==7 and set(added.Evidence)=={'Cause unavailable from these exports'}
    assert set(added['Earlier rank'])=={'Unavailable'}
    uploads['Later run']=b'broken'
    # Changing bytes without widget callbacks also invalidates the cached result.
    page.run()
    assert page.error and not page.dataframe and not page.table


def test_page_initial():
    page=AppTest.from_file(ROOT/'streamlit_app.py').run()
    assert not page.exception
    assert 'Upload two existing MarkerGeneFinderMini export bundles to compare.' in [x.value for x in page.info]


def test_identical_and_no_shared_page(monkeypatch):
    earlier=bundle(ROOT/'sample/A')
    uploads={'Earlier run':earlier,'Later run':earlier}
    def uploaded(label,**kwargs):
        item=io.BytesIO(uploads[label]);item.name=label+'.zip';return item
    monkeypatch.setattr(st,'file_uploader',uploaded)
    page=AppTest.from_file(ROOT/'streamlit_app.py').run(timeout=30)
    assert not page.exception
    assert any('Identical exported runs' in x.value for x in page.info)
    assert any('No exported shortlist changes' in x.value for x in page.success)
    csv=(ROOT/'sample/B/marker_gene_results.csv').read_text().replace('B cells','Other cluster')
    uploads['Later run']=bundle(ROOT/'sample/B',{'marker_gene_results.csv':csv})
    page.run()
    assert not page.exception
    assert any('no shared clusters' in x.value for x in page.warning)
    assert not page.selectbox


def test_plain_language_evidence_matches_engine():
    from upload import reason_label
    report,_=compare_uploads(bundle(ROOT/'sample/A'),bundle(ROOT/'extended/combined/B'))
    rows=[r for r in report['markers'] if r['change']!='unchanged']
    assert sum(reason_label(r)=='Below the new expression threshold' for r in rows)==7
    assert sum(reason_label(r)=='Outside the new top N' for r in rows)==3
    assert sum(reason_label(r)=='Rank changed after filtering earlier entries' for r in rows)==9
    mismatch,_=compare_uploads(bundle(ROOT/'sample/A'),bundle(ROOT/'extended/context_mismatch/B'))
    assert all(reason_label(r)=='Cause unavailable from these exports' for r in mismatch['markers'] if r['change']!='unchanged')


def test_pancreas_overview_and_layer(monkeypatch):
    paths={'Earlier run':ROOT/'pancreas/results/earlier.zip','Later run':ROOT/'pancreas/results/later.zip'}
    def uploaded(label,**kwargs):
        x=io.BytesIO(paths[label].read_bytes());x.name=label+'.zip';return x
    monkeypatch.setattr(st,'file_uploader',uploaded)
    page=AppTest.from_file(ROOT/'streamlit_app.py').run(timeout=30)
    assert not page.exception and not page.error
    overview=page.table[-2].value
    assert overview.Cluster.tolist()==['beta','duct','alpha']
    assert overview['Changed entries'].tolist()==[28,25,23]
    assert page.table[1].value.Layer.tolist()==['log1p_cp10k','log1p_cp10k']
    page.selectbox(key='shared_cluster').set_value('beta').run()
    assert not page.exception
    assert page.table[-1].value.iloc[0].to_dict()==dict(Added=9,Removed=9,Reranked=10,Unchanged=1)
    added=page.dataframe[0].value
    assert len(added)==9 and set(added.Evidence)=={'Cause unavailable from these exports'}
