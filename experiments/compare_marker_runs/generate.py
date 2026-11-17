"""Create two current-format PBMC exports and verify their comparison."""
import argparse
from dataclasses import replace
import hashlib
import json
from pathlib import Path
import tempfile
import zipfile

import numpy as np
import pandas as pd
import scanpy as sc

from marker_gene_finder import AnalysisConfig, analyze_file
from compare import compare,render

ROOT=Path(__file__).resolve().parent
SHA='0db367b991dd95809732b218539ede489bea99113807f62ebd7ccc970025fe38'


def main():
    p=argparse.ArgumentParser(__doc__);p.add_argument('input',type=Path);p.add_argument('--output',type=Path,default=ROOT/'sample')
    args=p.parse_args();out=args.output;out.mkdir(parents=True,exist_ok=True)
    assert hashlib.sha256(args.input.read_bytes()).hexdigest()==SHA
    config=AnalysisConfig(engine='scanpy',use_raw=True,cluster_column='louvain',target_clusters=('B cells',),top_n=20,min_pct=.05)
    for name,threshold in [('A',.05),('B',.10)]:
        directory=out/name;directory.mkdir(exist_ok=True)
        result=analyze_file(args.input,directory/'marker_gene_results.csv',replace(config,min_pct=threshold),hash_input_file=True)
        bundle=result.to_bundle()
        (out/(name+'.zip')).write_bytes(bundle)
        with zipfile.ZipFile(__import__('io').BytesIO(bundle)) as archive:
            for member in ['analysis_settings.json','analysis_config.json','provenance.json']:
                (directory/member).write_bytes(archive.read(member))
    report=compare(out/'A',out/'B');render(report,out/'report.html')
    (out/'comparison.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    assert report['config_changes']=={'min_pct':[.05,.1]} and report['threshold_only']
    removed=[r for r in report['markers'] if r['change']=='removed']
    assert removed and all(.05<=r['pct_in_A']<.1 and 'threshold crossing' in r['evidence'] for r in removed)
    identical=compare(out/'A',out/'A');render(identical,out/'identical.html')
    assert all(identical['counts'][k]==0 for k in ['added','removed','reranked'])
    assert not identical['config_changes']
    # Missing provenance must prevent causal attribution, even on the same CSV.
    with tempfile.TemporaryDirectory() as tmp:
        left=Path(tmp)
        for file in (out/'A').iterdir(): (left/file.name).write_bytes(file.read_bytes())
        provenance=json.loads((left/'provenance.json').read_text());provenance['input_sha256']=None
        (left/'provenance.json').write_text(json.dumps(provenance))
        unknown=compare(left,out/'B')
        assert not unknown['threshold_only']
        assert all(r['evidence']=='unresolved' for r in unknown['markers'] if r['change']=='removed')
        # A row without detection statistics cannot establish a threshold crossing.
        (left/'provenance.json').write_bytes((out/'A'/'provenance.json').read_bytes())
        missing=pd.read_csv(left/'marker_gene_results.csv').drop(columns=['pct_in'])
        missing.to_csv(left/'marker_gene_results.csv',index=False)
        unknown_stats=compare(left,out/'B')
        assert all(r['evidence']=='unresolved' for r in unknown_stats['markers'] if r['change']=='removed')
    # Independent full direct Scanpy results verify filtering and export truncation.
    data=sc.read_h5ad(args.input)
    sc.tl.rank_genes_groups(data,groupby='louvain',groups=['B cells'],reference='rest',use_raw=True,
        method='wilcoxon',tie_correct=True,corr_method='benjamini-hochberg',pts=True,n_genes=data.raw.n_vars)
    stats=sc.get.rank_genes_groups_df(data,group='B cells').set_index('names')
    mask=np.asarray(data.obs.louvain=='B cells')
    means=np.asarray(data.raw.X[mask].mean(0)).ravel();others=np.asarray(data.raw.X[~mask].mean(0)).ravel()
    stats['mean_difference']=pd.Series(means-others,index=data.raw.var_names).reindex(stats.index)
    full={}
    for name,t in [('A',.05),('B',.10)]:
        selected=stats[(stats.pct_nz_group>=t)&(stats.pct_nz_group>=stats.pct_nz_reference)&
          (stats.logfoldchanges>=0)&(stats.pvals_adj<=1)&(stats.mean_difference>=0)]
        full[name]=selected.reset_index().sort_values(['logfoldchanges','names'],ascending=[False,True]).reset_index(drop=True)
        exported=pd.read_csv(out/name/'marker_gene_results.csv')
        assert exported.gene.tolist()==full[name].names.head(20).tolist()
        for r in exported.itertuples():
            for a,b in [('pct_in','pct_nz_group'),('p_adj','pvals_adj'),('log2FC','logfoldchanges')]:
                np.testing.assert_allclose(getattr(r,a),stats.loc[r.gene,b],rtol=1e-6,atol=1e-12)
    added_verification=[]
    for row in report['markers']:
        if row['change']=='added':
            pos=full['A'].index[full['A'].names==row['gene']]
            assert len(pos)==1 and pos[0]>=20
            added_verification.append({'gene':row['gene'],'verified_full_rank_A':int(pos[0]+1),
                'report_claim':'A rank unavailable in exports; no causal inference',
                'independent_validation':'passed A filter but outside exported top20'})
    verification={'counts':report['counts'],'threshold_crossings_verified':len(removed),
        'identical_export_counts':identical['counts'],'missing_provenance_abstention_passed':True,'missing_statistics_abstention_passed':True,
        'direct_scanpy_matches':True,'added_genes_independent_check':added_verification,
        'input_sha256':SHA,'scanpy_version':__import__('importlib.metadata').metadata.version('scanpy')}
    (out/'verification.json').write_text(json.dumps(verification,indent=2)+'\n')
    print(json.dumps(verification,indent=2))


if __name__=='__main__':main()
