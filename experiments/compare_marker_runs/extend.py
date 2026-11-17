"""Extend PBMC export comparisons: top-N, two settings, and mismatched context."""
import argparse
from collections import Counter
from dataclasses import replace
import hashlib
import io
import json
from pathlib import Path
import zipfile

import numpy as np
import pandas as pd
import scanpy as sc

from marker_gene_finder import AnalysisConfig, analyze_file
from compare import compare, render
from generate import SHA

ROOT=Path(__file__).resolve().parent


def save(report,directory):
    directory.mkdir(parents=True,exist_ok=True)
    (directory/'comparison.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    render(report,directory/'report.html')


def main():
    p=argparse.ArgumentParser(__doc__)
    p.add_argument('input',type=Path)
    p.add_argument('--sample',type=Path,default=ROOT/'sample')
    p.add_argument('--output',type=Path,default=ROOT/'extended')
    args=p.parse_args();out=args.output;out.mkdir(parents=True,exist_ok=True)
    assert hashlib.sha256(args.input.read_bytes()).hexdigest()==SHA
    a=args.sample/'A'
    config=AnalysisConfig.from_json(a/'analysis_config.json')
    assert config.top_n==20 and config.min_pct==.05
    exports={}
    for name,cfg in [('top_n',replace(config,top_n=10)),('combined',replace(config,top_n=10,min_pct=.10))]:
        directory=out/name/'B';directory.mkdir(parents=True,exist_ok=True)
        result=analyze_file(args.input,directory/'marker_gene_results.csv',cfg,hash_input_file=True)
        with zipfile.ZipFile(io.BytesIO(result.to_bundle())) as bundle:
            (out/name/'B.zip').write_bytes(result.to_bundle())
            for member in ['analysis_config.json','analysis_settings.json','provenance.json']:
                (directory/member).write_bytes(bundle.read(member))
        report=compare(a,directory);save(report,out/name)
        assert report['provenance_comparable']
        assert set(report['config_changes'])==({'top_n'} if name=='top_n' else {'top_n','min_pct'})
        assert compare(a,out/name/'B.zip')==report
        exports[name]=directory
    # Current-format ZIP for the same baseline lets one command compare bundles.
    with zipfile.ZipFile(out/'A.zip','w',compression=zipfile.ZIP_DEFLATED) as z:
        for file in sorted(a.iterdir()):
            info=zipfile.ZipInfo(file.name,date_time=(1980,1,1,0,0,0));info.compress_type=zipfile.ZIP_DEFLATED
            z.writestr(info,file.read_bytes())
    assert compare(out/'A.zip',out/'top_n/B.zip')==compare(a,exports['top_n'])
    # Explicit synthetic metadata mismatch, NOT another biological run.
    mismatch=out/'context_mismatch'/'B';mismatch.mkdir(parents=True,exist_ok=True)
    for file in (args.sample/'B').iterdir():(mismatch/file.name).write_bytes(file.read_bytes())
    provenance=json.loads((mismatch/'provenance.json').read_text())
    provenance['input_sha256']='0'*64
    (mismatch/'provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')
    report=compare(a,mismatch);save(report,out/'context_mismatch')
    assert report['context_status']=='different' and 'input_sha256' in report['context_changes']
    assert all(not r['explainable'] and r['evidence']=='unresolved' for r in report['markers'] if r['change']!='unchanged')
    # Additional context checks: matrix, engine, grouping, target and missing metadata.
    variants=[('expression_source','AnnData.X'),('resolved_engine','native'),('engine_version','different')]
    for key,value in variants:
        changed=dict(provenance);changed['input_sha256']=json.loads((args.sample/'B/provenance.json').read_text())['input_sha256']
        changed[key]=value;(mismatch/'provenance.json').write_text(json.dumps(changed))
        check=compare(a,mismatch)
        assert not check['provenance_comparable'] and not any(r['explainable'] for r in check['markers'])
    (mismatch/'provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')
    # Direct full statistics independently verify the new top lists and every explanation.
    data=sc.read_h5ad(args.input)
    sc.tl.rank_genes_groups(data,groupby='louvain',groups=['B cells'],reference='rest',use_raw=True,
        method='wilcoxon',tie_correct=True,corr_method='benjamini-hochberg',pts=True,n_genes=data.raw.n_vars)
    stats=sc.get.rank_genes_groups_df(data,group='B cells').set_index('names')
    mask=np.asarray(data.obs.louvain=='B cells')
    means=np.asarray(data.raw.X[mask].mean(0)).ravel();others=np.asarray(data.raw.X[~mask].mean(0)).ravel()
    stats['mean_difference']=pd.Series(means-others,index=data.raw.var_names).reindex(stats.index)
    checks=[]
    for name,directory in exports.items():
        cfg=AnalysisConfig.from_json(directory/'analysis_config.json')
        passing=stats[(stats.pct_nz_group>=cfg.min_pct)&(stats.pct_nz_group>=stats.pct_nz_reference)&
            (stats.logfoldchanges>=0)&(stats.pvals_adj<=1)&(stats.mean_difference>=0)]
        ordered=passing.reset_index().sort_values(['logfoldchanges','names'],ascending=[False,True]).reset_index(drop=True)
        csv=pd.read_csv(directory/'marker_gene_results.csv')
        assert csv.gene.tolist()==ordered.names.head(cfg.top_n).tolist()
        report=compare(a,directory)
        for row in report['markers']:
            if row['evidence'].startswith('threshold crossing'):
                assert stats.loc[row['gene'],'pct_nz_group']<cfg.min_pct
            if row['rank_from_exported_prefix_B'] is not None:
                fullrank=int(ordered.index[ordered.names==row['gene']][0]+1)
                assert row['rank_from_exported_prefix_B']==fullrank
            if row['evidence'].startswith('top_n truncation'):
                assert row['rank_from_exported_prefix_B']>cfg.top_n
        changed=[r for r in report['markers'] if r['change']!='unchanged']
        checks.append({'comparison':name,'changes':report['counts'],'explainable':report['explainable_changes'],
            'changed_total':len(changed),'explainable_total':sum(r['explainable'] for r in changed),
            'evidence_types':dict(Counter(r['evidence'] for r in changed)),'direct_scanpy_match':True})
    original=compare(a,args.sample/'B')
    checks.insert(0,{'comparison':'original_min_pct','changes':original['counts'],'explainable':original['explainable_changes'],
        'changed_total':sum(original['counts'][k] for k in ['added','removed','reranked']),
        'explainable_total':sum(r['explainable'] for r in original['markers'])})
    verification={'comparisons':checks,'context_mismatch_suppresses_all_causes':True,
        'source_engine_version_mismatch_checks':True,'zip_directory_parity':True,'dataset_sha256':SHA}
    (out/'verification.json').write_text(json.dumps(verification,indent=2)+'\n')
    print(json.dumps(verification,indent=2))


if __name__=='__main__':main()
