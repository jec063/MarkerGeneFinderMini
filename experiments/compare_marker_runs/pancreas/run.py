"""Generate real multi-cluster pancreas exports and verify comparison evidence."""
import argparse
from dataclasses import replace
import hashlib
import io
import json
from pathlib import Path
import sys
import zipfile

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from compare import compare,render
from marker_gene_finder import AnalysisConfig,analyze_file

ROOT=Path(__file__).resolve().parent
URL='https://datasets.cellxgene.cziscience.com/ac56150b-add4-4336-9059-6d3d3ce17f3b.h5ad'
SHA='183673651cfa8c473a26641d42011d43be44eb2fea44e6e6ab8e2b0065d07483'
TARGETS=('alpha','beta','duct')


def main():
    p=argparse.ArgumentParser(__doc__);p.add_argument('input',type=Path)
    p.add_argument('--output',type=Path,default=ROOT/'results')
    args=p.parse_args();out=args.output;out.mkdir(parents=True,exist_ok=True)
    assert hashlib.sha256(args.input.read_bytes()).hexdigest()==SHA
    raw=ad.read_h5ad(args.input)
    assert raw.shape==(2126,15643)
    data=ad.AnnData(raw.raw.X.copy(),obs=raw.obs.copy(),var=raw.raw.var.copy())
    assert data.X.min()>=0
    np.testing.assert_array_equal(data.X.data,np.floor(data.X.data))
    data.var_names=pd.Index(data.var.feature_name.astype(str),name='gene_symbol')
    data.var_names_make_unique()
    sc.pp.normalize_total(data,target_sum=10000);sc.pp.log1p(data)
    # Store a named analysis layer; never infer preprocessing from raw.
    data.layers['log1p_cp10k']=data.X.copy()
    derived=args.input.parent/'analysis.h5ad';data.write_h5ad(derived)
    config=AnalysisConfig(engine='scanpy',cluster_column='cellular_classification',
        layer='log1p_cp10k',target_clusters=TARGETS,top_n=20,min_pct=.10)
    for name,threshold in [('earlier',.10),('later',.25)]:
        folder=out/name;folder.mkdir(exist_ok=True)
        result=analyze_file(derived,folder/'marker_gene_results.csv',replace(config,min_pct=threshold),hash_input_file=True)
        bundle=result.to_bundle();(out/(name+'.zip')).write_bytes(bundle)
        with zipfile.ZipFile(io.BytesIO(bundle)) as z:
            for member in ['analysis_settings.json','analysis_config.json','provenance.json']:
                (folder/member).write_bytes(z.read(member))
    report=compare(out/'earlier.zip',out/'later.zip')
    assert report['config_changes']=={'min_pct':[.10,.25]} and report['threshold_only']
    (out/'comparison.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    render(report,out/'report.html')
    # Independent full Scanpy result and filters for each annotated target.
    sc.tl.rank_genes_groups(data,groupby='cellular_classification',groups=list(TARGETS),
        reference='rest',layer='log1p_cp10k',use_raw=False,method='wilcoxon',tie_correct=True,
        corr_method='benjamini-hochberg',pts=True,n_genes=data.n_vars)
    full={};evidence=[];summary=[]
    for cluster in TARGETS:
        stats=sc.get.rank_genes_groups_df(data,group=cluster).set_index('names')
        mask=np.asarray(data.obs.cellular_classification==cluster)
        delta=np.asarray(data.X[mask].mean(0)).ravel()-np.asarray(data.X[~mask].mean(0)).ravel()
        stats['mean_difference']=pd.Series(delta,index=data.var_names).reindex(stats.index)
        for name,t in [('earlier',.10),('later',.25)]:
            passing=stats[(stats.pct_nz_group>=t)&(stats.pct_nz_group>=stats.pct_nz_reference)&
              (stats.logfoldchanges>=0)&(stats.pvals_adj<=1)&(stats.mean_difference>=0)]
            full[name]=passing.reset_index().sort_values(['logfoldchanges','names'],ascending=[False,True]).reset_index(drop=True)
            exported=pd.read_csv(out/name/'marker_gene_results.csv');exported=exported[exported.cluster==cluster]
            assert exported.gene.tolist()==full[name].names.head(20).tolist()
            for row in exported.itertuples():
                for field,scfield in [('pct_in','pct_nz_group'),('log2FC','logfoldchanges'),('p_adj','pvals_adj')]:
                    np.testing.assert_allclose(getattr(row,field),stats.loc[row.gene,scfield],rtol=1e-6,atol=1e-12)
        rows=[r for r in report['markers'] if r['cluster']==cluster]
        for r in rows:
            gene=r['gene']
            if r['evidence'].startswith('threshold crossing'):assert .1<=stats.loc[gene,'pct_nz_group']<.25
            if r['rank_from_exported_prefix_B'] is not None:
                assert r['rank_from_exported_prefix_B']==int(full['later'].index[full['later'].names==gene][0]+1)
            if r['change']=='added':
                pos=full['earlier'].index[full['earlier'].names==gene]
                assert len(pos)==1 and pos[0]>=20
                assert not r['explainable'] and r['rank_A'] is None
                evidence.append({'cluster':cluster,'gene':gene,'actual_earlier_postfilter_rank':int(pos[0]+1),
                    'report_knows_earlier_rank':False,'verification_only':'Earlier top20 truncation; not supplied to page'})
        counts={k:sum(r['change']==k for r in rows) for k in ['added','removed','reranked','unchanged']}
        summary.append(dict(cluster=cluster,cells=int(mask.sum()),**counts,
            membership_changes=counts['added']+counts['removed'],
            total_changed_entries=sum(counts[k] for k in ['added','removed','reranked']),
            explainable=sum(r['explainable'] for r in rows)))
    pd.DataFrame(summary).to_csv(out/'cluster_summary.csv',index=False)
    (out/'independent_verification.json').write_text(json.dumps({'source_url':URL,'source_sha256':SHA,
        'source_cells':data.n_obs,'source_features':data.n_vars,'all_source_cells_retained_as_comparisons':True,
        'groups':list(TARGETS),'layer':'log1p_cp10k','preprocessing':'source raw integer counts; normalize_total(10000), log1p; no scaling or imputation',
        'summary':summary,'independent_added_rank_checks':evidence,'direct_scanpy_statistics_and_order_match':True,
        'all_reported_threshold_and_rank_explanations_verified':True},indent=2)+'\n')
    print(pd.DataFrame(summary).to_string(index=False))


if __name__=='__main__':main()
