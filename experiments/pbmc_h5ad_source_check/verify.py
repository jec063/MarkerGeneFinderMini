"""Compare PBMC raw-source app exports, CLI and explicit direct Scanpy."""
import argparse
import hashlib
import json
import zipfile
from pathlib import Path
from importlib.metadata import version

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

from marker_gene_finder import AnalysisConfig
from marker_gene_finder.cli import main
from marker_gene_finder.expression_summary import summarize_expression

URL='https://raw.githubusercontent.com/chanzuckerberg/cellxgene/main/example-dataset/pbmc3k.h5ad'
SHA='0db367b991dd95809732b218539ede489bea99113807f62ebd7ccc970025fe38'
parser=argparse.ArgumentParser(__doc__)
parser.add_argument('directory',type=Path,help='Contains downloaded pbmc3k.h5ad and app exports')
args=parser.parse_args();root=args.directory
source=root/'pbmc3k.h5ad'
assert hashlib.sha256(source.read_bytes()).hexdigest()==SHA
x=ad.read_h5ad(source)
assert x.shape==(2638,1838) and x.raw.shape==(2638,13714)
assert not x.layers and x.obs.louvain.nunique()==8
assert x.X.min()<0 and x.raw.X.min()>=0
# Explicit raw selection; no implicit Scanpy default and no stored rankings reused.
sc.tl.rank_genes_groups(x,groupby='louvain',groups=['B cells'],reference='rest',
                       use_raw=True,method='wilcoxon',corr_method='benjamini-hochberg',
                       n_genes=x.raw.n_vars,pts=True,tie_correct=True)
r=sc.get.rank_genes_groups_df(x,group='B cells').set_index('names')
app=pd.read_csv(root/'app-markers.csv')
settings=json.loads((root/'app-settings.json').read_text())['settings']
with zipfile.ZipFile(root/'app-bundle.zip') as bundle:
    assert bundle.testzip() is None
    assert json.loads(bundle.read('analysis_settings.json'))['settings']==settings
assert settings['engine']=='scanpy' and settings['use_raw'] is True
assert settings['expression_source']=='AnnData.raw.X' and settings['layer'] is None
assert settings['target_clusters']==['B cells']
assert settings['rank_by']=='log2FC' and settings['top_n']==5
for row in app.itertuples():
 direct=r.loc[row.gene]
 for a,b in [('log2FC','logfoldchanges'),('p_value','pvals'),('p_adj','pvals_adj'),('pct_in','pct_nz_group'),('pct_out','pct_nz_reference')]:
  np.testing.assert_allclose(getattr(row,a),direct[b],rtol=1e-6,atol=1e-12)
# Apply actual permissive app filters, including mean and detection differences.
mask=np.asarray(x.obs.louvain=='B cells');m=x.raw.X
means=np.asarray(m[mask].mean(0)).ravel();outside=np.asarray(m[~mask].mean(0)).ravel()
rr=r.reindex(x.raw.var_names).copy()
rr['mean_difference']=means-outside
eligible=rr[(rr.logfoldchanges>=0)&(rr.pvals_adj<=1)&(rr.pct_nz_group>=rr.pct_nz_reference)&(rr.mean_difference>=0)]
expected=eligible.reset_index().rename(columns={'index':'gene','names':'gene'})
expected=expected.sort_values(['logfoldchanges','gene'],ascending=[False,True]).head(5)
assert app.gene.tolist()==expected.gene.tolist()
# Check exported dot/heatmap values directly on .raw, not scaled X.
summary=pd.read_csv(root/'app-summary.csv')
reused=summarize_expression(x,summary.gene.drop_duplicates(),cluster_column='louvain',use_raw=True)
for row in summary.itertuples():
 vals=x.raw.X[np.asarray(x.obs.louvain==row.cluster),x.raw.var_names.get_loc(row.gene)].toarray().ravel()
 np.testing.assert_allclose(row.mean_expression,vals.mean(),rtol=1e-6,atol=1e-8)
 np.testing.assert_allclose(row.fraction_expressing,(vals>0).mean(),atol=1e-12)
 one=reused[(reused.cluster==row.cluster)&(reused.gene==row.gene)].iloc[0]
 np.testing.assert_allclose(row.mean_expression,one.mean_expression,rtol=1e-6)
config=AnalysisConfig(engine='scanpy',cluster_column='louvain',use_raw=True,target_clusters=('B cells',))
config.to_json(root/'cli-config.json')
assert main(['analyze',str(source),str(root/'cli.csv'),'--config',str(root/'cli-config.json'),
             '--metadata',str(root/'cli-metadata.json')])==0
pd.testing.assert_frame_equal(app,pd.read_csv(root/'cli.csv'),check_dtype=False,rtol=1e-6,atol=1e-12)
report={'url':URL,'sha256':SHA,'shape_X':list(x.shape),'shape_raw':list(x.raw.shape),'layers':list(x.layers),
        'counts':x.obs.louvain.value_counts().to_dict(),'X_min':float(x.X.min()),'raw_min':float(x.raw.X.min()),
        'X_max':float(x.X.max()),'raw_max':float(x.raw.X.max()),'preprocessing_metadata':'No uns[log1p]; state unspecified by metadata; numerical ranges alone do not prove normalization.',
        'genes_verified':app.gene.tolist(),'plot_summary_rows_verified':len(summary),
        'scanpy_version':version('scanpy'),'anndata_version':version('anndata'),
        'direct_statistics_match':True,'filtered_rank_match':True,'cli_match':True,'source_settings_match':True}
(root/'verification.json').write_text(json.dumps(report,indent=2)+'\n')
print(app[['gene','log2FC','p_value','p_adj','pct_in','pct_out']].to_string(index=False))
print('All direct statistics, ranking, plot summaries, source settings and CLI checks passed.')
