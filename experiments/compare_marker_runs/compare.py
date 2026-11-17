"""Compare existing MarkerGeneFinderMini bundles or extracted export directories."""
import argparse
import hashlib
import html
import io
import json
from pathlib import Path
import zipfile

import pandas as pd


def load(path):
    path=Path(path)
    if path.is_dir():
        read=lambda name:(path/name).read_bytes()
        names={p.name for p in path.iterdir()}
    else:
        with zipfile.ZipFile(path) as archive:
            contents={name:archive.read(name) for name in archive.namelist()}
        read=lambda name:contents[name]
        names=set(contents)
    settings=json.loads(read('analysis_settings.json'))
    config=json.loads(read('analysis_config.json')) if 'analysis_config.json' in names else settings.get('analysis_config',settings.get('settings',{}))
    provenance=json.loads(read('provenance.json')) if 'provenance.json' in names else settings.get('provenance',{})
    markers=pd.read_csv(io.BytesIO(read('marker_gene_results.csv')))
    assert {'cluster','gene','rank'}.issubset(markers), 'Export needs cluster/gene/rank columns'
    assert not markers[['cluster','gene']].isna().any().any()
    markers['cluster']=markers.cluster.astype(str)
    assert not markers.duplicated(['cluster','gene']).any(), 'Duplicate cluster/gene keys'
    return dict(markers=markers,config=config,provenance=provenance,settings=settings,
                marker_sha256=hashlib.sha256(read('marker_gene_results.csv')).hexdigest())


def compare(left,right):
    a,b=load(left),load(right)
    differences={k:[a['config'].get(k),b['config'].get(k)]
                 for k in sorted(set(a['config'])|set(b['config'])) if a['config'].get(k)!=b['config'].get(k)}
    pa,pb=a['provenance'],b['provenance']
    comparable=bool(all(pa.get(k) and pb.get(k) for k in ['expression_source','resolved_engine','engine_version','dependency_versions'])
        and pa.get('input_sha256') and pa.get('input_sha256')==pb.get('input_sha256')
        and all(pa.get(k)==pb.get(k) for k in ['expression_source','layer','resolved_engine','engine_version','dependency_versions'])
        and pa.get('effective_config') and pb.get('effective_config'))
    effective_changes={k for k in set(pa.get('effective_config',{}))|set(pb.get('effective_config',{}))
                       if pa.get('effective_config',{}).get(k)!=pb.get('effective_config',{}).get(k)}
    threshold_only=(comparable and set(differences)=={'min_pct'} and effective_changes=={'min_pct'}
        and all(x['config'].get('min_pct')==x['provenance']['effective_config'].get('min_pct') for x in [a,b]))
    frames=[x['markers'].set_index(['cluster','gene']) for x in [a,b]]
    rows=[]
    for key in sorted(set(frames[0].index)|set(frames[1].index)):
        old=frames[0].loc[key] if key in frames[0].index else None
        new=frames[1].loc[key] if key in frames[1].index else None
        status='added' if old is None else 'removed' if new is None else 'reranked' if old['rank']!=new['rank'] else 'unchanged'
        reason='No membership/rank change.'
        evidence='direct export observation'
        statistic_changes=[]
        if old is not None and new is not None:
            for col in ['log2FC','p_value','p_adj','pct_in','pct_out','cluster_mean','other_mean']:
                if col in old and col in new and abs(old[col]-new[col])>1e-9+1e-6*abs(old[col]):statistic_changes.append(col)
        if status=='removed':
            if threshold_only and old is not None and 'pct_in' in old and pd.notna(old.pct_in) and b['config']['min_pct']>a['config']['min_pct'] and a['config']['min_pct']<=old.pct_in<b['config']['min_pct']:
                reason=f"Exported pct_in={old.pct_in:.6g} passes A min_pct={a['config']['min_pct']:g}, fails B min_pct={b['config']['min_pct']:g}."
                evidence='threshold crossing supported by same-input/source/engine/config provenance'
            else:
                reason='Insufficient export evidence to determine disappearance: could be filtering or top_n truncation; absence is not zero expression.'
                evidence='unresolved'
        elif status=='added':
            reason='Present only in B. B statistics are available, but A exports no row: cannot establish A rank/filter outcome from these exports. Top_n slot changes are possible.'
            evidence='unresolved'
        elif status=='reranked':
            reason='Shared marker moved in the exported ranked list.'
            if threshold_only and not statistic_changes:
                reason+=' Its exported statistics are unchanged; filtering can shift positions. Full pre-filter ordering is not exported.'
        if statistic_changes:reason+=' Changed exported statistics: '+', '.join(statistic_changes)+'.'
        rows.append(dict(cluster=key[0],gene=key[1],change=status,rank_A=None if old is None else int(old['rank']),
          rank_B=None if new is None else int(new['rank']),pct_in_A=None if old is None or 'pct_in' not in old else float(old.pct_in),
          pct_in_B=None if new is None or 'pct_in' not in new else float(new.pct_in),
          min_pct_A=a['config'].get('min_pct'),min_pct_B=b['config'].get('min_pct'),
          log2FC_A=None if old is None or 'log2FC' not in old else float(old.log2FC),
          log2FC_B=None if new is None or 'log2FC' not in new else float(new.log2FC),
          evidence=evidence,explanation=reason))
    counts={name:sum(r['change']==name for r in rows) for name in ['added','removed','reranked','unchanged']}
    return dict(counts=counts,config_changes=differences,provenance_comparable=comparable,
                threshold_only=threshold_only,runs=[{k:v for k,v in x.items() if k!='markers'} for x in [a,b]],markers=rows)


def render(report,path):
    def table(rows):return pd.DataFrame(rows).to_html(index=False,escape=True,float_format=lambda x:f'{x:.6g}',na_rep='not exported')
    body='<h1>CellCoPilot · Compare marker runs</h1><p>Existing exports only; no analysis is rerun to create this comparison. Missing rows are not zero expression.</p>'
    for name,run in zip(['A','B'],report['runs']):
        p=run['provenance'];c=run['config']
        body+=f'<h2>Run {name}</h2>'+table([{'matrix':p.get('expression_source','unspecified'),'engine':p.get('resolved_engine',c.get('engine','unspecified')),
             'engine_version':p.get('engine_version','unspecified'),'input_sha256':p.get('input_sha256','unspecified'),'marker_csv_sha256':run['marker_sha256'],'target_clusters':c.get('target_clusters'),
             'min_pct':c.get('min_pct'),'top_n':c.get('top_n'),'rank_by':c.get('rank_by')}])
        body+='<details><summary>Complete exported configuration/settings/provenance</summary><pre>'+html.escape(json.dumps(run,indent=2))+'</pre></details>'
    body+='<h2>Settings changed</h2>'+table([{'setting':k,'A':v[0],'B':v[1]} for k,v in report['config_changes'].items()])
    body+='<p>Comparable input/source/engine/dependency provenance: '+str(report['provenance_comparable'])+'. Confirmed sole analysis setting change is min_pct: '+str(report['threshold_only'])+'.</p>'
    body+='<h2>Marker changes</h2>'+table([report['counts']])
    if not any(report['counts'][k] for k in ['added','removed','reranked']):body+='<p>No membership or rank changes. Check exported configurations above for any non-marker differences.</p>'
    for name in ['removed','added','reranked','unchanged']:
        body+='<h3>'+name.title()+'</h3>'+table([r for r in report['markers'] if r['change']==name])
    body+='<p>min_pct and pct_in are fractions, not percentages. Log2FC is Scanpy’s approximate fold change. Means/detection use the selected matrix’s existing scale; raw is a storage location, not a preprocessing guarantee. A top_n CSV does not include failed genes or complete prior ranks.</p>'
    path.write_text('<!doctype html><html><head><meta charset="utf-8"><title>Compare marker runs</title><style>body{font:15px system-ui;margin:30px}table{border-collapse:collapse;font-size:12px}td,th{padding:7px;border:1px solid #ccc;max-width:360px;overflow-wrap:anywhere}pre{white-space:pre-wrap}h2{margin-top:30px}</style></head><body>'+body+'</body></html>')


def main():
    p=argparse.ArgumentParser(__doc__);p.add_argument('left',type=Path);p.add_argument('right',type=Path);p.add_argument('output',type=Path)
    args=p.parse_args();args.output.mkdir(parents=True,exist_ok=True)
    report=compare(args.left,args.right)
    (args.output/'comparison.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    render(report,args.output/'report.html')
    print(json.dumps(report['counts']))


if __name__=='__main__':main()
