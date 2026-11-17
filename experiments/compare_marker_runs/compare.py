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
    context_keys=['input_sha256','expression_source','layer','resolved_engine','engine_version','dependency_versions']
    context_changes={k:[pa.get(k),pb.get(k)] for k in context_keys if pa.get(k)!=pb.get(k)}
    for key in ['cluster_column','target_clusters','excluded_genes','excluded_prefixes','min_cluster_size','use_raw','layer']:
        if a['config'].get(key)!=b['config'].get(key):
            context_changes[key]=[a['config'].get(key),b['config'].get(key)]
    required=['input_sha256','expression_source','resolved_engine','engine_version','dependency_versions','effective_config']
    metadata_complete=all(pa.get(k) and pb.get(k) for k in required)
    comparable=bool(metadata_complete and not context_changes)
    effective_changes={k for k in set(pa.get('effective_config',{}))|set(pb.get('effective_config',{}))
                       if pa.get('effective_config',{}).get(k)!=pb.get('effective_config',{}).get(k)}
    controlled=(comparable and set(differences)<= {'min_pct','top_n'}
        and effective_changes==set(differences)
        and all(x['config'].get(k)==x['provenance']['effective_config'].get(k)
                for x in [a,b] for k in ['min_pct','top_n']))
    threshold_only=controlled and set(differences)=={'min_pct'}
    # Tightening a filter preserves ordering among survivors. An exported prefix
    # can therefore establish new ranks for old rows; never invent unseen rows.
    known_ranks={}
    monotone=controlled and b['config']['min_pct']>=a['config']['min_pct']
    if monotone:
        for cluster,old in a['markers'].groupby('cluster'):
            old=old.sort_values('rank')
            new=b['markers'][b['markers'].cluster==cluster].sort_values('rank')
            if 'pct_in' not in old or old.pct_in.isna().any():continue
            if old['rank'].tolist()!=list(range(1,len(old)+1)):continue
            if new['rank'].tolist()!=list(range(1,len(new)+1)):continue
            if len(old)>a['config']['top_n'] or len(new)>b['config']['top_n']:continue
            survivors=old[old.pct_in>=b['config']['min_pct']]
            expected=survivors.head(b['config']['top_n'])
            if expected.gene.tolist()!=new.gene.head(len(expected)).tolist():continue
            same_stats=True
            for gene in set(old.gene)&set(new.gene):
                oa=old[old.gene==gene].iloc[0];nb=new[new.gene==gene].iloc[0]
                for col in ['log2FC','p_value','p_adj','pct_in','pct_out','cluster_mean','other_mean']:
                    if col not in oa or col not in nb or not pd.notna(oa[col]) or not pd.notna(nb[col]):
                        same_stats=False;break
                    if abs(oa[col]-nb[col])>1e-9+1e-6*abs(oa[col]):same_stats=False
            if same_stats:
                for rank,gene in enumerate(survivors.gene,1):known_ranks[(cluster,gene)]=rank
    frames=[x['markers'].set_index(['cluster','gene']) for x in [a,b]]
    rows=[]
    for key in sorted(set(frames[0].index)|set(frames[1].index)):
        old=frames[0].loc[key] if key in frames[0].index else None
        new=frames[1].loc[key] if key in frames[1].index else None
        status='added' if old is None else 'removed' if new is None else 'reranked' if old['rank']!=new['rank'] else 'unchanged'
        reason='No membership/rank change.'
        evidence='direct export observation'
        explainable=False
        statistic_changes=[]
        if old is not None and new is not None:
            for col in ['log2FC','p_value','p_adj','pct_in','pct_out','cluster_mean','other_mean']:
                if col in old and col in new and abs(old[col]-new[col])>1e-9+1e-6*abs(old[col]):statistic_changes.append(col)
        if status=='removed':
            if controlled and old is not None and 'pct_in' in old and pd.notna(old.pct_in) and b['config']['min_pct']>a['config']['min_pct'] and a['config']['min_pct']<=old.pct_in<b['config']['min_pct']:
                reason=f"Exported pct_in={old.pct_in:.6g} passes A min_pct={a['config']['min_pct']:g}, fails B min_pct={b['config']['min_pct']:g}."
                evidence='threshold crossing supported by same-input/source/engine/config provenance'
                explainable=True
            elif key in known_ranks and known_ranks[key]>b['config']['top_n']:
                reason=f"Passes B min_pct={b['config']['min_pct']:g}; rank after filtering the exported A prefix is {known_ranks[key]}, outside B top_n={b['config']['top_n']}."
                evidence='top_n truncation established from exported prefix'
                explainable=True
            else:
                reason='Insufficient export evidence to determine disappearance: could be filtering or top_n truncation; absence is not zero expression.'
                evidence='unresolved'
        elif status=='added':
            reason='Present only in B. B statistics are available, but A exports no row: cannot establish A rank/filter outcome from these exports. Top_n slot changes are possible.'
            evidence='unresolved'
        elif status=='reranked':
            reason='Shared marker moved in the exported ranked list.'
            if key in known_ranks and known_ranks[key]==int(new['rank']) and not statistic_changes:
                reason+=f" Unchanged statistics; removing earlier A rows below B min_pct gives exported rank {known_ranks[key]}."
                evidence='rank shift established from exported prefix'
                explainable=True
        if status!='unchanged' and not comparable:
            reason+=' Analysis contexts differ or required metadata is missing; no causal explanation is assigned.'
            evidence='unresolved';explainable=False
        if statistic_changes:reason+=' Changed exported statistics: '+', '.join(statistic_changes)+'.'
        rows.append(dict(cluster=key[0],gene=key[1],change=status,rank_A=None if old is None else int(old['rank']),
          rank_B=None if new is None else int(new['rank']),pct_in_A=None if old is None or 'pct_in' not in old else float(old.pct_in),
          pct_in_B=None if new is None or 'pct_in' not in new else float(new.pct_in),
          min_pct_A=a['config'].get('min_pct'),min_pct_B=b['config'].get('min_pct'),
          log2FC_A=None if old is None or 'log2FC' not in old else float(old.log2FC),
          log2FC_B=None if new is None or 'log2FC' not in new else float(new.log2FC),
          explainable=explainable,rank_from_exported_prefix_B=known_ranks.get(key),evidence=evidence,explanation=reason))
    counts={name:sum(r['change']==name for r in rows) for name in ['added','removed','reranked','unchanged']}
    return dict(counts=counts,explainable_changes={k:sum(r['change']==k and r['explainable'] for r in rows) for k in ['added','removed','reranked']},config_changes=differences,provenance_comparable=comparable,
                context_status='different' if context_changes else 'comparable' if comparable else 'insufficient metadata',context_changes=context_changes,
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
    body+='<h2>Analysis context</h2><p>'+html.escape(report['context_status'])+'</p>'+table([{'field':k,'A':v[0],'B':v[1]} for k,v in report['context_changes'].items()])
    body+='<h2>Explainable shortlist changes</h2>'+table([report['explainable_changes']])
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
