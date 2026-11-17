"""Read-only experimental demo. Deploy separately from the stable app."""
from pathlib import Path
import hashlib
import json
import sys

import pandas as pd
import streamlit as st

ROOT=Path(__file__).resolve().parent
# Community Cloud and local launches need not share a working directory.
sys.path.insert(0,str(ROOT.parents[1]))
sys.path.insert(0,str(ROOT))
from upload import BundleError,compare_uploads,dataset_notice,reason_label,cluster_overview

st.set_page_config(page_title='CellCoPilot · Compare exported runs',layout='wide')
st.title('CellCoPilot · Compare exported runs')
st.write('Upload two MarkerGeneFinderMini export bundles to inspect changes in their exported marker shortlists and the evidence explaining them, or try the sample comparison below.')
st.caption('Experimental · read-only · no expression data upload, reanalysis, or annotation verdict')


def clear():
    for key in ['comparison','shared_clusters','upload_pair','shared_cluster']:
        st.session_state.pop(key,None)


def use_sample():
    clear()
    st.session_state['input_mode']='Muraro sample'


st.button('Try sample comparison',on_click=use_sample)
mode=st.radio('Comparison inputs',['Upload my bundles','Muraro sample'],
              key='input_mode',on_change=clear,horizontal=True)
if mode=='Muraro sample':
    st.info('Muraro pancreas sample: alpha, beta and duct versus the rest. Both runs use the same log1p_cp10k layer and Scanpy; min_pct changes from 0.10 to 0.25.')
    try:
        a,b=[(ROOT/'pancreas/results'/name).read_bytes() for name in ['earlier.zip','later.zip']]
    except OSError:
        clear()
        st.error('The bundled Muraro sample could not be loaded. Switch to Upload my bundles to compare your own exports.')
        st.stop()
else:
    earlier=st.file_uploader('Earlier run',type=['zip'],key='earlier',on_change=clear)
    later=st.file_uploader('Later run',type=['zip'],key='later',on_change=clear)
    if earlier is None or later is None:
        st.info('Upload two existing MarkerGeneFinderMini export bundles to compare.')
        st.stop()
    a,b=earlier.getvalue(),later.getvalue()
pair=(mode,hashlib.sha256(a).hexdigest(),hashlib.sha256(b).hexdigest())
if st.session_state.get('upload_pair')!=pair:
    clear()
    try:
        report,shared=compare_uploads(a,b)
    except BundleError as error:
        st.error(str(error));st.stop()
    st.session_state['comparison']=report
    st.session_state['shared_clusters']=shared
    # Deleting a widget key alone can restore its previous frontend value.
    # Initialize the new comparison explicitly before rendering the selector.
    if shared:st.session_state['shared_cluster']=shared[0]
    st.session_state['upload_pair']=pair
report=st.session_state['comparison']

# Changes precede shortlist counts; avoid A/B jargon in the visible headings.
st.header('Changed settings')
if report['config_changes']:
    st.table(pd.DataFrame([{'Setting':k,'Earlier run':json.dumps(v[0]),'Later run':json.dumps(v[1])} for k,v in report['config_changes'].items()]))
else:st.info('No exported configuration changes.')
st.subheader('Selected matrix and engine')
st.table(pd.DataFrame([{'Run':label,'Matrix':r['provenance'].get('expression_source','Unspecified'),
    'Layer':r['provenance'].get('layer') or 'None',
    'Engine':r['provenance'].get('resolved_engine','Unspecified'),'Engine version':r['provenance'].get('engine_version','Unspecified')}
    for label,r in zip(['Earlier run','Later run'],report['runs'])]))
notice=dataset_notice(report)
if 'cannot' in notice or 'differs' in notice:st.warning(notice)
else:st.caption(notice)
st.subheader('Analysis context')
if report['context_changes']:
    st.warning('Analysis contexts differ. Shortlist differences are shown without causal explanations.')
    st.table(pd.DataFrame([{'Field':k,'Earlier run':json.dumps(v[0]),'Later run':json.dumps(v[1])} for k,v in report['context_changes'].items()]))
elif not report['provenance_comparable']:
    st.warning('Required context metadata is incomplete. No causal explanations are assigned.')
else:st.caption('Recorded analysis contexts match. Only explanations supported by the comparison engine are shown.')
with st.expander('All exported settings and provenance'):
    st.json(report['runs'])

shared=st.session_state['shared_clusters']
if not shared:
    st.warning('These exported shortlists have no shared clusters. No cluster comparison can be shown.')
else:
    if len(shared)>1:
        overview=cluster_overview(report,shared)
        st.header('Changes across shared clusters')
        st.table(pd.DataFrame(overview))
        top=overview[0]['Changed entries']
        leaders=', '.join(r['Cluster'] for r in overview if r['Changed entries']==top)
        st.caption(f'Most changed by exported entries: {leaders} ({top}).')
        st.caption('Changed entries = added + removed + reranked; membership changes = added + removed. Neither measures biological change.')
    cluster=st.selectbox('Shared cluster',shared,key='shared_cluster')
    rows=[r for r in report['markers'] if r['cluster']==cluster]
    counts={k:sum(r['change']==k for r in rows) for k in ['added','removed','reranked','unchanged']}
    st.header('Exported shortlist changes')
    st.table(pd.DataFrame([{'Added':counts['added'],'Removed':counts['removed'],'Reranked':counts['reranked'],'Unchanged':counts['unchanged']}]))
    st.caption('Counts refer to exported shortlist membership and position, not expression absence or biological identity.')
    if not any(counts[k] for k in ['added','removed','reranked']):
        st.success('No exported shortlist changes for this cluster.')
    if a==b or report['runs'][0]==report['runs'][1]:st.info('Identical exported runs.')
    for kind,label in [('added','Added to exported shortlist'),('removed','Removed from exported shortlist'),('reranked','Reranked within exported shortlist')]:
        st.subheader(label)
        selected=[r for r in rows if r['change']==kind]
        if not selected:st.caption('None.');continue
        if kind=='added':
            st.caption('The Earlier run export has no row for these genes, so their earlier ranks and filter outcomes are unavailable. Added does not mean newly expressed.')
        elif kind=='removed':
            st.caption('The Later run export has no row for these genes. Removed means removed from its exported shortlist, not absence of expression.')
        display=[]
        for row in selected:
            display.append({'Gene':row['gene'],'Evidence':reason_label(row),
                'Earlier rank':str(row['rank_A']) if row['rank_A'] is not None else 'Unavailable',
                'Later rank':str(row['rank_B']) if row['rank_B'] is not None else 'Unavailable',
                'Earlier detection':row['pct_in_A'],'Later detection':row['pct_in_B'],
                'Earlier min_pct':row['min_pct_A'],'Later min_pct':row['min_pct_B'],
                'Later top N':report['runs'][1]['config'].get('top_n'),
                'Rank from earlier prefix after filtering':row['rank_from_exported_prefix_B'],
                'Explanation':row['explanation']})
        st.dataframe(pd.DataFrame(display),hide_index=True,width='stretch')
    st.caption('Detection/min_pct are fractions (0.10 = 10%). Missing ranks are unavailable, never zero. raw is a storage location; preprocessing depends on the file.')
st.download_button('Download comparison JSON',json.dumps(report,indent=2,allow_nan=False)+'\n',
                   file_name='comparison.json',mime='application/json',on_click='ignore')
st.caption('The download contains the complete existing diff for all clusters; the selector only filters this view.')
