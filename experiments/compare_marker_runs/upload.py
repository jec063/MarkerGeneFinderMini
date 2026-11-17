"""Bounded upload adapter for the existing export-only comparison engine."""
import io
import json
import re
import tempfile
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd

from marker_gene_finder import AnalysisConfig
from compare import compare, load

MAX_UPLOAD=20*1024*1024
MAX_EXPANDED=50*1024*1024
FILES={'marker_gene_results.csv','analysis_settings.json','analysis_config.json','provenance.json'}


class BundleError(ValueError):
    pass


def validate_bundle(content, directory):
    """Copy only bounded flat known members; reuse the existing bundle loader."""
    try:
        if len(content)>MAX_UPLOAD:raise BundleError('Bundle exceeds the 20 MiB upload limit.')
        with zipfile.ZipFile(io.BytesIO(content)) as archive:
            members=archive.infolist()
            if len(members)>32 or sum(m.file_size for m in members)>MAX_EXPANDED:
                raise BundleError('Bundle exceeds the member count or 50 MiB expanded limit.')
            names=[m.filename for m in members]
            if len(names)!=len(set(names)):raise BundleError('Bundle contains duplicate filenames.')
            for m in members:
                if '/' in m.filename or '\\' in m.filename or m.filename in {'.','..'}:
                    raise BundleError('Expected a flat MarkerGeneFinderMini bundle, without nested paths.')
                if m.flag_bits & 1 or (m.external_attr >> 16)&0o170000 == 0o120000:
                    raise BundleError('Encrypted files and symbolic links are not supported.')
            if not {'marker_gene_results.csv','analysis_settings.json'}.issubset(names):
                raise BundleError('Incomplete bundle: marker_gene_results.csv and analysis_settings.json are required.')
            for name in FILES&set(names):(directory/name).write_bytes(archive.read(name))
        settings=json.loads((directory/'analysis_settings.json').read_text())
        if not isinstance(settings,dict) or not isinstance(settings.get('settings'),dict):
            raise BundleError('Invalid analysis settings: expected the current versioned settings document.')
        if settings.get('schema_version')!=2:raise BundleError('Unsupported settings schema; expected version 2.')
        if not (directory/'analysis_config.json').exists() and not settings.get('analysis_config'):
            raise BundleError('Incomplete bundle: a replayable analysis_config.json or embedded analysis_config is required.')
        data=load(directory)  # existing CSV/config/provenance selection and key validation
        if not isinstance(data['config'],dict) or not isinstance(data['provenance'],dict):
            raise BundleError('Invalid configuration or provenance object.')
        AnalysisConfig.from_dict(data['config'])
        digest=data['provenance'].get('input_sha256')
        if digest is not None and (not isinstance(digest,str) or not re.fullmatch('[0-9a-fA-F]{64}',digest)):
            raise BundleError('Invalid recorded dataset SHA-256; expected 64 hexadecimal characters or no hash.')
        if data['provenance'].get('effective_config') is not None and not isinstance(data['provenance']['effective_config'],dict):
            raise BundleError('Invalid effective configuration in provenance.')
        frame=data['markers']
        if len(frame)>100000:raise BundleError('Marker table exceeds 100,000 exported rows.')
        if not pd.api.types.is_numeric_dtype(frame['rank']) or not np.isfinite(frame['rank']).all() or not ((frame['rank']>=1)&(frame['rank']%1==0)).all():
            raise BundleError('Marker ranks must be positive finite integers.')
        for column in ['pct_in','pct_out','log2FC','p_value','p_adj','cluster_mean','other_mean']:
            if column in frame:
                if not pd.api.types.is_numeric_dtype(frame[column]) or not np.isfinite(frame[column]).all():
                    raise BundleError(f'Marker {column} must contain finite numbers.')
                if column in ['pct_in','pct_out','p_value','p_adj'] and not frame[column].between(0,1).all():
                    raise BundleError(f'Marker {column} must be between 0 and 1.')
        return data
    except BundleError:raise
    except (ValueError,TypeError,KeyError,AssertionError,OSError,zipfile.BadZipFile,RuntimeError,NotImplementedError) as error:
        raise BundleError(f'Invalid or incomplete export bundle: {error}') from error


def compare_uploads(earlier,later):
    with tempfile.TemporaryDirectory() as tmp:
        root=Path(tmp);a=root/'earlier';b=root/'later';a.mkdir();b.mkdir()
        left=validate_bundle(earlier,a);right=validate_bundle(later,b)
        shared=sorted(set(left['markers'].cluster)&set(right['markers'].cluster))
        try:
            report=compare(a,b)  # exact existing machine-readable diff, no UI-specific inference
            json.dumps(report,allow_nan=False)
        except (ValueError,TypeError,KeyError,AssertionError,AttributeError) as error:
            raise BundleError(f'Invalid comparison metadata: {error}') from error
    return report,shared


def dataset_notice(report):
    a,b=[r['provenance'].get('input_sha256') for r in report['runs']]
    if not all(isinstance(v,str) and re.fullmatch('[0-9a-fA-F]{64}',v) for v in [a,b]):
        return 'Dataset identity cannot be verified: one or both bundles lack a valid recorded input SHA-256.'
    if a!=b:return 'Dataset identity differs: the bundles record different input SHA-256 hashes. No causal explanations are assigned.'
    return 'Dataset identity: matching recorded input SHA-256 hashes. Original datasets are not uploaded or independently rehashed here.'


def reason_label(row):
    if not row['explainable'] and row['change']!='unchanged':return 'Cause unavailable from these exports'
    if row['evidence'].startswith('threshold crossing'):return 'Below the new expression threshold'
    if row['evidence'].startswith('top_n truncation'):return 'Outside the new top N'
    if row['evidence'].startswith('rank shift'):return 'Rank changed after filtering earlier entries'
    return 'No shortlist change'


def cluster_overview(report,shared):
    """Summarize exported entries only; never interpret them as expression change."""
    rows=[]
    for cluster in shared:
        markers=[r for r in report['markers'] if r['cluster']==cluster]
        counts={k:sum(r['change']==k for r in markers) for k in ['added','removed','reranked']}
        rows.append({'Cluster':cluster,'Added':counts['added'],'Removed':counts['removed'],
            'Reranked':counts['reranked'],'Membership changes':counts['added']+counts['removed'],
            'Changed entries':sum(counts.values())})
    return sorted(rows,key=lambda r:(-r['Changed entries'],r['Cluster']))
