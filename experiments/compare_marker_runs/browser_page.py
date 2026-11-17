"""Browser parity/edge-case checks for the undeployed read-only page."""
import argparse
import io
import json
from pathlib import Path
import re
import zipfile

from playwright.sync_api import sync_playwright

from compare import compare

ROOT=Path(__file__).resolve().parent


def packed(directory,changes=None,omit=()):
    output=io.BytesIO()
    with zipfile.ZipFile(output,'w') as z:
        for p in directory.iterdir():
            if p.name not in omit:z.writestr(p.name,(changes or {}).get(p.name,p.read_bytes()))
    return output.getvalue()


def main():
    parser=argparse.ArgumentParser(__doc__)
    parser.add_argument('--url',default='http://127.0.0.1:8507')
    parser.add_argument('--output',type=Path,default=Path('/tmp/compare-page-verification'))
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=True)
    a=packed(ROOT/'sample/A')
    checks=[]
    with sync_playwright() as p:
        browser=p.chromium.launch(executable_path='/usr/bin/chromium',headless=True,args=['--no-sandbox'])
        page=browser.new_page(viewport={'width':1700,'height':1100},accept_downloads=True)
        page.goto(args.url)
        inputs=page.locator('input[type=file]')
        def upload(index,content,name):
            inputs.nth(index).set_input_files({'name':name,'mimeType':'application/zip','buffer':content})
        upload(0,a,'earlier.zip')
        scenarios=[('original',ROOT/'sample/B'),('top_n',ROOT/'extended/top_n/B'),
                   ('combined',ROOT/'extended/combined/B'),('identical',ROOT/'sample/A'),
                   ('context_mismatch',ROOT/'extended/context_mismatch/B')]
        for label,later in scenarios:
            upload(1,packed(later),label+'.zip')
            page.get_by_role('heading',name='Exported shortlist changes',exact=True).wait_for(timeout=30000)
            page.wait_for_timeout(700)
            expected=compare(ROOT/'sample/A',later)
            table=page.locator('[data-testid="stTable"]').last.inner_text()
            numbers=[int(s) for s in re.findall(r'\b\d+\b',table)]
            assert numbers==list(expected['counts'].values()),(label,table,numbers)
            if label=='identical':page.get_by_text('Identical exported runs.',exact=True).wait_for()
            if label=='context_mismatch':page.get_by_text('Analysis contexts differ.',exact=False).wait_for()
            with page.expect_download() as pending:page.get_by_text('Download comparison JSON',exact=True).click()
            download=pending.value;assert download.failure() is None
            file=args.output/(label+'.json');download.save_as(file)
            assert json.loads(file.read_text())==expected
            if label=='original':
                adds=[r for r in expected['markers'] if r['change']=='added']
                assert len(adds)==7 and all(not r['explainable'] for r in adds)
            page.screenshot(path=str(args.output/(label+'.png')),full_page=True)
            checks.append({'case':label,'counts':expected['counts'],'download_equals_standalone':True})
        upload(1,b'broken','invalid.zip')
        page.get_by_text('Invalid or incomplete export bundle:',exact=False).wait_for()
        assert page.get_by_text('Download comparison JSON',exact=True).count()==0
        assert page.get_by_role('heading',name='Exported shortlist changes',exact=True).count()==0
        upload(1,packed(ROOT/'sample/B',omit=['marker_gene_results.csv']),'incomplete.zip')
        page.get_by_text('Incomplete bundle:',exact=False).wait_for()
        assert page.get_by_text('Download comparison JSON',exact=True).count()==0
        # No shared exported labels is not a valid cluster comparison.
        csv=(ROOT/'sample/B/marker_gene_results.csv').read_text().replace('B cells','Other cluster')
        upload(1,packed(ROOT/'sample/B',{'marker_gene_results.csv':csv}),'no_shared.zip')
        page.get_by_text('These exported shortlists have no shared clusters.',exact=False).wait_for()
        assert page.get_by_role('combobox',name='Shared cluster',exact=True).count()==0
        provenance=json.loads((ROOT/'sample/B/provenance.json').read_text());provenance.pop('input_sha256')
        upload(1,packed(ROOT/'sample/B',{'provenance.json':json.dumps(provenance)}),'unknown_identity.zip')
        page.get_by_text('Dataset identity cannot be verified:',exact=False).wait_for()
        with page.expect_download() as pending:page.get_by_text('Download comparison JSON',exact=True).click()
        file=args.output/'unknown.json';pending.value.save_as(file)
        assert not any(r['explainable'] for r in json.loads(file.read_text())['markers'])
        # Changing the earlier upload also clears the previous successful comparison.
        upload(0,b'broken','invalid_earlier.zip')
        page.get_by_text('Invalid or incomplete export bundle:',exact=False).wait_for()
        assert page.get_by_text('Download comparison JSON',exact=True).count()==0
        browser.close()
    result={'scenarios':checks,'invalid_and_incomplete_feedback':True,'no_shared_clusters_feedback':True,
            'missing_identity_warning_and_suppression':True,'either_upload_clears_stale_results':True}
    (args.output/'browser-verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
