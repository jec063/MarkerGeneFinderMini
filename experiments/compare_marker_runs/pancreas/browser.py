"""Check real multi-cluster exports through the experimental page."""
import argparse
import json
from pathlib import Path
import re
import sys

from playwright.sync_api import sync_playwright

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from compare import compare

ROOT=Path(__file__).resolve().parent
parser=argparse.ArgumentParser(__doc__)
parser.add_argument('--url',default='http://127.0.0.1:8508')
parser.add_argument('--output',type=Path,default=Path('/tmp/compare-pancreas-browser'))
parser.add_argument('--expect-overview',action='store_true')
args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=True)
exports=ROOT/'results'
expected=compare(exports/'earlier.zip',exports/'later.zip')
checks=[]
with sync_playwright() as p:
    browser=p.chromium.launch(executable_path='/usr/bin/chromium',headless=True,args=['--no-sandbox'])
    overview=False
    for cluster in ['alpha','beta','duct']:
        page=browser.new_page(viewport={'width':1750,'height':1100},accept_downloads=True)
        page.goto(args.url)
        for i,label in enumerate(['earlier','later']):
            page.locator('input[type=file]').nth(i).set_input_files(str(exports/(label+'.zip')))
        page.get_by_role('heading',name='Exported shortlist changes',exact=True).wait_for(timeout=30000)
        page.wait_for_timeout(1200)
        settings=page.locator('[data-testid="stTable"]').first.inner_text()
        assert 'min_pct' in settings and '0.1' in settings and '0.25' in settings
        overview=bool(page.get_by_role('heading',name='Changes across shared clusters',exact=True).count())
        if args.expect_overview:
            assert overview
            page.get_by_text('Most changed by exported entries: beta',exact=False).wait_for()
        if cluster != 'alpha':
            selector=page.get_by_role('combobox',name='Shared cluster',exact=True)
            selector.scroll_into_view_if_needed()
            page.wait_for_timeout(300)
            selector.click()
            page.wait_for_timeout(500)
            page.get_by_role('option',name=cluster,exact=True).click()
            page.wait_for_timeout(1000)
        page.wait_for_timeout(700)
        rows=[r for r in expected['markers'] if r['cluster']==cluster]
        counts={k:sum(r['change']==k for r in rows) for k in ['added','removed','reranked','unchanged']}
        table=page.locator('[data-testid="stTable"]').last.inner_text()
        assert [int(s) for s in re.findall(r'\b\d+\b',table)]==list(counts.values()), (cluster,table,counts)
        with page.expect_download() as pending:page.get_by_text('Download comparison JSON',exact=True).click()
        download=pending.value;assert download.failure() is None
        path=args.output/(cluster+'.json');download.save_as(path)
        assert json.loads(path.read_text())==expected
        page.screenshot(path=str(args.output/(cluster+'.png')),full_page=True)
        checks.append({'cluster':cluster,'counts':counts,'matches_standalone':True})
        page.close()
    browser.close()
result={'changed_setting':{'min_pct':[.1,.25]},'overview_present':overview,'clusters':checks,
    'same_full_json_download_for_each_selected_cluster':True}
(args.output/'verification.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
