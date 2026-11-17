"""Verify sample/upload mode parity against the retained bundles in Chromium."""
import argparse
import json
from pathlib import Path

from playwright.sync_api import sync_playwright

from compare import compare

ROOT=Path(__file__).resolve().parent


def main():
    parser=argparse.ArgumentParser(__doc__)
    parser.add_argument('--url',default='http://127.0.0.1:8510')
    parser.add_argument('--output',type=Path,default=Path('/tmp/comparison-demo-browser'))
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=True)
    earlier=ROOT/'pancreas/results/earlier.zip';later=ROOT/'pancreas/results/later.zip'
    expected=compare(earlier,later)
    checks=[]
    with sync_playwright() as p:
        browser=p.chromium.launch(executable_path='/usr/bin/chromium',headless=True,args=['--no-sandbox'])
        page=browser.new_page(viewport={'width':1450,'height':1050},accept_downloads=True)
        page.goto(args.url)
        page.get_by_role('button',name='Try sample comparison',exact=True).wait_for(timeout=120000)
        assert page.locator('input[type=file]').count()==2
        assert 'inspect changes in their exported marker shortlists' in page.locator('body').inner_text()
        page.get_by_role('button',name='Try sample comparison',exact=True).click()
        page.get_by_role('heading',name='Changes across shared clusters',exact=True).wait_for()
        page.wait_for_timeout(800)
        assert page.get_by_role('radio',name='Muraro sample',exact=True).is_checked()
        assert page.locator('input[type=file]').count()==0
        settings=page.locator('[data-testid=stTable]').first.inner_text()
        assert 'min_pct' in settings and '0.1' in settings and '0.25' in settings
        page.get_by_text('Most changed by exported entries: beta (28).',exact=True).wait_for()
        def select(cluster):
            box=page.get_by_role('combobox',name='Shared cluster',exact=True)
            box.scroll_into_view_if_needed();page.wait_for_timeout(300);box.click();page.wait_for_timeout(500)
            page.get_by_role('option',name=cluster,exact=True).click();page.wait_for_timeout(1000)
        def download(name):
            with page.expect_download() as pending:page.get_by_text('Download comparison JSON',exact=True).click()
            path=args.output/(name+'.json');pending.value.save_as(path)
            assert json.loads(path.read_text())==expected
            checks.append({'mode':name,'download_matches_comparator':True})
        select('beta')
        assert [str(x) for x in [9,9,10,1]]==page.locator('[data-testid=stTable]').last.locator('td').all_text_contents()
        page.get_by_text('The Earlier run export has no row for these genes',exact=False).wait_for()
        download('sample-beta')
        page.get_by_role('heading',name='Added to exported shortlist',exact=True).scroll_into_view_if_needed()
        page.screenshot(path=str(args.output/'sample-beta.png'))
        select('duct')
        assert page.locator('[data-testid=stTable]').last.locator('td').all_text_contents()==['5','5','15','0']
        download('sample-duct')
        upload_label=page.get_by_text('Upload my bundles',exact=True)
        upload_label.scroll_into_view_if_needed();page.wait_for_timeout(300);upload_label.click()
        page.wait_for_timeout(800)
        assert page.get_by_role('radio',name='Upload my bundles',exact=True).is_checked()
        assert page.locator('input[type=file]').count()==2
        assert not page.get_by_role('heading',name='Exported shortlist changes',exact=True).count()
        assert not page.get_by_text('Download comparison JSON',exact=True).count()
        for i,path in enumerate([earlier,later]):page.locator('input[type=file]').nth(i).set_input_files(str(path))
        page.get_by_role('heading',name='Exported shortlist changes',exact=True).wait_for();page.wait_for_timeout(800)
        select('beta');download('uploaded-beta')
        # An invalid replacement must remove the uploaded result, not show the sample.
        page.locator('input[type=file]').nth(1).set_input_files({'name':'invalid.zip','mimeType':'application/zip','buffer':b'broken'})
        page.get_by_text('Invalid or incomplete export bundle:',exact=False).wait_for()
        assert not page.get_by_text('Download comparison JSON',exact=True).count()
        page.get_by_role('button',name='Try sample comparison',exact=True).click()
        page.get_by_role('heading',name='Exported shortlist changes',exact=True).wait_for();page.wait_for_timeout(800)
        assert page.get_by_role('combobox',name='Shared cluster',exact=True).input_value()=='alpha'
        assert not page.get_by_text('Invalid or incomplete export bundle:',exact=False).count()
        download('sample-after-invalid-upload')
        browser.close()
    result={'checks':checks,'setting_and_cluster_overview_correct':True,
            'sample_upload_switch_clears_visible_results':True,'invalid_upload_does_not_reuse_sample':True,
            'sample_action_recovers_and_resets_cluster':True}
    (args.output/'verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
