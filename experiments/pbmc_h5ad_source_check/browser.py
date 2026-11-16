"""Browser regression: upload, negative-source notice, raw Scanpy B cells, exports, source clearing."""
from pathlib import Path
import argparse
from playwright.sync_api import sync_playwright
parser=argparse.ArgumentParser(__doc__)
parser.add_argument('directory',type=Path)
parser.add_argument('--url',default='http://127.0.0.1:8506/')
args=parser.parse_args()
root=args.directory.resolve()
with sync_playwright() as p:
 browser=p.chromium.launch(executable_path='/usr/bin/chromium',headless=True,args=['--no-sandbox'])
 page=browser.new_page(viewport={'width':1450,'height':1100},accept_downloads=True)
 page.goto(args.url)
 page.locator('input[type=file]').first.set_input_files(str(root/'pbmc3k.h5ad'))
 page.get_by_text('Analysis settings',exact=True).wait_for(timeout=60000)
 page.wait_for_timeout(2000)
 page.get_by_role('combobox',name='Cluster column',exact=True).click()
 page.get_by_role('option',name='louvain',exact=True).click()
 page.wait_for_timeout(1000)
 assert page.get_by_text('The selected expression matrix contains negative values',exact=False).count()==1
 page.screenshot(path=str(root/'after-X-warning.png'),full_page=True)
 page.get_by_role('button',name='Find marker genes',exact=True).click()
 page.wait_for_timeout(2000)
 (root/'before-X.txt').write_text(page.locator('body').inner_text())
 print('X errors:',page.locator('[data-testid="stAlert"]').all_inner_texts(),flush=True)
 page.get_by_role('combobox',name='Analysis engine',exact=True).click()
 page.get_by_role('option',name='Scanpy',exact=True).click()
 page.wait_for_timeout(1000)
 page.get_by_role('combobox',name='Expression matrix',exact=True).click()
 print('sources',page.get_by_role('option').all_inner_texts(),flush=True)
 page.get_by_role('option',name='AnnData.raw',exact=True).click()
 page.wait_for_timeout(1000)
 # Remove seven tags by their accessible remove buttons; keep B cells.
 for label in ['CD14+ Monocytes','CD4 T cells','CD8 T cells','Dendritic cells','FCGR3A+ Monocytes','Megakaryocytes','NK cells']:
  page.get_by_role('button',name='Remove '+label,exact=True).click()
 page.get_by_role('button',name='Find marker genes',exact=True).click()
 page.get_by_text('Marker-gene results',exact=True).wait_for(timeout=120000)
 (root/'raw-result.txt').write_text(page.locator('body').inner_text())
 print(page.locator('body').inner_text()[-4500:],flush=True)
 for label,filename in [('Download results as CSV','app-markers.csv'),('Download analysis settings as JSON','app-settings.json'),('Download analysis bundle as ZIP','app-bundle.zip'),('Download expression summary as CSV','app-summary.csv')]:
  with page.expect_download(timeout=30000) as pending:page.get_by_text(label,exact=True).click()
  dl=pending.value;assert dl.failure() is None;dl.save_as(root/filename)
 page.screenshot(path=str(root/'raw-result.png'),full_page=True)
 page.get_by_role('combobox',name='Expression matrix',exact=True).click()
 page.get_by_role('option',name='AnnData.X (default)',exact=True).click()
 page.get_by_text('Marker-gene results',exact=True).wait_for(state='hidden',timeout=30000)
 assert page.get_by_text('Download results as CSV',exact=True).count()==0
 print('Source switch cleared results and downloads',flush=True)
 browser.close()
