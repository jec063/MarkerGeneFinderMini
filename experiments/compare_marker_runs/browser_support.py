"""Browser helpers for local pages and Community Cloud's framed app document."""
import time

from playwright.sync_api import Error, expect


def open_app(page,url,timeout=120000):
    """Find the rendered comparison app without depending on Cloud's iframe URL."""
    page.goto(url,wait_until='domcontentloaded',timeout=60000)
    deadline=time.monotonic()+timeout/1000
    while time.monotonic()<deadline:
        for frame in page.frames:
            try:
                if frame.get_by_role('button',name='Try sample comparison',exact=True).is_visible():
                    return frame
            except Error:
                if not frame.is_detached():raise
        page.wait_for_timeout(250)
    raise TimeoutError(f'Comparison app did not become ready at {url}')


def wait_counts(app,counts):
    """Wait for the requested rerun, rather than reading a stale/incomplete table."""
    cells=app.locator('[data-testid="stTable"]').last.locator('td')
    expect(cells).to_have_text([str(x) for x in counts],timeout=60000)
