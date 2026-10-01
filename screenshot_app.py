"""Capture a screenshot of a running Marker Gene Finder app."""

from __future__ import annotations

import argparse
from pathlib import Path

from playwright.sync_api import sync_playwright


def build_parser() -> argparse.ArgumentParser:
    """Create command-line options for screenshot capture."""
    parser = argparse.ArgumentParser(
        description="Capture a full-page screenshot of the Streamlit app."
    )
    parser.add_argument(
        "--url",
        default="http://localhost:8501",
        help="URL of the running app (default: %(default)s).",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("marker-gene-finder.png"),
        help="Screenshot path (default: %(default)s).",
    )
    return parser


def capture_screenshot(url: str, output: Path) -> None:
    """Open the app in Chromium and save a full-page PNG screenshot."""
    output.parent.mkdir(parents=True, exist_ok=True)
    with sync_playwright() as playwright:
        browser = playwright.chromium.launch()
        page = browser.new_page(viewport={"width": 1440, "height": 1000})
        page.goto(url, wait_until="domcontentloaded")
        page.get_by_text("CellCoPilot Marker Gene Finder", exact=False).first.wait_for(
            state="visible",
            timeout=30_000,
        )
        page.screenshot(path=str(output), full_page=True)
        browser.close()


def main() -> None:
    """Parse arguments and capture the screenshot."""
    args = build_parser().parse_args()
    capture_screenshot(args.url, args.output)
    print(f"Saved screenshot to {args.output}")


if __name__ == "__main__":
    main()
