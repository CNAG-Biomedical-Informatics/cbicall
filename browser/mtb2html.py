#!/usr/bin/env python3
"""Compatibility CLI for the reusable CBIcall mtDNA report module."""

import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_DIR = REPO_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from cbicall import mtdna_browser as _browser  # noqa: E402


BrowserError = _browser.BrowserError
build_download_buttons = _browser.build_download_buttons
load_payload = _browser.load_payload
render_html = _browser.render_html
render_report = _browser.render_report
resolve_json_path = _browser.resolve_json_path
_download_button = _browser._download_button
_json_for_script = _browser._json_for_script


def parse_args(argv=None):
    return _browser.build_html_parser().parse_args(argv)


def main(argv=None):
    return _browser.html_main(argv)


if __name__ == "__main__":
    raise SystemExit(main())
