#!/usr/bin/env python3
"""Compatibility CLI for the reusable CBIcall mtDNA report module."""

import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_DIR = REPO_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from cbicall import mtdna_browser as _browser  # noqa: E402


KEYS2REPORT = _browser.KEYS2REPORT
KEYS4HTML = _browser.KEYS4HTML
build_hash_out = _browser.build_hash_out
hash2array = _browser.hash2array
is_missing = _browser.is_missing
max_hf_any_sample = _browser.max_hf_any_sample
normalize_header = _browser.normalize_header
passes_hf_filter = _browser.passes_hf_filter
passes_maf_filter = _browser.passes_maf_filter
sort_sample_alphabetically = _browser.sort_sample_alphabetically
_html_cell = _browser._html_cell
_html_link = _browser._html_link
_is_http_url = _browser._is_http_url


def parse_args(argv=None):
    return _browser.build_json_parser().parse_args(argv)


def main(argv=None):
    return _browser.json_main(argv)


if __name__ == "__main__":
    raise SystemExit(main())
