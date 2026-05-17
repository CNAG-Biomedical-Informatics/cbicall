#!/usr/bin/env python3
"""
mtb2html.py - Generate a standalone CBIcall mtDNA HTML viewer for MToolBox output

Usage:
  mtb2json.py -i mit_prioritized_variants.txt -f json4html > mtDNA.json
  mtb2html.py --id peter --job-id 1234 --json mtDNA.json --out mtdna.html
"""

import argparse
import html as html_lib
import json
from pathlib import Path
import sys


COLUMNS = [
    "Sample",
    "Locus",
    "Variant allele",
    "Ref",
    "Alt",
    "AA change",
    "GT",
    "Depth",
    "Heteroplasmy",
    "tRNA annotation",
    "Disease score",
    "RNA predictions",
    "Mitomap disease",
    "Mitomap homoplasmy",
    "Mitomap heteroplasmy",
    "ClinVar",
    "OMIM",
    "dbSNP",
    "Mamit-tRNA",
    "AC/AN 1000G",
    "1000G homoplasmy",
    "1000G heteroplasmy",
]

DEFAULT_HIDDEN_COLUMNS = [9, 11, 13, 14, 18, 20, 21]
EVIDENCE_COLUMNS = [12, 15, 16, 17]
TEMPLATE_FILE = Path(__file__).resolve().parent / "templates" / "mtdna_report.html"


def parse_args():
    p = argparse.ArgumentParser(
        description="Generate a standalone CBIcall mtDNA HTML viewer"
    )
    p.add_argument("--id", required=True, help="Project ID to display in the HTML header")
    p.add_argument("--job-id", required=True, help="Job ID to display in the HTML header")
    p.add_argument(
        "--json",
        default="mtDNA.json",
        help="JSON file from mtb2json.py -f json4html",
    )
    p.add_argument("--out", "-o", default="mtdna.html", help="Output HTML file name")
    return p.parse_args()


def _download_button(href, label):
    return (
        '<a class="action-button" href="'
        + html_lib.escape(href, quote=True)
        + '" download>'
        + '<svg class="icon" aria-hidden="true"><use href="#icon-download"></use></svg>'
        + html_lib.escape(label)
        + "</a>"
    )


def build_download_buttons(base_dir="../01_mtoolbox"):
    files = [
        (base_dir + "/mit_prioritized_variants.txt", "Report"),
        (base_dir + "/mt_classification_best_results.csv", "Haplogroup"),
        (base_dir + "/VCF_file.vcf", "VCF"),
        (base_dir + "/mit.raw.json", "Raw JSON"),
    ]
    return "\n          ".join(_download_button(href, label) for href, label in files)


def resolve_json_path(json_file, html_out):
    json_path = Path(json_file)
    if json_path.is_absolute() or json_path.exists():
        return json_path

    out_dir = Path(html_out).resolve().parent if html_out else Path.cwd()
    colocated = out_dir / json_file
    if colocated.exists():
        return colocated

    return json_path


def load_payload(json_file, html_out=None):
    json_path = resolve_json_path(json_file, html_out)
    try:
        with open(json_path, "r", encoding="utf-8") as handle:
            payload = json.load(handle)
    except OSError as e:
        sys.exit("Cannot read {}: {}".format(json_path, e))
    except json.JSONDecodeError as e:
        sys.exit("Cannot parse {}: {}".format(json_path, e))

    if not isinstance(payload, dict) or not isinstance(payload.get("data"), list):
        sys.exit("{} must contain a top-level 'data' array".format(json_path))
    return payload


def _json_for_script(value):
    return json.dumps(value, ensure_ascii=False).replace("</", "<\\/")


def load_template():
    try:
        return TEMPLATE_FILE.read_text(encoding="utf-8")
    except OSError as e:
        sys.exit("Cannot read {}: {}".format(TEMPLATE_FILE, e))


def render_html(
    project_id,
    job_id,
    json_file="mtDNA.json",
    data_payload=None,
    disease_threshold=0.80,
):
    if data_payload is None:
        data_payload = {"data": []}

    html = load_template()
    html = html.replace("__PROJECT_ID__", html_lib.escape(project_id))
    html = html.replace("__JOB_ID__", html_lib.escape(job_id))
    html = html.replace("__SOURCE_JSON__", html_lib.escape(json_file))
    html = html.replace("__DOWNLOAD_BUTTONS__", build_download_buttons())
    html = html.replace("__DATA_JSON__", _json_for_script(data_payload))
    html = html.replace("__COLUMNS_JSON__", _json_for_script(COLUMNS))
    html = html.replace("__DEFAULT_HIDDEN_JSON__", json.dumps(DEFAULT_HIDDEN_COLUMNS))
    html = html.replace("__EVIDENCE_COLUMNS_JSON__", json.dumps(EVIDENCE_COLUMNS))
    html = html.replace("__DISEASE_THRESHOLD__", "{:.2f}".format(disease_threshold))
    return html


def main():
    args = parse_args()
    payload = load_payload(args.json, args.out)
    try:
        with open(args.out, "w", encoding="utf-8") as out:
            out.write(render_html(args.id, args.job_id, args.json, payload))
    except OSError as e:
        sys.exit("Cannot write {}: {}".format(args.out, e))
    return 0


if __name__ == "__main__":
    sys.exit(main())
