import importlib.util
import json
from pathlib import Path
import shutil
import subprocess
import sys


REPO_ROOT = Path(__file__).resolve().parents[1]


def load_browser_module(name):
    path = REPO_ROOT / "browser" / f"{name}.py"
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_mtb2json_html_cells_escape_values_and_links():
    mtb2json = load_browser_module("mtb2json")

    assert mtb2json._html_cell("Sample", "S-DNA_MIT<script>") == "S&lt;script&gt;"
    assert mtb2json._html_cell("HF", "0.4|<x>") == "0.4,<br />&lt;x&gt;"

    locus = mtb2json._html_cell("Locus", "MT-ND1")
    assert 'rel="noopener noreferrer"' in locus
    assert 'href="https://ghr.nlm.nih.gov/gene/MT-ND1#conditions"' in locus

    assert "<a " not in mtb2json._html_cell("OMIM_link", "javascript:alert(1)")


def test_mtb2html_builds_standalone_report_without_legacy_assets():
    mtb2html = load_browser_module("mtb2html")

    html = mtb2html.render_html(
        "P<1>",
        'J"2',
        'mit "x".json',
        {"data": [["S&lt;1&gt;", "MT-ND1", "123A", "A", "G"]]},
    )

    assert "P&lt;1&gt;" in html
    assert "J&quot;2" in html
    assert "mit &quot;x&quot;.json" in html
    assert 'id="mtdna-data"' in html
    assert "assets/" not in html
    assert "jquery" not in html.lower()
    assert 'return /[",\\n]/.test(text)' in html
    assert 'lines.join("\\n") + "\\n"' in html


def test_mtb2html_load_payload_resolves_json_next_to_output(tmp_path):
    mtb2html = load_browser_module("mtb2html")

    out = tmp_path / "browser" / "report.html"
    out.parent.mkdir()
    payload = {"data": [["sample"]]}
    (out.parent / "mit.json").write_text(json.dumps(payload), encoding="utf-8")

    assert mtb2html.load_payload("mit.json", out) == payload


def test_mtb2html_executes_in_headless_chromium_when_available():
    chromium = shutil.which("chromium") or shutil.which("chromium-browser")
    if chromium is None:
        return

    mtb2html = load_browser_module("mtb2html")
    smoke_dir = REPO_ROOT / ".pytest_cache" / "mtdna_browser"
    smoke_dir.mkdir(parents=True, exist_ok=True)
    html_path = smoke_dir / "mtdna.html"
    try:
        html_path.write_text(
            mtb2html.render_html(
                "project",
                "job",
                "mit.json",
                {"data": [["sample", "MT-ND1", "123A", "A", "G"]]},
            ),
            encoding="utf-8",
        )

        proc = subprocess.run(
            [
                chromium,
                "--headless",
                "--disable-gpu",
                "--no-sandbox",
                "--dump-dom",
                html_path.as_uri(),
            ],
            capture_output=True,
            text=True,
            check=True,
        )
    finally:
        html_path.unlink(missing_ok=True)

    assert 'id="metric-variants">1</div>' in proc.stdout
    assert "Showing 1-1 of 1 variants" in proc.stdout


def test_mtb2json_converts_reference_mtdna_report_for_browser():
    report = (
        REPO_ROOT
        / "examples/input/CNAG999_exome/CNAG99901P_ex"
        / "ref_cbicall_bash_mit_single_rsrs_gatk-3.5_649547582283533"
        / "01_mtoolbox/mit_prioritized_variants.txt"
    )

    proc = subprocess.run(
        [sys.executable, "browser/mtb2json.py", "-i", str(report), "-f", "json4html"],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
        check=True,
    )

    payload = json.loads(proc.stdout)
    assert "data" in payload
    assert isinstance(payload["data"], list)
    assert payload["data"]
