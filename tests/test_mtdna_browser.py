import importlib.util
import json
from pathlib import Path
import shutil
import subprocess
import sys

import pytest

from cbicall import mtdna_browser as browser


REPO_ROOT = Path(__file__).resolve().parents[1]


def load_wrapper(name):
    path = REPO_ROOT / "browser" / (name + ".py")
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def write_prioritized_report(path, rows=None):
    header = [
        "Variant Allele",
        "Sample",
        "Locus",
        "Aa Change",
        "REF",
        "ALT",
        "GT",
        "DP",
        "HF",
        "Disease Score",
        "Mitomap Associated Disease(s)",
        "ClinVar",
        "OMIM link",
        "dbSNP ID",
        "AC/AN 1000 Genomes",
    ]
    if rows is None:
        rows = [
            [
                "123A",
                "sample-DNA_MIT",
                "MT-ND1",
                "p.Val1Ala",
                "A",
                "G",
                "0/1",
                "100",
                "0.75",
                "0.91",
                "Example disease",
                "RCV0001",
                "https://omim.org/entry/1",
                "rs123",
                "0.001",
            ]
        ]
    lines = ["\t".join(header)] + ["\t".join(row) for row in rows]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def test_filter_helpers_cover_single_and_cohort_values():
    assert browser.is_missing(None)
    assert browser.is_missing("N/A")
    assert not browser.is_missing("0.4")
    assert browser.normalize_header([" Variant Allele "]) == ["Variant_Allele"]
    assert browser.sort_sample_alphabetically("") == ""
    assert browser.sort_sample_alphabetically("B,A") == "A,B"

    columns = {"HF": 0, "AC/AN_1000_Genomes": 1}
    assert browser.passes_maf_filter(["0.4", "0.5"], columns, 0)
    assert browser.passes_maf_filter(["0.4", ""], columns, 0.01)
    assert browser.passes_maf_filter(["0.4", "unknown"], columns, 0.01)
    assert not browser.passes_maf_filter(["0.4", "0.01"], columns, 0.01)
    assert browser.passes_maf_filter(["0.4"], {"HF": 0}, 0.01)

    assert browser.max_hf_any_sample("NA") is None
    assert browser.max_hf_any_sample("0.2,0.7") == 0.7
    assert browser.max_hf_any_sample("A:0.4,NA|B:bad,0.8") == 0.8
    assert browser.max_hf_any_sample("A:bad") is None
    assert browser.passes_hf_filter(["NA"], {"HF": 0}, 0.3)
    assert not browser.passes_hf_filter(["NA"], {"HF": 0}, 0.3, True)
    assert not browser.passes_hf_filter(["0.3"], {"HF": 0}, 0.3)
    assert browser.passes_hf_filter(["0.4"], {"HF": 0}, 0.3)
    assert browser.passes_hf_filter(["0.1"], {}, 0.3)
    assert browser.passes_hf_filter(["0.1"], {"HF": 0}, 0)


def test_parser_preserves_existing_filtering_and_json_shape(tmp_path, capsys):
    report = write_prioritized_report(
        tmp_path / "mit_prioritized_variants.txt",
        rows=[
            ["123A", "B-DNA_MIT,A-DNA_MIT", "MT-ND1", "p.Val1Ala", "A", "G", "0/1", "100", "A:0.75|B:0.5", "0.91", "Disease", "RCV1", "https://omim.org/1", "rs1", "0.001"],
            ["124A", "A-DNA_MIT", "MT-ND1", "syn", "A", "T", "0/1", "80", "0.8", "0.2", "", "", "", "", "0.001"],
            ["125A", "A-DNA_MIT", "MT-ND2", "p.X", "A", "C", "0/1", "80", "0.2", "0.2", "", "", "", "", "0.001"],
            ["126A", "A-DNA_MIT", "MT-ND3", "p.X", "A", "C", "0/1", "80", "0.8", "0.2", "", "", "", "", "0.1"],
        ],
    )

    parsed = browser.load_prioritized_variants(report)
    assert list(parsed) == ["A-DNA_MIT,B-DNA_MIT_MT-ND1_123A"]
    record = next(iter(parsed.values()))
    assert record["Sample"] == "A-DNA_MIT,B-DNA_MIT"
    assert record["HF"] == "A:0.75|B:0.5"

    assert browser.json_main(["-i", str(report), "-f", "json"]) == 0
    emitted = json.loads(capsys.readouterr().out)
    assert emitted == parsed

    assert browser.json_main(["-i", str(report), "-f", "json4html"]) == 0
    browser_payload = json.loads(capsys.readouterr().out)
    assert browser_payload["data"][0]["locus"] == "MT-ND1"

    assert browser.json_main(["-i", str(report), "-f", "tsv"]) == 0
    assert capsys.readouterr().out.startswith("Variant_Allele\tSample\tLocus")

    assert browser.json_main(["-i", str(report), "-f", "hash"]) == 0
    assert "MT-ND1" in capsys.readouterr().out


def test_build_hash_out_skips_short_and_incomplete_rows(capsys):
    header = browser.normalize_header(["Variant Allele", "Sample", "Locus", "HF"])
    result = browser.build_hash_out(
        [["short"], ["1A", "sample", "MT-ND1", "0.8"]],
        header,
        0.3,
        0,
        debug=True,
    )
    assert list(result) == ["sample_MT-ND1_1A"]
    assert "Skipping short line 2" in capsys.readouterr().err
    assert browser.build_hash_out([["1A"]], ["Variant_Allele"], 0, 0) == {}


def test_input_and_payload_validation_errors(tmp_path):
    empty = tmp_path / "empty.txt"
    empty.write_text("", encoding="utf-8")
    with pytest.raises(browser.BrowserError, match="Empty mtDNA input"):
        browser.load_prioritized_variants(empty)
    with pytest.raises(browser.BrowserError, match="Cannot read mtDNA input"):
        browser.load_prioritized_variants(tmp_path / "missing.txt")
    with pytest.raises(browser.BrowserError, match="top-level 'data' array"):
        browser.build_report_payload(
            {"data": None},
            project_id="p",
            job_id="j",
            source_name="bad.json",
        )

    payload = browser.build_report_payload(
        {
            "data": [
                {
                    "sample": "sample-DNA_MIT",
                    "locus": "MT-ND1",
                    "heteroplasmy": "0.5",
                    "diseaseScore": "0.9",
                    "dbsnp": "rs1",
                }
            ]
        },
        project_id="p",
        job_id="j",
        source_name="minimal.json",
    )
    assert payload["summary"]["samples"] == 1
    assert payload["summary"]["heteroplasmic"] == 1
    assert payload["summary"]["highDisease"] == 1


def test_structured_payload_has_mtdna_metrics_and_fields(tmp_path):
    report = write_prioritized_report(tmp_path / "report.txt")
    parsed = browser.load_prioritized_variants(report)
    rows = browser.hash_to_browser_rows(parsed)
    payload = browser.build_report_payload(
        {"data": rows},
        project_id="project",
        job_id="job",
        source_name=report.name,
    )

    assert payload["summary"] == {
        "variants": 1,
        "samples": 1,
        "evidence": 1,
        "highDisease": 1,
        "heteroplasmic": 1,
    }
    assert payload["samples"] == ["sample"]
    assert payload["loci"] == ["MT-ND1"]
    assert payload["rows"][0]["variantAllele"] == "123A"
    assert payload["rows"][0]["maxHeteroplasmy"] == 0.75


def test_legacy_browser_array_is_still_accepted():
    payload = browser.build_report_payload(
        {"data": [["S&lt;1&gt;", '<a href="x">MT-ND1</a>', "123A", "A", "G", "p.X", "1", "50", "0.5"]]},
        project_id="project",
        job_id="job",
        source_name="legacy.json",
    )
    assert payload["rows"][0]["sample"] == "S<1>"
    assert payload["rows"][0]["locus"] == "MT-ND1"
    assert payload["summary"]["heteroplasmic"] == 1


def test_legacy_html_helpers_escape_values_and_links():
    wrapper = load_wrapper("mtb2json")
    assert wrapper._html_cell("Sample", "S-DNA_MIT<script>") == "S&lt;script&gt;"
    assert wrapper._html_cell("HF", "0.4|<x>") == "0.4,<br />&lt;x&gt;"
    locus = wrapper._html_cell("Locus", "MT-ND1")
    assert 'rel="noopener noreferrer"' in locus
    assert 'href="https://ghr.nlm.nih.gov/gene/MT-ND1#conditions"' in locus
    assert "<a " not in wrapper._html_cell("OMIM_link", "javascript:alert(1)")
    assert wrapper._is_http_url("https://omim.org/1")
    assert "_plus_" in wrapper._html_cell("Mitomap_Associated_Disease(s)", "A+B")
    assert "ncbi.nlm.nih.gov/snp/rs1" in wrapper._html_cell("dbSNP_ID", "rs1")
    assert "omim.org/1" in wrapper._html_cell("OMIM_link", "https://omim.org/1")
    assert browser.hash2array({"key": {"Sample": "sample-DNA_MIT"}})[0][0] == "sample"


def test_report_embeds_tabulator_data_and_escapes_script_end(tmp_path):
    report = write_prioritized_report(tmp_path / "report.txt")
    parsed = browser.load_prioritized_variants(report)
    rows = browser.hash_to_browser_rows(parsed)
    rows[0]["mitomapDisease"] = "</script><script>alert(1)</script>"

    rendered = browser.render_html(
        "P<1>",
        'J"2',
        'mit "x".json',
        {"data": rows},
    )

    assert "P&lt;1&gt;" in rendered
    assert "J&quot;2" in rendered
    assert "mit &quot;x&quot;.json" in rendered
    assert "CBIcall mtDNA" in rendered
    assert 'id="detail-panel"' in rendered
    assert 'id="queue-tabs"' in rendered
    assert 'id="report-data"' in rendered
    assert "new Tabulator" in rendered
    assert "paginationSizeSelector: [25, 50, 100, 250]" in rendered
    assert "<\\/script><script>alert(1)<\\/script>" in rendered
    assert "Filtered JSON" in rendered
    assert "mit.filtered.json" in rendered
    assert "mit.raw.json" not in rendered
    assert "cdn.jsdelivr.net" not in rendered
    assert "BFF Tools Browser" not in rendered


def test_generate_report_writes_standalone_html_and_browser_json(tmp_path):
    report = write_prioritized_report(tmp_path / "report.txt")
    html_path = tmp_path / "browser" / "run.html"
    json_path = tmp_path / "browser" / "mit.json"

    summary = browser.generate_browser_report(
        report,
        html_path,
        project_id="project",
        job_id="run",
        browser_json_path=json_path,
    )

    assert summary["variants"] == 1
    assert html_path.is_file()
    assert json.loads(json_path.read_text(encoding="utf-8"))["data"][0]["locus"] == "MT-ND1"
    assert not (html_path.parent / "run.html.tmp").exists()


def test_report_asset_and_write_errors_are_wrapped(tmp_path, monkeypatch):
    monkeypatch.setattr(browser, "TEMPLATE_FILE", tmp_path / "missing-template.html")
    with pytest.raises(browser.BrowserError, match="Cannot read mtDNA browser assets"):
        browser.render_report({"projectId": "p", "jobId": "j", "source": "s"})

    def fail_write(self, content, encoding=None):
        raise OSError("read-only")

    monkeypatch.setattr(Path, "write_text", fail_write)
    with pytest.raises(browser.BrowserError, match="Cannot write"):
        browser._write_text_atomic(tmp_path / "out.html", "content")


def test_payload_loading_resolves_colocated_json_and_reports_errors(tmp_path):
    out = tmp_path / "browser" / "report.html"
    out.parent.mkdir()
    payload = {"data": [{"sample": "sample"}]}
    (out.parent / "mit.json").write_text(json.dumps(payload), encoding="utf-8")

    assert browser.load_payload("mit.json", str(out)) == payload
    with pytest.raises(browser.BrowserError, match="Cannot read"):
        browser.load_payload("missing.json", str(out))
    bad = out.parent / "bad.json"
    bad.write_text("not json", encoding="utf-8")
    with pytest.raises(browser.BrowserError, match="Cannot parse"):
        browser.load_payload(str(bad), str(out))
    wrong = out.parent / "wrong.json"
    wrong.write_text('{"rows": []}', encoding="utf-8")
    with pytest.raises(browser.BrowserError, match="top-level 'data' array"):
        browser.load_payload(str(wrong), str(out))
    with pytest.raises(browser.BrowserError, match="objects or legacy arrays"):
        browser.build_report_payload(
            {"data": ["bad"]},
            project_id="p",
            job_id="j",
            source_name="bad.json",
        )


def test_cli_wrappers_generate_json_and_html(tmp_path):
    report = write_prioritized_report(tmp_path / "report.txt")
    json_path = tmp_path / "mit.json"
    html_path = tmp_path / "report.html"

    with json_path.open("w", encoding="utf-8") as handle:
        subprocess.run(
            [sys.executable, "browser/mtb2json.py", "-i", str(report), "-f", "json4html"],
            cwd=REPO_ROOT,
            stdout=handle,
            text=True,
            check=True,
        )
    subprocess.run(
        [sys.executable, "browser/mtb2html.py", "--id", "project", "--job-id", "job", "--json", str(json_path), "--out", str(html_path)],
        cwd=REPO_ROOT,
        text=True,
        check=True,
    )
    assert "CBIcall mtDNA" in html_path.read_text(encoding="utf-8")


def test_in_process_cli_error_and_html_paths(tmp_path):
    with pytest.raises(SystemExit, match="Cannot read mtDNA input"):
        browser.json_main(["-i", str(tmp_path / "missing.txt"), "-f", "json"])

    payload_path = tmp_path / "mit.json"
    output_path = tmp_path / "report.html"
    payload_path.write_text('{"data": []}', encoding="utf-8")
    assert browser.html_main([
        "--id", "project",
        "--job-id", "job",
        "--json", str(payload_path),
        "--out", str(output_path),
    ]) == 0
    assert output_path.is_file()

    with pytest.raises(SystemExit, match="Cannot read"):
        browser.html_main([
            "--id", "project",
            "--job-id", "job",
            "--json", str(tmp_path / "missing.json"),
            "--out", str(output_path),
        ])


def test_report_executes_in_headless_chromium_when_available(tmp_path):
    chromium = shutil.which("chromium") or shutil.which("chromium-browser")
    if chromium is None:
        pytest.skip("Chromium is not installed")

    report = write_prioritized_report(tmp_path / "report.txt")
    parsed = browser.load_prioritized_variants(report)
    smoke_dir = REPO_ROOT / ".pytest_cache" / "mtdna_browser"
    smoke_dir.mkdir(parents=True, exist_ok=True)
    html_path = smoke_dir / "mtdna.html"
    try:
        html_path.write_text(
            browser.render_html(
                "project",
                "job",
                "mit.json",
                {"data": browser.hash_to_browser_rows(parsed)},
            ),
            encoding="utf-8",
        )

        proc = subprocess.run(
            [
                chromium,
                "--headless",
                "--disable-gpu",
                "--disable-dev-shm-usage",
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
    assert 'id="metric-variants">1</strong>' in proc.stdout
    assert 'class="tabulator"' in proc.stdout
    assert "1 of 1 variant" in proc.stdout
