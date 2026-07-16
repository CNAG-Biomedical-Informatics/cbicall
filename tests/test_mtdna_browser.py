import json
from pathlib import Path
import shutil
import subprocess
import sys

import pytest

from cbicall import mtdna_browser as browser


REPO_ROOT = Path(__file__).resolve().parents[1]

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

    report_payload = browser.build_report_payload(
        parsed,
        project_id="project",
        job_id="job",
        source_name=report.name,
    )
    assert report_payload["samples"] == ["A", "B"]

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
    payload = browser.build_report_payload(
        {
            "sample_MT-ND1_1A": {
                "Sample": "sample-DNA_MIT",
                "Locus": "MT-ND1",
                "Variant_Allele": "1A",
                "HF": "0.5",
                "Disease_Score": "0.9",
                "dbSNP_ID": "rs1",
            }
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
    payload = browser.build_report_payload(
        parsed,
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


def test_report_embeds_tabulator_data_and_escapes_script_end(tmp_path):
    report = write_prioritized_report(tmp_path / "report.txt")
    parsed = browser.load_prioritized_variants(report)
    next(iter(parsed.values()))["Mitomap_Associated_Disease(s)"] = "</script><script>alert(1)</script>"

    rendered = browser.render_html(
        "P<1>",
        'J"2',
        'mit.filtered "x".json',
        parsed,
    )

    assert "P&lt;1&gt;" in rendered
    assert "J&quot;2" in rendered
    assert "mit.filtered &quot;x&quot;.json" in rendered
    assert "CBIcall mtDNA" in rendered
    assert 'id="detail-panel"' in rendered
    assert 'id="view-tabs"' in rendered
    assert 'id="queue-tabs"' not in rendered
    assert 'id="view"' in rendered
    assert "All variants" in rendered
    assert "External evidence" in rendered
    assert "Heteroplasmic" in rendered
    assert "Review queue" not in rendered
    assert 'id="report-data"' in rendered
    assert "new Tabulator" in rendered
    assert "pagination: true" in rendered
    assert "paginationSize: 25" in rendered
    assert "paginationSizeSelector: [10, 25, 50, 100, 250]" in rendered
    assert 'id="table-pagination"' in rendered
    assert 'id="table-scroll-left"' in rendered
    assert 'id="table-scroll-right"' in rendered
    assert 'id="table-horizontal-scroll"' in rendered
    assert 'id="table-horizontal-scroll-content"' in rendered
    assert 'class="column-scroll-label">Columns</span>' in rendered
    assert "scrollTableColumns" in rendered
    assert "holder.scrollLeft = Math.max" in rendered
    assert "horizontalRail.scrollLeft = holder.scrollLeft" in rendered
    assert "holder.scrollLeft = horizontalRail.scrollLeft" in rendered
    assert 'addEventListener("scroll", handleTableHorizontalScroll' in rendered
    assert 'paginationElement: document.getElementById("table-pagination")' in rendered
    assert 'table.on("pageLoaded", updatePageCount)' in rendered
    assert "row._samples.indexOf(controls.sample.value)" in rendered
    assert "height: availableTableHeight()" in rendered
    assert "table.setHeight(availableTableHeight())" in rendered
    assert 'class="brand-mark"' not in rendered
    assert "<\\/script><script>alert(1)<\\/script>" in rendered
    assert "Filtered JSON" in rendered
    assert "mit.filtered.json" in rendered
    assert "mit.raw.json" not in rendered
    assert "cdn.jsdelivr.net" not in rendered
    assert "BFF Tools Browser" not in rendered


def test_generate_report_reads_canonical_json_and_writes_standalone_html(tmp_path):
    report = write_prioritized_report(tmp_path / "report.txt")
    filtered_json = tmp_path / "mtoolbox" / "mit.filtered.json"
    html_path = tmp_path / "browser" / "run.html"
    filtered_records = browser.write_filtered_json(report, filtered_json)

    summary = browser.generate_browser_report(
        filtered_json,
        html_path,
        project_id="project",
        job_id="run",
    )

    assert summary["variants"] == 1
    assert browser.load_filtered_variants(filtered_json) == filtered_records
    assert html_path.is_file()
    assert 'href="../01_mtoolbox/mit.filtered.json"' in html_path.read_text(encoding="utf-8")
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


def test_filtered_json_loading_reports_errors(tmp_path):
    with pytest.raises(browser.BrowserError, match="Cannot read"):
        browser.load_filtered_variants(tmp_path / "missing.json")
    bad = tmp_path / "bad.json"
    bad.write_text("not json", encoding="utf-8")
    with pytest.raises(browser.BrowserError, match="Cannot parse"):
        browser.load_filtered_variants(bad)
    wrong_shape = tmp_path / "wrong-shape.json"
    wrong_shape.write_text("[]", encoding="utf-8")
    with pytest.raises(browser.BrowserError, match="top-level object"):
        browser.load_filtered_variants(wrong_shape)
    wrong_record = tmp_path / "wrong-record.json"
    wrong_record.write_text('{"variant": []}', encoding="utf-8")
    with pytest.raises(browser.BrowserError, match="not an object"):
        browser.load_filtered_variants(wrong_record)


def test_cli_wrappers_generate_json_and_html(tmp_path):
    report = write_prioritized_report(tmp_path / "report.txt")
    json_path = tmp_path / "mit.filtered.json"
    html_path = tmp_path / "report.html"

    with json_path.open("w", encoding="utf-8") as handle:
        subprocess.run(
            [sys.executable, "browser/mtb2json.py", "-i", str(report), "-f", "json"],
            cwd=REPO_ROOT,
            stdout=handle,
            text=True,
            check=True,
        )
    subprocess.run(
        [sys.executable, "browser/mtb2html.py", "--id", "project", "--job-id", "job", "--filtered-json", str(json_path), "--out", str(html_path)],
        cwd=REPO_ROOT,
        text=True,
        check=True,
    )
    assert "CBIcall mtDNA" in html_path.read_text(encoding="utf-8")


def test_in_process_cli_error_and_html_paths(tmp_path):
    with pytest.raises(SystemExit, match="Cannot read mtDNA input"):
        browser.json_main(["-i", str(tmp_path / "missing.txt"), "-f", "json"])

    payload_path = tmp_path / "mit.filtered.json"
    output_path = tmp_path / "report.html"
    payload_path.write_text("{}", encoding="utf-8")
    assert browser.html_main([
        "--id", "project",
        "--job-id", "job",
        "--filtered-json", str(payload_path),
        "--out", str(output_path),
    ]) == 0
    assert output_path.is_file()

    with pytest.raises(SystemExit, match="Cannot read"):
        browser.html_main([
            "--id", "project",
            "--job-id", "job",
            "--filtered-json", str(tmp_path / "missing.json"),
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
                "mit.filtered.json",
                parsed,
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


def test_mtdna_workflows_use_only_canonical_filtered_json():
    for relative_path in (
        "workflows/bash/gatk-3.5/mit_single.sh",
        "workflows/bash/gatk-3.5/mit_cohort.sh",
    ):
        script = (REPO_ROOT / relative_path).read_text(encoding="utf-8")
        assert script.count('"$PYBINDIR/mtb2json.py"') == 1
        assert '-f json > "$mit_filtered_json"' in script
        assert '--filtered-json "$mit_filtered_json"' in script
