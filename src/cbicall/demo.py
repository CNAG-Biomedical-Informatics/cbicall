"""Generate a self-contained CBIcall demonstration from packaged artifacts."""

from __future__ import annotations

import hashlib
import json
import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List

from .html_reports import write_run_report_html
from .mtdna_browser import (
    BrowserError,
    generate_browser_report,
    load_prioritized_variants,
)


ASSET_DIR = Path(__file__).with_name("demo_assets")
PROJECT_ID = "CNAG99901P"


class DemoError(RuntimeError):
    """Raised when the packaged demonstration cannot be generated."""


@dataclass(frozen=True)
class DemoResult:
    output_dir: Path
    wes_report: Path
    wes_html: Path
    mtdna_html: Path
    mtdna_variants: int


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _load_wes_report(report_path: Path) -> Dict[str, object]:
    try:
        payload = json.loads(report_path.read_text(encoding="utf-8"))
    except OSError as exc:
        raise DemoError("Cannot read packaged WES report: {}".format(exc)) from exc
    except json.JSONDecodeError as exc:
        raise DemoError("Packaged WES report is invalid: {}".format(exc)) from exc
    if not isinstance(payload, dict):
        raise DemoError("Packaged WES report must be a JSON object")
    return payload


def _verify_wes_vcf(wes_dir: Path, report: Dict[str, object]) -> None:
    outputs = report.get("outputs")
    if not isinstance(outputs, dict):
        raise DemoError("Packaged WES report does not contain output metadata")
    reports = outputs.get("vcf_hash_reports")
    if not isinstance(reports, list) or not reports or not isinstance(reports[0], dict):
        raise DemoError("Packaged WES report does not contain a VCF fingerprint")
    expected = reports[0].get("raw_sha256")
    vcf_path = wes_dir / "02_varcall" / "CNAG99901P.hc.QC.vcf.gz"
    if not vcf_path.is_file():
        raise DemoError("Packaged WES VCF is missing: {}".format(vcf_path))
    observed = _sha256_file(vcf_path)
    if expected != observed:
        raise DemoError("Packaged WES VCF checksum does not match its audit report")


def _write_filtered_json(input_path: Path, output_path: Path) -> int:
    variants: List[Dict[str, object]] = load_prioritized_variants(input_path)
    output_path.write_text(
        json.dumps(variants, ensure_ascii=False, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return len(variants)


def _write_readme(result: DemoResult) -> None:
    (result.output_dir / "README.txt").write_text(
        "CBIcall demo\n"
        "============\n\n"
        "This directory is a resource-free tour generated from compact, "
        "precomputed outputs of the packaged CNAG99901P integration fixture. "
        "CBIcall did not execute BWA, GATK, or MToolBox during this demo. The "
        "fixture is test data and is not an analytical benchmark.\n\n"
        "Open these reports\n"
        "------------------\n"
        "- wes/run-report.html: native WES execution and audit summary\n"
        "- mtdna/02_browser/CNAG99901P.html: interactive mtDNA variant browser "
        "({} filtered variants)\n\n"
        "The underlying final WES VCF and QC summaries are under wes/. The "
        "mtDNA report, haplogroup file, VCF, and filtered JSON are under "
        "mtdna/01_mtoolbox/.\n\n"
        "To execute the corresponding workflows instead of viewing precomputed "
        "outputs, install the external resource bundle, set CBICALL_DATA, and run:\n\n"
        "  cbicall test --wes-bash --mit-bash -t 1\n"
        .format(result.mtdna_variants),
        encoding="utf-8",
    )


def run_demo(output_dir: Path) -> DemoResult:
    """Create WES and mtDNA demonstration reports in a new directory."""
    destination = output_dir.expanduser().resolve()
    if destination.exists():
        raise DemoError("Demo output already exists: {}".format(destination))

    wes_assets = ASSET_DIR / "wes"
    mtdna_assets = ASSET_DIR / "mtdna"
    if not wes_assets.is_dir() or not mtdna_assets.is_dir():
        raise DemoError("Packaged demo assets are missing; reinstall CBIcall")

    try:
        destination.mkdir(parents=True)
        wes_dir = destination / "wes"
        mtdna_dir = destination / "mtdna"
        shutil.copytree(wes_assets, wes_dir)
        shutil.copytree(mtdna_assets, mtdna_dir)

        wes_report = wes_dir / "run-report.json"
        wes_payload = _load_wes_report(wes_report)
        _verify_wes_vcf(wes_dir, wes_payload)
        wes_html = write_run_report_html(
            wes_report,
            wes_payload,
            html_path=wes_dir / "run-report.html",
        )

        mtoolbox_dir = mtdna_dir / "01_mtoolbox"
        browser_dir = mtdna_dir / "02_browser"
        prioritized = mtoolbox_dir / "mit_prioritized_variants.txt"
        filtered_json = mtoolbox_dir / "mit.filtered.json"
        mtdna_variants = _write_filtered_json(prioritized, filtered_json)
        mtdna_html = browser_dir / "CNAG99901P.html"
        summary = generate_browser_report(
            prioritized,
            mtdna_html,
            project_id=PROJECT_ID,
            job_id="precomputed-demo",
            browser_json_path=browser_dir / "mit.json",
        )
        if summary.get("variants") != mtdna_variants:
            raise DemoError("mtDNA browser and filtered JSON variant counts differ")

        result = DemoResult(
            output_dir=destination,
            wes_report=wes_report,
            wes_html=wes_html,
            mtdna_html=mtdna_html,
            mtdna_variants=mtdna_variants,
        )
        _write_readme(result)
        return result
    except DemoError:
        shutil.rmtree(destination, ignore_errors=True)
        raise
    except (BrowserError, OSError, TypeError, ValueError) as exc:
        shutil.rmtree(destination, ignore_errors=True)
        raise DemoError("Cannot generate demo: {}".format(exc)) from exc


__all__ = ["ASSET_DIR", "DemoError", "DemoResult", "run_demo"]
