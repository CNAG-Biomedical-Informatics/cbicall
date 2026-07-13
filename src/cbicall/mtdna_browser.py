from __future__ import annotations

import argparse
import csv
import html
import json
import re
import sys
from html.parser import HTMLParser
from pathlib import Path
from pprint import pprint
from typing import Any, Dict, Iterable, List, Optional, Sequence
from urllib.parse import quote, urlparse


ASSET_DIR = Path(__file__).with_name("mtdna_browser_assets")
TEMPLATE_FILE = ASSET_DIR / "report.html"
VENDOR_DIR = ASSET_DIR / "vendor"
TABULATOR_CSS = VENDOR_DIR / "tabulator-6.5.0.min.css"
TABULATOR_JS = VENDOR_DIR / "tabulator-6.5.0.min.js"

KEYS2REPORT = [
    "Variant_Allele",
    "Sample",
    "Locus",
    "Nt_Variability",
    "Codon_Position",
    "Aa_Change",
    "Aa_Variability",
    "REF",
    "ALT",
    "GT",
    "DP",
    "HF",
    "tRNA_Annotation",
    "Disease_Score",
    "RNA_predictions",
    "Mitomap_Associated_Disease(s)",
    "Mitomap_Homoplasmy",
    "Mitomap_Heteroplasmy",
    "Somatic_Mutations",
    "SM_Homoplasmy",
    "SM_Heteroplasmy",
    "ClinVar",
    "OMIM_link",
    "dbSNP_ID",
    "Mamit-tRNA_link",
    "PhastCons20Way",
    "PhyloP20Way",
    "AC/AN_1000_Genomes",
    "1000_Genomes_Homoplasmy",
    "1000_Genomes_Heteroplasmy",
]

KEYS4HTML = [
    "Sample",
    "Locus",
    "Variant_Allele",
    "REF",
    "ALT",
    "Aa_Change",
    "GT",
    "DP",
    "HF",
    "tRNA_Annotation",
    "Disease_Score",
    "RNA_predictions",
    "Mitomap_Associated_Disease(s)",
    "Mitomap_Homoplasmy",
    "Mitomap_Heteroplasmy",
    "ClinVar",
    "OMIM_link",
    "dbSNP_ID",
    "Mamit-tRNA_link",
    "AC/AN_1000_Genomes",
    "1000_Genomes_Homoplasmy",
    "1000_Genomes_Heteroplasmy",
]

BROWSER_FIELDS = [
    ("sample", "Sample", "Sample"),
    ("locus", "Locus", "Locus"),
    ("variantAllele", "Variant allele", "Variant_Allele"),
    ("ref", "Ref", "REF"),
    ("alt", "Alt", "ALT"),
    ("aaChange", "AA change", "Aa_Change"),
    ("gt", "GT", "GT"),
    ("depth", "Depth", "DP"),
    ("heteroplasmy", "Heteroplasmy", "HF"),
    ("trnaAnnotation", "tRNA annotation", "tRNA_Annotation"),
    ("diseaseScore", "Disease score", "Disease_Score"),
    ("rnaPredictions", "RNA predictions", "RNA_predictions"),
    ("mitomapDisease", "MITOMAP disease", "Mitomap_Associated_Disease(s)"),
    ("mitomapHomoplasmy", "MITOMAP homoplasmy", "Mitomap_Homoplasmy"),
    ("mitomapHeteroplasmy", "MITOMAP heteroplasmy", "Mitomap_Heteroplasmy"),
    ("clinvar", "ClinVar", "ClinVar"),
    ("omim", "OMIM", "OMIM_link"),
    ("dbsnp", "dbSNP", "dbSNP_ID"),
    ("mamitTrna", "Mamit-tRNA", "Mamit-tRNA_link"),
    ("populationFrequency", "AC/AN 1000G", "AC/AN_1000_Genomes"),
    ("populationHomoplasmy", "1000G homoplasmy", "1000_Genomes_Homoplasmy"),
    ("populationHeteroplasmy", "1000G heteroplasmy", "1000_Genomes_Heteroplasmy"),
    ("ntVariability", "Nucleotide variability", "Nt_Variability"),
    ("codonPosition", "Codon position", "Codon_Position"),
    ("aaVariability", "Amino-acid variability", "Aa_Variability"),
    ("somaticMutations", "Somatic mutations", "Somatic_Mutations"),
    ("somaticHomoplasmy", "Somatic homoplasmy", "SM_Homoplasmy"),
    ("somaticHeteroplasmy", "Somatic heteroplasmy", "SM_Heteroplasmy"),
    ("phastCons", "PhastCons20Way", "PhastCons20Way"),
    ("phyloP", "PhyloP20Way", "PhyloP20Way"),
]

DISEASE_THRESHOLD = 0.80
EVIDENCE_FIELDS = ("mitomapDisease", "clinvar", "omim", "dbsnp")


class BrowserError(RuntimeError):
    pass


class _TextExtractor(HTMLParser):
    def __init__(self) -> None:
        HTMLParser.__init__(self)
        self.parts = []  # type: List[str]

    def handle_data(self, data: str) -> None:
        self.parts.append(data)


def _strip_html(value: Any) -> str:
    parser = _TextExtractor()
    parser.feed(str(value or ""))
    return html.unescape("".join(parser.parts)).strip()


def is_missing(value: Optional[str]) -> bool:
    if value is None:
        return True
    return value.strip().upper() in ("", "NA", "N/A", ".")


def normalize_header(header_fields: Iterable[str]) -> List[str]:
    return [field.strip().replace(" ", "_") for field in header_fields]


def sort_sample_alphabetically(sample_str: str) -> str:
    parts = [sample for sample in sample_str.split(",") if sample]
    if not parts:
        return sample_str
    return ",".join(sorted(parts))


def passes_maf_filter(row: Sequence[str], col_idx: Dict[str, int], maf_cutoff: float) -> bool:
    if maf_cutoff <= 0:
        return True
    index = col_idx.get("AC/AN_1000_Genomes")
    if index is None:
        return True
    value = row[index].strip()
    if not value:
        return True
    try:
        return float(value) < maf_cutoff
    except ValueError:
        return True


def max_hf_any_sample(hf_str: str) -> Optional[float]:
    if is_missing(hf_str):
        return None

    values = []  # type: List[float]
    if "|" not in hf_str and ":" not in hf_str:
        blocks = [hf_str]
    else:
        blocks = [block.split(":", 1)[1] if ":" in block else block for block in hf_str.split("|")]

    for block in blocks:
        for item in block.split(","):
            item = item.strip()
            if is_missing(item):
                continue
            try:
                values.append(float(item))
            except ValueError:
                continue
    return max(values) if values else None


def passes_hf_filter(
    row: Sequence[str],
    col_idx: Dict[str, int],
    hf_cutoff: float,
    drop_missing_hf: bool = False,
) -> bool:
    if hf_cutoff <= 0:
        return True
    index = col_idx.get("HF")
    if index is None:
        return True
    maximum = max_hf_any_sample(row[index].strip())
    if maximum is None:
        return not drop_missing_hf
    return maximum > hf_cutoff


def build_hash_out(
    rows: Iterable[Sequence[str]],
    header: Sequence[str],
    hf_cutoff: float,
    maf_cutoff: float,
    drop_missing_hf: bool = False,
    debug: bool = False,
) -> Dict[str, Dict[str, str]]:
    col_idx = {name: index for index, name in enumerate(header)}
    hash_out = {}  # type: Dict[str, Dict[str, str]]

    for line_no, source_row in enumerate(rows, start=2):
        row = list(source_row)
        if len(row) < len(header):
            if debug:
                print("Skipping short line {}".format(line_no), file=sys.stderr)
            continue

        aa_index = col_idx.get("Aa_Change")
        if aa_index is not None and row[aa_index].strip() == "syn":
            continue
        if not passes_maf_filter(row, col_idx, maf_cutoff):
            continue
        if not passes_hf_filter(row, col_idx, hf_cutoff, drop_missing_hf):
            continue

        sample_index = col_idx.get("Sample")
        variant_index = col_idx.get("Variant_Allele")
        locus_index = col_idx.get("Locus")
        if sample_index is None or variant_index is None or locus_index is None:
            continue

        row[sample_index] = sort_sample_alphabetically(row[sample_index])
        index_key = "{}_{}_{}".format(
            row[sample_index].strip(),
            row[locus_index].strip(),
            row[variant_index].strip(),
        )
        hash_out[index_key] = {
            key: row[index].strip() if index is not None and index < len(row) else ""
            for key in KEYS2REPORT
            for index in [col_idx.get(key)]
        }
    return hash_out


def load_prioritized_variants(
    input_path: Path,
    hf_cutoff: float = 0.30,
    maf_cutoff: float = 0.01,
    keep_missing_hf: bool = False,
    debug: bool = False,
) -> Dict[str, Dict[str, str]]:
    try:
        with input_path.open("r", encoding="utf-8", errors="replace", newline="") as handle:
            reader = csv.reader(handle, delimiter="\t")
            try:
                header = normalize_header(next(reader))
            except StopIteration as exc:
                raise BrowserError("Empty mtDNA input file: <{}>".format(input_path)) from exc
            return build_hash_out(
                reader,
                header,
                hf_cutoff=hf_cutoff,
                maf_cutoff=maf_cutoff,
                drop_missing_hf=not keep_missing_hf,
                debug=debug,
            )
    except OSError as exc:
        raise BrowserError("Cannot read mtDNA input <{}>: {}".format(input_path, exc)) from exc


def _clean_sample(value: str) -> str:
    return value.replace("-DNA_MIT", "")


def _float_or_none(value: Any) -> Optional[float]:
    try:
        return float(str(value).strip())
    except (TypeError, ValueError):
        return None


def _sample_names(value: str) -> List[str]:
    names = []  # type: List[str]
    for item in re.split(r"[,|]", _clean_sample(value)):
        candidate = item.split(":", 1)[0].strip()
        if candidate and candidate not in names:
            names.append(candidate)
    return names


def hash_to_browser_rows(hash_out: Dict[str, Dict[str, str]]) -> List[Dict[str, Any]]:
    rows = []  # type: List[Dict[str, Any]]
    for index_key in sorted(hash_out):
        record = hash_out[index_key]
        row = {
            browser_key: _clean_sample(record.get(source_key, "")) if source_key == "Sample" else record.get(source_key, "")
            for browser_key, _label, source_key in BROWSER_FIELDS
        }
        row["id"] = index_key
        row["maxHeteroplasmy"] = max_hf_any_sample(record.get("HF", ""))
        row["diseaseScoreValue"] = _float_or_none(record.get("Disease_Score", ""))
        row["_samples"] = _sample_names(record.get("Sample", ""))
        row["_hasEvidence"] = any(str(row.get(field, "")).strip() for field in EVIDENCE_FIELDS)
        row["_highDisease"] = (
            row["diseaseScoreValue"] is not None
            and row["diseaseScoreValue"] >= DISEASE_THRESHOLD
        )
        row["_heteroplasmic"] = (
            row["maxHeteroplasmy"] is not None
            and 0 < row["maxHeteroplasmy"] < 1
        )
        rows.append(row)
    return rows


def _coerce_browser_rows(data: Any) -> List[Dict[str, Any]]:
    if not isinstance(data, list):
        raise BrowserError("Browser JSON must contain a top-level 'data' array")
    if all(isinstance(item, dict) for item in data):
        return [dict(item) for item in data]

    rows = []  # type: List[Dict[str, Any]]
    for row_index, item in enumerate(data):
        if not isinstance(item, list):
            raise BrowserError("Browser JSON data entries must be objects or legacy arrays")
        legacy = {
            source_key: _strip_html(item[index]) if index < len(item) else ""
            for index, source_key in enumerate(KEYS4HTML)
        }
        key = "{}_{}_{}".format(
            legacy.get("Sample", ""),
            legacy.get("Locus", ""),
            legacy.get("Variant_Allele", row_index),
        )
        rows.extend(hash_to_browser_rows({key: legacy}))
    return rows


def build_report_payload(
    data_payload: Dict[str, Any],
    *,
    project_id: str,
    job_id: str,
    source_name: str,
    disease_threshold: float = DISEASE_THRESHOLD,
) -> Dict[str, Any]:
    rows = _coerce_browser_rows(data_payload.get("data"))
    loci = sorted({str(row.get("locus", "")).strip() for row in rows if str(row.get("locus", "")).strip()})

    for row in rows:
        if not isinstance(row.get("_samples"), list):
            row["_samples"] = _sample_names(str(row.get("sample", "")))
        if row.get("maxHeteroplasmy") is None:
            row["maxHeteroplasmy"] = max_hf_any_sample(str(row.get("heteroplasmy", "")))
        score = row.get("diseaseScoreValue")
        if score is None:
            score = _float_or_none(row.get("diseaseScore"))
            row["diseaseScoreValue"] = score
        row["_hasEvidence"] = any(str(row.get(field, "")).strip() for field in EVIDENCE_FIELDS)
        row["_highDisease"] = score is not None and score >= disease_threshold
        row["_heteroplasmic"] = (
            row["maxHeteroplasmy"] is not None
            and 0 < row["maxHeteroplasmy"] < 1
        )

    samples = sorted({sample for row in rows for sample in row.get("_samples", [])})

    return {
        "projectId": project_id,
        "jobId": job_id,
        "source": source_name,
        "diseaseThreshold": disease_threshold,
        "columns": [
            {"key": key, "label": label}
            for key, label, _source_key in BROWSER_FIELDS
        ],
        "rows": rows,
        "samples": samples,
        "loci": loci,
        "summary": {
            "variants": len(rows),
            "samples": len(samples),
            "evidence": sum(1 for row in rows if row.get("_hasEvidence")),
            "highDisease": sum(1 for row in rows if row.get("_highDisease")),
            "heteroplasmic": sum(1 for row in rows if row.get("_heteroplasmic")),
        },
    }


def _json_for_script(value: Any) -> str:
    return json.dumps(value, ensure_ascii=False, separators=(",", ":")).replace("</", "<\\/")


def _download_button(href: str, label: str) -> str:
    return (
        '<a class="download-link" href="{}" download>'
        '<svg class="icon" aria-hidden="true"><use href="#icon-download"></use></svg>{}</a>'
    ).format(html.escape(href, quote=True), html.escape(label))


def build_download_buttons(base_dir: str = "../01_mtoolbox") -> str:
    files = [
        (base_dir + "/mit_prioritized_variants.txt", "Report"),
        (base_dir + "/mt_classification_best_results.csv", "Haplogroup"),
        (base_dir + "/VCF_file.vcf", "VCF"),
        (base_dir + "/mit.filtered.json", "Filtered JSON"),
    ]
    return "\n".join(_download_button(href, label) for href, label in files)


def render_report(payload: Dict[str, Any]) -> str:
    try:
        template = TEMPLATE_FILE.read_text(encoding="utf-8")
        tabulator_css = TABULATOR_CSS.read_text(encoding="utf-8")
        tabulator_js = TABULATOR_JS.read_text(encoding="utf-8")
    except OSError as exc:
        raise BrowserError("Cannot read mtDNA browser assets: {}".format(exc)) from exc

    return (
        template.replace("__TABULATOR_CSS__", tabulator_css)
        .replace("__TABULATOR_JS__", tabulator_js.replace("</script", "<\\/script"))
        .replace("__PROJECT_ID__", html.escape(str(payload["projectId"])))
        .replace("__JOB_ID__", html.escape(str(payload["jobId"])))
        .replace("__SOURCE_FILE__", html.escape(str(payload["source"])))
        .replace("__DOWNLOAD_BUTTONS__", build_download_buttons())
        .replace("__REPORT_DATA__", _json_for_script(payload))
    )


def resolve_json_path(json_file: str, html_out: Optional[str]) -> Path:
    json_path = Path(json_file)
    if json_path.is_absolute() or json_path.exists():
        return json_path
    out_dir = Path(html_out).resolve().parent if html_out else Path.cwd()
    colocated = out_dir / json_file
    return colocated if colocated.exists() else json_path


def load_payload(json_file: str, html_out: Optional[str] = None) -> Dict[str, Any]:
    json_path = resolve_json_path(json_file, html_out)
    try:
        with json_path.open("r", encoding="utf-8") as handle:
            payload = json.load(handle)
    except OSError as exc:
        raise BrowserError("Cannot read {}: {}".format(json_path, exc)) from exc
    except json.JSONDecodeError as exc:
        raise BrowserError("Cannot parse {}: {}".format(json_path, exc)) from exc
    if not isinstance(payload, dict) or not isinstance(payload.get("data"), list):
        raise BrowserError("{} must contain a top-level 'data' array".format(json_path))
    return payload


def render_html(
    project_id: str,
    job_id: str,
    json_file: str = "mit.json",
    data_payload: Optional[Dict[str, Any]] = None,
    disease_threshold: float = DISEASE_THRESHOLD,
) -> str:
    report = build_report_payload(
        data_payload or {"data": []},
        project_id=project_id,
        job_id=job_id,
        source_name=json_file,
        disease_threshold=disease_threshold,
    )
    return render_report(report)


def _write_text_atomic(output_path: Path, content: str) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    temporary = output_path.with_name(output_path.name + ".tmp")
    try:
        temporary.write_text(content, encoding="utf-8")
        temporary.replace(output_path)
    except OSError as exc:
        temporary.unlink(missing_ok=True)
        raise BrowserError("Cannot write {}: {}".format(output_path, exc)) from exc


def generate_browser_report(
    input_path: Path,
    output_path: Path,
    *,
    project_id: str,
    job_id: str,
    browser_json_path: Optional[Path] = None,
    hf_cutoff: float = 0.30,
    maf_cutoff: float = 0.01,
    keep_missing_hf: bool = False,
) -> Dict[str, int]:
    hash_out = load_prioritized_variants(
        input_path,
        hf_cutoff=hf_cutoff,
        maf_cutoff=maf_cutoff,
        keep_missing_hf=keep_missing_hf,
    )
    data_payload = {"data": hash_to_browser_rows(hash_out)}
    if browser_json_path is not None:
        _write_text_atomic(browser_json_path, json.dumps(data_payload, ensure_ascii=False) + "\n")
    payload = build_report_payload(
        data_payload,
        project_id=project_id,
        job_id=job_id,
        source_name=input_path.name,
    )
    _write_text_atomic(output_path, render_report(payload))
    return dict(payload["summary"])


def hash2array(hash_out: Dict[str, Dict[str, str]]) -> List[List[str]]:
    return [
        [_html_cell(key, hash_out[index_key].get(key, "")) for key in KEYS4HTML]
        for index_key in sorted(hash_out)
    ]


def _html_link(href: str, text: str) -> str:
    return '<a target="_blank" rel="noopener noreferrer" href="{}">{}</a>'.format(
        html.escape(href, quote=True), html.escape(text)
    )


def _is_http_url(value: str) -> bool:
    parsed = urlparse(value)
    return parsed.scheme in ("http", "https") and bool(parsed.netloc)


def _html_cell(key: str, value: str) -> str:
    temporary = value or ""
    if key == "Mitomap_Associated_Disease(s)" and temporary:
        temporary = temporary.replace("+", "_plus_")
    if key == "Sample" and temporary:
        temporary = _clean_sample(temporary)
    if key in ("GT", "DP", "HF") and temporary:
        return ",<br />".join(html.escape(part) for part in temporary.split("|"))
    if key == "Locus" and temporary:
        return _html_link("https://ghr.nlm.nih.gov/gene/" + quote(temporary, safe="") + "#conditions", temporary)
    if key == "dbSNP_ID" and temporary:
        return _html_link("https://www.ncbi.nlm.nih.gov/snp/" + quote(temporary, safe=""), temporary)
    if key == "OMIM_link" and temporary and _is_http_url(temporary):
        return _html_link(temporary, temporary)
    return html.escape(temporary)


def build_json_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Parse MToolBox prioritized variants and output hash/json/json4html/tsv"
    )
    parser.add_argument("-i", "--input", required=True, type=Path)
    parser.add_argument("-f", "--format", default="hash", choices=["hash", "json", "json4html", "tsv"])
    parser.add_argument("--HF", type=float, default=0.30, help="Heteroplasmic fraction cutoff [0.30]. Set to 0 to disable.")
    parser.add_argument("--MAF", type=float, default=0.01, help="1000G MAF cutoff [0.01]. Set to 0 to disable.")
    parser.add_argument("--denovo", action="store_true", help="Reserved; not implemented")
    parser.add_argument("--debug", action="store_true")
    parser.add_argument("--keep-missing-hf", action="store_true")
    return parser


def json_main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_json_parser().parse_args(argv)
    try:
        hash_out = load_prioritized_variants(
            args.input,
            hf_cutoff=args.HF,
            maf_cutoff=args.MAF,
            keep_missing_hf=args.keep_missing_hf,
            debug=args.debug,
        )
    except BrowserError as exc:
        raise SystemExit(str(exc)) from exc

    if args.format == "hash":
        pprint(hash_out)
    elif args.format == "json":
        print(json.dumps(hash_out, sort_keys=True))
    elif args.format == "json4html":
        print(json.dumps({"data": hash_to_browser_rows(hash_out)}, ensure_ascii=False))
    elif args.format == "tsv":
        writer = csv.writer(sys.stdout, delimiter="\t", lineterminator="\n")
        writer.writerow(KEYS2REPORT)
        for index_key in sorted(hash_out):
            writer.writerow([hash_out[index_key].get(key, "") for key in KEYS2REPORT])
    return 0


def build_html_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Generate a standalone CBIcall mtDNA report")
    parser.add_argument("--id", required=True, help="Project ID displayed in the report")
    parser.add_argument("--job-id", required=True, help="CBIcall run ID displayed in the report")
    parser.add_argument("--json", default="mit.json", help="Browser JSON from mtb2json.py -f json4html")
    parser.add_argument("--out", "-o", default="mtdna.html", help="Output HTML path")
    return parser


def html_main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_html_parser().parse_args(argv)
    try:
        payload = load_payload(args.json, args.out)
        _write_text_atomic(
            Path(args.out),
            render_html(args.id, args.job_id, args.json, payload),
        )
    except BrowserError as exc:
        raise SystemExit(str(exc)) from exc
    return 0


__all__ = [
    "BROWSER_FIELDS",
    "BrowserError",
    "DISEASE_THRESHOLD",
    "KEYS2REPORT",
    "KEYS4HTML",
    "build_download_buttons",
    "build_hash_out",
    "build_report_payload",
    "generate_browser_report",
    "hash2array",
    "hash_to_browser_rows",
    "html_main",
    "is_missing",
    "json_main",
    "load_payload",
    "load_prioritized_variants",
    "max_hf_any_sample",
    "normalize_header",
    "passes_hf_filter",
    "passes_maf_filter",
    "render_html",
    "render_report",
    "resolve_json_path",
    "sort_sample_alphabetically",
]
