#!/usr/bin/env python3
"""Download or install the external CBIcall data bundle.

The Google Drive download can be unreliable for large files. This script is
therefore state-based: users may let it download the files, or download them
manually into the target directory and then use the script for assembly,
checksum verification, extraction, and manifest creation.
"""

from __future__ import annotations

import argparse
import datetime as _dt
import hashlib
import json
import os
import re
import sys
import tarfile
from pathlib import Path
from typing import Iterable, List, Optional, Tuple


BUNDLE_ID = "cbicall-germline-resources"
BUNDLE_VERSION = "1"
CATALOG_VERSION = 1
ARCHIVE_NAME = "data.tar.gz"
CHECKSUM_NAME = "data.tar.gz.md5"
BUNDLE_IDENTIFIER_NAME = "cbicall-bundle-id.json"
DEFAULT_BUNDLE_KEY = "cbicall-germline-resources-v1"
CATALOG_NAME = "resources/cbicall-resource-catalog.json"
MANIFEST_NAME = "cbicall-resource-installation.json"
GOOGLE_DRIVE_FOLDER = "https://drive.google.com/drive/folders/13MqZk0MHN_MQdNyXwjz_QTjbl2Najkeg"

FILES = {
    CHECKSUM_NAME: "1MmeppkZ1xR9ODsBLuDWUh-2ob98Prxq6",
    "data.tar.gz.part-00": "12HQw_duxmlcASh7Yw8yL2-f-hyjqMim9",
    "data.tar.gz.part-01": "1Ejrn1oQ2WO3TSTmWSfoQ1IjYBOXivbRm",
    "data.tar.gz.part-02": "18CC1y2t5heV66-jOHCFIHqgj9GbT3QgF",
    "data.tar.gz.part-03": "1mq5t0U8GFk6ZdaMqBzWHD73LfYJ626vi",
    "data.tar.gz.part-04": "1v8rQA-qtgYgnPcOV8PzO-otw9dE9ePqt",
    "data.tar.gz.part-05": "1knI_DEdrlushYj5lOclnRNDdvSEmXQ5m",
}

PART_NAMES = sorted(name for name in FILES if name.startswith(f"{ARCHIVE_NAME}.part-"))
EXPECTED_TOP_LEVEL = ("Databases", "NGSutils")


def google_drive_url(file_id: str) -> str:
    return f"https://drive.google.com/uc?export=download&id={file_id}"


def repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


def default_catalog_path() -> Path:
    return repo_root() / CATALOG_NAME


def load_catalog_entry(catalog_path: Optional[Path], bundle_key: str) -> Optional[dict]:
    path = catalog_path or default_catalog_path()
    if not path.is_file():
        return None

    catalog = json.loads(path.read_text(encoding="utf-8"))
    bundles = catalog.get("bundles", {}) if isinstance(catalog, dict) else {}
    entry = bundles.get(bundle_key)
    if entry is None:
        raise SystemExit(f"Bundle {bundle_key!r} not found in registry: {path}")
    entry["_catalog_path"] = str(path)
    entry["key"] = bundle_key
    return entry


def default_bundle_metadata() -> dict:
    return {
        "schema_version": CATALOG_VERSION,
        "bundle_id": BUNDLE_ID,
        "bundle_version": BUNDLE_VERSION,
        "description": "CBIcall external tools and reference resources",
        "compatible_cbicall": ">=1.0,<2.0",
        "compatible_workflows": [
            "bash/wes/single/gatk-4.6",
            "bash/wes/cohort/gatk-4.6",
            "bash/wgs/single/gatk-4.6",
            "bash/wgs/cohort/gatk-4.6",
            "bash/mit/single/gatk-3.5",
            "bash/mit/cohort/gatk-3.5",
            "snakemake/wes/single/gatk-4.6",
            "snakemake/wes/cohort/gatk-4.6",
        ],
        "archive": {
            "source_name": ARCHIVE_NAME,
            "canonical_name": f"{BUNDLE_ID}-v{BUNDLE_VERSION}.tar.gz",
            "checksum_file": CHECKSUM_NAME,
            "parts": PART_NAMES,
        },
        "source": {
            "provider": "google_drive",
            "folder": GOOGLE_DRIVE_FOLDER,
            "files": FILES,
        },
        "remote_identifier": {
            "filename": BUNDLE_IDENTIFIER_NAME,
            "google_drive_file_id": None,
            "sha256": None,
            "expected": {
                "bundle": DEFAULT_BUNDLE_KEY,
            },
        },
        "layout": {
            "datadir": ".",
            "dbdir": "Databases",
            "ngsutils": "NGSutils",
            "expected_top_level": list(EXPECTED_TOP_LEVEL),
        },
        "env_bindings": {
            "DATADIR": ".",
            "DBDIR": "Databases",
            "NGSUTILS": "NGSutils",
        },
        "tools": {
            "gatk4": {
                "version": "4.6.2.0",
                "path_hint": "NGSutils/gatk/gatk-4.6.2.0/gatk",
            },
            "bwa": {
                "version": "0.7.18",
                "path_hint": "NGSutils/bwa-0.7.18/bwa",
            },
            "samtools": {
                "version": "0.1.19",
                "path_hint": "NGSutils/samtools-0.1.19/samtools",
            },
        },
        "resource_sets": {
            "b37": {
                "reference_fasta_hint": "Databases/GATK_bundle/b37/references_b37_Homo_sapiens_assembly19.fasta",
                "known_sites_hint": [
                    "Databases/dbSNP/human_9606_b144_GRCh37p13/All_20160408.vcf.gz",
                    "Databases/GATK_bundle/b37/b37_Mills_and_1000G_gold_standard.indels.b37.vcf.gz",
                ],
            },
            "hg38": {
                "reference_fasta_hint": "Databases/GATK_bundle/hg38/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta",
            },
            "rsrs": {
                "description": "Mitochondrial resources used by MToolBox workflows",
            },
        },
    }


def safe_slug(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "-", value).strip("-")


def bundle_slug(metadata: dict) -> str:
    bundle_id = safe_slug(str(metadata.get("bundle_id") or BUNDLE_ID))
    bundle_version = safe_slug(str(metadata.get("bundle_version") or BUNDLE_VERSION))
    if bundle_version.startswith("v"):
        return f"{bundle_id}-{bundle_version}"
    return f"{bundle_id}-v{bundle_version}"


def canonical_archive_name(metadata: dict) -> str:
    archive = metadata.get("archive") if isinstance(metadata.get("archive"), dict) else {}
    explicit_name = archive.get("canonical_name")
    if explicit_name:
        return Path(str(explicit_name)).name
    return f"{bundle_slug(metadata)}.tar.gz"


def bundle_key(metadata: dict) -> str:
    if metadata.get("key"):
        return str(metadata["key"])
    return f"{metadata.get('bundle_id', BUNDLE_ID)}-v{metadata.get('bundle_version', BUNDLE_VERSION)}"


def source_files(metadata: dict) -> dict:
    source = metadata.get("source") if isinstance(metadata.get("source"), dict) else {}
    files = source.get("files") if isinstance(source.get("files"), dict) else None
    return files or FILES


def archive_part_names(metadata: dict) -> List[str]:
    archive = metadata.get("archive") if isinstance(metadata.get("archive"), dict) else {}
    parts = archive.get("parts") if isinstance(archive.get("parts"), list) else None
    return list(parts or PART_NAMES)


def expected_top_level(metadata: dict) -> List[str]:
    layout = metadata.get("layout") if isinstance(metadata.get("layout"), dict) else {}
    expected = layout.get("expected_top_level") if isinstance(layout.get("expected_top_level"), list) else None
    return list(expected or EXPECTED_TOP_LEVEL)


def remote_identifier(metadata: dict) -> dict:
    identifier = metadata.get("remote_identifier") if isinstance(metadata.get("remote_identifier"), dict) else {}
    return identifier


def identifier_filename(metadata: dict) -> str:
    return Path(str(remote_identifier(metadata).get("filename") or BUNDLE_IDENTIFIER_NAME)).name


def identifier_file_id(metadata: dict) -> Optional[str]:
    return remote_identifier(metadata).get("google_drive_file_id")


def identifier_sha256(metadata: dict) -> Optional[str]:
    return remote_identifier(metadata).get("sha256")


def maybe_download_bundle_identifier(outdir: Path, metadata: dict, file_id: Optional[str], force: bool = False) -> None:
    if not file_id:
        return
    download_if_missing(identifier_filename(metadata), file_id, outdir, force=force)


def validate_bundle_identifier(outdir: Path, metadata: dict, expected_sha256: Optional[str] = None) -> Optional[dict]:
    path = outdir / identifier_filename(metadata)
    if not path.is_file():
        return None

    observed_sha256 = compute_sha256(path)
    expected_sha256 = expected_sha256 or identifier_sha256(metadata)
    if expected_sha256 and observed_sha256 != expected_sha256:
        raise SystemExit(
            "Bundle identifier SHA-256 verification failed.\n"
            f"  expected: {expected_sha256}\n"
            f"  observed: {observed_sha256}"
        )

    text = path.read_text(encoding="utf-8").strip()
    try:
        payload = json.loads(text)
    except json.JSONDecodeError:
        payload = {"bundle": text}
    if isinstance(payload, str):
        payload = {"bundle": payload}
    if not isinstance(payload, dict):
        raise SystemExit(f"Bundle identifier must be a string or JSON object: {path}")

    expected_bundle = remote_identifier(metadata).get("expected", {}).get("bundle") or bundle_key(metadata)
    observed_bundle = payload.get("bundle") or payload.get("bundle_key")
    if observed_bundle != expected_bundle:
        raise SystemExit(
            "Bundle identifier does not match the selected CBIcall catalog entry.\n"
            f"  expected: {expected_bundle}\n"
            f"  observed: {observed_bundle}"
        )

    payload["_source"] = str(path)
    payload["_sha256"] = observed_sha256
    return payload


def print_manual_download_instructions(outdir: Path, metadata: dict, bundle_id_file_id: Optional[str] = None) -> None:
    print("Manual download mode")
    print("====================")
    print()
    print(f"Download these files into: {outdir.resolve()}")
    print()
    file_id = bundle_id_file_id or identifier_file_id(metadata)
    if file_id:
        print(f"- {identifier_filename(metadata)}")
        print(f"  {google_drive_url(file_id)}")
    for filename, file_id in source_files(metadata).items():
        print(f"- {filename}")
        print(f"  {google_drive_url(file_id)}")
    print()
    print(f"Folder view: {GOOGLE_DRIVE_FOLDER}")
    print()
    print("After all files are present, run:")
    print(f"  python3 {Path(__file__).resolve()} --outdir {outdir.resolve()} --skip-download")


def download_if_missing(filename: str, file_id: str, outdir: Path, force: bool = False) -> None:
    output = outdir / filename
    if output.exists() and output.stat().st_size > 0 and not force:
        print(f"{filename} already exists; skipping download.")
        return

    try:
        import gdown
    except ImportError as exc:
        raise SystemExit(
            "Python package 'gdown' is required for automatic download.\n"
            "Install it with 'pip3 install gdown', or download the files manually "
            "and rerun with --skip-download."
        ) from exc

    url = google_drive_url(file_id)
    print(f"Downloading {filename}...")
    result = gdown.download(url, str(output), quiet=False)
    if result is None or not output.exists() or output.stat().st_size == 0:
        raise RuntimeError(
            f"Download did not produce a usable {filename}. "
            "Download the files manually and rerun with --skip-download."
        )


def archive_candidates(outdir: Path, metadata: dict) -> List[Path]:
    names = [canonical_archive_name(metadata), ARCHIVE_NAME]
    candidates = []
    for name in names:
        path = outdir / name
        if path not in candidates:
            candidates.append(path)
    return candidates


def existing_archive(outdir: Path, metadata: dict) -> Optional[Path]:
    for archive in archive_candidates(outdir, metadata):
        if archive.exists() and archive.stat().st_size > 0:
            return archive
    return None


def download_files(outdir: Path, metadata: dict, force: bool = False) -> None:
    archive = existing_archive(outdir, metadata)
    checksum_file = outdir / CHECKSUM_NAME
    if archive and checksum_file.exists() and not force:
        print(f"{archive.name} and {CHECKSUM_NAME} already exist; skipping source downloads.")
        return

    for filename, file_id in source_files(metadata).items():
        download_if_missing(filename, file_id, outdir, force=force)


def ensure_required_files(outdir: Path, metadata: dict) -> None:
    missing = [name for name in source_files(metadata) if not (outdir / name).is_file()]
    if missing:
        missing_text = "\n".join(f"  - {name}" for name in missing)
        raise SystemExit(
            "Missing required bundle files:\n"
            f"{missing_text}\n\n"
            "Either rerun without --skip-download, or download the missing files "
            "manually into the output directory."
        )


def assemble_archive(outdir: Path, metadata: dict, force: bool = False) -> Path:
    existing = existing_archive(outdir, metadata)
    if existing and not force:
        print(f"{existing.name} already exists; skipping assembly.")
        return existing

    archive = outdir / ARCHIVE_NAME
    if archive.exists() and archive.stat().st_size > 0 and not force:
        print(f"{ARCHIVE_NAME} already exists; skipping assembly.")
        return archive

    ensure_required_files(outdir, metadata)
    tmp_archive = outdir / f"{ARCHIVE_NAME}.tmp"

    print(f"Assembling {ARCHIVE_NAME} from split parts...")
    with tmp_archive.open("wb") as output:
        for part_name in archive_part_names(metadata):
            part = outdir / part_name
            print(f"  adding {part_name}")
            with part.open("rb") as handle:
                while True:
                    chunk = handle.read(1024 * 1024)
                    if not chunk:
                        break
                    output.write(chunk)

    tmp_archive.replace(archive)
    return archive


def parse_md5_file(path: Path) -> Tuple[str, Optional[str]]:
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line:
            continue
        match = re.search(r"\b([a-fA-F0-9]{32})\b", line)
        if not match:
            continue
        checksum = match.group(1).lower()
        remainder = line[match.end():].strip()
        target = remainder.lstrip("*").strip() or None
        return checksum, target
    raise ValueError(f"No MD5 checksum found in {path}")


def compute_md5(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        while True:
            chunk = handle.read(8 * 1024 * 1024)
            if not chunk:
                break
            digest.update(chunk)
    return digest.hexdigest()


def compute_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while True:
            chunk = handle.read(8 * 1024 * 1024)
            if not chunk:
                break
            digest.update(chunk)
    return digest.hexdigest()


def verify_archive(outdir: Path, archive: Optional[Path] = None) -> Tuple[str, str]:
    archive = archive or (outdir / ARCHIVE_NAME)
    checksum_file = outdir / CHECKSUM_NAME
    if not archive.is_file():
        raise SystemExit(f"Cannot verify checksum because {archive} does not exist.")
    if not checksum_file.is_file():
        raise SystemExit(f"Cannot verify checksum because {checksum_file} does not exist.")

    expected, target = parse_md5_file(checksum_file)
    if target and Path(target).name not in {ARCHIVE_NAME, archive.name}:
        print(f"Warning: checksum file refers to {target!r}, verifying {archive.name!r}.")

    print(f"Verifying {archive.name} with {CHECKSUM_NAME}...")
    actual = compute_md5(archive)
    if actual != expected:
        raise SystemExit(
            "Checksum verification failed.\n"
            f"  expected: {expected}\n"
            f"  observed: {actual}\n\n"
            "Remove the archive and rerun assembly, or download the parts again."
        )

    print("Checksum OK.")
    return expected, actual


def extraction_looks_done(outdir: Path, metadata: dict) -> bool:
    return all((outdir / name).exists() for name in expected_top_level(metadata))


def _safe_members(tar: tarfile.TarFile, destination: Path) -> Iterable[tarfile.TarInfo]:
    destination = destination.resolve()
    for member in tar.getmembers():
        member_path = (destination / member.name).resolve()
        if os.path.commonpath([str(destination), str(member_path)]) != str(destination):
            raise RuntimeError(f"Unsafe path in tar archive: {member.name}")
        yield member


def canonicalize_archive(outdir: Path, archive: Path, metadata: dict) -> Path:
    canonical = outdir / canonical_archive_name(metadata)
    if archive.resolve() == canonical.resolve():
        return archive
    if canonical.exists() and canonical.stat().st_size > 0:
        print(f"Using existing canonical archive: {canonical.name}")
        return canonical
    archive.replace(canonical)
    print(f"Renamed archive to {canonical.name}.")
    return canonical


def extract_archive(outdir: Path, archive: Path, metadata: dict, force: bool = False) -> bool:
    if extraction_looks_done(outdir, metadata) and not force:
        expected = ", ".join(expected_top_level(metadata))
        print(f"Extraction appears complete ({expected}); skipping extraction.")
        return False

    print(f"Extracting {archive.name} into {outdir.resolve()}...")
    with tarfile.open(archive, "r:gz") as tar:
        tar.extractall(outdir, members=_safe_members(tar, outdir))
    return True


def remove_parts(outdir: Path, metadata: dict) -> None:
    for part_name in archive_part_names(metadata):
        part = outdir / part_name
        if part.exists():
            part.unlink()
            print(f"Removed {part_name}.")


def write_manifest(
    outdir: Path,
    metadata: dict,
    archive: Path,
    expected_md5: Optional[str],
    observed_md5: Optional[str],
    extracted: bool,
    bundle_identifier: Optional[dict] = None,
    manifest_name: str = MANIFEST_NAME,
) -> Path:
    manifest = outdir / manifest_name
    public_metadata = {key: value for key, value in metadata.items() if not key.startswith("_")}
    payload = {
        "catalog_version": CATALOG_VERSION,
        "bundle_id": metadata.get("bundle_id", BUNDLE_ID),
        "bundle_version": metadata.get("bundle_version", BUNDLE_VERSION),
        "bundle_slug": bundle_slug(metadata),
        "catalog_entry": public_metadata,
        "catalog_source": metadata.get("_catalog_path"),
        "remote_identifier": bundle_identifier,
        "installed_at_utc": _dt.datetime.utcnow().replace(microsecond=0).isoformat() + "Z",
        "datadir": str(outdir.resolve()),
        "archive": {
            "filename": archive.name,
            "source_filename": ARCHIVE_NAME,
            "checksum_file": CHECKSUM_NAME,
            "expected_md5": expected_md5,
            "observed_md5": observed_md5,
            "verified": bool(expected_md5 and observed_md5 and expected_md5 == observed_md5),
        },
        "parts": [
            {
                "filename": part_name,
                "present": (outdir / part_name).exists(),
                "size_bytes": (outdir / part_name).stat().st_size if (outdir / part_name).exists() else None,
            }
            for part_name in archive_part_names(metadata)
        ],
        "expected_top_level": expected_top_level(metadata),
        "extracted": extracted or extraction_looks_done(outdir, metadata),
        "source": {
            "type": metadata.get("source", {}).get("provider", "google_drive"),
            "folder": metadata.get("source", {}).get("folder", GOOGLE_DRIVE_FOLDER),
            "file_ids": source_files(metadata),
        },
    }
    manifest.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"Wrote installation manifest: {manifest}")
    return manifest


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Download, assemble, verify, and extract the CBIcall external data bundle.",
    )
    parser.add_argument(
        "--outdir",
        default=".",
        help="Directory where the bundle parts, archive, and extracted data live (default: current directory).",
    )
    parser.add_argument(
        "--catalog",
        help=f"CBIcall resource catalog JSON (default: {CATALOG_NAME}).",
    )
    parser.add_argument(
        "--bundle",
        default=DEFAULT_BUNDLE_KEY,
        help=f"Bundle key from the CBIcall resource catalog (default: {DEFAULT_BUNDLE_KEY}).",
    )
    parser.add_argument(
        "--print-manual-download",
        action="store_true",
        help="Print file names and Google Drive URLs for manual download, then exit.",
    )
    parser.add_argument(
        "--skip-download",
        action="store_true",
        help="Do not use Python/gdown; expect all required files to already exist in --outdir.",
    )
    parser.add_argument(
        "--bundle-id-file-id",
        help="Google Drive file ID for the optional small bundle identifier file.",
    )
    parser.add_argument(
        "--expected-bundle-id-sha256",
        help="Expected SHA-256 of the optional small bundle identifier file.",
    )
    parser.add_argument(
        "--download-only",
        action="store_true",
        help="Download missing files, then stop before assembly and extraction.",
    )
    parser.add_argument(
        "--no-extract",
        action="store_true",
        help="Assemble and verify the archive, but do not extract it.",
    )
    parser.add_argument(
        "--remove-parts",
        action="store_true",
        help="Remove data.tar.gz.part-* files after the archive has been verified.",
    )
    parser.add_argument(
        "--force-download",
        action="store_true",
        help="Download files again even if they already exist.",
    )
    parser.add_argument(
        "--force-assemble",
        action="store_true",
        help="Recreate data.tar.gz even if it already exists.",
    )
    parser.add_argument(
        "--force-extract",
        action="store_true",
        help="Extract even if Databases and NGSutils already exist.",
    )
    parser.add_argument(
        "--manifest",
        default=MANIFEST_NAME,
        help=f"Installation manifest filename (default: {MANIFEST_NAME}).",
    )
    return parser


def main(argv: Optional[List[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    outdir = Path(args.outdir).expanduser()
    outdir.mkdir(parents=True, exist_ok=True)
    catalog_path = Path(args.catalog).expanduser() if args.catalog else None
    metadata = load_catalog_entry(catalog_path, args.bundle) or default_bundle_metadata()
    bundle_id_file_id = args.bundle_id_file_id or identifier_file_id(metadata)
    expected_bundle_id_sha256 = args.expected_bundle_id_sha256 or identifier_sha256(metadata)

    if args.print_manual_download:
        print_manual_download_instructions(outdir, metadata, bundle_id_file_id=bundle_id_file_id)
        return 0

    if not args.skip_download:
        maybe_download_bundle_identifier(outdir, metadata, bundle_id_file_id, force=args.force_download)

    bundle_identifier = validate_bundle_identifier(outdir, metadata, expected_sha256=expected_bundle_id_sha256)
    print(f"Bundle: {metadata.get('bundle_id', BUNDLE_ID)} v{metadata.get('bundle_version', BUNDLE_VERSION)}")

    if not args.skip_download:
        download_files(outdir, metadata, force=args.force_download)

    if args.download_only:
        print("Download-only mode complete.")
        archive = existing_archive(outdir, metadata) or (outdir / ARCHIVE_NAME)
        write_manifest(
            outdir,
            metadata,
            archive,
            None,
            None,
            extracted=False,
            bundle_identifier=bundle_identifier,
            manifest_name=args.manifest,
        )
        return 0

    archive = assemble_archive(outdir, metadata, force=args.force_assemble)
    expected_md5, observed_md5 = verify_archive(outdir, archive)
    archive = canonicalize_archive(outdir, archive, metadata)

    if args.remove_parts:
        remove_parts(outdir, metadata)

    extracted = False
    if not args.no_extract:
        extracted = extract_archive(outdir, archive, metadata, force=args.force_extract)

    write_manifest(
        outdir,
        metadata,
        archive,
        expected_md5,
        observed_md5,
        extracted=extracted,
        bundle_identifier=bundle_identifier,
        manifest_name=args.manifest,
    )

    print()
    print("Resource setup complete.")
    print(f"Set DATADIR to: {outdir.resolve()}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
