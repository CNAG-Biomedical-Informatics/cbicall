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
import urllib.request
from pathlib import Path
from typing import Iterable, List, Optional, Tuple


DEFAULT_BUNDLE_KEY = "cbicall-germline-resources-v1"
CATALOG_NAME = "resources/cbicall-resource-catalog.json"
DEFAULT_CATALOG_URL = (
    "https://raw.githubusercontent.com/mrueda/cbicall/refs/heads/main/"
    "resources/cbicall-resource-catalog.json"
)
MANIFEST_NAME = "cbicall-resource-installation.json"
CATALOG_VERSION = 1


def google_drive_url(file_id: str) -> str:
    return f"https://drive.google.com/uc?export=download&id={file_id}"


def repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


def default_catalog_path() -> Path:
    return repo_root() / CATALOG_NAME


def fetch_catalog(catalog_url: str) -> dict:
    try:
        with urllib.request.urlopen(catalog_url, timeout=30) as handle:
            return json.loads(handle.read().decode("utf-8"))
    except Exception as exc:
        raise SystemExit(
            "Could not load the CBIcall resource catalog.\n"
            f"  URL: {catalog_url}\n\n"
            "Provide a local catalog with --catalog, or check network access."
        ) from exc


def load_catalog(catalog_path: Optional[Path], catalog_url: str) -> Tuple[dict, str]:
    path = catalog_path or default_catalog_path()
    if path.is_file():
        return json.loads(path.read_text(encoding="utf-8")), str(path)

    if catalog_path is not None:
        raise SystemExit(f"Resource catalog not found: {path}")

    print(f"Resource catalog not found locally; fetching {catalog_url}")
    return fetch_catalog(catalog_url), catalog_url


def load_catalog_entry(catalog_path: Optional[Path], bundle_key: str, catalog_url: str = DEFAULT_CATALOG_URL) -> dict:
    catalog, catalog_source = load_catalog(catalog_path, catalog_url)
    bundles = catalog.get("bundles", {}) if isinstance(catalog, dict) else {}
    entry = bundles.get(bundle_key)
    if entry is None:
        raise SystemExit(f"Bundle {bundle_key!r} not found in resource catalog: {catalog_source}")
    entry = dict(entry)
    entry["_catalog_source"] = catalog_source
    if catalog_source.startswith("/"):
        entry["_catalog_path"] = catalog_source
    entry["key"] = bundle_key
    return entry


def safe_slug(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "-", value).strip("-")


def bundle_slug(metadata: dict) -> str:
    if not metadata.get("bundle_id") or not metadata.get("bundle_version"):
        return safe_slug(bundle_key(metadata))

    bundle_id = safe_slug(str(metadata["bundle_id"]))
    bundle_version = safe_slug(str(metadata["bundle_version"]))
    if bundle_version.startswith("v"):
        return f"{bundle_id}-{bundle_version}"
    return f"{bundle_id}-v{bundle_version}"


def canonical_archive_name(metadata: dict) -> str:
    archive = metadata.get("archive") if isinstance(metadata.get("archive"), dict) else {}
    explicit_name = archive.get("canonical_name")
    if explicit_name:
        return Path(str(explicit_name)).name
    return f"{safe_slug(bundle_key(metadata))}.tar.gz"


def archive_source_name(metadata: dict) -> str:
    archive = metadata.get("archive") if isinstance(metadata.get("archive"), dict) else {}
    if not archive.get("source_name"):
        raise SystemExit(f"Bundle {bundle_key(metadata)!r} is missing archive.source_name in the resource catalog.")
    return Path(str(archive["source_name"])).name


def checksum_file_name(metadata: dict) -> str:
    archive = metadata.get("archive") if isinstance(metadata.get("archive"), dict) else {}
    if not archive.get("checksum_file"):
        raise SystemExit(f"Bundle {bundle_key(metadata)!r} is missing archive.checksum_file in the resource catalog.")
    return Path(str(archive["checksum_file"])).name


def checksum_algorithm(metadata: dict) -> str:
    archive = metadata.get("archive") if isinstance(metadata.get("archive"), dict) else {}
    return str(archive.get("checksum_algorithm") or "md5").lower()


def bundle_key(metadata: dict) -> str:
    if metadata.get("key"):
        return str(metadata["key"])
    if metadata.get("bundle_id") and metadata.get("bundle_version"):
        return f"{metadata['bundle_id']}-v{metadata['bundle_version']}"
    return "unknown-bundle"


def source_files(metadata: dict) -> dict:
    source = metadata.get("source") if isinstance(metadata.get("source"), dict) else {}
    files = source.get("files") if isinstance(source.get("files"), dict) else None
    if files is None:
        raise SystemExit(f"Bundle {bundle_key(metadata)!r} is missing source.files in the resource catalog.")
    return dict(files)


def archive_part_names(metadata: dict) -> List[str]:
    archive = metadata.get("archive") if isinstance(metadata.get("archive"), dict) else {}
    parts = archive.get("parts") if isinstance(archive.get("parts"), list) else None
    if parts is not None:
        return list(parts)
    prefix = f"{archive_source_name(metadata)}.part-"
    return sorted(name for name in source_files(metadata) if name.startswith(prefix))


def expected_top_level(metadata: dict) -> List[str]:
    layout = metadata.get("layout") if isinstance(metadata.get("layout"), dict) else {}
    expected = layout.get("expected_top_level") if isinstance(layout.get("expected_top_level"), list) else None
    return list(expected or [])


def remote_identifier(metadata: dict) -> dict:
    identifier = metadata.get("remote_identifier") if isinstance(metadata.get("remote_identifier"), dict) else {}
    return identifier


def identifier_filename(metadata: dict) -> str:
    return Path(str(remote_identifier(metadata).get("filename") or "cbicall-bundle-id.json")).name


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


def require_bundle_identifier(outdir: Path, metadata: dict, bundle_identifier: Optional[dict]) -> dict:
    if bundle_identifier is not None:
        return bundle_identifier

    raise SystemExit(
        "Bundle identifier file was not found.\n"
        f"  expected file: {outdir / identifier_filename(metadata)}\n\n"
        "Either allow the script to download it, or place the file in --outdir "
        "and rerun with --skip-download."
    )


def print_bundle_identifier_summary(metadata: dict, bundle_identifier: dict) -> None:
    observed_bundle = bundle_identifier.get("bundle") or bundle_identifier.get("bundle_key")
    print("Bundle identifier OK.")
    print(f"  Bundle      : {observed_bundle}")
    print(f"  File        : {bundle_identifier.get('_source')}")
    print(f"  SHA-256     : {bundle_identifier.get('_sha256')}")
    print(f"  Catalog key : {bundle_key(metadata)}")


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
    folder = metadata.get("source", {}).get("folder") if isinstance(metadata.get("source"), dict) else None
    if folder:
        print()
        print(f"Folder view: {folder}")
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
    names = [canonical_archive_name(metadata), archive_source_name(metadata)]
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
    checksum_name = checksum_file_name(metadata)
    checksum_file = outdir / checksum_name
    if archive and checksum_file.exists() and not force:
        print(f"{archive.name} and {checksum_name} already exist; skipping source downloads.")
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

    source_name = archive_source_name(metadata)
    archive = outdir / source_name
    if archive.exists() and archive.stat().st_size > 0 and not force:
        print(f"{source_name} already exists; skipping assembly.")
        return archive

    ensure_required_files(outdir, metadata)
    tmp_archive = outdir / f"{source_name}.tmp"

    print(f"Assembling {source_name} from split parts...")
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
    return parse_md5_entries(path)[0]


def parse_md5_entries(path: Path) -> List[Tuple[str, Optional[str]]]:
    entries = []
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
        entries.append((checksum, target))
    if entries:
        return entries
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


def verify_file_md5(path: Path, expected: str) -> str:
    if not path.is_file():
        raise SystemExit(f"Cannot verify checksum because {path} does not exist.")

    actual = compute_md5(path)
    if actual != expected:
        raise SystemExit(
            "Checksum verification failed.\n"
            f"  file: {path}\n"
            f"  expected: {expected}\n"
            f"  observed: {actual}\n\n"
            "Remove the file and download it again."
        )
    return actual


def verify_bundle_checksums(outdir: Path, metadata: dict, archive: Optional[Path] = None) -> dict:
    if checksum_algorithm(metadata) != "md5":
        raise SystemExit(
            f"Unsupported archive checksum algorithm {checksum_algorithm(metadata)!r}; "
            "currently supported: md5."
        )

    checksum_name = checksum_file_name(metadata)
    checksum_file = outdir / checksum_name
    if not checksum_file.is_file():
        raise SystemExit(f"Cannot verify checksum because {checksum_file} does not exist.")

    entries = parse_md5_entries(checksum_file)
    part_names = set(archive_part_names(metadata))
    archive_names = {archive_source_name(metadata), canonical_archive_name(metadata)}
    if archive is not None:
        archive_names.add(archive.name)

    part_entries = [(expected, Path(target).name) for expected, target in entries if target and Path(target).name in part_names]
    archive_entries = [
        (expected, Path(target).name if target else None)
        for expected, target in entries
        if (target and Path(target).name in archive_names) or (target is None and len(entries) == 1 and archive is not None)
    ]

    if part_entries:
        expected_part_names = {name for _, name in part_entries}
        missing_entries = sorted(part_names - expected_part_names)
        missing_files = sorted(name for name in expected_part_names if not (outdir / name).is_file())
        if missing_entries or missing_files:
            details = []
            if missing_entries:
                details.append("missing checksum entries: " + ", ".join(missing_entries))
            if missing_files:
                details.append("missing part files: " + ", ".join(missing_files))
            raise SystemExit(
                "Checksum file covers split archive parts, but the part set is incomplete.\n"
                + "\n".join(f"  {detail}" for detail in details)
                + "\n\nCopy/download the split parts and rerun with --skip-download. "
                "A copied assembled data.tar.gz cannot be verified with this part-level checksum file."
            )

        print(f"Verifying split archive parts with {checksum_name}...")
        verified = []
        for expected, part_name in part_entries:
            actual = verify_file_md5(outdir / part_name, expected)
            verified.append({"filename": part_name, "expected_md5": expected, "observed_md5": actual, "verified": True})
        print("Checksum OK.")
        return {"algorithm": "md5", "scope": "parts", "checksum_file": checksum_name, "verified": True, "entries": verified}

    if archive_entries and archive is not None:
        expected, target_name = archive_entries[0]
        print(f"Verifying {archive.name} with {checksum_name}...")
        actual = verify_file_md5(archive, expected)
        print("Checksum OK.")
        return {
            "algorithm": "md5",
            "scope": "archive",
            "checksum_file": checksum_name,
            "verified": True,
            "entries": [{"filename": target_name or archive.name, "expected_md5": expected, "observed_md5": actual, "verified": True}],
        }

    if archive is not None:
        raise SystemExit(
            "Cannot verify the assembled archive with this checksum file.\n"
            f"  archive: {archive}\n"
            f"  checksum file: {checksum_file}\n\n"
            "The checksum entries do not refer to the assembled archive. "
            "Copy/download the split parts, or provide an archive-level checksum."
        )

    raise SystemExit(f"No usable checksum entries found in {checksum_file}.")


def extraction_looks_done(outdir: Path, metadata: dict) -> bool:
    expected = expected_top_level(metadata)
    return bool(expected) and all((outdir / name).exists() for name in expected)


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
    checksum_result: Optional[dict],
    extracted: bool,
    bundle_identifier: Optional[dict] = None,
    manifest_name: str = MANIFEST_NAME,
) -> Path:
    manifest = outdir / manifest_name
    public_metadata = {key: value for key, value in metadata.items() if not key.startswith("_")}
    payload = {
        "catalog_version": CATALOG_VERSION,
        "bundle_id": metadata.get("bundle_id") or bundle_key(metadata),
        "bundle_version": metadata.get("bundle_version"),
        "bundle_slug": bundle_slug(metadata),
        "catalog_entry": public_metadata,
        "catalog_source": metadata.get("_catalog_source") or metadata.get("_catalog_path"),
        "remote_identifier": bundle_identifier,
        "installed_at_utc": _dt.datetime.utcnow().replace(microsecond=0).isoformat() + "Z",
        "datadir": str(outdir.resolve()),
        "archive": {
            "filename": archive.name,
            "source_filename": archive_source_name(metadata),
            "checksum_file": checksum_file_name(metadata),
            "checksum_algorithm": checksum_algorithm(metadata),
            "verified": bool(checksum_result and checksum_result.get("verified")),
        },
        "checksum": checksum_result,
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
            "folder": metadata.get("source", {}).get("folder"),
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
        "--catalog-url",
        default=DEFAULT_CATALOG_URL,
        help="URL used to fetch the CBIcall resource catalog when the local catalog is not present.",
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
        "--verify-bundle-id-only",
        action="store_true",
        help="Download and verify only the small bundle identifier JSON, then exit.",
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
    metadata = load_catalog_entry(catalog_path, args.bundle, catalog_url=args.catalog_url)
    bundle_id_file_id = args.bundle_id_file_id or identifier_file_id(metadata)
    expected_bundle_id_sha256 = args.expected_bundle_id_sha256 or identifier_sha256(metadata)

    if args.print_manual_download:
        print_manual_download_instructions(outdir, metadata, bundle_id_file_id=bundle_id_file_id)
        return 0

    if not args.skip_download:
        maybe_download_bundle_identifier(outdir, metadata, bundle_id_file_id, force=args.force_download)

    bundle_identifier = validate_bundle_identifier(outdir, metadata, expected_sha256=expected_bundle_id_sha256)
    version_label = metadata.get("bundle_version")
    if version_label:
        print(f"Bundle: {metadata.get('bundle_id') or bundle_key(metadata)} v{version_label}")
    else:
        print(f"Bundle: {bundle_key(metadata)}")

    if args.verify_bundle_id_only:
        bundle_identifier = require_bundle_identifier(outdir, metadata, bundle_identifier)
        print_bundle_identifier_summary(metadata, bundle_identifier)
        return 0

    if not args.skip_download:
        download_files(outdir, metadata, force=args.force_download)

    if args.download_only:
        print("Download-only mode complete.")
        archive = existing_archive(outdir, metadata) or (outdir / archive_source_name(metadata))
        write_manifest(
            outdir,
            metadata,
            archive,
            None,
            extracted=False,
            bundle_identifier=bundle_identifier,
            manifest_name=args.manifest,
        )
        return 0

    archive = assemble_archive(outdir, metadata, force=args.force_assemble)
    checksum_result = verify_bundle_checksums(outdir, metadata, archive)
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
        checksum_result,
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
