import argparse
import gzip
import hashlib
import json
import os
import re
import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import yaml


class IntegrationTestError(RuntimeError):
    """Raised when an integration contract fails."""


@dataclass(frozen=True)
class TestSelection:
    key: str
    label: str
    fixture: str
    engine: Optional[str] = None
    optional_in_all: bool = False


TESTS: Dict[str, TestSelection] = {
    "wes-bash": TestSelection("wes-bash", "WES Bash", "native-wes-bash.yaml"),
    "wes-snakemake": TestSelection(
        "wes-snakemake",
        "WES Snakemake",
        "native-wes-snakemake.yaml",
        engine="snakemake",
        optional_in_all=True,
    ),
    "wes-nextflow": TestSelection(
        "wes-nextflow",
        "WES Nextflow",
        "native-wes-nextflow.yaml",
        engine="nextflow",
        optional_in_all=True,
    ),
    "mit-bash": TestSelection("mit-bash", "MIT Bash", "native-mit-bash.yaml"),
    "nf-core-demo": TestSelection(
        "nf-core-demo",
        "nf-core Demo",
        "nf-core-demo.yaml",
        engine="nextflow",
    ),
    "nf-core-sarek": TestSelection(
        "nf-core-sarek",
        "nf-core Sarek",
        "nf-core-sarek.yaml",
        engine="nextflow",
    ),
}


def fixture_dir(project_root: Path) -> Path:
    return project_root / "tests" / "fixtures" / "integration"


def load_contract(project_root: Path, selection: TestSelection) -> dict:
    path = fixture_dir(project_root) / selection.fixture
    if not path.is_file():
        raise FileNotFoundError(f"Integration contract not found: {path}")
    data = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    if not isinstance(data, dict):
        raise IntegrationTestError(f"Integration contract must be a mapping: {path}")
    data["_contract_path"] = str(path)
    return data


def _project_path(project_root: Path, value: str) -> Path:
    path = Path(value)
    if path.is_absolute():
        return path
    return project_root / path


def _work_path(workdir: Path, value: str) -> Path:
    path = Path(value)
    if path.is_absolute():
        return path
    return workdir / path


def list_run_dirs(base_dir: Path, run_glob: str) -> List[Path]:
    if not base_dir.is_dir():
        return []
    return sorted(path for path in base_dir.glob(run_glob) if path.is_dir())


def _parse_run_dir_from_stdout(stdout: str) -> Optional[Path]:
    for pattern in (
        r"^Working directory:\s*(.+)$",
        r"^\s*Report\s*=>\s*(.+/run-report\.json)\s*$",
        r"^\s*Log\s*=>\s*(.+)\s*$",
    ):
        matches = re.findall(pattern, stdout, flags=re.MULTILINE)
        if matches:
            value = matches[-1].strip()
            path = Path(value)
            if path.name == "run-report.json":
                return path.parent
            if path.suffix == ".log":
                return path.parent
            return path
    return None


def _new_run_dir(before: Iterable[Path], after: Iterable[Path]) -> Optional[Path]:
    before_set = {path.resolve() for path in before}
    new_dirs = [path for path in after if path.resolve() not in before_set]
    return sorted(new_dirs)[-1] if new_dirs else None


def _json_path(data, path: str):
    current = data
    for part in path.split("."):
        if isinstance(current, dict) and part in current:
            current = current[part]
        else:
            raise IntegrationTestError(f"Missing JSON path: {path}")
    return current


def _check_json_expectations(report_path: Path, expectations: List[dict]) -> None:
    data = json.loads(report_path.read_text(encoding="utf-8"))
    for item in expectations:
        path = item["path"]
        observed = _json_path(data, path)
        if "equals" in item and observed != item["equals"]:
            raise IntegrationTestError(
                f"JSON expectation failed for {path}: expected {item['equals']!r}, observed {observed!r}"
            )
        if "endswith" in item and not str(observed).endswith(str(item["endswith"])):
            raise IntegrationTestError(
                f"JSON expectation failed for {path}: expected suffix {item['endswith']!r}, observed {observed!r}"
            )
        if "contains" in item and str(item["contains"]) not in str(observed):
            raise IntegrationTestError(
                f"JSON expectation failed for {path}: expected to contain {item['contains']!r}, observed {observed!r}"
            )


def _read_text_lines(path: Path) -> List[str]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
        return [line.rstrip("\n") for line in handle]


def _normalized_text_hash(path: Path, pattern: str = "^#") -> Tuple[str, int]:
    regex = re.compile(pattern)
    lines = [line for line in _read_text_lines(path) if not regex.search(line)]
    lines = sorted(lines)
    payload = ("\n".join(lines) + "\n").encode("utf-8") if lines else b""
    return hashlib.sha256(payload).hexdigest(), len(lines)


def _canonical_json_hash(path: Path) -> str:
    payload = json.loads(path.read_text(encoding="utf-8"))
    canonical = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(canonical).hexdigest()


def _check_hashes(run_dir: Path, hashes: List[dict]) -> None:
    for item in hashes:
        path = _work_path(run_dir, item["path"])
        if not path.is_file():
            raise IntegrationTestError(f"Hash target does not exist: {path}")
        kind = item.get("type", "normalized_text")
        if kind in {"normalized_text", "normalized_vcf"}:
            observed, records = _normalized_text_hash(path, item.get("pattern", "^#"))
            expected_records = item.get("records")
            if expected_records is not None and int(expected_records) != records:
                raise IntegrationTestError(
                    f"Record count mismatch for {path}: expected {expected_records}, observed {records}"
                )
        elif kind == "canonical_json":
            observed = _canonical_json_hash(path)
        else:
            raise IntegrationTestError(f"Unsupported hash type in contract: {kind}")

        expected = item["sha256"]
        print(f"  Hash {item.get('name', path.name)}: {observed}")
        if observed != expected:
            raise IntegrationTestError(
                f"Hash mismatch for {path}: expected {expected}, observed {observed}"
            )


def validate_contract(run_dir: Path, contract: dict) -> None:
    errors: List[str] = []
    for rel_path in contract.get("required_files", []):
        path = _work_path(run_dir, rel_path)
        if not path.is_file():
            errors.append(f"missing file: {rel_path}")

    for pattern in contract.get("required_globs", []):
        matches = sorted(run_dir.glob(pattern))
        if not matches:
            errors.append(f"missing glob: {pattern}")

    report = run_dir / "run-report.json"
    if report.is_file():
        try:
            _check_json_expectations(report, contract.get("json_expectations", []))
        except IntegrationTestError as exc:
            errors.append(str(exc))
    elif contract.get("json_expectations"):
        errors.append("missing file: run-report.json")

    if errors:
        raise IntegrationTestError("\n".join(errors))

    _check_hashes(run_dir, contract.get("hashes", []))


def _write_inline_parameters(workdir: Path, key: str, payload: dict) -> Path:
    handle = tempfile.NamedTemporaryFile(
        "w",
        encoding="utf-8",
        prefix=f"cbicall-{key}.",
        suffix=".yaml",
        dir=workdir,
        delete=False,
    )
    with handle:
        yaml.safe_dump(payload, handle, sort_keys=False)
    return Path(handle.name)


def _parameter_file(project_root: Path, workdir: Path, key: str, contract: dict) -> Tuple[Path, Optional[Path]]:
    run = contract.get("run", {})
    if "parameter_file" in run:
        return _work_path(workdir, run["parameter_file"]), None
    if "parameters" in run:
        path = _write_inline_parameters(workdir, key, run["parameters"])
        return path, path
    raise IntegrationTestError("Contract run section must define parameter_file or parameters")


def _print_artifacts(label: str, run_dir: Path, contract: dict, launcher_log: Path) -> None:
    workflow_log = run_dir / contract["workflow_log"]
    print()
    print(f"{label} run artifacts:")
    for name, path in (
        ("Run dir", run_dir),
        ("Workflow log", workflow_log),
        ("Run report", run_dir / "run-report.json"),
        ("HTML report", run_dir / "run-report.html"),
        ("Launcher log", launcher_log),
        ("Contract", Path(contract["_contract_path"])),
    ):
        suffix = "" if path.exists() else " (not created)"
        print(f"  {name:<12} : {path}{suffix}")
    print()


def _cleanup(run_dir: Path, contract: dict, keep_external_work: bool) -> None:
    cleanup_paths = contract.get("cleanup_paths", [])
    if keep_external_work or not cleanup_paths:
        return
    print("Cleaning heavy execution state:")
    for rel_path in cleanup_paths:
        path = _work_path(run_dir, rel_path)
        if path.exists():
            if path.is_dir():
                shutil.rmtree(path)
            else:
                path.unlink()
            print(f"  removed {rel_path}")


def _run_one(
    *,
    project_root: Path,
    selection: TestSelection,
    threads: int,
    runtime_profile: str,
    skip_missing_optional: bool,
    keep_external_work: bool,
) -> Tuple[str, str, str]:
    if selection.engine and shutil.which(selection.engine) is None:
        detail = f"{selection.engine} not found"
        if skip_missing_optional and selection.optional_in_all:
            print(f"SKIP: {selection.label} requires {selection.engine} on PATH.")
            return selection.label, "skipped", detail
        raise IntegrationTestError(f"{selection.label} requires {selection.engine} on PATH.")

    contract = load_contract(project_root, selection)
    workdir = _project_path(project_root, contract.get("workdir", "examples/input"))
    base_dir = _work_path(workdir, contract.get("run", {}).get("base_dir", "."))
    run_glob = contract.get("run", {}).get("run_glob")
    if not run_glob:
        raise IntegrationTestError("Contract run section must define run_glob")

    param_file, temp_param = _parameter_file(project_root, workdir, selection.key, contract)
    launcher_fd, launcher_name = tempfile.mkstemp(
        prefix=f"cbicall-test-{selection.key}.",
        suffix=".log",
        dir=workdir,
        text=True,
    )
    os.close(launcher_fd)
    launcher_log = Path(launcher_name)

    before = list_run_dirs(base_dir, run_glob)
    cmd = [
        str(project_root / "bin" / "cbicall"),
        "run",
        "-p",
        str(param_file),
        "-t",
        str(threads),
        "--runtime-profile",
        str(runtime_profile),
    ]

    print("========================================")
    print(f"TEST: {selection.label}")
    print("========================================")
    print(f"Running {selection.label} integration test...")

    proc = subprocess.run(
        cmd,
        cwd=str(workdir),
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        check=False,
    )
    launcher_log.write_text(proc.stdout, encoding="utf-8")

    after = list_run_dirs(base_dir, run_glob)
    run_dir = _parse_run_dir_from_stdout(proc.stdout)
    if run_dir is not None and not run_dir.is_absolute():
        run_dir = workdir / run_dir
    if run_dir is None or not run_dir.exists():
        run_dir = _new_run_dir(before, after)

    try:
        if proc.returncode != 0:
            tail = "\n".join(proc.stdout.splitlines()[-40:])
            raise IntegrationTestError(
                f"{selection.label} cbicall command failed with return code {proc.returncode}.\n"
                f"Launcher log: {launcher_log}\n{tail}"
            )
        if run_dir is None:
            raise IntegrationTestError(f"{selection.label} finished but no run directory was found.")

        _print_artifacts(selection.label, run_dir, contract, launcher_log)
        print("Validating contract:")
        validate_contract(run_dir, contract)
        print("  OK required files")
        print("  OK required globs")
        print("  OK run-report fields")
        if contract.get("hashes"):
            print("  OK output hashes")
        _cleanup(run_dir, contract, keep_external_work)
        print(f"SUCCESS: {selection.label} integration contract passed.")
        return selection.label, "passed", str(run_dir)
    finally:
        if temp_param:
            temp_param.unlink(missing_ok=True)


def run_integration_tests(
    *,
    project_root: Path,
    selected: List[TestSelection],
    threads: int,
    runtime_profile: str,
    skip_missing_optional: bool = False,
    keep_external_work: bool = False,
) -> int:
    overall_status = 0
    summary: List[Tuple[str, str, str]] = []

    for selection in selected:
        try:
            summary.append(
                _run_one(
                    project_root=project_root,
                    selection=selection,
                    threads=threads,
                    runtime_profile=runtime_profile,
                    skip_missing_optional=skip_missing_optional,
                    keep_external_work=keep_external_work,
                )
            )
        except IntegrationTestError as exc:
            overall_status = 1
            print(f"ERROR: {exc}")
            summary.append((selection.label, "failed", str(exc).splitlines()[0]))

    print("========================================")
    print("Integration test summary")
    print("========================================")
    for label, status, detail in summary:
        suffix = f" ({detail})" if detail else ""
        print(f"{label}: {status}{suffix}")
    print("========================================")
    print("All requested tests finished.")
    print(f"Exit code: {overall_status}")
    print("========================================")
    return overall_status


def selected_tests_from_args(args: argparse.Namespace) -> List[TestSelection]:
    selected: List[TestSelection] = []
    if args.all or args.wes_bash:
        selected.append(TESTS["wes-bash"])
    if args.all or args.wes_snakemake:
        selected.append(TESTS["wes-snakemake"])
    if args.all or args.wes_nextflow:
        selected.append(TESTS["wes-nextflow"])
    if args.all or args.mit_bash:
        selected.append(TESTS["mit-bash"])
    if args.nf_core_demo:
        selected.append(TESTS["nf-core-demo"])
    if args.nf_core_sarek:
        selected.append(TESTS["nf-core-sarek"])
    return selected
