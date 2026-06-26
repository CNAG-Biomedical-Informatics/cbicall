import argparse
import gzip
import hashlib
import json
import os
import re
import shlex
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
    backend_executable: Optional[str] = None
    optional_in_all: bool = False


TESTS: Dict[str, TestSelection] = {
    "wes-bash": TestSelection("wes-bash", "WES Bash", "native-wes-bash.yaml"),
    "wes-cohort-bash": TestSelection(
        "wes-cohort-bash",
        "WES Cohort Bash",
        "native-wes-cohort-bash.yaml",
    ),
    "wes-cohort-bash-sharded": TestSelection(
        "wes-cohort-bash-sharded",
        "WES Cohort Bash Sharded",
        "native-wes-cohort-bash-sharded.yaml",
    ),
    "wes-bash-gatk35": TestSelection("wes-bash-gatk35", "WES Bash GATK 3.5", "native-wes-bash-gatk35.yaml"),
    "wes-snakemake": TestSelection(
        "wes-snakemake",
        "WES Snakemake",
        "native-wes-snakemake.yaml",
        backend_executable="snakemake",
        optional_in_all=True,
    ),
    "wes-nextflow": TestSelection(
        "wes-nextflow",
        "WES Nextflow",
        "native-wes-nextflow.yaml",
        backend_executable="nextflow",
        optional_in_all=True,
    ),
    "wes-cromwell": TestSelection(
        "wes-cromwell",
        "WES Cromwell",
        "native-wes-cromwell.yaml",
        optional_in_all=True,
    ),
    "mit-bash": TestSelection("mit-bash", "MIT Bash", "native-mit-bash.yaml"),
    "nf-core-demo": TestSelection(
        "nf-core-demo",
        "nf-core Demo",
        "nf-core-demo.yaml",
        backend_executable="nextflow",
    ),
    "nf-core-sarek": TestSelection(
        "nf-core-sarek",
        "nf-core Sarek",
        "nf-core-sarek.yaml",
        backend_executable="nextflow",
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


def _contract_env(project_root: Path, contract: dict) -> Dict[str, str]:
    values = contract.get("env") or {}
    if not isinstance(values, dict):
        raise IntegrationTestError("Integration contract env must be a mapping.")

    env = os.environ.copy()
    for key, value in values.items():
        if not isinstance(key, str) or not key:
            raise IntegrationTestError("Integration contract env keys must be non-empty strings.")
        if value is None:
            env.pop(key, None)
            continue

        text = str(value)
        path = Path(text)
        project_value = project_root / path
        if not path.is_absolute() and project_value.exists():
            text = str(project_value)
        env[key] = text
    return env


def _resolve_bcftools(project_root: Path, env: Dict[str, str]) -> str:
    bcftools = env.get("BCFTOOLS")
    if not bcftools:
        env_file = project_root / "workflows" / "bash" / "gatk-4.6" / "env.sh"
        script = f"source {shlex.quote(str(env_file))}; printf '%s\\n' \"${{BCFTOOLS:-bcftools}}\""
        proc = subprocess.run(
            ["bash", "-lc", script],
            cwd=str(project_root),
            env=env,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            check=False,
        )
        if proc.returncode != 0:
            raise IntegrationTestError(f"Could not resolve BCFTOOLS from {env_file}: {proc.stdout.strip()}")
        bcftools = proc.stdout.strip().splitlines()[-1] if proc.stdout.strip() else "bcftools"

    if os.sep in bcftools:
        if not Path(bcftools).is_file():
            raise IntegrationTestError(f"BCFTOOLS does not exist: {bcftools}")
        return bcftools
    if shutil.which(bcftools, path=env.get("PATH")) is None:
        raise IntegrationTestError(
            "bcftools is required only for the staged cohort gather/finalize integration test. "
            "CBIcall resolved BCFTOOLS to 'bcftools' from the environment/env.sh, but it is not on PATH. "
            "Load bcftools, add it to PATH, or set BCFTOOLS=/path/to/bcftools before running "
            "bin/cbicall test --wes-cohort-bash-sharded or --all."
        )
    return bcftools


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


def _check_text_expectations(run_dir: Path, expectations: List[dict]) -> None:
    for item in expectations:
        rel_path = item["path"]
        path = _work_path(run_dir, rel_path)
        if not path.is_file():
            raise IntegrationTestError(f"Text expectation target does not exist: {rel_path}")
        text = "\n".join(_read_text_lines(path))
        if "contains" in item and str(item["contains"]) not in text:
            raise IntegrationTestError(
                f"Text expectation failed for {rel_path}: expected to contain {item['contains']!r}"
            )
        if "startswith" in item and not text.startswith(str(item["startswith"])):
            raise IntegrationTestError(
                f"Text expectation failed for {rel_path}: expected prefix {item['startswith']!r}"
            )
        if "endswith" in item and not text.endswith(str(item["endswith"])):
            raise IntegrationTestError(
                f"Text expectation failed for {rel_path}: expected suffix {item['endswith']!r}"
            )
        if "matches" in item and re.search(str(item["matches"]), text, flags=re.MULTILINE) is None:
            raise IntegrationTestError(
                f"Text expectation failed for {rel_path}: expected regex {item['matches']!r}"
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


def _load_report(run_dir: Path) -> dict:
    report_path = run_dir / "run-report.json"
    if not report_path.is_file():
        raise IntegrationTestError(f"Missing run report: {report_path}")
    return json.loads(report_path.read_text(encoding="utf-8"))


def _report_value(report: dict, path: str):
    current = report
    for part in path.split("."):
        if isinstance(current, dict):
            current = current.get(part)
        else:
            return None
    return current


def _vcf_hash_from_report(report: dict) -> Optional[dict]:
    reports = _report_value(report, "outputs.vcf_hash_reports") or []
    if not isinstance(reports, list):
        return None
    usable = [item for item in reports if isinstance(item, dict) and item.get("normalized_sha256")]
    if not usable:
        return None
    return sorted(usable, key=lambda item: str(item.get("file") or item.get("path") or ""))[0]


def _vcf_hash_from_contract(run_dir: Path, contract: dict) -> Optional[dict]:
    for item in contract.get("hashes", []):
        if item.get("type") != "normalized_vcf":
            continue
        path = _work_path(run_dir, item["path"])
        if not path.is_file():
            continue
        observed, records = _normalized_text_hash(path, item.get("pattern", "^#"))
        return {
            "file": path.name,
            "path": str(path),
            "normalized_sha256": observed,
            "normalized_records": str(records),
            "source": "integration_contract",
        }
    return None


def _release_vcf_hash(run_dir: Path, contract: dict) -> dict:
    report = _load_report(run_dir)
    vcf_hash = _vcf_hash_from_report(report) or _vcf_hash_from_contract(run_dir, contract)
    if not vcf_hash:
        raise IntegrationTestError(f"No normalized VCF hash found for release comparison: {run_dir}")
    return vcf_hash


def _release_metadata(report: dict) -> dict:
    return {
        "pipeline": _report_value(report, "workflow.pipeline"),
        "mode": _report_value(report, "workflow.mode"),
        "genome": _report_value(report, "run.display_genome"),
        "software_stack": _report_value(report, "workflow.software_stack"),
        "resource": _report_value(report, "resources.bundle.key"),
    }


def _check_release_metadata(baseline_run: Path, candidate_run: Path) -> None:
    baseline = _release_metadata(_load_report(baseline_run))
    candidate = _release_metadata(_load_report(candidate_run))
    mismatches = [key for key, value in baseline.items() if candidate.get(key) != value]
    if mismatches:
        details = ", ".join(
            f"{key}: {baseline.get(key)!r} != {candidate.get(key)!r}" for key in mismatches
        )
        raise IntegrationTestError(f"Release metadata mismatch: {details}")


def _backend_is_available(selection: TestSelection) -> Tuple[bool, str]:
    if selection.key == "wes-cromwell":
        if os.environ.get("CROMWELL_JAR") or shutil.which("cromwell"):
            return True, ""
        return False, "CROMWELL_JAR or cromwell executable not found"
    if selection.backend_executable and shutil.which(selection.backend_executable) is None:
        return False, f"backend executable {selection.backend_executable} not found"
    return True, ""


def _short_hash(value: Optional[str]) -> str:
    if not value:
        return "(undef)"
    return value if len(value) <= 16 else f"{value[:12]}...{value[-8:]}"


def _release_hash_detail(hash_info: Optional[dict]) -> str:
    if not hash_info:
        return "(hash not available)"
    records = hash_info.get("normalized_records")
    record_text = f" | {records} records" if records is not None else ""
    return f"{_short_hash(hash_info.get('normalized_sha256'))}{record_text}"


def _status_marker(status: str) -> str:
    normalized = status.lower()
    if normalized in {"passed", "same final vcf"}:
        return "[PASS]"
    if normalized == "skipped":
        return "[SKIP]"
    if normalized == "failed":
        return "[FAIL]"
    return "[INFO]"


def _overall_marker(exit_code: int) -> str:
    return "[PASS]" if exit_code == 0 else "[FAIL]"


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

    try:
        _check_text_expectations(run_dir, contract.get("text_expectations", []))
    except IntegrationTestError as exc:
        errors.append(str(exc))

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


def _write_setup_files_from_runs(workdir: Path, contract: dict, setup_runs: Dict[str, Path]) -> List[Path]:
    paths: List[Path] = []
    for item in contract.get("setup_files_from_runs", []):
        run_key = item["run"]
        if run_key not in setup_runs:
            raise IntegrationTestError(f"Setup run was not executed: {run_key}")
        run_dir = setup_runs[run_key]
        source = _work_path(run_dir, item["source"])
        if not source.is_file():
            raise IntegrationTestError(f"Setup source does not exist: {source}")
        if item.get("source_index"):
            source_index = _work_path(run_dir, item["source_index"])
            if not source_index.is_file():
                raise IntegrationTestError(f"Setup source index does not exist: {source_index}")

        if item.get("sample_map"):
            path = _work_path(workdir, item["sample_map"])
            path.parent.mkdir(parents=True, exist_ok=True)
            lines = [f"{sample}\t{source}" for sample in item.get("samples", [])]
            if not lines:
                raise IntegrationTestError(f"Setup sample map requires at least one sample: {path}")
            path.write_text("\n".join(lines) + "\n", encoding="utf-8")
            paths.append(path)
    return paths


def _cleanup_setup_files(paths: List[Path], workdir: Path) -> None:
    for path in paths:
        path.unlink(missing_ok=True)
    parents = sorted({path.parent for path in paths}, key=lambda path: len(path.parts), reverse=True)
    for parent in parents:
        if parent == workdir:
            continue
        try:
            parent.rmdir()
        except OSError:
            pass


def _run_setup_runs(
    *,
    project_root: Path,
    contract: dict,
    threads: int,
    runtime_profile: str,
    keep_external_work: bool,
    stack: List[str],
) -> Dict[str, Path]:
    setup_runs: Dict[str, Path] = {}
    for key in contract.get("setup_runs", []):
        if key in stack:
            raise IntegrationTestError(f"Recursive integration setup dependency: {' -> '.join(stack + [key])}")
        if key not in TESTS:
            raise IntegrationTestError(f"Unknown integration setup dependency: {key}")
        print(f"Preparing setup run: {TESTS[key].label}")
        label, status, detail = _run_one(
            project_root=project_root,
            selection=TESTS[key],
            threads=threads,
            runtime_profile=runtime_profile,
            skip_missing_optional=False,
            keep_external_work=keep_external_work,
            stack=stack + [key],
        )
        if status != "passed":
            raise IntegrationTestError(f"Setup run did not pass: {label} ({status})")
        setup_runs[key] = Path(detail)
    return setup_runs


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


def _run_staged_finalize(
    *,
    project_root: Path,
    selection: TestSelection,
    workdir: Path,
    shard_run_dir: Path,
    contract: dict,
    threads: int,
    runtime_profile: str,
    keep_external_work: bool,
) -> Path:
    staged = contract.get("staged_finalize") or {}
    if not isinstance(staged, dict):
        raise IntegrationTestError("staged_finalize must be a mapping.")

    env = _contract_env(project_root, contract)
    bcftools = _resolve_bcftools(project_root, env)
    shard_vcfs = staged.get("shard_vcfs") or []
    if not shard_vcfs:
        raise IntegrationTestError("staged_finalize.shard_vcfs must list at least one raw shard VCF.")

    raw_vcfs = [_work_path(shard_run_dir, rel_path) for rel_path in shard_vcfs]
    for raw_vcf in raw_vcfs:
        if not raw_vcf.is_file():
            raise IntegrationTestError(f"Staged shard raw VCF does not exist: {raw_vcf}")
        if not Path(str(raw_vcf) + ".tbi").is_file():
            raise IntegrationTestError(f"Staged shard raw VCF index does not exist: {raw_vcf}.tbi")

    raw_vcfs_list = _work_path(workdir, staged.get("raw_vcfs_list", "raw_vcfs.list"))
    raw_vcfs_list.write_text("".join(f"{path}\n" for path in raw_vcfs), encoding="utf-8")
    gathered_vcf = _work_path(workdir, staged.get("gathered_vcf", "cohort.gathered.gv.raw.vcf.gz"))
    gathered_vcf.unlink(missing_ok=True)
    Path(str(gathered_vcf) + ".tbi").unlink(missing_ok=True)

    print("Gathering staged raw cohort VCFs:")
    print(f"  BCFTOOLS    : {bcftools}")
    print(f"  VCF list    : {raw_vcfs_list}")
    print(f"  Gathered VCF: {gathered_vcf}")

    for cmd in (
        [bcftools, "concat", "-f", str(raw_vcfs_list), "-Oz", "-o", str(gathered_vcf)],
        [bcftools, "index", "-t", str(gathered_vcf)],
    ):
        proc = subprocess.run(
            cmd,
            cwd=str(workdir),
            env=env,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            check=False,
        )
        if proc.returncode != 0:
            raise IntegrationTestError(
                "Staged cohort gather command failed with return code "
                f"{proc.returncode}: {' '.join(cmd)}\n{proc.stdout}"
            )

    finalize = staged.get("finalize") or {}
    if not isinstance(finalize, dict) or "parameters" not in finalize:
        raise IntegrationTestError("staged_finalize.finalize.parameters is required.")
    finalize_parameters = dict(finalize["parameters"] or {})
    finalize_parameters.setdefault("input_vcf", str(gathered_vcf))
    param_file = _write_inline_parameters(workdir, f"{selection.key}-finalize", finalize_parameters)

    launcher_fd, launcher_name = tempfile.mkstemp(
        prefix=f"cbicall-test-{selection.key}-finalize.",
        suffix=".log",
        dir=workdir,
        text=True,
    )
    os.close(launcher_fd)
    launcher_log = Path(launcher_name)

    base_dir = _work_path(workdir, finalize.get("base_dir", "."))
    run_glob = finalize.get("run_glob") or contract.get("run", {}).get("run_glob")
    if not run_glob:
        raise IntegrationTestError("staged_finalize.finalize.run_glob is required.")
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
    print(f"TEST: {selection.label} Finalize")
    print("========================================")
    print("Running staged cohort finalize integration test...")
    proc = subprocess.run(
        cmd,
        cwd=str(workdir),
        env=env,
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
                f"{selection.label} finalize cbicall command failed with return code {proc.returncode}.\n"
                f"Launcher log: {launcher_log}\n{tail}"
            )
        if run_dir is None:
            raise IntegrationTestError(f"{selection.label} finalize finished but no run directory was found.")

        finalize_contract = {
            "_contract_path": contract["_contract_path"],
            "workflow_log": finalize.get("workflow_log", contract.get("workflow_log", "workflow.log")),
            "required_files": finalize.get("required_files", []),
            "required_globs": finalize.get("required_globs", []),
            "text_expectations": finalize.get("text_expectations", []),
            "json_expectations": finalize.get("json_expectations", []),
            "hashes": finalize.get("hashes", []),
            "cleanup_paths": finalize.get("cleanup_paths", []),
        }
        _print_artifacts(f"{selection.label} Finalize", run_dir, finalize_contract, launcher_log)
        print("Validating staged finalize contract:")
        validate_contract(run_dir, finalize_contract)
        print("  OK finalize required files")
        print("  OK finalize run-report fields")
        _cleanup(run_dir, finalize_contract, keep_external_work)
        return run_dir
    finally:
        param_file.unlink(missing_ok=True)


def _run_one(
    *,
    project_root: Path,
    selection: TestSelection,
    threads: int,
    runtime_profile: str,
    skip_missing_optional: bool,
    keep_external_work: bool,
    stack: Optional[List[str]] = None,
) -> Tuple[str, str, str]:
    stack = stack or [selection.key]
    if selection.key == "wes-cromwell" and not (os.environ.get("CROMWELL_JAR") or shutil.which("cromwell")):
        detail = "CROMWELL_JAR or cromwell executable not found"
        if skip_missing_optional and selection.optional_in_all:
            print("SKIP: WES Cromwell requires CROMWELL_JAR or cromwell on PATH.")
            return selection.label, "skipped", detail
        raise IntegrationTestError("WES Cromwell requires CROMWELL_JAR or cromwell on PATH.")

    if selection.backend_executable and shutil.which(selection.backend_executable) is None:
        detail = f"backend executable {selection.backend_executable} not found"
        if skip_missing_optional and selection.optional_in_all:
            print(f"SKIP: {selection.label} requires backend executable {selection.backend_executable} on PATH.")
            return selection.label, "skipped", detail
        raise IntegrationTestError(f"{selection.label} requires backend executable {selection.backend_executable} on PATH.")

    contract = load_contract(project_root, selection)
    workdir = _project_path(project_root, contract.get("workdir", "examples/input"))
    base_dir = _work_path(workdir, contract.get("run", {}).get("base_dir", "."))
    run_glob = contract.get("run", {}).get("run_glob")
    if not run_glob:
        raise IntegrationTestError("Contract run section must define run_glob")

    setup_files: List[Path] = []
    temp_param: Optional[Path] = None
    setup_runs = _run_setup_runs(
        project_root=project_root,
        contract=contract,
        threads=threads,
        runtime_profile=runtime_profile,
        keep_external_work=keep_external_work,
        stack=stack,
    )
    setup_files.extend(_write_setup_files_from_runs(workdir, contract, setup_runs))
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
    for note in contract.get("notes", []):
        print(f"Note: {note}")

    proc = subprocess.run(
        cmd,
        cwd=str(workdir),
        env=_contract_env(project_root, contract),
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
        detail_run_dir = run_dir
        if contract.get("staged_finalize"):
            detail_run_dir = _run_staged_finalize(
                project_root=project_root,
                selection=selection,
                workdir=workdir,
                shard_run_dir=run_dir,
                contract=contract,
                threads=threads,
                runtime_profile=runtime_profile,
                keep_external_work=keep_external_work,
            )
        _cleanup(run_dir, contract, keep_external_work)
        print(f"SUCCESS: {selection.label} integration contract passed.")
        return selection.label, "passed", str(detail_run_dir)
    finally:
        if temp_param:
            temp_param.unlink(missing_ok=True)
        _cleanup_setup_files(setup_files, workdir)


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
    print(f"{'Result':<8} {'Test':<24} Details")
    print(f"{'-' * 8} {'-' * 24} {'-' * 40}")
    for label, status, detail in summary:
        details = detail or ""
        print(f"{_status_marker(status):<8} {label:<24} {details}")
    print("========================================")
    print("All requested tests finished.")
    print(f"Status: {_overall_marker(overall_status)} {'PASSED' if overall_status == 0 else 'FAILED'}")
    print(f"Exit code: {overall_status}")
    print("========================================")
    return overall_status


def run_release_equivalence_test(
    *,
    project_root: Path,
    threads: int,
    runtime_profile: str,
) -> int:
    baseline_selection = TESTS["wes-bash"]
    comparator_selections = [TESTS["wes-snakemake"], TESTS["wes-nextflow"], TESTS["wes-cromwell"]]
    results: List[Tuple[str, str, str]] = []
    comparison_hashes: Dict[str, dict] = {}
    overall_status = 0
    passed_comparisons = 0

    print("========================================")
    print("Backend equivalence")
    print("========================================")
    print("Running native WES backends and comparing normalized final VCF content.")

    try:
        baseline_label, baseline_status, baseline_detail = _run_one(
            project_root=project_root,
            selection=baseline_selection,
            threads=threads,
            runtime_profile=runtime_profile,
            skip_missing_optional=False,
            keep_external_work=False,
        )
        results.append((baseline_label, baseline_status, baseline_detail))
        baseline_run = Path(baseline_detail)
        baseline_contract = load_contract(project_root, baseline_selection)
        baseline_hash = _release_vcf_hash(baseline_run, baseline_contract)
        comparison_hashes[baseline_label] = baseline_hash
    except IntegrationTestError as exc:
        overall_status = 1
        results.append((baseline_selection.label, "failed", str(exc).splitlines()[0]))
        baseline_run = None
        baseline_hash = None

    if baseline_run is not None and baseline_hash is not None:
        for selection in comparator_selections:
            available, detail = _backend_is_available(selection)
            if not available:
                results.append((selection.label, "skipped", detail))
                continue

            try:
                label, status, run_detail = _run_one(
                    project_root=project_root,
                    selection=selection,
                    threads=threads,
                    runtime_profile=runtime_profile,
                    skip_missing_optional=False,
                    keep_external_work=False,
                )
                candidate_run = Path(run_detail)
                candidate_contract = load_contract(project_root, selection)
                _check_release_metadata(baseline_run, candidate_run)
                candidate_hash = _release_vcf_hash(candidate_run, candidate_contract)
                if candidate_hash.get("normalized_sha256") != baseline_hash.get("normalized_sha256"):
                    raise IntegrationTestError(
                        "normalized VCF hash differs from Bash baseline: "
                        f"{candidate_hash.get('normalized_sha256')} != {baseline_hash.get('normalized_sha256')}"
                    )
                baseline_records = baseline_hash.get("normalized_records")
                candidate_records = candidate_hash.get("normalized_records")
                if baseline_records is not None and candidate_records is not None and str(candidate_records) != str(baseline_records):
                    raise IntegrationTestError(
                        f"normalized VCF record count differs from Bash baseline: {candidate_records} != {baseline_records}"
                    )
                passed_comparisons += 1
                comparison_hashes[label] = candidate_hash
                results.append((label, "same final VCF", str(candidate_run)))
            except IntegrationTestError as exc:
                overall_status = 1
                results.append((selection.label, "failed", str(exc).splitlines()[0]))

    if baseline_run is not None and passed_comparisons == 0:
        overall_status = 1
        results.append(("Backend equivalence", "failed", "no non-Bash backend was available for comparison"))

    print("========================================")
    print("Backend equivalence summary")
    print("========================================")
    print("Baseline")
    for label, status, detail in results:
        if label != baseline_selection.label:
            continue
        hash_info = comparison_hashes.get(label)
        suffix = f" ({detail})" if detail else ""
        marker = _status_marker(status)
        if hash_info:
            print(f"  {marker:<8} {label:<13} => {status} | {_release_hash_detail(hash_info)}{suffix}")
        else:
            print(f"  {marker:<8} {label:<13} => {status}{suffix}")
    print()
    print("Backend equivalence")
    for label, status, detail in results:
        if label == baseline_selection.label:
            continue
        hash_info = comparison_hashes.get(label)
        suffix = f" ({detail})" if detail else ""
        marker = _status_marker(status)
        if hash_info:
            print(f"  {marker:<8} {label:<13} => {status} | {_release_hash_detail(hash_info)}{suffix}")
        else:
            print(f"  {marker:<8} {label:<13} => {status}{suffix}")
    print()
    print(f"Compared non-Bash backends: {passed_comparisons}")
    print(f"Status: {_overall_marker(overall_status)} {'PASSED' if overall_status == 0 else 'FAILED'}")
    print(f"Exit code: {overall_status}")
    print("========================================")
    return overall_status


def selected_tests_from_args(args: argparse.Namespace) -> List[TestSelection]:
    selected: List[TestSelection] = []
    if args.all or args.wes_bash:
        selected.append(TESTS["wes-bash"])
    # In --all, run MIT before cohort tests create additional WES setup dirs.
    # The MIT Bash workflow discovers the upstream WES BAM from the sample tree.
    if args.all or args.mit_bash:
        selected.append(TESTS["mit-bash"])
    if args.all or args.wes_cohort_bash:
        selected.append(TESTS["wes-cohort-bash"])
    if args.all or args.wes_cohort_bash_sharded:
        selected.append(TESTS["wes-cohort-bash-sharded"])
    if args.wes_bash_gatk35:
        selected.append(TESTS["wes-bash-gatk35"])
    if args.all or args.wes_snakemake:
        selected.append(TESTS["wes-snakemake"])
    if args.all or args.wes_nextflow:
        selected.append(TESTS["wes-nextflow"])
    if args.all or args.wes_cromwell:
        selected.append(TESTS["wes-cromwell"])
    if args.nf_core_demo:
        selected.append(TESTS["nf-core-demo"])
    if args.nf_core_sarek:
        selected.append(TESTS["nf-core-sarek"])
    return selected
