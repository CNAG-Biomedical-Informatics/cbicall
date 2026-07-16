"""Read-only inspection of workflow runtimes and configured Java executables."""

import hashlib
import os
import platform
import re
import shutil
import subprocess
import sys
from pathlib import Path
from typing import List, Optional

import yaml


def parse_version_text(text: str) -> Optional[str]:
    """Return the first dotted version string found in command output."""
    match = re.search(r"\b\d+(?:\.\d+)+(?:[-+._A-Za-z0-9]*)?\b", text)
    return match.group(0) if match else None


def command_version(command: str, args: List[str]) -> dict:
    """Run a bounded version probe and return structured status metadata."""
    path = shutil.which(command)
    report = {
        "name": command,
        "path": path,
        "command": " ".join([command, *args]),
        "status": "not_found" if path is None else "unknown",
        "version": None,
    }
    if path is None:
        return report

    try:
        completed = subprocess.run(
            [path, *args],
            capture_output=True,
            text=True,
            timeout=10,
            check=False,
        )
    except Exception as exc:
        report["status"] = "error"
        report["error"] = str(exc)
        return report

    output = "\n".join(part for part in [completed.stdout, completed.stderr] if part).strip()
    first_line = next((line.strip() for line in output.splitlines() if line.strip()), "")
    report["status"] = "ok" if completed.returncode == 0 else "nonzero_exit"
    report["returncode"] = completed.returncode
    report["version"] = parse_version_text(output)
    if first_line:
        report["detail"] = first_line[:200]
    return report


def architecture_key() -> str:
    """Return the architecture key used by native workflow configuration files."""
    machine = platform.machine().lower()
    if machine in {"x86_64", "amd64"}:
        return "amd64"
    if machine in {"aarch64", "arm64"}:
        return "aarch64"
    return machine


def source_bash_env_variable(env_file: str, variable: str) -> dict:
    """Read one variable after sourcing a registered Bash environment file."""
    report = {
        "source": str(env_file),
        "variable": variable,
        "status": "unknown",
        "path": None,
    }
    path = Path(str(env_file))
    if not path.is_file():
        report["status"] = "source_missing"
        return report
    try:
        completed = subprocess.run(
            [
                "bash",
                "-c",
                'source "$1" >/dev/null 2>&1; printf "%s" "${!2:-}"',
                "bash",
                str(path),
                variable,
            ],
            capture_output=True,
            text=True,
            timeout=10,
            check=False,
        )
    except Exception as exc:
        report["status"] = "error"
        report["error"] = str(exc)
        return report

    value = completed.stdout.strip()
    report["returncode"] = completed.returncode
    if completed.returncode != 0:
        report["status"] = "nonzero_exit"
        if completed.stderr:
            report["detail"] = completed.stderr.strip().splitlines()[0][:200]
        return report
    if not value:
        report["status"] = "not_set"
        return report
    report["status"] = "ok"
    report["path"] = value
    return report


def configured_native_java_report(workflow) -> Optional[dict]:
    """Inspect the Java executable selected by a native workflow definition."""
    if workflow is None or workflow.metadata.get("provider") == "nf-core":
        return None

    java_path = None
    source = None
    status = "not_configured"
    arch = architecture_key()

    if workflow.backend == "bash":
        env_file = workflow.helpers.get("env")
        if not env_file:
            return None
        env_report = source_bash_env_variable(env_file, "JAVA8")
        java_path = env_report.get("path")
        source = f"{env_file}:JAVA8"
        status = env_report.get("status") or status
        if not java_path:
            env_report["name"] = "configured_java"
            return env_report
    elif workflow.backend in {"snakemake", "nextflow"}:
        config_file = workflow.config_file
        if not config_file:
            return None
        config_path = Path(str(config_file))
        source = f"{config_file}:java.{arch}"
        if not config_path.is_file():
            return {
                "name": "configured_java",
                "source": source,
                "status": "source_missing",
                "path": None,
                "version": None,
            }
        try:
            config = yaml.safe_load(config_path.read_text(encoding="utf-8")) or {}
        except Exception as exc:
            return {
                "name": "configured_java",
                "source": source,
                "status": "error",
                "path": None,
                "version": None,
                "error": str(exc),
            }
        java_cfg = config.get("java") or {}
        java_path = java_cfg.get(arch)
        if not java_path:
            return {
                "name": "configured_java",
                "source": source,
                "status": "not_set",
                "path": None,
                "version": None,
            }
        status = "ok"
    else:
        return None

    version_report = command_version(str(java_path), ["-version"])
    version_report["name"] = "configured_java"
    version_report["source"] = source
    if status != "ok" and version_report.get("status") == "not_found":
        version_report["status"] = status
    return version_report


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def cromwell_jar_report(jar_path: Path, *, source: str) -> dict:
    """Inspect a Cromwell JAR and the Java launcher used to execute it."""
    java_cmd = os.environ.get("JAVA_CMD", "java")
    report = command_version(java_cmd, ["-jar", str(jar_path), "--version"])
    report["name"] = "cromwell"
    report["launcher"] = "java-jar"
    report["jar"] = str(jar_path)
    report["source"] = source
    if jar_path.is_file():
        report["jar_sha256"] = _sha256_file(jar_path)
    else:
        report["status"] = "not_found"
        report["error"] = f"Cromwell JAR does not exist: {jar_path}"
    return report


def cromwell_runtime_report() -> dict:
    """Inspect the Cromwell runtime selected by the execution layer."""
    jar = os.environ.get("CROMWELL_JAR")
    if jar:
        return cromwell_jar_report(Path(jar), source="CROMWELL_JAR")
    return command_version("cromwell", ["--version"])


def backend_runtime_report(backend: str) -> dict:
    """Inspect one workflow backend using the same command names as execution."""
    commands = {
        "bash": ("bash", ["--version"]),
        "snakemake": ("snakemake", ["--version"]),
        "nextflow": ("nextflow", ["-version"]),
    }
    if backend == "cromwell":
        return cromwell_runtime_report()
    if backend in commands:
        command, args = commands[backend]
        return command_version(command, args)
    return {
        "name": backend,
        "path": None,
        "status": "unsupported",
        "version": None,
    }


def runtime_report(backend: str, workflow=None) -> dict:
    """Build the runtime section stored in a CBIcall run report."""
    payload = {
        "python": {
            "executable": sys.executable,
            "version": sys.version.split()[0],
        },
        "java": command_version("java", ["-version"]),
        "backend": backend_runtime_report(backend),
    }
    configured_java = configured_native_java_report(workflow)
    if configured_java:
        payload["configured_java"] = configured_java
    return payload
