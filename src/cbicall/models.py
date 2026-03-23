from dataclasses import asdict, dataclass, field
from typing import Any, Dict, Mapping, Optional


@dataclass(frozen=True)
class InputsSpec:
    input_dir: Optional[str] = None
    sample_map: Optional[str] = None

    @classmethod
    def from_mapping(cls, data: Mapping[str, Any]) -> "InputsSpec":
        return cls(
            input_dir=data.get("input_dir"),
            sample_map=data.get("sample_map"),
        )

    def to_dict(self) -> Dict[str, Optional[str]]:
        return asdict(self)


@dataclass(frozen=True)
class WorkflowSpec:
    engine: str
    pipeline: str
    mode: str
    gatk_version: str
    entrypoint: Optional[str]
    config_file: Optional[str] = None
    helpers: Dict[str, str] = field(default_factory=dict)

    @classmethod
    def from_mapping(cls, data: Mapping[str, Any]) -> "WorkflowSpec":
        return cls(
            engine=str(data["engine"]),
            pipeline=str(data["pipeline"]),
            mode=str(data["mode"]),
            gatk_version=str(data["gatk_version"]),
            entrypoint=data.get("entrypoint"),
            config_file=data.get("config_file"),
            helpers=dict(data.get("helpers", {})),
        )

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class RunSettings:
    project_dir: str
    run_id: str
    threads: int
    debug: bool
    genome: Optional[str]
    cleanup_bam: bool
    inputs: InputsSpec
    workflow: WorkflowSpec
    workflow_rule: Optional[str] = None
    allow_partial_run: bool = False
    run_mode: str = "full"

    @classmethod
    def from_mapping(cls, data: Mapping[str, Any]) -> "RunSettings":
        project_dir = data["project_dir"]
        run_id = data["run_id"] if "run_id" in data else data.get("id")
        return cls(
            project_dir=str(project_dir),
            run_id=str(run_id),
            threads=int(data["threads"]),
            debug=bool(data["debug"]),
            genome=data.get("genome"),
            cleanup_bam=bool(data.get("cleanup_bam", False)),
            workflow_rule=data.get("workflow_rule"),
            allow_partial_run=bool(data.get("allow_partial_run", False)),
            run_mode=str(data.get("run_mode", "full")),
            inputs=InputsSpec.from_mapping(data.get("inputs", {})),
            workflow=WorkflowSpec.from_mapping(data["workflow"]),
        )

    def to_dict(self) -> Dict[str, Any]:
        return {
            "project_dir": self.project_dir,
            "run_id": self.run_id,
            "threads": self.threads,
            "debug": self.debug,
            "genome": self.genome,
            "cleanup_bam": self.cleanup_bam,
            "workflow_rule": self.workflow_rule,
            "allow_partial_run": self.allow_partial_run,
            "run_mode": self.run_mode,
            "inputs": self.inputs.to_dict(),
            "workflow": self.workflow.to_dict(),
        }


@dataclass(frozen=True)
class ResolvedConfig:
    user: str
    workflow_engine: str
    genome: Optional[str]
    pipeline: str
    mode: str
    gatk_version: str
    inputs: InputsSpec
    workflow: WorkflowSpec
    run_id: str
    date: str
    project_dir: str
    output_basename: Optional[str]
    hostname: str
    host_threads: int
    host_threads_minus_one: int
    compression_cmd: str
    workflow_rule: Optional[str] = None
    allow_partial_run: bool = False
    run_mode: str = "full"
    capture_label: Optional[str] = None
    arch: Optional[str] = None
    version: Optional[str] = None

    @classmethod
    def from_mapping(cls, data: Mapping[str, Any]) -> "ResolvedConfig":
        workflow_engine = data["workflow_engine"] if "workflow_engine" in data else data["workflow"]["engine"]
        pipeline = data["pipeline"] if "pipeline" in data else data["workflow"]["pipeline"]
        mode = data["mode"] if "mode" in data else data["workflow"]["mode"]
        gatk_version = data["gatk_version"] if "gatk_version" in data else data["workflow"]["gatk_version"]
        run_id = data["run_id"] if "run_id" in data else data["id"]
        project_dir = data["project_dir"]
        return cls(
            user=str(data.get("user", "")),
            workflow_engine=str(workflow_engine),
            genome=data.get("genome"),
            pipeline=str(pipeline),
            mode=str(mode),
            gatk_version=str(gatk_version),
            inputs=InputsSpec.from_mapping(data.get("inputs", {})),
            workflow=WorkflowSpec.from_mapping(data["workflow"]),
            run_id=str(run_id),
            date=str(data.get("date", "")),
            project_dir=str(project_dir),
            output_basename=data.get("output_basename"),
            hostname=str(data.get("hostname", "")),
            host_threads=int(data.get("host_threads", 0)),
            host_threads_minus_one=int(data.get("host_threads_minus_one", 0)),
            compression_cmd=str(data.get("compression_cmd", "")),
            workflow_rule=data.get("workflow_rule"),
            allow_partial_run=bool(data.get("allow_partial_run", False)),
            run_mode=str(data.get("run_mode", "full")),
            capture_label=data.get("capture_label"),
            arch=data.get("arch"),
            version=data.get("version"),
        )

    def to_dict(self) -> Dict[str, Any]:
        return {
            "user": self.user,
            "workflow_engine": self.workflow_engine,
            "genome": self.genome,
            "pipeline": self.pipeline,
            "mode": self.mode,
            "gatk_version": self.gatk_version,
            "inputs": self.inputs.to_dict(),
            "workflow": self.workflow.to_dict(),
            "run_id": self.run_id,
            "date": self.date,
            "project_dir": self.project_dir,
            "output_basename": self.output_basename,
            "hostname": self.hostname,
            "host_threads": self.host_threads,
            "host_threads_minus_one": self.host_threads_minus_one,
            "compression_cmd": self.compression_cmd,
            "workflow_rule": self.workflow_rule,
            "allow_partial_run": self.allow_partial_run,
            "run_mode": self.run_mode,
            "capture_label": self.capture_label,
            "arch": self.arch,
            "version": self.version,
        }
