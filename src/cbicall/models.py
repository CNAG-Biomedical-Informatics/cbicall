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
    pipeline_version: str
    entrypoint: Optional[str]
    config_file: Optional[str] = None
    helpers: Dict[str, str] = field(default_factory=dict)
    profiles: Dict[str, Dict[str, str]] = field(default_factory=dict)
    metadata: Dict[str, Any] = field(default_factory=dict)

    @classmethod
    def from_mapping(cls, data: Mapping[str, Any]) -> "WorkflowSpec":
        return cls(
            engine=str(data["engine"]),
            pipeline=str(data["pipeline"]),
            mode=str(data["mode"]),
            gatk_version=str(data["gatk_version"]),
            pipeline_version=str(data.get("pipeline_version", "legacy")),
            entrypoint=data.get("entrypoint"),
            config_file=data.get("config_file"),
            helpers=dict(data.get("helpers", {})),
            profiles={str(k): dict(v) for k, v in data.get("profiles", {}).items()},
            metadata=dict(data.get("metadata", {})),
        )

    def to_dict(self) -> Dict[str, Any]:
        return {
            "engine": self.engine,
            "pipeline": self.pipeline,
            "mode": self.mode,
            "gatk_version": self.gatk_version,
            "pipeline_version": self.pipeline_version,
            "entrypoint": self.entrypoint,
            "config_file": self.config_file,
            "helpers": dict(self.helpers),
            "metadata": dict(self.metadata),
        }


@dataclass(frozen=True)
class RunSettings:
    project_dir: str
    run_id: str
    threads: int
    debug: bool
    profile: str
    genome: Optional[str]
    cleanup_bam: bool
    inputs: InputsSpec
    workflow: WorkflowSpec
    snakemake_parameters: Dict[str, Any] = field(default_factory=dict)
    nextflow_parameters: Dict[str, Any] = field(default_factory=dict)
    nfcore_profile: Optional[str] = None
    nfcore_parameters: Dict[str, Any] = field(default_factory=dict)
    nfcore_singularity_cache_dir: Optional[str] = None
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
            profile=str(data.get("profile", "local")),
            snakemake_parameters=dict(data.get("snakemake_parameters", {})),
            nextflow_parameters=dict(data.get("nextflow_parameters", {})),
            nfcore_profile=data.get("nfcore_profile"),
            nfcore_parameters=dict(data.get("nfcore_parameters", {})),
            nfcore_singularity_cache_dir=data.get("nfcore_singularity_cache_dir"),
            genome=data.get("genome"),
            cleanup_bam=bool(data.get("cleanup_bam", False)),
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
            "profile": self.profile,
            "snakemake_parameters": dict(self.snakemake_parameters),
            "nextflow_parameters": dict(self.nextflow_parameters),
            "nfcore_profile": self.nfcore_profile,
            "nfcore_parameters": dict(self.nfcore_parameters),
            "nfcore_singularity_cache_dir": self.nfcore_singularity_cache_dir,
            "genome": self.genome,
            "cleanup_bam": self.cleanup_bam,
            "run_mode": self.run_mode,
            "inputs": self.inputs.to_dict(),
            "workflow": self.workflow.to_dict(),
        }


@dataclass(frozen=True)
class ResolvedConfig:
    user: str
    workflow_engine: str
    profile: str
    genome: Optional[str]
    pipeline: str
    mode: str
    gatk_version: str
    pipeline_version: str
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
    snakemake_parameters: Dict[str, Any] = field(default_factory=dict)
    nextflow_parameters: Dict[str, Any] = field(default_factory=dict)
    nfcore_profile: Optional[str] = None
    nfcore_parameters: Dict[str, Any] = field(default_factory=dict)
    nfcore_singularity_cache_dir: Optional[str] = None
    run_mode: str = "full"
    capture_label: Optional[str] = None
    arch: Optional[str] = None
    version: Optional[str] = None
    resources: Dict[str, Any] = field(default_factory=dict)

    @classmethod
    def from_mapping(cls, data: Mapping[str, Any]) -> "ResolvedConfig":
        workflow_engine = data["workflow_engine"] if "workflow_engine" in data else data["workflow"]["engine"]
        pipeline = data["pipeline"] if "pipeline" in data else data["workflow"]["pipeline"]
        mode = data["mode"] if "mode" in data else data["workflow"]["mode"]
        gatk_version = data["gatk_version"] if "gatk_version" in data else data["workflow"]["gatk_version"]
        pipeline_version = data.get("pipeline_version")
        if pipeline_version is None and "workflow" in data:
            pipeline_version = data["workflow"].get("pipeline_version")
        run_id = data["run_id"] if "run_id" in data else data["id"]
        project_dir = data["project_dir"]
        return cls(
            user=str(data.get("user", "")),
            workflow_engine=str(workflow_engine),
            profile=str(data.get("profile", "local")),
            snakemake_parameters=dict(data.get("snakemake_parameters", {})),
            nextflow_parameters=dict(data.get("nextflow_parameters", {})),
            nfcore_profile=data.get("nfcore_profile"),
            nfcore_parameters=dict(data.get("nfcore_parameters", {})),
            nfcore_singularity_cache_dir=data.get("nfcore_singularity_cache_dir"),
            genome=data.get("genome"),
            pipeline=str(pipeline),
            mode=str(mode),
            gatk_version=str(gatk_version),
            pipeline_version=str(pipeline_version) if pipeline_version is not None else "legacy",
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
            run_mode=str(data.get("run_mode", "full")),
            capture_label=data.get("capture_label"),
            arch=data.get("arch"),
            version=data.get("version"),
            resources=dict(data.get("resources", {})),
        )

    def to_dict(self) -> Dict[str, Any]:
        return {
            "user": self.user,
            "workflow_engine": self.workflow_engine,
            "profile": self.profile,
            "snakemake_parameters": dict(self.snakemake_parameters),
            "nextflow_parameters": dict(self.nextflow_parameters),
            "nfcore_profile": self.nfcore_profile,
            "nfcore_parameters": dict(self.nfcore_parameters),
            "nfcore_singularity_cache_dir": self.nfcore_singularity_cache_dir,
            "genome": self.genome,
            "pipeline": self.pipeline,
            "mode": self.mode,
            "gatk_version": self.gatk_version,
            "pipeline_version": self.pipeline_version,
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
            "run_mode": self.run_mode,
            "capture_label": self.capture_label,
            "arch": self.arch,
            "version": self.version,
            "resources": self.resources,
        }
