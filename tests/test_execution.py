import json

import pytest

from cbicall import execution
from cbicall.errors import WorkflowExecutionError, WorkflowResolutionError


def test_execution_builds_bash_command_with_flags(tmp_path, monkeypatch):
    recorded = {}

    def fake_run_cmd(cmd, cwd, log_path, env=None, backend=None):
        recorded.update(
            {"cmd": cmd, "cwd": cwd, "log_path": log_path, "env": env, "backend": backend}
        )

    monkeypatch.setattr(execution.WorkflowExecutor, "_run_cmd", staticmethod(fake_run_cmd))

    project_dir = tmp_path / "proj"
    project_dir.mkdir()
    script = str(tmp_path / "wes_single.sh")

    settings = {
        "project_dir": str(project_dir),
        "threads": 8,
        "run_id": "RID123",
        "debug": False,
        "cleanup_bam": True,
        "genome": "b37",
        "qc_coverage_region": "chr22",
        "inputs": {"input_dir": None, "sample_map": "map.txt"},
        "snakemake_parameters": {},
        "nextflow_parameters": {},
        "run_mode": "full",
        "workflow": {
            "backend": "bash",
            "pipeline": "wes",
            "mode": "single",
            "software_stack": "gatk-4.6",
            "registry_version": "v1",
            "entrypoint": script,
            "config_file": None,
            "helpers": {},
        },
    }

    obj = execution.WorkflowExecutor(settings)
    assert obj.run() is True

    assert recorded["cmd"][0] == script
    assert recorded["cmd"][1:3] == ["-t", "8"]
    assert "--pipeline" in recorded["cmd"]
    assert "--cleanup-bam" in recorded["cmd"]
    assert "--sample-map" in recorded["cmd"]
    assert recorded["env"]["GENOME"] == "b37"
    assert recorded["env"]["CBICALL_COVERAGE_REGION"] == "chr22"
    assert recorded["cwd"] == project_dir
    assert recorded["backend"] == "bash"
    contract = json.loads((project_dir / "cbicall-execution-contract.json").read_text(encoding="utf-8"))
    assert contract["kind"] == "cbicall_execution_contract"
    assert contract["workflow"]["key"] == "bash/wes/single/gatk-4.6/v1"
    assert contract["environment_overrides"] == {"CBICALL_COVERAGE_REGION": "chr22", "GENOME": "b37"}
    assert contract["run"]["qc_coverage_region"] == "chr22"
    assert contract["command"]["argv"] == recorded["cmd"]
    assert "PATH" not in contract["environment_overrides"]
    assert len(contract["fingerprint"]) == 64


def test_execution_passes_resolved_env_file_to_bash(tmp_path, monkeypatch):
    recorded = {}

    def fake_run_cmd(cmd, cwd, log_path, env=None, backend=None):
        recorded.update({"env": env})

    monkeypatch.setattr(execution.WorkflowExecutor, "_run_cmd", staticmethod(fake_run_cmd))

    project_dir = tmp_path / "proj"
    project_dir.mkdir()
    env_file = tmp_path / "cnag-hpc-env.sh"

    settings = {
        "project_dir": str(project_dir),
        "threads": 1,
        "run_id": "RIDENV",
        "debug": False,
        "genome": "b37",
        "inputs": {"input_dir": None, "sample_map": None},
        "snakemake_parameters": {},
        "nextflow_parameters": {},
        "run_mode": "full",
        "workflow": {
            "backend": "bash",
            "pipeline": "wes",
            "mode": "single",
            "software_stack": "gatk-4.6",
            "registry_version": "v1",
            "entrypoint": str(tmp_path / "wes_single.sh"),
            "config_file": None,
            "helpers": {"env": str(env_file)},
        },
    }

    assert execution.WorkflowExecutor(settings).run() is True
    assert recorded["env"]["CBICALL_ENV_FILE"] == str(env_file)
    assert recorded["env"]["CBICALL_COVERAGE_REGION"] == "chr1"


def test_execution_debug_prints(monkeypatch, tmp_path, capsys):
    def fake_run_cmd(cmd, cwd, log_path, env=None, backend=None):
        return None

    monkeypatch.setattr(execution.WorkflowExecutor, "_run_cmd", staticmethod(fake_run_cmd))

    project_dir = tmp_path / "proj"
    project_dir.mkdir()

    settings = {
        "project_dir": str(project_dir),
        "threads": 2,
        "run_id": "RIDDBG",
        "debug": True,
        "genome": "b37",
        "inputs": {"input_dir": None, "sample_map": None},
        "snakemake_parameters": {},
        "nextflow_parameters": {},
        "run_mode": "full",
        "workflow": {
            "backend": "bash",
            "pipeline": "wes",
            "mode": "single",
            "software_stack": "gatk-3.5",
            "registry_version": "v1",
            "entrypoint": str(tmp_path / "wes_single.sh"),
            "config_file": None,
            "helpers": {},
        },
    }

    assert execution.WorkflowExecutor(settings).run() is True
    out = capsys.readouterr().out
    assert "Log file:" in out
    assert "GENOME=b37" in out


def test_execution_builds_bash_command_gatk35_has_no_extra_flags(tmp_path, monkeypatch):
    recorded = {}

    def fake_run_cmd(cmd, cwd, log_path, env=None, backend=None):
        recorded.update({"cmd": cmd})

    monkeypatch.setattr(execution.WorkflowExecutor, "_run_cmd", staticmethod(fake_run_cmd))

    project_dir = tmp_path / "proj"
    project_dir.mkdir()
    script = str(tmp_path / "wes_single.sh")

    settings = {
        "project_dir": str(project_dir),
        "threads": 4,
        "run_id": "RID999",
        "debug": False,
        "cleanup_bam": True,
        "genome": "b37",
        "inputs": {"input_dir": None, "sample_map": "map.txt"},
        "snakemake_parameters": {},
        "nextflow_parameters": {},
        "run_mode": "full",
        "workflow": {
            "backend": "bash",
            "pipeline": "wes",
            "mode": "single",
            "software_stack": "gatk-3.5",
            "registry_version": "v1",
            "entrypoint": script,
            "config_file": None,
            "helpers": {},
        },
    }

    assert execution.WorkflowExecutor(settings).run() is True
    assert recorded["cmd"] == [script, "-t", "4"]


def test_execution_builds_snakemake_command_and_config(tmp_path, monkeypatch):
    recorded = {}

    def fake_run_cmd(cmd, cwd, log_path, env=None, backend=None):
        recorded.update({"cmd": cmd, "backend": backend})

    monkeypatch.setattr(execution.WorkflowExecutor, "_run_cmd", staticmethod(fake_run_cmd))

    project_dir = tmp_path / "proj"
    project_dir.mkdir()
    snakefile = str(tmp_path / "wes_single.smk")

    settings = {
        "project_dir": str(project_dir),
        "threads": 12,
        "run_id": "RID777",
        "debug": False,
        "genome": "hg38",
        "qc_coverage_region": "chr22",
        "inputs": {"input_dir": None, "sample_map": None},
        "snakemake_parameters": {},
        "nextflow_parameters": {},
        "run_mode": "full",
        "workflow": {
            "backend": "snakemake",
            "pipeline": "wes",
            "mode": "single",
            "software_stack": "gatk-4.6",
            "registry_version": "v1",
            "entrypoint": snakefile,
            "config_file": str(tmp_path / "config.yaml"),
            "helpers": {},
        },
    }

    assert execution.WorkflowExecutor(settings).run() is True
    cmd = recorded["cmd"]
    assert cmd[:3] == ["snakemake", "--forceall", "all"]
    assert "--configfile" in cmd
    assert str(tmp_path / "config.yaml") in cmd
    assert "genome=hg38" in cmd
    assert "qc_coverage_region=chr22" in cmd
    assert "pipeline=wes" in cmd
    assert recorded["backend"] == "snakemake"


def test_execution_builds_snakemake_partial_rule_command(tmp_path, monkeypatch):
    recorded = {}

    def fake_run_cmd(cmd, cwd, log_path, env=None, backend=None):
        recorded.update({"cmd": cmd, "backend": backend})

    monkeypatch.setattr(execution.WorkflowExecutor, "_run_cmd", staticmethod(fake_run_cmd))

    project_dir = tmp_path / "proj"
    project_dir.mkdir()
    snakefile = str(tmp_path / "wes_single.smk")

    settings = {
        "project_dir": str(project_dir),
        "threads": 6,
        "run_id": "RIDPART",
        "debug": False,
        "genome": "b37",
        "inputs": {"input_dir": None, "sample_map": None},
        "snakemake_parameters": {"target": "call_variants"},
        "nextflow_parameters": {},
        "run_mode": "partial",
        "workflow": {
            "backend": "snakemake",
            "pipeline": "wes",
            "mode": "single",
            "software_stack": "gatk-4.6",
            "registry_version": "v1",
            "entrypoint": snakefile,
            "config_file": str(tmp_path / "config.yaml"),
            "helpers": {},
        },
    }

    assert execution.WorkflowExecutor(settings).run() is True
    assert recorded["cmd"][:3] == ["snakemake", "--forceall", "call_variants"]


def test_execution_builds_nextflow_command_and_helpers(tmp_path, monkeypatch):
    recorded = {}

    def fake_run_cmd(cmd, cwd, log_path, env=None, backend=None):
        recorded.update({"cmd": cmd, "backend": backend})

    monkeypatch.setattr(execution.WorkflowExecutor, "_run_cmd", staticmethod(fake_run_cmd))

    project_dir = tmp_path / "proj"
    project_dir.mkdir()
    workflow = str(tmp_path / "wgs_single.nf")
    config = str(tmp_path / "config.yaml")

    settings = {
        "project_dir": str(project_dir),
        "threads": 8,
        "run_id": "RIDNF",
        "debug": False,
        "genome": "hg38",
        "qc_coverage_region": "chr22",
        "inputs": {"input_dir": None, "sample_map": None},
        "snakemake_parameters": {},
        "nextflow_parameters": {"emit_report": True, "scatter_count": 2},
        "run_mode": "full",
        "cleanup_bam": True,
        "workflow": {
            "backend": "nextflow",
            "pipeline": "wgs",
            "mode": "single",
            "software_stack": "gatk-4.6",
            "registry_version": "v1",
            "entrypoint": workflow,
            "config_file": config,
            "helpers": {
                "coverage": str(tmp_path / "coverage.sh"),
                "vcf2sex": str(tmp_path / "vcf2sex.sh"),
                "vcf2hash": str(tmp_path / "vcf2hash.sh"),
            },
        },
    }

    assert execution.WorkflowExecutor(settings).run() is True
    cmd = recorded["cmd"]
    assert cmd[:3] == ["nextflow", "run", workflow]
    assert "-params-file" in cmd
    assert config in cmd
    assert "--pipeline" in cmd and "wgs" in cmd
    assert "--genome" in cmd and "hg38" in cmd
    assert "--cleanup_bam" in cmd and "true" in cmd
    assert "--qc_coverage_region" in cmd and "chr22" in cmd
    assert "--vcf2hash_script" in cmd
    assert "--emit_report" in cmd
    assert "--scatter_count" in cmd and "2" in cmd
    assert recorded["backend"] == "nextflow"
    contract = json.loads((project_dir / "cbicall-execution-contract.json").read_text(encoding="utf-8"))
    assert contract["workflow"]["backend"] == "nextflow"
    assert contract["workflow"]["provider"] == "cbicall"
    assert contract["command"]["normalized_sha256"]


def test_execution_builds_nextflow_cohort_command_with_sample_map(tmp_path, monkeypatch):
    recorded = {}

    def fake_run_cmd(cmd, cwd, log_path, env=None, backend=None):
        recorded.update({"cmd": cmd, "backend": backend})

    monkeypatch.setattr(execution.WorkflowExecutor, "_run_cmd", staticmethod(fake_run_cmd))

    project_dir = tmp_path / "proj"
    project_dir.mkdir()
    sample_map = tmp_path / "sample_map.tsv"
    sample_map.write_text("S1\t/s1.g.vcf.gz\n", encoding="utf-8")

    settings = {
        "project_dir": str(project_dir),
        "threads": 4,
        "run_id": "RIDNFCOHORT",
        "debug": False,
        "genome": "b37",
        "inputs": {"input_dir": None, "sample_map": str(sample_map)},
        "snakemake_parameters": {},
        "nextflow_parameters": {},
        "run_mode": "full",
        "cleanup_bam": False,
        "workflow": {
            "backend": "nextflow",
            "pipeline": "wes",
            "mode": "cohort",
            "software_stack": "gatk-4.6",
            "registry_version": "v1",
            "entrypoint": str(tmp_path / "wes_cohort.nf"),
            "config_file": str(tmp_path / "config.yaml"),
            "helpers": {"vcf2hash": str(tmp_path / "vcf2hash.sh")},
        },
    }

    assert execution.WorkflowExecutor(settings).run() is True
    cmd = recorded["cmd"]
    assert "--sample_map" in cmd and str(sample_map) in cmd
    assert "--workspace" in cmd and "cohort.genomicsdb.RIDNFCOHORT" in cmd
    assert "--vcf2hash_script" in cmd


def test_execution_nextflow_cohort_requires_sample_map(tmp_path):
    project_dir = tmp_path / "proj"
    project_dir.mkdir()

    settings = {
        "project_dir": str(project_dir),
        "threads": 4,
        "run_id": "RIDNFCOHORT",
        "debug": False,
        "genome": "b37",
        "inputs": {"input_dir": None, "sample_map": None},
        "snakemake_parameters": {},
        "nextflow_parameters": {},
        "run_mode": "full",
        "workflow": {
            "backend": "nextflow",
            "pipeline": "wes",
            "mode": "cohort",
            "software_stack": "gatk-4.6",
            "registry_version": "v1",
            "entrypoint": str(tmp_path / "wes_cohort.nf"),
            "config_file": str(tmp_path / "config.yaml"),
            "helpers": {},
        },
    }

    with pytest.raises(WorkflowResolutionError, match="sample_map is required"):
        execution.WorkflowExecutor(settings).run()


def test_execution_builds_nfcore_sarek_command_and_params_file(tmp_path, monkeypatch):
    recorded = {}
    monkeypatch.setattr(execution.platform, "machine", lambda: "x86_64")

    def fake_run_cmd(cmd, cwd, log_path, env=None, backend=None):
        recorded.update({"cmd": cmd, "cwd": cwd, "log_path": log_path, "backend": backend, "env": env})

    monkeypatch.setattr(execution.WorkflowExecutor, "_run_cmd", staticmethod(fake_run_cmd))

    project_dir = tmp_path / "proj"
    project_dir.mkdir()
    sample_map = tmp_path / "samplesheet.csv"
    sample_map.write_text("patient,sample,lane,fastq_1,fastq_2\n", encoding="utf-8")

    settings = {
        "project_dir": str(project_dir),
        "threads": 2,
        "run_id": "RIDSAREK",
        "debug": False,
        "genome": "external",
        "nfcore_profile": "docker",
        "nfcore_singularity_cache_dir": str(tmp_path / "nxf-cache"),
        "nfcore_parameters": {
            "input": str(sample_map),
            "genome": "GATK.GRCh38",
            "tools": "haplotypecaller",
            "max_memory": "30.GB",
        },
        "inputs": {"input_dir": None, "sample_map": None},
        "snakemake_parameters": {},
        "nextflow_parameters": {},
        "run_mode": "full",
        "workflow": {
            "backend": "nextflow",
            "pipeline": "sarek",
            "mode": "cohort",
            "software_stack": "nf-core",
            "registry_version": "v1",
            "entrypoint": "nf-core/sarek",
            "config_file": None,
            "helpers": {},
            "metadata": {
                "provider": "nf-core",
                "source": "nf-core/sarek",
                "release": "3.8.1",
                "default_outdir": "sarek",
            },
        },
    }

    assert execution.WorkflowExecutor(settings).run() is True
    cmd = recorded["cmd"]
    assert cmd[:5] == ["nextflow", "run", "nf-core/sarek", "-r", "3.8.1"]
    assert "-profile" in cmd and "docker" in cmd
    assert "-c" in cmd and str(project_dir / "cbicall_external_nextflow.config") in cmd
    assert "-work-dir" in cmd and str(project_dir / "work") in cmd
    params_file = project_dir / "cbicall_external_nextflow.params.yaml"
    config_file = project_dir / "cbicall_external_nextflow.config"
    assert params_file.is_file()
    assert config_file.is_file()
    params_text = params_file.read_text(encoding="utf-8")
    assert f"input: {sample_map}" in params_text
    assert f"outdir: {project_dir / 'sarek'}" in params_text
    assert "genome: GATK.GRCh38" in params_text
    assert "max_cpus: 2" in params_text
    assert "max_memory: 30.GB" in params_text
    assert "tools: haplotypecaller" in params_text
    config_text = config_file.read_text(encoding="utf-8")
    assert "resourceLimits = [ cpus: 2, memory: 30.GB ]" in config_text
    assert "singularity {" in config_text
    assert "apptainer {" in config_text
    assert f"cacheDir = '{tmp_path / 'nxf-cache'}'" in config_text
    assert f"libraryDir = '{tmp_path / 'nxf-cache'}'" in config_text
    assert "NXF_SINGULARITY_CACHEDIR" not in recorded["env"]
    assert "NXF_SINGULARITY_LIBRARYDIR" not in recorded["env"]
    assert "NXF_APPTAINER_CACHEDIR" not in recorded["env"]
    assert "NXF_APPTAINER_LIBRARYDIR" not in recorded["env"]
    assert recorded["cwd"] == project_dir
    assert recorded["log_path"] == project_dir / "nf-core_sarek_cohort.log"
    assert recorded["backend"] == "nextflow"
    contract = json.loads((project_dir / "cbicall-execution-contract.json").read_text(encoding="utf-8"))
    generated = {item["role"]: item for item in contract["generated_files"]}
    assert set(generated) == {"nf-core:params", "nf-core:config"}
    assert generated["nf-core:params"]["status"] == "present"
    assert generated["nf-core:config"]["normalized_sha256"]
    assert contract["workflow"]["provider"] == "nf-core"


def test_execution_nfcore_sarek_pins_amd64_docker_platform_on_arm64(tmp_path, monkeypatch):
    recorded = {}
    monkeypatch.setattr(execution.platform, "machine", lambda: "aarch64")

    def fake_run_cmd(cmd, cwd, log_path, env=None, backend=None):
        recorded.update({"cmd": cmd, "cwd": cwd, "backend": backend})

    monkeypatch.setattr(execution.WorkflowExecutor, "_run_cmd", staticmethod(fake_run_cmd))

    project_dir = tmp_path / "proj"
    project_dir.mkdir()
    sample_map = tmp_path / "samplesheet.csv"
    sample_map.write_text("patient,sample,lane,fastq_1,fastq_2\n", encoding="utf-8")

    settings = {
        "project_dir": str(project_dir),
        "threads": 6,
        "run_id": "RIDSAREKARM",
        "debug": False,
        "genome": "external",
        "nfcore_profile": "docker",
        "nfcore_parameters": {
            "input": str(sample_map),
            "genome": "GATK.GRCh38",
            "tools": "haplotypecaller",
        },
        "inputs": {"input_dir": None, "sample_map": None},
        "snakemake_parameters": {},
        "nextflow_parameters": {},
        "run_mode": "full",
        "workflow": {
            "backend": "nextflow",
            "pipeline": "sarek",
            "mode": "cohort",
            "software_stack": "nf-core",
            "registry_version": "v1",
            "entrypoint": "nf-core/sarek",
            "config_file": None,
            "helpers": {},
            "metadata": {
                "provider": "nf-core",
                "source": "nf-core/sarek",
                "release": "3.8.1",
                "default_outdir": "sarek",
            },
        },
    }

    assert execution.WorkflowExecutor(settings).run() is True
    config_text = (project_dir / "cbicall_external_nextflow.config").read_text(encoding="utf-8")
    assert "resourceLimits = [ cpus: 6 ]" in config_text
    assert "runOptions = '--platform linux/amd64'" in config_text


def test_snakemake_missing_snakefile_raises(tmp_path):
    project_dir = tmp_path / "proj"
    project_dir.mkdir()

    settings = {
        "project_dir": str(project_dir),
        "threads": 2,
        "run_id": "RIDSMK",
        "debug": False,
        "genome": "b37",
        "inputs": {"input_dir": None, "sample_map": None},
        "snakemake_parameters": {},
        "nextflow_parameters": {},
        "run_mode": "full",
        "workflow": {
            "backend": "snakemake",
            "pipeline": "wes",
            "mode": "single",
            "software_stack": "gatk-4.6",
            "registry_version": "v1",
            "entrypoint": None,
            "config_file": None,
            "helpers": {},
        },
    }

    with pytest.raises(RuntimeError, match="Missing Snakefile"):
        execution.WorkflowExecutor(settings).run()


def test_execution_raises_if_projectdir_missing(tmp_path):
    settings = {
        "project_dir": str(tmp_path / "does_not_exist"),
        "threads": 2,
        "run_id": "RID0",
        "debug": False,
        "inputs": {"input_dir": None, "sample_map": None},
        "snakemake_parameters": {},
        "nextflow_parameters": {},
        "run_mode": "full",
        "workflow": {
            "backend": "bash",
            "pipeline": "wes",
            "mode": "single",
            "software_stack": "gatk-3.5",
            "registry_version": "v1",
            "entrypoint": str(tmp_path / "wes_single.sh"),
            "config_file": None,
            "helpers": {},
        },
    }
    with pytest.raises(RuntimeError, match="Project directory does not exist"):
        execution.WorkflowExecutor(settings).run()


def _write_cromwell_config(tmp_path):
    config = tmp_path / "config.yaml"
    config.write_text(
        "datadir: /data\n"
        "dbdir: '{datadir}/Databases'\n"
        "ngsutils: '{datadir}/NGSutils'\n"
        "tmpdir: '{datadir}/tmp'\n"
        "mem: 8G\n"
        "mem_genotype: 64G\n"
        "gatk4_cmd: '{ngsutils}/gatk/gatk --java-options -Xmx{mem}'\n"
        "tools:\n"
        "  amd64:\n"
        "    bwa: '{ngsutils}/bwa/bwa'\n"
        "    samtools: '{ngsutils}/samtools/samtools'\n"
        "  aarch64:\n"
        "    bwa: '{ngsutils}/bwa-arm/bwa'\n"
        "    samtools: '{ngsutils}/samtools-arm/samtools'\n"
        "resources:\n"
        "  b37:\n"
        "    bundle: '{dbdir}/GATK_bundle/b37'\n"
        "    ref: '{bundle}/ref.fasta'\n"
        "    refgz: '{bundle}/ref.fasta.gz'\n"
        "    dbsnp: '{dbdir}/dbsnp.vcf.gz'\n"
        "    mills_indels: '{bundle}/mills.vcf.gz'\n"
        "    kg_indels: '{bundle}/kg.vcf.gz'\n"
        "    hapmap: '{bundle}/hapmap.vcf.gz'\n"
        "    omni: '{bundle}/omni.vcf.gz'\n"
        "    interval_list: '{bundle}/targets.interval_list'\n"
        "    snp_res: '-resource:dbsnp,known=true {dbsnp}'\n"
        "    indel_res: '-resource:mills,known=true {mills_indels}'\n",
        encoding="utf-8",
    )
    return config


def test_execution_builds_and_promotes_cromwell_wes_single(tmp_path, monkeypatch):
    recorded = {}
    project_dir = tmp_path / "CNAG99901P_ex" / "cbicall_cromwell_run"
    project_dir.mkdir(parents=True)
    input_dir = project_dir.parent
    (input_dir / "CNAG99901P_ex_S2_L001_R1_001.fastq.gz").write_text("r1\n", encoding="utf-8")
    (input_dir / "CNAG99901P_ex_S2_L001_R2_001.fastq.gz").write_text("r2\n", encoding="utf-8")
    jar = tmp_path / "cromwell.jar"
    jar.write_text("jar\n", encoding="utf-8")
    monkeypatch.setenv("CROMWELL_JAR", str(jar))
    monkeypatch.setenv("JAVA_CMD", "/usr/bin/java")
    monkeypatch.setattr(execution.platform, "machine", lambda: "x86_64")

    def fake_run_cmd(cmd, cwd, log_path, env=None, backend=None):
        recorded.update({"cmd": cmd, "cwd": cwd, "backend": backend})
        log_path.write_text("cromwell log\n", encoding="utf-8")
        task = tmp_path / "task-output"
        (task / "02_varcall").mkdir(parents=True)
        (task / "03_stats").mkdir()
        (task / "logs").mkdir()
        outputs = {}
        for name in ("hc.g.vcf.gz", "hc.raw.vcf.gz", "hc.QC.vcf.gz"):
            path = task / "02_varcall" / f"CNAG99901P.{name}"
            path.write_text(name + "\n", encoding="utf-8")
            key = {"hc.g.vcf.gz": "gvcf", "hc.raw.vcf.gz": "raw_vcf", "hc.QC.vcf.gz": "qc_vcf"}[name]
            outputs[f"CBIcallWesSingle.{key}"] = str(path)
        for key, filename in {
            "coverage": "CNAG99901P.coverage.txt",
            "sex": "CNAG99901P.sex.txt",
            "vcf_hash": "CNAG99901P.vcf.sha256.txt",
        }.items():
            path = task / "03_stats" / filename
            path.write_text(key + "\n", encoding="utf-8")
            outputs[f"CBIcallWesSingle.{key}"] = str(path)
        log = task / "logs" / "CNAG99901P.log"
        log.write_text("task log\n", encoding="utf-8")
        outputs["CBIcallWesSingle.logs"] = [str(log)]
        (cwd / "cbicall_cromwell.metadata.json").write_text(
            __import__("json").dumps({"outputs": outputs}),
            encoding="utf-8",
        )

    monkeypatch.setattr(execution.WorkflowExecutor, "_run_cmd", staticmethod(fake_run_cmd))
    wdl = tmp_path / "wes_single.wdl"
    wdl.write_text("version 1.0\n", encoding="utf-8")
    config = _write_cromwell_config(tmp_path)

    settings = {
        "project_dir": str(project_dir),
        "threads": 2,
        "run_id": "RIDCROM",
        "debug": False,
        "profile": "local",
        "genome": "b37",
        "qc_coverage_region": "chr22",
        "inputs": {"input_dir": str(input_dir), "sample_map": None},
        "snakemake_parameters": {},
        "nextflow_parameters": {},
        "cromwell_parameters": {"extra_label": "audit"},
        "run_mode": "full",
        "cleanup_bam": False,
        "workflow": {
            "backend": "cromwell",
            "pipeline": "wes",
            "mode": "single",
            "software_stack": "gatk-4.6",
            "registry_version": "v1",
            "entrypoint": str(wdl),
            "config_file": str(config),
            "helpers": {
                "coverage": str(tmp_path / "coverage.sh"),
                "vcf2sex": str(tmp_path / "vcf2sex.sh"),
                "vcf2hash": str(tmp_path / "vcf2hash.sh"),
            },
        },
    }

    assert execution.WorkflowExecutor(settings).run() is True
    cmd = recorded["cmd"]
    assert cmd[:4] == ["/usr/bin/java", "-jar", str(jar), "run"]
    assert "--inputs" in cmd and str(project_dir / "cbicall_cromwell.inputs.json") in cmd
    assert "--options" in cmd and str(project_dir / "cbicall_cromwell.options.json") in cmd
    assert "--metadata-output" in cmd
    assert recorded["backend"] == "cromwell"
    inputs = __import__("json").loads((project_dir / "cbicall_cromwell.inputs.json").read_text(encoding="utf-8"))
    assert inputs["CBIcallWesSingle.id"] == "CNAG99901P"
    assert inputs["CBIcallWesSingle.bwa"] == "/data/NGSutils/bwa/bwa"
    assert inputs["CBIcallWesSingle.ref"] == "/data/Databases/GATK_bundle/b37/ref.fasta"
    assert inputs["CBIcallWesSingle.qc_coverage_region"] == "chr22"
    assert inputs["CBIcallWesSingle.extra_label"] == "audit"
    assert (project_dir / "cbicall_cromwell.fastq_pairs.tsv").read_text(encoding="utf-8").count("\n") == 1
    assert (project_dir / "02_varcall" / "CNAG99901P.hc.QC.vcf.gz").is_file()
    assert (project_dir / "03_stats" / "CNAG99901P.vcf.sha256.txt").is_file()
    assert (project_dir / "logs" / "CNAG99901P.log").is_file()
    contract = json.loads((project_dir / "cbicall-execution-contract.json").read_text(encoding="utf-8"))
    generated = {item["role"]: item for item in contract["generated_files"]}
    assert set(generated) == {"cromwell:inputs", "cromwell:options", "cromwell:fastq_pairs"}
    assert generated["cromwell:inputs"]["status"] == "present"
    assert generated["cromwell:options"]["sha256"]
    assert contract["backend_parameters"]["cromwell_parameters"] == {"extra_label": "audit"}


def test_execution_builds_and_promotes_cromwell_wgs_cohort(tmp_path, monkeypatch):
    recorded = {}
    project_dir = tmp_path / "cohort" / "cbicall_cromwell_run"
    project_dir.mkdir(parents=True)
    sample_map = tmp_path / "sample_map.tsv"
    sample_map.write_text("S1\t/path/S1.g.vcf.gz\nS2\t/path/S2.g.vcf.gz\n", encoding="utf-8")
    jar = tmp_path / "cromwell.jar"
    jar.write_text("jar\n", encoding="utf-8")
    monkeypatch.setenv("CROMWELL_JAR", str(jar))
    monkeypatch.setattr(execution.platform, "machine", lambda: "x86_64")

    def fake_run_cmd(cmd, cwd, log_path, env=None, backend=None):
        recorded.update({"cmd": cmd, "cwd": cwd, "backend": backend})
        log_path.write_text("cromwell log\n", encoding="utf-8")
        task = tmp_path / "cohort-task-output"
        (task / "02_varcall").mkdir(parents=True)
        (task / "03_stats").mkdir()
        (task / "logs").mkdir()
        outputs = {}
        for key, filename in {
            "raw_vcf": "cohort.gv.raw.vcf.gz",
            "qc_vcf": "cohort.gv.QC.vcf.gz",
        }.items():
            path = task / "02_varcall" / filename
            path.write_text(key + "\n", encoding="utf-8")
            outputs[f"CBIcallCohort.{key}"] = str(path)
        hash_path = task / "03_stats" / "cohort.gv.QC.vcf.sha256.txt"
        hash_path.write_text("hash\n", encoding="utf-8")
        outputs["CBIcallCohort.vcf_hash"] = str(hash_path)
        log = task / "logs" / "cohort_joint_genotyping.log"
        log.write_text("task log\n", encoding="utf-8")
        outputs["CBIcallCohort.logs"] = [str(log)]
        (cwd / "cbicall_cromwell.metadata.json").write_text(
            __import__("json").dumps({"outputs": outputs}),
            encoding="utf-8",
        )

    monkeypatch.setattr(execution.WorkflowExecutor, "_run_cmd", staticmethod(fake_run_cmd))
    wdl = tmp_path / "wgs_cohort.wdl"
    wdl.write_text("version 1.0\n", encoding="utf-8")
    config = _write_cromwell_config(tmp_path)

    settings = {
        "project_dir": str(project_dir),
        "threads": 3,
        "run_id": "RIDCOHORT",
        "debug": False,
        "profile": "local",
        "genome": "b37",
        "inputs": {"input_dir": None, "sample_map": str(sample_map)},
        "snakemake_parameters": {},
        "nextflow_parameters": {},
        "cromwell_parameters": {"min_snp_for_vqsr": 10},
        "run_mode": "full",
        "cleanup_bam": False,
        "workflow": {
            "backend": "cromwell",
            "pipeline": "wgs",
            "mode": "cohort",
            "software_stack": "gatk-4.6",
            "registry_version": "v1",
            "entrypoint": str(wdl),
            "config_file": str(config),
            "helpers": {
                "coverage": str(tmp_path / "coverage.sh"),
                "vcf2sex": str(tmp_path / "vcf2sex.sh"),
                "vcf2hash": str(tmp_path / "vcf2hash.sh"),
            },
        },
    }

    assert execution.WorkflowExecutor(settings).run() is True
    inputs = __import__("json").loads((project_dir / "cbicall_cromwell.inputs.json").read_text(encoding="utf-8"))
    assert inputs["CBIcallCohort.pipeline"] == "wgs"
    assert inputs["CBIcallCohort.sample_map"] == str(sample_map.resolve())
    assert inputs["CBIcallCohort.workspace"] == "cohort.genomicsdb.RIDCOHORT"
    assert inputs["CBIcallCohort.gatk4_cmd"].endswith("-Xmx64G")
    assert inputs["CBIcallCohort.min_snp_for_vqsr"] == 10
    assert not (project_dir / "cbicall_cromwell.fastq_pairs.tsv").exists()
    assert (project_dir / "02_varcall" / "cohort.gv.QC.vcf.gz").is_file()
    assert (project_dir / "03_stats" / "cohort.gv.QC.vcf.sha256.txt").is_file()
    assert (project_dir / "logs" / "cohort_joint_genotyping.log").is_file()


def test_execution_cromwell_requires_runtime(tmp_path, monkeypatch):
    monkeypatch.delenv("CROMWELL_JAR", raising=False)
    monkeypatch.setenv("PATH", str(tmp_path / "empty-path"))
    monkeypatch.setattr(execution.shutil, "which", lambda name: None)
    project_dir = tmp_path / "CNAG99901P_ex" / "run"
    project_dir.mkdir(parents=True)
    (project_dir.parent / "S_R1_001.fastq.gz").write_text("r1\n", encoding="utf-8")
    (project_dir.parent / "S_R2_001.fastq.gz").write_text("r2\n", encoding="utf-8")
    wdl = tmp_path / "wes_single.wdl"
    wdl.write_text("version 1.0\n", encoding="utf-8")
    config = _write_cromwell_config(tmp_path)
    settings = {
        "project_dir": str(project_dir),
        "threads": 1,
        "run_id": "RID",
        "debug": False,
        "profile": "local",
        "genome": "b37",
        "inputs": {"input_dir": str(project_dir.parent), "sample_map": None},
        "snakemake_parameters": {},
        "nextflow_parameters": {},
        "cromwell_parameters": {},
        "cleanup_bam": False,
        "workflow": {
            "backend": "cromwell",
            "pipeline": "wes",
            "mode": "single",
            "software_stack": "gatk-4.6",
            "registry_version": "v1",
            "entrypoint": str(wdl),
            "config_file": str(config),
            "helpers": {"coverage": "coverage.sh", "vcf2sex": "vcf2sex.sh", "vcf2hash": "vcf2hash.sh"},
        },
    }
    with pytest.raises(WorkflowResolutionError, match="Cromwell is not available"):
        execution.WorkflowExecutor(settings).run()

def test_run_raises_on_invalid_backend(tmp_path):
    project_dir = tmp_path / "proj"
    project_dir.mkdir()
    settings = {
        "project_dir": str(project_dir),
        "threads": 2,
        "run_id": "RIDBAD",
        "debug": False,
        "inputs": {"input_dir": None, "sample_map": None},
        "snakemake_parameters": {},
        "nextflow_parameters": {},
        "run_mode": "full",
        "workflow": {
            "backend": "nope",
            "pipeline": "wes",
            "mode": "single",
            "software_stack": "gatk-3.5",
            "registry_version": "v1",
            "entrypoint": None,
            "config_file": None,
            "helpers": {},
        },
    }
    with pytest.raises(WorkflowResolutionError, match="Invalid workflow_backend"):
        execution.WorkflowExecutor(settings).run()


def test_run_cmd_nonzero_return_raises(monkeypatch, tmp_path):
    class P:
        returncode = 1

    def fake_run(cmd, cwd, env, stdout, stderr, check):
        return P()

    monkeypatch.setattr(execution.subprocess, "run", fake_run)
    with pytest.raises(WorkflowExecutionError, match="returncode=1") as excinfo:
        execution.WorkflowExecutor._run_cmd(
            cmd=["false"],
            cwd=tmp_path,
            log_path=tmp_path / "log.txt",
            env=None,
            backend="bash",
        )
    msg = str(excinfo.value)
    assert "backend=bash" in msg
    assert "Command: false" in msg
    assert f"Working directory: {tmp_path}" in msg


def test_run_cmd_exception_raises(monkeypatch, tmp_path):
    def fake_run(cmd, cwd, env, stdout, stderr, check):
        raise OSError("boom")

    monkeypatch.setattr(execution.subprocess, "run", fake_run)
    with pytest.raises(WorkflowExecutionError, match="could not start command") as excinfo:
        execution.WorkflowExecutor._run_cmd(
            cmd=["echo", "ok"],
            cwd=tmp_path,
            log_path=tmp_path / "log.txt",
            env=None,
            backend="snakemake",
        )
    msg = str(excinfo.value)
    assert "backend=snakemake" in msg
    assert "Command: echo ok" in msg
