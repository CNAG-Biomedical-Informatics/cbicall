# Configuration Reference

CBIcall runs from a YAML parameters file plus CLI runtime settings.

```bash
bin/cbicall run -p parameters.yaml -t 4
```

Unknown YAML keys are rejected, so misspellings fail early instead of being ignored.
Analysis configuration is defined in YAML. Runtime controls such as thread count,
color output, validation commands, and the CBIcall native runtime profile are
selected on the CLI.

:::tip[Minimal WES single-sample run]
```yaml
mode:            single
pipeline:        wes
workflow_backend: bash
software_stack:    gatk-4.6
input_dir:       CNAG999_exome/CNAG99901P_ex
genome:          b37
```
:::

## Core Keys

| Key | Default | Values | Use |
| --- | --- | --- | --- |
| `mode` | `single` | `single`, `cohort` | Selects one-sample processing or cohort-level processing. |
| `pipeline` | `wes` | `wes`, `wgs`, `mit`; external names are registry-defined | Selects the analysis type. For `workflow_provider: nf-core`, the value is resolved through the workflow registry. |
| `workflow_backend` | `bash` | `bash`, `snakemake`, `nextflow` | Selects the execution backend supported by the current workflows. |
| `software_stack` | `gatk-3.5` | `gatk-3.5`, `gatk-4.6` | Selects the GATK release for CBIcall-native workflows. Use `gatk-4.6` for current bundled WES/WGS workflows. |
| `workflow_provider` | `cbicall` | `cbicall`, `nf-core` | Selects whether the workflow is a CBIcall-maintained implementation or an external nf-core workflow. Use `workflow_provider: nf-core` for external nf-core workflows. |
| `resource` | `cbicall-germline-resources-v1` | resource key | Selects one entry from `resources/cbicall-resource-catalog.json`. |
| `genome` | inferred | `b37`, `hg38`, `rsrs`, `external` | Reference genome. If omitted, CBIcall uses `b37` for WES/WGS, `rsrs` for mtDNA, and `external` for nf-core/Sarek. |
| `input_dir` | `null` | path | Input sample or project directory. Relative paths are resolved from the YAML file location. |
| `sample_map` | `null` | path | Cohort-mode TSV containing sample IDs and gVCF paths. Relative paths are resolved from the YAML file location. |
| `project_dir` | `cbicall` | path or prefix | Prefix for the generated run directory. |
| `cleanup_bam` | `false` | `true`, `false` | Deletes intermediate BAM and BAI files after successful WES/WGS single-sample runs. |

The resource catalog is the inventory of selectable resource entries and their
workflow compatibility metadata.

## Compatibility Matrix

| Workflow | Supported |
| --- | --- |
| `gatk-4.6` + `bash` + `wes single/cohort` | Yes |
| `gatk-4.6` + `bash` + `wgs single/cohort` | Yes |
| `gatk-4.6` + `snakemake` + `wes single/cohort` | Yes |
| `gatk-4.6` + `snakemake` + `wgs single/cohort` | Yes |
| `gatk-4.6` + `nextflow` + `wes single/cohort` | Yes |
| `gatk-4.6` + `nextflow` + `wgs single/cohort` | Yes |
| `nf-core` + `nextflow` + `demo single` | External nf-core smoke test |
| `nf-core` + `nextflow` + `sarek cohort` | External nf-core/Sarek workflow |
| `gatk-3.5` + `bash` + `wes single/cohort` | Legacy |
| `gatk-3.5` + `bash` + `mit single/cohort` | Yes, x86_64 only |
| `mit` + `snakemake` or `nextflow` | No |
| `gatk-3.5` + `snakemake` or `nextflow` | No |

:::info[Genome rules]
- `pipeline: mit` always uses `genome: rsrs`.
- External nf-core pipelines use `genome: external`; the reference is selected by the nf-core parameters in `nfcore_parameters`.
- `genome: hg38` is supported only with `pipeline: wgs`.
- `pipeline: wes` currently uses `b37`.
:::

## Input Rules

### Single-Sample WES/WGS

Use `input_dir` pointing to the sample directory containing paired FASTQ files.

```yaml
mode:            single
pipeline:        wes
workflow_backend: bash
software_stack:    gatk-4.6
input_dir:       CNAG999_exome/CNAG99901P_ex
genome:          b37
```

### Cohort WES/WGS

Use `sample_map` pointing to a TSV with sample identifiers and gVCF paths.

```yaml
mode:            cohort
pipeline:        wes
workflow_backend: bash
software_stack:    gatk-4.6
genome:          b37
sample_map:      ./sample_map.tsv
```

### mtDNA

mtDNA workflows consume BAMs from previous WES/WGS runs. They do not start from FASTQ files.

```yaml
mode:            single
pipeline:        mit
workflow_backend: bash
software_stack:    gatk-3.5
input_dir:       CNAG999_exome/CNAG99901P_ex
```

### nf-core/demo

`nf-core/demo` is useful for testing CBIcall's external nf-core support without
modeling a full biological workflow. It uses the nf-core `test` profile.

```yaml
mode:             single
pipeline:         demo
workflow_backend:  nextflow
workflow_provider: nf-core
resource:         nf-core-demo-managed-resources-v1
nfcore_profile: test,singularity
```

Use the checked-in `test,singularity` profile on HPC. Here, `test` supplies
nf-core's built-in demo inputs and smoke-test settings, while `singularity`
selects the Singularity/Apptainer runtime. On an x86_64 Docker workstation,
`test,docker` is also possible. If `nfcore_parameters.input` is set, it overrides
the input supplied by the `test` profile. For workstation and cluster runs, see
[nf-core External Workflows](../backends/nf-core).

### nf-core/Sarek

Sarek is launched as an external nf-core Nextflow workflow. CBIcall validates the
YAML, pins the registered nf-core release, writes a small params file in the run
directory, and leaves Sarek outputs in their native layout under `sarek/`.

```yaml
mode:             cohort
pipeline:         sarek
workflow_backend:  nextflow
workflow_provider: nf-core
resource:         nf-core-sarek-managed-resources-v1
nfcore_profile: singularity
# nfcore_singularity_cache_dir: nxf-singularity-cache
nfcore_parameters:
  input: sarek_samplesheet.csv
  genome: GATK.GRCh38
  tools: haplotypecaller
  skip_tools: haplotypecaller_filter
  wes: true
  intervals: ../../workflows/nf-core/sarek/grch38_chr22_test.bed
  max_memory: 30.GB
```

CBIcall does not interpret Sarek-specific parameters. Values under
`nfcore_parameters` are passed to the generated nf-core parameters file. Use the
samplesheet format and parameter names expected by the selected Sarek release.

For nf-core/Sarek, the CLI thread value is written to the generated params file
as `max_cpus`. For example, `bin/cbicall run -p nf-core-sarek.yaml -t 6` passes
`max_cpus: 6` to Sarek and writes a small Nextflow config with
`process.resourceLimits` so individual processes do not request more than six
CPUs. Memory caps stay in `nfcore_parameters`, for example `max_memory: 30.GB`;
CBIcall writes the same value to Nextflow `process.resourceLimits.memory`.
On HPC, set `nfcore_singularity_cache_dir` to a user- or project-owned
directory so the generated Nextflow config points away from unreadable
site-level container libraries. If the HPC module exports `NXF_*` variables,
keep those exports in the shell or SLURM bootstrap before invoking CBIcall.
For the tiny chr22 smoke test, `skip_tools: haplotypecaller_filter` avoids a
GATK `FilterVariantTranches` failure caused by too few overlapping resource
variants. Remove it for production Sarek runs if you want Sarek's default
HaplotypeCaller filtering.

## Bundle Provenance

`resource` selects the external tools and reference data expected for the run.
CBIcall checks that the selected resource is compatible with the resolved
workflow and records resource key, version, and fingerprint provenance in
`log.json` and `run-report.json`.

Use [Resource Validation](../usage/resource-validation) for resource checks and
[Run Comparison](../usage/run-comparison) to compare repeated runs.

## Pipeline Implementation Version

Each workflow registry entry has a CBIcall pipeline implementation version,
currently `v1` for the bundled workflows. Normal YAML files do not need to set
this; the registry provides the default.

Set `pipeline_version` only when a registry entry exposes more than one
implementation and a run must pin a non-default one.

## Runtime Profiles

CBIcall native profiles select the environment mapping used by a workflow. The
default profile is `local`; additional profiles can be declared in the workflow
registry when the same workflow needs more than one `env.sh` layout, for example
on a shared HPC system.

Select a non-default CBIcall profile on the CLI:

```bash
bin/cbicall run -p parameters.yaml -t 4 --runtime-profile cnag-hpc
```

Validate the parameters YAML and resolved setup without starting the workflow:

```bash
bin/cbicall validate-parameters -p parameters.yaml --runtime-profile cnag-hpc
```

The `profile` key is not accepted in the parameters YAML. During a real run, the
resolved profile and selected environment file are written to `log.json`.
`validate-parameters` prints the same resolved values without creating a run
directory or log file.

## Command Utilities

| Command | Use |
| --- | --- |
| `bin/cbicall run -p parameters.yaml -t 4 [--runtime-profile cnag-hpc]` | Execute a normal analysis run. |
| `bin/cbicall validate-parameters -p parameters.yaml [--runtime-profile cnag-hpc]` | Dry-run preflight for one concrete run. It validates the parameters YAML, workflow, runtime profile env file, and selected resource without launching the workflow. |
| `bin/cbicall validate-resources` | Check the resource catalog and, optionally, one resource key. |
| `bin/cbicall compare-runs RUN_A RUN_B [RUN_C ...]` | Compare two or more run directories or `run-report.json` files. |
| `bin/cbicall render-report RUN_DIR` | Regenerate `run-report.html` from an existing `run-report.json` without rerunning the workflow. |
| `bin/cbicall test --wes-bash [--runtime-profile cnag-hpc]`, `--wes-snakemake`, `--wes-nextflow`, `--mit-bash`, `--nf-core-demo`, `--nf-core-sarek`, or `--all` | Runs contract-based integration examples. `--runtime-profile` is forwarded to the internal `cbicall run` calls. `--wes-bash` is the required native baseline test; Snakemake, Nextflow, and nf-core tests require their backends on `PATH`. |

For a higher-level explanation of included pipelines versus execution backends,
see [Included Pipelines](../pipelines/overview) and
[Native Backends](../backends/native).

## Advanced Keys

| Key | Default | Use |
| --- | --- | --- |
| `pipeline_version` | Registry default, currently `v1` | Advanced pin for a specific CBIcall pipeline implementation. Leave unset for normal runs. |
| `snakemake_parameters` | `{}` | Snakemake-specific options. `target` selects a Snakemake target instead of the default `all`; other keys are passed through as extra `--config key=value` entries after CBIcall-managed config values. |
| `nextflow_parameters` | `{}` | Native CBIcall Nextflow parameters passed as `--key value`. CBIcall blocks keys it owns, such as `pipeline`, `genome`, `threads`, helper scripts, and cohort workspace settings. |
| `nfcore_profile` | `null` | nf-core profile passed to external nf-core workflows, for example `docker`, `singularity`, or `test,singularity`. |
| `nfcore_parameters` | `{}` | Pass-through nf-core parameters written to the generated params file. CBIcall controls `outdir` and `max_cpus`. |
| `nfcore_singularity_cache_dir` | `null` | Optional Singularity/Apptainer image cache directory for external nf-core workflows. CBIcall writes it to the generated Nextflow config as cache and library directories. |
| `organism` | `Homo sapiens` | Metadata field. |
| `technology` | `Illumina HiSeq` | Metadata field. |

:::note[Backend-specific parameters]
Use `snakemake_parameters`, `nextflow_parameters`, and `nfcore_parameters` only for parameters owned by that backend or external workflow. CBIcall still owns the compatibility contract and blocks overrides of core values it resolves itself.
:::

## Output Directory Naming

Every run gets a generated directory:

```text
<project_dir>_<workflow_backend>_<software_stack>_<pipeline>_<mode>_<genome>_<run-id>/
```

External nf-core workflows use `software_stack: nf-core`; the displayed genome
label is inferred from `nfcore_parameters.genome` when present:

```text
<project_dir>_nextflow_nf-core_<pipeline>_<mode>_<display-genome>_<run-id>/
```

When `input_dir` is set, this directory is created inside `input_dir`.
See [Outputs](outputs) for the files produced by each workflow.
