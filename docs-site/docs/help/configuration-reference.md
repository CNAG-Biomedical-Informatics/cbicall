# Configuration Reference

CBIcall runs from a YAML parameters file plus the CLI thread setting.

```bash
bin/cbicall -p parameters.yaml -t 4
```

Unknown YAML keys are rejected, so misspellings fail early instead of being ignored.
Run configuration is defined in YAML. The CLI supplies runtime controls such as the parameter file, thread count, color output, and validation/test commands; it does not override YAML analysis keys.

:::tip[Minimal WES single-sample run]
```yaml
mode:            single
pipeline:        wes
workflow_engine: bash
gatk_version:    gatk-4.6
input_dir:       CNAG999_exome/CNAG99901P_ex
genome:          b37
```
:::

## Core Keys

| Key | Default | Values | Use |
| --- | --- | --- | --- |
| `mode` | `single` | `single`, `cohort` | Selects one-sample processing or cohort-level processing. |
| `pipeline` | `wes` | `wes`, `wgs`, `mit` | Selects the analysis type. |
| `workflow_engine` | `bash` | `bash`, `snakemake` | Selects the execution backend supported by the current workflows. |
| `profile` | `local` | `local`, `cnag-hpc` | Selects the runtime environment file. `cnag-hpc` uses `cnag-hpc-env.sh` instead of the default `env.sh` for Bash workflows. |
| `gatk_version` | `gatk-3.5` | `gatk-3.5`, `gatk-4.6` | Selects the workflow version. Use `gatk-4.6` for current WES/WGS workflows. |
| `resource_bundle` | `cbicall-germline-resources-v1` | local catalog key | Selects the external dependency bundle declared in `resources/cbicall-resource-catalog.json`. |
| `genome` | inferred | `b37`, `hg38`, `rsrs` | Reference genome. If omitted, CBIcall uses `b37` for WES/WGS and `rsrs` for mtDNA. |
| `input_dir` | `null` | path | Input sample or project directory. Relative paths are resolved from the YAML file location. |
| `sample_map` | `null` | path | Cohort-mode TSV containing sample IDs and gVCF paths. Relative paths are resolved from the YAML file location. |
| `project_dir` | `cbicall` | path or prefix | Prefix for the generated run directory. |
| `cleanup_bam` | `false` | `true`, `false` | Deletes intermediate BAM and BAI files after successful WES/WGS single-sample runs. |

## Compatibility Matrix

| Workflow | Supported |
| --- | --- |
| `gatk-4.6` + `bash` + `wes single/cohort` | Yes |
| `gatk-4.6` + `bash` + `wgs single/cohort` | Yes |
| `gatk-4.6` + `snakemake` + `wes single/cohort` | Yes |
| `gatk-4.6` + `snakemake` + `wgs single/cohort` | Yes |
| `gatk-3.5` + `bash` + `wes single/cohort` | Legacy |
| `gatk-3.5` + `bash` + `mit single/cohort` | Yes, x86_64 only |
| `mit` + `snakemake` | No |
| `gatk-3.5` + `snakemake` | No |

:::info[Genome rules]
- `pipeline: mit` always uses `genome: rsrs`.
- `genome: hg38` is supported only with `pipeline: wgs`.
- `pipeline: wes` currently uses `b37`.
:::

## Input Rules

### Single-Sample WES/WGS

Use `input_dir` pointing to the sample directory containing paired FASTQ files.

```yaml
mode:            single
pipeline:        wes
workflow_engine: bash
gatk_version:    gatk-4.6
input_dir:       CNAG999_exome/CNAG99901P_ex
genome:          b37
```

### Cohort WES/WGS

Use `sample_map` pointing to a TSV with sample identifiers and gVCF paths.

```yaml
mode:            cohort
pipeline:        wes
workflow_engine: bash
gatk_version:    gatk-4.6
genome:          b37
sample_map:      ./sample_map.tsv
```

### mtDNA

mtDNA workflows consume BAMs from previous WES/WGS runs. They do not start from FASTQ files.

```yaml
mode:            single
pipeline:        mit
workflow_engine: bash
gatk_version:    gatk-3.5
input_dir:       CNAG999_exome/CNAG99901P_ex
```

## Resource Bundle Provenance

`resource_bundle` pins the external tools and reference data expected for the run. CBIcall resolves this readable bundle key against `resources/cbicall-resource-catalog.json`, checks compatibility with the resolved workflow and pipeline implementation version, and writes a compact `resource_bundle` provenance block to `log.json`.

Two runs used the same declared external dependency set when their `resource_bundle.fingerprint` values match. The fingerprint is computed from the local catalog entry, so changes to tool versions, reference paths, archive naming, or compatibility metadata create a different value.

At preflight/runtime, CBIcall also resolves `DATADIR` from the selected `env.sh` or Snakemake config. If `cbicall-resource-installation.json` or `cbicall-bundle-id.json` exists in that directory, CBIcall checks that it matches the selected `resource_bundle`. The identifier file is verified by SHA-256 when the catalog pins one, and the installation manifest is checked against the local catalog fingerprint. Older manual installs without these metadata files are reported as `metadata_not_found`.

## Pipeline Implementation Version

Each registry entry also has a CBIcall pipeline implementation version, currently `v1` for the bundled workflows. This is different from `gatk_version`: `gatk_version` selects the tool family and workflow directory, while the pipeline implementation version identifies the CBIcall script/Snakefile revision selected inside that registry entry.

Normal YAML files do not need to set this. The workflow registry provides the default, currently `v1`.

Set `pipeline_version` only when a registry entry exposes more than one implementation and a run must pin a non-default one. The resolved value is written to `log.json` and `run-report.json`.

## Runtime Profiles

Profiles select the environment mapping used by a workflow. The default profile is `local`; additional profiles can be declared in the workflow registry when the same workflow needs more than one `env.sh` layout, for example on a shared HPC system.

Select a non-default profile in YAML:

```yaml
profile: cnag-hpc
```

Check the resolved setup without starting the workflow:

```bash
bin/cbicall doctor -p parameters.yaml
```

During a real run, the resolved `profile` and selected environment file are written to `log.json`. `doctor` prints the same resolved values without creating a run directory or log file.

## Command Utilities

| Command | Use |
| --- | --- |
| `bin/cbicall doctor -p parameters.yaml` | Dry-run preflight for one concrete run. It resolves the parameter YAML, workflow, profile env file, and resource bundle without launching the workflow. |
| `bin/cbicall validate-registry` | Developer check for `workflows/registry/workflows.yaml` against `workflows/schema/workflows.schema.json`. |
| `bin/cbicall validate-resources` | Developer check for `resources/cbicall-resource-catalog.json` and its workflow compatibility keys. |
| `bin/cbicall test --wes`, `--mit`, or `--all` | Runs the bundled integration examples without remembering the script path. |

```bash
bin/cbicall doctor -p parameters.yaml
bin/cbicall validate-registry
bin/cbicall validate-resources
bin/cbicall test --wes -t 1
bin/cbicall test --all -t 1
```

## Advanced Keys

| Key | Default | Use |
| --- | --- | --- |
| `pipeline_version` | Registry default, currently `v1` | Advanced pin for a specific CBIcall pipeline implementation. Leave unset for normal runs. |
| `workflow_rule` | `null` | Snakemake target for a partial run. Leave unset for normal full runs. |
| `allow_partial_run` | `false` | Must be `true` when `workflow_rule` is set. This prevents accidental partial starts. |
| `organism` | `Homo sapiens` | Metadata field. |
| `technology` | `Illumina HiSeq` | Metadata field. |

:::warning[Partial runs]
Partial runs are intended for targeted Snakemake execution and restarts. If `workflow_rule` is set without `allow_partial_run: true`, CBIcall refuses to start.
:::

## Output Directory Naming

Every run gets a generated directory:

```text
<project_dir>_<workflow_engine>_<pipeline>_<mode>_<genome>_<gatk_version>_<run-id>/
```

When `input_dir` is set, this directory is created inside `input_dir`.
See [Outputs](outputs) for the files produced by each workflow.
