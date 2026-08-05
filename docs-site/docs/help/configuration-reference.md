# Configuration Reference

CBIcall runs from a YAML parameters file plus CLI runtime settings.

```bash
cbicall run -p parameters.yaml -t 4
```

Unknown YAML keys are rejected, so misspellings fail early instead of being ignored.
Analysis configuration is defined in YAML. Runtime controls such as thread count,
color output, validation commands, and the Bash runtime profile are
selected on the CLI.

A **YAML contract** is the parameters YAML after CBIcall has validated and
resolved it against the workflow registry and resource catalog. `run` and
`validate-parameters` use the same validation and resolution path;
`validate-parameters` stops before launching the workflow.

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
| `workflow_backend` | `bash` | `bash`, `snakemake`, `nextflow`, `cromwell` | Selects the execution backend supported by the current workflows. |
| `software_stack` | `gatk-3.5` | `gatk-3.5`, `gatk-4.6` | Selects the GATK release for bundled workflows. Use `gatk-4.6` for current WES/WGS workflows. |
| `workflow_provider` | `cbicall-core` | `cbicall-core`, `nf-core` | Selects the registered workflow collection. `cbicall-core` contains CBIcall's bundled native implementations; use `nf-core` for registered external nf-core workflows. |
| `resource` | `cbicall-germline-resources-v1` | resource key | Selects one entry from `resources/cbicall-resource-catalog.json`. |
| `genome` | inferred | `b37`, `hg38`, `rsrs`, `external` | Reference genome. If omitted, CBIcall uses `b37` for WES/WGS, `rsrs` for mtDNA, and `external` for nf-core/Sarek. |
| `input_dir` | `null` | path | Input sample or project directory. Relative paths are resolved from the YAML file location. |
| `sample_map` | `null` | path | Cohort-mode TSV containing sample IDs and gVCF paths. Relative paths are resolved from the YAML file location. |
| `project_dir` | `cbicall` | path or prefix | Prefix for the generated run directory. |
| `cleanup_bam` | `false` | `true`, `false` | Deletes intermediate BAM and BAI files after successful WES/WGS single-sample runs. |
| `qc_coverage_region` | `chr1` | contig name | Contig used only for the lightweight coverage summary. It does not change variant-calling intervals. |

The resource catalog is the inventory of selectable resource entries and their
workflow compatibility metadata.

## Supported Combinations

See the [compatibility matrix](../pipelines/overview#compatibility-matrix) for
the pipeline, mode, software stack, genome, and backend combinations currently
supported by CBIcall.

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

#### Staged Cohort Runs

Bundled `cbicall-core` GATK 4.6 cohort runs can be split into shard jobs and one finalize job.
This is useful when running chromosomes in parallel on a scheduler.

| Key | Default | Use |
| --- | --- | --- |
| `cohort_stage` | `all` | `all` runs the complete workflow, `shard` stops after one raw interval VCF, and `finalize` filters a gathered raw VCF. |
| `interval_shard` | `null` | Contig or interval-list label required by `cohort_stage: shard`. |
| `input_vcf` | `null` | Gathered raw VCF required by `cohort_stage: finalize`. |
| `output_basename` | `null` | Output stem, such as `cohort.chr1` for a shard or `cohort` for the final callset. |

Shard one contig:

```yaml
mode:            cohort
pipeline:        wgs
workflow_backend: bash
software_stack:    gatk-4.6
genome:          hg38
sample_map:      ./sample_map.tsv
cohort_stage:    shard
interval_shard:  chr1
output_basename: cohort.chr1
```

:::note[GenomicsDB workspace]
There is no user-facing `workspace` key in the parameters YAML. CBIcall controls
GenomicsDB workspace names and creates one workspace per run under
`01_genomicsdb/cohort.genomicsdb.<run-id>`. Use `output_basename` for shard-specific
VCF names.
:::

After all raw shard VCFs have been concatenated and indexed, run final filtering:

```yaml
mode:            cohort
pipeline:        wgs
workflow_backend: bash
software_stack:    gatk-4.6
genome:          hg38
cohort_stage:    finalize
input_vcf:       ./cohort.gathered.gv.raw.vcf.gz
output_basename: cohort
```

Staged cohort keys are currently supported only with bundled
`software_stack: gatk-4.6`, `mode: cohort`, and `workflow_backend` set to
`bash`, `snakemake`, `nextflow`, or `cromwell`. See the WES/WGS cohort page for
a GNU parallel chromosome-sharding example.

### mtDNA

mtDNA workflows consume BAMs from previous bundled Bash WES/WGS runs. They do
not start from FASTQ files.

```yaml
mode:            single
pipeline:        mit
workflow_backend: bash
software_stack:    gatk-3.5
input_dir:       CNAG999_exome/CNAG99901P_ex
```

### nf-core/demo

`nf-core/demo` tests the external-provider path without the CBIcall resource
bundle:

```yaml
mode:             single
pipeline:         demo
workflow_backend:  nextflow
workflow_provider: nf-core
resource:         nf-core-demo-managed-resources-v1
nfcore_profile: test,singularity
```

Use `test,singularity` on HPC or `test,docker` on an x86_64 Docker workstation.

### nf-core/Sarek

Sarek is launched through the registered external nf-core provider:

```yaml
mode:             cohort
pipeline:         sarek
workflow_backend:  nextflow
workflow_provider: nf-core
resource:         nf-core-sarek-managed-resources-v1
nfcore_profile: singularity
nfcore_parameters:
  input: sarek_samplesheet.csv
  genome: GATK.GRCh38
  tools: haplotypecaller
  wes: true
```

CBIcall passes `nfcore_parameters` to the selected nf-core release and controls
the output directory and CPU limit. Use the upstream pipeline's parameter names
and samplesheet format. See [External nf-core Workflows](../backends/nf-core)
for profiles, container caches, memory limits, and the packaged smoke tests.

## Bundle Provenance

`resource` selects the external tools and reference data expected for the run.
CBIcall checks that the selected resource is compatible with the resolved
workflow and records resource key, version, and fingerprint provenance in
`log.json` and `run-report.json`.

Use [Resource Validation](../usage/resource-validation) for resource checks and
[Run Comparison](run-comparison) to compare repeated runs.

## Runtime Profiles

`--runtime-profile` is a CLI-only setting for bundled Bash environment mappings.
The default is `local`; a site can register another mapping such as `cnag-hpc`.

Select a non-default CBIcall profile on the CLI:

```bash
cbicall run -p parameters.yaml -t 4 --runtime-profile cnag-hpc
```

Validate the parameters YAML and resolved setup without starting the workflow:

```bash
cbicall validate-parameters -p parameters.yaml --runtime-profile cnag-hpc
```

The `profile` key is not accepted in the parameters YAML. During a real run, the
resolved profile and selected environment file are written to `log.json`.
`validate-parameters` prints the same resolved values without creating a run
directory or log file.

<details>
<summary>How backend profiles differ</summary>

For Bash, CBIcall resolves the selected environment file and passes it through
`CBICALL_ENV_FILE`. Snakemake, Nextflow, Cromwell, and nf-core use their own
configuration or profile mechanisms and do not use this Bash switch.

</details>

## Related commands

This page documents YAML configuration. See [General Usage](../usage) for
CLI syntax, [Integration Tests](../validation/integration-tests) for shipped
workflow checks, [Resource Validation](../usage/resource-validation) for the
resource catalog, and [Run Comparison](run-comparison) for audit reports.

## Advanced Keys

| Key | Default | Use |
| --- | --- | --- |
| `registry_version` | Registry default, currently `v1` | Advanced pin for a specific CBIcall registry version. Leave unset for normal runs. |
| `snakemake_parameters` | `{}` | Snakemake-specific options. `target` selects a Snakemake target instead of the default `all`; other keys are passed through as extra `--config key=value` entries after CBIcall-managed config values. |
| `nextflow_parameters` | `{}` | Extra parameters for bundled Nextflow workflows. CBIcall blocks keys it owns, such as `pipeline`, `genome`, threads, helper scripts, and cohort workspace settings. |
| `cromwell_parameters` | `{}` | Extra WDL inputs for bundled Cromwell workflows. CBIcall blocks overrides of tool paths, references, sample identity, pipeline, coverage region, and thread count. |
| `nfcore_profile` | `null` | nf-core profile passed to external nf-core workflows, for example `docker`, `singularity`, or `test,singularity`. |
| `nfcore_parameters` | `{}` | Pass-through nf-core parameters written to the generated params file. CBIcall controls `outdir` and `max_cpus`. |
| `nfcore_singularity_cache_dir` | `null` | Optional Singularity/Apptainer image cache directory for external nf-core workflows. CBIcall writes it to the generated Nextflow config as cache and library directories. |
| `organism` | `Homo sapiens` | Metadata field. |
| `technology` | `Illumina HiSeq` | Metadata field. |

:::note[Backend-specific parameters]
Use `snakemake_parameters`, `nextflow_parameters`, `cromwell_parameters`, and `nfcore_parameters` only for parameters owned by that backend or external workflow. CBIcall still owns the compatibility contract and blocks overrides of core values it resolves itself.
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
