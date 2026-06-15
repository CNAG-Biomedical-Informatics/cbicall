import WorkflowCompatibilityMatrix from '@site/src/components/WorkflowCompatibilityMatrix.mdx';

# Included Pipelines

This page summarizes the analysis pipelines CBIcall can launch and the backend combinations currently supported.

| Concept | YAML key | Meaning |
| --- | --- | --- |
| Pipeline | `pipeline` | The analysis family: `wes`, `wgs`, `mit`, or a registered external pipeline such as `sarek`. |
| Mode | `mode` | The run shape, usually `single` or `cohort`. |
| Backend | `workflow_backend` | The execution technology: Bash, Snakemake, Nextflow, or Cromwell. |
| Provider | `workflow_provider` | Who supplies the implementation: `cbicall` for native CBIcall pipelines, `nf-core` for registered external provider entries. |

A pipeline implementation is **native** when it produces the CBIcall output contract: standard
run directory, logs, reports, output inventory, and final-output fingerprints
when available. External provider entries can still be launched and audited by CBIcall,
but they keep their upstream output layout and runtime assumptions.

## Pipeline Guides

| Goal | Page |
| --- | --- |
| Process one WES/WGS sample from FASTQ | [WES/WGS Single-Sample](wes-wgs-single) |
| Joint-genotype a WES/WGS cohort from gVCFs | [WES/WGS Cohort](wes-wgs-cohort) |
| Run mitochondrial variant calling | [mtDNA](mtdna) |
| Run selected external nf-core provider entries | [External nf-core](../backends/nf-core) |

For a backend-centered view, see [Native Backends](../backends/native).

:::info[Backend-native validation]
CBIcall validates the parameters YAML against the workflow registry and resource
catalog, then writes a `cbicall-execution-contract.json` for the concrete command
it launches. Syntax or semantic validation inside each workflow language remains
backend-native: use `bash -n`, Snakemake lint/dry-run checks, Nextflow validation,
or `womtool validate` for WDL/Cromwell workflows.
:::

## Compatibility Matrix

<WorkflowCompatibilityMatrix />
