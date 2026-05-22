# Native Backends

CBIcall separates **what analysis is run** from **how it is executed**.
[Included Pipelines](../pipelines/overview) describes the shipped WES, WGS, and
mtDNA analyses. This page describes the native execution backends maintained in
the CBIcall repository.

The workflow backend controls how the selected pipeline is launched.

| Native backend | Role | Typical use |
| --- | --- | --- |
| Bash | Runs CBIcall-maintained shell workflows directly | Production runs with the bundled WES, WGS, and mtDNA scripts |
| Snakemake | Runs CBIcall-maintained Snakemake workflows | Rule-based execution and partial workflow targets |
| Nextflow | Runs CBIcall-maintained Nextflow workflows | Alternative workflow-backend implementation for WES/WGS |

### Bash

The Bash backend runs CBIcall-maintained shell scripts from the `workflows/bash`
directory. It is the broadest backend in CBIcall and includes WES, WGS, and
mtDNA workflows.

### Snakemake

The Snakemake backend runs CBIcall-maintained Snakefiles from
`workflows/snakemake`. It supports WES/WGS workflows and is useful when rule
selection or partial execution is needed.

### Nextflow

The native Nextflow backend runs CBIcall-maintained Nextflow workflows from
`workflows/nextflow`. It supports WES/WGS workflows and follows CBIcall's run
directory and provenance model.

CBIcall can also launch registered external nf-core workflows through Nextflow.
Those workflows are documented separately because they keep their native output
layout and runtime assumptions. See [nf-core](nf-core).

## Backend Selection

For CBIcall-native workflows, select the backend with `workflow_backend`:

```yaml
pipeline:        wes
workflow_backend: bash
software_stack:    gatk-4.6
```

```yaml
pipeline:        wes
workflow_backend: snakemake
software_stack:    gatk-4.6
```

```yaml
pipeline:        wgs
workflow_backend: nextflow
software_stack:    gatk-4.6
```

## What to Read Next

| Goal | Page |
| --- | --- |
| Run WES/WGS single-sample workflows | [WES/WGS Single-Sample](../pipelines/wes-wgs-single) |
| Run WES/WGS cohort workflows | [WES/WGS Cohort](../pipelines/wes-wgs-cohort) |
| Run mtDNA workflows | [mtDNA](../pipelines/mtdna) |
| Run registered external nf-core workflows | [nf-core](nf-core) |
| Add a new workflow entry | [Adding a Pipeline](../technical-details/adding-a-pipeline) |
