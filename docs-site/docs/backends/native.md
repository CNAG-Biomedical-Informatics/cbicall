# Native Backends

CBIcall separates **what analysis is run** from **how it is executed**. The
analysis is selected with `pipeline` and `mode`; the execution technology is
selected with `workflow_backend`.

A pipeline implementation is **CBIcall-native** when it produces the CBIcall output contract:
the generated `cbicall_*` run directory, standard logs, audit reports, output
inventory, and final-output fingerprints when available. Native status is about
that output contract, not the workflow language.

## Supported Native Backends

| Backend | Workflow form | Required runner |
| --- | --- | --- |
| **Bash** | Shell scripts | Bash |
| **Snakemake** | Snakefile rules | Snakemake |
| **Nextflow** | Nextflow processes | Nextflow and Java |
| **Cromwell** | WDL workflows and tasks | Cromwell and Java |

:::tip[What `cbicall-core` means]
CBIcall is the framework. `cbicall-core` is the ready-to-run workflow collection
distributed and maintained with it.
:::

External nf-core provider entries are different: they are launched through the Nextflow
backend, but they keep their upstream nf-core output layout and runtime
assumptions. See [External nf-core](nf-core).

Use `cbicall doctor` to see which backend runners are available in the current
installation.

## Selecting a Backend

```yaml
pipeline: wes
mode: single
workflow_provider: cbicall-core
workflow_backend: bash
software_stack: gatk-4.6
```

For these bundled workflows, `workflow_provider: cbicall-core` is the default
and can be omitted. Change `workflow_backend` to `snakemake`, `nextflow`, or
`cromwell` when the selected pipeline, mode, and software stack support that
backend.

Use [Included Pipelines](../pipelines/overview) for the compatibility matrix and
pipeline-specific guides.
