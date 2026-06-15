# Native Backends

CBIcall separates **what analysis is run** from **how it is executed**. The
analysis is selected with `pipeline` and `mode`; the execution technology is
selected with `workflow_backend`.

A pipeline implementation is **CBIcall-native** when it produces the CBIcall output contract:
the generated `cbicall_*` run directory, standard logs, audit reports, output
inventory, and final-output fingerprints when available. Native status is about
that output contract, not the workflow language.

## Supported Native Backends

| Native backend | Current native role |
| --- | --- |
| Bash | Broadest native backend; includes WES, WGS, and mtDNA pipelines. |
| Snakemake | Native WES/WGS implementation for rule-based execution and partial targets. |
| Nextflow | Native WES/WGS implementation that follows the CBIcall run layout. |
| Cromwell | Native WDL implementation for GATK 4.6 WES/WGS pipelines. |

External nf-core provider entries are different: they are launched through the Nextflow
backend, but they keep their upstream nf-core output layout and runtime
assumptions. See [External nf-core](nf-core).

## Selecting a Backend

```yaml
pipeline: wes
mode: single
workflow_provider: cbicall
workflow_backend: bash
software_stack: gatk-4.6
```

For native pipeline implementations, `workflow_provider: cbicall` is the default. Change
`workflow_backend` to `snakemake`, `nextflow`, or `cromwell` when the selected
pipeline/mode/software stack supports that backend.

Use [Included Pipelines](../pipelines/overview) for the compatibility matrix and
pipeline-specific guides.
