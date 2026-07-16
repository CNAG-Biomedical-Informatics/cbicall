# Quickstart

Use this page to inspect CBIcall's reports immediately, then choose the workflow
family you want to execute first.

If you still need to choose an installation method or workflow, start with the
[Overview](../overview).

:::info[YAML contract]
In CBIcall, the parameters YAML becomes a **YAML contract** after CBIcall has
validated it and resolved it against the workflow registry and resource catalog.
Both `validate-parameters` and `run` perform this validation; `validate-parameters`
stops before launching the workflow.
:::

## 1. Confirm the CLI

```bash
cbicall --help
cbicall --version
cbicall doctor
```

You should see the command help, the installed CBIcall version, and a concise
installation report. `doctor` checks the packaged contracts, `CBICALL_DATA`
bundle metadata, and available workflow backends. Missing optional backends are
reported as warnings.

## 2. Explore the Reports

Generate a WES audit report and an interactive mtDNA browser from packaged
example outputs:

```bash
cbicall demo
```

CBIcall writes both reports under `cbicall-demo/` and prints their paths. This
command needs no external resource bundle, workflow backend, Java installation,
or container runtime.

:::note[Precomputed demonstration]
The demo uses sanitized, precomputed outputs from the packaged CNAG99901P
integration fixture. It exercises CBIcall's reporting code but does not execute
BWA, GATK, or MToolBox, and it is not an analytical benchmark.
:::

Use a different empty destination when needed:

```bash
cbicall demo --output-dir my-cbicall-demo
```

## 3. Choose an Execution Test

For actual workflow execution, choose one of these paths:

| Path | CBIcall bundle | Other requirements | Use when |
| --- | --- | --- | --- |
| nf-core provider | No | Nextflow and the selected container runtime | You want to test external-provider orchestration. |
| Native CBIcall WES/mtDNA | Yes | Tools supplied by the bundle | You want to execute packaged native workflows against CBIcall integration contracts. |

For nf-core, CBIcall validates the YAML and records provenance, while nf-core and
Nextflow manage the workflow's own test data, containers, and references.

## 4. Option A: Run nf-core Without the CBIcall Bundle

From `examples/input`, run the lightweight nf-core demo example:

```bash
cd examples/input
cbicall validate-parameters -p nf-core-demo.yaml --no-color
cbicall run -p nf-core-demo.yaml -t 4 --no-color
```

This does not require the CBIcall germline resource bundle or `DATADIR`. It does
require Nextflow and the container/runtime profile selected in the YAML, for
example `test,singularity` on HPC or `test,docker` on a Docker workstation.

## 5. Option B: Run the Native WES Example Test

Point CBIcall to the installed external bundle, then run the test:

```bash
export CBICALL_DATA=/absolute/path/to/cbicall-data
cbicall test --wes-bash -t 1
```

This runs the bundled Bash WES workflow and validates the generated VCF against
the expected normalized hash declared by the integration contract. It requires
the CBIcall-provided resource bundle to be installed. `CBICALL_DATA` is applied
consistently to the Bash, Snakemake, Nextflow, and Cromwell native backends.

:::tip[Need deeper checks?]
Use [Integration Tests](../validation/integration-tests) for bundled WES/mtDNA test details,
[Configuration Reference](../help/configuration-reference) for
`validate-parameters`, and [Resource Validation](resource-validation) for checking
the selected resource entry.
:::

## 6. Optional: Run the mtDNA Test

Run the WES and mtDNA integration contracts together so that the WES test
produces the BAM consumed by MToolBox:

```bash
cbicall test --wes-bash --mit-bash -t 1
```

:::info[Architecture]
The mtDNA workflow uses MToolBox and is x86_64-only. If you are on ARM / aarch64, run WES/WGS workflows there but move mtDNA runs to an x86_64 host.
:::

## 7. Run With Your Own YAML

Once the integration tests work, use the normal invocation:

```bash
cbicall run -p parameters.yaml -t 4
```

| Option | Meaning |
| --- | --- |
| `-p` | YAML parameters file. |
| `-t` | Threads passed to the workflow. |

For most WES/WGS runs, start with 4 threads and adjust after checking
[Performance](../help/performance).

## Next Steps

| Goal | Page |
| --- | --- |
| Run nf-core workflows | [nf-core Provider](../backends/nf-core) |
| Run real WES/WGS data | [End-to-end Example: WES](end-to-end-example-wes) |
| Run mtDNA analysis | [End-to-end Example: mtDNA](end-to-end-example-mit) |
| Check reproducibility | [Run Comparison](../help/run-comparison) |
| Understand generated files | [Outputs](../help/outputs) |
