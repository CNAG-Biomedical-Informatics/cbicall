# Quickstart

Use this page to check that CBIcall runs, then choose the workflow family you
want to exercise first.

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
bin/cbicall --help
bin/cbicall --version
```

You should see the command help and the installed CBIcall version.

## 2. Choose a First Test

CBIcall now has two useful first-run paths:

| Path | Needs the CBIcall resource bundle? | Use when |
| --- | --- | --- |
| nf-core provider | No | You want to test CBIcall orchestration immediately with Nextflow/nf-core. |
| Native CBIcall WES/mtDNA | Yes | You want to test the bundled Bash/Snakemake/Nextflow/Cromwell workflows against CBIcall's integration contracts. |

For nf-core, CBIcall validates the YAML and records provenance, while nf-core and
Nextflow manage the workflow's own test data, containers, and references.

## 3. Option A: Run nf-core Without the CBIcall Bundle

From `examples/input`, run the lightweight nf-core demo example:

```bash
cd examples/input
../../bin/cbicall validate-parameters -p nf-core-demo.yaml --no-color
../../bin/cbicall run -p nf-core-demo.yaml -t 4 --no-color
```

This does not require the CBIcall germline resource bundle or `DATADIR`. It does
require Nextflow and the container/runtime profile selected in the YAML, for
example `test,singularity` on HPC or `test,docker` on a Docker workstation.

## 4. Option B: Run the Native WES Example Test

From the repository root:

```bash
bin/cbicall test --wes-bash -t 1
```

This runs the bundled Bash WES workflow and validates the generated VCF against
the expected normalized hash declared by the integration contract. It requires
the CBIcall-provided resource bundle to be installed and configured.

:::tip[Need deeper checks?]
Use [Integration Tests](../validation/integration-tests) for bundled WES/mtDNA test details,
[Configuration Reference](../help/configuration-reference) for
`validate-parameters`, and [Resource Validation](resource-validation) for checking
the selected resource entry.
:::

## 5. Optional: Run the mtDNA Test

Run this after the WES test works:

```bash
bin/cbicall test --mit-bash -t 1
```

:::info[Architecture]
The mtDNA workflow uses MToolBox and is x86_64-only. If you are on ARM / aarch64, run WES/WGS workflows there but move mtDNA runs to an x86_64 host.
:::

## 6. Run With Your Own YAML

Once the integration tests work, use the normal invocation:

```bash
bin/cbicall run -p parameters.yaml -t 4
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
| Check reproducibility | [Run Comparison](../validation/run-comparison) |
| Understand generated files | [Outputs](../help/outputs) |
