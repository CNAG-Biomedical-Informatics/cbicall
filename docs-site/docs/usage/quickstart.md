# Quickstart

Use this page to check that CBIcall runs on the shipped example data, then run
your own YAML file.

If you still need to choose an installation method or workflow, start with the
[Overview](../overview).

## 1. Confirm the CLI

```bash
bin/cbicall --help
bin/cbicall --version
```

You should see the command help and the installed CBIcall version.

## 2. Run the WES Example Test

From the repository root:

```bash
bin/cbicall test --wes-bash -t 1
```

This runs the required bundled Bash WES example and compares the generated VCF
with the shipped reference output.

:::tip[Need deeper checks?]
Use [Integration Tests](integration-tests) for bundled WES/mtDNA test details,
[Configuration Reference](../help/configuration-reference) for
`validate-parameters`, and [Resource Validation](resource-validation) for checking
the selected resource entry.
:::

## 3. Optional: Run the mtDNA Test

Run this after the WES test works:

```bash
bin/cbicall test --mit-bash -t 1
```

:::info[Architecture]
The mtDNA workflow uses MToolBox and is x86_64-only. If you are on ARM / aarch64, run WES/WGS workflows there but move mtDNA runs to an x86_64 host.
:::

## 4. Run With Your Own YAML

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
| Run real WES/WGS data | [End-to-end Example: WES](end-to-end-example-wes) |
| Run mtDNA analysis | [End-to-end Example: mtDNA](end-to-end-example-mit) |
| Check reproducibility | [Run Comparison](run-comparison) |
| Understand generated files | [Outputs](../help/outputs) |
