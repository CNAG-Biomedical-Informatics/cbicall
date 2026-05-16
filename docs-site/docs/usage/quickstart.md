# Quickstart

Use this page when you want to verify that CBIcall runs correctly on the shipped example data.

:::tip[Fast path]
Run the WES integration test first. If it passes, the CLI, example layout, reference paths, and core GATK workflow are working together.
:::

If you still need to choose an installation method or workflow, start with [Choose Your Path](choose-your-path).

## 1. Confirm the CLI

```bash
bin/cbicall --help
bin/cbicall --version
```

You should see the command help and the installed CBIcall version.

## 2. Run the Minimal Reproducibility Check

From the repository root, run the bundled WES check:

```bash
bin/cbicall validate-registry
bin/cbicall validate-resources
bin/cbicall doctor -p examples/input/param.yaml
bin/cbicall test --wes -t 1
```

This is a small framework-level reproducibility check on shipped example data.
It validates the workflow registry and resource catalog, resolves the
user-facing YAML without launching GATK, then runs the example workflow and
compares the generated VCF against the shipped reference output.

This test:

| Check | What it confirms |
| --- | --- |
| `validate-registry` | The workflow registry conforms to its JSON Schema. |
| `validate-resources` | The resource catalog is well formed and compatible workflow keys exist in the registry. |
| `doctor` | The YAML resolves to a declared workflow, profile, pipeline implementation version, and selected bundle. |
| `test --wes` | The CLI can execute the bundled WES workflow and reproduce the shipped reference VCF under deterministic comparison rules. |

When it succeeds, the run directory looks like this:

```text
CNAG999_exome/CNAG99901P_ex/cbicall_bash_wes_single_b37_gatk-4.6_*/
  01_bam/
  02_varcall/
  03_stats/
  logs/
```

The main file to inspect is:

```text
02_varcall/CNAG99901P.hc.QC.vcf.gz
```

The test command also prints the exact run directory, workflow log, `run-report.json`, launcher log, and output file used for comparison.

## 3. Run the mtDNA Integration Test

Run this after the WES test works. The mtDNA workflow depends on the WES/WGS-style project structure and consumes BAMs from previous WES/WGS runs.

```bash
bin/cbicall test --mit -t 1
```

When it succeeds, the run directory looks like this:

```text
CNAG999_exome/CNAG99901P_ex/cbicall_bash_mit_single_rsrs_gatk-3.5_*/
  01_mtoolbox/
  02_browser/
```

Start with:

```text
02_browser/README.txt
```

The final summary reports whether the mtDNA text and JSON outputs matched the shipped references, plus the exact run directory used for the comparison.

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

For most WES/WGS runs, start with 4 threads and adjust after checking [Performance](../help/performance).

## Next Steps

| Goal | Page |
| --- | --- |
| Pick the right workflow | [Choose Your Path](choose-your-path) |
| Run real WES/WGS data | [End-to-end Example: WES](end-to-end-example-wes) |
| Run mtDNA analysis | [End-to-end Example: mtDNA](end-to-end-example-mit) |
| Understand generated files | [Outputs](../help/outputs) |
| Edit the YAML safely | [Configuration Reference](../help/configuration-reference) |
