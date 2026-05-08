# Quickstart

Use this page when you want to verify that CBIcall runs correctly on the shipped example data.

:::tip[Fast path]
Run the WES smoke test first. If it passes, the CLI, example layout, reference paths, and core GATK workflow are working together.
:::

If you still need to choose an installation method or workflow, start with [Choose Your Path](choose-your-path).

## 1. Confirm the CLI

```bash
bin/cbicall --help
bin/cbicall --version
```

You should see the command help and the installed CBIcall version.

## 2. Run the WES Smoke Test

```bash
cd examples/input
./run_tests.sh --wes
```

This test:

| Check | What it confirms |
| --- | --- |
| Starts `cbicall` | The CLI and Python package are usable. |
| Runs the bundled WES example | The example input layout is valid. |
| Compares the output VCF | The result matches the shipped reference output. |

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

## 3. Run the mtDNA Smoke Test

Run this after the WES test works. The mtDNA workflow depends on the WES/WGS-style project structure and consumes BAMs from previous WES/WGS runs.

```bash
cd examples/input
./run_tests.sh --mit
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

:::info[Architecture]
The mtDNA workflow uses MToolBox and is x86_64-only. If you are on ARM / aarch64, run WES/WGS workflows there but move mtDNA runs to an x86_64 host.
:::

## 4. Run With Your Own YAML

Once the smoke tests work, use the normal invocation:

```bash
bin/cbicall -p parameters.yaml -t 4
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

