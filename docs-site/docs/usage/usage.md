# General Usage

CBIcall is normally run with a YAML parameters file and a thread count:

```bash
bin/cbicall run -p parameters.yaml -t 4
```

Use this page for command syntax and common execution patterns. For YAML keys and supported combinations, see [Configuration Reference](../help/configuration-reference).

:::tip[Typical workflow]
Choose or edit a YAML file, run `bin/cbicall run -p ... -t ...`, then inspect the generated run directory and `log.json`.
:::

## Input Model

For WES/WGS single-sample runs, `input_dir` points to a sample directory containing paired FASTQ files:

```text
CNAG999_exome/CNAG99901P_ex/
  *_R1_*.fastq.gz
  *_R2_*.fastq.gz
```

For WES/WGS cohort runs, `sample_map` points to a TSV containing sample IDs and gVCF paths.

For mtDNA runs, CBIcall expects BAM files from previous WES/WGS single-sample runs. mtDNA workflows do not start from FASTQ files.

## CLI Synopsis

```text
cbicall run -p <parameters_file.yaml> -t <n_threads> [options]
```

Use `run` for normal analysis execution. The other subcommands perform preflight checks, validation, or bundled integration tests.

| Argument | Meaning |
| --- | --- |
| `-p`, `--param` | YAML parameters file. |
| `-t`, `--threads` | Number of CPUs/cores/threads passed to the workflow. |

| Option | Meaning |
| --- | --- |
| `-verbose` | Enable verbose output. |
| `-debug`, `--debug` | Debugging level. |
| `-nc`, `--no-color` | Disable ANSI colors explicitly. Colors are automatically disabled when output is redirected. |
| `-v`, `--version` | Show version information. |
| `-h`, `--help` | Show brief help. |
| `-man` | Show full command-line documentation. |

## Common Commands

```bash
bin/cbicall run -p wes_single.yaml -t 4
bin/cbicall run -p wes_cohort.yaml -t 8 -verbose
bin/cbicall run -p mit_single.yaml -t 4 > run.log 2>&1
nohup bin/cbicall run -p parameters.yaml -t 4 > run.log 2>&1 &
```

| Pattern | Use when |
| --- | --- |
| `bin/cbicall run -p wes_single.yaml -t 4` | Normal foreground run. |
| `-verbose` | You want more CLI output while the workflow starts. |
| `> run.log 2>&1` | You are saving terminal output to a file. ANSI colors are disabled automatically. |
| `nohup ... &` | You need a simple long-running background job outside a scheduler. |

## Integration Tests

Use the built-in test command from the repository root. It runs the bundled example workflow and compares outputs against the shipped reference files.

These tests are intended as framework-level integration validation: they exercise
the normal CLI, parameter YAML resolution, workflow registry dispatch,
pipeline-version provenance, resource-bundle reporting, workflow execution, and
deterministic output comparison. They do not replace biological or clinical
validation of the underlying variant-calling methods.

```bash
bin/cbicall test --wes -t 1
bin/cbicall test --mit -t 1
bin/cbicall test --all -t 1
bin/cbicall validate-resources
bin/cbicall compare-runs run_a/ run_b/ run_c/ --output compare-report.txt --html compare-report.html
```

| Command | Use |
| --- | --- |
| `bin/cbicall test --wes -t 1` | Fast WES integration test. Run this first. |
| `bin/cbicall test --mit -t 1` | mtDNA integration test after the WES path is working. |
| `bin/cbicall test --all -t 1` | Run all bundled integration examples. |
| `bin/cbicall validate-resources` | Validate the resource catalog and its workflow compatibility keys. Add `--bundle <key>` to check one bundle entry. |
| `bin/cbicall compare-runs run_a/ run_b/ run_c/ --output compare-report.txt --html compare-report.html` | Compare run reports for workflow, resource, and output fingerprint differences. See [Run Comparison](run-comparison). |

:::tip[Thread choice]
For most WES/WGS runs, start with **4 threads per task**. See [Performance](../help/performance) for the benchmark and scaling guidance.
:::

## Minimal YAML Examples

### WES Single-Sample

```yaml
mode:            single
pipeline:        wes
workflow_engine: bash
gatk_version:    gatk-4.6
input_dir:       CNAG999_exome/CNAG99901P_ex
genome:          b37
```

### WES Cohort

```yaml
mode:            cohort
pipeline:        wes
workflow_engine: bash
gatk_version:    gatk-4.6
genome:          b37
sample_map:      ./sample_map.tsv
```

### mtDNA Single-Sample

```yaml
mode:            single
pipeline:        mit
workflow_engine: bash
gatk_version:    gatk-3.5
input_dir:       CNAG999_exome/CNAG99901P_ex
```

## Partial Runs

Partial runs are intended for targeted Snakemake execution and restarts. Leave `workflow_rule` unset for a normal full run.

```yaml
workflow_engine: snakemake
workflow_rule: call_variants
allow_partial_run: true
```

Behavior:

- if `workflow_rule` is unset, the full workflow runs
- if `workflow_rule` is set, the run is marked as partial in metadata and CLI output
- if `workflow_rule` is set without `allow_partial_run: true`, CBIcall refuses to start

## Outputs and Logs

Each run creates a directory containing:

- the main workflow log
- `log.json` with CLI arguments, resolved configuration, and YAML parameters
- workflow outputs such as VCFs, BAMs, statistics, and browser reports

See [Outputs](../help/outputs) for the file reference.

:::info[Reproducibility]
Keep `log.json` with the workflow outputs. It records the CLI arguments, resolved configuration, and parameters used for the run.
:::

## Recommended Specifications

CBIcall is designed for multi-core Linux desktop, workstation, server, and HPC environments.

| Resource | Recommendation |
| --- | --- |
| CPU | At least 4 cores |
| RAM | At least 8 GB for single-sample runs; more for cohort joint genotyping |
| Disk | At least 250 GB, depending on input size and whether BAM cleanup is enabled |
| Source installation | Python 3.8+, Java 8, Snakemake when using Snakemake workflows |

## Next Steps

- First run: [Quickstart](quickstart)
- Decide workflow choices: [Choose Your Path](choose-your-path)
- Full YAML reference: [Configuration Reference](../help/configuration-reference)
- Output files: [Outputs](../help/outputs)
