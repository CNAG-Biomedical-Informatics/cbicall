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

For external nf-core workflows, workflow-specific inputs are passed in
`nfcore_parameters`. For example, Sarek expects its samplesheet under
`nfcore_parameters.input`. CBIcall launches the registered Nextflow workflow and
leaves its native output directory layout unchanged.

## CLI Synopsis

```text
cbicall run -p <parameters_file.yaml> -t <n_threads> [options]
```

Use `run` for normal analysis execution. Other subcommands are documented on
their dedicated pages.

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

:::tip[Thread choice]
For most WES/WGS runs, start with **4 threads per task**. See [Performance](../help/performance) for the benchmark and scaling guidance.
:::

## Minimal YAML Examples

### WES Single-Sample

```yaml
mode:            single
pipeline:        wes
workflow_backend: bash
software_stack:    gatk-4.6
input_dir:       CNAG999_exome/CNAG99901P_ex
genome:          b37
```

### WES Cohort

```yaml
mode:            cohort
pipeline:        wes
workflow_backend: bash
software_stack:    gatk-4.6
genome:          b37
sample_map:      ./sample_map.tsv
```

### WGS Single-Sample With Nextflow

```yaml
mode:            single
pipeline:        wgs
workflow_backend: nextflow
software_stack:    gatk-4.6
input_dir:       CNAG999_genome/CNAG99901P_wg
genome:          hg38
```

### WES Single-Sample With Cromwell

The checked-in example is `examples/input/wes_cromwell.yaml`. Set
`CROMWELL_JAR`, put `cromwell` on `PATH`, or put `cromwell*.jar` on `PATH` before running it.

```yaml
mode:            single
pipeline:        wes
workflow_backend: cromwell
software_stack:    gatk-4.6
input_dir:       CNAG999_exome/CNAG99901P_ex
genome:          b37
```

### mtDNA Single-Sample

```yaml
mode:            single
pipeline:        mit
workflow_backend: bash
software_stack:    gatk-3.5
input_dir:       CNAG999_exome/CNAG99901P_ex
```

### nf-core/Sarek Through CBIcall

The checked-in example is `examples/input/nf-core-sarek.yaml`.

```yaml
mode:             cohort
pipeline:         sarek
workflow_backend:  nextflow
workflow_provider: nf-core
resource:         nf-core-sarek-managed-resources-v1
nfcore_profile: singularity
# nfcore_singularity_cache_dir: nxf-singularity-cache
nfcore_parameters:
  input: sarek_samplesheet.csv
  genome: GATK.GRCh38
  tools: haplotypecaller
  skip_tools: haplotypecaller_filter
  wes: true
  intervals: ../../workflows/nextflow/nf-core/sarek/grch38_chr22_test.bed
  max_memory: 30.GB
```

### nf-core/demo Smoke Test

The checked-in example is `examples/input/nf-core-demo.yaml`.

```yaml
mode:             single
pipeline:         demo
workflow_backend:  nextflow
workflow_provider: nf-core
resource:         nf-core-demo-managed-resources-v1
nfcore_profile: test,singularity
```

This example uses nf-core's built-in test profile with Singularity/Apptainer,
which is the recommended profile on HPC. On an x86_64 Docker workstation,
`test,docker` is also possible.

For workstation and cluster runs with nf-core workflows, see
[nf-core Workflows](../backends/nf-core).

## Backend-Specific Parameters

Use backend-specific parameter blocks for values that belong to an execution backend rather than to CBIcall's global analysis contract.

For targeted Snakemake execution, set a Snakemake target:

```yaml
workflow_backend: snakemake
snakemake_parameters:
  target: call_variants
```

For native CBIcall Nextflow workflows, pass extra workflow parameters with:

```yaml
workflow_backend: nextflow
nextflow_parameters:
  publish_debug_files: true
```

For native CBIcall Cromwell workflows, pass advanced WDL inputs with
`cromwell_parameters`. CBIcall still owns tool paths, reference paths, sample
identity, genome, pipeline, and thread count.

CBIcall still owns core values such as `pipeline`, `genome`, `threads`, helper script paths, and cohort workspace names, and refuses backend parameter blocks that try to override them.

## Outputs and Logs

Each run creates a directory containing:

- the main workflow log
- `log.json` with CLI arguments, resolved configuration, and YAML parameters
- workflow outputs such as VCFs, BAMs, statistics, and browser reports

External workflows can keep their own native output layout. For nf-core/Sarek,
CBIcall writes `cbicall_external_nextflow.params.yaml` and Sarek writes under
`sarek/` inside the CBIcall run directory.

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
| Source installation | Python 3.8+, Java 8, plus Snakemake, Nextflow, or Cromwell when using those workflow backends |

## Next Steps

- First run: [Quickstart](quickstart)
- Full YAML reference: [Configuration Reference](../help/configuration-reference)
- Output files: [Outputs](../help/outputs)
