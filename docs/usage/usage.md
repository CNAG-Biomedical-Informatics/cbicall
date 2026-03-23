# General Usage

CBIcall execution requires:

## Input Files

A folder containing paired-end FASTQ files, for example:

```text
CNAG999_exome/CNAG99901P_ex//*{R1,R2}*fastq.gz
```

Example input data is available under `examples/input/`.

## Parameters File

CBIcall is configured through a YAML parameters file.

Below are the main parameters and their default values.

## Essential Parameters

```yaml
mode:              single
pipeline:          wes
input_dir:         null
sample_map:        null
workflow_engine:   bash
gatk_version:      gatk-4.6
genome:            b37
cleanup_bam:       false
workflow_rule:     null
allow_partial_run: false
project_dir:       cbicall
```

## Optional Parameters

```yaml
organism:    Homo sapiens
technology:  Illumina HiSeq
```

CBIcall creates a dedicated project directory to store analysis outputs. Multiple runs can coexist without modifying the original inputs.

## Parameter Reference

### `cleanup_bam`

Set to `true` to delete `01_bam/*.{bam,bai}` after the workflow completes.

### `gatk_version`

Supported values:

- `gatk-3.5`
- `gatk-4.6`

### `genome`

Supported values:

- `b37` (default for non-mtDNA runs)
- `hg38`
- `rsrs` (forced for `pipeline: mit`)

### `mode`

Supported values:

- `single`
- `cohort`

### `pipeline`

Supported values:

- `wes`
- `wgs`
- `mit`

For `cohort` analysis, complete the corresponding `single` analysis for each sample first.

### `project_dir`

Prefix used for the generated run directory, for example `cancer_sample_001`.
It can also contain a path such as `foo/cancer_sample_001`.

When `input_dir` is defined, the final run directory is created below `input_dir`, with a unique run identifier appended automatically.

### `input_dir`

Path, relative or absolute, to the directory containing FASTQ inputs for analysis.

Example:

```text
examples/input/CNAG999_exome/CNAG99901P_ex
```

### `sample_map`

Cohort-mode parameter pointing to the TSV file containing sample identifiers and GVCF paths.

Example:

- [examples/input/sample_map.tsv](../../examples/input/sample_map.tsv)

### `workflow_engine`

Supported workflow engines:

- `bash`
- `snakemake`

### `workflow_rule`

Optional workflow target used to start a partial run. Leave it unset for a normal full workflow execution.

Partial runs are currently supported for the `snakemake` workflow engine.

### `allow_partial_run`

Safety switch for partial runs. Set it to `true` together with `workflow_rule` to explicitly allow a partial workflow start.

If `workflow_rule` is set but `allow_partial_run` is not `true`, CBIcall refuses to start the run.

## CLI Synopsis

```text
cbicall -p <parameters_file.yaml> -t <n_threads> [options]
```

Arguments:

- `-p`, `--param`: Parameters input file (YAML)
- `-t`, `--threads`: Number of CPUs/Cores/Threads

Options:

- `-debug`, `--debug`: Debugging level
- `-h`, `--help`: Brief help message
- `-man`: Full documentation
- `-v`, `--version`: Show version information
- `-verbose`: Enable verbose output
- `-nc`, `--no-color`: Disable ANSI colors in STDOUT

## Example Commands

```bash
bin/cbicall -p param_file.yaml -t 8
bin/cbicall -p param_file.yaml -t 4 -verbose
bin/cbicall -p param_file.yaml -t 16 > log 2>&1
python3 -u ../../bin/cbicall -p param.yaml -t 8 --no-color > log
$path_to_cbicall/bin/cbicall -p param_file.yaml -t 8 -debug 5
nohup bin/cbicall -p param_file.yaml -t 4 &
```

## Partial Runs

Partial runs are intended for workflow-layer restarts or targeted execution, without changing the output directory naming scheme.

Example:

```yaml
workflow_engine: snakemake
workflow_rule: call_variants
allow_partial_run: true
```

Behavior:

- if `workflow_rule` is unset, the full workflow runs
- if `workflow_rule` is set, the run is marked as partial in metadata and CLI output
- if `workflow_rule` is set without `allow_partial_run: true`, the run is rejected

## Outputs and Logging

Each run creates a project directory containing:

- the workflow log file
- `log.json` with CLI arguments, resolved configuration, and YAML parameters
- workflow-generated outputs

The project directory name is not changed when `workflow_rule` is used.

## Recommended Specifications

CBIcall is optimized for multi-core Linux desktop, workstation, or server environments.

- Works on amd64 and arm64 architectures, including Apple Silicon
- Debian-based distributions are preferred, though others may work
- At least 8 GB RAM
- At least 4 CPU cores
- At least 250 GB disk space

## Additional Requirements for Source Installation

These are required only when installing CBIcall from source:

- Python 3.8 or newer
- Java 8
- Snakemake
